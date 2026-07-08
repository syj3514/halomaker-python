#!/usr/bin/env python
import argparse
import os
import sys
from types import SimpleNamespace
from pathlib import Path

# HDF5 file locking fails (errno 11, "Resource temporarily unavailable") on some
# network/parallel filesystems (Lustre/NFS). GasMaker is the sole writer of its
# output file, so disabling the lock is safe; must be set before any h5py import.
# Export HDF5_USE_FILE_LOCKING to override.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

_ROOT = Path(__file__).resolve().parent
_SRC = str(_ROOT / "src")
if _SRC not in sys.path:
    sys.path.append(_SRC)

from gasmaker import GasMaker
from gasmaker.progress import MODES as PROGRESS_MODES, Progress
from gasmaker.config import (
    INPUTFILES_FILE,
    PARAM_FILE,
    read_inputfiles,
    read_params,
    validate_unique_outputs,
)
# The rur adapter is the default snapshot reader. It is imported lazily (and
# `rur` itself only when a reader is constructed) so the GasMaker core does not
# depend on rur. To use a different simulation/format, implement the
# gasmaker.readers.base.CellReader interface and pass it to GasMaker(...).
from gasmaker.readers.rur import RurCellReader


def _parse_root_ids(value):
    root_ids = []
    for item in value.split(","):
        item = item.strip()
        if item:
            root_ids.append(int(item))
    if not root_ids:
        raise argparse.ArgumentTypeError("expected at least one root id")
    return root_ids


def _select_root_ids(maker, positional_root_id, roots, root_ids, max_roots):
    if root_ids is not None:
        selected = root_ids
    elif positional_root_id is not None:
        selected = [positional_root_id]
    elif roots == "all":
        rows = maker.catalog.root_rows()
        selected = [int(root_id) for root_id in maker.catalog.halos["id"][rows]]
    else:
        raise ValueError(f"Unsupported root selection: {roots}")
    if max_roots is not None:
        selected = selected[:max_roots]
    return selected


def _build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("catalog")
    parser.add_argument("repo")
    parser.add_argument("iout", type=int)
    parser.add_argument("root_id", type=int, nargs="?")
    parser.add_argument("--mode", default="nh2")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="output HDF5 path (default: gas_bricks{iout:05d}.h5)",
    )
    parser.add_argument("--roots", choices=("all",), default="all")
    parser.add_argument("--root-ids", type=_parse_root_ids)
    parser.add_argument("--max-roots", type=int)
    parser.add_argument("--overwrite", "--force", action="store_true")
    parser.add_argument(
        "--stop-after-roots",
        type=int,
        help="stop after processing this many new roots; useful for restart tests",
    )
    parser.add_argument(
        "--simulate-crash-after-summary-root",
        type=int,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--radius-field", default="r")
    parser.add_argument("--padding", type=float, default=1.1)
    parser.add_argument("--nthread", type=int, default=1)
    parser.add_argument("--overlap-depth", type=int, default=8)
    parser.add_argument("--overlap-tolerance", type=float, default=0.05)
    parser.add_argument("--overlap-threshold", type=float, default=0.1)
    parser.add_argument("--read-grav", action="store_true")
    parser.add_argument(
        "--progress",
        choices=PROGRESS_MODES,
        default="auto",
        help="run progress display: auto (tty->bar, else plain), bar, plain, quiet",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1,
        help="in plain mode, print one line every N completed roots",
    )
    parser.add_argument(
        "--rur-path",
        default=None,
        help="path to a rur checkout (default: use installed rur, or $RUR_PATH)",
    )
    return parser


def _describe_root_policy(args):
    if args.root_id is not None:
        policy = f"single root {args.root_id}"
    elif args.root_ids:
        policy = f"root_ids={','.join(str(item) for item in args.root_ids)}"
    else:
        policy = "all level-1 roots"
    if args.max_roots is not None:
        policy += f" (max_roots={args.max_roots})"
    return policy


def _run_job(args):
    import time

    if args.output is None:
        args.output = Path(f"gas_bricks{args.iout:05d}.h5")

    # Config-mode provenance paths (set by _run_config_mode; "" in CLI mode).
    config_param_path = getattr(args, "_config_param_path", "")
    config_inputfiles_path = getattr(args, "_config_inputfiles_path", "")

    progress = Progress(
        mode=getattr(args, "progress", "auto"),
        every=getattr(args, "progress_every", 1),
    )
    progress.banner(
        catalog=args.catalog,
        repo=args.repo,
        iout=args.iout,
        mode=args.mode,
        output=args.output,
        nthread=args.nthread,
        policy=_describe_root_policy(args),
    )

    stage_t0 = time.perf_counter()
    try:
        reader = RurCellReader(
            args.repo, args.iout, mode=args.mode, rur_path=args.rur_path
        )
    except ImportError as exc:
        raise SystemExit(
            f"Could not import the rur RAMSES reader from '{args.rur_path}'. "
            f"Install rur or set --rur-path, or supply another reader that "
            f"implements gasmaker.readers.base.CellReader. ({exc})"
        )
    progress.stage(
        f"snapshot reader ready (rur mode={args.mode}, iout={args.iout}) — "
        f"{time.perf_counter() - stage_t0:.1f}s"
    )

    stage_t0 = time.perf_counter()
    maker = GasMaker(
        args.catalog,
        reader,
        radius_field=args.radius_field,
        padding=args.padding,
        overlap_depth=args.overlap_depth,
        overlap_tolerance=args.overlap_tolerance,
        overlap_threshold=args.overlap_threshold,
    )
    n_halos = len(maker.catalog.halos)
    n_roots = len(maker.catalog.root_rows())
    progress.stage(
        f"catalog loaded: {n_halos} halos, {n_roots} level-1 roots — "
        f"{time.perf_counter() - stage_t0:.1f}s"
    )
    try:
        root_ids = _select_root_ids(
            maker,
            args.root_id,
            args.roots,
            args.root_ids,
            args.max_roots,
        )
        status = maker.process_roots(
            root_ids,
            args.output,
            config_param_path=config_param_path,
            config_inputfiles_path=config_inputfiles_path,
            overwrite=args.overwrite,
            read_grav=args.read_grav,
            nthread=args.nthread,
            stop_after_roots=args.stop_after_roots,
            simulate_crash_after_summary_root=(
                args.simulate_crash_after_summary_root
            ),
            progress=progress,
        )
    finally:
        reader.close()

    progress.finish(status)
    print(f"requested_roots={','.join(str(item) for item in status['requested'])}")
    print(f"processed_roots={','.join(str(item) for item in status['processed'])}")
    print(f"skipped_roots={','.join(str(item) for item in status['skipped'])}")
    print(f"remaining_roots={','.join(str(item) for item in status['remaining'])}")
    print(f"output={args.output}")


def _run_config_mode(parser):
    if not (Path(PARAM_FILE).exists() and Path(INPUTFILES_FILE).exists()):
        parser.error(
            f"missing positional arguments or config files; run without arguments "
            f"only when {PARAM_FILE} and {INPUTFILES_FILE} exist in the current "
            f"directory"
        )

    try:
        params = read_params(PARAM_FILE)
        jobs = read_inputfiles(INPUTFILES_FILE)
        validate_unique_outputs(jobs)
    except (OSError, ValueError) as exc:
        raise SystemExit(str(exc))

    defaults = parser.parse_args(["__catalog__", "__repo__", "0"])
    for job in jobs:
        args = SimpleNamespace(**vars(defaults))
        # Pass absolute paths of the config files used, for header provenance.
        args._config_param_path = str(Path(PARAM_FILE).resolve())
        args._config_inputfiles_path = str(Path(INPUTFILES_FILE).resolve())
        args.catalog = job.catalog
        args.repo = job.repo
        args.iout = job.iout
        args.root_id = None
        args.output = job.output
        for key, value in params.items():
            setattr(args, key, value)
        _run_job(args)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    parser = _build_parser()
    if argv:
        _run_job(parser.parse_args(argv))
    else:
        _run_config_mode(parser)


if __name__ == "__main__":
    main()
