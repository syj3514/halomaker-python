#!/usr/bin/env python
import argparse
from pathlib import Path

from gasmaker import GasMaker
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


def main():
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
        "--rur-path",
        default=None,
        help="path to a rur checkout (default: use installed rur, or $RUR_PATH)",
    )
    args = parser.parse_args()
    if args.output is None:
        args.output = Path(f"gas_bricks{args.iout:05d}.h5")

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
    maker = GasMaker(
        args.catalog,
        reader,
        radius_field=args.radius_field,
        padding=args.padding,
        overlap_depth=args.overlap_depth,
        overlap_tolerance=args.overlap_tolerance,
        overlap_threshold=args.overlap_threshold,
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
            overwrite=args.overwrite,
            read_grav=args.read_grav,
            nthread=args.nthread,
            stop_after_roots=args.stop_after_roots,
            simulate_crash_after_summary_root=(
                args.simulate_crash_after_summary_root
            ),
        )
    finally:
        reader.close()

    print(f"requested_roots={','.join(str(item) for item in status['requested'])}")
    print(f"processed_roots={','.join(str(item) for item in status['processed'])}")
    print(f"skipped_roots={','.join(str(item) for item in status['skipped'])}")
    print(f"remaining_roots={','.join(str(item) for item in status['remaining'])}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
