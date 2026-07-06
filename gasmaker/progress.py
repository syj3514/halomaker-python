"""Run-wide progress reporting for GasMaker.

A single small object owns every user-facing progress line, so the pipeline
and CLI stay free of presentation logic. Output is strictly cosmetic: nothing
here touches physics, schema, or the HDF5 output.

All progress output goes to **stderr** (adopted from the Codex TASK-23 entry):
stdout carries only the final machine-readable summary lines, so
``python GasMaker.py > result.log`` keeps the log parseable while progress
stays visible on the terminal, and ``> run.log 2>&1`` captures the full
timeline.

Modes
-----
- ``auto``  (default): ``bar`` when stderr is an interactive terminal,
  otherwise ``plain`` (so ``> run.log 2>&1`` never records ``\r`` spam).
- ``bar``   : tqdm progress bar over the root loop; stage lines are routed
  through ``tqdm.write`` so they do not break the bar.
- ``plain`` : one timestamped line per stage and per root (or every
  ``every`` roots), with elapsed / ETA — designed for log files.
- ``quiet`` : stage banner and final summary only.
"""

from __future__ import annotations

import sys
import time
from datetime import datetime

MODES = ("auto", "bar", "plain", "quiet")


def _hms(seconds):
    seconds = max(0, int(round(seconds)))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m:02d}:{s:02d}"


class Progress:
    def __init__(self, mode="auto", every=1, stream=None):
        if mode not in MODES:
            raise ValueError(f"progress expects one of {MODES}, got {mode!r}")
        self.stream = stream if stream is not None else sys.stderr
        if mode == "auto":
            mode = "bar" if sys.stderr.isatty() else "plain"
        self.mode = mode
        self.every = max(1, int(every))
        self.t0 = time.perf_counter()
        self._bar = None
        self._todo = 0
        self._done_now = 0
        self._read_total = 0.0
        self._compute_total = 0.0
        self._slowest = None  # (seconds, root_id)

    # -- low-level ---------------------------------------------------------
    def _emit(self, text):
        if self._bar is not None:
            from tqdm import tqdm
            tqdm.write(text, file=self.stream)
        else:
            print(text, file=self.stream, flush=True)

    def _line(self, text):
        self._emit(f"[{datetime.now():%H:%M:%S}] {text}")

    # -- stages ------------------------------------------------------------
    def banner(self, *, catalog, repo, iout, mode, output, nthread, policy):
        if self.mode == "quiet":
            self._line(f"GasMaker start: iout={iout} output={output}")
            return
        bar = "=" * 66
        self._emit(bar)
        self._emit(f" GasMaker  catalog : {catalog}")
        self._emit(f"           repo    : {repo}  (mode={mode}, iout={iout})")
        self._emit(f"           output  : {output}   nthread={nthread}")
        self._emit(f"           roots   : {policy}")
        self._emit(bar)

    def stage(self, text):
        if self.mode == "quiet":
            return
        self._line(text)

    # -- root loop ---------------------------------------------------------
    def start_roots(self, requested, already_done, todo):
        self._todo = todo
        self._loop_t0 = time.perf_counter()
        self.stage(
            f"roots: {requested} requested · {already_done} already done "
            f"(resume-skip) · {todo} to process"
        )
        if self.mode == "bar" and todo:
            from tqdm import tqdm
            self._bar = tqdm(
                total=todo, unit="root", desc="roots", ncols=100,
                file=sys.stderr, dynamic_ncols=False,
            )

    def root_done(self, metrics):
        self._done_now += 1
        self._read_total += metrics.read_seconds
        self._compute_total += metrics.compute_seconds
        total = metrics.read_seconds + metrics.compute_seconds
        if self._slowest is None or total > self._slowest[0]:
            self._slowest = (total, metrics.root_id)

        if self._bar is not None:
            self._bar.set_postfix_str(
                f"root={metrics.root_id} halos={metrics.descendant_count} "
                f"read={metrics.read_seconds:.1f}s comp={metrics.compute_seconds:.1f}s",
                refresh=False,
            )
            self._bar.update(1)
            return
        if self.mode == "quiet":
            return
        if self._done_now % self.every and self._done_now != self._todo:
            return
        elapsed = time.perf_counter() - self._loop_t0
        rate = elapsed / self._done_now
        eta = rate * (self._todo - self._done_now)
        self._line(
            f"[{self._done_now}/{self._todo}] root={metrics.root_id} "
            f"halos={metrics.descendant_count} cells={metrics.cells_read} "
            f"read={metrics.read_seconds:.1f}s comp={metrics.compute_seconds:.1f}s "
            f"elapsed={_hms(elapsed)} eta={_hms(eta)}"
        )

    # -- end ---------------------------------------------------------------
    def finish(self, status):
        if self._bar is not None:
            self._bar.close()
            self._bar = None
        wall = time.perf_counter() - self.t0
        self._line(
            f"done: {len(status['processed'])} processed · "
            f"{len(status['skipped'])} skipped (resume) · "
            f"{len(status['remaining'])} remaining · wall {_hms(wall)}"
        )
        if self._done_now:
            slow_s, slow_id = self._slowest
            self._line(
                f"totals: read {self._read_total:.1f}s · "
                f"compute {self._compute_total:.1f}s · "
                f"slowest root {slow_id} ({slow_s:.1f}s)"
            )
