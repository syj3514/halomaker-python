"""Pluggable snapshot readers for GasMaker.

`base` defines the reader contract; concrete adapters (e.g. `rur`) are imported
lazily so the GasMaker core does not hard-depend on any one I/O library.
"""
from .base import CellReader, ReadResult

__all__ = ["CellReader", "ReadResult"]
