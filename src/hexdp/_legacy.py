"""Support for code migrated from experimentell/: the repository root, which the old modules use
for their default input and output paths (data/, output/). Only meaningful in a source checkout
(editable install); new code takes paths as arguments."""

from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
