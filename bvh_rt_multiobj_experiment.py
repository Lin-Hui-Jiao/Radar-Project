#!/usr/bin/env python3
"""Compatibility entry for tools/bvh_rt_multiobj_experiment.py."""

from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    script = Path(__file__).resolve().parent / "tools" / "bvh_rt_multiobj_experiment.py"
    runpy.run_path(str(script), run_name="__main__")
