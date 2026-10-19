"""Reproducible diagnostics using the existing committed NEID solar metadata."""
from __future__ import annotations
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from research.neid_variability import load_accepted  # noqa: E402
from research.neid_robustness import (  # noqa: E402
    profile_intrinsic_jitter, night_cluster_bootstrap, leave_one_night_out,
)


def main() -> None:
    source = ROOT / "public/data/neid-solar-2459700-2459713.json"
    data = load_accepted(source)
    result = {
        "data_kind": "existing public NEID Level-2 solar RV metadata",
        "source_file": str(source.relative_to(ROOT)),
        "selection": "existing analysisPass filter",
        "profiled_white_jitter": profile_intrinsic_jitter(data),
        "night_cluster_bootstrap": night_cluster_bootstrap(data),
        "leave_one_night_out": leave_one_night_out(data),
        "limitations": [
            "BJD-day groups are not guaranteed to equal local observing nights.",
            "White jitter combines stellar, instrument and unmodelled correlated terms.",
            "The short window cannot establish solar rotation or planetary signals.",
            "The night bootstrap does not account for long-term nonstationarity.",
        ],
    }
    target = ROOT / "research/neid-robustness-summary.json"
    target.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Wrote {target} for {len(data['bjd'])} accepted exposures")


if __name__ == "__main__":
    main()
