"""Reproduce independent night-block diagnostics from retained public NEID JSON."""
from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from research.neid_variability import load_accepted  # noqa: E402
from research.night_block_audit import (  # noqa: E402
    quality_cut_decomposition, leave_one_night_out, night_cluster_bootstrap,
)


def main():
    path = ROOT / "public/data/neid-solar-2459700-2459713.json"
    data = load_accepted(path)
    result = {
        "provenance": str(path.relative_to(ROOT)),
        "accepted_exposures": int(data["bjd"].size),
        "analyses": {
            "snr_night_decomposition": quality_cut_decomposition(data),
            "leave_one_night_out": leave_one_night_out(data),
            "night_block_bootstrap": night_cluster_bootstrap(data),
        },
        "limitations": "Four observed nights; block-bootstrap percentiles are exploratory, "
        "not calibrated confidence intervals. Quality cuts are nonrandom and not activity corrections.",
    }
    output = ROOT / "research/night-block-audit-results.json"
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = result["analyses"]["snr_night_decomposition"]
    valid = [r for r in rows if r["valid"]]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot([r["snr_cut"] for r in valid],
            [r["total_rms_mps"] for r in valid], "o-", label="Weighted total RMS")
    ax.plot([r["snr_cut"] for r in valid],
            [r["within_night_rms_mps"] for r in valid], "s-", label="Within-night RMS")
    ax.set(xlabel="S/N minimum", ylabel="Weighted RV RMS (m/s)",
           title="NEID quality-cut / night-composition sensitivity")
    ax.legend()
    fig.tight_layout()
    fig.savefig(ROOT / "research/night-block-audit.png", dpi=150)
    plt.close(fig)
    print(f"Night-block audit: {result['accepted_exposures']} exposures; {output}")


if __name__ == "__main__":
    main()
