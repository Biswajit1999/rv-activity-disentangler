"""Query a reproducible public NEID Level-2 solar RV interval through TAP."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
from io import StringIO
import json
from pathlib import Path

import numpy as np
import requests

TAP_URL = "https://neid.ipac.caltech.edu/TAP/sync"
QUERY = """SELECT TOP 2000 obsdate,ccfjdsum,ccfrvmod,dvrms,extsnr,sunagl,flagged,l2filename,l2checksum,swversion
FROM neidsolarl2 WHERE ccfjdsum BETWEEN 2459700 AND 2459713 ORDER BY ccfjdsum"""


def main() -> None:
    response = requests.post(TAP_URL, data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": QUERY}, timeout=120)
    response.raise_for_status()
    rows = list(csv.DictReader(StringIO(response.text)))
    if len(rows) != 2000:
        raise ValueError(f"Expected a complete 2,000-row benchmark window, received {len(rows)}")

    rv = np.array([float(row["ccfrvmod"]) * 1000 for row in rows])
    center = float(np.nanmedian(rv))
    records = []
    for row, value in zip(rows, rv):
        uncertainty = float(row["dvrms"]) * 1000
        snr = float(row["extsnr"])
        archive_flag = int(row["flagged"])
        records.append({
            "timestampUtc": row["obsdate"], "bjd": float(row["ccfjdsum"]),
            "rvMetersPerSecond": round(value, 4), "rvRelativeMetersPerSecond": round(value-center, 4),
            "uncertaintyMetersPerSecond": round(uncertainty, 4), "snr": round(snr, 4),
            "sunAltitudeDegrees": round(float(row["sunagl"]), 4), "archiveFlag": archive_flag,
            "analysisPass": archive_flag == 0 and uncertainty < 10 and snr >= 10,
            "filename": row["l2filename"], "archiveMd5": row["l2checksum"], "pipelineVersion": row["swversion"],
        })

    accepted = [record for record in records if record["analysisPass"]]
    accepted_rv = np.array([record["rvRelativeMetersPerSecond"] for record in accepted])
    accepted_snr = np.array([record["snr"] for record in accepted])
    accepted_altitude = np.array([record["sunAltitudeDegrees"] for record in accepted])
    payload = {
        "schema": "rv-activity.neid-solar/1", "generatedAtUtc": datetime.now(timezone.utc).isoformat(),
        "source": {"tapUrl": TAP_URL, "query": QUERY, "querySha256": hashlib.sha256(QUERY.encode()).hexdigest(), "archiveHelp": "https://neid.ipac.caltech.edu/help_solar.php", "credit": "NASA Exoplanet Science Institute / NEID"},
        "selection": {"bjdStart": 2459700, "bjdStop": 2459713, "returned": len(records), "accepted": len(accepted), "archiveFlagged": sum(record["archiveFlag"] == 1 for record in records), "rvCenterMetersPerSecond": center, "rule": "archiveFlag=0 AND dvrms<10 m/s AND extsnr>=10"},
        "diagnostics": {"acceptedRmsMetersPerSecond": float(np.sqrt(np.mean(accepted_rv**2))), "rvSnrPearson": float(np.corrcoef(accepted_rv, accepted_snr)[0,1]), "rvSunAltitudePearson": float(np.corrcoef(accepted_rv, accepted_altitude)[0,1])},
        "records": records,
        "warning": "Archive flags and Level-2 formal RV errors are quality metadata, not a complete activity model. Correlation with S/N or solar altitude is diagnostic and does not establish causation.",
    }
    output = Path("public/data/neid-solar-2459700-2459713.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    print(json.dumps({"output": str(output), "rows": len(records), "accepted": len(accepted), "rms": payload["diagnostics"]["acceptedRmsMetersPerSecond"]}, indent=2))


if __name__ == "__main__": main()
