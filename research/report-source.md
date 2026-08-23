# Internal research synthesis and claim ledger

## Decision

Make a public NEID solar RV quality audit the primary product. Keep the known-truth injection model as a unit-test layer.

## Primary evidence

1. NExScI documents public NEID/EXPRES solar Level-2 data and the RV, uncertainty, S/N, solar altitude, and flag columns: https://neid.ipac.caltech.edu/help_solar.php
2. pyNEID documents unauthenticated public TAP access: https://pyneid.readthedocs.io/en/latest/
3. The official service endpoint is https://neid.ipac.caltech.edu/TAP/sync

## Claim ledger

| Claim | Evidence | Confidence | Design response |
|---|---|---:|---|
| rows are real public NEID Level-2 metadata | official TAP response | high | version exact ADQL and SHA-256 |
| archive flags matter | NExScI column contract | high | preserve and reveal rejected rows |
| formal error alone defines science quality | unsupported | high | add transparent S/N rule and non-claim |
| S/N/altitude correlations identify stellar activity | unsupported | high | label Pearson values diagnostic, not causal |
| metadata-only analysis removes activity | false | high | name missing spectroscopic/solar proxies explicitly |
