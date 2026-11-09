# MERGE edge case: target_abolished

Synthetic test only; not an approved taxonomic proposal.

Rows merge Ascovirus hvav3a → Ascovirus sfav1a, then abolish sfav1a. Expect MERGE.DEST_NOT_IN_NEW. Both species are absent, and tnav2a survives. Detection does not roll back either operation.

Recovered byte-for-byte from the seven scratch cases documented in the 009F timeline (second round, 2026-11-05). The first three use the edited Marseilleviridae workbook; the other four use Ascovirus.

Run from the repository root:

```sh
bash regression_test_og.sh merge_edge_target_abolished
```

The harness compares the QC report and log with the baselines; it does not compare `msl.tsv` or `msl_load.sql`. Expected ERROR/WARNING diagnostics are part of these negative tests.
