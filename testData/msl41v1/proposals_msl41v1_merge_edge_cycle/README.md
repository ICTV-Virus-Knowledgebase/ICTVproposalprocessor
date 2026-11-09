# MERGE edge case: cycle

Synthetic test only; not an approved taxonomic proposal.

Rows merge Ascovirus hvav3a → Ascovirus sfav1a, then sfav1a → hvav3a. Expect MERGE.DEST_IS_SRC; sfav1a survives and hvav3a remains merged into it. The ordering pass also reports CHAIN.HIDDEN.

Recovered byte-for-byte from the seven scratch cases documented in the 009F timeline (second round, 2026-11-05). The first three use the edited Marseilleviridae workbook; the other four use Ascovirus.

Run from the repository root:

```sh
bash regression_test_og.sh merge_edge_cycle
```

The harness compares the QC report and log with the baselines; it does not compare `msl.tsv` or `msl_load.sql`. Expected ERROR/WARNING diagnostics are part of these negative tests.
