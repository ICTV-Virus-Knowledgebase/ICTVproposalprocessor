# MERGE edge case: chain_forward

Synthetic test only; not an approved taxonomic proposal.

Rows 5–6 merge Ascovirus hvav3a → Ascovirus sfav1a → Ascovirus tnav2a. Both source species disappear; both merge targets resolve to tnav2a; expect MERGE.DEST_MERGED.

Recovered byte-for-byte from the seven scratch cases documented in the 009F timeline (second round, 2026-11-05). The first three use the edited Marseilleviridae workbook; the other four use Ascovirus.

Run from the repository root:

```sh
bash regression_test_og.sh merge_edge_chain_forward
```

The harness compares the QC report and log with the baselines; it does not compare `msl.tsv` or `msl_load.sql`. Expected ERROR/WARNING diagnostics are part of these negative tests.
