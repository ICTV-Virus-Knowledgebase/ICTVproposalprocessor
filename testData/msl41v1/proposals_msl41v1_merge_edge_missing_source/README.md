# MERGE edge case: missing_source

Synthetic test only; not an approved taxonomic proposal.

Row 12 names nonexistent Marseillevirus massiliensis as its source. Expect MERGE.DEST_PENDING_RENAME_FAILED on row 13 before the merge is applied; senegalense remains.

Recovered byte-for-byte from the seven scratch cases documented in the 009F timeline (second round, 2026-11-05). The first three use the edited Marseilleviridae workbook; the other four use Ascovirus.

Run from the repository root:

```sh
bash regression_test_og.sh merge_edge_missing_source
```

The harness compares the QC report and log with the baselines; it does not compare `msl.tsv` or `msl_load.sql`. Expected ERROR/WARNING diagnostics are part of these negative tests.
