# MERGE edge case: chain_reverse

Synthetic test only; not an approved taxonomic proposal.

The same chain as chain_forward, with spreadsheet rows reversed. The processor orders the dependent merges. Expect the same surviving species and final targets, with MERGE.DEST_MERGED.

Recovered byte-for-byte from the seven scratch cases documented in the 009F timeline (second round, 2026-11-05). The first three use the edited Marseilleviridae workbook; the other four use Ascovirus.

Run from the repository root:

```sh
bash regression_test_og.sh merge_edge_chain_reverse
```

The harness compares the QC report and log with the baselines; it does not compare `msl.tsv` or `msl_load.sql`. Expected ERROR/WARNING diagnostics are part of these negative tests.
