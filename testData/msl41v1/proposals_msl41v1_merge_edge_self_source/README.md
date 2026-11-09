# MERGE edge case: self_source

Synthetic test only; not an approved taxonomic proposal.

MERGE source is also the source of the pending RENAME. Expect MERGE.DEST_IS_SRC on row 13; row 12 still creates Alphamarseillevirus massiliense. Marseillevirus senegalense remains and consequently has a genus/binomial error.

Recovered byte-for-byte from the seven scratch cases documented in the 009F timeline (second round, 2026-11-05). The first three use the edited Marseilleviridae workbook; the other four use Ascovirus.

Run from the repository root:

```sh
bash regression_test_og.sh merge_edge_self_source
```

The harness compares the QC report and log with the baselines; it does not compare `msl.tsv` or `msl_load.sql`. Expected ERROR/WARNING diagnostics are part of these negative tests.
