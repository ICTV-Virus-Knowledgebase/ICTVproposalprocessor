# MERGE edge case: failed_rename

Synthetic test only; not an approved taxonomic proposal.

Genus rename row 5 is blank. Species rename row 12 fails binomial validation. Expect MERGE.DEST_PENDING_RENAME_FAILED on row 13. This is the mock proposal demonstrating detection without recovery: senegalense is already removed and its merge still points to Marseillevirus massiliense.

Recovered byte-for-byte from the seven scratch cases documented in the 009F timeline (second round, 2026-11-05). The first three use the edited Marseilleviridae workbook; the other four use Ascovirus.

Run from the repository root:

```sh
bash regression_test_og.sh merge_edge_failed_rename
```

The harness compares the QC report and log with the baselines; it does not compare `msl.tsv` or `msl_load.sql`. Expected ERROR/WARNING diagnostics are part of these negative tests.

## Walkthrough: the failed-rename recovery question

Open `2026.902F.N.v1.rename_fails.xlsx`, worksheet **Proposal Template**.
This is the original scratch mock; it is not a proposed change to the live input.

| Step | Workbook row | What happens |
| --- | --- | --- |
| Initial state | — | Marseillevirus contains massiliense and senegalense. Alphamarseillevirus does not exist. |
| Missing prerequisite | 5 (blank) | The genus rename to Alphamarseillevirus was deliberately removed. |
| Merge runs first | 13 | senegalense is merged into the existing massiliense taxon, anticipating row 12. senegalense disappears from the working new MSL. |
| Dependent rename fails | 12 | The proposed species name Alphamarseillevirus massiliense does not match its still-existing genus Marseillevirus. RENAME.SPECIES_BINOMIAL_MISMATCH is reported. |
| Final check | 13 | MERGE.DEST_PENDING_RENAME_FAILED is reported. The merge already applied is still present. |

In the generated `testResults/msl41v1/proposals_msl41v1_merge_edge_failed_rename/msl.tsv`,
Marseillevirus massiliense remains, senegalense is absent, and Alphamarseillevirus
massiliense is absent. In `msl_load.sql`, old taxnode 202503803 still has
`out_change='merge'` pointing to Marseillevirus massiliense. This is an intentionally
invalid partial result, useful for reviewing the QC errors, not a load-ready proposal.

## Recovery options for discussion only

No recovery behavior was added by the 2026-11-06 fix.

1. **Correct the proposal and rerun from the original MSL.** Restore the missing
   genus rename, or otherwise correct the intended taxonomy, then regenerate the
   outputs. This is the simplest operational recovery before any database load.
2. **Defer the merge until its dependency succeeds.** If row 12 fails, row 13
   never removes senegalense. This requires scheduling dependent merges before
   final checks and accounting for other operations that depend on those merges.
3. **Roll back the merge and dependent operations.** Restoring only the deleted
   species is insufficient: the previous-MSL out_change/out_target fields and
   internal merge links must also be restored, and later operations may depend on
   the partial result. A proposal-level transaction or replay would be more reliable
   than manually undoing one row.
4. **Block load-ready export when errors remain.** This does not restore the working
   taxonomy, but can prevent a partial result from being loaded. Diagnostic exports
   could remain available. This would change the processor/export contract.

Recovery after loading incorrect SQL is a separate database/VMR repair problem;
none of these tests executes a database load.
