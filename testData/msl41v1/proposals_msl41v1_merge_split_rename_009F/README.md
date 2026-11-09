# proposals_msl41v1_merge_split_rename_009F

Regression test for two `merge_proposal_zips.R` defects found while processing MSL42 proposal `2026.009F.N.v2.Marseilleviridae_reorg`. The fixes were made on 2026-11-05.

Run it against `current_msl/msl41v1`:

    bash ./regression_test_og.sh merge_split_rename_009F

## Input

`2026.009F.N.v2.Marseilleviridae_reorg.xlsx` is the edited MSL42 working copy from `MSL42/` (edited 2026-10-27). It is not the original submission. Its edits are described in `MSL42_Results/edited_proposals/2026.009F.N.v2.Marseilleviridae_reorg.md`. The DOCX is not included, because only the XLSX is needed to exercise the code.

In MSL41, family *Marseilleviridae* has two genera, each with two species. The proposal:

- renames genus *Marseillevirus* to *Alphamarseillevirus* (row 5), and species *Marseillevirus massiliense* to *Alphamarseillevirus massiliense* (row 12);
- merges *Marseillevirus senegalense* into *Alphamarseillevirus massiliense* (row 13);
- renames genus *Losannavirus* to *Betamarseillevirus* (row 6);
- splits *Losannavirus* into *Gammamarseillevirus* and *Deltamarseillevirus* (rows 10–11), with no split row that keeps the name *Losannavirus*;
- moves *Losannavirus tunisense* to *Gammamarseillevirus tunisense* (row 15), and creates 3 genera and 12 species.

## What it tests

### Merge into a taxon that the same proposal renames

Merges run before renames, so when row 13 runs, *Alphamarseillevirus massiliense* and its genus do not exist yet.

- **Before the fix:** row 13 failed with `MERGE_MOVE_RENAME.PARENT_NOT_FOUND_2`. *Marseillevirus senegalense* stayed in the renamed genus, so row 5 also reported `RENAME_GENUS.SPECIES_BINOMIAL_MISMATCH`.
- **After the fix:** row 13 finds row 12, which renames *Marseillevirus massiliense* to that name, and merges into that species. It logs INFO `MERGE.DEST_PENDING_RENAME`, then `MERGE_INTO.OK`.

The merge records its target (`out_target`) when it is applied, which is before the rename. After all changes, the processor now replaces that value with the target's final name if the target's lineage changed. The database load matches `out_target` to a name or lineage in the new MSL. In `msl_load.sql`, the update for *Marseillevirus senegalense* (taxnode 202503803) should read `out_change='merge'`, `out_target='Alphamarseillevirus massiliense'`.

### Split source that the same proposal renames

The split rows mark *Losannavirus* for removal at the end of the run, because no split row keeps its name. Row 6 renames that same taxon to *Betamarseillevirus*.

- **Before the fix:** the processor tried to remove *Betamarseillevirus*, and reported `SPLIT.IMPLICIT_ABOLISH_WITH_KIDS` (3 rows, one per species) plus INFO `SPLIT.IMPLICIT_ABOLISH`.
- **After the fix:** a split source that the same proposal renames, moves, promotes or demotes counts as kept, so it is not removed. Baseline `msl38/proposals2023.009F` covers the promote case.

## Expected result

The QC report has no ERROR or WARNING rows. The new MSL has 7 genera and 15 species in *Marseilleviridae*, as the proposal DOCX describes.

Known issue, not fixed yet: in `msl.tsv`, the `lineage` text for *Gammamarseillevirus tunisense* still reads `…;Betamarseillevirus;Losannavirus tunisense`. MOVE does not update a moved taxon's own lineage text. Its `parent_id` is correct (*Gammamarseillevirus*), and `msl_load.sql` does not load the lineage column. Fixing this changes the QC notes of every baseline that has a successful move, promote or demote, so it was deferred.
