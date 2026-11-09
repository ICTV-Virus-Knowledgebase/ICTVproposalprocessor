# proposals_msl41v1_name_case_sw5_075B

Regression test for the case-insensitive taxon-name check in `merge_proposal_zips.R` (commit 8824d74, 2026-11-04). The defect was found while processing MSL42 proposal `2026.075B.N.v1.Herellevirales_1no_6nf_3nsf_14_ng_211ns`. The baseline was generated on 2026-11-09 with `regression_test_og.sh`.

Run it against `current_msl/msl41v1`:

    bash ./regression_test_og.sh name_case_sw5

## Input

`2026.075B.N.v1.Herellevirales_1no_6nf_3nsf_14_ng_211ns.xlsx` is the original, unedited submission from `MSL42/`. Its investigation is in `MSL42_Results/error_proposals_not_edited/README.md`. The proposal makes 265 changes, most of them creates, and reorganizes family *Herelleviridae* into the new order *Herellevirales*. MSL42 has no DOCX for this proposal. The two similarity-matrix workbooks that came with it (7.5 MB, supplementary files) are not included, because only the XLSX is needed to exercise the code.

## What it tests

Row 224 creates species *Kochikohdavirus sw5*, with exemplar OX463805 and virus name Enterococcus phage Sw5. MSL41 already has *Kochikohdavirus Sw5* (ON286976.1) in the same genus, for a virus with the same name. The two names differ only in the capital S.

- **Before the fix:** the taxon-name checks compared names exactly, so `sw5` did not match `Sw5`. Row 224 got `CREATE.OK`, and `msl_load.sql` created a second species next to *Kochikohdavirus Sw5*. The duplicate-accession check did not catch it either, because the accessions are different.
- **After the fix:** row 224 gets ERROR `DEST.IN_CUR_CASE` ("Change=CREATE NEW, but taxon name already exists with different capitalization"), and the change is not applied.

## Expected result

The QC report has one ERROR row, row 224 `DEST.IN_CUR_CASE`, and no WARNING rows. The other 264 changes are applied; the log reports `ERROR:1 INFO:263 SUCCESS:264`. *Kochikohdavirus sw5* is in neither `msl.tsv` nor `msl_load.sql`, and *Kochikohdavirus Sw5* is unchanged.

## Not covered here

This case only exercises a create that collides with a current taxon. It does not exercise:

- a split from *Sw5* to *sw5*, which used to create *sw5* as a new species;
- collisions with a historical name (`DEST.IN_OLD_CASE`) or with a name created earlier in the same run (`DEST.IN_NEW_CASE`);
- collisions at higher ranks;
- a rename that only changes capitalization, which is still allowed.

Known gap, not fixed yet: a "Move; rename" change (`.action` = `move`) skips every destination-name check, exact or case-only, so it can still rename a taxon to an existing name.

This is a full proposal with many moves and promotes. Changes to unrelated QC notes, such as the deferred MOVE lineage-text fix described in `proposals_msl41v1_merge_split_rename_009F/README.md`, will also change this baseline. Only a change at row 224 means this fix has regressed.
