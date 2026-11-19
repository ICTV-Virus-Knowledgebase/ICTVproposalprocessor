# proposals_msl41v1_accession_trailing_period_010M

Regression test for removing a period that ends an accession, in both `merge_proposal_zips.R` and `merge_proposal_zips.py` (2026-11-19). The defect was found in MSL42 proposal `2026.010M.N.v2.orthobunyavirus_22nsp`. The baselines were generated on 2026-11-19 with `regression_test_og.sh` (R) and `regression_test_python.sh` (Python).

Run it against `current_msl/msl41v1`:

    bash ./regression_test_og.sh accession_trailing_period
    bash ./regression_test_python.sh accession_trailing_period

## Input

`2026.010M.N.v2.orthobunyavirus_22nsp.xlsx` is the MSL42 submission from `MSL42/` (identical to the original in `2026/`), with one cell edited for this test. It creates 22 species in genus *Orthobunyavirus* (rows 5 to 26). The DOCX that came with it isn't included, because only the XLSX is needed to exercise the code.

**Edited cell (2026-11-19, test copy only; `MSL42/` and `2026/` are unchanged):** AE24 (*Orthobunyavirus Tweentantanense*, exemplar accession) was `MK896515.1`, the L segment only. It is now `S: MK896513.; M: MK896514. ; L: MK896515.1`. The S and M accessions come from the proposal's DOCX, which gives S MK896513.1, M MK896514.1 and L MK896515.1 for isolate HB7-451; their version numbers were replaced with periods to make the test. The edit was made by patching the sheet XML (a new shared string, and the cell's index changed). A cell-by-cell comparison of every sheet with `readxl` showed only AE24 changed. Only `xl/worksheets/sheet2.xml` and `xl/sharedStrings.xml` differ in the zip.

## What it tests

A rule on `exemplarAccession` removes a period that ends an accession, at the end of the cell or before a `;`, and logs WARNING `XLSX.ACCESSION.TRAILING_PERIOD_REMOVED`. Before this fix, both scripts kept the period in `msl.tsv` and `msl_load.sql`.

- **Row 16, end of the cell** (unedited, as submitted): species *Orthobunyavirus ipeanense*, accession `MG029287.`, is stored as `MG029287`. Note: `exemplarAccession : MG029287[.]`.
- **Row 24, before a `;`** (the edited cell): `S: MK896513.; M: MK896514. ; L: MK896515.1` is stored as `S: MK896513; M: MK896514; L: MK896515.1`. A period right before `;` and one followed by a space are both removed, and the versioned `MK896515.1` at the end is kept. Note: `exemplarAccession : S: MK896513[.;] M: MK896514[. ;] L: MK896515.1`. The brackets show the whole match, `;` included; the `;` is put back.

## Expected result

- Rows 16 and 24 each have the WARNING above. Row 16 then has ERROR `CREATE.DUP_ACC` (`accession=MG029287`, the same accession as *Orthobunyavirus nepuyoi* in MSL41) and `CREATE.OK`. `CREATE.DUP_ACC` is a soft error, so the species is still created.
- Rows 6, 8, 10, 11, 13, 21, 22 and 23 also have `CREATE.DUP_ACC`. Their accessions are already in MSL41, which has nothing to do with this fix.
- All 22 species are created, and the log reports `ERROR:9 WARNING:2 SUCCESS:22`. Python also reports `INFO:1`, an `XLSX.CODE_MISS` for cell E1 (`PYTHON_ONLY` `modern_header_code`) that R doesn't have.
- None of `MG029287.`, `MK896513.` or `MK896514.` appears in `msl.tsv` or `msl_load.sql`. R's and Python's `msl_load.sql` are identical. Their `msl.tsv` files differ only in `in_change` (`PYTHON_ONLY` `in_change_text`).
