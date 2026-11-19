# proposals_msl41v1_accession_trailing_period_010M

Regression test for removing a period that ends an accession, in both `merge_proposal_zips.R` and `merge_proposal_zips.py` (2026-11-19). The defect was found in MSL42 proposal `2026.010M.N.v2.orthobunyavirus_22nsp`. The baselines were generated on 2026-11-19 with `regression_test_og.sh` (R) and `regression_test_python.sh` (Python).

Run it against `current_msl/msl41v1`:

    bash ./regression_test_og.sh accession_trailing_period
    bash ./regression_test_python.sh accession_trailing_period

## Input

`2026.010M.N.v2.orthobunyavirus_22nsp.xlsx` is the original, unedited submission from `MSL42/` (byte-identical to the copy in `2026/`). It creates 22 species in genus *Orthobunyavirus* (rows 5 to 26). The DOCX that came with it isn't included, because only the XLSX is needed to exercise the code.

## What it tests

Row 16 creates species *Orthobunyavirus ipeanense* with exemplar accession `MG029287.`, which ends with a period.

- **Before the fix:** both scripts kept the period, so `msl.tsv` and `msl_load.sql` stored `MG029287.`, and the row's `CREATE.DUP_ACC` note read `accession=MG029287.`.
- **After the fix:** a new value-validation rule on `exemplarAccession` removes a period at the end of the cell or before a `;`. Row 16 gets WARNING `XLSX.ACCESSION.TRAILING_PERIOD_REMOVED` (note `exemplarAccession : MG029287[.]`), and the species is created with accession `MG029287`.

## Expected result

- Row 16 has the WARNING above, then ERROR `CREATE.DUP_ACC` (`accession=MG029287`; *Orthobunyavirus nepuyoi* in MSL41 has the same accession) and `CREATE.OK`. `CREATE.DUP_ACC` is a soft error: the species is still created.
- Rows 6, 8, 10, 11, 13, 21, 22 and 23 also have `CREATE.DUP_ACC`. These accessions are already in MSL41, and they have nothing to do with this fix.
- All 22 species are created; the log reports `ERROR:9 WARNING:1 SUCCESS:22` (Python also `INFO:1`). `MG029287.` (with the period) appears in neither `msl.tsv` nor `msl_load.sql`.
- Python also reports INFO `XLSX.CODE_MISS` for cell E1 (`PYTHON_ONLY` `modern_header_code`); R doesn't.

## Not covered here

The rule also removes a period before a `;` in a list (`S: X1.; M: X2.` becomes `S: X1; M: X2`). This proposal has no such list, so that path isn't exercised.
