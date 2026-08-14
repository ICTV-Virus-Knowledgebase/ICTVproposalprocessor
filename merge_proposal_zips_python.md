# Python proposal validator (first conversion)

`merge_proposal_zips.py` is a native Python first-version conversion of
`merge_proposal_zips.R`. It keeps the R CLI option names and output filenames.
Correctness is prioritized over performance in this version.

## Conda

```bash
conda env create -f environment.python.yml
conda activate ictv-proposal-validator-python
./merge_proposal_zips.py \
  --refDir=current_msl/msl41v1 \
  --proposalsDir=testData/msl41v1/proposals_msl41v1_rename \
  --outDir=results \
  --msl
```

An existing bash invocation can replace:

```bash
Rscript merge_proposal_zips.R [flags]
```

with:

```bash
./merge_proposal_zips.py [the same flags]
```

The four cache flags remain accepted so callers do not fail, but emit a warning
and do not load or save a cache.

## Container

```bash
docker build -f Dockerfile.python -t ictv-proposal-validator-python .
docker run --rm \
  -v "$PWD/testData:/testData:ro" \
  -v "$PWD/testResults:/testResults:rw" \
  ictv-proposal-validator-python \
  --refDir=/app/current_msl/msl41v1 \
  --proposalsDir=/testData/msl41v1/proposals_msl41v1_rename \
  --outDir=/testResults/msl41v1/proposals_msl41v1_rename \
  --msl
```

## Intentional first-version boundaries

- Accession reuse is checked only against `species_isolates.utf8.txt`,
  `abolished_accessions.tsv`, and other changes in the current proposal run.
  `taxonomy_node_export.utf8.txt` is not used as an accession source.
- Both old literal-`NULL` reference files and MySQL-style quoted TSV files
  using unquoted `\N`, `""`, and backslash escapes are accepted.
- Accessions in numeric ranges are expanded and checked one by one, with a
  safety ceiling of 10,000 members per range.
- `msl_load.sql` calls `VMR_update_from_MSL()` before `QC_run_modules(NULL)`.
  The repository does not include that stored procedure's definition; if the
  deployed procedure has parameters, adjust that one generated call.
- Proposal and related filenames containing spaces are renamed with underscores
  when collision-free. A QC warning records the rename or a collision.
- Taxonomy helper values are recomputed after cell cleanup, so surrounding
  whitespace cannot change source, parent, suffix, case, or binomial checks.
- Genome Coverage is required only for species rows whose action is
  `Create new`; existing taxa retain it across other actions.
- Same-proposal dependencies are applied before their consumers, including
  created parents and descendant moves needed before an old branch is removed.

## Outputs

With `--msl`, the conversion writes `msl.tsv`, `msl_load.sql`, and
`QC.docx_summary.tsv`, in addition to the QC TSV and XLSX summaries. `-D` also
writes `changes.tsv`.

The consolidated Excel reports are `QC.pretty_summary.all.xlsx` (all QC
messages) and `QC.pretty_summary.issues.xlsx` (ERROR and WARNING).
