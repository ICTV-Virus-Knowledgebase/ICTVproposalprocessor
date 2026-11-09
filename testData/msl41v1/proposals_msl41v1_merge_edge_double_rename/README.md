# MERGE dependency followed by another rename: double_rename

Synthetic test only; not an approved proposal.

| Row | Operation | Current species | Proposed species |
| --- | --- | --- | --- |
| 5 | Merge | Ascovirus hvav3a | Ascovirus interim |
| 6 | Rename | Ascovirus sfav1a | Ascovirus interim |
| 7 | Rename | Ascovirus interim | Ascovirus final |

The row-6 operation succeeds on the original sfav1a taxon, which row 7 renames again. Expect no ERROR or WARNING, one MERGE_INTO.OK, and a final merge target of Ascovirus final. Both hvav3a and sfav1a names disappear; Ascovirus final and Ascovirus tnav2a survive. Before the 2026-11-06 fix, row 5 incorrectly reported MERGE.DEST_PENDING_RENAME_FAILED.

Run `bash regression_test_og.sh merge_edge_double_rename` from the repository root. The harness compares the QC report and log with the baselines; it does not compare `msl.tsv` or `msl_load.sql`.
