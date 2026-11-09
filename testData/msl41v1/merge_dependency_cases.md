# MERGE dependency regression cases

All workbooks below are synthetic tests against `current_msl/msl41v1`.
The seven original scratch workbooks from the 2026-11-05 second-round investigation
were recovered byte-for-byte. Two success cases were added on 2026-11-06.

| Directory suffix (`proposals_msl41v1_merge_edge_…`) | Coverage |
| --- | --- |
| [self_source](proposals_msl41v1_merge_edge_self_source/README.md) | Pending rename points to merge source; refuse merge and preserve source for rename |
| [failed_rename](proposals_msl41v1_merge_edge_failed_rename/README.md) | Mock failed dependency, with a step-by-step explanation and recovery options |
| [missing_source](proposals_msl41v1_merge_edge_missing_source/README.md) | Pending rename source cannot be found; refuse merge |
| [chain_forward](proposals_msl41v1_merge_edge_chain_forward/README.md) | A → B → C resolves both merges to C |
| [chain_reverse](proposals_msl41v1_merge_edge_chain_reverse/README.md) | Same chain with spreadsheet rows reversed |
| [cycle](proposals_msl41v1_merge_edge_cycle/README.md) | Refuse merging the survivor back into a taxon already merged away |
| [target_abolished](proposals_msl41v1_merge_edge_target_abolished/README.md) | Report unresolved destination after its abolition |
| [double_rename](proposals_msl41v1_merge_edge_double_rename/README.md) | MERGE A → C, RENAME B → C, RENAME C → D; no false dependency error |
| [move_then_rename](proposals_msl41v1_merge_edge_move_then_rename/README.md) | Same success case with MOVE/RENAME B → C |

Run from the repository root:

```sh
bash regression_test_og.sh merge_edge
```

The harness compares the generated R QC report and log with the matching baselines
under `testResults/msl41v1`. It does not compare `msl.tsv` or `msl_load.sql`, so the
surviving taxa and SQL merge destinations described in the case READMEs are not
checked automatically. Negative cases intentionally contain QC errors; a passing
run means those expected diagnostics are unchanged.
