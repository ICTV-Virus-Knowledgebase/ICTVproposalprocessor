# IVK-406 Extra file crash

If a non-proposal file, that does NOT include '_Suppl.' in the filename is included, the program will crash. 

Before crashing, it will output QC.summary.tsv, which will include the error: 
  * NON_PROPOSAL_NOT_SUPPL
    * Note a proposal XLSX: All supplemental XLSX files should have '_Suppl.' in the name!
    * testData/msl40v1/proposals_msl41_min_extra_file/2024.035B.Uc.v3.Tubulavirales_sim_matrix.xlsx

