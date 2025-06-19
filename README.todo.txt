## IVK-XX1 binomial species name ##

accidentally commited many files with binomial species name fix - should just be merge.R and testResults.
[main 923af18] IVK-XXX fix binomial species name: only allow 2 words, not mulitple
 12 files changed, 16 insertions(+), 2 deletions(-)
 create mode 100644 testData/msl39v4/proposals_msl40_abolish_family_with_kids/README.md
 create mode 100644 testData/msl39v4/proposals_msl40_accession_dash/2024.001P.A.v2.Fimoviridae_1nsp.docx
 create mode 100644 testData/msl39v4/proposals_msl40_accession_dash/2024.001P.A.v2.Fimoviridae_1nsp.xlsx
 create mode 100644 testData/msl39v4/proposals_msl40_accession_dash/README.md
 create mode 100644 testData/msl39v4/proposals_msl40_accession_range/2024.001P.A.v2.Fimoviridae_1nsp.docx
 create mode 100644 testData/msl39v4/proposals_msl40_accession_range/2024.001P.A.v2.Fimoviridae_1nsp.xlsx
 create mode 100644 testData/msl39v4/proposals_msl40_accession_range/2024.005M.A.v2.Cardoreovirus_1nsp.docx
 create mode 100644 testData/msl39v4/proposals_msl40_accession_range/2024.005M.A.v2.Cardoreovirus_1nsp.xlsx
 create mode 100644 testData/msl39v4/proposals_msl40_duplicate_accession/2024.010P.Geminiviridae_Begomovirus_19nsp.docx
 create mode 100644 testData/msl39v4/proposals_msl40_duplicate_accession/2024.010P.Geminiviridae_Begomovirus_19nsp.xlsx
 create mode 100644 testData/msl39v4/proposals_msl40_duplicate_accession/README.md

## IVK-XX2 duplicate accessions ##

First switch to loading species_isolates.txt, instead of VMR.xls
Then parse accessions to make a list
Remember to strip .1's
Parse incoming accessions, and compare

