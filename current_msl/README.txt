#
# reference data for proposal_validator
#
#
# excel files
#

are downloaded from https://ictv.global
	TP_Template_Excel_module_2023_v1.xlsx
	VMR_MSL38_v1.xlsx

#
# txt (TSV) files are extract from MSSQL on Windows
# then copied down.
#
        # extract
	RDC "AWS ICTV Dev"
	export_current_msl_tables.bat
	# download
	taxonomy_level.utf8.txt
	taxonomy_molecule.utf8.txt
	taxonomy_node_export.utf8.txt
	virus_isolates.utf8.txt
	vmr.utf8.txt
#
# the script can then load these in R and save a cached 
# version as 
#
	.RData
	
