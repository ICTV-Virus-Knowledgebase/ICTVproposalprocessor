@REM
@REM export tables to tsv needed for proposal_validator's current_msl/ cache
@REM
@REM 20230614 MSL38
@REM
@REM basic method from 
@REM   https://stackoverflow.com/questions/1355876/export-table-to-file-with-column-headers-column-names-using-the-bcp-utility-an
@REM utf-8 (-f o:65001 ) from 
@REM   https://stackoverflow.com/questions/41561658/i-need-dump-table-from-sql-server-to-csv-in-utf-8
@REM

sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline38].[dbo].[taxonomy_level]"| findstr /v /c:"-" /b > "taxonomy_level.utf8.txt"
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline38].[dbo].[taxonomy_molecule]"| findstr /v /c:"-" /b > "taxonomy_molecule.utf8.txt"
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline38].[dbo].[taxonomy_node_export]"| findstr /v /c:"-" /b > "taxonomy_node_export.utf8.txt"

sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline38].[dbo].[vmr]"| findstr /v /c:"-" /b > "vmr.utf8.txt"
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline38].[dbo].[virus_isolates]"| findstr /v /c:"-" /b > "virus_isolates.utf8.txt"
