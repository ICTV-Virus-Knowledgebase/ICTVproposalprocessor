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


@REM primary data tables
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline39].[dbo].[taxonomy_toc]"| findstr /v /c:"-" /b > "taxonomy_toc.utf8.txt"
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline39].[dbo].[taxonomy_node_export]"| findstr /v /c:"-" /b > "taxonomy_node_export.utf8.txt"
@REM replaced with species_isolates
@REM sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline39].[dbo].[vmr]"| findstr /v /c:"-" /b > "vmr.utf8.txt"
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline39].[dbo].[species_isolates]"| findstr /v /c:"-" /b > "species_isolates.utf8.txt"

@REM CV tables
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline39].[dbo].[taxonomy_level]"| findstr /v /c:"-" /b > "taxonomy_level.utf8.txt"
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline39].[dbo].[taxonomy_molecule]"| findstr /v /c:"-" /b > "taxonomy_molecule.utf8.txt"
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline39].[dbo].[taxonomy_host_source]"| findstr /v /c:"-" /b > "taxonomy_host_source.utf8.txt"
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline39].[dbo].[taxonomy_genome_coverage]"| findstr /v /c:"-" /b > "taxonomy_genome_coverage.utf8.txt"
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline39].[dbo].[taxonomy_change_in]"| findstr /v /c:"-" /b > "taxonomy_change_in.utf8.txt"
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline39].[dbo].[taxonomy_change_out]"| findstr /v /c:"-" /b > "taxonomy_change_out.utf8.txt"



@REM cache tables
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline39].[dbo].[taxonomy_node_delta]"| findstr /v /c:"-" /b > "taxonomy_node_delta.utf8.txt"
sqlcmd -s"	" -f o:65001 -W -Q "set nocount on; select * from [ICTVonline39].[dbo].[taxonomy_node_merge_split]"| findstr /v /c:"-" /b > "taxonomy_node_merge_split.utf8.txt"

@REM exports

@REM NCBI_linkout is special
@REM each row is a single column with embedded char(10)+char(13) to make multiple lines
@REM Strip column name and "underline" formatting: findstr /v /c:"#" /b | findstr /v /c:"------------" /x
@REM Strip "Warning: Null value is eliminated by an aggregate": findstr /v /c:"Changed" /b
sqlcmd -s"" -f o:65001 -W -Q "use ICTVonline39; set nocount on; DECLARE @nl varchar(10); SET @nl=char(13)+char(10);exec [NCBI_linkout_ft_export] NULL, @nl" | findstr /v /c:"#" /b | findstr /v /c:"------------" /x | findstr /v /c:"Changed" /b | findstr /v /c:Warning /b > "linkout.utf8.txt"


@REM convenience views

@REM copy back to laptop
copy /Y  *.txt \\tsclient\ICTV\xfer\prod\export_msl
copy /Y  *.bat \\tsclient\ICTV\xfer\prod\export_msl
