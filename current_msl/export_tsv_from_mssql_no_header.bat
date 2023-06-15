@REM dump a table to a TSV file
@REM -t = field delimiter
@REM -c = character mode (as opposed to "native")
@REM -T = connects to SQL Server with a trusted connection using integrated security
@REM       if -T not specified, you need to specify -U and -P to login
@REM 
@REM Easiest to use a bat file to get the -t"[tab_char]" to work - can't type TAB in to CMD
@REM well, actually, in CMD, you do ALT-009 (which is command-009 on a mac via RDC) to enter a TAB

bcp "select * from ICTVonline38.dbo.taxonomy_toc" queryout output.tsv -t"	" -c -T