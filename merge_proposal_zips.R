#!/usr/bin/env Rscript
#### SUMMARY ####
#
# Load and QC ICTV taxonomic proposals
# apply proposed changes to current MSL
# Create SQL to load new MSL into MSSQL DB
#
#
# WARNING: we use data.TABLE instead of data.FRAME
#          (much higher performance)
#
#### TO DO ####
#
# IMPLEMENT MERGE 
#
#### CODE OVERVIEW ####
#
# load_proposal_cache() 
# or 
#   scan_for_proposals
#     scan for files (docx, xlsx) & merge -> proposalsDf
#   load & qc proposals
#     for(code) {
#       load docx 
#       load xlsx
#       qc_xlsx
#     }
# merge ->  allChangesDf
# find and link "chained" changes
# order allChangesDf
#
# for(change) {
#   apply_change
#
# write allErrorsDf - > QC.pretty.xlsx, QC.regression.tsv
# write curMSL & newMSL -> new_msl.sql
#    newMSL -> inserts
#    curMSL -> updates (out_target, out_change, out_filename, out_notes)

# This script 
#   * load previous MSL data into a data frame
#   * scans params$ ./proposals for *.(xlsx|docx) (proposals)
#   * iterates over each proposal .xlsx
#     * reads the proposal.xlsx file into a data table 
#     * merge each row (change)
#       * match up to existing taxon
#       * QC data, etc.
#       * implement change in new MSL data frame
#       * record errors, warnings, etc
#       * create merged change set in one merged dataframe
#   * write merged proposal data frame to a Unicode(UTF-16LE) TSV file (params$merged) that can be loaded into MSSQL on Windows using "Import Data...." 
#   * write a status sheet listing parsing and QC success/fail status for each proposal (params$status)
#   * write new MSL load & updates to prev_msl.out_*
#
#### LIBRARIES ####
# this allows modification in place of a (data) passed
# to a subroutine (pass-by-reference feature)
suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(yaml))
#suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))  # also part of tidyverse
suppressPackageStartupMessages(library(dplyr))  # also part of tidyverse
suppressPackageStartupMessages(library(readr))  # also part of tidyverse

suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(writexl) )# another option library(openxlsx)
# this package is difficult to install for Dockerfile, use other methods
#suppressPackageStartupMessages(library(DescTools)) # AscToChar
suppressPackageStartupMessages(library(qdapTools)) # read docx
suppressPackageStartupMessages(library(tools)) # file_ext()


#### GLOBAL VARS ####
#   vector      .GlobalEnv$params 
#   data.table  .GlobalEnv$allErrorDf
#   data.frame  .GlobalEnv$proposalsDf  
#     code      (row.name)
#     xlsx      (xlsx basename + .xlsx)
#     xlsxpath  (xlsx path)
#     docx      (docx basename +.docx)
#     docxpath  (docx path)
#     basename  (basename for proposal/.zip: from xlsx|docx)
#     scAbbrev  (single letter study section abbrev)
#     cleanbase (strip off version, workflow status and .fix, to get final, production filename)
#
# BUILT FROM REF DIR
#   list        .GlobalEnv$cvList            of vectors
#   list        .GlobalEnv$dbCvMapList       of vectors
#   data.table  .GlobalEnv$oldMSLs
#   data.table  .GlobalEnv$curMSL
#   data.table  .GlobalEnv$newMSL
#   vector      .GlobalEnv$scAbbrevNameMap
#   data.frame  .GlobalEnv$vmrDf
#
# WORKING VARS
#   int         .GlobalEnv$actionOrder
#
##### INIT #####
.GlobalEnv$actionOrder = 0
.GlobalEnv$cvList      = list()
##### create allErrorDf   ##### 
# 
# error reporting data frame
# 
.GlobalEnv$allErrorDf = data.table(
  "subcommittee" = character(),
  "code" = character(),
  "docx" = character(),
  "xlsx" = character(),
  "row"  = integer(),
  "change" = character(),
  "rank" = character(),
  "taxon" = character(),
  "level" = factor(levels=c("ERROR","WARNING","INFO")),
  "error" = character(),
  "message" = character(),
  "notes" = character(),
  "scAbbrev" = character(),
  "validator_version" = character(),
  "order"  = integer()
)


#### Parse Args ####
suppressPackageStartupMessages(library("optparse"))
option_list <- list( 
  # verbose/quiet/debug
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output"),
  make_option(c("-t", "--tmi"), action="store_true", default=FALSE,
              help="Print lots of extra output (Too Much Information)"),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose", 
              help="Print no output"),
  make_option(c("-d", "--debug"), action="store_true", default=FALSE, dest="debug_on_error", 
              help="Call browser() on data error [default \"%default\"]"),
  make_option(c("--noInfo"), action="store_false", default=TRUE, dest="show.xlsx.code_miss", 
              help="Supress INFO level warnings from output [default \"%default\"]"),
  make_option(c("-c", "--useCache"), action="store_true", default=FALSE, dest="use_cache", 
              help="Load .RData cache in refDir [default \"%default\"]"),
  make_option(c("-u", "--updateCache"), action="store_true", default=FALSE, dest="update_cache", 
              help="Write .RData cache in refDir after processings all proposals in proposalDir [default \"%default\"]"),
  make_option(c("-C", "--loadProposalCache"), action="store_true", default=FALSE, dest="load_proposal_cache", 
              help="Load .RData cache in proposalDir [default \"%default\"]"),
  make_option(c("-U", "--saveProposalCache"), action="store_true", default=FALSE, dest="save_proposal_cache", 
              help="Save .RData cache in proposalDir [default \"%default\"]"),
  make_option(c("--msl"), action="store_true", default=FALSE, dest="export_msl", 
              help="Write resulting MSL to load_msl.sql and msl.tsv [default \"%default\"]"),
  make_option(c("-D","--outputDeltas"), action="store_true", default=FALSE, dest="output_change_report", 
              help="Write taxa change list to a 2nd sheet in QC report [default \"%default\"]"),
  
  
  # in/out/ref directories
  make_option(c("-i","--proposalsDir"), 
              default="proposalsTest", 
              #default="testData", 
              #default="proposalsTestEmpty", 
              #default="proposalTest1Doc", 
              #default="proposalTestNonProXlsx", 
              dest="proposals_dir",
              help = "Directory to scan for YYYY.###SC.*.xlsx proposal files [default \"%default\"]"),
  make_option(c("-o","--outDir"), 
	      default="results", 
	      #default="testResults", 
	      dest="out_dir",
              help = "Directory to write outputs to [default \"%default\"]"),
  make_option(c("-r","--refDir"), default="current_msl/msl39v4", dest="ref_dir", 
              help="Directory from which read current MSL and CV data from [default \"%default\"]"),
  
  # out filenames
  make_option(c("--mslTsv"), default="msl.tsv", dest="msl_tsv",
              help="Output file, in outDir, for resulting MSL w/o IDs (for diff) [default \"%default\"]"),
  make_option(c("--docxSummary"), default="QC.docx_summary.tsv", dest="proposals_meta",
              help="Output file, in outDir, for meta-data from proposal docx's [default \"%default\"]"),
  
  
  make_option(c("--mslLoadSql"), default="msl_load.sql", dest="sql_load_filename",
              help="Output file, in outDir, for SQL to load new MSL into db [default \"%default\"]"),
  make_option(c("--sqlInsertBatch"), default="200", dest="sql_insert_batch_size",
              help="Number of rows of data per SQL INSERT line in outDir/mslLoadSql file [default \"%default\"]"),
  make_option(c("--taxnodeIdDelta"), default="100000", dest="taxnode_delta", 
              help="Amount to increment taxnode_ids by to create a new MSL [default \"%default\"]"),
  make_option(c("--qcPrettySummary"), default="QC.pretty_summary.all.xlsx", dest="qc_summary_fname", 
              help="Primary error/warning output xlsx file, returned to user [default \"%default\"]"),
  make_option(c("--qcTsvSummary"), default="QC.summary.tsv", dest="qc_summary_tsv_fname", 
              help="Primary error/warning output TSV file, parsed by web app [default \"%default\"]"),
  make_option(c("--qcTsvRegression"), default="QC.regression.new.tsv", dest="qc_regression_tsv_fname", 
              help="Primary error/warning output TSV file, used by regression testing (no version column) [default \"%default\"]"),
  
  # ref filenames
  make_option(c("--dbSourceHost"), default="taxonomy_host_source.utf8.txt", dest="db_host_source_fname",
              help="Reference file listing allowable source/host values [default \"%default\"]"),
  make_option(c("--dbRanks"), default="taxonomy_level.utf8.txt", dest="db_rank_fname",
              help="Reference file listing allowable ranks [default \"%default\"]"),
  make_option(c("--dbMolecules"), default="taxonomy_molecule.utf8.txt", dest="db_molecule_fname",
              help="Reference filelisting allowable genomic molecules [default \"%default\"]"),
  make_option(c("--dbTaxa"), default="taxonomy_node_export.utf8.txt", dest="db_taxonomy_node_fname",
              help="Reference file listing all historic taxa [default \"%default\"]"),
  make_option(c("--cvTemplate"), default="TP_Template.xlsx", dest="template_xlsx_fname",
              help="Template proposal xlsx, used to load CVs [default \"%default\"]"),
  make_option(c("--cvTemplateSheet"), default="Menu Items (Do not change)", dest="template_xlsx_sheet",
              help="Template proposal xlsx, used to load CVs [default \"%default\"]"),
  make_option(c("--vmr"), default="VMR.xlsx", dest="vmr_fname",
              help="VMR to check for accession re-use [default \"%default\"]"),
  make_option(c("--templateURL"), default="https://ictv.global/taxonomy/templates", dest="template_url",
              help="URL for out-of-date template message [default \"%default\"]"),
  make_option(c("--cacheFile"), default=".RData", dest="cache_fname",
              help="Filename for refDir/cache [default \"%default\"]"),
  
  # input filenames
  make_option(c("--infileSupplPat"), default="_Suppl.", dest="infile_suppl_pat",
              help="Input files (in proposalDir/) containing this pattern in filename are ignored [default \"%default\"]"),
  make_option(c("--mode"), default="validate", dest="processing_mode",
              help="Stringency Mode: 'validate' (any .doc|docx|xls|xlsx), 'draft' (YYYY.###A.v#.'), or 'final' (no '.v#') [default \"%default\"]"),
  make_option(c("--version_file"), default="version_git.txt", dest="version_file",
              help="File to use to report version. v#.#.###### (last part is git hash); run ./version_git.sh to update [default \"%default\"]"),
  make_option(c("--version"), dest="version",
              help="code version  [default \"%default\"]"),
  
  # metadata for new MSL
  make_option(c("--newMslName"), default="YYYY", dest="msl_name",
              help="Root node name for new MSL [default \"%default\"]"),
  make_option(c("--newMslNotes"), default="Provisional EC ##, Online meeting, July YYYY; Email ratification March YYYY (MSL ###)", dest="msl_notes",
              help="Description for the new MSL for database loading [default \"%default\"]")
  
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
params <- parse_args(OptionParser(option_list=option_list))
#debuggingState(FALSE)
#cat("# DebuggingState()=",debuggingState(),"\n")

##### debug overrides #####
if( interactive() ) {
  print("!!!!||||||||||||||||||||||||!!!!")
  print("!!!! DEBUG OVERIDES ENGAGED !!!!")
  print("!!!!||||||||||||||||||||||||!!!!")
  options(error=recover)
  #options(error=browser)
  params$debug=F
  # WARNING - this will store all the other debug settings into the proposalDir/.RData file!
  params$load_proposal_cache = F
  params$save_proposal_cache = F
  # defeat auto-caching when debugging
  #rm(docxList,xlsxList,changeList)
  params$verbose = F
  params$tmi = F
  params$debug_on_error = F
  params$processing_mode = 'final'
  #params$output_change_report = F
  params$export_msl = F
#  params$test_case_dir = "crash"
#  params$test_case_dir = "proposalsEC55.1"
  # MSL39v2
  params$test_case_dir = 'proposals_msl39v2'
  params$proposals_dir = paste0("testData/msl39/",params$test_case_dir)
  params$out_dir       = paste0("testResults/msl39/",params$test_case_dir)

  # fast debugging merge/chain
  params$processing_mode ="final"
#  params$proposals_dir = "./MSL39v6/proposalsFinal"
#  params$out_dir       = "./MSL39v6/results/proposalsFinal"
  #params$proposals_dir = "EC55"
  #params$out_dir       = "EC55_results"
  params$qc_regression_tsv_fname = "QC.regression.new.tsv"
  cat("!! VERBOSE: ",params$versbose, "\n")
  cat("!! TMI:     ",params$tmi, "\n")
  cat("!! MODE:    ",params$processing_mode, "\n")
  cat("!! SRC_DIR: ",params$proposals_dir, "\n")
  cat("!! OUT_DIR: ",params$out_dir, "\n")
}
#
# args trickery - TMI implies verbose
#
if(params$tmi) { params$verbose=T}

#### IN/OUT DIRS ####

if(params$verbose){ cat(paste0("PROPOSALS_DIR:  ",params$proposals_dir,"/\n")) }
if(params$verbose){ cat(paste0("OUT_DIR:        ",params$out_dir,"/\n")) }
dir.create(params$out_dir,recursive=T,showWarnings = F)

if(params$verbose){ cat(paste0("REF_DIR:        ",params$ref_dir,"/\n"))}


#### UTILITY FUNCTIONS ####

# log an error or status to global list
#
# logs to global allErrorDf
log_change_error = function(curChangeDf,levelStr, errorCode, errorStr, notes="") {
  if( is.null(curChangeDf) && interactive() ) {browser()}
  return(log_error(code=curChangeDf$.code, linenum=curChangeDf$.linenum, 
                   action=curChangeDf$change,   # original text from proposal.xlsx
                   #action=curChangeDf$.action, # processed internal verb
                   actionOrder=curChangeDf$.order,
                   rank=curChangeDf$rank, # ? ????
                   taxon=curChangeDf$.changeTaxon, # ?????
                   levelStr=levelStr, errorCode=errorCode, errorStr=errorStr, notes=notes))
}
log_error=function(code,linenum,action,rank,taxon,levelStr,errorCode,errorStr,notes="",actionOrder="") {
  #cat("log_error(",code,",linenum=",linenum,"...,action=",action,",rank=",rank,",taxon=",taxon,"...)\n")
  #
  # handle non-numeric linenum, gracefully
  # sometimes linenum is empty, or is "code:linenum" 
  # sometimes it's just a numeric, like it should be.
  #
  real_linenum = NA
  if( class(linenum) == "character" ) {
    linenum_parts = unlist(strsplit(linenum,":"))
    if( length(linenum_parts)>1 ) { 
      real_linenum = as.numeric(linenum_parts[2])
    }
  }
  nextErrorDf = data.frame(
    subcommittee = proposalsDf[code,]$subcommittee,
    code = code,
    row = as.numeric(linenum),
    # overide global actionOrder, if needed
    order = ifelse(!is.na(actionOrder) && actionOrder != "",actionOrder,.GlobalEnv$actionOrder),
    change = action, 
    rank=rank,
    taxon = taxon,
    docx = ifelse(is.na(proposalsDf[code,]$docx),"MISSING",proposalsDf[code,]$docx),
    xlsx = ifelse(is.na(proposalsDf[code,]$xlsx),"MISSING",proposalsDf[code,]$xlsx),
    level = levelStr,
    error = errorCode,
    message = errorStr, 
    notes = notes)
  # if ERROR, debug, if desired
  if( params$debug_on_error && interactive() && levelStr=="ERROR") { browser()}
  # appended to GLOBAL error list
  # (fill=TRUE works because first arg of rbind is a data.table)
  .GlobalEnv$allErrorDf=rbind(.GlobalEnv$allErrorDf,nextErrorDf,fill=TRUE)
} # log_error()

#
# write error output file
# 
# we do this frequently, so there will be an artifact to return to the user
# in case we crash part way through
#
write_error_summary = function(errorDf,final=FALSE) {
 
  # make sure version is correct
  errorDf$validator_version = params[["version"]]

  # sort by row of worksheet
  errorSortCols = errorDf[,c("code","row","order")]
  errorSortCols$row_n = as.integer(errorSortCols$row)
  errorSortCols$order_n = as.integer(errorSortCols$order)
  errorsSorted = do.call(order,errorSortCols[,c("code","row_n","order_n")])
  
  # ..........PRETTY FORMAT ..............
  # write current error list
  # TODO: add line breaks, bold
  # openxlsx
  # r2excel: (install fails?) http://www.sthda.com/english/wiki/r2excel-read-write-and-format-easily-excel-files-using-r-software#install-and-load-r2excel-package
  # xlsx (requires java)
  # XLConnect (requires rJava)
  
  prettyErrorDf = data.frame(errorDf %>% filter(FALSE))
  if( nrow(errorDf) > 0) {
    prettyRow = 0
    prevCode = ""
    
    for(i in seq(1,nrow(errorDf)) ) {
      row=errorsSorted[i]
      # add blank line and header when document changes
      rowCode = errorDf[row,"code"]
      if(!is.na(rowCode) && rowCode != prevCode) { 
        prevCode = errorDf[row,"code"]
        
        prettyErrorDf[prettyRow,c("subcommittee")] = c(errorDf[row,"subcommittee"])
        prettyRow=prettyRow+1
        
        prettyErrorDf[prettyRow,c("subcommittee","code","xlsx")] = c(
          errorDf[row,"subcommittee"],
          errorDf[row,"code"],
          ifelse(is.na(allErrorDf[row,"xlsx"]),
                 errorDf[row,"docx"],
                 errorDf[row,"xlsx"]
          )
        )
        prettyRow=prettyRow+1
      }
      
      # copy other lines as-is
      prettyErrorDf[prettyRow,]=errorDf[row,]
      prettyRow=prettyRow+1
    }
  }
  prettyCols = grep(names(prettyErrorDf),pattern="(code|docx|subcommittee)",invert=T,value=T)
  
  # if( params$output_change_list ) {
  #   #
  #   # create change list report
  #   #
  #   print(names(.GlobalEnv$curMSL %>% filter(!is.na(.actionOrder))))
  #   print(names(.GlobalEnv$newMSL %>% filter(!is.na(.actionOrder))))
  #   print(intersect( names(.GlobalEnv$newMSL %>% filter(!is.na(.actionOrder))),
  #                    names(.GlobalEnv$curMSL %>% filter(!is.na(.actionOrder)))
  #   ))
  #   finalChangesColumns = c(".actionOrder","taxnode_id.old","lineage.old","rank.old","name.old","change","name.new","rank.new","lineage.new","taxnode_id.new")
  #   finalChangesDf=merge(.GlobalEnv$curMSL %>% filter(!is.na(.actionOrder)),
  #                        .GlobalEnv$newMSL %>% filter(!is.na(.actionOrder)),
  #                        by=c(".actionOrder"),
  #                        suffixes=c(".old",".new"),
  #                        all=TRUE)
  #   # only one change should be !NA.
  #   finalChangesDf$change = paste0(
  #     na.omit(finalChangesDf$out_change.old), 
  #     na.omit(finalChangesDf$in_change.new)
  #   )
  #   finalChangesOrder = order(finalChangesDf$.actionOrder)
  # }
  
  #
  # if final report, also save by SC
  #
  if( final ) {
    # this should get chunked into files/worksheets by "subcommittee"
    for(committee in levels(as.factor(prettyErrorDf$subcommittee)) ) {
      subcommitteeFilename = str_replace_all(str_replace(committee,pattern="(.*) \\(([A-Z])\\) .*","\\2 \\1")," ","_")
      filename = file.path(params$out_dir,
                           paste0("QC.pretty_summary.",subcommitteeFilename,".xlsx")
      )
      
      select = (prettyErrorDf$subcommittee == committee)
      write_xlsx( x=prettyErrorDf[select,prettyCols],path=filename)
      if(params$verbose) {cat("Wrote: ", filename, " (",nrow(prettyErrorDf[select,]),"rows )\n")}
    }
  }
  
  #..... WRITE FILES ....
  # XLS version (for end human user)
  fname = file.path(params$out_dir,params$qc_summary_fname)
  xlsxSheets = list("QC_Report"=as.data.frame(prettyErrorDf)[,prettyCols])
#  if( params$output_change_report ) {
#    xlsxSheets[["Change_List"]] = as.data.frame(finalChangesDf)[finalChangesOrder,finalChangesColumns]
#    }
  write_xlsx( x=xlsxSheets,path=fname)
  if(params$verbose || params$tmi) {
    cat(paste0("Wrote: ", fname, " with ",length(xlsxSheets)," sheets\n"))
    cat(paste0("     sheet '", names(xlsxSheets)[1]," ",nrow(xlsxSheets[[1]])," rows)\n"))
#    if( params$output_change_report ) {    
#      cat(paste0("     sheet '", names(xlsxSheets)[2]," ",nrow(xlsxSheets[[2]])," rows)\n"))
#    }
  }
  # TSV version for web app parsing
  fnameTsv = file.path(params$out_dir,params$qc_summary_tsv_fname)
  write_delim(x=errorDf[errorsSorted],file=fnameTsv, delim="\t")
  if(params$verbose || params$tmi) {cat(paste0("Wrote: ", fnameTsv, " (",nrow(errorDf)," rows)\n"))}
  else {cat(paste0("Wrote: OUTDIR/", params$qc_summary_tsv_fname, " (",nrow(errorDf)," rows)\n"))}
  # TSV version for regression testing (no version column)
  fnameTsv = file.path(params$out_dir,params$qc_regression_tsv_fname)
  nonVersionCols=grep(names(errorDf),pattern="version",invert=T)
  write_delim(x=errorDf[errorsSorted,..nonVersionCols],file=fnameTsv, delim="\t")
  if(params$verbose || params$tmi) {cat(paste0("Wrote: ", fnameTsv, " (",nrow(errorDf)," rows)\n"))}
}

#
# load version file
#
load_version = function() {
  # 
  # does it exist? 
  # 
  if( !file.exists(params$version_file) ) {
    # try to create it
    cat(paste0("BUILD_VERSION_GIT: ",params$version_file," does not exists. Rebuilding with ./version_git.sh","\n"))
    system("./version_git.sh")
  }
  if( !file.exists(params$version_file) ) {
    cat(paste0("BUILD_VERSION_GIT: ",params$version_file," still does not exist, after running ./version_git.sh","\n"))
    stop(1)
  }
  #
  # load version number, store in params
  #
  .GlobalEnv$params[["version"]] = read.table("version_git.txt")[1,1]
  #
  # print
  # 
  cat( paste0("VERSION: ",params$version,"\n"))
}

#### load version ####

load_version()

#### SECTION: LOAD REF  ####

#........... function to copy prev MSL to new MSL......................................
#
# this uses the DATABASE schema for the taxonomy_node table (ie, naming convention)
#
# extra admin fields added, which wont be saved
#   .prev_taxnode_id
#   prev_proposals
#
create_new_msl = function(curMSL,prev_msl,dest_msl,taxnode_delta) {
  # curMSL = curMSL; prev_msl=max(curMSL$msl_release_num); dest_msl=max(curMSL$msl_release_num)+1; taxnode_delta=as.integer(params$taxnode_delta)
  
  #
  # copy previous years MSL rows to create a basis for this years
  #
  # copies the rows having msl_release_num=prev_msl and offsets all their taxnode_id (and FK(taxnode_id) columns by (taxnode_delta)
  #
  # copy the rows
  newMSL=copy(subset(curMSL, msl_release_num==prev_msl))
  
  #
  # add some admin rows
  #
  
  # copy original taxnode_id's into a new column
  newMSL[,".prev_taxnode_id"] = newMSL[,"taxnode_id"]
  
  # CSV of proposals that have already modified this taxon in this MSL
  newMSL[,"prev_proposals"] = NA_character_
  
  # add missing columns
  for( name in setdiff(c("exemplar_name","genome_coverage"),names(newMSL))) {
    newMSL[,name] = NA_character_
  }
  #
  # update all ids (taxnode_id and all FK's to taxnnode_id) for new MSL
  #
  
  # update the node ids
  fkTaxnodeIdCols = c("taxnode_id","parent_id","tree_id"
                      # ,paste0(tolower(
                      #   c("Realm", "Subrealm", "Kingdom", "Subkingdom", 
                      #     "Phylum", "Subphylum", "Class", "Subclass", "Order", "Suborder", 
                      #     "Family", "Subfamily", "Genus", "Subgenus", "Species")),"_id")
  )
  newMSL[,(fkTaxnodeIdCols) := lapply(.SD, function(id) id+as.integer(taxnode_delta)),.SDcols=fkTaxnodeIdCols]
  
  # update the MSL number
  newMSL[,msl_release_num:=dest_msl]
  
  # make row ids match
  #rownames(newMSL) = newMSL$taxnode_id
  
  # clear in/out changes
  clearCols = c(grep(names(newMSL),pattern="^(in|out)_",value=T))
  newMSL[,(clearCols) := NA]
  
  # add admin columns
  newMSL[,".emptyReported"] = NA_character_
  newMSL[,".otherLineage"] = NA_character_
  newMSL[,".otherLineageProposal"] = NA_character_
  newMSL[,".otherLineageAction"] = NA_character_
  
  # rename release
  newMSL[newMSL$rank=="tree","name"]  = params$msl_name
  newMSL[newMSL$rank=="tree","notes"] = params$msl_notes
  
  # add the new rows
  #return(rbind(taxonomy_node, newMSL))
  return(newMSL)
} # create_new_msl()


# ............... load_reference() ..................
#
# creates data structures
#   list        .GlobalEnv$cvList            of vectors
#   data.table  .GlobalEnv$oldMSLs
#   data.table  .GlobalEnv$curMSL
#   data.table  .GlobalEnv$newMSL
#   vector      .GlobalEnv$scAbbrevNameMap
#   list        .GlobalEnv$dbCvMapList       of vectors
#   data.frame  .GlobalEnv$vmrDf
refGlobalObjs = c("cvList","oldMSLs","curMSL","newMSL","scAbbrevNameMap","dbCvMapList","vmrDf")
#
# internal 
#   vector      taxonomy_node_names
#   data.table  taxonomyDt 
#
# MSL37 loads in 0.22 sec (nThreads=1), slower with more threads

load_reference_cache=function() {
  cacheFilename = file.path(params$ref_dir,params$cache_fname)
  if(params$verbose){ cat("LOAD: cache file", cacheFilename, "...\n")}
  load(file=cacheFilename,verbose=T,envir=.GlobalEnv)
  
}
save_reference_cache=function() {
  cacheFilename = file.path(params$ref_dir,params$cache_fname)
  if(params$verbose){ cat("WRITE: cache file ", cacheFilename, "\n")}
  save(file=cacheFilename,list=refGlobalObjs)

}
load_reference=function() {
  #
  # this uses the DATABASE schema for the taxonomy_node table (ie, naming convention)
  #
  
  #### load taxonomy_node ####
  
  #taxonomy_node_seq = "select  * from taxonomy_node_export"
  taxonomy_node_names = c(
    taxnode_id="integer",parent_id="integer",tree_id="integer",msl_release_num="integer",level_id="integer",
    name="character",
    ictv_id="integer",
    molecule_id="integer",
    abbrev_csv="character",
    genbank_accession_csv="character",
    genbank_refseq_accession_csv="character",
    refseq_accession_csv="character",
    isolate_csv="character",
    notes="character",
    is_ref="factor",is_official="factor",is_hidden="factor",is_deleted="factor",
    is_deleted_next_year="factor",is_typo="factor",is_renamed_next_year="factor",
    is_obsolete="factor",
    in_change="factor",in_target="character",in_filename="character",in_notes="character",
    out_change="factor",out_target="character",out_filename="character",out_notes="character",
    # start_num_sort="integer",
    # row_num="integer",filename="character",xref="character",
    # realm_id="integer",realm_kid_ct="integer",realm_desc_ct="integer",
    # subrealm_id="integer",subrealm_kid_ct="integer",subrealm_desc_ct="integer",
    # kingdom_id="integer",kingdom_kid_ct="integer",kingdom_desc_ct="integer",
    # subkingdom_id="integer",subkingdom_kid_ct="integer",subkingdom_desc_ct="integer",
    # phylum_id="integer",phylum_kid_ct="integer",phylum_desc_ct="integer",
    # subphylum_id="integer",subphylum_kid_ct="integer",subphylum_desc_ct="integer",
    # class_id="integer",class_kid_ct="integer",class_desc_ct="integer",
    # subclass_id="integer",subclass_kid_ct="integer",subclass_desc_ct="integer",
    # order_id="integer",order_kid_ct="integer",order_desc_ct="integer",
    # suborder_id="integer",suborder_kid_ct="integer",suborder_desc_ct="integer",
    # family_id="integer",family_kid_ct="integer",family_desc_ct="integer",
    # subfamily_id="integer",subfamily_kid_ct="integer",subfamily_desc_ct="integer",
    # genus_id="integer",genus_kid_ct="integer",genus_desc_ct="integer",
    # subgenus_id="integer",subgenus_kid_ct="integer",subgenus_desc_ct="integer",
    # species_id="integer",species_kid_ct="integer",species_desc_ct="integer",
    # taxa_kid_cts="character",taxa_desc_cts="character",
    # inher_molecule_id="integer",left_idx="integer",right_idx="integer",node_depth="integer",
    lineage="character",
    cleaned_name="character",
    # cleaned_problem="character",flags="character",
    # "_numKids"="integer","_out_target_parent"="factor","_out_target_name"="factor",
    rank="factor",
    #tree="factor",
    #  realm="factor",subrealm="factor",kingdom="factor",subkingdom="factor",phylum="factor",subphylum="factor",class="factor",subclass="factor",order="factor",suborder="factor",family="factor",subfamily="factor",genus="factor",subgenus="factor",species="factor",
    molecule="factor"
    #inher_molecule="factor"
  )
  #taxonomyDf=read.delim(file=params$prev_taxa_fname,header=FALSE,col.names=names(taxonomy_node_names),stringsAsFactors=FALSE,na.strings="NULL")
  dbTaxonomyNodeFilename=file.path(params$ref_dir, params$db_taxonomy_node_fname)
  taxonomyDt=withCallingHandlers(
    warning = function(cnd) {
      stop(paste0("ERROR: bad export in ",dbTaxonomyNodeFilename,": ",cnd$message ))
    },
    fread(file=dbTaxonomyNodeFilename,
                     header=TRUE,#header=FALSE,col.names=names(taxonomy_node_names),
                     colClasses=as.character(taxonomy_node_names),
                     stringsAsFactors=FALSE,na.strings=c("","NULL"),
                     key=c("taxnode_id"), index=c("name","msl_release_num","parent_id"),
                     nThread = 1
    )
  )
  cat("Previous taxa:",dim(taxonomyDt), " from ",dbTaxonomyNodeFilename,"\n")
  if( !("host_source" %in% names(taxonomyDt)) ) {
    cat("WARNING: no host_source column in taxonomy_node dump!!! (Adding)\n")
    taxonomyDt[,"host_source"] = NA_character_
  }
  
  #### subset taxonomy: last, old, & cur  #### 
  
  lastMSL = max(as.integer(taxonomyDt$msl_release_num))
  .GlobalEnv$oldMSLs = subset(taxonomyDt, msl_release_num < lastMSL)
  .GlobalEnv$curMSL = subset(taxonomyDt, msl_release_num==lastMSL)
  # add accounting columns
  .GlobalEnv$curMSL[,"out_updated"] = FALSE
  
  #
  # copy prev MSL to make new MSL,
  # to which we will try and apply these edits
  #
  .GlobalEnv$newMSL=create_new_msl(.GlobalEnv$curMSL,lastMSL, lastMSL+1, params$taxnode_delta)
  
  
  #### load CVs  #### 
  
  # DB CV loads
  .GlobalEnv$dbCvList = list()
  .GlobalEnv$dbCvMapList = list()
  
  ##### CVs from db dumps #####
  dbRankFilename=file.path(params$ref_dir, params$db_rank_fname)
  rankCV = read.delim(file=dbRankFilename,header=TRUE, stringsAsFactors=TRUE,na.strings=c("NULL"))
  rownames(rankCV) = rankCV$id
  .GlobalEnv$dbCvList[["rank"]] = rankCV
  .GlobalEnv$dbCvMapList[["rank"]] = rankCV$id
  names(.GlobalEnv$dbCvMapList[["rank"]]) = rankCV$name
  if(params$verbose) {cat("RankCV: ", dim(rankCV), " from ",dbRankFilename,"\n")}
  
  
  dbMoleculeFilename=file.path(params$ref_dir, params$db_molecule_fname)
  moleculeCV = read.delim(file=dbMoleculeFilename,header=TRUE, stringsAsFactors=TRUE,na.strings=c("NULL"))
  rownames(moleculeCV) = moleculeCV$id
  .GlobalEnv$dbCvList[["molecule"]] = moleculeCV
  .GlobalEnv$dbCvMapList[["molecule"]] = moleculeCV$id
  names(.GlobalEnv$dbCvMapList[["molecule"]]) = moleculeCV$abbrev
  if( params$verbose) {cat("MoleculeCV: ", dim(moleculeCV), " from ",dbMoleculeFilename,"\n")}
  
  ##### CVs from Proposals Template #####
  
  #
  # this uses the PROPOSAL.XLSX schema (naming convention)
  #
  # we should probably load CVs from a series of DB dumps, rather than using
  # the excel file as the reference, in order to make them easier to update. 
  #
  
  refProposalTemplateFilename=file.path(params$ref_dir, params$template_xlsx_fname)
  templateProposalCV = suppressMessages(data.frame(read_excel(refProposalTemplateFilename,sheet = params$template_xlsx_sheet,col_names = FALSE)))
  
  #cvDf = data.f rame(trib[,])  # remove "select one" line
  if( params$verbose) {cat("ProposalTemplate[",params$template_xlsx_sheet,"]: ", dim(templateProposalCV), " from ",refProposalTemplateFilename,"\n")}
  # patch up studysection, which is a map, not just a CV list
  templateProposalCV[1,6]="Subcommittee Abbrev"
  templateProposalCV[1,7]="Subcommittee Name"
  # fix any \r, \n, \r, etc in column headers, and convert punctuation to spaces, too
  templateProposalCV[1,]=gsub("[^[:alpha:]]+"," ",templateProposalCV[1,])
  for(cv_col in 1:ncol(templateProposalCV)) {
    cv_name = templateProposalCV[1,cv_col]
    cv = templateProposalCV[,cv_col][-1]
    # clean UTF8-NB_space and other unprintable whitespaces
    cvClean = gsub("\u00A0"," ",cv) # UTF8-NBSP
    cvClean = gsub("[^[:alnum:][:punct:]]+"," ",cvClean) # anything else
    .GlobalEnv$cvList[[cv_name]]=c(cvClean[!is.na(cvClean)],NA)
    if(params$tmi) {cat("ProposalTemplateCV[",cv_name,"]: ", length(cvList[[cv_name]]), " from ",refProposalTemplateFilename,":",params$template_xlsx_sheet,"\n")}
    
  }
  
  # map to actual input xlsx column names
  cvNameMap = c(
    "Genome coverage"=    "genomeCoverage",
    "Genome composition"= "molecule",
    "Host Source"=        "hostSource",
    "Change"=             "change",
    "Proposed Rank"=      "rank",
    "Subcommittee Abbrev"="scAbbrev",
    "Subcommittee Name"=  "scName"
  )
  names(.GlobalEnv$cvList)=cvNameMap[names(cvList)]
  
  # action mappings - map XLSX's "change" to internal "action"
  change2actionMap = c(
    "create" = "new"
    ,"create new" = "new"
    ,"rename" = "rename"
    ,"abolish" = "abolish"
    # types of moves
    ,"move" = "move"
    ,"move; rename" = "move"
    ,"promote" = "promote"
    ,"demote" = "demote"
    ,"merge" = "merge"
    # split is treated as a move or a create
    # depending on whether the name changes
    ,"split" = "split"
  )
  .GlobalEnv$cvList[[".change2action"]]=change2actionMap
  
  # remove ("Please select",NA) from "change" & "rank" CVs - that is a required field
  for( cv in c("change","rank") ) { #},"scAbbrev","scName","hostSource") ) {
    isRemoveTerm = tolower(gsub("[^[:alnum:]]","",cvList[[cv]])) %in% c(
      "pleaseselect", # 2023
      NA
    ) 
    .GlobalEnv$cvList[[cv]] = cvList[[cv]][!isRemoveTerm]
  }
  
  # clean UTF8-NB_space (mac) 
  # pattern=paste0(AscToChar(194),AscToChar(160))
  for( cv in c("change","rank","scAbbrev","scName") ) {
    for( i in seq(1:length(cvList[[cv]])) ) {
      term = cvList[[cv]][i]
      if( is.na(term) || regexpr(text=term, pattern="\xC2\xA0") > 0 ) {
        term = gsub(paste0("(\xC2\xA0)+")," ",term)
        .GlobalEnv$cvList[[cv]][i] = term
      }
    }
  }
  
  #
  ##### build Subcommittee Map   ##### 
  #
  .GlobalEnv$scAbbrevNameMap = cvList[["scName"]]
  names(.GlobalEnv$scAbbrevNameMap) = cvList[["scAbbrev"]]
  # remove them from "official" cvList, as that list 
  # also servers as a list of require columns for the proposal xlsx
  .GlobalEnv$cvList = cvList[!names(cvList) %in% c("scName","scAbbrev")]

  #
  # add XLS terms be aliases to DB ID maps
  #
  cv = "molecule"
  isOk = 
    # things that match after removing spaces
    gsub(pattern=" ",replacement="",x=cvList[[cv]]) %in% names(dbCvMapList[[cv]]) &
    # things that DON'T match exactly
    !cvList[[cv]] %in% names(dbCvMapList[[cv]])
  
  aliases = dbCvMapList[["molecule"]][gsub(pattern=" ",replacement="",x=cvList[["molecule"]])[isOk]]
  names(aliases) = cvList[[cv]][isOk]
  # add aliases
  dbCvMapList[[cv]] = c(
    # original map
    dbCvMapList[[cv]],
    # additional aliases
    aliases,
    # hack for "multiple" which was in xlsx, but doesn't map to the db
    "multiple"=NA
  )

  # 
  ##### replace rank CV from db   ##### 
  #
  #cvList[["rank"]] = [["rank"]]
  
  # 
  ##### augment hostSource CV from db   ##### 
  #
  # this is just names, no IDs
  #
  dbHostSourceFilename=file.path(params$ref_dir, params$db_host_source_fname)
  hostSourceCV = read.delim(file=dbHostSourceFilename,header=TRUE, stringsAsFactors=TRUE,na.strings=c("NULL"))
  #rownames(hostSourceCV) = hostSourceCV$id
  .GlobalEnv$dbCvList[["hostSource"]] = hostSourceCV$host_source
  .GlobalEnv$dbCvMapList[["hostSource"]] = hostSourceCV$host_source
  names(.GlobalEnv$dbCvMapList[["hostSource"]]) = hostSourceCV$host_source
  if(params$verbose) {cat("HostSourceCV: ", dim(hostSourceCV), " from ",dbHostSourceFilename,"\n")}
  .GlobalEnv$cvList[["hostSource"]] = union(cvList[["hostSource"]],dbCvList[["hostSource"]])
  #
  ##### CVs from VMR xlsx  #####
  #
  vmrFilename = file.path(params$ref_dir, params$vmr_fname)
  .GlobalEnv$vmrDf = data.frame(
    # suppress errors about column names
    # https://github.com/tidyverse/readxl/issues/580#issuecomment-519804058
    suppressMessages(
      read_excel(
        path    = vmrFilename,
        #sheet   = "Terms",
        trim_ws = TRUE,
        na      = c("Please select","[Please select]","[Please\u00A0select]"),
        skip    = 0,
        range   = cell_cols("A:AO"),
        col_names = TRUE
      )
    )
  )
} # load_reference()



#### load ref ####

if(params$use_cache && !params$update_cache && file.exists(cacheFilename)) {
  load_reference_cache()
} else {
  load_reference()
}
if( params$update_cache ) {
  save_reference_cache()
}

#### SCAN FOR PROPOSALS ####

scan_for_proposals = function() {
  if( params$debug ) {cat("SCAN_FOR_PROPOSALS(",params$proposals_dir,")\n")}
    # 
  ##### filename regex #####
  #
  if( params$processing_mode == "final") {
    filenameFormatRegex="^[0-9][0-9][0-9][0-9]\\.[0-9][0-9][0-9][A-Z]X*\\.[^ ]*"
    filenameFormatMsg="final:####[A-Z].###[A-Z][X].____"
  } else if( params$processing_mode == "draft") {
    filenameFormatRegex="^[0-9][0-9][0-9][0-9]\\.[0-9][0-9][0-9][A-Z]X*\\.[A-Za-z]+\\.v[0-9]+\\.[^ ]*"
    filenameFormatMsg="draft:####[A-Z].###[A-Z][X].[A-Z]+*.v#.____"
  } else if( params$processing_mode == "validate") {
    # allow any doc/xls file
    filenameFormatRegex="^.*"
    filenameFormatMsg="validate:*.____"
  } else {
    # unsupported format
    cat(paste0("ERROR: --mode='",params$processing_mode,"' is not a valid option: validate, draft, or final\n"))
    quit(save="no", status=1)
  }
  if(params$verbose){ cat("(PROCESSING) MODE      :",params$processing_mode,"\n") }
  
  ##### scan proposal_dir ##### 
  
  if( file_test("-d",params$proposals_dir) ) {
    # recursively scqan directory for .doc/.xls files 
    inputFiles = data.frame(docpath=list.files(path=params$proposals_dir,
                                               pattern="^[^~.].*\\.(doc|xls)x*$", 
                                               recursive=T, full.names=TRUE)
    )
  } else if( file_test("-f",params$proposals_dir) ) {
    # the "path" is actually a file
    inputFiles = data.frame(docpath=c(params$proposals_dir))
  } else {
    errorDf = data.frame(notes=paste0("Input folder='",params$proposals_dir,"/'"))
    errorDf$level = "ERROR"
    errorDf$error = "PROPOSAL_DIR_NO_EXIST"
    errorDf$message = "path not found"
    .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
    write_error_summary(.GlobalEnv$allErrorDf)
    cat("# ERROR: ",errorDf$error, ":", params$proposals_dir, "\n")
    quit(save="no", status=1)
  }
  # debugging - only process the regex-matching files
  #inputFiles = inputFiles[grep(inputFiles$docpath,pattern="2023.013P"),,drop=FALSE]
  
  # list files found, if in TMI mode
  if(params$tmi) { 
    cat("# xls|doc(x) files found: N=",nrow(inputFiles),"\n")
    for( f in inputFiles$docpath ) {
      cat("\t",f,"\n")
    }
  }
  
  # check if we found any files at all
  if( nrow(inputFiles) == 0 ) {
    errorDf = data.frame(notes=paste0("Input folder='",params$proposals_dir,"/'"))
    errorDf$level = "ERROR"
    errorDf$error = "NO_INPUT_FILES"
    errorDf$message = "found no .xls[x] or .doc[x] files"
    .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
    write_error_summary(.GlobalEnv$allErrorDf)
    cat("# ERROR: NO_INPUT_FILES found in ",params$proposals_dir, "\n")
    quit(save="no", status=1)
  }
  
  ##### parse file paths/names #####
  inputFiles$path         = dirname(inputFiles$docpath)
  inputFiles$file         = basename(inputFiles$docpath)
  inputFiles$basename     = gsub("(.*).(doc|xls)x*$","\\1",inputFiles$file)
  # code may not be in early versions of filenames
  parseableFilenames = grep(inputFiles$basename,pattern="^[0-9][0-9][0-9][0-9]\\.[0-9][0-9][0-9][A-Z]X*\\.")
  inputFiles[parseableFilenames,"code"]          = sub("^([0-9]+\\.[0-9]+[A-Z]X*).*","\\1",inputFiles[parseableFilenames,]$basename)
  inputFiles[parseableFilenames,"scAbbrev"]      = sub("^[0-9][0-9][0-9][0-9]\\.[0-9][0-9][0-9]([A-Z])X*$","\\1",inputFiles[parseableFilenames,"code"] )
  
  ##### filter filenames #####
  # remove all *.Ud.* files, when in draft/final modes
  if( params$processing_mode %in% c("draft","final") ) {
    filtered = grep(inputFiles$file, pattern=".*\\.Ud\\..*")
    if( length(filtered)) {
      errorDf = inputFiles[filtered,]
      errorDf$level = "INFO"
      errorDf$error = "IGNORE_FNAME_.Ud."
      errorDf$message = "All files with '.Ud.' in filename are ignored"
      .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
      write_error_summary(.GlobalEnv$allErrorDf)
      if(params$tmi) { cat(paste0("# .Ud. files filtered out: N=",length(filtered),"\n"))}
    }
    inputFiles = inputFiles[grep(inputFiles$file, pattern=".*\\.Ud\\..*", invert=T),]
    if(params$tmi) { cat("# xls|doc(x) files after .Ud. removal: N=",nrow(inputFiles),"\n")}
  }
  
  # ignore "Suppl" files 
  filtered = grep(inputFiles$file,pattern=params$infile_suppl_pat)
  if( length(filtered) ) {
    errorDf = inputFiles[filtered,]
    errorDf$level = "INFO"
    errorDf$error = "IGNORE_FNAME_SUPPL"
    errorDf$message = paste0("All files with '",params$infile_suppl_pat,"' in filename are ignored")
    .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
    write_error_summary(.GlobalEnv$allErrorDf)
    if(params$tmi) { cat(paste0("# Suppl files filtered out: N=",length(filtered),"\n"))}
  }
  inputFiles = inputFiles[grep(inputFiles$file,pattern=params$infile_suppl_pat, invert=T),]
  if(params$tmi) { cat("# xls|doc(x) files after Suppl/appendix removal: N=",nrow(inputFiles),"\n")}
  
  # ignore "appendix" files 
  filtered = grep(inputFiles$file,pattern="appendix",ignore.case = T)
  if( length(filtered) ) {
    errorDf = inputFiles[filtered,]
    errorDf$level = "INFO"
    errorDf$error = "IGNORE_FNAME_APPENDIXL"
    errorDf$message = paste0("All files with '","appendix","' in filename are ignored")
    .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
    write_error_summary(.GlobalEnv$allErrorDf)
    if(params$tmi) { cat(paste0("# Appendix files filtered out: N=",length(filtered),"\n"))}
  }
  inputFiles = inputFiles[grep(inputFiles$file,pattern="appendix",ignore.case = T, invert=T),]
  if(params$tmi) { cat("# xls|doc(x) files after appendix removal: N=",nrow(inputFiles),"\n")}
  
  ##### SECTION scan DOCX #####
  #
  # break filename into proposal code, filename and basename
  #
  docxs          = inputFiles[grep(inputFiles$file,pattern="\\.docx*$"),,drop=FALSE]
  colnames(docxs)[which(colnames(docxs) %in% c("docpath","file") )] <- c("docxpath","docx")
  #names(docxs)   = c("docxpath","path","docx","basename","code","scAbbrev")
  
  if(params$tmi) { cat("# doc(x) files found: N=",nrow(docxs),"\n")}
  
  # QQQ use log_error()!!!
  
  # check for duplicate Proposal IDs
  dups = duplicated(docxs$code)
  allDups =docxs$code %in% docxs$code[dups]
  if(sum(dups) > 0) {
    errorDf = docxs[allDups, c("scAbbrev", "code", "docx")]
    errorDf$level = "ERROR"
    errorDf$error = "DUPCODE.DOCX"
    errorDf$message = "duplicate proposal ID"
    .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
    # output: @DD make this web-app friedly
    #kable(errorDf,caption = paste0("QC01: ERRORS: dupliate docx proposal IDs"))
    write_xlsx(x=errorDf,path=file.path(params$out_dir,"QC01.docx_duplicate_ids.xlsx"))
    write_error_summary(.GlobalEnv$allErrorDf)
  }
  # use codes for row names, but fall back to ordinals, if no code in filename
  rownames(docxs)=ifelse(is.na(docxs$code),rownames(docxs),docxs$code)
  docxs[,"code"] = rownames(docxs)
  
  # spaces in filenames - move this out of DOCX/XLSX and into common
  spacedOut = grep(pattern=" ",docxs$docx)
  if(sum(spacedOut) > 0) {
    errorDf = docxs[spacedOut, c("scAbbrev", "code", "docx")]
    errorDf$level = "WARNING"
    errorDf$error = "DOCX_FILENAME_SPACES"
    errorDf$message = "filename contains a space: please replace with _ or -"
    errorDf$notes = gsub("( )","[\\1]",docxs[spacedOut,]$docx)
    .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
  }
  
  docxBadFnameFormats = grep(docxs$docx, pattern=paste0(filenameFormatRegex,".docx$"), invert=T)
  if( length(docxBadFnameFormats) > 0 ) {
    if( params$processing_mode %in% c("draft","final") ) {
      # supress picky error message
      errorDf= docxs[docxBadFnameFormats,c("scAbbrev","code","docx")]
      errorDf$level= "WARNING"
      errorDf$error = "DOCX.BAD_FILENAME_FORMAT"
      errorDf$message = paste0("Should be '",filenameFormatMsg,".docx'")
      errorDf$notes= "####[A-Z]=year/study_section, ###[A-Z]=index/type, [A-Z]+=status, v#=version"
      .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
    }
  }
  #
  ##### SECTION: scan XLSX ##### 
  #
  
  xlsxs          = inputFiles[grep(inputFiles$file,pattern="\\.xlsx*$"),,drop=FALSE]
  colnames(xlsxs)[which(colnames(xlsxs) %in% c("docpath","file") )] <- c("xlsxpath","xlsx")
  #names(xlsxs)   = c("xlsx path","path","xlsx","basename","code","scAbbrev")
  
  if(params$tmi) { cat("# xls(x) files found: N=",nrow(xlsxs),"\n")}
  
  # QC for duplicate codes
  dups = duplicated(xlsxs$code) & !is.na(xlsxs$code)
  allDups =xlsxs$code %in% xlsxs$code[dups]
  if( sum(dups) > 0 ) {
    # track if all dups were resolved by being non-template files
    dups_ok = TRUE 
    filter_rows = NULL
    
    if(params$tmi) { cat("START: can not proceed with duplicated proposal IDs:", paste(levels(as.factor(xlsxs[dups,"code"])),collapse=","),"\n") }
    # error details
    errorDf = xlsxs[allDups, c("scAbbrev", "code", "xlsx")]
    errorDf$level = "ERROR"
    errorDf$error = "XLSX_DUPCODE"
    for( code in xlsxs$code[dups] ) {
      if(params$tmi) {
        print("====================================================================================" )
        print(paste("====", code, ":  DUP XLSXS ===="))
        print("====================================================================================" )
      }
      for( row in rownames(xlsxs[xlsxs$code==code,]) ) {
        if(params$tmi) { 
          print("------------------------------------------------------------------------------------" )
          print(paste("D5UP_XLSX_FORMAT_CHECK: ", xlsxs[row,"xlsx"] ))
        }
        load_check = load_proposal(code,  xlsxs[row,"xlsxpath"] )
        if(is.null(load_check[["proposalDf"]])) {
          # report the non-proposal non-suppl-named xlsx file
          if(params$tmi) { print(paste0("LOAD_CHECK_RESULT:  NOT PROPOSAL: ", xlsxs[row,"xlsx"])) }
          subErrorDf = xlsxs[row,]
          subErrorDf[,c("docpath","file")] = subErrorDf[,c("xlsxpath","xlsx")]
          subErrorDf$level = "WARNING"
          subErrorDf$error = "NON_PROPOSAL_NOT_SUPPL"
          subErrorDf$message = paste0("Note a proposal XLSX: All supplemental XLSX files should have '_Suppl.' in the name!")
          .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, subErrorDf),fill=TRUE)
          write_error_summary(.GlobalEnv$allErrorDf)
          if(params$tmi) { cat(paste0("# DUP file filtered out (not a template): code=",code,", file=",xlsxs[row,"xlsx"],"\n"))}
          # remove it from dups list
          allDups[row]= FALSE
          filter_rows=c(filter_rows,row)
        } else {
          if(params$tmi) { print(paste0("LOAD_CHECK_RESULT:  *IS PROPOSAL: ", xlsxs[row,"xlsx"])) }
        }
      } # for each dup for this code
      
      # check if there was only one parsable xlsx
      code_template_count = nrow(xlsxs[xlsxs$code==code && allDups,])
      if( code_template_count > 1 ) {
          dups_ok = FALSE
      }
      # build list of all xlsxs using duplicate codes
      errorDf[errorDf$code==code,"message"] = 
        paste0("duplicate proposal ID: ",
               paste(xlsxs$xlsx[xlsxs$code==code],collapse=","))
    } # for each code with dups
    #
    # finalize
    #
    # append to global error list
    .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
    # output
    write_error_summary(.GlobalEnv$allErrorDf)
    
    if( !dups_ok ) {
      # terminate
      cat("ERROR: can not proceed with duplicated proposal IDs:", paste(levels(as.factor(xlsxs[dups,"code"])),collapse=","),"\n")
      stop(1)
    } else {
      # remove the xlsxs rows that got filtered 
      xlsxs = xlsxs[!(rownames(xlsxs) %in% filter_rows),]
    }
  }
  rownames(xlsxs) = ifelse(is.na(xlsxs$code),rownames(xlsxs),xlsxs$code)
  xlsxs[,"code"] = rownames(xlsxs)
  #
  # check that xlsx names match format - move this out of DOCX/XLSX and into common
  #
  # production format (no version)
  xlsxBadFnameFormats = grep(xlsxs$xlsx, pattern=paste0(filenameFormatRegex,".xlsx*$"), invert=T)
  if( length(xlsxBadFnameFormats) > 0 ) {
    if( params$processing_mode %in% c("draft","final") ) {
      # supress picky error message
      errorDf= xlsxs[xlsxBadFnameFormats,c("scAbbrev","code","xlsx")]
      errorDf$level= "WARNING"
      errorDf$error = "XLSX.BAD_FILENAME_FORMAT"
      errorDf$message = paste0("Should be '",filenameFormatMsg,".xlsx'")
      errorDf$notes= "####[A-Z]=year/study_section, ###[A-Z]=index/type, [A-Z]+=status, v#=version"
      .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
    }
    write_error_summary(.GlobalEnv$allErrorDf)  
  }
  # spaces in XLSX filenames
  spacedOut = grep(pattern=" ",xlsxs$xlsx)
  if(sum(spacedOut) > 0) {
    errorDf = xlsxs[spacedOut, c("scAbbrev", "code", "xlsx")]
    errorDf$level = "WARNING"
    errorDf$error = "XLSX_FILENAME_SPACES"
    errorDf$message = "filename contains a space: please replace with _ or -"
    errorDf$notes = gsub("( )","[\\1]",xlsxs[spacedOut,]$xlsx)
    .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
    write_error_summary(.GlobalEnv$allErrorDf)  
  }
  #
  # check if any filenames passed QC
  #
  # check if we found any files at all
  if( nrow(xlsxs) == 0 ) {
    errorDf = data.frame(notes=paste0("Input folder='",params$proposals_dir,"/'"))
    errorDf$level = "ERROR"
    errorDf$error = "NO_INPUT_FILES"
    errorDf$message = "found no .xls[x] files that passed filename QC"
    .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
    write_error_summary(.GlobalEnv$allErrorDf)
    cat("# ERROR: NO_INPUT_FILES found in ",params$proposals_dir, "\n")
    quit(save="no", status=1)
  }
  
  #
  ##### merge XLSX list into DOCX list #####
  #
  
  proposalsDf = merge(xlsxs, docxs, by=c("code","path","scAbbrev"), 
                      all=T, suffixes = c(".xlsx",".docx"))
  # prioritize the xlsx basename for zip file naming and reporting.
  proposalsDf$basename = ifelse(!is.na(proposalsDf$basename.xlsx),proposalsDf$basename.xlsx, proposalsDf$basename.docx)
  # remove xlsx/docx specific basenames
  proposalsDf = subset(proposalsDf, select = -c(basename.xlsx, basename.docx) ) 
  
  rownames(proposalsDf)=proposalsDf$code
  
  # strip off version, workflow status and .fix, to get final, production filename
  proposalsDf$cleanbase= gsub("^([0-9]+\\.[0-9]+[A-Z]X*)(\\.[A-Z]+)(\\.v[0-9]+)*(\\.fix)*(\\..*)$","\\1\\5",proposalsDf$basename)
  
  # QC  - missing xlsx file
  missing= is.na(proposalsDf$xlsx)
  if( sum(missing) > 0 ) {
    errorDf= proposalsDf[missing,c("code","docx")]
    # suggest possible matches based on ID
    errorDf$xlsx= NA 
    errorDf$row = NA
    errorDf$level= "ERROR"
    errorDf$error = "XLSX.MISSING"
    errorDf$message = "DOCX has no matching XLSX"
    errorDf$notes= "Suggestions: contact corresponding author"
    for(row in rownames(errorDf) ) {
      guesses = sum(xlsxs$code==errorDf[row,"code"],na.rm=TRUE)
      if( guesses == 1 ) { 
        # if we have a unique guess, then make it a WARNING, and use that guess
        proposalsDf[row,"xlsx"]=xlsxs[xlsxs$code==errorDf[row,"code"],"xlsx"]
        proposalsDf[row,"xlsxpath"]=xlsxs[xlsxs$code==errorDf[row,"code"],"xlsxpath"]
        errorDf[row,]$xlsx = proposalsDf[row,"xlsx"]
        errorDf[row,]$row = NA
        errorDf[row,]$level = "WARNING"
        errorDf[row,]$error = "xlsx.TYPO"
        errorDf[row,]$notes = paste("Using best guess:",proposalsDf[row,"xlsx"])
      } else if( guesses > 1 ) {
        # just list all the options we found
        errorDf[row,]$error = "xlsx.MULTIPLE"
        errorDf[row,]$error = paste0("Multiple .xlsx files start with code '",errorDf[row,"code"],"', and aren't marked '",params$infile_suppl_pat,"'")
        errorDf[row,"notes"] = paste("SUGGESTIONS:",paste(paste0(xlsxs[xlsxs$xlsxID==errorDf[row,"docxID"],"basename"],".xlsx"),collapse=", "))
      } else {
        # 
        # can't find the xlsx
        #
        # let the error pass through
        proposalsDf[]
      }
    }
    # suppress name-mismatch warnings
    loadErrorDfFilt = errorDf
    if( params$processing_mode == "validate" ) {
      # filter out some errors
      loadErrorDfFilt = errorDf %>% filter(error != "xlsx.TYPO")
    }
    if(nrow(loadErrorDfFilt) > 0) {
      # append to global list
      .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
    }
    write_error_summary(.GlobalEnv$allErrorDf)  
  }
  # QQQ don't check for missing docx files? 
  
  #
  # get SC names from last letter of code
  #
  proposalsDf$scAbbrev = NA
  parsableCode = grep(proposalsDf$code, pattern="[0-9][0-9][0-9][0-9]\\.[0-9][0-9][0-9]([A-Z])X*" )
  if( length(parsableCode) >0 ) {
    proposalsDf[parsableCode,"scAbbrev"] =  gsub("[0-9][0-9][0-9][0-9]\\.[0-9][0-9][0-9]([A-Z])X*","\\1",proposalsDf$code[parsableCode])
  }
  # QC
  badProposalAbbrevs = !(proposalsDf$scAbbrev %in% names(scAbbrevNameMap))
  if(sum(badProposalAbbrevs)>0) {
    if( params$processing_mode %in% c("draft","final") ) {
      # supress picky error message
      errorDf = proposalsDf[badProposalAbbrevs, c("scAbbrev", "code", "xlsx", "docx")]
      errorDf$level = "WARNING"
      errorDf$error = "CODE_BAD_SC_ABBREV"
      errorDf$message = "Last letter of CODE not a valid ICTV Subcommittee letter"
      errorDf$notes = paste0("'",proposalsDf[badProposalAbbrevs,"scAbbrev"],"' not in [",
                             paste0(names(scAbbrevNameMap),collapse=","),"]")
      .GlobalEnv$allErrorDf = rbindlist(list(.GlobalEnv$allErrorDf, errorDf),fill=TRUE)
    }
    write_error_summary(.GlobalEnv$allErrorDf)  
  }
  proposalsDf$subcommittee = ifelse(badProposalAbbrevs,
                                    ifelse(is.na(proposalsDf$scAbbrev),"unspecified",paste0("unknown-",proposalsDf$scAbbrev)),
                                    scAbbrevNameMap[proposalsDf$scAbbrev])

  # move results to global scope
  .GlobalEnv$proposalsDf = proposalsDf
  return(proposalsDf)
} # scan_for_proposals()

####  LOAD/QC #### 

##### proposal structure/mapping tables #####
#
# this uses the PROPOSAL.XLSX schema (naming convention)
#

# dput(unname(as.vector(df[1,])))
xlsx_v1_row2=c(
  "CURRENT TAXONOMY", NA_character_, NA_character_, NA_character_, NA_character_, 
  NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, 
  NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, 
  NA_character_, "PROPOSED TAXONOMY", NA_character_, NA_character_, NA_character_, 
  NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, 
  NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, 
  NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, 
  NA_character_, NA_character_, NA_character_, "SPECIFY PROPOSED CHANGE", NA_character_, 
  "COMMENTS", NA, NA, NA, NA_character_, 
  NA_character_, NA_character_, NA_character_, NA_character_
)

xlsx_v2_row2=c(
  "CURRENT TAXONOMY", NA_character_, NA_character_, NA_character_, NA_character_, 
  NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, 
  NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, 
  "PROPOSED TAXONOMY", NA_character_, NA_character_, NA_character_, 
  NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, 
  NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, 
  NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, 
  NA_character_, NA_character_, NA_character_, "SPECIFY PROPOSED CHANGE", NA_character_, 
  "COMMENTS"
)

xlsx_2023_row3 = c("CURRENT TAXONOMY" , NA_character_, NA_character_ , NA_character_, NA_character_,
                   NA_character_, NA_character_, NA_character_, NA_character_, NA_character_,
                   NA_character_, NA_character_, NA_character_, NA_character_, NA_character_,
                   "PROPOSED TAXONOMY", NA_character_, NA_character_, NA_character_, NA_character_,
                   NA_character_, NA_character_, NA_character_, NA_character_, NA_character_,
                   NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, 
                   "DESCRIPTIVES", NA_character_, NA_character_, NA_character_, NA_character_,
                   NA_character_, NA_character_, "ACTION", NA_character_, "COMMENTS"
)
# dput(unname(as.vector(df[2,])))
xlsx_v1_row3=c(
  "Realm", "Subrealm", "Kingdom", "Subkingdom", 
  "Phylum", "Subphylum", "Class", "Subclass", "Order", "Suborder", 
  "Family", "Subfamily", "Genus", "Subgenus", "Species",
  # V1 only, removed in V2
  "Exemplar GenBank Accession Number", 
  "Realm", "Subrealm", "Kingdom", "Subkingdom", "Phylum", "Subphylum", 
  "Class", "Subclass", "Order", "Suborder", "Family", "Subfamily", 
  "Genus", "Subgenus", "Species",
  "Exemplar GenBank Accession Number", 
  "Exemplar \r\nvirus name",
  "Virus name abbrevn", 
  "Exemplar\r\nisolate designation", 
  "Genome coverage", 
  "Genome composition", 
  "Host/Source", 
  "Change", 
  "Rank"
)
xlsx_v2_row3=c(
  # V2 added column
  Code=NA_character_, 
  "Realm", "Subrealm", "Kingdom", "Subkingdom",
  "Phylum","Subphylum", "Class", "Subclass", "Order", "Suborder", 
  "Family", "Subfamily", "Genus", "Subgenus", "Species",
  "Realm", "Subrealm", "Kingdom", "Subkingdom", "Phylum","Subphylum",
  "Class", "Subclass", "Order", "Suborder", "Family", "Subfamily",
  "Genus", "Subgenus", "Species",
  "Exemplar GenBank Accession Number", 
  "Exemplar \r\nvirus name", 
  "Virus name abbreviation", 
  "Exemplar\r\nisolate designation", 
  "Genome coverage", 
  "Genome composition", 
  "Host/Source", 
  "Change", 
  "Rank",
  # V2 added column
  "Comments"
)

xlsx_2023_row4=c(
  "Realm", "Subrealm", "Kingdom", "Subkingdom", 
  "Phylum", "Subphylum", "Class", "Subclass", "Order", "Suborder", 
  "Family", "Subfamily", "Genus", "Subgenus", "Species",
  "Realm", "Subrealm", "Kingdom", "Subkingdom", "Phylum", "Subphylum", 
  "Class", "Subclass", "Order", "Suborder", "Family", "Subfamily", 
  "Genus", "Subgenus", "Species",
  "Exemplar GenBank Accession Number", 
  "Exemplar virus name",
  "Virus name abbreviation", 
  "Exemplar isolate designation", 
  "Genome coverage", 
  "Genome composition", 
  "Host/Source", 
  "Change", 
  "Proposed Rank",
  "Comments"
)
#
# changeDf column names
#
# this is the mapped subset of columns pulled from the proposal
# normalized across xlsx template versions
#
xlsx_change_ranks=c("Realm", "Subrealm", "Kingdom", "Subkingdom", 
                    "Phylum", "Subphylum", "Class", "Subclass", "Order", "Suborder", 
                    "Family", "Subfamily", "Genus", "Subgenus", "Species"
)

# record mapping of V2 column names to generic v1/v2 merged, cleaned names
xlsx_change_other=c( 
  "Exemplar GenBank Accession Number" = "exemplarAccession", 
  "Exemplar \r\nvirus name"           = "exemplarName", 
  "Virus name abbreviation"           = "Abbrev", 
  "Exemplar\r\nisolate designation"   = "exemplarIsolate", 
  "Genome coverage"                   = "genomeCoverage", 
  "Genome composition"                = "molecule",   # Genome Composition", 
  "Host/Source"                       = "hostSource",
  "Change"                            = "change",
  "Rank"                              = "rank",
  "Comments"                          = "comments"
)
xlsx_change_colnames = c(  
  paste0("src",xlsx_change_ranks), # 1:15
  tolower(xlsx_change_ranks),      # 16:30
  xlsx_change_other )              # 31:39

xlsx_change_srcCols = c(1:15)
xlsx_change_destCols = c(16:30)

xlsx_change_srcDest_colnames = xlsx_change_colnames[c(xlsx_change_srcCols,xlsx_change_destCols)]

#
# map xlsx VERSION columns to normalized "Change" columns
#
xlsx_v1_change_cols = c(
  1:15,  # srcRanks
  17:31, # dest Ranks
  32:41  # other
)
xlsx_v2_change_cols = c(
  2:16,  # srcRanks
  17:31, # dest Ranks
  32:41  # other
)
xlsx_2023_change_cols =  c(
  1:15,  # srcRanks
  16:30, # dest Ranks
  31:40  # other
)

##### regex list value QC lists #####
#
# regex patterns to recognize illegal cell content
#
#xlsx_col_info = data.frame(row.names=xlsx_change_colnames,
#                           pattern = rep(NA_character_,length(xlsx_change_colnames)),
#                           pat_warn = rep(NA_character_,length(xlsx_change_colnames))
#                           )

value_validation = data.frame(
  pat_name = "remove-non-alpha-numeric",
  # non-species taxa
  col = grep(xlsx_change_colnames[c(xlsx_change_srcCols,xlsx_change_destCols)],pattern="pecies",invert = T, value=T),
  regex = "([^[:alnum:]-]+)",
  type = "replace",
  class = "INFO",
  code = "XLSX.NON_ALPHA-NUMERIC",
  warn = "non-(AlphaNumeric or hyphen) characters removed",
  replace = ""
)

value_validation = rbind(value_validation, data.frame(
  pat_name = "remove-non-alpha-numeric-space",
  col = c("srcSpecies","species"),
  regex = "([^[:alnum:] -]+)",
  type = "replace",
  class = "INFO",
  code = "XLSX.NON_ALPHA-NUMERIC-SPACE",
  warn = "non-(AlphaNumeric or hyphen) characters removed",
  replace = ""
))

# Species name must be "word[space]word"
value_validation = rbind(value_validation, data.frame(
  pat_name = "species 2 words with a space",
  col = c("species"),
  regex = "^([[:alnum:]-]+ [[:alnum:] -]+)$",
  type = "required",
  class = "ERROR",
  code = "XLSX.SPECIES_BAD_NAME",
  warn = "Species name must be 'genus[space]species' binomial naming",
  replace = ""
))

# allow internal spaces on species and columns with lists/isolate names
value_validation = rbind(value_validation, data.frame(
  pat_name = "change.remove-non-alpha-numeric-space-semi",
  col = c("change"),
  regex = "([^[:alnum:] ;-]+)",
  type = "replace",
  class = "INFO",
  code = "XLSX.CHANGE_CV_REMOVE_NON_ALPHA-NUMERIC-SPACE-SEMI",
  warn = "Change col:non-(AlphaNumeric,hyphen,space,dash) characters removed",
  replace = " "
))

# Abbrev: semi-colon separated lists
value_validation = rbind(value_validation, data.frame(
  pat_name = "abbrev.remove-non-semicolon-sep-list",
  col = c("Abbrev"),
  regex = "([^[:alnum:];,_/. -]+)",
  type = "replace",
  class = "INFO",
  code = "XLSX.ABBREV.REMOVE_ILLEGAL_CHARS",
  warn = "Should be semicolon-separated list",
  replace = ""
))

# Accession:  convert " and " to ";" in lists
value_validation = rbind(value_validation, data.frame(
  pat_name = "accession.and-to-semi",
  col = c("exemplarAccession"),
  regex = "( +and +)",
  type = "replace",
  class = "INFO",
  code = "XLSX.ACCESSION.REPLACE_AND_WITH_SEMI",
  warn = "Should be semicolon-separated list, with optional colon-separated name prefixes",
  replace = ";"
))

value_validation = rbind(value_validation, data.frame(
  pat_name = "accession.acc-paren-to-name-acc",
  col = c("exemplarAccession"),
  regex = "([A-Z0-9.]+) *\\(([^\\)]+) segment\\)",
  type = "replace",
  class = "INFO",
  code = "XLSX.ACCESSION.REPLACE_PAREN_WITH_PREFIX",
  warn = "Should be semicolon-separated list, with optional colon-separated name prefixes",
  replace = "\\2:\\1"
))

# Accession:  semi-colon separated list, with colon-separated labels
value_validation = rbind(value_validation, data.frame(
  pat_name = "accession.remove-non-semi-colon-sep-list",
  col = c("exemplarAccession"),
  regex = "([^[:alnum:]:;_/. -]+)",
  type = "replace",
  class = "INFO",
  code = "XLSX.ACCESSION.REMOVED_ILLEGAL_CHARS",
  warn = "Should be semicolon-separated list, with optional colon-separated name prefixes",
  replace = ""
))


# Accession:  semi-colon separated list, with colon-separated labels
value_validation = rbind(value_validation, data.frame(
  pat_name = "accession.remove-non-semi-colon-sep-list",
  col = c("exemplarName","exemplarIsolate","genomeCoverage","molecule","hostSource","comments"),
  regex = "([^[:print:]]+)",
  type = "replace",
  class = "INFO",
  code = "XLSX.NON_PRINTABLE_REMOVED",
  warn = "unprintable/non-ASCII removed",
  replace = " "
))

# exemplarIsolate and exemplarName are anything printable. 
#xlsx_col_info[c("exemplarName","exemplarIsolate"),"pattern"] = "([^[:print:]]+)"
#xlsx_col_info[c("exemplarName","exemplarIsolate"),"pat_warn"] = "unprintable/non-ASCII"
#xlsx_col_info[c("exemplarName","exemplarIsolate"),"pat_code"] = "XLSX.NON_PRINTABLE"
# most taxa can only be alphanumeric plus a hyphen
# xlsx_col_info[c(xlsx_change_colnames[c(xlsx_change_srcCols,xlsx_change_destCols)],"rank"),"pattern"] = "([^[:alnum:]-]+)"
# xlsx_col_info[c(xlsx_change_colnames[c(xlsx_change_srcCols,xlsx_change_destCols)],"rank"),"pat_warn"] = "non-(AlphaNumeric or hyphen)"
# xlsx_col_info[c(xlsx_change_colnames[c(xlsx_change_srcCols,xlsx_change_destCols)],"rank"),"pat_code"] = "XLSX.NON_ALPHA-NUMERIC"
# xlsx_col_info[c(xlsx_change_colnames[c(xlsx_change_srcCols,xlsx_change_destCols)],"rank"),"pat_replace"] = ""
# # allow internal spaces on species and columns with lists/isolate names
# xlsx_col_info[c("srcSpecies","species"),"pattern"] = "([^[:alnum:]]$+)"
# xlsx_col_info[c("srcSpecies","species"),"pat_warn"] = "non-(AlphaNumeric,hyphen,space)"
# xlsx_col_info[c("srcSpecies","species"),"pat_code"] = "XLSX.NON_ALPHA-NUMERIC-SPACE"
## allow internal spaces on species and columns with lists/isolate names
#xlsx_col_info[c("change"),"pattern"] = "([^[:alnum:] ;-]+)"
#xlsx_col_info[c("change"),"pat_warn"] = "non-(AlphaNumeric,hyphen,space,semicolon)"
#xlsx_col_info[c("change"),"pat_code"] = "XLSX.NON_ALPHA-NUMERIC-SPACE-SEMI"
# semi-colon separated lists
#xlsx_col_info[c( "Abbrev"),"pattern"] = "([^[:alnum:];_/ . -]+)"
#xlsx_col_info[c( "Abbrev"),"pat_warn"] = "Should be semicolon-separated list"
#xlsx_col_info[c( "Abbrev"),"pat_code"] = "XLSX.SEMICOLON-SEP-LIST"
#xlsx_col_info[c( "exemplarAccession"),"pattern"] = "([^[:alnum:]:;_/ .-]+)"
#xlsx_col_info[c( "exemplarAccession"),"pat_warn"] = "Should be semicolon-separated list of [name:]accession"
#xlsx_col_info[c( "exemplarAccession"),"pat_code"] = "XLSX.SEMICOLON-SEP-LIST-NAMED"
# exemplarIsolate and exemplarName are anything printable. 
#xlsx_col_info[c("exemplarName","exemplarIsolate"),"pattern"] = "([^[:print:]]+)"
#xlsx_col_info[c("exemplarName","exemplarIsolate"),"pat_warn"] = "unprintable/non-ASCII"
#xlsx_col_info[c("exemplarName","exemplarIsolate"),"pat_code"] = "XLSX.NON_PRINTABLE"

# others just must be printable
#otherCols = is.na(xlsx_col_info$pattern)
#xlsx_col_info[otherCols,"pattern"] = "([^[:print:]]+)"
#xlsx_col_info[otherCols,"pat_warn"] = "unprintable/non-ASCII"
#xlsx_col_info[otherCols,"pat_code"] = "XLSX.NON_PRINTABLE"

# ```
# # process changes function
# 
# ```{r process_changes_helper_functions, echo=FALSE}
#
# this uses the DATABASE schema for the taxonomy_node table (ie, naming convention)
#

realmSpeciesCols = 1:15

##### taxon disection functions #####

taxon_get_name = function(realmSpecies) {
  # 
  # get the NAME of the lowest ranked taxon
  #
  nonNAs = !is.na(realmSpecies[realmSpeciesCols])
  if( sum(nonNAs) > 0 ) {
    return(realmSpecies[max(which(nonNAs))])
  } else {
    # no taxon - could be the src of a create, or the dest of abolish
    #if( params$debug_on_error) { print("INTERNAL ERROR: taxon_get_name()=NA"); browser()}
    return(NA)
  }
}
taxon_get_parent_name = function(realmSpecies) {
  #
  # get the NAME of the penultimate low ranked taxon
  #
  nonNAs = !is.na(realmSpecies[realmSpeciesCols])
  if( sum(nonNAs) > 1 ) {
    # get 2nd to last rank 
    return(realmSpecies[sort(which(nonNAs),decreasing=T)[2]])
  } else if( sum(nonNAs) == 1 ) {
    # only one name means the parent is the tree root
    return(params$msl_name)
  } else {
    # no taxon - could be the src of a create, or the dest of abolish
    #if( params$debug_on_error) { print("INTERNAL ERROR: taxon_get_parent_name()=NA"); browser()}
    return(NA)
  }
}
taxon_get_lineage = function(realmSpecies,masks=c()) {
  #   
  # concat with ; all non-NA ranks, except those in masks
  #
  nonNAs = !is.na(realmSpecies[realmSpeciesCols])
  # is a matrix, rows=taxon, columns=ranks
  if(!is.null(masks) ) {
    # change to FALSE all the listed ranks
    nonNAs[tolower(xlsx_change_ranks) %in% tolower(masks)] = FALSE
  }
  if( sum(nonNAs) > 0 ) {
    return(paste0(realmSpecies[nonNAs],collapse=";"))
  } else {
    # no taxon - could be the src of a create, or the dest of abolish
    #if( params$debug_on_error) { print("INTERNAL ERROR: taxon_get_lineage()=NA"); browser()}
    return(NA)
  }
}
taxon_get_parent_lineage = function(realmSpecies) {
  nonNAs = !is.na(realmSpecies[realmSpeciesCols])
  if( sum(nonNAs) > 0 ) {
    rank=taxon_get_rank(realmSpecies)
    return(taxon_get_lineage(realmSpecies,masks=c(rank)))
  } else {
    # no taxon - could be the src of a create, or the dest of abolish
    #if( params$debug_on_error) { print("INTERNAL ERROR: taxon_get_parent_name()=NA"); browser()}
    return(NA)
  }
}
taxon_get_rank = function(realmSpecies) {
  nonNAs = !is.na(realmSpecies[realmSpeciesCols])
  if( sum(nonNAs) > 0 ) {
    # 
    # get the rank of the lowest-rank taxon 
    #
    return(as.character(dbCvList[["rank"]]$name[max(which(nonNAs))+1]))
  } else {
    # no taxon - could be the src of a create, or the dest of abolish
    #if( params$debug_on_error) { print("INTERNAL ERROR: taxon_get_parent_name()=NA"); browser()}
    return(NA)
  }
}

#
# error formating functions
#

# formatting
sdiffOpen = "["
sdiffClose = "]"
sdiffSep = "//"

diff_strings = function(s1, s2) {
  #cat("diff_strings(",s1,",",s2,")\n")
  # find left prefix
  lmatch = 0
  while( str_trunc(s1,lmatch+1,"right","")==str_trunc(s2,lmatch+1,"right","")) { 
    #cat("[",lmatch,"] ",str_trunc(s1,lmatch,"right",""),"==",str_trunc(s2,lmatch,"right",""),"\n");
    lmatch = lmatch + 1
  }
  # find right prefix
  rmatch = 0
  while( str_trunc(s1,rmatch+1,"left","")==str_trunc(s2,rmatch+1,"left","")) { 
    #cat("[",rmatch,"] ",str_trunc(s1,rmatch,"left",""),"==",str_trunc(s2,rmatch,"left",""),"\n");
    rmatch = rmatch + 1
  }
  # print substring analysis
  if((lmatch + rmatch) < 4) {
    # not very similar - just print whole words
    return(paste0(
      sdiffOpen, s1, sdiffSep, s2, sdiffClose
    ))
  }
  return(paste0(
    str_trunc(s1,lmatch,"right",""),
    sdiffOpen, 
    substr(s1, lmatch+1,str_length(s1)-rmatch),
    sdiffSep,
    substr(s2, lmatch+1,str_length(s2)-rmatch),
    sdiffClose,
    str_trunc(s1,rmatch,"left","")
  ))
  # debug
  # s1 = c("one2three","four5six"); s2 =c("onetwothree", "fourfivesix")
  # diff_strings(s1[1], s2[1])
  # diff_strings(s1[2], s2[2])
}

diff_lineages=function(lin1, lin2) {
  lin1ex=unlist(strsplit(as.character(lin1),";"))
  lin2ex=unlist(strsplit(as.character(lin2),";"))
  
  # match lengths - would be better if we knew the actual ranks, and could align them
  if(length(lin1ex) < length(lin2ex)) {
    lin1ex=c(lin1ex,rep("",length(lin2ex)-length(lin1ex)))
  } else if( length(lin1ex) > length(lin2ex)) {
    lin2ex=c(lin2ex,rep("",length(lin1ex)-length(lin2ex)))
  }
  
  diffs=lin1ex==lin2ex
  tups=c()
  for( i in 1:length(diffs)) {
    tups[i] = ifelse(diffs[i],lin1ex[i],diff_strings(lin1ex[i],lin2ex[i]))
  }
  #tups=ifelse(diffs,lin1ex,paste0('[[',lin1ex,'//',lin2ex,']]'))
  #tups=ifelse(diffs,lin1ex,diff_strings(lin1ex,lin2ex))
  
  return(paste0(
    ifelse(sum(!diffs)==0,"[identical]",""),
    paste0(tups,collapse=";")))
  
  # test
  # diff_lineages(paste0(s1,collapse=";"), paste0(s2,collapse=";"))
}

#  
##### load/QC functions #####

#
# load .xlsx file into DF
#
# NOTE on global variables (ick):
# when called from scan_proposal_dir:
#    ARG: xlsxpath_override will be set
#    GLOBAL: proposalsDf wont exist
#    GLOBAL: 
# 
load_proposal = function(code,xlsxpath_override=NA) {
  cat("# LOAD_PROPOSAL(",code,",",xlsxpath_override,")\n")
  
  # get xlsxpath, if not specified
  xlsxpath = ifelse(is.na(xlsxpath_override),
                    proposalsDf[code,"xlsxpath"],
                    xlsxpath_override)
                    
  # return data
  proposalDf = NULL
  
  # get worksheet names
  proposalsSheetNameRegex = "proposal*"
  sheetNames = suppressMessages(
    excel_sheets(xlsxpath)
  )
  #
  # compare to expected sets of sheet names, warn if they added extra sheets
  #
  sheets2022=c("Proposals Template","Menu Items (Do not change)")
  sheets2023=c("Instructions","Proposal Template","Menu Items (Do not change)")
  extraSheets=c()
  templateSheets=c()
  if( sum(sheets2022 %in% sheetNames) == length(sheets2022) ) {
    templateSheets = sheets2022
    extraSheets = setdiff(sheetNames,sheets2022)
  } else if( sum(sheets2023 %in% sheetNames) == length(sheets2023) ) {
    templateSheets = sheets2023
    extraSheets = setdiff(sheetNames,sheets2023)
  } 
  if( length(extraSheets > 1) ) {
    cat(paste0("code ", code, ": Extra worksheets found in xlsx: '",paste0(extraSheets,collapse="','"),"'\n"))
    if(is.na(xlsxpath_override)) {
      log_error(code,linenum=0,action="OPEN_XLSX",actionOrder=actionOrder, 
              rank="",taxon="",
              levelStr="WARNING",errorCode="XLSX_EXTRA_SHEETS",
              errorStr="XLS file has additional sheets not present in template",
              notes=paste0("Worksheet(s) named '",paste0(extraSheets,collapse="','"),"' were added ")
     )
    }
  }
  
  # Try and find the one sheet that contains the proposal data
  proposalsSheetNames = grep(proposalsSheetNameRegex,sheetNames,ignore.case = T, value=T)
  if( is.na(proposalsSheetNames[1]) ) {
    # ERROR if can't find it. 
    cat(paste0("code ", code, ": A worksheet name matching '",proposalsSheetNameRegex,"' was not found: ","logging error","\n"))
    if(is.na(xlsxpath_override)) {
      log_error(code,linenum=0,action="OPEN_XLSX",actionOrder=actionOrder, rank="",taxon="",
              levelStr="ERROR","XLSX_NOT_PROPOSAL",
              errorStr="XLS no acceptably named 'proposal' sheet",
              notes=paste0("A worksheet name matching '",proposalsSheetNameRegex,"' was not found, only: ",
                     "'",paste0(sheetNames,collapse="','"),"'"
              )
      )
    }
  } else if( length(proposalsSheetNames) > 1 ) {
    # ERROR if find more than one
    cat(paste0("code ", code, ": ERROR: More than 1 worksheet named like '",proposalsSheetNameRegex,"' found: '",paste0(proposalsSheetNames,collapse="','"),"'\n"))
    if(is.na(xlsxpath_override)) {
      log_error(code,linenum=0,action="OPEN_XLSX",actionOrder=actionOrder, rank="",taxon="",
        levelStr="ERROR","XLSX_MULTI_PROPOSAL", 
        errorStr="XLS file has more than one sheet named 'proposal'",
        notes=paste0("More than 1 worksheet named '",proposalsSheetNameRegex,"' were found: '",paste0(proposalsSheetNames,collapse="','"),"'")
      )
    }
  } else {
    # OK: found one and only one proposal sheet - read that one
    
    # read XLSX file, no column names
    proposalDf = data.frame(
      # suppress errors about column names
      # https://github.com/tidyverse/readxl/issues/580#issuecomment-519804058
      suppressMessages(
        read_excel(
          xlsxpath,
          sheet=proposalsSheetNames[1],
          trim_ws = TRUE,
          na = c("Please select","[Please select]","[Please\u00A0select]"),
          skip = 2,
          range = cell_cols("A:AO"),
          col_names = FALSE
        )
      )
    )
    proposalDf[,"code"]=code
  }
  
  # return value
  return(list(proposalDf=proposalDf))
  
}

#
# Load  .docx into DF
#
load_proposal_docx = function(code) {
  cat("# LOAD_PROPOSAL_DOCX(",code,")\n")

  # return data
  metaDf = data.frame(code=c(code))

  # check if we even HAVE a docx file!
  if( is.na(proposalsDf[code,"docxpath"]) ) {
    if( params$processing_mode %in% c("draft","final") ) {
      # supress picky error message
      log_error(code,linenum=0,action="OPEN_DOCX",actionOrder=actionOrder,
                rank="",taxon="",
                levelStr = "WARNING",
                errorCode = "DOCX_MISSING",
                errorStr = ".docx file not available")
    }
    return(list(metaDf=metaDf))
  }
  
  # scanner will find .doc files, but we only support .docx right now.
  docFileExt = file_ext(proposalsDf[code,"docxpath"])
  if( docFileExt != "docx" ) {
    log_error(code,linenum=0,action="OPEN_DOCX",actionOrder=actionOrder,
              rank="",taxon="",
              levelStr = "WARNING",
              errorCode = "DOCX_BAD_EXT",
              errorStr = paste0("File type .",docFileExt," not (yet) supported. Only .docx"),
              notes=paste0("DOCX_FILENAME=",proposalsDf[code,"docx"])
    )
    return(list(metaDf=metaDf))
  } 

  # read DOCX file, no column names
  txt = read_docx(proposalsDf[code,"docxpath"])
  
  
  # get title
  titlePattern="^Short title:"
  titleIdx = grep(titlePattern, txt, ignore.case = T)
  if( length(titleIdx) > 0 ) {
    metaDf$title=
      gsub(pattern="\t+"," ",
           gsub(pattern = "\n+","; ",
                strsplit(txt[titleIdx],split="[sS]hort title *: *")[[1]][2]
           )
      )
  } else {
    log_error(code,linenum=0,action="OPEN_DOCX",actionOrder=actionOrder, 
              rank="",taxon="",
              levelStr  = "WARNING",
              errorCode = "DOCX_TITLE_MISSING",
              errorStr  =  "Title not found",
              notes     = paste0("Line containing expression '",titlePattern,"' not found")
    )
  }
  
  # authors
  authorsPattern="^Author.* and email address.*"
  authorsIdx=grep(pattern=authorsPattern, txt, ignore.case = T)
  if(length(authorsIdx) > 0 ) {
    metaDf$authorsEmails = gsub(pattern="\t+"," ",
                                gsub(pattern = "\n+","; ",
                                     txt[authorsIdx+1]
                                )
    )
  } else {
    log_error(code,linenum=0,action="OPEN_DOCX",actionOrder=actionOrder, 
              rank="",taxon="",
              levelStr  = "WARNING",
              errorCode = "DOCX_AUTHORS_MISSING",
              errorStr  =  "Authors list not found",
              notes     = paste0("Line containing expression '",authorsPattern,"' not found")
    )
  }
  
  # corresponding author
  corrAuthorPattern="^Corresponding author"
  corrAuthIdx=grep(pattern=corrAuthorPattern, txt, ignore.case = T)
  if(length(corrAuthIdx)>0) {
    metaDf$correspondingAuthor = gsub(pattern="\t+"," ",
                                      gsub(pattern = "\n+","; ",txt[corrAuthIdx+1]
                                      )
    )
  } else {
    log_error(code,linenum=0,action="OPEN_DOCX",actionOrder=actionOrder, 
              rank="",taxon="",
              levelStr  = "WARNING",
              errorCode = "DOCX_CORR_AUTHOR_MISSING",
              errorStr  =  "Corresponding Author not found",
              notes     = paste0("Line containing expression '",corrAuthorPattern,"' not found")
    )
   }
  
  # abstract
  abstractPattern="^Abstract"
  abstractIdx=grep(pattern=abstractPattern, txt, ignore.case = T)
  if(length(abstractIdx) >0) {
    metaDf$abstract = gsub(pattern="\t+"," ",
                           gsub(pattern = "\n+","; ",txt[abstractIdx+1]
                           )
    )
  } else {
    log_error(code,linenum=0,action="OPEN_DOCX",actionOrder=actionOrder, 
              rank="",taxon="",
              levelStr  = "WARNING",
              errorCode = "DOCX_ABSTRACT_MISSING",
              errorStr  =  "Abstract not found",
              notes     = paste0("Line containing expression '",abstractPattern,"' not found")
    )
  }

  # return value
  return(list(metaDf=metaDf))
}


# id="2022.003S"; proposalDf=load_proposal(id) # V1 template
# id="2022.002S"; proposalDf=load_proposal(id) # V2 template
# id="2022.015M"; proposalDf=load_proposal(id) # line 5 hostSource invalid chars? 
# code="2022.005D"; proposalDf=load_proposal(code); code; # line 5 genomeComposition missing spaces! 
# code="2022.003M"; proposalDf=load_proposal(code); code; # line 5 hostSource specify="river sediment"

qc_proposal = function(code, proposalDf) {
  # set up error tracking
  errorDf = allErrorDf[FALSE,]
  templateVersion = "error"
  xlxs_colnames = toupper(c(letters,paste0("a",letters))) 
  

  #### guess template version ####

  
  # check row 3, cell 1 for 2023 and later version numbers
  if(!is.null(proposalDf[2,1]) && (
    (substring(proposalDf[2,1],1,13) == "version 2023.")
    || (substring(proposalDf[2,1],1,13) == "version 2024.")
    )) {
    
    # "YYYY.##."
    templateVersionX = substring(proposalDf[2,1],9) 
    # "YYYY."
    templateVersion  = substring(templateVersionX,1,5)

    if( templateVersion == "2024." ) {
      # 2024.1 is the same as 2023.# series
      if( templateVersionX == "2024.1") { templateVersion = "2023." } # QQQQ debug fix after test
      else {
        # new, unspported version
        proposalsDf[code,"templateVersion"]="error"   
        log_error(code,linenum=2,action="OPEN_XLSX",actionOrder=actionOrder, 
                  rank="",taxon="",
                  levelStr="ERROR",errorCode="XLSX.TEMPLATE_UNK",errorStr="XLSX template version",
                  notes=paste0("Cell A2='",proposalDf[2,1],"'")
        )
        return(list())
      }
    }
      
    # 
    # process 2023-style template layout
    #
    
    # complain about any formatting changes
    row3match = (
      tolower(gsub("[^[:alnum:]]","",proposalDf[3, seq(from = 1, to = min(length(xlsx_2023_row3), ncol(proposalDf)))]))
      == 
      tolower(gsub("[^[:alnum:]]","",xlsx_2023_row3))
    )
    #row3match[is.na(row3match)]=FALSE
    row3mismatchCt = sum(!row3match, na.rm = T)
    row3Error= paste("row3 mismatches template at: ", 
                              paste(paste0("column ",
                                           toupper(c(letters,paste0("a",letters)))[which(!row3match)],
                                           "='",
                                           proposalDf[2,which(!row3match)],
                                           "' instead of '",
                                           xlsx_2023_row3[which(!row3match)],"'"),
                                    collapse="; ")
    )
                              
    row4match = (
      tolower(gsub("[^[:alnum:]]","",proposalDf[4, seq(from = 1, to = min(length(xlsx_2023_row4), ncol(proposalDf)))]))
      == 
      tolower(gsub("[^[:alnum:]]","",xlsx_2023_row4))
    )
    row4match[is.na(row4match)]=FALSE
    row4mismatchCt = sum(!row4match, na.rm = T)
    row4Error= paste("row4 mismatches template at: ", 
                     paste(paste0("column ",
                                  toupper(c(letters,paste0("a",letters)))[which(!row4match)],
                                  "='",
                                  proposalDf[2,which(!row4match)],
                                  "' instead of '",
                                  xlsx_2023_row4[which(!row4match)],"'"),
                           collapse="; ")
    )
    # report error if header lines have un-expected cell values
    if( row3mismatchCt > 0 || row4mismatchCt > 0 ) { 
      proposalsDf[code,"templateVersion"]="error"   
      if( row3mismatchCt > 0 ) {
        log_error(code,linenum=2,action="OPEN_XLSX",actionOrder=actionOrder, 
                  rank="",taxon="",
                  levelStr="ERROR",errorCode="XLSX.TEMPLATE_LINE3",errorStr="XLSX template version",
                  notes=row3Error
        )
      }
      if( row4mismatchCt > 0 ) {
        log_error(code,linenum=2,action="OPEN_XLSX",actionOrder=actionOrder, 
                  rank="",taxon="",
                  levelStr="ERROR",errorCode="XLSX.TEMPLATE_LINE4",errorStr="XLSX template version",
                  notes=row4Error
        )
      }
      # stop processing this file
      return(list())
    } # line3 or line4 error
  } else {
    #
    # figure out v1 or v2 templates from contents of rows 2 & 3
    # 
    
    # check row 2 to validate version of templates
    row2v1match = proposalDf[2, seq(from = 1, to = min(length(xlsx_v1_row2), ncol(proposalDf)))] == xlsx_v1_row2
    row2v1matchCt = sum(row2v1match, na.rm = T)
    row2v2match = proposalDf[2, seq(from = 2, to = min(length(xlsx_v2_row2) +
                                                         1, ncol(proposalDf)))] == xlsx_v2_row2
    row2v2matchCt = sum(row2v2match, na.rm = T)
    
    templateLine2Error = ""
    if (row2v1matchCt == sum(!is.na(xlsx_v1_row2))) {
      templateLine2Version = "v1"
    } else if (row2v2matchCt == sum(!is.na(xlsx_v2_row2))) {
      templateLine2Version = "v2"
    } else {
      templateLine2Version = "unrecognized"
      # report the unexpected values
      if( row2v1matchCt > row2v2matchCt ) { 
        templateLine2Error= paste("similar to v1, but with ", 
                                  paste(paste0("column ",toupper(c(letters,paste0("a",letters)))[which(!row2v1match)],"='",proposalDf[2,which(!row2v1match)],"' instead of '",xlsx_v1_row2[which(!row2v1match)],"'"),collapse="; ")
        )
      } 
      if( row2v1matchCt < row2v2matchCt ) { 
        templateLine2Error= paste("similar to v2, but with ", 
                                  paste(paste0("column ",toupper(c(letters,paste0("a",letters)))[which(!row2v2match)],"='",proposalDf[2,which(!row2v2match)],"' instead of '",xlsx_v1_row2[which(!row2v2match)],"'"),collapse="; ")
        )
      } 
    }
    if(params$verbose>1) { cat("     ", code, " XLSX.Row 2: ",templateLine2Version," ",templateLine2Error,"\n")}
    
    # check row 3 to validate version of templates
    row3v1match = proposalDf[3, seq(from = 1, to = min(length(xlsx_v1_row3), ncol(proposalDf)))] == xlsx_v1_row3
    row3v1matchCt = sum(row3v1match, na.rm=T)
    row3v2match = proposalDf[3, seq(from = 1, to = min(length(xlsx_v2_row3), ncol(proposalDf)))] == xlsx_v2_row3
    row3v2matchCt = sum(row3v2match, na.rm=T)
    templateLine3Error = ""
    if (sum(row3v1match,na.rm=T) == sum(!is.na(xlsx_v1_row3))) {
      templateLine3Version = "v1"
    } else if (sum(row3v2match, na.rm=T) == sum(!is.na(xlsx_v2_row3))) {
      templateLine3Version = "v2"
    } else {
      templateLine3Version = "unrecognized"
      # report the unexpected values
      if( row3v1matchCt > row3v2matchCt ) { 
        templateLine3Error= paste("similar to v1, but with ", 
                                  paste(paste0("column ",toupper(c(letters,paste0("a",letters)))[which(!row3v1match)],"='",proposalDf[3,which(!row3v1match)],"' instead of '",xlsx_v1_row3[which(!row3v1match)],"'"),collapse="; ")
        )
      } 
      if( row3v1matchCt < row3v2matchCt ) { 
        templateLine3Error= paste("similar to v2, but with ", 
                                  paste(paste0("column ",toupper(c(letters,paste0("a",letters)))[which(!row3v2match)],"='",proposalDf[3,which(!row3v2match)],"' instead of '",xlsx_v2_row3[which(!row3v2match)],"'"),collapse="; ")
        )
      } 
    }
    if(params$verbose>1){cat("     ", code, " XLSX.Row 3: ",templateLine3Version," ",templateLine3Error,"\n")}
    
    # check for row2/row3 mismatch or errors
    if (templateLine2Version != templateLine2Version || "unrecognized" %in% c(templateLine2Version,templateLine3Version)) {
      
      proposalsDf[code,"templateVersion"]="error"   
      log_error(code,linenum=3,action="OPEN_XLSX",actionOrder=actionOrder, 
                rank="",taxon="",
                levelStr="ERROR",errorCode="XLSX.TEMPLATE_UNK",errorStr="XLSX template version",
                notes=paste0("ROW2 is ", templateLine2Version, " (",templateLine2Error,")",
                             ", ROW3 is ", templateLine3Version," (",templateLine3Error,")")
      )
      return(list())
    } else {
      templateVersion = templateLine2Version
    }
    # finish up templateVersion
    proposalsDf[code,"templateVersion"]=templateVersion
    if( templateVersion == "v1") {
      # WARNING for outdated templates
      errorDf=log_error(code,linenum=3,action="OPEN_XLSX",actionOrder=actionOrder, 
                        rank="",taxon="",
                        levelStr="INFO",errorCode = "XLSX.OLD_TEMPLATE_V1",errorStr = "XLSX template version",
                        notes=paste0("You are using version ",templateVersion,". Please get the latest version from ",params$templateURL)
      )
    }
    #
    # QC Proposal ID
    #
    codeValue="missing" 
    codeCell="undefined"
    codeRow = NA
    if( templateVersion == "v1" ) { codeValue= proposalDf[1,1]; codeCell="A1"; codeRow=1 }
    if( templateVersion == "v2" ) { codeValue= proposalDf[3,1]; codeCell="A3"; codeRow=3 }
    if( templateVersion == "2023." ) { codeValue= proposalDf[1,5]; codeCell="E1"; codeRow=1 }
    if( templateVersion == "2024." ) { codeValue= proposalDf[1,5]; codeCell="E1"; codeRow=1 }
    if( codeValue != code ) {
      if( str_starts(codeValue,"Code") ) {
        if(params$show.xlsx.code_miss) {
          log_error(code,linenum=codeRow,action="OPEN_XLSX",actionOrder=actionOrder, 
                    rank="",taxon="",
                    levelStr="INFO","XLSX.CODE_MISS", "XLSX code missing", 
                    paste0("XLSX cell ",codeCell,
                           " is ", "'",codeValue, 
                           "'; replace with the actual code: '",code,"'")
          )
        }
      } else {
        log_error(code,linenum=codeRow,action="OPEN_XLSX",actionOrder=actionOrder, 
                  rank="",taxon="",
                  levelStr="WARNING", errorCode="XLSX.CODE_BAD",errorStr="XLSX code wrong", 
                  notes=paste0("XLSX cell ",codeCell,
                               " does not match proposal code from filename: ", 
                               "'",codeValue, "' should be '", code,"' ")
        )
      }
    }
  }
  
  if(params$verbose){cat("     ", code, " XLSX template ",templateVersion,"\n")}
    
  
  #
  # extract & standardize
  # cols&rows with change data 
  #
  changeDf = data.frame()

  # map columns
  firstDataRow=4
  if(templateVersion=="v1") { firstDataRow=4; changeDf = proposalDf[firstDataRow:nrow(proposalDf),xlsx_v1_change_cols] }
  else if(templateVersion=="v2") { firstDataRow=4; changeDf = proposalDf[firstDataRow:nrow(proposalDf),xlsx_v2_change_cols] }
  else if(templateVersion=="2023." ) {firstDataRow=5; changeDf = proposalDf[firstDataRow:nrow(proposalDf),xlsx_2023_change_cols]}
  else if(templateVersion=="2024." ) {firstDataRow=5; changeDf = proposalDf[firstDataRow:nrow(proposalDf),xlsx_2023_change_cols]}
  else {
    log_error(code,linenum=firstDataRow,action="OPEN_XLSX",actionOrder=actionOrder, 
              rank="",taxon="",
              levelStr="ERROR", errorCode="XLSX.EMPTY",errorStr="XLSX no change rows found")
    return(list())
    
    }
  colnames(changeDf) = xlsx_change_colnames
  # a flag to exclude rows with irrecoverable errors
  changeDf[,".noErrors"] = TRUE
  changeDf[,".errors"] = NA_character_
  changeDf[,".code"] = code
  changeDf[,".linenum"] = rownames(changeDf)
  
  # find non-empty rows
  hasData = apply(changeDf[,xlsx_change_srcDest_colnames],1,function(row){return(sum(is.na(row))!=length(row))})
  changeDf=changeDf[hasData,]
  if( nrow(changeDf) == 0 ) {
    log_error(code,linenum=firstDataRow,action="OPEN_XLSX",actionOrder=actionOrder, 
              rank="",taxon="",
              levelStr="ERROR", errorCode="XLSX.EMPTY",errorStr="XLSX no change rows found")
    return(list())
  } else {
    proposalsDf[code,"nChanges"] = nrow(changeDf)
  }

  #
  #### extract src/dest taxon names for error reporting ####
  #
  # must be before QC spaces/quotes, etc, as error reporting there
  # uses .changeTaxon!
  changeDf$.srcTaxon =     apply(changeDf[,xlsx_change_srcCols], 1,taxon_get_name
                                 )
  changeDf$.srcRank =      apply(changeDf[,xlsx_change_srcCols], 1,taxon_get_rank)
  changeDf$.srcLineage =   apply(changeDf[,xlsx_change_srcCols], 1,taxon_get_lineage)

  changeDf$.destTaxon =         apply(changeDf[,xlsx_change_destCols],1,taxon_get_name
                                     )
  changeDf$.destRank =          apply(changeDf[,xlsx_change_destCols],1,taxon_get_rank)
  changeDf$.destLineage =       apply(changeDf[,xlsx_change_destCols], 1,taxon_get_lineage)
  changeDf$.destParentName =    apply(changeDf[,xlsx_change_destCols], 1,taxon_get_parent_name)
  changeDf$.destParentRank =    apply(changeDf[,xlsx_change_destCols], 1,function(x) {names(taxon_get_parent_name(x)[1])})
  changeDf$.destParentLineage = apply(changeDf[,xlsx_change_destCols], 1,taxon_get_parent_lineage)
  
  # active taxon: src, or, if missing, dest
  changeDf$.changeTaxon = with(changeDf,ifelse(!is.na(.srcTaxon),.srcTaxon,.destTaxon))
   
  #
  #### QC spaces, quotes, etc ####
  # 
  # TODO: remove leading/trailing spaces, quotes
  dataColumns = grep(names(changeDf),pattern='^\\.', value=T, invert=T)
  for( col in dataColumns ) {
    # col = "srcSpecies" # debug
    #  col = "hostSource" # debug
    
    #
    # silent fixes
    # 

    # non-breaking spaces
    # Mac: AscToChar(202)
    pattern="0xCA"; pat_warn="non-breaking space character"; pat_replace=" "
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(params$tmi && length(qc.matches)>0) { 
      cat("TMI:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") 
    }
    changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    # linux & mac
    pattern="\u00A0"; pat_warn="non-breaking space character"; pat_replace=" "
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(params$tmi && length(qc.matches)>0) { 
      cat("TMI:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") 
    }
    changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    
    # long-dashes:  "–" and "—"
    #pattern=paste0("([\u2012\u2013\u2014\u2015",AscToChar(207),AscToChar(208),"]+)"); pat_warn="long-dash[en/em dashes]";pat_replace="-"
    pattern=paste0("([\u2012\u2013\u2014\u2015]+)"); pat_warn="long-dash[en/em dashes]";pat_replace="-"
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(params$tmi && length(qc.matches)>0) { 
      cat("TMI:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") 
    }
    changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    
    # curvy quotes
    # Mac: AscToChar(210)AscToChar(211)
    # Linux: AscToChar(c(226, 128, 156)) & AscToChar(c(226, 128, 157))
    # AscToChar(211) and SacToChar(210) caused problem son mac, but only when run from the command-line, not inside RStudio!
    #pattern=paste0("([“”",AscToChar(210),AscToChar(211),AscToChar(c(226, 128, 156)),AscToChar(c(226, 128, 157)),"]+)"); pat_warn="curvy quotes";pat_replace='"'
    #pattern=paste0("([“”",AscToChar(c(226, 128, 156)),AscToChar(c(226, 128, 157)),"]+)"); pat_warn="curvy quotes";pat_replace='"'
    pattern=paste0("([\u201C\u201D\u201E\u201F]+)"); pat_warn="curvy quotes";pat_replace='"'
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(params$tmi && length(qc.matches)>0) { 
      cat("TMI:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") 
    }
    changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])

    #
    # curvy single-quotes
    # AscToChar(212)AscToChar(213)
    #pattern=paste0("([‘’]+)"); pat_warn="curvy single quotes";pat_replace="'"
    pattern=paste0("([\u2018\u2019\u201A\u201B]+)"); pat_warn="curvy single quotes";pat_replace="'"	
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(params$tmi && length(qc.matches)>0) { 
      cat("TMI:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") 
    }
    changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])

    
    # newline
    pattern="([\r\n]+)"; pat_warn="newline character(s)";pat_replace=";"
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(params$tmi && length(qc.matches)>0) { 
      cat("TMI:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") 
    }
    changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    
    # leading white space
    pattern="^([ \t]+)"; pat_warn="leading whitespace"; pat_replace=""
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(params$tmi && length(qc.matches)>0) { 
      cat("TMI:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") 
    }
    changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    
    # trailing white space
    pattern="([ \t]+)$"; pat_warn="trailing whitespace"; pat_replace=""
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(params$tmi && length(qc.matches)>0) { 
      cat("TMI:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") 
    }
    changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    
    #
    # repeated spaces
    #
    pattern='( [ ]+)'; pat_warn="repeated_spaces"; pat_replace=" "
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(params$tmi && length(qc.matches)>0) { 
      cat("TMI:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") 
    }
    changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    
    #
    # quotes
    #
    if( col != "comments" ) {
      pattern='(["]+)'; pat_warn="quote"; pat_replace=""
      qc.matches =grep(changeDf[,col],pattern=pattern)
      if(length(qc.matches)>0) { 
        if(params$verbose) { cat("INFO:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") }
        log_error(code,linenum=rownames(changeDf)[qc.matches],
                  action=changeDf$change[qc.matches],actionOrder=actionOrder, 
                  rank=changeDf$rank[qc.matches],taxon=changeDf$.changeTaxon[qc.matches],
                  levelStr="INFO",errorCode="XLSX.QUOTES_REMOVED", errorStr=paste("XLSX has",pat_warn),
                  notes=paste0(paste(col,gsub(pattern,"[\\1]",changeDf[qc.matches,col]),sep=":")," (replacing with '",pat_replace,"')")
        )
        # backup original values
        changeDf[qc.matches,paste0(col,"_orig")] == changeDf[qc.matches,col]
        # remove non-ascii chars
        changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
      } 
    }
  }
   
  #
  # check regex's for each column
  for(i in rownames(value_validation) ) {
    # 
    # more rigorous for some columns
    #
    #pattern=xlsx_col_info[col,"pattern"]; pat_warn=xlsx_col_info[col,"pat_warn"]
    col = value_validation[i,]$col
    error = value_validation[i,]$code
    if(  params$tmi ) {
        cat(with(value_validation[i,],paste0("\tTMI: column_validation(",col,",",code,"):s/",regex,"/",replace,"/\n")))
    }
    if( value_validation[i,]$type == "replace") {
      # 
      # check for presence of regex's to be replaced
      #
      qc.matches =grep(changeDf[,col],pattern=value_validation[i,]$regex)
      if(length(qc.matches)>0) { 
        if( params$verbose) { cat(paste0(value_validation[i,]$class,":",error," has ",length(qc.matches)," cells with ",value_validation[i,]$warn," in column '",col,"'\n")) }
        charEncodings = ""
        if(params$tmi) { 
          # show code points in error message, mostly for I18N problems
          charEncodings = paste0("; charEncodings",
                                 paste(unlist(strsplit(changeDf[qc.matches,col],"")),
                                       CharToAsc(unlist(strsplit(changeDf[qc.matches,col],""))), 
                                       sep="=", collapse=","))
        }
        log_error(code,linenum=rownames(changeDf)[qc.matches],
                  action=changeDf$change[qc.matches],actionOrder=actionOrder, 
                  rank=changeDf$rank[qc.matches], taxon=changeDf$.changeTaxon[qc.matches],
                  levelStr=value_validation[i,]$class,errorCode=error, errorStr=paste("XLSX has",value_validation[i,]$warn),
                  notes=paste(col,":",gsub(value_validation[i,]$regex,"[\\1]",changeDf[qc.matches,col]),
                              charEncodings)
        )
        # backup original values
        changeDf[qc.matches,paste0(col,"_orig")] == changeDf[qc.matches,col]
        # remove non-ascii chars
        changeDf[,col] = gsub(value_validation[i,]$regex,value_validation[i,]$replace,changeDf[,col])
      }
    } else if( value_validation[i,]$type == "required" ) {
      #
      # check required regex matches
      #
      # find values failing regex, and also not NA
      qc.matches = 
        !grepl(changeDf[,col],pattern=value_validation[i,]$regex) &
        !is.na(changeDf[,col])
      
      if(sum(qc.matches)>0) { 
        if(params$verbose) { cat(value_validation[i,]$class,":",error,"has",sum(qc.matches),"cells with",value_validation[i,]$warn,"in column",col,
                                 " in ", code, " on line(s) ", paste(rownames(changeDf)[qc.matches],collapse=","), "\n") 
        }
        # .action is not yet defined, because "change" column has not been QC'ed. Use original "change" text from xlsx
        log_error(code,linenum=rownames(changeDf)[qc.matches],
                  action=changeDf$change[qc.matches],actionOrder=actionOrder, 
                  rank=changeDf$rank[qc.matches], taxon=changeDf$.changeTaxon[qc.matches],
                  levelStr=value_validation[i,]$class,errorCode=error,errorStr=value_validation[i,]$warn,
                  notes=paste(col,gsub(value_validation[i,]$regex,"[\\1]",changeDf[qc.matches,col]),sep=":"))
        # backup original values
        #changeDf[qc.matches,paste0(col,"_orig")] == changeDf[qc.matches,col]
        # remove non-ascii chars
        #changeDf[,col] = gsub(value_validation[i,]$regex,value_validation[i,]$replace,changeDf[,col])
      }
    }  # for col
  } # for regex 
  
  
  #
  #### QC controlled vocabularies ####
  #

  # check that all columns are present
  missingCvNames = grep(
    # cvList names not in changeDf columns
    x=names(cvList)[!names(cvList) %in% names(changeDf)],
    # except those with "."
    pattern="\\.",value=T,invert=T)
                       
  if(length(missingCvNames)>0) {
    if(params$verbose) { cat("ERROR:",code,"missing CV columns [",paste(missingCvNames, collapse=","),"] on row 3\n") }
    log_error(code,linenum=3,action=NA, actionOrder=actionOrder, 
              rank=NA, taxon=NA,
              levelStr="ERROR",errorCode="XLSX.MISSING_COLUMN",errorStr=paste("XLSX missing column"),
              notes=paste("XLSX missing header [",missingCvNames,"]")) 
    return(list())
  }
  
  #
  # convert "please select" to NA
  #
  isPleaseSelect = regexpr(text=changeDf, pattern=".*please.*select.*", ignore.case = T)>0
  isPleaseSelect[is.na(isPleaseSelect)]=FALSE
  if(sum(isPleaseSelect,na.rm=T)>0) {
    # replace "[Please Select]" with NA
    changeDf[isPleaseSelect]=NA
  }
  
  #
  # check that all the terms are legal
  #
  # remove lines containing invalid terms
  #
  correctionsDf = changeDf[FALSE,]
  for(cv in grep(x=names(cvList),pattern="\\.",value=T,invert=T) ) {
    # cv="molecule" # debug
    
    # compare to CV terms, removing all spaces and capitalization
    isTermPerfect = changeDf[,cv] %in% cvList[[cv]]
    isTermClose = 
      tolower(gsub(pattern="[^[:alpha:];+-]",replacement="",changeDf[,cv])) %in% 
      tolower(gsub(pattern="[^[:alpha:];+-]",replacement="",cvList[[cv]]))

    # find the term we were close to
    for( row in rownames(changeDf)[!isTermPerfect & isTermClose] ) {
       correctedTermIdx = which(
        tolower(gsub(pattern="[^[:alpha:];+-]",replacement="",changeDf[row,cv]))
        ==
        tolower(gsub(pattern="[^[:alpha:];+-]",replacement="",cvList[[cv]]))
      )
      # check for multiple close matches
      if(length(correctedTermIdx) > 1) {
        # matching multiple things solves nothin
        isTermClose[which(row==rownames(changeDf))]=FALSE
      } else {
        # one, and only one approximate match: then use it!
       correctionsDf[row,cv] = cvList[[cv]][correctedTermIdx]
      }
    }
    # 
    # reject unknown terms
    #
    badTerms=changeDf[!isTermClose,]
    if(nrow(badTerms)>0) {
      # complain
      if(params$verbose) { cat("ERROR:",code,"illegal term in CV column '",cv,"':",
                               "[",paste(paste0(rownames(badTerms),":",badTerms[,cv]), collapse=","),"] on rows (",paste(rownames(badTerms),collapse=","),")\n") }
      # for optional CVs, NA the term, for obligatory ones, flag the line
      changeDf[rownames(badTerms),".errors"] = paste("CV '",cv,"' incorrect value [",badTerms[,cv],"]. Valid terms: [",paste0(cvList[[cv]],collapse=","),"]")
      levelStr="UNKNOWN"
      if(cv %in% c("genomeCoverage","molecule","hostSource") ) {
        # remove value, issue WARNING
        changeDf[rownames(badTerms),cv] = NA
        levelStr="WARNING"
      } else {
        # skip line, issue ERROR
        changeDf[rownames(badTerms),".noErrors"] = FALSE
        levelStr="ERROR"
      }
      # issue error/warning
      log_error(code,linenum=rownames(badTerms),action=badTerms$change,actionOrder=actionOrder,
                rank=badTerms$rank, taxon=badTerms$.changeTaxon,
                levelStr=levelStr,errorCode="XLSX.INVALID_TERM",errorStr=paste("XLSX incorrect term in column",cv),
                notes=paste("XLSX incorrect value [",badTerms[,cv],"]. Valid terms: [",paste0(cvList[[cv]],collapse=","),"]")) 
      if(levelStr == "ERROR") {
        next # skip this line
      }
    } # if badTerms
    #
    # whine about terms with spacing/caps problems
    #
    flawedTerms=changeDf[isTermClose & !isTermPerfect,] #cv, drop=FALSE]
    if(nrow(flawedTerms)>0) {
      # complain, never mind, just auto-correct
      if( params$tmi ) {
        if(params$verbose || params$tmi ) { 
          cat("TMI:",code,"typo in CV column '",cv,"':",
              "[",paste(paste0(rownames(flawedTerms),":",flawedTerms[,cv]), collapse=","),"] on rows (",paste(rownames(flawedTerms),collapse=","),")\n") 
        }
        if( params$tmi ) { 
          # only warn about terms with spacing issues if we're in TMI mode
          log_error(code,linenum=rownames(flawedTerms),action=flawedTerms$change, actionOrder = actionOrder,
                    rank=flawedTerms$rank, taxon=flawedTerms$.changeTaxon,
                    levelStr="TMI",errorCode="XLSX.TYPO_TERM",errorStr=paste("fixed term with typo (space,caps) in column ",cv),
                    notes=paste0("Term '",flawedTerms[,cv],"' replaced with '",correctionsDf[rownames(flawedTerms),cv],"'. Valid terms: [",paste0(cvList[[cv]],collapse=","),"]")) 
        }
      }
      # correct
      changeDf[rownames(flawedTerms),cv] = correctionsDf[rownames(flawedTerms),cv]
    } # if flawedTerms
  } # for cvList
  
  # clean change -> action
  changeDf$.action = cvList[[".change2action"]][tolower(changeDf$change)]
  badRows = changeDf[is.na(changeDf$.action),]
  if( nrow(badRows) > 0 ) { 
    # mark bad rows
    changeDf[rownames(badRows),".noErrors"] = FALSE
    log_error(code,linenum=rownames(badRows),action=badRows$change,actionOrder=actionOrder, 
              rank=badRows$rank,taxon=badRows$.changeTaxon,
              levelStr="ERROR",errorCode="ACTION.UNK",errorStr=paste("XLSX incorrect term in column","Change"),
              notes=paste("XLSX incorrect value [",badRows[,"change"],"]. Valid terms: [",paste0(names(cvList[[".change2action"]]),collapse=","),"]")
    ) 
  }
  
  # return warnings, if any
  return(list(changeDf=changeDf))
} # qc_proposal()
# ```
# 
# 
# # Load and QC
# 
# ```{r load_and_qc}

##### LOAD xlsx/DOCX, QC ##### 


load_and_qc_proposals = function(proposalsDf,changeList) {
  #
  # this uses the PROPOSAL.XLSX schema (naming convention)
  #
  
  #
  # extracted changes
  #
  
  # make this re-entrant - only create/delete lists if they don't exist
  if(!exists("docxList"))   { docxList   = list() }
  if(!exists("xlsxList"))   { xlsxList   = list() }
  if(!exists("changeList")) { changeList = list() }
  codes = rownames(proposalsDf) # all proposals
  #codes = c("2022.003S","2022.002S")  # templates: V1, V2
  #codes = c(codes,"2022.007F", "2022.006P") # missing xlsx, typo in xlsx
  # codes = c("2022.085B") # MOVE
  for( code in codes ) {
    # code = codes[1]; code
    # code = codes[which(code==codes)+1]; code; # next
    # code = "2022.003S" # V1
    # code = "2022.002S" # V2
    #  code =  "2022.016P" # space in docx filename
    #
    # load
    #
    
    if(!is.null(changeList[[code]])) {
      # already loaded
    } else {
      # LOAD
      xlsx_fname=proposalsDf[code,"xlsx"]
      if(is.na(xlsx_fname)) {
        cat("# SKIPPED: ",code," (no .xlsx)\n")
      } else {
        #
        # try to load proposal metdata from DOCX
        #
        if(is.null(docxList[[code]])){
          # load raw docx file
          # results=list(metaDf,errorDf)
          results= load_proposal_docx(code)
          docxList[[code]] = results[["metaDf"]] 
          errorDf = .GlobalEnv$allErrorDf %>% filter(code == code)
          cat("# LOADED: ",code," DOCX with ",nrow(errorDf)," errors/warnings\n")
        } else {
          cat("#  DOCX FROM CACHE: ",code,"\n")
        }
        .GlobalEnv$proposalsDf[code,names(docxList[[code]])] = docxList[[code]]
        #
        # try to load proposal from XLSX
        #
        if(is.null(xlsxList[[code]])){
          # load raw xlsx file
          results =  load_proposal(code)
          changeDf = results[["proposalDf"]]
          errorDf  =  .GlobalEnv$allErrorDf %>% filter(code == code)
          cat("# LOADED: ",code," XLS with ",nrow(errorDf)," errors/warnings\n")
          xlsxList[[code]] = results[["proposalDf"]]
          
          # if load had errors/warnings, update error file
          if( nrow(errorDf)>0  ){
            write_error_summary(.GlobalEnv$allErrorDf)          
          }
          cat("# LOADED: ",code,"\n")
          
        } else {
          cat("# XLSX FROM CACHE: ",code,"\n")
        }
        proposalDf = xlsxList[[code]]
        
        #
        # qc
        #
        if( is.null(proposalDf) ) {
          # load failed, move on
          cat("SKIP: ", code, ": proposal could not be loaded\n")
        } else {
          # successful load
          cat(paste0("# QC start: ", code, ": proposal loaded",
                     " (",which(codes==code)," out of ",length(codes),")\n"))
          
          results = qc_proposal(code,proposalDf)
          errorDf = .GlobalEnv$allErrorDf %>% filter(code == code)
          cat("# QCed:     ",code," with ",nrow(errorDf)," errors/warnings\n")
          changeDf = results[["changeDf"]]
          if(!is.null(changeDf)) {
            changeList[[code]]=changeDf
          }
        } # loadable 
      } # LOAD if exists(xlsx)
    } # LOAD if not in changeList
  }
  # load summary
  cat("changeList: ",paste(names(changeList)),"\n")
  
  # error summary
  write_error_summary(.GlobalEnv$allErrorDf)
  #kable(allErrorDf,caption = paste0("QC: Summary of ERRORs and WARNINGs"))
  
  return(changeList)
} # load_and_qc_proposals()

#####  PROPOSAL CACHE #####

###### LOAD PROPOSAL CACHE #####

# Save loaded proposals, so we can work on applying them
loadCacheFilename = file.path(params$proposals_dir,".RData")

if(!exists("changeList")) { changeList = list() }

if( params$load_proposal_cache && file.exists(loadCacheFilename) ) {
  cat("LOADing *PROPOSALS* from ", loadCacheFilename,"...\n")
  load(file=loadCacheFilename,verbose = T)
  cat("LOAD DONE","\n")
} else { 
  .GlobalEnv$proposalsDf = scan_for_proposals()
  .GlobalEnv$changeList = load_and_qc_proposals(proposalsDf,changeList)
}

###### SAVE PROPOSAL CACHE #####

if( params$save_proposal_cache ) {
  save(file=loadCacheFilename,list=c(
    "proposalsDf",
    "changeList"
  ))
  cat("WROTE loaded and merged proposals to ",loadCacheFilename,"\n")
}


#
#### PRE-PROCESS MERGED ACTIONS ####
#

##### merge proposals ######

#
# iterate through proposalsDf
#   merge into allChangesDf
#
skipped=0
merged=0
codes = rownames(proposalsDf)
allChangeDf = NULL

#codes = c('2022.085B'); code=codes[1]; code; # debug
for( code in codes) {
  # code = rownames(proposalsDf)[1]; code # debug 
  # code = rownames(proposalsDf)[2]; code # debug 
  # code = rownames(proposalsDf)[3]; code # debug 
  # code = rownames(proposalsDf)[4]; code # debug 
  # code = rownames(proposalsDf)[8]; code # debug - 20221013 current crash: around "Suid alphaherpesvirus 1"
  # make sure we were able to load the changesDf
  # code='2022.001A' # create Itzamnaviridae in [tree root]
  # code='2022.085B'; code; # debug - MOVE test
  # code = rownames(proposalsDf)[which(code==rownames(proposalsDf))+1];code; # DEBUG NEXT
  # code = codes[which(codes==code)+1] # DEBUG NEXT
  
  changeDf = changeList[[code]]
  if( is.null(changeDf) ) {
    skipped = skipped+1
    if(params$verbose>1) {cat("SKIP ",code,": no changeDf loaded\n")}
  } else {
    # remove error lines
    cleanChangeDf = changeDf[changeDf$.noErrors == TRUE,]
    #
    # iterate over changes
    #
    cat(paste0("# MERGE PROC: ",code," with ", nrow(cleanChangeDf), " changes ",
               " (",which(codes==code)," out of ",length(codes),")\n"))
    
    # append 
    if( is.null(allChangeDf) ) {
      allChangeDf = cleanChangeDf
    } else {
      # append 
      allChangeDf = as.data.frame(rbindlist(list(allChangeDf, cleanChangeDf), fill=T))
    }
    
    merged=merged+1
  } # if have proposal XLSX
} # for proposal code

# if there are any changes, add index by code:linenum
tot_changes = 0
if(!is.null(allChangeDf)) {
  allChangeDf$.codeLine = paste0(allChangeDf$.code,":",allChangeDf$.linenum)
  rownames(allChangeDf) = allChangeDf$.codeLine 
  tot_changes = nrow(allChangeDf)
}

cat("# DONE. MERGED",merged," proposal(s), skipped",skipped,"; total changes:",tot_changes,"\n")
#


##### sort merged actions #####
#
#
# sort proposal by action/rank
#
# first all "abolish" & "merge" from bottom to top rank
# than all other actions from top rank to bottom
#

# action sort
# -1  abolish, merge
#  1  create,  split
#  2  rename,  move
actionOrderDf = rbind(
  # destructive: run first, species->root by rank
  "abolish"  = c(actionOrder = -1, subRankOrder= 0),
  "merge"    = c(actionOrder = -1, subRankOrder= 0),
  # constructive: root -> species by rank
  # first create (w/in a rank)
  "new"      = c(actionOrder =  1, subRankOrder= 1),
  "promote"  = c(actionOrder =  1, subRankOrder= 1),
  "demote"   = c(actionOrder =  1, subRankOrder= 1), 
  "split"    = c(actionOrder =  1, subRankOrder= 1),
  # then move (w/in a rank)
  "rename"   = c(actionOrder =  1, subRankOrder= 2),
  "move"     = c(actionOrder =  1, subRankOrder= 2)
)

allChangeDf$.actionOrder = actionOrderDf[allChangeDf$.action,"actionOrder"]
allChangeDf$.subRankOrder = actionOrderDf[allChangeDf$.action,"subRankOrder"]

# rank sort - opposite directions depending on actionOrder
rankMap = dbCvList[["rank"]]$id
names(rankMap) = dbCvList[["rank"]]$name
if( !is.null(nrow(allChangeDf)) && nrow(allChangeDf) > 0 ) {
  allChangeDf[,".rankOrder"] = dbCvMapList[["rank"]][allChangeDf$rank] * sign(allChangeDf$.actionOrder)
  
  # tentative sort 
  allChangeDfOrder = NULL
  allChangeDfOrder = order(
    allChangeDf$.actionOrder, 
    allChangeDf$.rankOrder, 
    allChangeDf$.subRankOrder)

  # stash that
  allChangeDf[allChangeDfOrder,".changeOrder"] = seq(1,nrow(allChangeDf))
  
}
##
##### detect chaining #####
##
# we find series of changes that affect the same taxon
# we then order them so that they will get applied correctly. 
# 
# for example, [B->C], A->B must be re-ordered [A->B, B>C]
#
# previously, actions were sorted and .changeOrder was set
# we set .chainOrder, using decimals added to the .changeOrder of the first link
# After we finish, we re-set the .changeOrder variable. 
# original .changeOrder is stored in 
##

if( !is.null(nrow(allChangeDf)) && nrow(allChangeDf) > 0 ) {
  # by default, everything is not chained
  .GlobalEnv$allChangeDf$.chainOrder = 0
  .GlobalEnv$allChangeDf$.chainStart = ""
  .GlobalEnv$allChangeDf$.originalOrder = .GlobalEnv$allChangeDf$.changeOrder
  .GlobalEnv$allChangeDf$.srcRankTaxon  = paste0(allChangeDf$.srcRank,":",allChangeDf$.srcTaxon)
  .GlobalEnv$allChangeDf$.destRankTaxon = paste0(allChangeDf$.destRank,":",allChangeDf$.destTaxon)
  
  # scan for chained changes: prev.destTaxon = next.srcTaxon (and not self)
  # note that moves, merges and splits can all have srcTaxon=destTaxon!
  #
  link_cols= c(".codeLine",".srcRankTaxon",".destRankTaxon","change",".changeOrder",".chainOrder",".chainStart")
  chainedTaxaDf = merge(
    x=(allChangeDf[,link_cols] %>% filter(.destRankTaxon != 'NA:NA')),
    y=(allChangeDf[,link_cols]  %>% filter(.srcRankTaxon != 'NA:NA')),
    by.x = ".destRankTaxon", 
    by.y=".srcRankTaxon", 
    suffixes = c("Prev","Next")
  ) %>% filter(.codeLinePrev != .codeLineNext)
  
  if( nrow(chainedTaxaDf) == 0) {
    # no chains to re-order
    if(params$verbose){cat("# CHAIN: NO chained changes found\n")}
  } else {
    
    if(params$verbose){cat("# CHAIN: Found", nrow(chainedTaxaDf), "links!\n")}

    # index df
    rownames(chainedTaxaDf) = paste0(chainedTaxaDf$.codeLinePrev,":",chainedTaxaDf$.codeLineNext)
    chainTaxaCols = c(".codeLinePrev",".changeOrderPrev",".srcRankTaxon","changePrev",".destRankTaxon","changeNext",".destRankTaxonNext",".codeLineNext",".changeOrderNext")
    # View(chainedTaxaDf[,chainTaxaCols])
    
    if(params$tmi){with(chainedTaxaDf, 
                        cat("",paste0("    # LINK :        ",.codeLinePrev," ",.srcRankTaxon,
                                      " [",changePrev,"] ",
                                      .destRankTaxon," ", .codeLineNext, 
                                      "\n")))}
    chainDepth = 1
    
    #
    # mark first link in each chain
    #
    # links where srcRankTaxon does NOT appear in the DEST of another link
    chainStartsSimple =  !(chainedTaxaDf$.srcRankTaxon %in% chainedTaxaDf$.destRankTaxon)
    # links [A,A] (src=dest), where A is not in the DEST of any other link, except the  link itself
    chainStartsNoOp = chainedTaxaDf$.srcRankTaxon == chainedTaxaDf$.destRankTaxon & !(chainedTaxaDf$.srcRankTaxon %in% chainedTaxaDf$.destRankTaxon[chainedTaxaDf$.srcRankTaxon != chainedTaxaDf$.destRankTaxon])
    # merge start lists
    chainStartLinks = chainStartsSimple | chainStartsNoOp
    
    
    # mark order in  chain list (.chainOrder = [.changeOrder].1)
    chainedTaxaDf[chainStartLinks,".chainOrderPrev"] = paste0(chainedTaxaDf[chainStartLinks,".changeOrderPrev"],".",chainDepth)
    chainedTaxaDf[chainStartLinks,".chainStartPrev"] = chainedTaxaDf[chainStartLinks,".codeLinePrev"]
    
    # isolate the active links we're building chains on
    activeLinksDf = chainedTaxaDf[chainStartLinks,] 

    if(params$verbose){cat("# CHAIN: level ", chainDepth, " has ", length(unique(activeLinksDf$.codeLinePrev)), " startLinks in ",nrow(activeLinksDf), " chains","\n")}
    if(params$tmi){with(activeLinksDf, 
                        cat("",paste0("    #",
                                      " CHAIN ", .chainStartPrev, 
                                      " LINK :        ",.codeLinePrev," ",.srcRankTaxon,
                                      " [",changePrev,"] ",
                                      .destRankTaxon," ", .codeLineNext, 
                                      "\n")))}
    
    # copy .chainOrder to master change list
    # use fact rowname(allChangeDf)=.codeLine
    .GlobalEnv$allChangeDf[activeLinksDf$.codeLinePrev,".chainOrder"] = activeLinksDf$.chainOrderPrev
    .GlobalEnv$allChangeDf[activeLinksDf$.codeLinePrev,".chainStart"] = activeLinksDf$.chainStartPrev
    
    if(params$verbose){with(activeLinksDf,
                            cat("",paste0("     #   SET prevChainOrder:  ",
                                          .codeLinePrev, " = ", .chainOrderPrev, 
                                          " which is", .srcRankTaxon, 
                                          " [", changePrev, "] ", 
                                          .destRankTaxon, 
                                          "\n")
                            ))}

    #
    # get and mark next link in each chain
    #
    while( nrow(activeLinksDf)>0 ) {
      
      if(params$tmi){cat("",paste0("    # CHAIN: level ", chainDepth, " process link.Next\n"))}

      # print active links
      chainDepth = chainDepth +1
   
      # 
      # set depth of 2nd change in the link
      #
      # use the greater of
      #   the chainHead.actionOrder + chainDepth/10
      #   the  nextLink.actionOrder + chainDepth/10
      #
      # this minimizes re-order, at the expense of less of an audit trail (not all orders have the head as prefix)
      chainOrderNext = as.integer(chainedTaxaDf[rownames(activeLinksDf),".chainOrderPrev"])+chainDepth/10
      orderNext =      as.integer(chainedTaxaDf[rownames(activeLinksDf),".changeOrderNext"])+chainDepth/10
      chainedTaxaDf[rownames(activeLinksDf),".chainOrderNext"] = ifelse(chainOrderNext > orderNext, chainOrderNext, orderNext)
 
      # carry over chain name
      chainedTaxaDf[rownames(activeLinksDf),".chainStartNext"] = chainedTaxaDf[rownames(activeLinksDf),".chainStartPrev"]
      
      # copy to global
      .GlobalEnv$allChangeDf[chainedTaxaDf[rownames(activeLinksDf),".codeLineNext"],".chainOrder"] = chainedTaxaDf[rownames(activeLinksDf),".chainOrderNext"] 
      .GlobalEnv$allChangeDf[chainedTaxaDf[rownames(activeLinksDf),".codeLineNext"],".chainStart"] = chainedTaxaDf[rownames(activeLinksDf),".chainStartNext"] 
      
      if(params$tmi){with(chainedTaxaDf[rownames(activeLinksDf),],
                          cat("",paste0("     #   SET nextChainOrder:  ",
                                        .codeLineNext, " = ",.chainOrderNext, 
                                        " which is ", .destRankTaxon, 
                                        " [", changeNext, "] ", 
                                        .destRankTaxonNext, 
                                        "\n")
                          ))}
      
      # copy to "in" side of any following links
      for( activeRow in rownames(activeLinksDf) ) {
        # activeRow = rownames(activeLinksDf)[1] # debug
        activeRowDf = chainedTaxaDf[activeRow,]
        # copy outgoing (Next) .chainOrder to incoming (Prev) side of any following links
        onwardLinks=chainedTaxaDf$.codeLinePrev == activeRowDf$.codeLineNext
        if( sum(onwardLinks) > 0) {
          # copy that chainOrder from the Next of this link, to Prev of next link
          chainedTaxaDf[onwardLinks,".chainOrderPrev"] = activeRowDf$.chainOrderNext
          chainedTaxaDf[onwardLinks,".chainStartPrev"] = activeRowDf$.chainStartNext
          
          if(params$tmi){with(chainedTaxaDf[onwardLinks,],
                              cat("",paste0("     #   SET NEXT.prevChainOrder:  ",
                                            .codeLinePrev, " = ",.chainOrderPrev, 
                                            " which is ", .srcRankTaxon, 
                                            " [", changePrev, "] ", 
                                            .destRankTaxon, 
                                            "\n")
                              ))}
        } # if onward links

      } # for each active link/chain
      
      # find next links in this chain
      activeLinksDf = chainedTaxaDf[
        # links the connect to the Next of activeLinks
        chainedTaxaDf$.codeLinePrev %in% activeLinksDf[,".codeLineNext"] &
          # AND haven't been traversed already
          !(chainedTaxaDf$.chainStartPrev %in% activeLinksDf[,".chainStartNext"])
        ,]
      
      if(params$verbose){cat("# CHAIN: level ", chainDepth, " has ", length(unique(activeLinksDf$.codeLinePrev)), " startLinks in ",nrow(activeLinksDf), " chains","\n")}
      if(params$tmi){with(activeLinksDf, 
                          cat("",paste0("    #",
                                        " CHAIN ", .chainStartPrev, 
                                        " LINK :        ",.codeLinePrev," ",.srcRankTaxon,
                                        " [",changePrev,"] ",
                                        .destRankTaxon," ", .codeLineNext, 
                                        "\n")))}
      
        
    } # while more links to traverse
    #
    # record final ordering
    #
    .GlobalEnv$allChangeDf$.chainedOrder=ifelse(allChangeDf$.chainOrder == "0", allChangeDf$.changeOrder, allChangeDf$.chainOrder)
    allChangeDfOrder = order(as.numeric(.GlobalEnv$allChangeDf$.chainedOrder))
    .GlobalEnv$allChangeDf[allChangeDfOrder,".changeOrder"] = seq(1,nrow(allChangeDf))
    #
    # errors for any links that didn't get traversed
    #
    badLinksDf = chainedTaxaDf %>% filter(.chainOrderPrev==0 | .chainOrderNext==0)
    if( nrow(badLinksDf)) {
      # these changes were in chains, but didn't get traversed. 
      # could be chains w/o heads, etc
      for( badIdx in rownames(badLinksDf)) {
        prevTaxon = .GlobalEnv$allChangeDf[badLinksDf[badIdx,".codeLinePrev"],]
        badLinkDf = badLinksDf[badIdx,]
        
        
        log_error(prevTaxon$.code,linenum=prevTaxon$.linenum,action=prevTaxon$change,actionOrder=prevTaxon$actionOrder, 
                  rank=prevTaxon$rank,taxon=prevTaxon$.changeTaxon,
                  levelStr="ERROR",errorCode="CHAIN.HIDDEN",errorStr="Chained change order not resolved",
                  notes=paste("Change links with ",badLinkDf$.codeLineNext,": ",
                              badLinkDf$.destRankTaxon," [",badLinkDf$changeNext,"] ",badLinkDf$.destRankTaxonNext)
            
        )
      } # for hidden link
      # write any linking errors out
      write_error_summary(.GlobalEnv$allErrorDf)
    } # if any hidden links
  } # if there are CHAINS
} # if there are changes
#### APPLY CHANGES functions ####
#
# This uses changes (in the PROPOSAL.XLSX schema) to make changes (one per line)
# to the curMSL (in the database taxonomy_node table schema)
#

# excel                                # R.changeDf                             # db.load_next_msl          # db.taxonomy_node          # R.newMSL
#  "Exemplar GenBank Accession Number" = "exemplarAccession",                   "exemplarAccession"         "genbank_accession_csv"     "genbank_accession_csv"
#  "Exemplar \r\nvirus name"           = "exemplarName",                        "exemplarName"              ["exemplar_name"]           "exemplar_name"
# "Virus name abbreviation"            = "Abbrev",                              "Abbrev"                    "abbrev_csv"                "abbrev_csv"
# "Exemplar\r\nisolate designation"    = "exemplarIsolate",                     "exemplarIsolate"           "isolate_csv"               "isolate_csv"
# "Genome coverage"                    = "genomeCoverage",                      "isComplete"                ["genome_coverage"]         "genome_coverage"
# "Genome composition"                 = "molecule",                            "molecule"                  "molecule_id"               "molecule_id"
# "Host/Source"                        = "hostSource",                          "hostSource"                ["host_source"]             "host_source"
# "Change"                             = "change",                              "change","_action"          "in_change","out_change"    "in_change", "out_change"(prevMSL)
# "Rank"                               = "rank",                                "rank","_dest_taxon_rank"   "level_id"                  "level_id"
# "Comments"                           = "comments"                             "comments"                  NA                          "comments"

# map XLSX column names to DB column names
xlsx2dbMap = c(
      "exemplarAccession"="genbank_accession_csv"
      # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
      ,"exemplarName" = "exemplar_name"
      ,"Abbrev"="abbrev_csv"
      # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
      ,"exemplarIsolate"="isolate_csv"
      # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
      ,"genomeCoverage"="genome_coverage"
      # genomeComposition = molecule_id 
      ,"molecule" = "molecule_id"
      # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
      ,"hostSource" = "host_source"
      ,"rank" = "level_id"
      ,"comments" = "notes"
  )
# ............................................................................
#
# update_lineage (.GlobalEnv$newMSL; recursive)
# 
# ............................................................................
# INPUTs:
#   taxnode_id = node whose children need updating
#   destLineage = lineage of taxnode_id
#
update_lineage = function(parent_id,parent_lineage, 
                          parent_otherLineage, parent_otherLineageProposal, 
                          parent_otherLineageAction){
  
  # update parent
  parent_row = which(.GlobalEnv$newMSL$taxnode_id == parent_id)
  .GlobalEnv$newMSL[parent_row,".otherLineage"]         = parent_otherLineage # should append? 
  .GlobalEnv$newMSL[parent_row,".otherLineageProposal"] = ifelse(is.na(.GlobalEnv$newMSL[parent_row,".otherLineageProposal"]),parent_otherLineageProposal,
                                                                 paste(.GlobalEnv$newMSL[parent_row,".otherLineageProposal"] ,parent_otherLineageProposal,sep=";"))
  .GlobalEnv$newMSL[parent_row,".otherLineageAction"]   = ifelse(is.na(.GlobalEnv$newMSL[parent_row,".otherLineageAction"]),parent_otherLineageAction,
                                                                 paste(.GlobalEnv$newMSL[parent_row,".otherLineageAction"] ,parent_otherLineageAction,sep=";"))
  
  # find kids, get their taxnode_is
  ct = 0
  kid_rows= which(.GlobalEnv$newMSL$parent_id == parent_id)
  if(length(kid_rows)>0) {
    # update kids' lineage
    .GlobalEnv$newMSL[kid_rows,"lineage"] = paste(parent_lineage,.GlobalEnv$newMSL[kid_rows,]$name,sep=";")
    if(!is.na(parent_otherLineage)) {
      .GlobalEnv$newMSL[kid_rows,".otherLineage"]   = paste(parent_otherLineage,.GlobalEnv$newMSL[kid_rows,]$name,sep=";") # should append?
    }
    
    # recurse
    for(kid_row in kid_rows ) {
      if(params$verbose>1) {
        cat("update_lineage: ",
            .GlobalEnv$newMSL$taxnode_id[kid_row]," ",
            .GlobalEnv$newMSL$lineage[kid_row],
            " (",as.character(.GlobalEnv$newMSL$rank[kid_row]),")",
            ifelse(is.na(.GlobalEnv$newMSL$.otherLineage[kid_row]),"",paste0("[AKA ",.GlobalEnv$newMSL$.otherLineage[kid_row])), 
            "\n")
      }
      ct = ct + update_lineage(.GlobalEnv$newMSL$taxnode_id[kid_row],.GlobalEnv$newMSL$lineage[kid_row],
                               .GlobalEnv$newMSL$.otherLineage[kid_row], parent_otherLineageProposal, 
                               parent_otherLineageAction )
    }
  }
  return(ct)
}

# ............................................................................
#
# apply_changes
#
# ............................................................................
apply_changes = function(changesDf) {
  #
  # iterate over changes
  # 
  errorDf = allErrorDf[FALSE,]
  
  # genera to scan for binomial issues when we're done
  # QQQQ: should replace this with an admin column on newMSL!
  renamedGenera = data.frame(name=c("genusName"),taxnode_id=0, code="2023.000A", linenum="0", actionOrder=0)
  
  #### admin columns ####
  .GlobalEnv$curMSL[,c(".split",".split_kept")] = FALSE
  .GlobalEnv$curMSL[,c(".split_code",".split_linenum",".split_actionOrder",".split_acc_used")] = NA_character_

  .GlobalEnv$newMSL[,c(".split",".split_kept")] = FALSE
  .GlobalEnv$newMSL[,c(".split_code",".split_linenum",".split_actionOrder",".split_acc_used")] = NA_character_

  rowsWoErrors = rownames(changesDf[changesDf$.noErrors,])
  actionOrder = 0
  #### FOR(row) ####
  for(row in rowsWoErrors) {
    # row = rownames(changesDf)[1] ; row; changesDf[row,"change"];# debug
    # row = rownames(changesDf)[63] ; row; changesDf[row,"change"];# debug
    # row = "12"; linenum; changesDf[row,];
    # row = rownames(changesDf)[which(row==rownames(changesDf))+1]; row # debug: next row
    # .........................................................................
    #
    #
    ##### change definition #####
    #
    # .........................................................................
    
    # track order of actions
    actionOrder = actionOrder + 1 

    # get this change's line
    curChangeDf = changesDf[row,]
    action = curChangeDf$change
    code   = curChangeDf$.code
    linenum = curChangeDf$.linenum
    proposalBasename=proposalsDf[code,"basename"]
    proposalZip = paste0(proposalBasename,".zip")
    curChangeDf$.order = actionOrder
    
    #
    # check CVs
    #

    # rank
    isRankOk = tolower(curChangeDf$rank) %in% dbCvList[["rank"]]$name
    if( !isRankOk ) {
        log_change_error(changeDf=change,levelStr="WARNING", errorCode="RANK.UNK",errorStr="Rank name unknown", 
                  notes=paste0(", rank=", curChangeDf$rank, 
                               "; knownRanks=", paste(dbCvList[["rank"]]$name[-1],collapse = ","))
        )
    } 
    else {
      rankClean = tolower(curChangeDf$rank)
    }
    
    # linenum = rownames(changeDf)[1] # debug
    srcRealmSpecies = curChangeDf[xlsx_change_srcCols]
    destRealmSpecies = curChangeDf[xlsx_change_destCols]
    
    # should these be a vector operation, adding columns to changeDF? 
    # if we do that, we need to separate srcTaxon's label (rank) and value (name)
    # QQQQ: use the columns already in changeDf [.srcTaxon,.srcRank, etc]
    srcTaxonName  =taxon_get_name(srcRealmSpecies)
    srcTaxonRank  =taxon_get_rank(srcRealmSpecies)
    srcLineage    =taxon_get_lineage(srcRealmSpecies)
    destTaxonName =taxon_get_name(destRealmSpecies)
    destTaxonRank =taxon_get_rank(destRealmSpecies)
    destParentName=taxon_get_parent_name(destRealmSpecies)
    destLineage   =taxon_get_lineage(destRealmSpecies)
    destParentLineage =taxon_get_lineage(destRealmSpecies,masks=c(destTaxonRank))

    # check where this name may or may not exist already
    srcNewMatches= (.GlobalEnv$newMSL$name==as.character(srcTaxonName))
    destNewMatches=(.GlobalEnv$newMSL$name==as.character(destTaxonName))
    destCurMatches=(.GlobalEnv$curMSL$name==as.character(destTaxonName))
    destOldMatches=(.GlobalEnv$oldMSLs$name==as.character(destTaxonName))
    
    #
    # check that at least one of dest/src taxa were given
    #
    if(is.na(srcTaxonName) && is.na(destTaxonName)) {
      log_change_error(changeDf=change,levelStr="ERROR", errorCode="XLSX.NO_SRC_DEST", 
                       errorStr=paste0("Change=",toupper(action),", but neither 'current taxonomy', nor 'proposed taxonomy', were specified")
      )
      next;
      
    }
    #
    # check if lineage not correct in "current" 
    #
    if( sum(srcNewMatches,na.rm=T)==1 ) {
      srcNewTaxon = .GlobalEnv$newMSL[srcNewMatches,]
      if(curChangeDf$.srcLineage != srcNewTaxon$lineage) {
         if(!is.na(srcNewTaxon$.otherLineage) ) {
          # this is a know change, done by another proposal
           log_change_error(curChangeDf=curChangeDf,levelStr="INFO", errorCode="SRC.LINEAGE_CHANGED", 
                            errorStr=paste0("CURRENT lineage changed by another proposal"), 
                            notes=paste0("proposal(s)=",srcNewTaxon$.otherLineageProposal, 
                                         " did a ",srcNewTaxon$.otherLineageAction, 
                                         " PROPOSAL//NEW_MSL=", diff_lineages(curChangeDf$.srcLineage,srcNewTaxon$lineage)
                            )
           )
        } 
        else {
          # unexplained lineage mismatch
          log_change_error(curChangeDf=curChangeDf,levelStr="WARNING", errorCode="SRC.LINEAGE_WRONG", 
                           errorStr=paste0("'current taxonomy' in proposal doesn't match current MSL. Typo?"), 
                           notes=paste0("MSL//PROPOSED=", diff_lineages(curChangeDf$.srcLineage,curChangeDf$.destLineage),
                                        ", PROPOSAL=", curChangeDf$.srcLineage, 
                                        "; MSL=", .GlobalEnv$newMSL[srcNewMatches,"lineage"])
          )
        } # unexplained lineage mismatch
      }
    } # one match
    #
    # make sure new taxon does not exist already
    #
    # only check if there IS a destTaxonName
    # move/merge/split are not checked:
    #    move is lineageB/taxonA to lineageC/taxonA
    #    merge usually is taxonA into existing_taxonB
    #   split can be to the same name, if src=dest
    if(!is.na(destTaxonName) && !(curChangeDf$.action %in% c("move","merge"))
       && !(curChangeDf$.action == 'split' && !is.na(srcTaxonName) && destTaxonName == srcTaxonName)) {
      #
      # verify we're not re-creating something that exists, or existed
      #
      
      # check if exists in the current MSL
      if(sum(destCurMatches)>0) {
        matchDf = .GlobalEnv$curMSL[destCurMatches,]
        matchDf = matchDf[order(matchDf$msl_release_num,decreasing = T)[1],]
        matchLineage = paste0( "MSL", matchDf$msl_release_num,":",matchDf$lineage)
        log_change_error(curChangeDf=curChangeDf,levelStr="ERROR", errorCode="DEST.IN_CUR", 
                         errorStr=paste0("Change=",toupper(action),", but taxon name already exists"), 
                         notes=paste0("proposed ",destTaxonRank,"=", destTaxonName, 
                                      "; existing ",matchDf$rank,"=", matchLineage)
        )
        next;
      }
      # check if exists in a historical MLS (not the most recent)
      if(sum(destOldMatches)>0) {
        matchDf = .GlobalEnv$oldMSLs[destOldMatches,]
        matchDf = matchDf[order(matchDf$msl_release_num,decreasing = T)[1],]
        matchLineage = paste0( "MSL", matchDf$msl_release_num,":",matchDf$lineage)
        log_change_error(curChangeDf,levelStr="ERROR", errorCode="DEST.IN_OLD", 
                         errorStr=paste0("Change=",toupper(action),",  but taxon name existed historically"), 
                         notes=paste0("proposed ",destTaxonRank,"=", destTaxonName, 
                                      "; existing ",matchDf$rank,"=", matchLineage)
        )
        next;
      }
      # check if already created in newMSL
      if(sum(destNewMatches)>0 ) {
        matchDf = .GlobalEnv$newMSL[destNewMatches,]
        matchDf = matchDf[order(matchDf$msl_release_num,decreasing = T)[1],]
        matchLineage = paste0( "MSL", matchDf$msl_release_num,":",matchDf$lineage)
        log_change_error(curChangeDf, "ERROR", "DEST.IN_NEW", 
                         errorStr=paste0("Change=",toupper(action),", but taxon name already created in new MSL"), 
                         notes=paste0("proposed ",destTaxonRank,"=", destTaxonName, 
                                      "; existing ",matchDf$rank,"=", matchLineage,
                                      "; otherProposal=", matchDf[1,"in_filename"]
                         )
        )
        next;
      }
      # WARN: check if taxon rank matches change rank
      if ( tolower(destTaxonRank) != tolower(curChangeDf$rank) ) {
        log_change_error(curChangeDf, "WARNING", "DEST.RANK_MISMATCH", 
                         errorStr=paste0("Change=",toupper(action)," but proposed taxon rank does not match [rank] column"), 
                         notes=paste0("rankColumn=", curChangeDf$rank,
                                      ", proposedTaxonRank=",destTaxonRank,
                                      ", proposedTaxonomy=", destLineage)
        )
      }
    } # if(!is.na(destTaxonName))
    
    # debug
    #if(code=="2022.003M" && interactive()) {browser()}
    
    #
    #### CREATE/SPLIT ####
    
    #
    #
    # Note: doesn't fix left_idx/right_idx, etc
    #
    # SPLIT (srcName != destName)
    #     split (same name) is handled by "move" coe
    #     can be new name, new species
    #     can be in any order
    # if there isn't one with same name or same accession by end of proposal, 
    # then delete original name in newMSL (must remember that!) 
    #   curMSL$.split_kept = T means the MOVE code saw a srcName = destName split
    #   curMSL$.split = T means this code saw a srcName!=destName split
    #
    if(    curChangeDf$.action %in% c("new")
       || (curChangeDf$.action %in% c("split") && srcTaxonName != destTaxonName) 
    ) {

      # index of our "current"/src taxon,  in the newMSL
      srcNewTaxonIdx = 0 # not yet found
      srcCurTaxonIdx = 0 
      
      #
      # check SRC taxon
      #
      if(curChangeDf$.action == "new" ) {
        # should have NO source taxon
        if( !is.na(srcTaxonName)) {
          # if srcTaxon was specified in xlsx (shouldn't be for NEW)  (curChangeDf$.destTaxon)
          # override .changeTaxon: set = .destTaxon
          # by default, it's .srcTaxon, if not empty, but the warning here is that it IS not empty when 
          # it should be, so the .destTaxon is more informative, if this is actually a split
          curChangeDf$.changeTaxon = curChangeDf$.destTaxon
          log_change_error(curChangeDf, "WARNING", "CREATE.W_SRC", 
                           errorStr="Change=CREATE, but 'current taxonomy' is not empty; perhaps you meant SPLIT", 
                           notes=paste0("CURRENT//PROPOSED=", diff_lineages(curChangeDf$.srcLineage,curChangeDf$.destLineage))
          )
        }
      } 
      else if(curChangeDf$.action == "split" ) {
#DEBUG        if(srcTaxonName == "Erythroparvovirus pinniped1"){browser()} #DEBUG
        if( is.na(srcTaxonName) ) {
          # split MUST have a source taxon (curChangeDf$.destTaxon)
          log_change_error(curChangeDf, "ERROR", "SPLIT.NO_SRC", 
                           errorStr="Change=SPLIT, but 'current taxonomy' columns are empty", 
                           notes=paste0("proposedTaxonomy=",curChangeDf$.destLineage)
          )
          next;
        } 
        else {
          # has a source taxon, let's find it
          srcCurTaxonIdx=(.GlobalEnv$curMSL$name==as.character(srcTaxonName))
          if( sum(srcCurTaxonIdx) == 0 ) {
            # split MUST have a source taxon  (curChangeDf$.destTaxon)
            log_change_error(curChangeDf, "ERROR", "SPLIT.SRC_NO_EXIST", 
                             errorStr="Change=SPLIT, but 'current taxonomy' doesn't exist", 
                             notes=paste0("currentTaxonomy=",curChangeDf$.srcLineage)
            )
            next;
          }
          # remember this taxon has been split - may need to 
          # obsolete it, later, if no split directive keeps the old name
          
          # mark original, just in case that's usefule
          .GlobalEnv$curMSL[srcCurTaxonIdx,".split"] = TRUE
          .GlobalEnv$curMSL[srcCurTaxonIdx,".split_code"] = code
          .GlobalEnv$curMSL[srcCurTaxonIdx,".split_linenum"] = linenum
          .GlobalEnv$curMSL[srcCurTaxonIdx,".split_actionOrder"] = actionOrder
          .GlobalEnv$curMSL[srcCurTaxonIdx,".actionOrder"] = actionOrder
          
          # mark the newMSL entry we may delete, if not already marked
          if(!.GlobalEnv$curMSL[srcCurTaxonIdx,".split_kept"]) {
            # NOT already seen by SPLIT=, so WE must mark
            #srcNewTaxonIdx=(.GlobalEnv$newMSL$name==as.character(srcTaxonName))
            srcNewTaxonIdx=(.GlobalEnv$newMSL$.prev_taxnode_id==.GlobalEnv$curMSL$taxnode_id[srcCurTaxonIdx])
            if( sum(srcNewTaxonIdx,na.rm=T)>0) {
              # mark for possible deletion, if not also flagged for preservation
              .GlobalEnv$newMSL[srcNewTaxonIdx,".split"] = TRUE 
              .GlobalEnv$newMSL[srcNewTaxonIdx,".split_code"] = code
              .GlobalEnv$newMSL[srcNewTaxonIdx,".split_linenum"] = linenum
              .GlobalEnv$newMSL[srcNewTaxonIdx,".split_actionOrder"] = actionOrder
            }
            else {
              # ERROR can't find src taxon in newMSL
              log_change_error(curChangeDf, "ERROR", "SPLIT.SRC_NOT_IN_NEW", 
                               errorStr="Change=SPLIT, but 'current taxonomy' can not be found in newMSL", 
                               notes=paste0("currentTaxonomy=",curChangeDf$.srcLineage)
              )
              next;
            }
          } # hasn't been marked by SPLIT=
        } # has srcTaxon
      } # split
     
      # 
      # check DEST taxon specified
      #
      if(is.na(destTaxonName)) {
        log_change_error(curChangeDf, "ERROR", "CREATE.NO_DEST",
                         errorStr=paste0("Change=",toupper(curChangeDf$.action),", but 'proposed taxonomy' columns are empty"), 
                         notes=paste0("CURRENT//PROPOSED=", diff_lineages(curChangeDf$.srcLineage, curChangeDf$.destLineage),
                                      ",CURRENT=", curChangeDf$.srcLineage,
                                      ", PROPOSED=",curChangeDf$.destLineage)
        )
        next;
      }
      
      #
      # check if same accession number already exists
      #
      isDupAccession = (.GlobalEnv$newMSL$genbank_accession_csv == curChangeDf$exemplarAccession)
      if(sum(isDupAccession, na.rm=TRUE)>0) {
        
        # unless this is a SPLIT doing a rename, but keeping the isolate/accession
        if( curChangeDf$.action == 'split' &&
            .GlobalEnv$curMSL$.split[srcCurTaxonIdx] &&
            !is.na(.GlobalEnv$curMSL$genbank_accession_csv[srcCurTaxonIdx]) &&
            .GlobalEnv$curMSL$genbank_accession_csv[srcCurTaxonIdx] == curChangeDf$exemplarAccession) {
  
          # this is ok for a split to re-use an accession under a new or same name, but just once
          if( !is.na(.GlobalEnv$curMSL$.split_acc_used[srcCurTaxonIdx]) ) {
             # can't use it more than once
            log_change_error(curChangeDf, "WARNING", "SPLIT.REUSE_ACC", 
                             errorStr=paste0("Change=",toupper(curChangeDf$.action),", accession re-used more than once in a split"), 
                             notes=paste0("accession=", curChangeDf$exemplarAccession,
                                    "; already reused on ", .GlobalEnv$curMSL$.split_acc_used[srcCurTaxonIdx])
            )
            next;
          } # illegal accession reuse
          else {
            # record that a split has used this accession
            .GlobalEnv$curMSL[srcCurTaxonIdx,".split_acc_used"] =paste0(code,":",linenum)
            if(params$tmi) {
              log_change_error(curChangeDf, "INFO", "SPLIT.REUSE_ACC", 
                               errorStr=paste0("Change=",toupper(curChangeDf$.action),", accession re-used in a split"), 
                               notes=paste0("accession=", curChangeDf$exemplarAccession,
                                      "; currentName=", srcTaxonName, "; proposedName=", destTaxonName)
              )
            } # tmi - accession reuse ok
          }
        } # acc re-use in split
        else {
          # not in a split
          
          # Pending/dup is hard error in FINAL mode, otherwise WARNING
          errLevel = ifelse(params$processing_mode=="final","ERROR","WARNING")
          
          # this is a warning/error depending on mode
          if(  curChangeDf$exemplarAccession %in% c("pending","Pending") ) { 
            log_change_error(curChangeDf, errLevel, "CREATE.PENDING_ACC", 
                             errorStr=paste0("Change=",toupper(curChangeDf$.action),", accession number is 'pending'"), 
                             notes=paste0("accession=", curChangeDf$exemplarAccession)
            )
          } 
          else {
            # hard error for each other species with this accession
            ## QQQ what proposal created this? (from this round? historically? Need a function: last_modified())
            log_change_error(curChangeDf, errLevel, "CREATE.DUP_ACC", 
                             errorStr=paste0("Change=",toupper(curChangeDf$.action),", a species with this accession number already exists"), 
                             notes=paste0("accession=", curChangeDf$exemplarAccession, ", existingSpecies=",.GlobalEnv$newMSL[isDupAccession,]$lineage)
            )
          }     
          # hard error in FINAL mode
          if(params$processing_mode=="final") { next }
        } # acc re-use not in split
         
      } # acc re-use
      
      # 
      # verify that PARENT taxon exists already in newMSL
      #
      #if(curChangeDf$.linenum == 8) {print("CREATE Orthodiscovirus oryzae: incorrect CREATE.PARENT_LINEAGE WARNING"); browser()}
      destParentNameAlias =""
      if(is.na(destParentName)) {
        # no parent - use root node
        parentDestNewMatches=(.GlobalEnv$newMSL$taxnode_id==.GlobalEnv$newMSL$tree_id)
      } 
      else {
        # newMSL: find parent by name (must be unique) 
        parentDestNewMatches=(.GlobalEnv$newMSL$name==as.character(destParentName))
        if(sum(parentDestNewMatches)==0 ) {
          # find parent by other means: look for name in old (curMSL), and assume it got renamed already in newMSL
          # and find it's new entry in newMSL
          prevDestParent = curMSL %>% filter(name==as.character(destParentName))
          if(nrow(prevDestParent)==1) {
            # get what that older name became
            parentDestNewMatches = (.GlobalEnv$newMSL$taxnode_id==prevDestParent$.out_taxnode_id)
            parentDestNew = .GlobalEnv$newMSL[parentDestNewMatches,]
            destParentNameAlias =    paste0(" [AKA ",parentDestNew$name,"]")
            # INFO message
            log_change_error(curChangeDf, "INFO", "CREATE.PARENT_ALREADY_CHANGED", 
                             errorStr="PROPOSED parent taxon already modified", 
                             notes=paste0("proposal(s)=", parentDestNew$.otherLineageProposal, " did a '", parentDestNew$.otherLineageAction,"' ")
            )
          } 
          else {
            # changeDf: look for another proposal that already changed the name in newMSL?
            changedParent = .GlobalEnv$allChangeDf %>% filter(.srcTaxon==as.character(destParentName))
            if(nrow(changedParent)==1) {
              # get what that older name became
              parentDestNewTaxon  = (.GlobalEnv$newMSL %>% filter(name==changedParent$.destTaxon))
 
              # INFO message
              log_change_error(curChangeDf, "INFO", "CREATE.PARENT_ALREADY_CHANGED_2",
                               errorStr="PROPOSED parent taxon already modified", 
                               notes=paste0("proposal(s)=", parentDestNewTaxon$.otherLineageProposal, " did a '", parentDestNewTaxon$.otherLineageAction,"' ")
              )
            } 
            else if(sum(parentDestNewMatches)>1) {
              log_change_error(curChangeDf, "ERROR", "CREATE.PARENT_MULTIPLE_MATCH", 
                                 errorStr="Probably parent was split by another proposal", 
                                 notes=paste0("proposal(s)=", changedParent$.code, " did a '",changedParent$change,"' ",
                                        "of '",changedParent$.srcTaxon, "' to '", changedParent$.destTaxon,"'")
                )
            } 
            else if(sum(parentDestNewMatches)==0) {
              # just missing
              log_change_error(curChangeDf, "ERROR", "CREATE.PARENT_NO_EXIST", 
                               errorStr=paste0("Change=",toupper(curChangeDf$.action),", but parent rank taxon does not exist"), 
                               notes=paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage)
              )
            } # just missing
          } # search changeDf
        } # search curMSL
      } # search newMSL
        
      if(params$verbose) {
        labelStr = ifelse(curChangeDf$.action=="split",
                          "CREATE(SPLIT):     ",
                          "CREATE:            ")
        cat(paste0(labelStr,toupper(curChangeDf$rank)," code:",code," line:",linenum," '",destTaxonName, "' findParent(",destParentName,destParentNameAlias,")=",sum(parentDestNewMatches)),"\n")
      }

      # no parent found: skip record, otherwise, get parent
      if(sum(parentDestNewMatches)!=1) {
        next;
      }
      # get actual parent
      destParentTaxon = .GlobalEnv$newMSL[parentDestNewMatches,]
      
      # WARN: check if taxon rank matches change rank
      if ( tolower(destTaxonRank) != tolower(curChangeDf$rank) ) {
        log_change_error(curChangeDf, "WARNING", "CREATE.RANK_MISMATCH", 
                         errorStr=paste0("Change=",toupper(curChangeDf$.action),", but proposed taxon rank does not match [rank] column"), 
                         notes=paste0("ankColumn=", curChangeDf$rank,
                                      ", proposedTaxonRank=",destTaxonRank,
                                      ", proposedTaxonomy=", destLineage)
        )
      }
      
      ###### binomial check ######
      # if new SPECIES, check binomial prefix matches parent genus
      #
      if( destTaxonRank=="species" ) {
        
        # 
        # check hostSource
        #
        if( is.na(curChangeDf$hostSource) || grepl("please.*select",curChangeDf$hostSource, ignore.case=T) ) {
          log_change_error(curChangeDf, "WARNING", "CREATE.SPECIES_NO_HOST_SOURCE", 
                           errorStr=paste0("Change=",toupper(curChangeDf$.action),", but proposed species must have a host/source value")
          )
        }
        
        # check if parent is a genus
        if(destParentTaxon$rank=="genus") {
          destParentGenusTaxon = destParentTaxon
        } 
        else {
          if(destParentTaxon$rank == "subgenus" ) {
            # check up one rank
            destParentGenusTaxon = .GlobalEnv$newMSL %>% filter(taxnode_id == destParentTaxon$parent_id)
            if(destParentGenusTaxon$rank != "genus") {
              # can't find genus parent
              log_change_error(curChangeDf, "ERROR", "CREATE.SPECIES_SUBGENUS_NO_GENUS", 
                               errorStr=paste0("Change=",toupper(curChangeDf$.action),", can not find parent genus"), 
                               notes=paste0("subgenus '",destParentTaxon$name,"' is not in a genus;",
                                            " it's parent '",destParentGenusTaxon$name,"' is a ",destParentGenusTaxon$rank )
              )
              next;
            }
          } 
          else {
            # parent isn't subgenus or genus
            log_change_error(curChangeDf, "ERROR", "CREATE.SPECIES_NO_GENUS", 
                             errorStr=paste0("Change=",toupper(curChangeDf$.action),", can not find parent genus"), 
                             notes=paste0("parent '",destParentTaxon$name,"' is a ",destParentTaxon$rank )
            )
            next;
          }
        } 
        # have genus (parent or grandparent)
        # check binomial naming
        if( str_detect(destTaxonName,paste0(destParentGenusTaxon$name," ")) != TRUE ) {
          log_change_error(curChangeDf, "ERROR", "CREATE.SPECIES_BINOMIAL_MISMATCH", 
                            errorStr=paste0("Change=",toupper(curChangeDf$.action),", but proposed species names does not start with 'genus[space]' per binomial naming convention"), 
                            notes=paste0("parent genus name =", destParentGenusTaxon$name)
          )
          next;
        }
        
      }
      
      # 
      #####  create new taxon #####
      #
    
      # get parent
      newTaxon = .GlobalEnv$newMSL[parentDestNewMatches,]
  
      # WARN if PARENT_LINEAGE is not expected, AND USE PARENT LINEAGE
      if( !is.na(destParentLineage) && (newTaxon$lineage != destParentLineage) ) {
        # check if parent was known to be changed
        if( !is.na(destParentTaxon$.otherLineage)  ) {
          if( curChangeDf$.destParentName == destParentTaxon$name ) {
            # cases where parent name itself changed would have been reported above as CREATE.PARENT_ALREADY_CHANGED
            # this is grandparent and above changes
            log_change_error(curChangeDf, "INFO", "CREATE.PARENT_RENAMED",
                             errorStr=paste0("Change=",toupper(curChangeDf$.action),", proposed parent taxon exists, but not with expected name/lineage, using observed lineage"),
                             notes=paste0("otherProposal(s)=", destParentTaxon$.otherLineageProposal, " did a ", destParentTaxon$.otherLineageAction,
                                          ", PROPOSED//OBSERVED=", diff_lineages(destParentLineage, newTaxon$lineage)
                                          #", PROPOSED=", destParentLineage,
                                          #", OBSERVED=",newTaxon$lineage,
                                          
                             )
            )
          }
        } 
        else {
          # WARNING: wasn't renamed, that we know of
          log_change_error(curChangeDf, "WARNING", "CREATE.PARENT_LINEAGE", 
                           errorStr=paste0("Change=",toupper(curChangeDf$.action),", proposed parent taxon exists, but not with expected lineage, using observed lineage"), 
                           notes=paste0("PROPOSED//OBSERVED=", diff_lineages(destParentLineage, newTaxon$lineage),
                                        ", otherProposals=",newTaxon$prev_proposals)
          )
        }
        destLineage=paste0(newTaxon$lineage,";",destTaxonName)
      }
      
      # add new info - primary columns
      newTaxon[1,"in_change"]   = curChangeDf$.action
      newTaxon[1,"in_filename"] = 
        ifelse(is.na(newTaxon[1,"in_filename"]),proposalZip,
               paste0(newTaxon[1,"in_filename"],";",proposalZip))
      newTaxon[1,"in_notes"]    = paste0("xlsx_row=",linenum)
      newTaxon[1,"in_target"]   = destLineage
      
      newTaxon[1,"name"]        = destTaxonName
      newTaxon[1,"cleaned_name"]= destTaxonName
      newTaxon[1,"level_id"]    = dbCvList[["rank"]]$id[dbCvList[["rank"]]$name==destTaxonRank]
      newTaxon[1,"rank"]        = destTaxonRank
      newTaxon[1,"parent_id"]   = newTaxon[1,"taxnode_id"]
      newTaxon[1,"taxnode_id"]  = max(.GlobalEnv$newMSL$taxnode_id)+1
      newTaxon[1,"ictv_id"]     = newTaxon$taxnode_id
      
      # propagate alternate names/lineages
      if( !is.na(destParentTaxon$.otherLineage) ) {
        newTaxon[1,".otherLineage"] =  paste0(destParentTaxon$.otherLineage,";",destTaxonName)
        newTaxon[1,".otherLineageProposal"] = destParentTaxon$.otherLineageProposal
        newTaxon[1,".otherLineageAction"]   = destParentTaxon$.otherLineageAction
      } 
      
      # genomeComposition = molecule_id 
      newTaxon[1,xlsx2dbMap["molecule"]] = ifelse(is.na(curChangeDf$molecule),NA,dbCvMapList[["molecule"]][curChangeDf$molecule])
      
      # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
      newTaxon[1,xlsx2dbMap["hostSource"]] = curChangeDf[1,"hostSource"] 
      
      # comments
      newTaxon[1,xlsx2dbMap["comments"]]= curChangeDf[1,"comments"]
      
      # change tracking
      newTaxon[1,".actionOrder"] = actionOrder
      #
      # if split, do admin for split
      #
      if( curChangeDf$.action == "split") {
        #
        # if the list of split directives does NOT include the current name, 
        # then we will have to (after the proposal is finished) abolish the copy
        # of the current name in the new MSL. Here, we mark that we've seen the 
        # split directive with the same name, so we should NOT abolish this one.
        # admin; mark we kept the original name in the split
        #.GlobalEnv$curMSL[srcPrevTarget,".split_kept"] = TRUE # old copy, just to know
        newTaxon[1,".split_kept"] = TRUE  # new copy - don't delete!
        newTaxon[1,".split_code"] = code  
        newTaxon[1,".split_actionOrder"] = actionOrder  
        newTaxon[1,".split_linenum"] = linenum  # new copy - mark which line saved
        
        # set IN change for SPLIT
        newTaxon[1,"in_change"] = curChangeDf$.action
        newTaxon[1,"in_filename"] = 
          ifelse(is.na(newTaxon[1,"in_filename"]),proposalZip,
                 paste0(newTaxon[1,"in_filename"],";",proposalZip))
        newTaxon[1,"in_target"] = srcLineage
        newTaxon[1,"in_notes"] = paste0("linenum=",linenum) # add comments?
      }
      
      #
      ##### create.species #####  
      #
      if( destTaxonRank == "species" ) {
        #
        ## for species only
        #
        
        # "genbank_accession_csv"
        newTaxon[1,xlsx2dbMap["exemplarAccession"]] = curChangeDf[1,"exemplarAccession"] 
        # exemplar_name
        # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
        newTaxon[1,xlsx2dbMap["exemplarName"]] = curChangeDf[1,"exemplarName"] 
        # "abbrev_csv"
        newTaxon[1,xlsx2dbMap["Abbrev"]] = curChangeDf[1,"Abbrev"] 
        # "isolate_csv"
        newTaxon[1,xlsx2dbMap["exemplarIsolate"]] = curChangeDf[1,"exemplarIsolate"] 
        # genome_coverage
        # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
        newTaxon[1,xlsx2dbMap["genomeCoverage"]]= curChangeDf[1,"genomeCoverage"] 
        # genome_coverage
        # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
        newTaxon[1,xlsx2dbMap["genomeCoverage"]]= curChangeDf[1,"genomeCoverage"] 
        # molecule_id is above (not species-specific)
      }
      
      # info-only columns - wont be saved to DB
      newTaxon[1,"rank"]       = destTaxonRank
      newTaxon[1,"lineage"]    = destLineage
      
      # clear some columns inherited from parent
      newTaxon[1,".prev_taxnode_id"] = NA
      newTaxon[1,"prev_proposals"] = paste0( proposalsDf[code,"basename"],".zip:",linenum) 
      
      # add new taxon to newMSL
      #if(params$tmi) {print(paste0("rbindlist(newMSL,",newTaxon[1,"taxnode_id"],":",destLineage,")"))}
      .GlobalEnv$newMSL <- rbindlist(list(.GlobalEnv$newMSL,newTaxon),fill=TRUE)
      
      # SUCCESS message
      errorCodeCreateMap=c("new"="CREATE.OK","split"="SPLIT.OK")
      #curChangeDf,levelStr, errorCode, errorStr, notes
      
      log_change_error(curChangeDf, 
                       levelStr="SUCCESS", errorCode=toupper(errorCodeCreateMap[curChangeDf$.action]), 
                       errorStr=paste0("Change=",toupper(curChangeDf$.action),", applied successfully"), 
                       notes=paste0("Create ", newTaxon$rank," of '",newTaxon$lineage,"'")
      )
    } # create new taxon (new, split!=)
    else if(curChangeDf$.action %in% c("rename") ) {
    
      # ........................................................................
      #
      #### RENAME ####  
      #
      # ........................................................................
      
      # check if srcTaxon was specified in xlsx (required)
      if(is.na(srcTaxonName)) {
        log_change_error(curChangeDf, "ERROR", "RENAME.WO_SRC", 
                         errorStr="Change=RENAME, but 'current taxonomy' columns are empty", 
                         notes=paste0("currentTaxonomy=", srcLineage, ", destTaxonomy=", destLineage)
        )
        next;
      }
      
      # check if renamed taxon has the same name (TP:Error-R1)
      if(srcTaxonName == destTaxonName) {
        log_change_error(curChangeDf, "WARNING", "RENAME.SAME_NAME", 
                         errorStr="Change=RENAME, but new name is the same; might this be a move?", 
                         notes=paste0("currentVsProposed=",diff_lineages(srcLineage, destLineage),
                                      ", currentTaxonomy=", srcLineage, ", destTaxonomy=", destLineage)
        )
      }
      
      
      #
      # find the target to rename in NEW
      #
      srcNewTarget=(.GlobalEnv$newMSL$name==as.character(srcTaxonName))
      if(params$verbose) {cat(paste0("RENAME:            ",toupper(curChangeDf$rank)," code:",code," line:",linenum," '",destTaxonName, "' findTarget(",srcTaxonName,")=",sum(srcNewTarget)),"\n")}
    
      if(sum(srcNewTarget)==0) {
        # check if someone else already modified it
        prevSrcTaxon = curMSL %>% filter(name==as.character(srcTaxonName))
        if(nrow(prevSrcTaxon)==0) {
          # just not found
          log_change_error(curChangeDf, "ERROR", "RENAME.NO_EXIST", "Change=RENAME, but taxon does not exist", 
                           errorStr=paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName)
          )
           next;
        } 
        else {
          # previous record - must have been modified already this MSL
          log_change_error(curChangeDf, "ERROR", "RENAME.SRC_ALREADY_CHANGED", "Change=RENAME, but taxon already modified", 
                           errorStr=paste0("proposal(s)=", prevSrcTaxon$out_filename, " also did a '",prevSrcTaxon$out_change,"' to '", prevSrcTaxon$out_target,"'")
          )
          next;
        }
      } 
      else if(sum(srcNewTarget)>1) {
        log_change_error(curChangeDf, "ERROR", "RENAME.MANY", "Change=RENAME, multiple taxa exist with parent name", 
                         errorStr=paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName
                                         , " matches ", sum(srcNewTarget) ))
         next;
      } 
      else {
        
        # find original taxon in prevMSL being renamed
        srcPrevTarget=(.GlobalEnv$curMSL$taxnode_id==.GlobalEnv$newMSL[srcNewTarget,]$.prev_taxnode_id)
        srcNewParent= (.GlobalEnv$newMSL$taxnode_id==.GlobalEnv$newMSL[srcNewTarget,]$parent_id )
        # print(paste("srcPrevTarget=",sum(srcPrevTarget),"srcNewParent=",sum(srcNewParent)))
        
        ###### binomial check ######
        
        # if rename SPECIES, check binomial prefix matches parent genus
        #
        if( destTaxonRank=="species" ) {
          # check if parent is a genus
          destParentTaxon = as.data.frame(newMSL[srcNewParent,])
          if(destParentTaxon$rank=="genus") {
            destParentGenusTaxon = destParentTaxon
          } 
          else {
            if(destParentTaxon$rank == "subgenus" ) {
              # check up one rank
              destParentGenusTaxon = .GlobalEnv$newMSL %>% filter(taxnode_id == destParentTaxon$parent_id)
              if(destParentGenusTaxon$rank != "genus") {
                # can't find genus parent
                log_change_error(curChangeDf, "ERROR", "RENAME.SPECIES_SUBGENUS_NO_GENUS", 
                                 errorStr="Change=RENAME, can not find parent genus", 
                                 notes=paste0("subgenus '",destParentTaxon$name,"' is not in a genus;",
                                              " it's parent '",destParentGenusTaxon$name,"' is a ",destParentGenusTaxon$rank )
                )
                next;
              }
            } 
            else {
              # parent isn't subgenus or genus
              log_change_error(curChangeDf, "ERROR", "CREATE.SPECIES_NO_GENUS", 
                               errorStr="Change=CREATE, can not find parent genus", 
                               notes=paste0("parent '",destParentTaxon$name,"' is a ",destParentTaxon$rank )
              )
              next;
            }
          } 
          # have genus (parent or grandparent)
          # check binomial naming
          if( str_detect(destTaxonName,paste0(destParentGenusTaxon$name," ")) != TRUE ) {
            log_change_error(curChangeDf, "ERROR", "RENAME.SPECIES_BINOMIAL_MISMATCH", 
                             errorStr="Change=RENAME, but proposed species names does not start with 'genus[space]' per binomial naming convention", 
                             notes=paste0("parent genus name =", destParentGenusTaxon$name)
            )
            next;
          }
          
        }
        
        ##### rename apply changes  #####
        
        # change the name and update lineage
        srcNewTargetOldLineage = .GlobalEnv$newMSL[srcNewTarget,"lineage"]
        .GlobalEnv$newMSL[srcNewTarget,"name"] = destTaxonName
        .GlobalEnv$newMSL[srcNewTarget,"lineage"] = paste(newMSL[srcNewParent,"lineage"],destTaxonName,sep=";") # should we do this? 
        .GlobalEnv$newMSL[srcNewTarget,".actionOrder"] = actionOrder
        
        # RECURSE TO SET LINEAGE OF dest's KIDS (RENAME)
        update_lineage(.GlobalEnv$newMSL$taxnode_id[srcNewTarget],.GlobalEnv$newMSL[srcNewTarget,"lineage"], 
                       # otherLineage
                       srcNewTargetOldLineage,
                       # otherLineageProposal
                       paste0(proposalZip, ":" , curChangeDf$.linenum), 
                       # otherLineageAction
                       paste0(curChangeDf$.action, " ", curChangeDf$.srcRank," ", curChangeDf$.srcTaxon, " to ", curChangeDf$.destTaxon) )
        
        # add to list of genera to check for binomial problems after proposal is complete
       if(  destTaxonRank == "genus" ) {
          renamedGenera = rbind(renamedGenera, 
                                data.frame(name =c(destTaxonName$genus), 
                                           taxnode_id=.GlobalEnv$newMSL$taxnode_id[srcNewTarget],
                                           code=code, 
                                           linenum= linenum,
                                           actionOrder = actionOrder)
          )
        } 
        
        # check if that is the expected lineage
        # RENAME: WARN that PARENT_LINEAGE is not expected, AND USE MSL/OBSERVED PARENT LINEAGE
        # DON'T NEED THIS - would already be reported above as SRC. LINEAGE_CHANGED
        # if( curChangeDf$.destParentName != .GlobalEnv$newMSL[srcNewParent,"name"] ) {
        #   parentRank = .GlobalEnv$newMSL[srcNewParent,]$rank
        #   errorDf=log_error(curChangeDf$.code,curChangeDf$.linenum,curChangeDf$.order,curChangeDf$change,curChangeDf$rank,curChangeDf$.changeTaxon,
        #                    "WARNING", "RENAME.PARENT_NAME", 
        #                    "parent name not as expected; was this meant to be a move? Is there a typo?", 
        #                    paste0("parent in this proposal: '", curChangeDf$.destParentName,"'[", curChangeDf$.destParentRank, "], ",
        #                          "but MSL parent is '",.GlobalEnv$newMSL[srcNewParent,"name"],"'[", .GlobalEnv$newMSL[srcNewParent,]$rank, "]; ",
        #                          "other new proposals: ",.GlobalEnv$newMSL[srcNewTarget,"prev_proposals"])
        #   )
        # }
        
        # append proposal
        .GlobalEnv$newMSL[srcNewTarget,"prev_proposals"] = paste0(
          ifelse(is.na(.GlobalEnv$newMSL[srcNewTarget,"prev_proposals"]),"",paste0(.GlobalEnv$newMSL[srcNewTarget,"prev_proposals"],",")),
          proposalZip,":",linenum)
        
        # put out_* changes on curMSL
        .GlobalEnv$curMSL[srcPrevTarget,"out_updated"] = TRUE  # admin; mark this to save to db
        .GlobalEnv$curMSL[srcPrevTarget,"out_change"] = "rename"
        .GlobalEnv$curMSL[srcPrevTarget,"out_filename"] = 
            ifelse(is.na(.GlobalEnv$curMSL[srcPrevTarget,"out_filename"]),proposalZip,
                   paste0(.GlobalEnv$curMSL[srcPrevTarget,"out_filename"],";",proposalZip))
        .GlobalEnv$curMSL[srcPrevTarget,"out_target"] = destTaxonName
        .GlobalEnv$curMSL[srcPrevTarget,"out_notes"] = paste0("linenum=",linenum)
        .GlobalEnv$curMSL[srcPrevTarget,".out_taxnode_id"] = .GlobalEnv$newMSL[srcNewTarget,"taxnode_id"]
        .GlobalEnv$curMSL[srcPrevTarget,".actionOrder"] = actionOrder
        
        # success note
        errorDf=log_change_error(curChangeDf, "SUCCESS", "RENAME.OK", "Change=RENAME, applied successfully", 
                                 paste0("RENAME ",srcTaxonRank, 
                                        " from ", "'",srcTaxonName,"'",
                                        " to ",   "'", .GlobalEnv$newMSL[srcNewTarget,"name"], "'",
                                        " in ",   "'", .GlobalEnv$newMSL[srcNewParent,"lineage"], "'"
                                 )
        )
        
      } # target found
    } # RENAME
    else if(curChangeDf$.action %in% c("abolish") ) { 
      # ........................................................................
      #
      ##### ABOLISH #####  
      #
      # ........................................................................
      # check if srcTaxon was specified in xlsx (required)
      if(is.na(srcTaxonName)) {
        errorDf=log_error(code=curChangeDf$.code,linenum=curChangeDf$.linenum,action=curChangeDf$change,actionOrder=curChangeDf$.order,
                          curChangeDf$rank,curChangeDf$.changeTaxon,
                         "ERROR", "ABOLISH.WO_SRC", "Change=ABOLISH, but 'current taxonomy' columns are empty", 
                         paste0("destTaxonomy=", destLineage)
        )
        next;
      }
      
      #
      # find the target to abolish in NEW and CUR
      #
      srcNewTarget=(.GlobalEnv$newMSL$name==as.character(srcTaxonName))
      srcPrevTarget=(.GlobalEnv$curMSL$taxnode_id==.GlobalEnv$newMSL[srcNewTarget,]$.prev_taxnode_id)

      if(params$verbose) {cat(paste0("ABOLISH:           ",code," line ",linenum," '",srcTaxonName, "' findTarget(",srcTaxonName,")=",sum(srcNewTarget),"/",sum(srcPrevTarget)),"\n")}
      
      if(sum(srcNewTarget,na.rm=TRUE)==0) {
        # check if someone else already modified it
        prevSrcTaxon = curMSL %>% filter(name==as.character(srcTaxonName))
        if(nrow(prevSrcTaxon)==0 ) {
          # not modified, just missing
          log_change_error(curChangeDf, "ERROR", "ABOLISH.NO_EXIST", "Change=ABOLISH, but taxon does not exist", 
                           notes=paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName))
          next;
        } 
        else {
          # previous record - must have been modified already this MSL
          log_change_error(curChangeDf, "ERROR", "ABOLISH.SRC_ALREADY_CHANGED", "Change=ABOLISH, but taxon already modified", 
                            notes=paste0("proposal(s)=", prevSrcTaxon$out_filename, " also did a '",prevSrcTaxon$out_change,"' to '", prevSrcTaxon$out_target,"'")
          )
          next;
        }
      } 
      else if(sum(srcNewTarget,na.rm=TRUE)>1) {
         # this should never happen
        log_change_error(curChangeDf, "ERROR", "ABOLISH.MANY", "Change=ABOLISH, multiple taxa exist with name", 
                         notes=paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName
                                      , " matches ", sum(srcNewTarget) ))
        next;
      }
      #
      # check for kids - can't abolish before kids
      #
      srcKids=(.GlobalEnv$newMSL$parent_id==.GlobalEnv$newMSL[srcNewTarget,"taxnode_id"])
      if(sum(srcKids,na.rm=TRUE)>0) {
        log_change_error(curChangeDf, "ERROR", "ABOLISH.KIDS", "Change=ABOLISH, taxon still has un-abolished children", 
                         notes=paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", kids: N=",sum(srcKids,na.rm=TRUE),
                                      ", NAMES=[", paste(.GlobalEnv$newMSL$rank[srcKids],.GlobalEnv$newMSL$name[srcKids],sep=":"),"]" ) 
        )
        next;
      }     
      
      #
      # remove taxon from NEW
      #
      .GlobalEnv$newMSL = .GlobalEnv$newMSL[!srcNewTarget,]
      
      # put out_* changes on curMSL
      .GlobalEnv$curMSL[srcPrevTarget,"out_updated"] = TRUE  # admin; mark this to save to db
      .GlobalEnv$curMSL[srcPrevTarget,"out_change"] = "abolish"
      .GlobalEnv$curMSL[srcPrevTarget,"out_filename"] = 
        ifelse(is.na(.GlobalEnv$curMSL[srcPrevTarget,"out_filename"]),proposalZip,
               paste0(.GlobalEnv$curMSL[srcPrevTarget,"out_filename"],";",proposalZip))
      .GlobalEnv$curMSL[srcPrevTarget,"out_target"] = destTaxonName
      .GlobalEnv$curMSL[srcPrevTarget,"out_notes"] = paste0("linenum=",linenum) # add comments?
      .GlobalEnv$curMSL[srcPrevTarget,".actionOrder"] = actionOrder
      
      # success note
      log_change_error(curChangeDf, "SUCCESS", "ABOLISH.OK", "Change=ABOLISH, applied successfully", 
                       notes=paste0("ABOLISH ",srcTaxonRank, " named ",srcTaxonName)
      )
      
    } # ABOLISH
    else if(curChangeDf$.action %in% c("move","promote","demote") 
            || ( curChangeDf$.action =='split' && srcTaxonName == destTaxonName) 
            || ( curChangeDf$.action =='merge' && srcTaxonName == destTaxonName) 
    ) { 
      # ........................................................................
      #
      ##### MOVE,PROMOTE/DEMOTE,SPLIT=,MERGE= ##### 
      #
      # Only SPLIT where srcName = destName
      # MOVE w/ and w/o rename
      #
      # these all require a SRC and DEST, and change structure
      #
      # ........................................................................
 
      # check if srcTaxon was specified in xlsx (required)
      if(is.na(srcTaxonName) && curChangeDf$.action =='move') {
        # they didn't specify the source, but  it's a MOVE
        # we should be able to find it
        guessSource = (.GlobalEnv$newMSL$name==as.character(destTaxonName))
        if( sum(guessSource,na.rm=TRUE) == 1) {
          # found it, so WARNING
          log_change_error(curChangeDf, "WARNING", "MOVE.NO_SRC_CAN_GUESS", "Change=MOVE, but 'current taxonomy' columns are empty", 
                           notes=paste0("destTaxonomy=", destLineage,
                                        "; we guess you meant=", .GlobalEnv$newMSL$lineage[guessSource])
          )
          srcTaxonName=.GlobalEnv$newMSL$name[guessSource]
          srcTaxonRank=.GlobalEnv$newMSL$rank[guessSource]
          srcLineage  =.GlobalEnv$newMSL$lineage[guessSource]
        } 
        else {
          # didn't find it, or found multiple, so ERROR
          log_change_error(curChangeDf, "WARNING", "MOVE.NO_SRC", "Change=MOVE, but 'current taxonomy' columns are empty", 
                           notes=paste0("destTaxonomy=", destLineage,
                                        "; no existing taxon '",destTaxonName,"' found")
          )
          next
        }
      } 
      
      ##### validate rank #####
      #   1. consistent
      #   2. move does not change
      #   3. promote/demote does
      #
      if( destTaxonRank != rankClean ) {
        # 1. consistent: src rank doesn't match action rank
        log_change_error(curChangeDf, "ERROR", "MOVE.DEST_WRONG_RANK", 
                         errorStr=paste0("Change=",toupper(curChangeDf$.action),", but 'proposed taxonomy' ranks doesn't match 'proposed rank' column"), 
                         notes=paste0("rank=", rankClean, ", srcTaxonomy=[", srcTaxonRank,"]",srcLineage)
        )
        next;
      }
      
      if( curChangeDf$.action %in% c("move") && srcTaxonRank != destTaxonRank ) {
        # 2. MOVE: src/dest ranks NOT the same
        log_change_error(curChangeDf, "ERROR", "MOVE.DIFF_RANK", 
                         errorStr="Change=MOVE, but 'current taxonomy' and 'proposed taxonomy' ranks differ; use promote or demote, not move, to change rank", 
                         notes=paste0("CUR//PROPOSED=[",srcTaxonRank,"//",destTaxonRank,"] ", diff_lineages(srcLineage,destLineage)
                         )
        )
        next;
      } 
      else if( curChangeDf$.action %in% c("promote","demote") && srcTaxonRank == destTaxonRank ) {
        # 3. PROMOTE/DEMOTE src/dest ranks ARE the same 
        log_change_error(curChangeDf, "ERROR", "CHANGE_RANK.SAME_RANK", 
                         errorStr=paste0("Change=",toupper(curChangeDf$.action),", but 'current taxonomy' and 'proposed taxonomy' ranks are the same; use MOVE"), 
                         notes=paste0("rank=", rankClean, 
                                      ", CUR//PROPOSED=",diff_lineages(srcLineage,destLineage))
                         
        )
        next;
      }
      
      #
      # find the SRC target to MOVE in NEW and CUR
      #
      srcNewTarget =(.GlobalEnv$newMSL$name==as.character(srcTaxonName))
      srcPrevTarget=(.GlobalEnv$curMSL$taxnode_id==.GlobalEnv$newMSL[srcNewTarget,]$.prev_taxnode_id)

      if(params$verbose) {cat(paste0("MOVE:              ",toupper(curChangeDf$rank)," code:",code," line:",linenum," '",srcTaxonName, "' findTarget(",srcTaxonName,")=",sum(srcNewTarget),"/",sum(srcPrevTarget)),"\n")}
    
      if(sum(srcNewTarget,na.rm=TRUE)==0) {
        # check prevMSL, to see if something else already moved it
        prevSrcTaxon = curMSL %>% filter(name==as.character(srcTaxonName))
        if(nrow(prevSrcTaxon)==0 ) {
          # no previous record, just doesn't exist
          log_change_error(curChangeDf, "ERROR", "MOVE.NO_EXIST", 
                           errorStr=paste0("Change=",toupper(curChangeDf$.action),", but taxon does not exist"), 
                           notes=paste0("taxon=", srcTaxonName, ",CUR//PROPOSED=",diff_lineages(srcLineage,destLineage),
                                        ", CUR=",srcLineage,", PROPOSED=", destLineage))
          next;
        } 
        else {
          # previous record - must have been modified already this MSL
          log_change_error(curChangeDf, "ERROR", "MOVE.SRC_ALREADY_CHANGED", 
                           errorStr=paste0("Change=",toupper(curChangeDf$.action),", but taxon already modified"), 
                           notes=paste0("proposal(s)=", prevSrcTaxon$out_filename, 
                                        " also did a '",prevSrcTaxon$out_change,"'",
                                        " to '", prevSrcTaxon$out_target,"'")
          )
          next;
        }
      } 
      else if(sum(srcNewTarget,na.rm=TRUE)>1) {
        log_change_error(curChangeDf, "ERROR", "MOVE.MANY", 
                         errorStr=paste0("Change=",toupper(curChangeDf$.action),", multiple taxa exist with name"), 
                         notes=paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName
                                      , " matches ", sum(srcNewTarget) ))
         next;
      }

      #
      # check if same accession number already exists
      #
      #  !srcNewTarget - prevents checking against ourself
      #
      isDupAccession = (.GlobalEnv$newMSL$genbank_accession_csv == curChangeDf$exemplarAccession & !srcNewTarget)
      if(sum(isDupAccession, na.rm=TRUE)>0) {
        log_change_error(curChangeDf, "ERROR", "MOVE.DUP_ACC", 
                         errorStr=paste0("Change=",toupper(curChangeDf$.action),", a species with this accession number already exists"), 
                         notes=paste0("accession=", curChangeDf$exemplarAccession, ", existingSpecies=",.GlobalEnv$newMSL[isDupAccession,]$lineage)
        )
        ## QQQ what proposal created this? (from this round? historically? Need a function: last_modified())
        next;
      }
      
      # 
      # verify that PARENT taxon exists already in newMSL
      #
      parentDestNewMatches = c()
      if(is.na(destParentName)) {
        # no parent - use root node
        parentDestNewMatches=(.GlobalEnv$newMSL$taxnode_id==.GlobalEnv$newMSL$tree_id)
      } 
      else {
        # find parent by name (must be unique)
        parentDestNewMatches=(.GlobalEnv$newMSL$name==as.character(destParentName))
      }
      
      if(params$verbose) {cat(paste0("MOVE:              ",toupper(curChangeDf$rank)," code:",code," line ",linenum," '",destTaxonName, "' findParent(",destParentName,")=",sum(parentDestNewMatches)),"\n")}
    
      if(sum(parentDestNewMatches,na.rm=TRUE)==0) {
        # check if it used to exist
        prevDestParent = curMSL %>% filter(name==as.character(destParentName))
        if(nrow(prevDestParent)==0) {
          # just missing
          log_change_error(curChangeDf, "ERROR", "MOVE.PARENT_NO_EXIST", 
                           errorStr=paste0("Change=",toupper(curChangeDf$.action),", but parent rank taxon does not exist"), 
                           notes=paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage)
          )
          next;
        } 
        else {
          # someone already modified it - tell me who!
          log_change_error(curChangeDf, "ERROR", "MOVE.PARENT_ALREADY_CHANGED", 
                           errorStr=paste0("Change=",toupper(curChangeDf$.action),", but parent taxon already modified"), 
                           notes=paste0("proposal(s)=", prevDestParent$out_filename, " also did a '",prevDestParent$out_change,"' to '", prevDestParent$out_target,"'")
          )
          next;
        }
      } 
      else if(sum(parentDestNewMatches,na.rm=TRUE)>1) {
        # this should never happen
        log_change_error(curChangeDf, "ERROR", "MOVE.PARENT_MANY", 
                         errorStr=paste0("Change=",toupper(curChangeDf$.action),", multiple taxa exist with parent name"), 
                         notes=paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage)
        )
        next;
      } 
      
      # remember to check later if a genus name change broke binomial
      if(  destTaxonRank == "genus" ) {
        renamedGenera = rbind(renamedGenera, 
                              data.frame(name =c(destTaxonName$genus), 
                                         taxnode_id=.GlobalEnv$newMSL$taxnode_id[srcNewTarget],
                                         code=code, 
                                         linenum= linenum,
                                         actionOrder = actionOrder)
        )
       } 
       
      #
      # check that promote/demote change name
      #
      if( curChangeDf$.action %in% c("promote","demote") ) {
          if(srcTaxonName == destTaxonName) {
            log_change_error(curChangeDf, "WARNING", "CHANGE_RANK.SAME_NAME", 
                             errorStr=paste0("Change=",toupper(curChangeDf$.action),", but current and proposed names are the same."), 
                             notes=paste0("current name=[",srcTaxonRank,"]", srcTaxonNameprevDestParent$out_filename, " and proposed name [",destTaxonRank,"]",destTaxonName)
            )
          }
      }
      
      # 
      # MOVE the taxon
      #

      # get new parent
      destParentTaxon = .GlobalEnv$newMSL[parentDestNewMatches,]
       # WARN if PARENT_LINEAGE is not expected
      if( !is.na(destParentLineage) && (destParentTaxon$lineage != destParentLineage) ) {
        # check if this is a known change
        
        if( !is.na(destParentTaxon$.otherLineage)  ) { 
          # a parent of the destParent was renamed, moved, etc.
          log_change_error(curChangeDf, "INFO", "MOVE.PROPOSED_PARENT_LINEAGE_CHANGED", 
                           errorStr=paste0("PROPOSED parent lineage already modified"), 
                           notes=paste0("proposal(s)=",destParentTaxon$.otherLineageProposal, 
                                  " did a ", destParentTaxon$.otherLineageAction, 
                                  "; OBSERVED//PROPOSED=", diff_lineages(destParentTaxon$lineage, destParentLineage)
                           )
          )
        } else if( exists("prevDestParent") && nrow(prevDestParent) && (prevDestParent$lineage != destParentLineage) ) {
          # lineage mismatch not explained by curMSL lineage
          log_change_error(curChangeDf,  "WARNING", "MOVE.PARENT_LINEAGE", 
                           errorStr=paste0("PROPOSED parent taxon exists, but not with expected lineage"), 
                           notes=paste0(",PROPOSED//CUR=",diff_lineages(destParentLineage,destParentTaxon$lineage),
                                        ", PROPOSED=", destParentLineage,
                                        ", CUR=",destParentTaxon$lineage,
                                        ", otherProposals=",destParentTaxon$prev_proposals)
          )
        } else {
          # destParentLineage matches original (curMSL) lineage. 
          log_change_error(curChangeDf, "INFO", "MOVE.PARENT_LINEAGE_UPDATED2", 
                           errorStr=paste0("PROPOSED parent taxon exists, but with updated lineage"), 
                           notes=paste0(",PROPOSED//CUR=",diff_lineages(destParentLineage,destParentTaxon$lineage),
                                        ", PROPOSED=", destParentLineage,
                                        ", CUR=",destParentTaxon$lineage
                                        #", otherProposals=",destParentTaxon$prev_proposals)
                           )
          )
          
        } # dest lineage known update
        
      } # dest lineage unexpected
      
      ##### binomial check (species) #####
      #
      # if move SPECIES, check binomial prefix matches parent genus
      #
      if( destTaxonRank=="species" ) {
        # check if parent is a genus
        if(destParentTaxon$rank=="genus") {
          destParentGenusTaxon = destParentTaxon
        } else {
          if(destParentTaxon$rank == "subgenus" ) {
            # check up one rank
            destParentGenusTaxon = .GlobalEnv$newMSL %>% filter(taxnode_id == destParentTaxon$parent_id)
            if(destParentGenusTaxon$rank != "genus") {
              # can't find genus parent
              log_change_error(curChangeDf, "ERROR", "MOVE.SPECIES_SUBGENUS_NO_GENUS", 
                               errorStr=paste0("Change=",toupper(curChangeDf$.action),", can not find parent genus"), 
                               notes=paste0("subgenus '",destParentTaxon$name,"' is not in a genus;",
                                            " it's parent '",destParentGenusTaxon$name,"' is a ",destParentGenusTaxon$rank )
              )
              next;
            }
          } else {
            # parent isn't subgenus or genus
            log_change_error(curChangeDf, "ERROR", "MOVE.SPECIES_NO_GENUS", 
                             errorStr=paste0("Change=",toupper(curChangeDf$.action),", can not find parent genus"), 
                             notes=paste0("parent '",destParentTaxon$name,"' is a ",destParentTaxon$rank )
            )
            next;
          }
        } 
        # have genus (parent or grandparent)
        # check binomial naming
        if( str_detect(destTaxonName,paste0(destParentGenusTaxon$name," ")) != TRUE ) {
          log_change_error(curChangeDf, "ERROR", "MOVE.SPECIES_BINOMIAL_MISMATCH", 
                           errorStr=paste0("Change=",toupper(curChangeDf$.action),", but proposed species names does not start with 'genus[space]' per binomial naming convention"), 
                           notes=paste0("parent genus name =", destParentGenusTaxon$name)
          )
          next;
        }
      }
      
      
      # add new info - primary columns
      .GlobalEnv$newMSL[srcNewTarget,"name"]       = destTaxonName
      # this part shouldn't be a change unless this is promote/demote....
      .GlobalEnv$newMSL[srcNewTarget,"level_id"]   = dbCvList[["rank"]]$id[dbCvList[["rank"]]$name==destTaxonRank]
      .GlobalEnv$newMSL[srcNewTarget,"rank"]       = destTaxonRank
      .GlobalEnv$newMSL[srcNewTarget,"parent_id"]  = destParentTaxon[1,"taxnode_id"]
      .GlobalEnv$newMSL[srcNewTarget,".actionOrder"]  = actionOrder
      
      # genomeComposition = molecule_id 
      if(!is.na(curChangeDf[1,"molecule"])) {
        .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["molecule"]] = dbCvMapList[["molecule"]][curChangeDf[1,"molecule"] ]
      }
      
      # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
      if(!is.na(curChangeDf[1,"hostSource"])) {
        .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["hostSource"]] = curChangeDf[1,"hostSource"] 
      }
      # comments
      if(!is.na(curChangeDf[1,"comments"])) {
        .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["comments"]]= curChangeDf[1,"comments"]
      }
      
      ## for species only
      if( destTaxonRank == "species" ) {
        # "genbank_accession_csv"
        if(!is.na(curChangeDf[1,"exemplarAccession"])) {
          .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["exemplarAccession"]] = curChangeDf[1,"exemplarAccession"] 
        }
        # exemplar_name
        # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
        if(!is.na(curChangeDf[1,"exemplarName"])) {
          .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["exemplarName"]] = curChangeDf[1,"exemplarName"] 
        }
        # "abbrev_csv"
        if(!is.na(curChangeDf[1,"Abbrev"])) {
          .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["Abbrev"]] = curChangeDf[1,"Abbrev"] 
        }
        # "isolate_csv"
        if(!is.na(curChangeDf[1,"exemplarIsolate"])) {
          .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["exemplarIsolate"]] = curChangeDf[1,"exemplarIsolate"] 
        }
        # genome_coverage
        # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
        if(!is.na(curChangeDf[1,"genomeCoverage"])) {
          .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["genomeCoverage"]]= curChangeDf[1,"genomeCoverage"] 
        }
      }
      
      #  RECURSE TO SET LINEAGE OF destPARENT's KIDS (MOVE)
      moveDiff = diff_lineages(curChangeDf$.srcLineage,.GlobalEnv$newMSL$lineage[srcNewTarget])
      update_lineage(.GlobalEnv$newMSL$taxnode_id[srcNewTarget],.GlobalEnv$newMSL$lineage[srcNewTarget],
                     paste0(destParentTaxon$.otherLineage,";",.GlobalEnv$newMSL$name[srcNewTarget]),
                     # otherLineageProposal
                     paste0(proposalZip, ":" , curChangeDf$.linenum), 
                     # otherLineageAction
                     paste0(curChangeDf$.action, " ", curChangeDf$.srcRank," ", moveDiff )
      )
      
      #
      # update curMSL [out_*] (for non-splits)
      #
      if( curChangeDf$.action == 'split' ) {
        # SPLIT=
        # if the list of split directives does NOT include the current name, 
        # then we will have to (after the proposal is finished) abolish the copy
        # of the current name in the new MSL. Here, we mark that we've seen the 
        # split directive with the same name, so we should NOT abolish this one.
        # admin; mark we kept the original name in the split
        .GlobalEnv$curMSL[srcPrevTarget,".split_kept"] = TRUE # old copy, just to know
        .GlobalEnv$newMSL[srcNewTarget,".split_kept"] = TRUE  # new copy - don't delete!
        .GlobalEnv$newMSL[srcNewTarget,".split_code"] = code  
        .GlobalEnv$newMSL[srcNewTarget,".split_actionOrder"] = actionOrder  
        .GlobalEnv$newMSL[srcNewTarget,".split_linenum"] = curChangeDf$.linenum  # new copy - mark which line saved
        
        # set IN change for SPLIT
        .GlobalEnv$newMSL[srcNewTarget,"in_change"] = curChangeDf$.action
        .GlobalEnv$newMSL[srcNewTarget,"in_filename"] = 
            ifelse(is.na(.GlobalEnv$newMSL[srcNewTarget,"in_filename"] ),proposalZip,
            paste0(.GlobalEnv$newMSL[srcNewTarget,"in_filename"],";",proposalZip))
        .GlobalEnv$newMSL[srcNewTarget,"in_target"] = srcLineage
        .GlobalEnv$newMSL[srcNewTarget,"in_notes"] = paste0("linenum=",linenum) # add comments?
      } 
      else if( curChangeDf$.action == 'merge' ) {
        # MERGE=
        # if the list of split directives does NOT include the current name, 
        # then we will have to (after the proposal is finished) abolish the copy
        # of the current name in the new MSL. Here, we mark that we've seen the 
        # split directive with the same name, so we should NOT abolish this one.
        # admin; mark we kept the original name in the split
        .GlobalEnv$curMSL[srcPrevTarget,".merge_kept"] = TRUE # old copy, just to know
        .GlobalEnv$newMSL[srcNewTarget, ".merge_kept"] = TRUE  # new copy - don't delete!
        .GlobalEnv$newMSL[srcNewTarget, ".merge_code"] = code  
        .GlobalEnv$newMSL[srcNewTarget, ".merge_actionOrder"] = actionOrder  
        .GlobalEnv$newMSL[srcNewTarget, ".merge_linenum"] = curChangeDf$.linenum  # new copy - mark which line saved
        
        # set OUT change for MERGE
        .GlobalEnv$curMSL[srcPrevTarget,"out_change"] = curChangeDf$.action
        .GlobalEnv$curMSL[srcPrevTarget,"out_filename"] = 
          ifelse(is.na(.GlobalEnv$curMSL[srcPrevTarget,"out_filename"]),proposalZip,
                 paste0(.GlobalEnv$curMSL[srcPrevTarget,"out_filename"],";",proposalZip))
        .GlobalEnv$curMSL[srcPrevTarget,"out_target"] = destLineage
        .GlobalEnv$curMSL[srcPrevTarget,"out_notes"] = paste0("linenum=",linenum) # add comments?
        .GlobalEnv$curMSL[srcPrevTarget,".actionOrder"] = actionOrder
      } 
      else {
        # set OUT change for all others
        .GlobalEnv$curMSL[srcPrevTarget,"out_updated"] = TRUE  # admin; mark this to save to db
        .GlobalEnv$curMSL[srcPrevTarget,"out_change"] = curChangeDf$.action
        .GlobalEnv$curMSL[srcPrevTarget,"out_filename"] = 
          ifelse(is.na(.GlobalEnv$curMSL[srcPrevTarget,"out_filename"]),proposalZip,
                 paste0(.GlobalEnv$curMSL[srcPrevTarget,"out_filename"],";",proposalZip))
        .GlobalEnv$curMSL[srcPrevTarget,"out_target"] = destLineage
        .GlobalEnv$curMSL[srcPrevTarget,"out_notes"] = paste0("linenum=",linenum) # add comments?
        .GlobalEnv$curMSL[srcPrevTarget,".out_taxnode_id"] = .GlobalEnv$newMSL[srcNewTarget,"taxnode_id"] 
        .GlobalEnv$curMSL[srcPrevTarget,".actionOrder"] = actionOrder
      }
      
      # success note
      errorCodeMoveMap=c("move"="MOVE.OK","promote"="PROMOTE.OK","demote"="DEMOTE.OK","split"="SPLIT=.OK","merge"="MERGE=.OK")
      log_change_error(curChangeDf, "SUCCESS",  errorCodeMoveMap[curChangeDf$.action], 
                       errorStr=paste0("Change=",toupper(curChangeDf$.action),", applied successfully"), 
                       notes=paste0(toupper(curChangeDf$.action)," ",
                                    .GlobalEnv$curMSL$rank[srcPrevTarget], " named '", .GlobalEnv$curMSL[srcPrevTarget,"name"],    "'", 
                                    " to ", .GlobalEnv$newMSL$rank[srcNewTarget], " named '", .GlobalEnv$newMSL[srcNewTarget,"name"],    "'", 
                                    " PROPOSED//CUR=",  diff_lineages(
                                      .GlobalEnv$newMSL[srcNewTarget, "lineage"],
                                      .GlobalEnv$curMSL[srcPrevTarget,"lineage"]
                                    )
                                    
                       )
      )
    } # MOVE/PROMOTE/DEMOTE/SPLIT=/MERGE=
    else if(curChangeDf$.action %in% c("merge") ) { 
      # ........................................................................
      #
      ###### MERGE != ######
      #
      # this can be 3 cases
      # curMSL    newMSL  case
      #   A   ->   A      (A exists; just set out_*)  [implemneted with move & split=]
      #   A   ->   B      (B exists; just set out_*)
      #   A   ->   C      (C not exits; like rename: set out_*)
      # 
      # also may need to delete A from newMSL. Mark and check afterwards.
      # ........................................................................
      
      # check if srcTaxon was specified in xlsx (required)
      if(is.na(srcTaxonName)) {
        log_change_error(curChangeDf, "ERROR", "MERGE.WO_SRC", 
                         errorStr="Change=MERGE, but 'current taxonomy' columns are empty", 
                         notes=paste0("currentTaxonomy=", srcLineage, ", destTaxonomy=", destLineage)
        )
        next;
      }
      
      # check if MERGED taxon has the same name 
      # that should have been handled above with SPLIT= and MOVE
      if(srcTaxonName == destTaxonName) {
        log_change_error(curChangeDf, "WARNING", "MERGE.INTERNAL_SAME_NAME", 
                         errorStr="Change=MERGE, current=proposed, internal error", 
                         notes=paste0("MERGE(current=proposed) should be handled elsewhere;",
                                      " programming error")
        )
      }
      
      
      #
      ##### find SRC #####
      #
      # first look in curMSL
      srcCurTaxon = curMSL %>% filter(name==as.character(srcTaxonName))
      if(params$verbose) {cat(paste0("MERGE:             ",toupper(curChangeDf$rank)," code:",code," line:",linenum," src.findTarget(cur,",srcTaxonName,")=",nrow(srcCurTaxon),"\n"))}
      if(nrow(srcCurTaxon)==0) {
        # not found in curMSL - see if it's been renamed, look in newMSL
        srcNewTaxon= newMSL %>% filter(name==as.character(srcTaxonName))
        if(params$verbose) {cat(paste0("MERGE:             ",toupper(curChangeDf$rank)," code:",code," line:",linenum," '", "' src.findTarget(new,",srcTaxonName,")=",nrow(srcNewTaxon),"\n"))}
        if(nrow(srcNewTaxon)==0 ) {
          # not found anywhere
          log_change_error(curChangeDf, "ERROR", "MERGE.SRC_NOT_FOUND", "Change=MERGE, but taxon not found", 
                           errorStr=paste0("Looked for=",srcTaxonName,
                                           "; merge(CURRENT//PROPOSED)=", diff_lineages(srcLineage,destLineage))
          )
          next;
        }
        else { 
          # found in newMSL, but not in curMSL; go find original curMSL name
          srcCurTaxon = curMSL %>% filter(taxnode_id==srcNewTaxon$.prev_taxnode_id)
          # warn about discrepancy
          log_change_error(curChangeDf, "WARN", "MERGE.SRC_ALREADY_CHANGED", "Change=MERGE, but current taxon already changed", 
                           errorStr=paste0("Looked for=",srcTaxonName,
                                           "; ORIGNAL//FOUND=", diff_lineages(srcCurTaxon$lineage,srcNewTaxon$lineage),
                                           "; proposal=",," already did=",)
          )
        }
      }
      
      ##### tag SRC in newMSL #####

      # create selectors for update (should re-write to use data.table update) 
      srcCurSelect = curMSL$taxnode_id == srcCurTaxon$taxnode_id
      srcNewSelect = newMSL$taxnode_id == srcNewTaxon$taxnode_id

      # mark original, just in case that's usefule
      .GlobalEnv$curMSL[srcCurSelect,".merge"] = TRUE
      .GlobalEnv$curMSL[srcCurSelect,".merge_code"] = code
      .GlobalEnv$curMSL[srcCurSelect,".merge_linenum"] = linenum
      .GlobalEnv$curMSL[srcCurSelect,".merge_actionOrder"] = actionOrder
      .GlobalEnv$curMSL[srcCurSelect,".actionOrder"] = actionOrder
      
      # mark the newMSL entry for possible, if no MERGE= elsewhere
      .GlobalEnv$newMSL[srcNewSelect,".merge"] = TRUE 
      .GlobalEnv$newMSL[srcNewSelect,".merge_code"] = code
      .GlobalEnv$newMSL[srcNewSelect,".merge_linenum"] = linenum
      .GlobalEnv$newMSL[srcNewSelect,".merge_actionOrder"] = actionOrder
     
      ##### look for DEST #####
      
      foundMergeTarget=FALSE
      destNewTaxon=.GlobalEnv$newMSL %>% filter(name==as.character(destTaxonName))
      if(params$verbose) {cat(paste0("MERGE:             ",toupper(curChangeDf$rank)," code:",code," line:",linenum," dest.findTarget(new,",destTaxonName,")=",nrow(destNewTaxon),"\n"))}
      if(nrow(destNewTaxon)==1) {
        # great, found it
        foundMergeTarget = TRUE
      } 
      else if(nrow(destNewTaxon)==0) {
        # check if someone else already modified it
        destCurTaxon = curMSL %>% filter(name==as.character(destTaxonName))
        if(nrow(destCurTaxon)==0) {
          # just not there - we'll be creating it!
          foundMergeTarget = FALSE
        } 
        else if(nrow(destCurTaxon)==1 ) {
          # previous record - must have been modified already this MSL
          # find current matching record in newMSL
          destNewTaxon = (.GlobalEnv$newMSL$.prev_taxnode_id==destCurTaxon$taxnode_id)
          if(nrow(destNewTaxon)==1) {
            # great, found it
            foundMergeTarget = TRUE
            # warn about it (possibly contradictory) change
            log_change_error(curChangeDf, "WARNING", "MERGE.DEST_ALREADY_CHANGED", "Change=MERGE, but proposed taxon already modified", 
                             errorStr=paste0("proposal(s)=", destCurTaxon$out_filename, " also did a '",destCurTaxon$out_change,"' to '", destCurTaxon$out_target,"'",
                                             "; PROPOSED//OBSERVED=",diff_lineage(destLineage,destCurTaxon$lineage))
            )
          } 
          else if(nrow(destNewTaxon) == 0) { 
            # just not there - we'll be creating it!
            foundMergeTarget = FALSE
          } 
          else if(nrow(destNewTaxon)>1) { 
            # crazy - perhaps it was split
            log_change_error(curChangeDf, "ERROR", "MERGE.DEST_MANY", "Change=MERGE, but proposed taxon already split", 
                             errorStr=paste0("proposal(s)=", destNewTaxon$in_filename[1], " already did a SPLIT? ")
            )
            next;
          } # multiple matches
        } # dest found in curMSL
        else if(nrow(destCurTaxon)>1 ) {
          # crazy - perhaps it was split
          log_change_error(curChangeDf, "ERROR", "MERGE.DEST_MANY_2", "Change=MERGE, but multiple taxa", 
                           errorStr=paste0("proposal(s)=", destCurTaxon$out_filename[1], " already did a MERGED? ")
          )
          next;
        } # dest multiple in curMSL
      } # dest not in newMSL
      
      ##### rename/move or merge #####  
      if( foundMergeTarget )  {
        ###### . merge into existing ######
        if(params$verbose) {cat(paste0("MERGE_INTO:        ", srcTaxonName, " into ", destTaxonName," (",destNewTaxon$lineage,")\n"))}
        
        # set out_* in curMSL for SRC
        .GlobalEnv$curMSL[srcCurSelect,"out_updated"] = TRUE  # admin; mark this to save to db
        .GlobalEnv$curMSL[srcCurSelect,"out_change"] = "merge"
        .GlobalEnv$curMSL[srcCurSelect,"out_filename"] = proposalZip
        .GlobalEnv$curMSL[srcCurSelect,"out_target"] = destNewTaxon$lineage
        .GlobalEnv$curMSL[srcCurSelect,"out_notes"] = paste0("linenum=",linenum)
        .GlobalEnv$curMSL[srcCurSelect,".out_taxnode_id"] = destNewTaxon$taxnode_id
        .GlobalEnv$curMSL[srcCurSelect,".actionOrder"] = actionOrder
        
        # set admin columns in newMSL for DEST
        destNewSelect = newMSL$taxnode_id == destNewTaxon$taxnode_id
        .GlobalEnv$newMSL[destNewSelect,".actionOrder"] = actionOrder
        
        #
        # remove the now abandon record for src in newMSL (because it is linked
        # to dest, instead)
        #
        
        # check for kids
        srcNewKidsSelect=(.GlobalEnv$newMSL$parent_id==.GlobalEnv$newMSL[srcNewSelect,"taxnode_id"])
        if(sum(srcNewKidsSelect,na.rm=TRUE)>0) {
          log_change_error(curChangeDf, "ERROR", "MERGE_INTO.KIDS.UNIMP", "Change=MERGE, merging child taxa not implemented (yet) ", 
                           notes=paste0("taxon='", srcTaxonName, "', lineage='",srcLineage,"', kids: N=",sum(srcNewKidsSelect,na.rm=TRUE),
                                        ", NAMES=[", paste(.GlobalEnv$newMSL$rank[srcNewKidsSelect],.GlobalEnv$newMSL$name[srcNewKidsSelect],sep=":"),"]" ) 
          )
          next;
        }     
        
        # remove taxon from NEW
        .GlobalEnv$newMSL = .GlobalEnv$newMSL[!srcNewSelect,]
        
        # message setup
        logNotes =paste0(toupper(curChangeDf$.action)," ",
                         srcCurTaxon$rank, " named '", srcCurTaxon$name, "'", 
                         " into ", destNewTaxon$rank, " named '", destNewTaxon$name, "'", 
                         " CUR//PROPOSED=",  
                         diff_lineages(srcCurTaxon$lineage,destNewTaxon$lineage)
        )
        
        # we don't implement merging child nodes, so only SPECIES implemented
        if(srcCurTaxon$rank != "species") {
          log_change_error(curChangeDf, "ERROR",  "MERGE_INTO.SRC_NOT_SPECIES.UNIMP", 
                           errorStr=paste0("Change=MERGE only species merging has been implemented"),
                           notes=logNotes
          )
          
        }
        # success note
        log_change_error(curChangeDf, "SUCCESS",  "MERGE_INTO.OK", 
                         errorStr=paste0("Change=",toupper(curChangeDf$.action),", applied successfully"), 
                         notes=logNotes
        )
      } # MERGE_INTO
      else { # MERGE_MOVE_RENAME
        ###### . move/rename taxon ######
        if(params$verbose) {cat(paste0("MERGE_MOVE_RENAME: ", srcTaxonName, " into ", destTaxonName," (",destLineage,")\n"))}
        
        # find parent in newMSL
        # 
        # verify that PARENT taxon exists already in newMSL
        #
        parentDestNameAlias =""
        parentDestNewSelect = NA
        parentDestNewTaxon = NA
        parentDestNameAlias = ""
        if(is.na(destParentName)) {
          # no parent - use root node
          parentDestNewSelect=(.GlobalEnv$newMSL$taxnode_id==.GlobalEnv$newMSL$tree_id)
        } 
        else {
          # newMSL: find parent by name (must be unique) 
          parentDestNewSelect=(.GlobalEnv$newMSL$name==as.character(destParentName))
        }
        parentDestNewTaxon = .GlobalEnv$newMSL[parentDestNewSelect,]
  
        if(params$verbose) {cat(paste0("MERGE_MOVE_RENAME: ",toupper(curChangeDf$rank)," code:",code," line:",linenum," '",destTaxonName, "' findParent(",destParentName,parentDestNameAlias,")=",sum(parentDestNewSelect)),"\n")}
          
        if(sum(parentDestNewSelect)==0 ) {
          # find parent by other means: look for name in old (curMSL), and assume it got renamed already in newMSL
          # and find it's new entry in newMSL
          parentDestCurTaxon = curMSL %>% filter(name==as.character(destParentName))
          if(nrow(parentDestCurTaxon)==1) {
            # get what that older name became
            parentDestNewSelect = (.GlobalEnv$newMSL$.prev_taxnode_id==parentDestCurTaxon$taxnode_id)
            parentDestNewTaxon = .GlobalEnv$newMSL[parentDestNewSelect,]
            parentDestNameAlias =    paste0(" [AKA ",parentDestNewTaxon$name,"]")
            # INFO message
            log_change_error(curChangeDf, "INFO", "MERGE_MOVE_RENAME.PARENT_ALREADY_CHANGED", 
                             errorStr="PROPOSED parent taxon already modified", 
                             notes=paste0("proposal(s)=", parentDestNewTaxon$.otherLineageProposal, " did a '", parentDestNewTaxon$.otherLineageAction,"' ")
            )
          } 
          else if(nrow(parentDestCurTaxon)==0 ) {
            # changeDf: look for another proposal that already changed the name in newMSL?
            changedParent = .GlobalEnv$allChangeDf %>% filter(.srcTaxon==as.character(destParentName))
            if(nrow(changedParent)==1) {
              # get what that older name became
              parentDestNewTaxon  = (.GlobalEnv$newMSL %>% filter(name==changedParent$.destTaxon))
              
              # INFO message
              log_change_error(curChangeDf, "INFO", "MERGE_MOVE_RENAME.PARENT_ALREADY_CHANGED_2",
                               errorStr="PROPOSED parent taxon already modified", 
                               notes=paste0("proposal(s)=", parentDestNewTaxon$.otherLineageProposal, " did a '", parentDestNewTaxon$.otherLineageAction,"' ")
              )
            } 
            else {
              log_change_error(curChangeDf, "ERROR", "MERGE_MOVE_RENAME.PARENT_NOT_FOUND_2",
                               errorStr="PROPOSED parent taxon not found", 
                               notes=paste0("'",destParentName,"' : ",destParentLineage)
              )
              next
            }
          }
          else if(nrow(parentDestCurTaxon)>1) {
            log_change_error(curChangeDf, "ERROR", "MERGE_MOVE_RENAME.PARENT_MULTIPLE_MATCH", 
                             errorStr="Probably parent was split by another proposal", 
                             notes=paste0("proposal(s)=", changedParent$.code, " did a '",changedParent$change,"' ",
                                          "of '",changedParent$.srcTaxon, "' to '", changedParent$.destTaxon,"'")
            )
            next
          } 
        }
        else if(sum(parentDestNewSelect)==0) {
          # just missing
          log_change_error(curChangeDf, "ERROR", "MERGE_MOVE_RENAME.PARENT_NO_EXIST", 
                           errorStr=paste0("Change=",toupper(curChangeDf$.action),", but parent rank taxon does not exist"), 
                           notes=paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage)
          )
          next
        } # just missing
        
        ###### . create binomial check ######
        
        # if  SPECIES, check binomial prefix matches parent genus
        #
        if( destTaxonRank=="species" ) {
          # check if parent is a genus
          if(parentDestNewTaxon$rank=="genus") {
            destParentGenusTaxon = parentDestNewTaxon
          } 
          else {
            if(parentDestNewTaxon$rank == "subgenus" ) {
              # check up one rank
              destParentGenusTaxon = .GlobalEnv$newMSL %>% filter(taxnode_id == parentDestNewTaxon$parent_id)
              if(destParentGenusTaxon$rank != "genus") {
                # can't find genus parent
                log_change_error(curChangeDf, "ERROR", "MERGE_CREATE.SPECIES_SUBGENUS_NO_GENUS", 
                                 errorStr="Change=MERGE_CREATE, can not find parent genus", 
                                 notes=paste0("subgenus '",parentDestNewTaxon$name,"' is not in a genus;",
                                              " it's parent '",destParentGenusTaxon$name,"' is a ",destParentGenusTaxon$rank )
                )
                next;
              }
            } 
            else {
              # parent isn't subgenus or genus
              log_change_error(curChangeDf, "ERROR", "MERGE_CREATE.SPECIES_NO_GENUS", 
                               errorStr="Change=MERGE_CREATE, can not find parent genus", 
                               notes=paste0("parent '",destParentGenusTaxon$name,"' is a ",destParentGenusTaxon$rank )
              )
              next;
            }
          } 
          # have genus (parent or grandparent)
          # check binomial naming
          if( str_detect(destTaxonName,paste0(destParentGenusTaxon$name," ")) != TRUE ) {
            errLevel = ifelse(params$processing_mode %in% c("draft","final"),"ERROR","WARNING")
            log_change_error(curChangeDf, errLevel, "CREATE_RENAME.SPECIES_BINOMIAL_MISMATCH", 
                             errorStr="Change=RENAME, but proposed species names does not start with 'genus[space]' per binomial naming convention", 
                             notes=paste0("GENUS//SPECIES =", diff_strings(destParentGenusTaxon$name,destTaxonName))
            )
            if( errLevel=="ERROR") {next}
          }
        } # species (check binomial)
        
        # add to list of genera to check for binomial problems after proposal is complete
        if(  destTaxonRank == "genus" ) {
          renamedGenera = rbind(renamedGenera, 
                                data.frame(name =c(destTaxonName$genus), 
                                           taxnode_id=.GlobalEnv$newMSL$taxnode_id[srcNewTarget],
                                           code=code, 
                                           linenum= linenum,
                                           actionOrder = actionOrder)
          )
        } 
        
        # check if that is the expected lineage
        # RENAME: WARN that PARENT_LINEAGE is not expected, AND USE MSL/OBSERVED PARENT LINEAGE
        # DON'T NEED THIS - would already be reported above as SRC. LINEAGE_CHANGED
        # if( curChangeDf$.destParentName != .GlobalEnv$newMSL[srcNewParent,"name"] ) {
        #   parentRank = .GlobalEnv$newMSL[srcNewParent,]$rank
        #   errorDf=log_error(curChangeDf$.code,curChangeDf$.linenum,curChangeDf$.order,curChangeDf$change,curChangeDf$rank,curChangeDf$.changeTaxon,
        #                    "WARNING", "RENAME.PARENT_NAME", 
        #                    "parent name not as expected; was this meant to be a move? Is there a typo?", 
        #                    paste0("parent in this proposal: '", curChangeDf$.destParentName,"'[", curChangeDf$.destParentRank, "], ",
        #                          "but MSL parent is '",.GlobalEnv$newMSL[srcNewParent,"name"],"'[", .GlobalEnv$newMSL[srcNewParent,]$rank, "]; ",
        #                          "other new proposals: ",.GlobalEnv$newMSL[srcNewTarget,"prev_proposals"])
        #   )
        # }
        
        # append proposal
        destNewSelect = srcNewSelect
        
        .GlobalEnv$newMSL[destNewSelect,"prev_proposals"] = paste0(
          ifelse(is.na(.GlobalEnv$newMSL[destNewSelect,"prev_proposals"]),"",paste0(.GlobalEnv$newMSL[destNewSelect,"prev_proposals"],",")),
          proposalZip,":",linenum)
        
        # update (rename/move) dest in newMSL
        .GlobalEnv$newMSL[destNewSelect,"name"] = destTaxonName
        .GlobalEnv$newMSL[destNewSelect,"lineage"] = paste0(parentDestNewTaxon$lineage,";",destTaxonName)
        .GlobalEnv$newMSL[destNewSelect,".actionOrder"] = actionOrder
        
        # put out_* changes on curMSL
        .GlobalEnv$curMSL[srcCurSelect,"out_updated"] = TRUE  # admin; mark this to save to db
        .GlobalEnv$curMSL[srcCurSelect,"out_change"] = "merge"
        .GlobalEnv$curMSL[srcCurSelect,"out_filename"] = proposalZip
        .GlobalEnv$curMSL[srcCurSelect,"out_target"] = newMSL[destNewSelect,"lineage"]
        .GlobalEnv$curMSL[srcCurSelect,"out_notes"] = paste0("linenum=",linenum)
        .GlobalEnv$curMSL[srcCurSelect,".out_taxnode_id"] = newMSL[destNewSelect,"taxnode_id"]
        .GlobalEnv$curMSL[srcCurSelect,".actionOrder"] = actionOrder
        
        # success note
        errorDf=log_change_error(curChangeDf, "SUCCESS", "MERGE_MOVE_RENAME.OK", "Change=MERGE, applied successfully", 
                                 paste0("MERGE ",srcTaxonRank, 
                                        " from ", "'",srcTaxonName,"'",
                                        " to ",   "'", newMSL[destNewSelect,"name"], "'",
                                        " in ",   "'", parentDestNewTaxon$lineage, "'",
                                        ", CURRENT//PROPOSED=", diff_lineages(curMSL[srcCurSelect,"lineage"], newMSL[destNewSelect,"lineage"] )
                                 )
        )
        
      } # MERGE_MOVE_RENAME
    } # MERGE
    else {
      # ........................................................................
      #
      ##### UNKNOWN ACTION not implemented. #####
      #
      # ........................................................................
      log_change_error(curChangeDf, "ERROR", "CHANGE.UNIMP", 
                       errorStr=paste0("Change=",toupper(action)," is NOT (yet) implemented"), 
                       notes=paste0("lineageChange=", diffLineageString(srcLineage,destLineage))
      )
      next;
      
    } # ACTION UN-IMP
    
  } # for each line in XLSX
  
  # ........................................................................
  #
  #### POST-QC #### 
  #
  # ........................................................................
  .GlobalEnv$actionOrder = .GlobalEnv$actionOrder+1
  
  #
  ##### empty taxa #####
  # did proposal create empty (non-species) taxa? (global scan)
  #
  kidCounts = apply(.GlobalEnv$newMSL[,"taxnode_id"],1,
                    function(parent_id,MSL) { sum(MSL$parent_id==parent_id)},
                    MSL=.GlobalEnv$newMSL)
  # scan for taxa with no kids, not species, not yet reported.
  # convert to data.frame, as the loop will work nicely
  # which is not true of the default, un-key'ed data.table.
  emptyTaxa = as.data.frame(
      .GlobalEnv$newMSL[kidCounts == 0,] 
      %>% filter(level_id != 600) 
      %>% filter(is.na(.emptyReported))
      )
  for( emptyTaxaRow in rownames(emptyTaxa) ) {
    # get a single taxon
    emptyTaxon = emptyTaxa[emptyTaxaRow,]
    
    # default error report
    # (need to improve change tracking to do better!)
    # NOTE we don't even try to track species moved out of genera
    errLineNum     = ""
    errChange      = ""
    errText        = "check if lower rank taxa were moved elsewhere?"
    # split
    if ( !is.na(emptyTaxon$.split_linenum) ) { 
      errLineNum=emptyTaxon$.split_linenum
      errChange = "Split"
      errText = paste0(proposalsDf[emptyTaxon$.split_code,]$basename,":",emptyTaxon$.split_linenum," did a split")
    } else if( !is.na(emptyTaxon$prev_proposals) ) {
      #errLineNum=""
      errChange="previous"
      errText = emptyTaxon$prev_proposals
    } else if( !is.na(emptyTaxon$.otherLineageProposal) ) {
      #errLineNum=emptyTaxon$.otherLineageLineNum
      errChange ="previous"
      errText = paste0(emptyTaxon$.otherLineageProposal," ",emptyTaxon$.otherLineageAction)
    }
    # mostly from create - ops we need code, not filename :-( )
    #if( !is.na(emptyTaxon$in_filename) ) {
    #  errLinenum = emptyTaxon$in_notes
    #  err
    #}
    log_error(code,linenum=errLineNum,action=errChange,actionOrder=actionOrder,
              rank=emptyTaxon$rank,taxon=emptyTaxon$name,
              levelStr="ERROR", errorCode="PROPOSAL.EMPTY_TAXA", 
              errorStr=paste0("Proposal created empty (non-species) taxa"), 
              notes=paste0("the ",emptyTaxon$rank," '",emptyTaxon$name,"' is empty - it does not contain any lower rank taxons; ", 
                           errText)
    )
    # mark empty taxa so we don't re-report them
    .GlobalEnv$newMSL[.GlobalEnv$newMSL$name %in% emptyTaxon$name,".emptyReported"] = proposalsDf[code,"xlsx"]
  }
  
  ##### renamed genera - binomial #####
  # skip initial row - just a placeholder for dataframe structure
  for(  genusCheckRow in rownames(renamedGenera[-1,]) )  {
    # count these as actions
    .GlobalEnv$actionOrder = .GlobalEnv$actionOrder+1
    
    # get genus 
    renamedGenusTaxon = .GlobalEnv$newMSL %>% filter(name == renamedGenera[genusCheckRow,"name"] )
    # subgenera
    renamedGeneraSubgeneraTaxons = .GlobalEnv$newMSL %>% 
      filter(parent_id %in% c(renamedGenusTaxon$taxnode_id)) %>%
      filter(level_id == 550 )
    # species
    renamedGeneraSpeciesTaxons = .GlobalEnv$newMSL %>% 
      filter(parent_id %in% c(renamedGenusTaxon$taxnode_id,renamedGeneraSubgeneraTaxons$taxnode_id)) %>%
      filter(level_id==600)
    
    # check names
    binomialViolations =str_detect(renamedGeneraSpeciesTaxons$name,paste0("^",renamedGenusTaxon$name," "), negate=TRUE)
    if( sum(binomialViolations) > 0 ) {
      # get offending species list
      speciesDf = as.data.frame(renamedGeneraSpeciesTaxons)[binomialViolations,]
      speciesDf = speciesDf[order(speciesDf$name),]
      log_error(code=renamedGenera[genusCheckRow,"code"],linenum=renamedGenera[genusCheckRow,"linenum"],
                action="rename_genus",actionOrder=renamedGenera[genusCheckRow,"actionOrder"],
                rank="species",taxon=speciesDf$name,
                # override default order
                levelStr="ERROR", errorCode="RENAME_GENUS.SPECIES_BINOMIAL_MISMATCH", 
                errorStr="Change=RENAME_GENUS, but species name does not start with 'genus[space]' per binomial naming convention", 
                notes=paste0("new genus name: ", renamedGenusTaxon$name )
      )
    }
      
  }
  ##### split taxa - delete #####
  #
  # if a taxon is split, and no split keeps it's original name
  # then quietly delete the original name from the new MSL
  # (it would have been created when the curMSL was copied to create
  # the new one)
  splitDeleteIdx = .GlobalEnv$newMSL$.split & !.GlobalEnv$newMSL$.split_kept
  splitKeepIdx   = .GlobalEnv$newMSL$.split & .GlobalEnv$newMSL$.split_kept
  if( sum(splitDeleteIdx) > 0 ) {
    # count these as actions
    .GlobalEnv$actionOrder = .GlobalEnv$actionOrder+1

        # log that we're removing them
    # function(errorDf,code,row,change,rank,taxon,levelStr,errorCode,errorStr,notes
    log_error(code=.GlobalEnv$newMSL$.split_code[splitDeleteIdx],linenum=.GlobalEnv$newMSL$.split_linenum[splitDeleteIdx],
              action="split_abolish",actionOrder=actionOrder,
              rank=as.character(.GlobalEnv$newMSL$rank[splitDeleteIdx]), 
              taxon=.GlobalEnv$newMSL$name[splitDeleteIdx],
              levelStr="INFO", errorCode="SPLIT.IMPLICIT_ABOLISH", 
              errorStr="Change=SPLIT, but no split line kept original name, so remove original name from new MSL", 
              notes=paste0("removed ",.GlobalEnv$newMSL$rank[splitDeleteIdx]," ",.GlobalEnv$newMSL$name[splitDeleteIdx],
                           " from ",.GlobalEnv$newMSL$lineage[splitDeleteIdx] )
    )
    
    # check for kids - work back to front for rowIdx stability across row deletions
    for(rowIdx in rev(which(splitDeleteIdx)) ) {
      # count these as actions
      .GlobalEnv$actionOrder = .GlobalEnv$actionOrder+1
      
      # get kid count for that taxon
      srcKids=(.GlobalEnv$newMSL$parent_id==.GlobalEnv$newMSL$taxnode_id[rowIdx])
      
      if(sum(srcKids,na.rm=TRUE)>0) {
        # error: can't abolish something with kids
        log_error(code=.GlobalEnv$newMSL[rowIdx,]$.split_code,
                  linenum=.GlobalEnv$newMSL[rowIdx,]$.split_linenum,
                  action="split_abolish",actionOrder=actionOrder,
                  rank=as.character(.GlobalEnv$newMSL[rowIdx,]$rank),.GlobalEnv$newMSL[rowIdx,]$name,
                  levelStr= "ERROR", errorCode="SPLIT.IMPLICIT_ABOLISH_WITH_KIDS", 
                  errorStr="Change=ABOLISH, taxon still has un-abolished/moved children", 
                  notes=paste0("taxon=", .GlobalEnv$newMSL[rowIdx,]$name, ", lineage=",.GlobalEnv$newMSL[rowIdx,]$lineage,", kids: N=",sum(srcKids,na.rm=TRUE),
                               ", NAMES=[",paste(.GlobalEnv$newMSL$rank[srcKids],.GlobalEnv$newMSL$name[srcKids],sep=":"),"]"
                  ) 
        )
      } else {
        # no kids, nuke it
        # remove the rows from newMSL
        .GlobalEnv$newMSL = .GlobalEnv$newMSL[-rowIdx,]
        # and from the keep list
        splitKeepIdx=splitKeepIdx[-rowIdx]
      }
    } # for taxon to split_abolish (implicit)
  }
    
  if( sum(splitKeepIdx) > 0 && params$verbose ) {
    # log that we're NOT removing them
    # function(errorDf,code,row,change,rank,taxon,levelStr,errorCode,errorStr,notes
    log_error(code=.GlobalEnv$newMSL$.split_code[splitKeepIdx],
              linenum=.GlobalEnv$newMSL$.split_linenum[splitKeepIdx],
              action="split_keep",actionOrder=actionOrder,
              rank=.GlobalEnv$newMSL$rank[splitKeepIdx], 
              taxon=.GlobalEnv$newMSL$name[splitKeepIdx],
              levelStr="INFO",  errorCode="SPLIT.KEEP", 
              errorStr="Change=SPLIT, but one split directive kept the original name", 
              notes=paste0("split and keep ",.GlobalEnv$newMSL$rank[splitKeepIdx]," '",.GlobalEnv$newMSL$name[splitKeepIdx],
                           "' from ",.GlobalEnv$newMSL$lineage[splitKeepIdx] )
    )
  }
  return(list(errorDf=errorDf))
} # apply_changes()
# ```
# 
# # process changes
# 
# RE-RUN FROM HERE!
# 
# ```{r process_changes}
#### SECTION apply changes #####
#
# GLOBAL VARS
#
#
# debug 
#cat("allErrorDf:",tracemem(allErrorDf),"\n");
if( params$tmi ) {
    cat("Summary of allErrorDf$code:\n" )
    summary(as.factor(.GlobalEnv$allErrorDf$code))
}

if( !is.null(nrow(allChangeDf)) && nrow(allChangeDf) > 0 ) {
  
  
  # Apply the Change
  apply_changes(allChangeDf[allChangeDfOrder,])
  
  # Internal sanity check
  if(!(20070000 %in% .GlobalEnv$newMSL$ictv_id)) { 
    cat("!!!!!!!!!!!!!!!! LOST ROOT NODE in ",code," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    if(interactive()) {browser()}
  }
  #
  # report on error counts from apply
  # (row=0 is load)
  applyErrorDf = .GlobalEnv$allErrorDf %>% filter(row > 0)
  # per proposal
  for( code in levels(as.factor(allChangeDf$.code)) ) {
    changeMatches = allChangeDf$.code == code
    errorMatches  = applyErrorDf$code == code
    levelSummary = summary(applyErrorDf[errorMatches,]$level)
    levelSummary = levelSummary[levelSummary>0]
    cat("# APPLIED:", sprintf("%3d",sum(changeMatches)),"changes from",proposalsDf[code,"basename"]," with",paste(names(levelSummary),levelSummary,sep=":"),"\n")
    
  }
  # global summary
  levelSummary = summary(.GlobalEnv$allErrorDf$level)
  levelSummary = levelSummary[levelSummary>0]
  cat("# TOTAL  :", sprintf("%3d",nrow(allChangeDf)),"changes ",
      "from",length(levels(as.factor(allChangeDf$.code))),"proposals ",
      "with",paste(names(levelSummary),levelSummary,sep=":"),"\n")

  
} # > 0 proposals to process

#### SECTION write error list #####
#
# set "final" to be TRUE (which produces per-subcommittee error files) only if not in "validate" mode
#print(dim(.GlobalEnv$allErrorDf))
#print(dim(.GlobalEnv$curMSL))
#print(dim(.GlobalEnv$newMSL))

#print("write_error_summary(final)")
write_error_summary(.GlobalEnv$allErrorDf,params$processing_mode != "validate")

#### EXPORT MSL ####
if(params$export_msl) {
  # 
  # # Export MSL to a TSV that can be easily diff'ed
  # - no taxnode_ids
  # - no dates
  # - no filenames
  #
  # 
  # ```{r tsv_export}
  #### SECTION export to TSV #####
  tsvColList = c(#"taxnode_id",
    #"parent_id",
    #"tree_id",
    "msl_release_num",
    "level_id",
    "name",
    #"ictv_id",
    "molecule_id",
    "abbrev_csv",
    "genbank_accession_csv",
    "genbank_refseq_accession_csv",
    "refseq_accession_csv",
    "isolate_csv",
    "notes",
    # historic columns we don't update
    #"is_ref",
    #"is_official",
    #"is_hidden",
    #"is_deleted",
    #"is_deleted_next_year",
    #"is_typo",
    #"is_renamed_next_year",
    #"is_obsolete",
    "in_change",
    "in_target",
    "in_filename",  # will be stripped to just  code
    "in_notes",
    "out_change",
    "out_target",
    "out_filename", # will be stripped to just  code
    "out_notes",
    #"lineage", # computed by trigger in db
    #"cleaned_name", # computed by trigger in db
    #"rank", # should be in level_id
    # "molecule", # should be in molecule_id
    # program admin columns - not in db
    #"out_updated",
    #".prev_taxnode_id",
    #"prev_proposals",
    "host_source", 
    "exemplar_name",
    "genome_coverage",
    "notes" 
    #"lineage"  # computed by trigger in db
  )
  #prepare data for export
  tsvDf = data.frame(newMSL[,..tsvColList])
  tsvDf$in_filename = gsub("^([0-9]+\\.[0-9A-Z]+)\\..*", "\\1...",newMSL$in_filename)
  tsvDf$out_filename = gsub("^([0-9]+\\.[0-9A-Z]+)\\..*", "\\1...",newMSL$out_filename)
  
  newTsvFilename = file.path(params$out_dir,"msl.tsv")
  tsvout=file(newTsvFilename,"wt")
  cat("Writing ",newTsvFilename,"\n")
  # header
  cat(paste0(
    tsvColList,
    collapse="\t"
  ),
  "\n",
  file=tsvout)
  # body
  for(i in seq(1,nrow(tsvDf))) {
    cat(paste0(
      tsvDf[i,tsvColList],
      collapse="\t"
    ),
    "\n",
    file=tsvout)
  }
  close(tsvout)
  cat(paste0("WROTE   ", newTsvFilename, " (",nrow(tsvDf)," rows)\n"))
  
  # ........................................................................
  # 
  # # Export proposal summary to TSV
  #
  # metadata from DOCX
  #
  #### SECTION export DOCX meta to TSV #####
  tsvDocxColList = c(
    "subcommittee",
    "code",
    "docx",
    "xlsx",
    "title",
    "authorsEmails",
    "correspondingAuthor",
    "abstract"
  )
  #
  #prepare data for export
  # add missing columns (happens when there are no docx file to parse)
  for(col in tsvDocxColList) {
    if(!(col %in% names(proposalsDf))) {
      proposalsDf[,col] = ""
    }
  }
  #
  # subset columns
  tsvDocxDf = data.frame(proposalsDf[,tsvDocxColList])
  
  docxTsvFilename = file.path(params$out_dir,params$proposals_meta)
  docxMetaOut=file(docxTsvFilename,"wt")
  cat("Writing ",docxTsvFilename,"\n")
  # header
  cat(paste0(
    tsvDocxColList,
    collapse="\t"
  ),
  "\n",
  file=docxMetaOut)
  # body
  for(i in seq(1,nrow(tsvDocxDf))) {
    cat(paste0(
      tsvDocxDf[i,tsvDocxColList],
      collapse="\t"
    ),
    "\n",
    file=docxMetaOut)
  }
  close(docxMetaOut)
  cat(paste0("WROTE   ", docxTsvFilename," (",nrow(tsvDocxDf), ")\n"))
  
  #### SECTION export to SQL #####
  # ........................................................................
  #
  # SQL to insert new MSL
  #
  # column names and data codings in newMSL should match the taxonomy_node table. 
  #
  # convert newMSL to a data.frame, so we can convert factor columns to character, 
  # then generate SQL
  #
  # remember to 
  #   1. don't quote NULLs (is.na)
  #   2. escape apostrophies
  #   3. single-quot evalues
  #
  # ........................................................................
  # cat(paste0('"',paste(names(newMSL),collapse='",\n"'),'"'))
  sqlColList = c("taxnode_id",
              "parent_id",
              "tree_id",
              "msl_release_num",
              "level_id",
              "name",
              "ictv_id",
              "molecule_id",
              "abbrev_csv",
              "genbank_accession_csv",
              "genbank_refseq_accession_csv",
              "refseq_accession_csv",
              "isolate_csv",
              "notes",
              # historic columns we don't update
              #"is_ref",
              #"is_official",
              "is_hidden",
              #"is_deleted",
              #"is_deleted_next_year",
              #"is_typo",
              #"is_renamed_next_year",
              #"is_obsolete",
              "in_change",
              "in_target",
              "in_filename",
              "in_notes",
              "out_change",
              "out_target",
              "out_filename",
              "out_notes",
              #"lineage", # computed by trigger in db
              #"cleaned_name", # computed by trigger in db
              #"rank", # should be in level_id
              # "molecule", # should be in molecule_id
              # program admin columns - not in db
              #"out_updated",
              #".prev_taxnode_id",
              #"prev_proposals",
              "host_source",
              "exemplar_name",
              "genome_coverage"
              #,"notes" 
              #"lineage"  # computed by trigger in db
              )
  # get MSL num
  sql_msl_num = .GlobalEnv$newMSL$msl_release_num[1]
  
  newSqlFilename = file.path(params$out_dir,params$sql_load_filename)
  sqlout=file(newSqlFilename,"wt",encoding = "UTF-8")
  cat("Writing ",newSqlFilename,"\n")
  
  #
  # output start transaction
  # 
  cat("-- begin transaction\n", file=sqlout)
  cat("-- rollback transaction\n", file=sqlout)
  
  # ........................................................................
  #
  ##### SQL to add new MSL to [taxonomy_toc]   ##### 
  #
  # ........................................................................
  cat("insert into [taxonomy_toc] ([tree_id],[msl_release_num],[comments]) ",
      "values (", paste(
       sort(levels(as.factor(newMSL$tree_id)),decreasing = T)[1],
        sql_msl_num,
        "NULL",
      sep=","),
      ")\n",
      file=sqlout
  )
  
  # ........................................................................
  #
  ##### SQL to insert newMSL into [taxonomy_node] #####
  #
  # convert factors to character
  # (was easier to do on a data.frame, because of data.table FAQ 1.1)
  #
  newMslStr = as.data.frame(newMSL)
  for( col in sqlColList)  {
    if( class(newMslStr[,col]) == 'factor' ) {
      cat("factor: ",col, "\n")
    }
    newMslStr[,col] = as.character(newMslStr[,col])
  }
  
  #
  # add taxa to [taxonomy_node]
  #
  # output an insert statement for each row
  #
  rowCount = 1
  for(  row in order(newMSL$level_id) )  {
    #row=head(order(newMSL$level_id),n=1) # debug
    
    # insert several rows per batch insert  statement
    # much faster - fewer trigger calls
    # but makes localizing errors harder
    if( rowCount %% as.integer(params$sql_insert_batch_size) == 1 ) {
      cat(paste0("insert into [taxonomy_node] ",
                 "([",
                 paste0(sqlColList,collapse="],["),
                 "])",
                 "\n",
                 " values ",
                 "\n"
      ),file=sqlout)
    } else {
      # separate each value "(x,y,..)" in the batch by a comma
      cat(",",file=sqlout)
    }
    #
    # actual row values
    #
    cat(paste0("(",
               paste0(
                 # convert NA to NULL (not 'NA')
                 ifelse(is.na(newMslStr[row, sqlColList]), "NULL",
                        paste0("'",
                               # escape apostrophies as double-appostrophies (MSSQL)
                               gsub(
                                 "'", "''", newMslStr[row, sqlColList]
                               )
                               , "'")),
                 collapse = ","
               )
               , ")"),
        " -- lineage=",
        as.character(newMslStr[row, "lineage"]) ,
        "\n",
        file = sqlout
      )
    #
    # count rows
    # 
    rowCount=rowCount+1
  }
  
  # QC queries
  cat("-- QC queries \n", file=sqlout)
  cat(paste("
  select level_id,
      new_ct=count(case when in_change like 'new' then 1 end), 
    	other_ct=count(case when in_change<>'new' then 1 end),
    	null_ct=count(case when in_change is null then 1 end),
    	total_ct=count(*)
  from taxonomy_node
  where msl_release_num = ",sql_msl_num,"
  group by level_id
  order by level_id
  "), file=sqlout)
  
  # QC query
  cat(paste("-- QC Query
  select 
       level_id,in_change, ct=count(*)
  from taxonomy_node
  where msl_release_num = ",sql_msl_num,"
  group by level_id,in_change
  order by level_id,in_change
  
  "), file=sqlout)
  # ........................................................................
  #
  ##### SQL to update prevMSL.out in [taxonomy_node]   ##### 
  #
  # column names and data codings in newMSL should match the taxonomy_node table. 
  #
  # convert newMSL to a data.frame, so we can convert factor columns to character, 
  # then generate SQL
  #
  # ........................................................................
  
  cat("Writing out_* updates for prevMSL\n")
  #
  # only update these columns
  #
  curMslColList = c("out_change","out_target","out_filename","out_notes")
  #
  # convert factors to character
  # (was easier to do on a data.frame, because of data.table FAQ 1.1)
  #
  curMslStr = as.data.frame(curMSL %>% filter(!is.na(out_change)) )
  for( col in curMslColList)  {
    if( class(curMslStr[,col]) == 'factor' ) {
      cat("factor: ",col, "\n")
    }
    curMslStr[,col] = as.character(curMslStr[,col])
  }
  
  #
  # output start transaction
  # 
  cat("-- begin transaction\n", file=sqlout)
  cat("-- rollback transaction\n", file=sqlout)
  
  #
  # output an insert statement for each row
  #
  for(  row in rownames(curMslStr) )  {
    #row=head(rownames(curMslStr),n=1) # debug
    cat(paste0("update [taxonomy_node] set ",
               
               paste0(
                 paste0(
                   # column names
                   "[",curMslColList,"]=",
                   # convert NA to NULL (not 'NA')
                   ifelse(is.na(curMslStr[row,curMslColList]),"NULL",
                          paste0("'",
                                 # escape apostrophes as double-apostrophes (MSSQL)
                                 gsub("'","''", curMslStr[row,curMslColList])
                                 ,"'"))
                 ,collapse=","),
               " where [taxnode_id]=",
               curMslStr[row,"taxnode_id"]
               ) 
    ),
    "\n",
    file=sqlout)
  }
  # QC SQL
  cat("\n-- build deltas (~7min) \n
  EXEC rebuild_delta_nodes_2 NULL
  
  -- rebuild merge/split (seconds) \n
  EXEC [dbo].[rebuild_node_merge_split] 
  
  ", file=sqlout)
  cat(paste("
  -- NOW check if all newMSL have delta in, and prevMSL have delta out
  select 'prevMSL w/o delta to new', count(*) from taxonomy_node
  where msl_release_num = ",sql_msl_num-1,"
  and taxnode_id not in (select prev_taxid from taxonomy_node_delta)
  union all
  select  'newMSL w/o delta to prev',  count(*) from taxonomy_node
  where msl_release_num = ",sql_msl_num,"
  and taxnode_id not in (select new_taxid from taxonomy_node_delta)
  "),
  file=sqlout)
  
  # QC query
  cat(paste("-- QC Query
  select 
       level_id,out_change, ct=count(*)
  from taxonomy_node
  where msl_release_num = ",sql_msl_num-1,"
  group by level_id,out_change
  order by level_id,out_change
  
  "), file=sqlout)
  #
  # close/flush file
  #
  close(sqlout)
  cat("WROTE   ", newSqlFilename, "\n")
}

if( FALSE ) {
  #
  # stats to match in db
  #
  summary(as.factor(newMSL$in_change))
  # new      split  NA's 
  # 1043     0      13653 
  summary(as.factor(curMSL$out_change))
  #abolish  demote   merge    move promote  rename    type    NA's 
  #     20       0       0      20       0    1662       0   11971
  
  # summarize taxa by rank/change in new MSL
  data.frame(in_change=summary(as.factor(paste0(newMSL$rank,".",newMSL$in_change))))
  data.frame(in_change=summary(as.factor(paste0(newMSL$in_change))))
  
  # summarize taxa by rank/change in new MSL
  data.frame(out_change=summary(as.factor(paste0(curMSL$rank,".",curMSL$out_change))))
  data.frame(out_change=summary(as.factor(paste0(curMSL$out_change))))
}
if( FALSE ) {
  rdataFilename = paste0(params$proposals_dir,"/.RData")
  save.image(file=rdataFilename)
  cat("WROTE",rdataFilename,"\n")
}

cat("# COMPLETED.\n")
