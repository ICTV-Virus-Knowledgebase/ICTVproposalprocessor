#!/usr/bin/env Rscript
#
# TODO: 
#  1. add arg parsing
#  2. move from params to opt? or keep params for Rmd compatibility?
#  3. make all filenames into dir+filename
#  4. update VMR - to table dump? Factor our accessions? 
#  5. add 2023v1 template support
# 
# Validate ICTV proposal(s).xlsx
#
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
  make_option(c("--msl"), action="store_true", default=FALSE, dest="export_msl", 
              help="Write resulting MSL to load_msl.sql and msl.tsv [default \"%default\"]"),
  
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
  make_option(c("-r","--refDir"), default="current_msl", dest="ref_dir", 
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
  make_option(c("--cvTemplate"), default="TP_Template_Excel_module_2023_v1.xlsx", dest="template_xlsx_fname",
              help="Template proposal xlsx, used to load CVs [default \"%default\"]"),
  make_option(c("--cvTemplateSheet"), default="Menu Items (Do not change)", dest="template_xlsx_sheet",
              help="Template proposal xlsx, used to load CVs [default \"%default\"]"),
  make_option(c("--vmr"), default="VMR_MSL38_v1.xlsx", dest="vmr_fname",
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
  # defeat auto-caching when debugging
  #rm(docxList,xlsxList,changeList)
  params$verbose = T
  params$tmi = T
  params$debug_on_error = T
  params$mode = 'draft'
  params$export_msl = T
  params$proposals_dir = "./MSL39dbg"
  params$out_dir       = "./MSL39dbg_results"
  #params$proposals_dir = "EC55"
  #params$out_dir       = "EC55_results"
  params$qc_regression_tsv_fname = "QC.regression.new.tsv"
  cat("!! VERBOSE: ",params$versbose, "\n")
  cat("!! TMI:     ",params$tmi, "\n")
  cat("!! MODE:    ",params$mode, "\n")
  cat("!! SRC_DIR: ",params$proposals_dir, "\n")
  cat("!! OUT_DIR: ",params$out_dir, "\n")
}
#
# args trickery - TMI implies verbose
#
if(params$tmi) { params$verbose=T}

#
# WARNING: we use data.TABLE instead of data.FRAME
#
# this allows modification in place of a (data) passed
# to a subroutine (pass-by-reference feature)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(writexl) )# another option library(openxlsx)
suppressPackageStartupMessages(library(DescTools)) # AscToChar
suppressPackageStartupMessages(library(qdapTools)) # read docx

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
  errorSortCols = errorDf[,c("code","row")]
  errorSortCols$row_n = as.integer(errorSortCols$row)
  errorsSorted = do.call(order,errorSortCols[,c("code","row_n")])
  
  #### pretty format ####
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
  
  #### write files ####
  # XLS version (for end human user)
  fname = file.path(params$out_dir,params$qc_summary_fname)
  write_xlsx( x=as.data.frame(prettyErrorDf)[,prettyCols][,prettyCols],path=fname)
  if(params$verbose || params$tmi) {cat(paste0("Wrote: ", fname, " (",nrow(errorDf)," rows)\n"))}
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

#### load cache (ref/.Rdata) ####
if(params$verbose){ cat(paste0("REF_DIR:        ",params$ref_dir,"/\n"))}
cacheFilename=paste0(params$ref_dir,"/",params$cache_fname)
if(params$use_cache && !params$update_cache && file.exists(cacheFilename)) {
  if(params$verbose){ cat("Loading image ", cacheFilename, "...\n")}
  load(file=cacheFilename)
  #rm("changeList") # to re-run only QC, delete this cache
  #cat("RM(changeList) # re-run QC")
} else {
  if(params$verbose){ cat("SKIP: loading cache file ", cacheFilename, "\n")}
  
  
  # 
  # final vs draft proposal filename patterns.
  #
  if( params$processing_mode == "final") {
    filenameFormatRegex="^[0-9][0-9][0-9][0-9]\\.[0-9][0-9][0-9][A-Z]\\.[A-Za-z]+\\.[^ ]*"
    filenameFormatMsg="####[A-Z].###[A-Z].[A-Z]+.____"
  } else if( params$processing_mode == "draft") {
    filenameFormatRegex="^[0-9][0-9][0-9][0-9]\\.[0-9][0-9][0-9][A-Z]\\.[A-Za-z]+\\.v[0-9]+\\.[^ ]*"
    filenameFormatMsg="####[A-Z].###[A-Z].[A-Z]+.v#.____"
  } else if( params$processing_mode == "validate") {
    # allow any doc/xls file
    filenameFormatRegex="^.*"
    filenameFormatMsg="*.____"
  } else {
    # unsupported format
    cat(paste0("ERROR: --mode='",params$processing_mode,"' is not a valid option: validate, draft, or final\n"))
    quit(save="no", status=1)
  }
  if(params$verbose){ cat("(PROCESSING) MODE      :",params$processing_mode,"\n") }
  
  # --------------------------------------------------------
  # << NEXT >>
  # --------------------------------------------------------
  # * implement MOVE / re_lineage_tree()
  #      * test: 2022.001B lines 5:6
  #   * RE_ORG to separate scan/load from process
  #      * make it easy to re-process - move globals to process section
  #   * check for UN-IMP actions
  #     * split
  #     * promote
  #     * demote
  #     * merge
  #   * SQL export
  #   * error/warning report by directory
  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  # CURRENT STATUS:
  # summary(as.factor(paste(allErrorDf$level, allErrorDf$error)))
  #              ERROR ACTION.UNK          ERROR CREATE.DUP_ACC  ERROR CREATE.PARENT_NO_EXIST 
  #                            21                            21                            29 
  #            ERROR CREATE.W_SRC             ERROR DEST.IN_CUR             ERROR DEST.IN_NEW 
  #                            52                            14                             6 
  #         ERROR RENAME.NO_EXIST           ERROR RENAME.WO_SRC WARNING CREATE.PARENT_LINEAGE 
  #                            25                             1                           214 
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
  
  #```{r setup, include=FALSE}
  
  
  # rm("xlxsList","changeList") # to re-run xlsx file loading & QC, delete these caches
  # rm("changeList") # to re-run only QC, delete this cache
  
  #```
  
  
  # function to copy prev MSL to new MSL
  
  #```{r function_copy_prev_msl_to_new_msl, echo=FALSE}
  
  #
  # this uses the DATABASE schema for the taxonomy_node table (ie, naming convention)
  #
  # extra admin fields added, which wont be saved
  #   prev_taxnode_id
  #   prev_proposals
  #
  createNewMSL = function(curMSL,prev_msl,dest_msl,taxnode_delta) {
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
    newMSL[,"prev_taxnode_id"] = newMSL[,"taxnode_id"]
    
    # CSV of proposals that have already modified this taxon in this MSL
    newMSL[,"prev_proposals"] = NA_character_
    
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
    
    # rename release
    newMSL[newMSL$rank=="tree","name"]  = params$msl_name
    newMSL[newMSL$rank=="tree","notes"] = params$msl_notes
    
    # add the new rows
    #return(rbind(taxonomy_node, newMSL))
    return(newMSL)
  }
  
  #
  # get next available taxnode_id in a given MSL
  #
  assignNextTaxnodeId = function(taxonomy_node,msl) {
    return(max(subset(taxonomy_node,msl_release_num==msl)$taxnode_id)+1)
    
  }
  #
  ##### load version #####
 
   load_version()
  
  # ```
  # # load previous MSLs
  # ```{r load prev MSL, echo=F}
  
  #
  # this uses the DATABASE schema for the taxonomy_node table (ie, naming convention)
  #
  
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
  taxonomyDt=fread(file=dbTaxonomyNodeFilename,
                   header=TRUE,#header=FALSE,col.names=names(taxonomy_node_names),
                   colClasses=as.character(taxonomy_node_names),
                   stringsAsFactors=FALSE,na.strings=c("","NULL"))
  cat("Previous taxa:",dim(taxonomyDt), " from ",dbTaxonomyNodeFilename,"\n")
  if( !("host_source" %in% names(taxonomyDt)) ) {
    cat("WARNING: no host_source column in taxonomy_node dump!!! (Adding)\n")
    taxonomyDt[,"host_source"] = NA_character_
  }
  
  # ```
  # 
  # # load db CV tables 
  # ```{r load db cv tables, echo=F}
  
  #### SECTION load CVs and preMSLs #### 
  
  #
  # this uses the DATABASE schema for the taxonomy_node table (ie, naming convention)
  #
  
  # DB CV loads
  dbCvList = list()
  dbCvMapList = list()

  ##### CVs from db dumps #####
  dbRankFilename=file.path(params$ref_dir, params$db_rank_fname)
  rankCV = read.delim(file=dbRankFilename,header=TRUE, stringsAsFactors=TRUE,na.strings=c("NULL"))
  rownames(rankCV) = rankCV$id
  dbCvList[["rank"]] = rankCV
  dbCvMapList[["rank"]] = rankCV$id
  names(dbCvMapList[["rank"]]) = rankCV$name
  if(params$verbose) {cat("RankCV: ", dim(rankCV), " from ",dbRankFilename,"\n")}
  
  
  dbMoleculeFilename=file.path(params$ref_dir, params$db_molecule_fname)
  moleculeCV = read.delim(file=dbMoleculeFilename,header=TRUE, stringsAsFactors=TRUE,na.strings=c("NULL"))
  rownames(moleculeCV) = moleculeCV$id
  dbCvList[["molecule"]] = moleculeCV
  dbCvMapList[["molecule"]] = moleculeCV$id
  names(dbCvMapList[["molecule"]]) = moleculeCV$abbrev
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
  cvList=list()
  for(cv_col in 1:ncol(templateProposalCV)) {
    cv_name = templateProposalCV[1,cv_col]
    cv = templateProposalCV[,cv_col][-1]
    # clean UTF8-NB_space and other unprintable whitespaces
    cvClean = gsub("\u00A0"," ",cv) # UTF8-NBSP
    cvClean = gsub("[^[:alnum:][:punct:]]+"," ",cvClean) # anything else
    cvList[[cv_name]]=c(cvClean[!is.na(cvClean)],NA)
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
  names(cvList)=cvNameMap[names(cvList)]
  
  # remove ("Please select",NA) from "change" & "rank" CVs - that is a required field
  for( cv in c("change","rank") ) { #},"scAbbrev","scName","hostSource") ) {
    isRemoveTerm = tolower(gsub("[^[:alnum:]]","",cvList[[cv]])) %in% c(
      "pleaseselect", # 2023
        NA
    ) 
    cvList[[cv]] = cvList[[cv]][!isRemoveTerm]
  }

    # clean UTF8-NB_space (mac) 
  for( cv in c("change","rank","scAbbrev","scName") ) {
    for( i in seq(1:length(cvList[[cv]])) ) {
        term = cvList[[cv]][i]
        if( is.na(term) || regexpr(text=term, pattern=paste0(AscToChar(194),AscToChar(160))) > 0 ) {
          term = gsub(paste0("(",AscToChar(194),AscToChar(160),")+")," ",term)
          cvList[[cv]][i] = term
        }
    }
  }
  
  #
  ##### build Subcommittee Map   ##### 
  #
  scAbbrevNameMap = cvList[["scName"]]
  names(scAbbrevNameMap) = cvList[["scAbbrev"]]
  # remove them from "official" cvList, as that list 
  # also servers as a list of require columns for the proposal xlsx
  cvList = cvList[!names(cvList) %in% c("scName","scAbbrev")]
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
  ##### augment hostSource CV from db   ##### 
  #
  # this is just names, no IDs
  #
  dbHostSourceFilename=file.path(params$ref_dir, params$db_host_source_fname)
  hostSourceCV = read.delim(file=dbHostSourceFilename,header=TRUE, stringsAsFactors=TRUE,na.strings=c("NULL"))
  #rownames(hostSourceCV) = hostSourceCV$id
  dbCvList[["hostSource"]] = hostSourceCV$host_source
  dbCvMapList[["hostSource"]] = hostSourceCV$host_source
  names(dbCvMapList[["hostSource"]]) = hostSourceCV$host_source
  if(params$verbose) {cat("HostSourceCV: ", dim(hostSourceCV), " from ",dbHostSourceFilename,"\n")}
  cvList[["hostSource"]] = union(cvList[["hostSource"]],dbCvList[["hostSource"]])
  #
  ##### CVs from VMR xlsx  #####
  #
  vmrFilename = file.path(params$ref_dir, params$vmr_fname)
  vmrDf = data.frame(
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
  
  #
  
  # 
  #### create allErrorDf   #### 
  # 
  # error reporting data frame
  # 
  allErrorDf = data.table(
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
    "validator_version" = character()
  )
  .GlobalEnv$loadErrorDf = allErrorDf %>% filter(FALSE)

  
  #
  #### setup taxonomies: last,old, & cur  #### 
  #
  lastMSL = max(as.integer(taxonomyDt$msl_release_num))
  .GlobalEnv$oldMSLs = subset(taxonomyDt, msl_release_num < lastMSL)
  .GlobalEnv$curMSL = subset(taxonomyDt, msl_release_num==lastMSL)
  # add accounting columns
  .GlobalEnv$curMSL[,"out_updated"] = FALSE
  
  #
  # copy prev MSL to make new MSL,
  # to which we will try and apply these edits
  #
  .GlobalEnv$newMSL=createNewMSL(.GlobalEnv$curMSL,lastMSL, lastMSL+1, params$taxnode_delta)
  
  }
#
#
#### SAVE REF CACHE ####
#
#
if( params$update_cache ) {
  if(params$verbose){ cat("WRITE: cache file ", cacheFilename, "\n")}
  params$update_cache=FALSE
  save.image(file=cacheFilename)
}
if(params$verbose){ cat(paste0("PROPOSALS_DIR:  ",params$proposals_dir,"/\n")) }
if(params$verbose){ cat(paste0("OUT_DIR:        ",params$out_dir,"/\n")) }
dir.create(params$out_dir,recursive=T,showWarnings = F)
#
#
#### SECTION scan for files ####
#
#
inputFiles = data.frame(docpath=list.files(path=params$proposals_dir,
                                           pattern="^[^~.].*\\.(doc|xls)x*$", 
                                        recursive=T, full.names=TRUE)
)
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
  .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
  write_error_summary(.GlobalEnv$loadErrorDf)
  cat("# ERROR: NO_INPUT_FILES found in ",params$proposals_dir, "\n")
  quit(save="no", status=1)
}

# parse file paths/names
inputFiles$path         = dirname(inputFiles$docpath)
inputFiles$file         = basename(inputFiles$docpath)
inputFiles$basename     = gsub("(.*).(doc|xls)x*$","\\1",inputFiles$file)
# code may not be in early versions of filenames
parseableFilenames = grep(inputFiles$basename,pattern="^[0-9][0-9][0-9][0-9]\\.[0-9][0-9][0-9][A-Z]\\.")
inputFiles[parseableFilenames,"code"]          = sub("^([0-9]+\\.[0-9]+[A-Z]).*","\\1",inputFiles[parseableFilenames,]$basename)
inputFiles[parseableFilenames,"scAbbrev"]      = sub("^[0-9][0-9][0-9][0-9]\\.[0-9][0-9][0-9]([A-Z])\\.","\\1",inputFiles[parseableFilenames,"code"] )


# remove all *.Ud.* files, when in draft/final modes
if( params$processing_mode %in% c("draft","final") ) {
  filtered = grep(inputFiles$file, pattern=".*\\.Ud\\..*")
  if( length(filtered)) {
    errorDf = inputFiles[filtered,]
    errorDf$level = "INFO"
    errorDf$error = "IGNORE_FNAME_.Ud."
    errorDf$message = "All files with '.Ud.' in filename are ignored"
    .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
    write_error_summary(.GlobalEnv$loadErrorDf)
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
  .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
  write_error_summary(.GlobalEnv$loadErrorDf)
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
  .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
  write_error_summary(.GlobalEnv$loadErrorDf)
  if(params$tmi) { cat(paste0("# Appendix files filtered out: N=",length(filtered),"\n"))}
}
inputFiles = inputFiles[grep(inputFiles$file,pattern="appendix",ignore.case = T, invert=T),]
if(params$tmi) { cat("# xls|doc(x) files after appendix removal: N=",nrow(inputFiles),"\n")}

##### SECTION scan DOCX #####
#
# break filename into proposal code, filename and basename
#
docxs          = inputFiles[grep(inputFiles$file,pattern="\\.docx*$"),,drop=FALSE]
names(docxs)   = c("docxpath","path","docx","basename","code","scAbbrev")

if(params$tmi) { cat("# doc(x) files found: N=",nrow(docxs),"\n")}


# check for duplicate Proposal IDs
dups = duplicated(docxs$code)
allDups =docxs$code %in% docxs$code[dups]
if(sum(dups) > 0) {
  errorDf = docxs[allDups, c("scAbbrev", "code", "docx")]
  errorDf$level = "ERROR"
  errorDf$error = "DUPCODE.DOCX"
  errorDf$message = "duplicate proposal ID"
  .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
  # output: @DD make this web-app friedly
  #kable(errorDf,caption = paste0("QC01: ERRORS: dupliate docx proposal IDs"))
  write_xlsx(x=errorDf,path=file.path(params$out_dir,"QC01.docx_duplicate_ids.xlsx"))
  write_error_summary(.GlobalEnv$loadErrorDf)
}
# use codes for row names, but fall back to ordinals, if no code in filename
rownames(docxs)=ifelse(is.na(docxs$code),rownames(docxs),docxs$code)
docxs[,"code"] = rownames(docxs)

# spaces in filenames
spacedOut = grep(pattern=" ",docxs$docx)
if(sum(spacedOut) > 0) {
  errorDf = docxs[spacedOut, c("scAbbrev", "code", "docx")]
  errorDf$level = "WARNING"
  errorDf$error = "DOCX_FILENAME_SPACES"
  errorDf$message = "filename contains a space: please replace with _ or -"
  errorDf$notes = gsub("( )","[\\1]",docxs[spacedOut,]$docx)
  .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
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
    .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
  }
}
#
##### SECTION: scan XLSX ##### 
#

xlsxs          = inputFiles[grep(inputFiles$file,pattern="\\.xlsx*$"),,drop=FALSE]
names(xlsxs)   = c("xlsxpath","path","xlsx","basename","code","scAbbrev")

if(params$tmi) { cat("# xls(x) files found: N=",nrow(xlsxs),"\n")}

# QC for duplicate codes
dups = duplicated(xlsxs$code) & !is.na(xlsxs$code)
allDups =xlsxs$code %in% xlsxs$code[dups]
if( sum(dups) > 0 ) {
  # error details
  errorDf = xlsxs[allDups, c("scAbbrev", "code", "xlsx")]
  errorDf$level = "ERROR"
  errorDf$error = "XLSX_DUPCODE"
  for( code in xlsxs$code[dups] ) {
    # build list of all xlsxs using duplicate codes
    errorDf[errorDf$code==code,"message"] = 
      paste0("duplicate proposal ID: ",
             paste(xlsxs$xlsx[xlsxs$code==code],collapse=","))
    
  }
  # append to global list
  .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
  # output
  write_error_summary(.GlobalEnv$loadErrorDf)
  
  # terminate
  cat("ERROR: can not proceed with duplicated proposal IDs:", paste(levels(as.factor(xlsxs[dups,"code"])),collapse=","),"\n")
  stop(1)
}
rownames(xlsxs) = ifelse(is.na(xlsxs$code),rownames(xlsxs),xlsxs$code)
xlsxs[,"code"] = rownames(xlsxs)
#
# check that xlsx names match format
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
    .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
  }
  write_error_summary(.GlobalEnv$loadErrorDf)  
}
# spaces in XLSX filenames
spacedOut = grep(pattern=" ",xlsxs$xlsx)
if(sum(spacedOut) > 0) {
  errorDf = xlsxs[spacedOut, c("scAbbrev", "code", "xlsx")]
  errorDf$level = "WARNING"
  errorDf$error = "XLSX_FILENAME_SPACES"
  errorDf$message = "filename contains a space: please replace with _ or -"
  errorDf$notes = gsub("( )","[\\1]",xlsxs[spacedOut,]$xlsx)
  .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
  write_error_summary(.GlobalEnv$loadErrorDf)  
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
  .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
  write_error_summary(.GlobalEnv$loadErrorDf)
  cat("# ERROR: NO_INPUT_FILES found in ",params$proposals_dir, "\n")
  quit(save="no", status=1)
}

#
#### merge XLSX list into DOCX list, verify ####
#

proposals = data.frame(
  row.names=c(union(rownames(xlsxs),rownames(docxs)))
  )
proposals$code = ifelse(is.na(xlsxs[rownames(proposals),"code"]),rownames(proposals),xlsxs[rownames(proposals),"code"])
# XLSX only fields
proposals$xlsx     = xlsxs[rownames(proposals),"xlsx"]
proposals$xlsxpath = xlsxs[rownames(proposals),"xlsxpath"]
# DOCS only fields
proposals$docx = docxs[rownames(proposals),"docx"]
proposals$docxpath = docxs[rownames(proposals),"docxpath"]
# MERGE 
proposals$basename = ifelse(!is.na(xlsxs[proposals$code,"basename"]),
                            xlsxs[proposals$code,"basename"],
                            docxs[proposals$code,"basename"])
proposals$scAbbrev = ifelse(!is.na(xlsxs[proposals$code,"scAbbrev"]),
                            xlsxs[proposals$code,"scAbbrev"],
                            docxs[proposals$code,"scAbbrev"])
# strip off version, workflow status and .fix, to get final, production filename
proposals$cleanbase= gsub("^([0-9]+\\.[0-9]+[A-Z])(\\.[A-Z]+)(\\.v[0-9]+)*(\\.fix)*(\\..*)$","\\1\\5",proposals$basename)

# QC  - missing xlsx file
missing= is.na(proposals$xlsx)
if( sum(missing) > 0 ) {
  errorDf= proposals[missing,c("code","docx")]
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
      proposals[row,"xlsx"]=xlsxs[xlsxs$code==errorDf[row,"code"],"xlsx"]
      proposals[row,"xlsxpath"]=xlsxs[xlsxs$code==errorDf[row,"code"],"xlsxpath"]
      errorDf[row,]$xlsx = proposals[row,"xlsx"]
      errorDf[row,]$row = NA
      errorDf[row,]$level = "WARNING"
      errorDf[row,]$error = "XSLX.TYPO"
      errorDf[row,]$notes = paste("Using best guess:",proposals[row,"xlsx"])
    } else if( guesses > 1 ) {
      # just list all the options we found
      errorDf[row,]$error = "XSLX.MULTIPLE"
      errorDf[row,]$error = paste0("Multiple .xlsx files start with code '",errorDf[row,"code"],"', and aren't marked '",params$infile_suppl_pat,"'")
      errorDf[row,"notes"] = paste("SUGGESTIONS:",paste(paste0(xlsxs[xlsxs$xlsxID==errorDf[row,"docxID"],"basename"],".xslx"),collapse=", "))
    } else {
      # 
      # can't find the xlsx
      #
      # let the error pass through
      proposals[]
    }
  }
  # suppress name-mismatch warnings
  loadErrorDfFilt = errorDf
  if( params$processing_mode == "validate" ) {
      # filter out some errors
      loadErrorDfFilt = errorDf %>% filter(error != "XSLX.TYPO")
  }
  if(nrow(loadErrorDfFilt) > 0) {
    # append to global list
    .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
   }
  write_error_summary(.GlobalEnv$loadErrorDf)  
}
# QQQ don't check for missing docx files? 

#
# get SC names from last letter of code
#
proposals$scAbbrev = NA
parsableCode = grep(proposals$code, pattern="[0-9][0-9][0-9][0-9]\\.[0-9][0-9][0-9]([A-Z])" )
if( length(parsableCode) >0 ) {
  proposals[parsableCode,"scAbbrev"] =  gsub("[0-9][0-9][0-9][0-9]\\.[0-9][0-9][0-9]([A-Z])","\\1",proposals$code[parsableCode])
}
# QC
badProposalAbbrevs = !(proposals$scAbbrev %in% names(scAbbrevNameMap))
if(sum(badProposalAbbrevs)>0) {
  if( params$processing_mode %in% c("draft","final") ) {
    # supress picky error message
    errorDf = proposals[badProposalAbbrevs, c("scAbbrev", "code", "xlsx", "docx")]
    errorDf$level = "WARNING"
    errorDf$error = "CODE_BAD_SC_ABBREV"
    errorDf$message = "Last letter of CODE not a valid ICTV Subcommittee letter"
    errorDf$notes = paste0("'",proposals[badProposalAbbrevs,"scAbbrev"],"' not in [",
                           paste0(names(scAbbrevNameMap),collapse=","),"]")
    .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
  }
  write_error_summary(.GlobalEnv$loadErrorDf)  
}
proposals$subcommittee = ifelse(badProposalAbbrevs,
                          ifelse(is.na(proposals$scAbbrev),"unspecified",paste0("unknown-",proposals$scAbbrev)),
                          scAbbrevNameMap[proposals$scAbbrev])


# ```
# 
# # Load summary
# 
# ```{r load_summary,echo=FALSE}
# kable(caption=paste0("SUMMARY: Proposal .xlsx file found in ",params$proposals_dir,"/"),
#       x=data.frame(nProposals=summary(as.factor(proposals$subcommittee))),)
# ```
# 
# ## QC setup
# ```{r qc_setup_template_versions, echo=TRUE}

####  SECTION: QC setup #### 
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
# normalized across XSLX template versions
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

# ```
# # process changes function
# 
# ```{r process_changes_helper_functions, echo=FALSE}
#
# this uses the DATABASE schema for the taxonomy_node table (ie, naming convention)
#

realmSpeciesCols = 1:15

getTaxon = function(realmSpecies) {
  # 
  # get the NAME of the lowest ranked taxon
  #
  nonNAs = !is.na(realmSpecies[realmSpeciesCols])
  if( sum(nonNAs) > 0 ) {
    return(realmSpecies[max(which(nonNAs))])
  } else {
    # no taxon - could be the src of a create, or the dest of abolish
    #if( params$debug_on_error) { print("INTERNAL ERROR: getTaxon()=NA"); browser()}
    return(NA)
  }
}
getParentTaxon = function(realmSpecies) {
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
    #if( params$debug_on_error) { print("INTERNAL ERROR: getParentTaxon()=NA"); browser()}
    return(NA)
  }
}
getLineage = function(realmSpecies,masks=c()) {
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
    #if( params$debug_on_error) { print("INTERNAL ERROR: getLineage()=NA"); browser()}
    return(NA)
  }
}
getParentLineage = function(realmSpecies) {
  nonNAs = !is.na(realmSpecies[realmSpeciesCols])
  if( sum(nonNAs) > 0 ) {
    rank=getTaxonRank(realmSpecies)
    return(getLineage(realmSpecies,masks=c(rank)))
  } else {
    # no taxon - could be the src of a create, or the dest of abolish
    #if( params$debug_on_error) { print("INTERNAL ERROR: getParentTaxon()=NA"); browser()}
    return(NA)
  }
}
getTaxonRank = function(realmSpecies) {
  nonNAs = !is.na(realmSpecies[realmSpeciesCols])
  if( sum(nonNAs) > 0 ) {
    # 
    # get the rank of the lowest-rank taxon 
    #
    return(as.character(rankCV$name[max(which(nonNAs))+1]))
  } else {
    # no taxon - could be the src of a create, or the dest of abolish
    #if( params$debug_on_error) { print("INTERNAL ERROR: getParentTaxon()=NA"); browser()}
    return(NA)
  }
}

#
# error formating functions
#

diffLineageStrings=function(lin1, lin2) {
  lin1ex=unlist(strsplit(as.character(lin1),";"))
  lin2ex=unlist(strsplit(as.character(lin2),";"))
  
  # match lengths - would be better if we knew the actual ranks, and could align them
  if(length(lin1ex) < length(lin2ex)) {
    lin1ex=c(lin1ex,rep("",length(lin2ex)-length(lin1ex)))
  } else if( length(lin1ex) > length(lin2ex)) {
    lin2ex=c(lin2ex,rep("",length(lin1ex)-length(lin2ex)))
  }
 
  diffs=lin1ex==lin2ex
  tups=ifelse(diffs,lin1ex,paste0('[',lin1ex,'//',lin2ex,']'))

  return(paste0(
    ifelse(sum(!diffs)==0,"[identical]",""),
    paste0(tups,collapse=";")))
}
# 
# # QC functions
# 
# ```{r qc_functions_defined, echo=FALSE}

#
# this uses the PROPOSAL.XLSX schema (naming convention)
#

#
# load .xlsx file into DF
#
load_proposal = function(code) {
  cat("# LOAD_PROPOSAL(",code,")\n")

  # return data
  errorDf = .GlobalEnv$allErrorDf[FALSE,]
  proposalDf = NULL
  
  # get worksheet names
  proposalsSheetNameRegex = "proposal*"
  sheetNames = suppressMessages(
    excel_sheets(proposals[code,"xlsxpath"])
  )
  #
  # comapre to expected sets of sheet names, warn if they added extra sheets
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
    errorDf=addError(errorDf,
                     code,"","","","",
                     "WARNING","XLSX_EXTRA_SHEETS",
                     "XLS file has additional sheets not present in template",
                     paste0("Worksheet(s) named '",paste0(extraSheets,collapse="','"),"' were added ")
    )
  }
  
  # Try and find the one sheet that contains the proposal data
  proposalsSheetNames = grep(proposalsSheetNameRegex,sheetNames,ignore.case = T, value=T)
  if( is.na(proposalsSheetNames[1]) ) {
    # ERROR if can't find it. 
    cat(paste0("code ", code, ": A worksheet name matching '",proposalsSheetNameRegex,"' was not found","logging error","\n"))
    errorDf=addError(errorDf,
                     code,"","","","",
                     "ERROR","XLSX_NOT_PROPOSAL",
                     "XLS does not an appropriately named 'proposal' sheet",
                     paste0("A worksheet name matching '",proposalsSheetNameRegex,"' was not found, only: ",
                            "'",paste0(sheetNames,collapse="','"),"'"
                            )
                     )
  } else if( length(proposalsSheetNames) > 1 ) {
    # ERROR if find more than one
    cat(paste0("code ", code, ": ERROR: More than 1 worksheet named like '",proposalsSheetNameRegex,"' found: '",paste0(proposalsSheetNames,collapse="','"),"'\n"))
    errorDf=addError(errorDf,
                     code,"","","","",
                     "ERROR","XLSX_MULTI_PROPOSAL", 
                     "XLS file has more than one sheet named 'proposal'",
                     paste0("More than 1 worksheet named  '",proposalsSheetNameRegex,"' were found: '",paste0(proposalsSheetNames,collapse="','"),"'")

    )
  } else {
    # fOK: ound one and only one proposal sheet - read that one
    
    # read XLSX file, no column names
    proposalDf = data.frame(
      # suppress errors about column names
      # https://github.com/tidyverse/readxl/issues/580#issuecomment-519804058
      suppressMessages(
        read_excel(
          proposals[code,"xlsxpath"],
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
  return(list(proposalDf=proposalDf,errorDf=errorDf))
  
}

#
# Load  .docx into DF
#
load_proposal_docx = function(code) {
  cat("# LOAD_PROPOSAL_DOCX(",code,")\n")

  # return data
  errorDf = allErrorDf[FALSE,]
  metaDf = data.frame(code=c(code))

  # check if we even HAVE a docx file!
  if( is.na(proposals[code,"docxpath"]) ) {
    if( params$processing_mode %in% c("draft","final") ) {
      # supress picky error message
      errorDf = proposals[code, c("code", "xlsx")]
      errorDf$level = "WARNING"
      errorDf$error = "DOCX_MISSING"
      errorDf$message = ".docx file not available"
      .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
    }
    return(list(metaDf=metaDf,errorDf=errorDf))
  }
  # read DOCX file, no column names
  txt = read_docx(proposals[code,"docxpath"])
  
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
    errorDf = proposals[code, c("subcommittee", "code", "xlsx", "docx")]
    errorDf$level = "WARNING"
    errorDf$error = "DOCX_TITLE_MISSING"
    errorDf$message = "Title not found"
    errorDf$notes = paste0("Line containing expression '",authorsPattern,"' not found")
    .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
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
    errorDf = proposals[code, c("subcommittee", "code", "xlsx", "docx")]
    errorDf$level = "WARNING"
    errorDf$error = "DOCX_AUTHORS_MISSING"
    errorDf$message = "Authors list not found"
    errorDf$notes = paste0("Line containing expression '",authorsPattern,"' not found")
    .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
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
    errorDf = proposals[code, c("subcommittee", "code", "xlsx", "docx")]
    errorDf$level = "WARNING"
    errorDf$error = "DOCX_CORR_AUTHOR_MISSING"
    errorDf$message = "Corresponding Author not found"
    errorDf$notes = paste0("Line containing expression '",corrAuthorPattern,"' not found")
    .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
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
    errorDf = proposals[code, c("subcommittee", "code", "xlsx", "docx")]
    errorDf$level = "WARNING"
    errorDf$error = "DOCX_ABSTRACT_MISSING"
    errorDf$message = "Abstract not found"
    errorDf$notes = paste0("Line containing expression '",abstractPattern,"' not found")
    .GlobalEnv$loadErrorDf = rbindlist(list(.GlobalEnv$loadErrorDf, errorDf),fill=TRUE)
  }

  # return value
  return(list(metaDf=metaDf,errorDf=errorDf))
}

# QC the .xlsx data at the sheet level
#
# returns: errorDf
addError=function(errorDf,code,row,change,rank,taxon,levelStr,errorCode,errorStr,notes) {
  #browser()
  nextErrorDf = data.frame(
    subcommittee = proposals[code,]$subcommittee,
    code = code,
    row = as.numeric(row),
    change = change, 
    rank=rank,
    taxon = taxon,
    docx = ifelse(is.na(proposals[code,]$docx),"MISSING",proposals[code,]$docx),
    xlsx = ifelse(is.na(proposals[code,]$xlsx),"MISSING",proposals[code,]$xlsx),
    level = levelStr,
    error = errorCode,
    message = errorStr, 
    notes = notes)
  # if ERROR, debug, if desired
  if( params$debug_on_error && interactive() && levelStr=="ERROR") { browser()}
  # return appended error list
  return(rbind(errorDf,nextErrorDf,fill=TRUE))
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
  if(!is.null(proposalDf[2,1]) && (substring(proposalDf[2,1],1,13) == "version 2023.")) {
    
    templateVersion = substring(proposalDf[2,1],9,13)
    # 
    # process 2023 template layout
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
      
      proposals[code,"templateVersion"]="error"   
      errorDf=addError(errorDf,code,3,NA,NA,NA,"ERROR","XLSX.TEMPLATE_UNK","XLSX template version",
                              paste0("ROW2 is ", templateLine2Version, " (",templateLine2Error,")",
                                ", ROW3 is ", templateLine3Version," (",templateLine3Error,")")
      )
      return(list(errorDf=errorDf))
    } else {
      templateVersion = templateLine2Version
    }
    # finish up templateVersion
    proposals[code,"templateVersion"]=templateVersion
    if( templateVersion == "v1") {
      # WARNING for outdated templates
      errorDf=addError(errorDf,code,3,NA,NA,NA,"INFO","XLSX.OLD_TEMPLATE_V1", "XLSX template version",
                       paste0("You are using version ",templateVersion,". Please get the latest version from ",params$templateURL)
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
    if( codeValue != code ) {
      if( str_starts(codeValue,"Code") ) {
        if(params$show.xlsx.code_miss) {
          errorDf=addError(errorDf,code,codeRow,NA,NA,NA,"INFO","XLSX.CODE_MISS", "XLSX code missing", 
                           paste0("XLSX cell ",codeCell,
                                  " is ", "'",codeValue, 
                                  "'; replace with the actual code: '",code,"'")
          )
        }
      } else {
        errorDf=addError(errorDf,code,codeRow,NA,NA,NA,"WARNING", "XLSX.CODE_BAD","XLSX code wrong", 
                         paste0("XLSX cell ",codeCell,
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
  if(templateVersion=="v2") { firstDataRow=4; changeDf = proposalDf[firstDataRow:nrow(proposalDf),xlsx_v2_change_cols] }
  if(templateVersion=="2023." ) {firstDataRow=5; changeDf = proposalDf[firstDataRow:nrow(proposalDf),xlsx_2023_change_cols]}
  colnames(changeDf) = xlsx_change_colnames
  # a flag to exclude rows with irrecoverable errors
  changeDf[,".noErrors"] = TRUE
  changeDf[,".errors"] = NA_character_
  
  # find non-empty rows
  hasData = apply(changeDf[,xlsx_change_srcDest_colnames],1,function(row){return(sum(is.na(row))!=length(row))})
  changeDf=changeDf[hasData,]
  if( nrow(changeDf) == 0 ) {
    errorDf=addError(errorDf,code,firstDataRow,NA,NA,NA,"ERROR", "XLSX.EMPTY","XLSX no change rows found", "")
    return(list(errorDf=errorDf))
  } else {
    proposals[code,"nChanges"] = nrow(changeDf)
  }

  #
  #### extract src/dest taxon names for error reporting ####
  #
  # must be before QC spaces/quotes, etc, as error reporting there
  # uses .changeTaxon!
  changeDf$.srcTaxon =     apply(changeDf[,xlsx_change_srcCols], 1,getTaxon)
  changeDf$.srcRank =      apply(changeDf[,xlsx_change_srcCols], 1,getTaxonRank)
  changeDf$.srcLineage =   apply(changeDf[,xlsx_change_srcCols], 1,getLineage)

  changeDf$.destTaxon =         apply(changeDf[,xlsx_change_destCols],1,getTaxon)
  changeDf$.destRank =          apply(changeDf[,xlsx_change_destCols],1,getTaxonRank)
  changeDf$.destLineage =       apply(changeDf[,xlsx_change_destCols], 1,getLineage)
  changeDf$.destParentName =    apply(changeDf[,xlsx_change_destCols], 1,getParentTaxon)
  changeDf$.destParentLineage = apply(changeDf[,xlsx_change_destCols], 1,getParentLineage)
  
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
    pattern=AscToChar(202); pat_warn="non-breaking space character"; pat_replace=" "
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
        errorDf=addError(errorDf,code,rownames(changeDf)[qc.matches],
                         changeDf$change[qc.matches],changeDf$rank[qc.matches],changeDf$.changeTaxon[qc.matches],
                         "INFO","XLSX.QUOTES_REMOVED", paste("XLSX has",pat_warn),
                         paste0(paste(col,gsub(pattern,"[\\1]",changeDf[qc.matches,col]),sep=":")," (replacing with '",pat_replace,"')")
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
    if( params$tmi ) {
        cat(with(value_validation[i,],paste0("\tTMI: column_validation(",col,",",code,"):s/",regex,"/",replace,"/\n")))
    }
    if( value_validation[i,]$type == "replace") {
      # 
      # check for presence of regex's to be replaced
      #
      qc.matches =grep(changeDf[,col],pattern=value_validation[i,]$regex)
      if(length(qc.matches)>0) { 
        if(params$verbose) { cat(paste0(value_validation[i,]$class,":",error," has ",length(qc.matches)," cells with ",value_validation[i,]$warn," in column '",col,"'\n")) }
        #browser()
        charEncodings = ""
        if(params$tmi) { 
          # show code points in error message, mostly for I18N problems
          charEncodings = paste0("; charEncodings",
                                 paste(unlist(strsplit(changeDf[qc.matches,col],"")),
                                       CharToAsc(unlist(strsplit(changeDf[qc.matches,col],""))), 
                                       sep="=", collapse=","))
        }
        errorDf=addError(errorDf,code,rownames(changeDf)[qc.matches],
                         changeDf$change[qc.matches],changeDf$rank[qc.matches], changeDf$.changeTaxon[qc.matches],
                         value_validation[i,]$class,error, paste("XLSX has",value_validation[i,]$warn),
                         paste(col,":",gsub(value_validation[i,]$regex,"[\\1]",changeDf[qc.matches,col]),
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
        if(params$verbose) { cat(value_validation[i,]$class,":",error,"has",sum(qc.matches),"cells with",value_validation[i,]$warn,"in column",col,"\n") }
        errorDf=addError(errorDf,code,rownames(changeDf)[qc.matches],
                         changeDf$change[qc.matches],changeDf$rank[qc.matches], changeDf$.changeTaxon[qc.matches],
                         value_validation[i,]$class,error,value_validation[i,]$warn,
                         paste(col,gsub(value_validation[i,]$regex,"[\\1]",changeDf[qc.matches,col]),sep=":"))
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
  missingCvNames = names(cvList)[!names(cvList) %in% names(changeDf)]
  if(length(missingCvNames)>0) {
    if(params$verbose) { cat("ERROR:",code,"missing CV columns [",paste(missingCvNames, collapse=","),"] on row 3\n") }
    errorDf=addError(errorDf,code,3,
                     NA,NA,NA,
                     "ERROR","XLSX.MISSING_COLUMN",paste("XLSX missing column"),
                     paste("XLSX missing header [",missingCvNames,"]")) 
    return(list(errorDf=errorDf))
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
  for(cv in names(cvList)) {
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
      if(cv %in% c("genomeCoverage","molecule","hostSource") ) {
          # remove value, issue WARNING
          changeDf[rownames(badTerms),cv] = NA
          errorDf=addError(errorDf,code,rownames(badTerms),
                           badTerms$change,badTerms$rank, badTerms$.changeTaxon,
                           "WARNING","XLSX.INVALID_TERM",paste("XLSX incorrect term in column",cv),
                           paste("XLSX incorrect value [",badTerms[,cv],"]. Valid terms: [",paste0(cvList[[cv]],collapse=","),"]")) 
      } else {
         changeDf[rownames(badTerms),".noErrors"] = FALSE
         # issue error and skip this line
         errorDf=addError(errorDf,code,rownames(badTerms),
                          badTerms$change,badTerms$rank, badTerms$.changeTaxon,
                          "ERROR","XLSX.INVALID_TERM",paste("XLSX incorrect term in column",cv),
                          paste("XLSX incorrect value [",badTerms[,cv],"]. Valid terms: [",paste0(cvList[[cv]],collapse=","),"]")) 
      }
      #browser()
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
          errorDf=addError(errorDf,code,rownames(flawedTerms),
                           flawedTerms$change, flawedTerms$rank, flawedTerms$.changeTaxon,
                           "TMI","XLSX.TYPO_TERM",paste("fixed term with typo (space,caps) in column ",cv),
                           paste0("Term '",flawedTerms[,cv],"' replaced with '",correctionsDf[rownames(flawedTerms),cv],"'. Valid terms: [",paste0(cvList[[cv]],collapse=","),"]")) 
        }
      }
      # correct
      changeDf[rownames(flawedTerms),cv] = correctionsDf[rownames(flawedTerms),cv]
    } # if flawedTerms
  } # for cvList
  
  #
  # sort proposal by action/rank
  #
  # first all "abolish" & "merge" from bottom to top rank
  # than all other actions from top rank to bottom
  #
  
  # action sort
  destroyRows = grep(changeDf$change, pattern="(abolish|merge)", ignore.case = T)
  createRows  = grep(changeDf$change, pattern="(abolish|merge)", ignore.case = T, invert = T)
  changeDf[destroyRows,".sortAction"] = 1
  changeDf[createRows, ".sortAction"] = 2
  
  # rank sort - opposite directions depending on .sortAction
  rankMap = dbCvList[["rank"]]$id
  names(rankMap) = dbCvList[["rank"]]$name
  changeDf[,".sortRank"] = dbCvMapList[["rank"]][changeDf$rank] * ifelse(changeDf$.sortAction==1,-1,1)
  
  # final sort
  actionRankOrder = order(changeDf$.sortAction, changeDf$.sortRank)

  # return warnings, if any
  return(list(errorDf=errorDf,changeDf=changeDf[actionRankOrder,]))
} # qc_proposal()
# ```
# 
# 
# # Load and QC
# 
# ```{r load_and_qc}

#### SECTION: load and QC xlsx's #### 
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
codes = rownames(proposals) # all proposals
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
    xlsx_fname=proposals[code,"xlsx"]
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
        errorDf = results[["errorDf"]]
        cat("# LOADED: ",code," DOCX with ",nrow(errorDf)," errors/warnings\n")
        if(nrow(errorDf)>0){.GlobalEnv$loadErrorDf=rbindlist(list(.GlobalEnv$loadErrorDf,errorDf),fill=TRUE)}
      } else {
        cat("#  DOCX FROM CACHE: ",code,"\n")
      }
      proposals[code,names(docxList[[code]])] = docxList[[code]]
      #
      # try to load proposal from XLSX
      #
      if(is.null(xlsxList[[code]])){
        # load raw xlsx file
        results =  load_proposal(code)
        changeDf = results[["proposalDf"]]
        errorDf  = results[["errorDf"]]
        cat("# LOADED: ",code," XLS with ",nrow(errorDf)," errors/warnings\n")
        if(nrow(errorDf)>0){.GlobalEnv$loadErrorDf=rbindlist(list(.GlobalEnv$loadErrorDf,errorDf),fill=TRUE)}
        xlsxList[[code]] = results[["proposalDf"]]
        
        # if load had errors/warnings, update error file
        if( nrow(errorDf)>0  ){
          write_error_summary(.GlobalEnv$loadErrorDf)          
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
        errorDf = results[["errorDf"]]
        cat("# QCed:     ",code," with ",nrow(errorDf)," errors/warnings\n")
        if(nrow(errorDf)>0){.GlobalEnv$loadErrorDf=rbindlist(list(.GlobalEnv$loadErrorDf,errorDf),fill=TRUE)}
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
write_error_summary(.GlobalEnv$loadErrorDf)
#kable(allErrorDf,caption = paste0("QC: Summary of ERRORs and WARNINGs"))
 
# DEBUG: break to RStudio debugger after file load/QC
#if(interactive()){browser()}

# ```
# 
# # recursive lineage setter
# 
# ``` {r recursive_lineage, echo=FALSE}
#
# INPUTs:
#   taxnode_id = node whose children need updating
#   destLineage = lineage of taxnode_id
#
update_lineage = function(parent_id,parent_lineage){
  # 
  ct = 0
  # find kids, get their taxnode_is
  kid_rows= which(.GlobalEnv$newMSL$parent_id == parent_id)
  if(length(kid_rows)>0) {
    # update kids' lineage
    .GlobalEnv$newMSL[kid_rows,"lineage"] = paste(parent_lineage,.GlobalEnv$newMSL[kid_rows,]$name,sep=";")

    # recurse
    for(kid_row in kid_rows ) {
      if(params$verbose>1) {
        cat("update_lineage: ",
            .GlobalEnv$newMSL$taxnode_id[kid_row]," ",
            .GlobalEnv$newMSL$lineage[kid_row],
            " (",as.character(.GlobalEnv$newMSL$rank[kid_row]),")\n")
      }
      ct = ct + update_lineage(.GlobalEnv$newMSL$taxnode_id[kid_row],.GlobalEnv$newMSL$lineage[kid_row])
    }
  }
  return(ct)
}

# test
if( FALSE ) {
  oldLineage = curMSL$lineage
  # Poxviridae
  parent_id=202204737; parent_lineage=.GlobalEnv$newMSL$lineage[which(.GlobalEnv$newMSL$taxnode_id==parent_id)]
  ct = update_lineage(parent_id, parent_lineage)
  print(ct)
  sum(oldLineage!=.GlobalEnv$newMSL$lineage,na.rm=T)
}
# ```

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

# ```
# 
# ```{r process_changes_function, echo=FALSE}

apply_changes = function(code,proposalBasename,changeDf) {
  # proposalBasename=proposals[code,"basename"]; proposalBasename # debug
  proposalZip = paste0(proposalBasename,".zip")
  #
  # iterate over changes
  # 
  errorDf = allErrorDf[FALSE,]
  
  # genera to scan for binomial issues when we're done
  # QQQQ: should replace this with an admin column on newMSL!
  renamedGenera = c()
  
  #### admin columns ####
  .GlobalEnv$curMSL[,c(".split",".split_kept")] = FALSE
  .GlobalEnv$curMSL[,c(".split_code",".split_linenum",".split_acc_used")] = NA_character_

  .GlobalEnv$newMSL[,c(".split",".split_kept")] = FALSE
  .GlobalEnv$newMSL[,c(".split_code",".split_linenum",".split_acc_used")] = NA_character_

    # action mappings - externalize or move up
  actionCV = c(
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
  linesWoErrors = rownames(changeDf[changeDf$.noErrors,])
  for(linenum in linesWoErrors) {
    # linenum = rownames(changeDf)[1] ; linenum; changeDf[linenum,"change"];# debug
    # linenum = rownames(changeDf)[63] ; linenum; changeDf[linenum,"change"];# debug
    # linenum = "12"; linenum; changeDf[linenum,];
    # linenum = rownames(changeDf)[which(linenum==rownames(changeDf))+1]; linenum # debug: next row

    # get this change's line
    change = changeDf[linenum,]
    action = change$change
    
    #
    # check CVs
    #

    # clean action
    actionClean = actionCV[tolower(action)]
    if( is.na(actionClean) ) { 
      errorDf = with(changeDf[linenum,],addError(errorDf,code,linenum,change,rank,.changeTaxon,
                                                 "ERROR","ACTION.UNK","Unknown action/change",action))
      next;
    }
    
   
    # rank
    isRankOk = tolower(change$rank) %in% rankCV$name
    if( !isRankOk ) {
        errorDf=addError(errorDf,code,linenum,change$change, change$rank, change$.changeTaxon,
                         "WARNING", "RANK.UNK", "Rank name unknown", 
                         paste0(", rank=", change$rank, 
                                "; knownRanks=", paste(rankCV$name[-1],collapse = ","))
        )
    } else {
      rankClean = tolower(change$rank)
    }
    
    # linenum = rownames(changeDf)[1] # debug
    srcRealmSpecies = change[xlsx_change_srcCols]
    destRealmSpecies = change[xlsx_change_destCols]
    
    # should these be a vector operation, adding columns to changeDF? 
    # if we do that, we need to separate srcTaxon's label (rank) and value (name)
    # QQQQ: use the columns already in changeDf [.srcTaxon,.srcRank, etc]
    srcTaxonName  =getTaxon(srcRealmSpecies)
    srcTaxonRank  =getTaxonRank(srcRealmSpecies)
    srcLineage    =getLineage(srcRealmSpecies)
    destTaxonName =getTaxon(destRealmSpecies)
    destTaxonRank =getTaxonRank(destRealmSpecies)
    destParentName=getParentTaxon(destRealmSpecies)
    destLineage   =getLineage(destRealmSpecies)
    destParentLineage =getLineage(destRealmSpecies,masks=c(destTaxonRank))

    # check where this name may or may not exist already
    srcNewMatches= (.GlobalEnv$newMSL$name==as.character(srcTaxonName))
    destNewMatches=(.GlobalEnv$newMSL$name==as.character(destTaxonName))
    destCurMatches=(.GlobalEnv$curMSL$name==as.character(destTaxonName))
    destOldMatches=(.GlobalEnv$oldMSLs$name==as.character(destTaxonName))

    #
    # check that at least one of dest/src taxa were given
    #
    if(is.na(srcTaxonName) && is.na(destTaxonName)) {
         errorDf=addError(errorDf,code,linenum,change$change,change$rank,change$.changeTaxon, 
                          "ERROR", "XLSX.NO_SRC_DEST", 
                         paste0("Change=",toupper(action),", but neither 'current taxonomy', nor 'proposed taxonomy', were specified"), 
                         ""
        )
        next;
     
    }
    #
    # check if lineage not correct in "current" 
    #
    if( sum(srcNewMatches,na.rm=T)==1 ) {
      if(change$.srcLineage != .GlobalEnv$newMSL[srcNewMatches,"lineage"]) {
        
        errorDf=addError(errorDf,code,linenum,change$change,change$rank,change$.changeTaxon, 
                         "WARNING", "SRC.LINEAGE_WRONG", 
                         paste0("'current taxonomy' in proposal doesn't match current MSL. Typo?"), 
                         paste0("proposal 'current taxonomy'=", change$.srcLineage, 
                                "; MSL lineage=", .GlobalEnv$newMSL[srcNewMatches,"lineage"])
        )
      }
    }
    #
    # make sure new taxon does not exist already
    #
    if(!is.na(destTaxonName) && !(actionClean %in% c("move"))
       && !(actionClean == 'split' && !is.na(srcTaxonName) && destTaxonName == srcTaxonName))  {
      # check if exists in the current MSL
      if(sum(destCurMatches)>0) {
        matchDf = .GlobalEnv$curMSL[destCurMatches,]
        matchDf = matchDf[order(matchDf$msl_release_num,decreasing = T)[1],]
        matchLineage = paste0( "MSL", matchDf$msl_release_num,":",matchDf$lineage)
        errorDf=addError(errorDf,code,linenum,change$change,change$rank,change$.changeTaxon,
                         "ERROR", "DEST.IN_CUR", 
                         paste0("Change=",toupper(action),", but taxon name already exists"), 
                         paste0("proposed ",destTaxonRank,"=", destTaxonName, 
                                "; existing ",matchDf$rank,"=", matchLineage)
        )
        next;
      }
      # check if exists in a historical MLS (not the most recent)
      if(sum(destOldMatches)>0) {
        matchDf = .GlobalEnv$oldMSLs[destOldMatches,]
        matchDf = matchDf[order(matchDf$msl_release_num,decreasing = T)[1],]
        matchLineage = paste0( "MSL", matchDf$msl_release_num,":",matchDf$lineage)
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "ERROR", "DEST.IN_OLD", 
                         paste0("Change=",toupper(action),",  but taxon name existed historically"), 
                         paste0("proposed ",destTaxonRank,"=", destTaxonName, 
                                "; existing ",matchDf$rank,"=", matchLineage)
        )
        next;
      }
      # check if already created in newMSL
      if(sum(destNewMatches)>0) {
        matchDf = .GlobalEnv$newMSL[destNewMatches,]
        matchDf = matchDf[order(matchDf$msl_release_num,decreasing = T)[1],]
        matchLineage = paste0( "MSL", matchDf$msl_release_num,":",matchDf$lineage)
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "ERROR", "DEST.IN_NEW", 
                         paste0("Change=",toupper(action),", but taxon name already created in new MSL"), 
                         paste0("proposed ",destTaxonRank,"=", destTaxonName, 
                                "; existing ",matchDf$rank,"=", matchLineage,
                                "; otherProposal=", matchDf[1,"in_filename"]
                              )
        )
        next;
      }
      # WARN: check if taxon rank matches change rank
      if ( tolower(destTaxonRank) != tolower(change$rank) ) {
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "WARNING", "DEST.RANK_MISMATCH", 
                         paste0("Change=",toupper(action)," but proposed taxon rank does not match [rank] column"), 
                         paste0("rankColumn=", change$rank,
                                ", proposedTaxonRank=",destTaxonRank,
                                ", proposedTaxonomy=", destLineage)
        )
      }
    } # if(!is.na(destTaxonName))
    
    # debug
    #if(code=="2022.003M" && interactive()) {browser()}
    
    #  -------------------------------------------------------------------------
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
    #  -------------------------------------------------------------------------
    if(    actionClean %in% c("new")
       || (actionClean %in% c("split") && srcTaxonName != destTaxonName) 
       ) {

      # index of our "current"/src taxon,  in the newMSL
      srcNewTaxonIdx = 0 # not yet found
      srcCurTaxonIdx = 0 
      
      #
      # check SRC taxon
      #
      if(actionClean == "new" ) {
        # should have NO source taxon
        if( !is.na(srcTaxonName)) {
          # if srcTaxon was specified in xlsx (shouldn't be for NEW)
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                           "WARNING", "CREATE.W_SRC", "Change=CREATE, but 'current taxonomy' is not empty; perhaps you meant SPLIT", 
                           paste0("currentVsProposed=", diffLineageStrings(change$.srcLineage,change$.destLineage),
                                  ", currentTaxonomy=", change$.srcLineage,
                                  ", proposedTaxonomy=",change$.destLineage)
          )
        }
      } else if(actionClean == "split" ) {
        if( is.na(srcTaxonName) ) {
          # split MUST have a source taxon
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                           "ERROR", "SPLIT.NO_SRC", "Change=SPLIT, but 'current taxonomy' columns are empty", 
                           paste0("proposedTaxonomy=",change$.destLineage)
          )
          next;
        } else {
          # has a source taxon, let's find it
          srcCurTaxonIdx=(.GlobalEnv$curMSL$name==as.character(srcTaxonName))
          if( sum(srcCurTaxonIdx) == 0 ) {
            # split MUST have a source taxon
            errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                             "ERROR", "SPLIT.SRC_NO_EXIST", "Change=SPLIT, but 'current taxonomy' doesn't exist", 
                             paste0("currentTaxonomy=",change$.srcLineage)
            )
            next;
          }
          # remember this taxon has been split - may need to 
          # obsolete it, later, if no split directive keeps the old name
          
          # mark original, just in case that's usefule
          .GlobalEnv$curMSL[srcCurTaxonIdx,".split"] = TRUE
          .GlobalEnv$curMSL[srcCurTaxonIdx,".split_code"] = code
          .GlobalEnv$curMSL[srcCurTaxonIdx,".split_linenum"] = linenum
          
          # mark the newMSL entry we may delete
          srcNewTaxonIdx=(.GlobalEnv$newMSL$name==as.character(srcTaxonName))
          if( sum(srcNewTaxonIdx)>0) {
            # mark for possible deletion, if not also flagged for preservation
            .GlobalEnv$newMSL[srcNewTaxonIdx,".split"] = TRUE 
            .GlobalEnv$newMSL[srcNewTaxonIdx,".split_code"] = code
            .GlobalEnv$newMSL[srcNewTaxonIdx,".split_linenum"] = linenum
          }
        }
      }
     
      # 
      # check DEST taxon specified
      #
      if(is.na(destTaxonName)) {
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                         "ERROR", "CREATE.NO_DEST", paste0("Change=",toupper(actionClean),", but 'proposed taxonomy' columns are empty"), 
                         paste0("currentVsProposed=", diffLineageStrings(change$.srcLineage, change$.destLineage),
                                ",currentTaxonomy=", change$.srcLineage,
                                ", proposedTaxonomy=",change$.destLineage)
        )
        next;
      }
      
      #
      # check if same accession number already exists
      #
      isDupAccession = (.GlobalEnv$newMSL$genbank_accession_csv == change$exemplarAccession)
      if(sum(isDupAccession, na.rm=TRUE)>0) {
        
        # unless this is a SPLIT doing a rename, but keeping the isolate/accession
        if( actionClean == 'split' &&
            .GlobalEnv$curMSL$.split[srcCurTaxonIdx] &&
            !is.na(.GlobalEnv$curMSL$genbank_accession_csv[srcCurTaxonIdx]) &&
            .GlobalEnv$curMSL$genbank_accession_csv[srcCurTaxonIdx] == change$exemplarAccession) {
  
          # this is ok for a split to re-use an accession under a new or same name, but just once
          if( !is.na(.GlobalEnv$curMSL$.split_acc_used[srcCurTaxonIdx]) ) {
             # can't use it more than once
            errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                             "WARNING", "SPLIT.REUSE_ACC", paste0("Change=",toupper(actionClean),", accession re-used more than once in a split"), 
                             paste0("accession=", change$exemplarAccession,
                                    "; already reused on ", .GlobalEnv$curMSL$.split_acc_used[srcCurTaxonIdx])
            )
            next;
          } # illegal accession reuse
          else {
            # record that a split has used this accession
            .GlobalEnv$curMSL[srcCurTaxonIdx,".split_acc_used"] =paste0(code,":",linenum)
            if(params$tmi) {
              errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                               "INFO", "SPLIT.REUSE_ACC", paste0("Change=",toupper(actionClean),", accession re-used in a split"), 
                               paste0("accession=", change$exemplarAccession,
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
          if(  change$exemplarAccession %in% c("pending","Pending") ) { 
             errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                             errLevel, "CREATE.PENDING_ACC", paste0("Change=",toupper(actionClean),", accession number is 'pending'"), 
                             paste0("accession=", change$exemplarAccession)
            )
          } else {
            # hard error for each other species with this accession
            ## QQQ what proposal created this? (from this round? historically? Need a function: last_modified())
            errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                             errLevel, "CREATE.DUP_ACC", paste0("Change=",toupper(actionClean),", a species with this accession number already exists"), 
                             paste0("accession=", change$exemplarAccession, ", existingSpecies=",.GlobalEnv$newMSL[isDupAccession,]$lineage)
            )
          }     
          # hard error in FINAL mode
          if(params$processing_mode=="final") { next }
        } # acc re-use not in split
         
      } # acc re-use
      
      # 
      # verify that PARENT taxon exists already in newMSL
      #
      if(is.na(destParentName)) {
        # no parent - use root node
        parentDestNewMatches=(.GlobalEnv$newMSL$taxnode_id==.GlobalEnv$newMSL$tree_id)
      } else {
        # find parent by name (must be unique)
        parentDestNewMatches=(.GlobalEnv$newMSL$name==as.character(destParentName))
     }
 
      if(params$verbose) {print(paste0("CREATE: ",code," line ",linenum," '",destTaxonName, "' findParent(",destParentName,")=",sum(parentDestNewMatches)))}
    
      if(sum(parentDestNewMatches)==0) {
        # check if it used to exist
        prevDestParent = curMSL %>% filter(name==as.character(destParentName))
        if(nrow(prevDestParent)==0) {
          # just missing
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                            "ERROR", "CREATE.PARENT_NO_EXIST", paste0("Change=",toupper(actionClean),", but parent rank taxon does not exist"), 
                            paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage)
           )
           next;
        } else {
          # someone already modified it - tell me who!
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "WARNING", "CREATE.PARENT_ALREADY_CHANGED", paste0("Change=",toupper(actionClean),", but parent taxon already modified"), 
                           paste0("proposal(s)=", prevDestParent$out_filename, " did a '",prevDestParent$out_change,"' ",
                                  "of '",prevDestParent$name, "' to '", prevDestParent$out_target,"'")
          )
          # find what it became
          parentDestNewMatches=(.GlobalEnv$newMSL$lineage==as.character(prevDestParent$out_target)|.GlobalEnv$newMSL$name==as.character(prevDestParent$out_target))
          if(sum(parentDestNewMatches) != 1) {
            # couldn't find new version, what the heck?!!?
            errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                             "ERROR", "CREATE.PARENT_CHANGED_BUT_MISSING", paste0("Change=",toupper(actionClean),", but parent taxon already modified"), 
                             paste0("can't find  '", prevDestParent$out_target,"' in new MSL: did someone else further change it?")
            )
            next;
          }
          
                 }
      } else if(sum(parentDestNewMatches)>1) {
         # this should never happen
         errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                          "ERROR", "CREATE.PARENT_MANY", paste0("Change=",toupper(actionClean),", multiple taxa exist with parent name"), 
                          paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage)
         )
         next;
      }
      # WARN: check if taxon rank matches change rank
      if ( tolower(destTaxonRank) != tolower(change$rank) ) {
        errorDf=addError(errorDf,code,linenum,change$change,change$rank,change$.destTaxon,
                         "WARNING", "CREATE.RANK_MISMATCH", 
                         paste0("Change=",toupper(actionClean),", but proposed taxon rank does not match [rank] column"), 
                         paste0("ankColumn=", change$rank,
                                ", proposedTaxonRank=",destTaxonRank,
                                ", proposedTaxonomy=", destLineage)
        )
      }
      
      #### binomial check ####
      # if new SPECIES, check binomial prefix matches parent genus
      #
      if( destTaxonRank=="species" ) {
        
        # 
        # check hostSource
        #
        if( is.na(change$hostSource) || grepl("please.*select",change$hostSource, ignore.case=T) ) {
          errorDf=addError(errorDf,code,linenum,change$change,change$rank,change$.destTaxon,
                           "WARNING", "CREATE.SPECIES_NO_HOST_SOURCE", 
                           paste0("Change=",toupper(actionClean),", but proposed species must have a host/source value"), 
                           ""
          )
        }
        
        # check if parent is a genus
        destParentTaxon = .GlobalEnv$newMSL[parentDestNewMatches,]
        if(destParentTaxon$rank=="genus") {
          destParentGenusTaxon = destParentTaxon
        } else {
          if(destParentTaxon$rank == "subgenus" ) {
            # check up one rank
            destParentGenusTaxon = .GlobalEnv$newMSL %>% filter(taxnode_id == destParentTaxon$parent_id)
            if(destParentGenusTaxon$rank != "genus") {
              # can't find genus parent
              errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                               "ERROR", "CREATE.SPECIES_SUBGENUS_NO_GENUS", 
                               paste0("Change=",toupper(actionClean),", can not find parent genus"), 
                               paste0("subgenus '",destParentTaxon$name,"' is not in a genus;",
                                      " it's parent '",destParentGenusTaxon$name,"' is a ",destParentGenusTaxon$rank )
              )
              next;
            }
          } else {
            # parent isn't subgenus or genus
            errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                             "ERROR", "CREATE.SPECIES_NO_GENUS", 
                             paste0("Change=",toupper(actionClean),", can not find parent genus"), 
                             paste0("parent '",destParentTaxon$name,"' is a ",destParentTaxon$rank )
            )
            next;
          }
        } 
        # have genus (parent or grandparent)
        # check binomial naming
        if( str_detect(destTaxonName,paste0(destParentGenusTaxon$name," ")) != TRUE ) {
          errorDf=addError(errorDf,code,linenum,change$change,change$rank,change$.destTaxon,
                           "ERROR", "CREATE.SPECIES_BINOMIAL_MISMATCH", 
                           paste0("Change=",toupper(actionClean),", but proposed species names does not start with 'genus[space]' per binomial naming convention"), 
                           paste0("parent genus name =", destParentGenusTaxon$name)
          )
          next;
        }
        
      }
      
      # 
      #####  create new taxon  ##### 
      #
      
      
      # get parent
      newTaxon = .GlobalEnv$newMSL[parentDestNewMatches,]
      
      # WARN if PARENT_LINEAGE is not expected, AND USE PARENT LINEAGE
      if( !is.na(destParentLineage) && (newTaxon$lineage != destParentLineage) ) {
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                         "WARNING", "CREATE.PARENT_LINEAGE", 
                         paste0("Change=",toupper(actionClean),", proposed parent taxon exists, but not with expected lineage, using observed lineage"), 
                         paste0("proposedVsObserved=", diffLineageStrings(destParentLineage, newTaxon$lineage),
                                ", proposedParentLineage=", destParentLineage,
                                ", observedParentLineage=",newTaxon$lineage,
                                ", otherProposals=",newTaxon$prev_proposals)
        )
        destLineage=destParentLineage
      }
      
      # add new info - primary columns
      newTaxon[1,"in_change"]   = actionClean
      newTaxon[1,"in_filename"] = proposalZip
      newTaxon[1,"in_notes"]    = paste0("xlsx_row=",linenum)
      newTaxon[1,"in_target"]   = destLineage
      
      newTaxon[1,"name"]        = destTaxonName
      newTaxon[1,"cleaned_name"]= destTaxonName
      newTaxon[1,"level_id"]    = rankCV$id[rankCV$name==destTaxonRank]
      newTaxon[1,"rank"]        = destTaxonRank
      newTaxon[1,"parent_id"]   = newTaxon[1,"taxnode_id"]
      newTaxon[1,"taxnode_id"]  = max(.GlobalEnv$newMSL$taxnode_id)+1
      newTaxon[1,"ictv_id"]     = newTaxon$taxnode_id
      
      # genomeComposition = molecule_id 
      newTaxon[1,xlsx2dbMap["molecule"]] = ifelse(is.na(change$molecule),NA,dbCvMapList[["molecule"]][change$molecule])
      
      # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
      newTaxon[1,xlsx2dbMap["hostSource"]] = change[1,"hostSource"] 
      
      # comments
      newTaxon[1,xlsx2dbMap["comments"]]= change[1,"comments"]
      
      
      #
      #### CREATE.SPECIES ####  
      #
      if( destTaxonRank == "species" ) {
        #
        ## for species only
        #
        
        # "genbank_accession_csv"
        newTaxon[1,xlsx2dbMap["exemplarAccession"]] = change[1,"exemplarAccession"] 
        # exemplar_name
        # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
        newTaxon[1,xlsx2dbMap["exemplarName"]] = change[1,"exemplarName"] 
        # "abbrev_csv"
        newTaxon[1,xlsx2dbMap["Abbrev"]] = change[1,"Abbrev"] 
        # "isolate_csv"
        newTaxon[1,xlsx2dbMap["exemplarIsolate"]] = change[1,"exemplarIsolate"] 
        # genome_coverage
        # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
        newTaxon[1,xlsx2dbMap["genomeCoverage"]]= change[1,"genomeCoverage"] 
        # genome_coverage
        # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
        newTaxon[1,xlsx2dbMap["genomeCoverage"]]= change[1,"genomeCoverage"] 
        # molecule_id is above (not species-specific)
      }
      
      # info-only columns - wont be saved to DB
      newTaxon[1,"rank"]       = destTaxonRank
      newTaxon[1,"lineage"]    = destLineage
      
      # clear some columns inherited from parent
      newTaxon[1,"prev_taxnode_id"] = NA
      newTaxon[1,"prev_proposals"] = paste0(code,":",linenum) 
      
      # add new taxon to newMSL
      if(params$verbose) {print(paste0("rbindlist(newMSL,",newTaxon[1,"taxnode_id"],":",destLineage,")"))}
      .GlobalEnv$newMSL <- rbindlist(list(.GlobalEnv$newMSL,newTaxon),fill=TRUE)
      
      # SUCCESS message
      errorCodeCreateMap=c("new"="CREATE.OK","split"="SPLIT.OK")
      errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                       "SUCCESS", errorCodeCreateMap[actionClean], 
                       paste0("Change=",toupper(actionClean),", applied successfully"), 
                       paste0("Create ", newTaxon$rank," of '",newTaxon$lineage,"'")
      )
    } # create new taxon
    else  if(actionClean %in% c("rename") ) {
    
      #  -------------------------------------------------------------------------
      #
      #### RENAME ####
      #
      #  -------------------------------------------------------------------------

      # check if srcTaxon was specified in xlsx (required)
      if(is.na(srcTaxonName)) {
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "ERROR", "RENAME.WO_SRC", "Change=RENAME, but 'current taxonomy' columns are empty", 
                         paste0("currentTaxonomy=", srcLineage, ", destTaxonomy=", destLineage)
        )
        next;
      }
      
      # check if renamed taxon has the same name (TP:Error-R1)
      if(srcTaxonName == destTaxonName) {
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "WARNING", "RENAME.SAME_NAME", "Change=RENAME, but new name is the same; might this be a move?", 
                         paste0("currentVsProposed=",diffLineageStrings(srcLineage, destLineage),
                                ", currentTaxonomy=", srcLineage, ", destTaxonomy=", destLineage)
        )
      }
      
      
      #
      # find the target to rename in NEW
      #
      srcNewTarget=(.GlobalEnv$newMSL$name==as.character(srcTaxonName))
      if(params$verbose) {print(paste0("RENAME: ",code," line ",linenum," '",destTaxonName, "' findTarget(",srcTaxonName,")=",sum(srcNewTarget)))}
    
      if(sum(srcNewTarget)==0) {
        # check if someone else already modified it
        prevSrcTaxon = curMSL %>% filter(name==as.character(srcTaxonName))
        #browser()
        if(nrow(prevSrcTaxon)==0) {
          # just not found
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                            "ERROR", "RENAME.NO_EXIST", "Change=RENAME, but taxon does not exist", 
                            paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName))
           next;
        } else {
          # previous record - must have been modified already this MSL
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "ERROR", "RENAME.SRC_ALREADY_CHANGED", "Change=RENAME, but taxon already modified", 
                           paste0("proposal(s)=", prevSrcTaxon$out_filename, " also did a '",prevSrcTaxon$out_change,"' to '", prevSrcTaxon$out_target,"'")
          )
          next;
        }
      } else if(sum(srcNewTarget)>1) {
         errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                          "ERROR", "RENAME.MANY", "Change=RENAME, multiple taxa exist with parent name", 
                          paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName
                                , " matches ", sum(srcNewTarget) ))
         next;
      } else {
        
        # find original taxon in prevMSL being renamed
        srcPrevTarget=(.GlobalEnv$curMSL$taxnode_id==.GlobalEnv$newMSL[srcNewTarget,]$prev_taxnode_id)
        srcNewParent= (.GlobalEnv$newMSL$taxnode_id==.GlobalEnv$newMSL[srcNewTarget,]$parent_id )
        # print(paste("srcPrevTarget=",sum(srcPrevTarget),"srcNewParent=",sum(srcNewParent)))
        
        ##### binomial check #####
        
        # if rename SPECIES, check binomial prefix matches parent genus
        #
        if( destTaxonRank=="species" ) {
          # check if parent is a genus
          destParentTaxon = as.data.frame(newMSL[srcNewParent,])
          if(destParentTaxon$rank=="genus") {
            destParentGenusTaxon = destParentTaxon
          } else {
            if(destParentTaxon$rank == "subgenus" ) {
              # check up one rank
              destParentGenusTaxon = .GlobalEnv$newMSL %>% filter(taxnode_id == destParentTaxon$parent_id)
              if(destParentGenusTaxon$rank != "genus") {
                # can't find genus parent
                errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                                 "ERROR", "RENAME.SPECIES_SUBGENUS_NO_GENUS", 
                                 "Change=RENAME, can not find parent genus", 
                                 paste0("subgenus '",destParentTaxon$name,"' is not in a genus;",
                                        " it's parent '",destParentGenusTaxon$name,"' is a ",destParentGenusTaxon$rank )
                )
                next;
              }
            } else {
              # parent isn't subgenus or genus
              errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                               "ERROR", "CREATE.SPECIES_NO_GENUS", 
                               "Change=CREATE, can not find parent genus", 
                               paste0("parent '",destParentTaxon$name,"' is a ",destParentTaxon$rank )
              )
              next;
            }
          } 
          # have genus (parent or grandparent)
          # check binomial naming
          if( str_detect(destTaxonName,paste0(destParentGenusTaxon$name," ")) != TRUE ) {
            errorDf=addError(errorDf,code,linenum,change$change,change$rank,change$.destTaxon,
                             "ERROR", "RENAME.SPECIES_BINOMIAL_MISMATCH", 
                             "Change=RENAME, but proposed species names does not start with 'genus[space]' per binomial naming convention", 
                             paste0("parent genus name =", destParentGenusTaxon$name)
            )
            next;
          }
          
        }
        
        ##### apply changes #####
        # change the name and update lineage
        .GlobalEnv$newMSL[srcNewTarget,"name"] = destTaxonName
        .GlobalEnv$newMSL[srcNewTarget,"lineage"] = paste(newMSL[srcNewParent,"lineage"],destTaxonName,sep=";") # should we do this? 
        
        # add to list of genera to check for binomial problems after proposal is complete
       if(  destTaxonRank == "genus" ) {
          renamedGenera[destTaxonName$genus] = linenum
        } 
        
        # check if that is the expected lineage
        # WARN that PARENT_LINEAGE is not expected, AND USE MSL/OBSERVED PARENT LINEAGE
        if( change$.destParentName != .GlobalEnv$newMSL[srcNewParent,"name"] ) {
          parentRank = .GlobalEnv$newMSL[srcNewParent,]$rank
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "WARNING", "RENAME.PARENT_NAME", 
                           "parent name not as expected; was this meant to be a move? Is there a typo?", 
                           paste0("parent ", parentRank, " in this proposal: '", change$.destParentName,"', ",
                                  "but MSL parent is ", parentRank, " '",.GlobalEnv$newMSL[srcNewParent,"name"],"'; ",
                                  "other new proposals: ",.GlobalEnv$newMSL[srcNewTarget,"prev_proposals"])
          )
         }
 
        # append proposal
        .GlobalEnv$newMSL[srcNewTarget,"prev_proposals"] = paste0(
          ifelse(is.na(.GlobalEnv$newMSL[srcNewTarget,"prev_proposals"]),"",paste0(.GlobalEnv$newMSL[srcNewTarget,"prev_proposals"],",")),
          proposalZip,":",linenum)
        
        # put out_* changes on curMSL
        .GlobalEnv$curMSL[srcPrevTarget,"out_updated"] = TRUE  # admin; mark this to save to db
        .GlobalEnv$curMSL[srcPrevTarget,"out_change"] = "rename"
        .GlobalEnv$curMSL[srcPrevTarget,"out_filename"] = proposalZip
        .GlobalEnv$curMSL[srcPrevTarget,"out_target"] = destTaxonName
        .GlobalEnv$curMSL[srcPrevTarget,"out_notes"] = paste0("linenum=",linenum)
        
        # success note
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "SUCCESS", "RENAME.OK", "Change=RENAME, applied successfully", 
                         paste0("RENAME ",srcTaxonRank, 
                                " from ", "'",srcTaxonName,"'",
                                " to ",   "'", .GlobalEnv$newMSL[srcNewTarget,"name"], "'",
                                " in ",   "'", .GlobalEnv$newMSL[srcNewTarget,"lineage"], "'"
                         )
        )
        
        ## ZZZ add comments to out_notes?
      } # target found
      # RENAME
    } else if(actionClean %in% c("abolish") ) { 
      #  -------------------------------------------------------------------------
      #
      #####  ABOLISH       #### 
      #
      #  -------------------------------------------------------------------------
      # check if srcTaxon was specified in xlsx (required)
      if(is.na(srcTaxonName)) {
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "ERROR", "ABOLISH.WO_SRC", "Change=ABOLISH, but 'current taxonomy' columns are empty", 
                         paste0("destTaxonomy=", destLineage)
        )
        next;
      }
      
      #
      # find the target to abolish in NEW and CUR
      #
      srcNewTarget=(.GlobalEnv$newMSL$name==as.character(srcTaxonName))
      srcPrevTarget=(.GlobalEnv$curMSL$taxnode_id==.GlobalEnv$newMSL[srcNewTarget,]$prev_taxnode_id)

      if(params$verbose) {print(paste0("ABOLISH: ",code," line ",linenum," '",srcTaxonName, "' findTarget(",srcTaxonName,")=",sum(srcNewTarget),"/",sum(srcPrevTarget)))}
      
      if(sum(srcNewTarget,na.rm=TRUE)==0) {
        # check if someone else already modified it
        prevSrcTaxon = curMSL %>% filter(name==as.character(srcTaxonName))
        #browser()
        if(nrow(prevSrcTaxon)==0 ) {
          # not modified, just missing
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                          "ERROR", "ABOLISH.NO_EXIST", "Change=ABOLISH, but taxon does not exist", 
                          paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName))
          next;
        } else {
          # previous record - must have been modified already this MSL
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "ERROR", "ABOLISH.SRC_ALREADY_CHANGED", "Change=ABOLISH, but taxon already modified", 
                           paste0("proposal(s)=", prevSrcTaxon$out_filename, " also did a '",prevSrcTaxon$out_change,"' to '", prevSrcTaxon$out_target,"'")
          )
          next;
        }
      } else if(sum(srcNewTarget,na.rm=TRUE)>1) {
         # this should never happen
         errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                          "ERROR", "ABOLISH.MANY", "Change=ABOLISH, multiple taxa exist with name", 
                          paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName
                                , " matches ", sum(srcNewTarget) ))
         next;
      }
      #
      # check for kids - can't abolish before kids
      #
      srcKids=(.GlobalEnv$newMSL$parent_id==.GlobalEnv$newMSL[srcNewTarget,"taxnode_id"])
      if(sum(srcKids,na.rm=TRUE)>0) {
         errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                          "ERROR", "ABOLISH.KIDS", "Change=ABOLISH, taxon still has un-abolished children", 
                         paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", kids: N=",sum(srcKids,na.rm=TRUE),
                                ", NAMES=[", join(.GlobalEnv$newMSL[srcKids,"rank"],.GlobalEnv$newMSL[srcKids,"name"],sep=":") ) )
         next;
      }     
      
      #
      # remove taxon from NEW
      #
      .GlobalEnv$newMSL = .GlobalEnv$newMSL[!srcNewTarget,]
      
      # put out_* changes on curMSL
      .GlobalEnv$curMSL[srcPrevTarget,"out_updated"] = TRUE  # admin; mark this to save to db
      .GlobalEnv$curMSL[srcPrevTarget,"out_change"] = "abolish"
      .GlobalEnv$curMSL[srcPrevTarget,"out_filename"] = proposalZip
      .GlobalEnv$curMSL[srcPrevTarget,"out_target"] = destTaxonName
      .GlobalEnv$curMSL[srcPrevTarget,"out_notes"] = paste0("linenum=",linenum) # add comments?
      
      # success note
      errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                       "SUCCESS", "ABOLISH.OK", "Change=ABOLISH, applied successfully", 
                       paste0("ABOLISH ",srcTaxonRank, " named ",srcTaxonName)
      )
      
    } else if(actionClean %in% c("move","promote","demote","merge") 
              || ( actionClean =='split' && srcTaxonName == destTaxonName) ) { 
      #  -------------------------------------------------------------------------
      #
      #### MOVE, PROMOTE/DEMOTE, MERGE, SPLIT #### 
      #
      # Only SPLIT where srcName = destName
      # MOVE w/ and w/o rename
      #
      # these all require a SRC and DEST, and change structure
      #
      #  -------------------------------------------------------------------------
      # check if srcTaxon was specified in xlsx (required)
      if(is.na(srcTaxonName) && actionClean =='move') {
        # they didn't specify the source, but  it's a MOVE
        # we should be able to find it
        if(actionClean %in% c("move") ) {
          guessSource = (.GlobalEnv$newMSL$name==as.character(destTaxonName))
          if( sum(guessSource,na.rm=TRUE) == 1) {
            # found it, so WARNING
            errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                             "WARNING", "MOVE.NO_SRC_CAN_GUESS", "Change=MOVE, but 'current taxonomy' columns are empty", 
                             paste0("destTaxonomy=", destLineage,
                                    "; we guess you meant=", .GlobalEnv$newMSL$lineage[guessSource])
            )
            srcTaxonName=.GlobalEnv$newMSL$name[guessSource]
            srcTaxonRank=.GlobalEnv$newMSL$rank[guessSource]
            srcLineage  =.GlobalEnv$newMSL$lineage[guessSource]
          } else {
            # didn't find it, or found multiple, so ERROR
            errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                             "WARNING", "MOVE.NO_SRC", "Change=MOVE, but 'current taxonomy' columns are empty", 
                             paste0("destTaxonomy=", destLineage,
                                    "; no existing taxon '",destTaxonName,"' found")
            )
            next
          }
        } else { 
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "WARNING", "CHANGE_RANK.NO_SRC", 
                           paste0("Change=",toupper(actionClean),", but 'current taxonomy' columns are empty"), 
                           paste0("destTaxonomy=", destLineage)
          )
          next
          }
      }
      
      #
      ##### validate rank #####
      #   1. consistent
      #   2. move does not change
      #   3. promote/demote does
      #
      if( destTaxonRank != rankClean ) {
          # 1. consistent: src rank doesn't match action rank
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "ERROR", "MOVE.DEST_WRONG_RANK", 
                           paste0("Change=",toupper(actionClean),", but 'proposed taxonomy' ranks doesn't match 'proposed rank' column"), 
                           paste0("rank=", rankClean, ", srcTaxonomy=[", srcTaxonRank,"]",srcLineage)
          )
      }
      if( actionClean %in% c("move") && srcTaxonRank != destTaxonRank ) {
        # 2. MOVE: src/dest ranks NOT the same
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "ERROR", "MOVE.DIFF_RANK", 
                         "Change=MOVE, but 'current taxonomy' and 'proposed taxonomy' ranks differ; use promote or demote, not move, to change rank", 
                         paste0("srcVsDest=[",srcTaxonRank,"//",destTaxonRank,"] ", diffLineageStrings(srcLineage,destLineage),
                                ", srcTaxonomy=[", srcTaxonRank,"]",srcLineage,
                                ", destTaxonomy=[", destTaxonRank,"]", destLineage)
        )
        next;
      } else if( actionClean %in% c("promote","demote") && srcTaxonRank == destTaxonRank ) {
        # 3. PROMOTE/DEMOTE src/dest ranks ARE the same 
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "ERROR", "CHANGE_RANK.SAME_RANK", 
                         paste0("Change=",toupper(actionClean),", but 'current taxonomy' and 'proposed taxonomy' ranks are the same; use MOVE"), 
                         paste0("rank=", rankClean, 
                                ", srcVsDest=",diffLineageStrings(srcLineage,destLineage),
                                ", srcTaxonomy=[", srcTaxonRank,"]", srcLineage,
                                ", destTaxonomy=[", destTaxonRank,"]", destLineage)
        )
      }
      
      #
      # find the SRC target to MOVE in NEW and CUR
      #
      srcNewTarget =(.GlobalEnv$newMSL$name==as.character(srcTaxonName))
      srcPrevTarget=(.GlobalEnv$curMSL$taxnode_id==.GlobalEnv$newMSL[srcNewTarget,]$prev_taxnode_id)

      if(params$verbose) {print(paste0("MOVE: ",code," line ",linenum," '",srcTaxonName, "' findTarget(",srcTaxonName,")=",sum(srcNewTarget),"/",sum(srcPrevTarget)))}
    
      if(sum(srcNewTarget,na.rm=TRUE)==0) {
        # check prevMSL, to see if something else already moved it
        prevSrcTaxon = curMSL %>% filter(name==as.character(srcTaxonName))
        #browser()
        if(nrow(prevSrcTaxon)==0 ) {
          # no previous record, just doesn't exist
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                            "ERROR", "MOVE.NO_EXIST", paste0("Change=",toupper(actionClean),", but taxon does not exist"), 
                           paste0("taxon=", srcTaxonName, ",srcVsDest=",diffLineageStrings(srcLineage,destLineage),
                                  ", lineage=",srcLineage,", proposedLineage=", destLineage))
          next;
        } else {
          # previous record - must have been modified already this MSL
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "ERROR", "MOVE.SRC_ALREADY_CHANGED", paste0("Change=",toupper(actionClean),", but taxon already modified"), 
                           paste0("proposal(s)=", prevSrcTaxon$out_filename, " also did a '",prevSrcTaxon$out_change,"' to '", prevSrcTaxon$out_target,"'")
                           )
          next;
        }
      } else if(sum(srcNewTarget,na.rm=TRUE)>1) {
         errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                          "ERROR", "MOVE.MANY", paste0("Change=",toupper(actionClean),", multiple taxa exist with name"), 
                         paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName
                                , " matches ", sum(srcNewTarget) ))
         next;
      }

      #
      # check if same accession number already exists
      #
      #  !srcNewTarget - prevents checking against ourself
      #
      isDupAccession = (.GlobalEnv$newMSL$genbank_accession_csv == change$exemplarAccession & !srcNewTarget)
      if(sum(isDupAccession, na.rm=TRUE)>0) {
        errorDf=addError(errorDf,code,linenum,change$change,change$rank,change$.changeTaxon,
                         "ERROR", "MOVE.DUP_ACC", paste0("Change=",toupper(actionClean),", a species with this accession number already exists"), 
                         paste0("accession=", change$exemplarAccession, ", existingSpecies=",.GlobalEnv$newMSL[isDupAccession,]$lineage)
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
      } else {
        # find parent by name (must be unique)
        parentDestNewMatches=(.GlobalEnv$newMSL$name==as.character(destParentName))
      }
      
      if(params$verbose) {print(paste0("MOVE: ",code," line ",linenum," '",destTaxonName, "' findParent(",destParentName,")=",sum(parentDestNewMatches)))}
    
      if(sum(parentDestNewMatches,na.rm=TRUE)==0) {
        # check if it used to exist
        prevDestParent = curMSL %>% filter(name==as.character(destParentName))
        if(nrow(prevDestParent)==0) {
          # just missing
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "ERROR", "MOVE.PARENT_NO_EXIST", paste0("Change=",toupper(actionClean),", but parent rank taxon does not exist"), 
                           paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage)
          )
          next;
        } else {
          # someone already modified it - tell me who!
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "ERROR", "MOVE.PARENT_ALREADY_CHANGED", paste0("Change=",toupper(actionClean),", but parent taxon already modified"), 
                           paste0("proposal(s)=", prevDestParent$out_filename, " also did a '",prevDestParent$out_change,"' to '", prevDestParent$out_target,"'")
          )
          next;
        }
      } else if(sum(parentDestNewMatches,na.rm=TRUE)>1) {
        # this should never happen
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "ERROR", "MOVE.PARENT_MANY", paste0("Change=",toupper(actionClean),", multiple taxa exist with parent name"), 
                         paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage)
        )
        next;
      } 
      
      # remember to check later if a genus name change broke binomial
      if(  destTaxonRank == "genus" ) {
        renamedGenera[destTaxonName$genus] = linenum
      } 
       
      #
      # check that promote/demote change name
      #
      if( actionClean %in% c("promote","demote") ) {
          if(srcTaxonName == destTaxonName) {
            errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                             "WARNING", "CHANGE_RANK.SAME_NAME", paste0("Change=",toupper(actionClean),", but current and proposed names are the same."), 
                             paste0("current name=[",srcTaxonRank,"]", srcTaxonNameprevDestParent$out_filename, " and proposed name [",destTaxonRank,"]",destTaxonName)
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
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "WARNING", "MOVE.PARENT_LINEAGE", 
                         paste0("Change=",toupper(actionClean),", proposed parent taxon exists, but not with expected lineage"), 
                         paste0(",proposedVsObserved=",diffLineageStrings(destParentLineage,destParentTaxon$lineage),
                                ", proposedParentLineage=", destParentLineage,
                                ", observedParentLineage=",destParentTaxon$lineage,
                                ", otherProposals=",destParentTaxon$prev_proposals)
        )
      }
      
      #### binomial check (species) ####
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
              errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                               "ERROR", "MOVE.SPECIES_SUBGENUS_NO_GENUS", 
                               paste0("Change=",toupper(actionClean),", can not find parent genus"), 
                               paste0("subgenus '",destParentTaxon$name,"' is not in a genus;",
                                      " it's parent '",destParentGenusTaxon$name,"' is a ",destParentGenusTaxon$rank )
              )
              next;
            }
          } else {
            # parent isn't subgenus or genus
            errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                             "ERROR", "MOVE.SPECIES_NO_GENUS", 
                             paste0("Change=",toupper(actionClean),", can not find parent genus"), 
                             paste0("parent '",destParentTaxon$name,"' is a ",destParentTaxon$rank )
            )
            next;
          }
        } 
        # have genus (parent or grandparent)
        # check binomial naming
        if( str_detect(destTaxonName,paste0(destParentGenusTaxon$name," ")) != TRUE ) {
          errorDf=addError(errorDf,code,linenum,change$change,change$rank,change$.destTaxon,
                           "ERROR", "MOVE.SPECIES_BINOMIAL_MISMATCH", 
                           paste0("Change=",toupper(actionClean),", but proposed species names does not start with 'genus[space]' per binomial naming convention"), 
                           paste0("parent genus name =", destParentGenusTaxon$name)
          )
          next;
        }
      }
      
      
      # add new info - primary columns
      .GlobalEnv$newMSL[srcNewTarget,"name"]       = destTaxonName
      # this part shouldn't be a change unless this is promote/demote....
      .GlobalEnv$newMSL[srcNewTarget,"level_id"]   = rankCV$id[rankCV$name==destTaxonRank]
      .GlobalEnv$newMSL[srcNewTarget,"rank"]       = destTaxonRank
      .GlobalEnv$newMSL[srcNewTarget,"parent_id"]  = destParentTaxon[1,"taxnode_id"]
      
      
      # genomeComposition = molecule_id 
      if(!is.na(change[1,"molecule"])) {
        .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["molecule"]] = dbCvMapList[["molecule"]][change[1,"molecule"] ]
      }
      
      # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
      if(!is.na(change[1,"hostSource"])) {
        .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["hostSource"]] = change[1,"hostSource"] 
      }
      # comments
      if(!is.na(change[1,"comments"])) {
        .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["comments"]]= change[1,"comments"]
      }
      
      ## for species only
      if( destTaxonRank == "species" ) {
        # "genbank_accession_csv"
        if(!is.na(change[1,"exemplarAccession"])) {
          .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["exemplarAccession"]] = change[1,"exemplarAccession"] 
        }
        # exemplar_name
        # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
        if(!is.na(change[1,"exemplarName"])) {
          .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["exemplarName"]] = change[1,"exemplarName"] 
        }
        # "abbrev_csv"
        if(!is.na(change[1,"Abbrev"])) {
          .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["Abbrev"]] = change[1,"Abbrev"] 
        }
        # "isolate_csv"
        if(!is.na(change[1,"exemplarIsolate"])) {
          .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["exemplarIsolate"]] = change[1,"exemplarIsolate"] 
        }
        # genome_coverage
        # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
        if(!is.na(change[1,"genomeCoverage"])) {
          .GlobalEnv$newMSL[srcNewTarget,xlsx2dbMap["genomeCoverage"]]= change[1,"genomeCoverage"] 
        }
      }
      
      #  RECURSE TO SET LINEAGE OF destPARENT's KIDS
      update_lineage(destParentTaxon$taxnode_id,destParentTaxon$lineage)
      
      #
      # update curMSL [out_*] (for non-splits)
      #
      if( actionClean == 'split' ) {
        #
        # if the list of split directives does NOT include the current name, 
        # then we will have to (after the proposal is finished) abolish the copy
        # of the current name in the new MSL. Here, we mark that we've seen the 
        # split directive with the same name, so we should NOT abolish this one.
        # admin; mark we kept the original name in the split
        .GlobalEnv$curMSL[srcPrevTarget,".split_kept"] = TRUE # old copy, just to know
        .GlobalEnv$newMSL[srcNewTarget,".split_kept"] = TRUE  # new copy - don't delete!
        .GlobalEnv$newMSL[srcNewTarget,".split_code"] = code  
        .GlobalEnv$newMSL[srcNewTarget,".split_linenum"] = linenum  # new copy - mark which line saved
        
        # set IN change for SPLIT
        .GlobalEnv$newMSL[srcNewTarget,"in_change"] = actionClean
        .GlobalEnv$newMSL[srcNewTarget,"in_filename"] = proposalZip
        .GlobalEnv$newMSL[srcNewTarget,"in_target"] = srcLineage
        .GlobalEnv$newMSL[srcNewTarget,"in_notes"] = paste0("linenum=",linenum) # add comments?
        
      } else {
        # set OUT change for all others
        .GlobalEnv$curMSL[srcPrevTarget,"out_updated"] = TRUE  # admin; mark this to save to db
        .GlobalEnv$curMSL[srcPrevTarget,"out_change"] = actionClean
        .GlobalEnv$curMSL[srcPrevTarget,"out_filename"] = proposalZip
        .GlobalEnv$curMSL[srcPrevTarget,"out_target"] = destLineage
        .GlobalEnv$curMSL[srcPrevTarget,"out_notes"] = paste0("linenum=",linenum) # add comments?
      }
      
      # success note
      errorCodeMoveMap=c("move"="MOVE.OK","promote"="PROMOTE.OK","demote"="DEMOTE.OK","merge"="MERGE.OK","split"="SPLIT=.OK")
      errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                       "SUCCESS",  errorCodeMoveMap[actionClean], 
                       paste0("Change=",toupper(actionClean),", applied successfully"), 
                       paste0(toupper(actionClean)," ",
                              .GlobalEnv$curMSL$rank[srcPrevTarget], " named '", .GlobalEnv$curMSL[srcPrevTarget,"name"],    "'", 
                              " to ", .GlobalEnv$newMSL$rank[srcNewTarget], " named '", .GlobalEnv$newMSL[srcNewTarget,"name"],    "'", 
                              " change=",  diffLineageStrings(.GlobalEnv$curMSL[srcPrevTarget,"lineage"],.GlobalEnv$newMSL[srcNewTarget, "lineage"])
                              
                       )
      )
    } else {
      #  -------------------------------------------------------------------------
      #
      #### UNKNOWN ACTION not implemented. ####
      #
      #  -------------------------------------------------------------------------
      errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                       "ERROR", "CHANGE.UNIMP", 
                       paste0("Change=",toupper(action)," is NOT (yet) implemented"), 
                       paste0("lineageChange=", diffLineageString(srcLineage,destLineage))
                       )
      next;
  
    } # ACTION UN-IMP
 
  } # for each line in XLSX

    #  -------------------------------------------------------------------------
  #
  #### POST-QC #### 
  #
  #  -------------------------------------------------------------------------
  
  #
  ##### empty taxa #####
  # did proposal create empty (non-species) taxa? (global scan)
  #
  kidCounts = apply(.GlobalEnv$newMSL[,"taxnode_id"],1,
                    function(parent_id,MSL) { sum(MSL$parent_id==parent_id)},
                    MSL=.GlobalEnv$newMSL)
  # scan for taxa with no kids, not species, not yet reported.
  emptyTaxa = .GlobalEnv$newMSL[kidCounts == 0,]%>% filter(level_id != 600) %>% filter(is.na(.emptyReported))
  if( nrow(emptyTaxa) > 0 ) {
    errorDf=addError(errorDf,code,NA, NA,emptyTaxa$rank,emptyTaxa$name,
                     "ERROR", "PROPOSAL.EMPTY_TAXA", 
                     paste0("Proposal created empty (non-species) taxa"), 
                     paste0("the ",emptyTaxa$rank," '",emptyTaxa$name,"' is empty - it does not contain any lower rank taxons")
    )
    # mark empty taxa so we don't re-report them
    .GlobalEnv$newMSL[.GlobalEnv$newMSL$name %in% emptyTaxa$name,".emptyReported"] = xlsxs[code,"xlsx"]
  }
  
  ##### renamed genera - binomial #####
  for( destGenusName in names(renamedGenera) )  {
    renamedGeneraTaxIds = .GlobalEnv$newMSL %>% filter(name == destGenusName )
    # subgenera
    renamedGeneraSubgeneraTaxIds = .GlobalEnv$newMSL %>% 
      filter(parent_id %in% c(renamedGeneraTaxIds$taxnode_id)) %>%
      filter(level_id == 550 )
    # species
    renamedGeneraSpecies = .GlobalEnv$newMSL %>% 
      filter(parent_id %in% c(renamedGeneraTaxIds$taxnode_id,renamedGeneraSubgeneraTaxIds$taxnode_id)) %>%
      filter(level_id==600)
    
    # check names
    binomialViolations =str_detect(renamedGeneraSpecies$name,paste0("^",destGenusName," "), negate=TRUE)
    if( sum(binomialViolations) > 0 ) {
      # get offending species list
      speciesDf = as.data.frame(renamedGeneraSpecies)[binomialViolations,]
      
      errorDf=addError(errorDf,code,renamedGenera[destGenusName],"rename_genus","species",speciesDf$name,
                       "ERROR", "RENAME_GENUS.SPECIES_BINOMIAL_MISMATCH", 
                       "Change=RENAME_GENUS, but species name does not start with 'genus[space]' per binomial naming convention", 
                       paste0("new genus name: ", destGenusName )
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
    # log that we're removing them
    # function(errorDf,code,row,change,rank,taxon,levelStr,errorCode,errorStr,notes
    errorDf=addError(errorDf,.GlobalEnv$newMSL$.split_code[splitDeleteIdx],.GlobalEnv$newMSL$.split_linenum[splitDeleteIdx],
                     "split_abolish",as.character(.GlobalEnv$newMSL$rank[splitDeleteIdx]), .GlobalEnv$newMSL$name[splitDeleteIdx],
                     "INFO", "SPLIT.IMPLICIT_ABOLISH", 
                     "Change=SPLIT, but no split line kept original name, so remove original name from new MSL", 
                     paste0("removed ",.GlobalEnv$newMSL$rank[splitDeleteIdx]," ",.GlobalEnv$newMSL$name[splitDeleteIdx],
                            " from ",.GlobalEnv$newMSL$lineage[splitDeleteIdx] )
    )
    
    # check for kids - work back to front for rowIdx stability across row deletions
    for(rowIdx in rev(which(splitDeleteIdx)) ) {
      srcKids=(.GlobalEnv$newMSL$parent_id==.GlobalEnv$newMSL$taxnode_id[rowIdx])
      
      if(sum(srcKids,na.rm=TRUE)>0) {
        # error: can't abolish something with kids
        errorDf=addError(errorDf,.GlobalEnv$newMSL[row,]$.split_code,.GlobalEnv$newMSL[row,]$.split_linenum, 
                         "split_abolish",as.character(.GlobalEnv$newMSL[row,]$rank),.GlobalEnv$newMSL[row,]$name,
                         "ERROR", "SPLIT.IMPLICIT_ABOLISH_WITH_KIDS", "Change=ABOLISH, taxon still has un-abolished/moved children", 
                         paste0("taxon=", .GlobalEnv$newMSL[row,]$name, ", lineage=",.GlobalEnv$newMSL[row,]$lineage,", kids: N=",sum(srcKids,na.rm=TRUE),
                                ", NAMES=[", join(.GlobalEnv$newMSL[srcKids,"rank"],.GlobalEnv$newMSL[srcKids,"name"],sep=":") 
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
    errorDf=addError(errorDf,.GlobalEnv$newMSL$.split_code[splitKeepIdx],.GlobalEnv$newMSL$.split_linenum[splitKeepIdx],
                     "split_keep",.GlobalEnv$newMSL$rank[splitKeepIdx], .GlobalEnv$newMSL$name[splitKeepIdx],
                     "INFO", "SPLIT.KEEP", 
                     "Change=SPLIT, but one split directive kept the original name", 
                     paste0("split and keep ",.GlobalEnv$newMSL$rank[splitKeepIdx]," '",.GlobalEnv$newMSL$name[splitKeepIdx],
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
# ------------------------------------------------------------------------------
#
# GLOBAL VARS
# 
# ------------------------------------------------------------------------------

#
# DEBUG ONLY - delete all load/QC errors in allErrorDf
#.GlobalEnv$allErrorDf = .GlobalEnv$loadErrorDf %>% filter(TRUE)
.GlobalEnv$allErrorDf = .GlobalEnv$loadErrorDf
# debug 
#cat("allErrorDf:",tracemem(allErrorDf),"\n");
if( params$tmi ) {
    cat("Summary of errorDf$code:" )
    summary(as.factor(errorDf$code))
    cat("Summary of allErrorDf$code:" )
    summary(as.factor(allErrorDf$code))
}
#
# iterate through proposals
#   iterate through changes
#
skipped=0
processed=0
codes = rownames(proposals)
#codes = c('2022.085B'); code=codes[1]; code; # debug
for( code in codes) {
  # code = rownames(proposals)[1]; code # debug 
  # code = rownames(proposals)[2]; code # debug 
  # code = rownames(proposals)[3]; code # debug 
  # code = rownames(proposals)[4]; code # debug 
  # code = rownames(proposals)[8]; code # debug - 20221013 current crash: around "Suid alphaherpesvirus 1"
  # make sure we were able to load the changesDf
  # code='2022.001A' # create Itzamnaviridae in [tree root]
  # code='2022.085B'; code; # debug - MOVE test
  # code = rownames(proposals)[which(code==rownames(proposals))+1];code; # DEBUG NEXT
  # code = codes[which(codes==code)+1] # DEBUG NEXT
  
  changeDf = changeList[[code]]
  if( is.null(changeDf) ) {
    skipped = skipped+1
    if(params$verbose>1) {cat("SKIP ",code,": no changeDf loaded\n")}
  } else {
    #
    # iterate over changes
    #
    cat(paste0("# START PROC: ",code," with ", nrow(changeDf), " changes ",
        " (",which(codes==code)," out of ",length(codes),")\n"))
    
    # debug if code of intereif( code %in% c("2022.001S") ) { cat("DEBUG ",code,"\n"); if(interactive()) {browser() }}
    
    
    # Apply the Change
    results=apply_changes(code,proposals[code,"cleanbase"],changeDf)
    
    # Internal sanity check
    if(!(20070000 %in% .GlobalEnv$newMSL$ictv_id)) { 
      cat("!!!!!!!!!!!!!!!! LOST ROOT NODE in ",code," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      if(interactive()) {browser()}
    }
    errorDf = results[["errorDf"]]
    levelSummary = summary(errorDf$level)
    levelSummary = levelSummary[levelSummary>0]
    cat("# PROCCESSED:",code,"with",paste(names(levelSummary),levelSummary,sep=":"),"\n")
    if(nrow(errorDf)>0){.GlobalEnv$allErrorDf=rbindlist(list(.GlobalEnv$allErrorDf,errorDf),fill=TRUE)}
    
    processed=processed+1
  } # if have proposal XLSX
} # proposal code
cat("# DONE. Found",length(rownames(proposals)),"proposals; Processed",processed,", skipped",skipped,"\n")

#### SECTION write error list #####
#
# set "final" to be TRUE (which produces per-subcommittee error files) only if not in "validate" mode

write_error_summary(.GlobalEnv$allErrorDf,params$processing_mode != "validate")

#### EXPORT MSL ####
if(params$export_msl) {
  # ------------------------------------------------------------------------------
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
    #"prev_taxnode_id",
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
  cat("WROTE   ", newTsvFilename, "\n")
  
  # ------------------------------------------------------------------------------
  # 
  # # Export proposal summary to TSV
  #
  # metadata from DOCX
  #
  # ```{r tsv_export_docx_meta}
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
    if(!(col %in% names(proposals))) {
      proposals[,col] = ""
    }
  }
  #
  # subset columns
  tsvDocxDf = data.frame(proposals[,tsvDocxColList])
  
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
  cat("WROTE   ", docxTsvFilename, "\n")
  
  # ```
  
  # ------------------------------------------------------------------------------
  # 
  # # Export SQL to update the db
  # 
  # ```{r sql_export}
  #### SECTION export to SQL #####
  #### open issues #####
  # 1. why is there a "molecule" column?
  # 2. notes vs comments columns? 
  #
  # ------------------------------------------------------------------------------
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
  # ------------------------------------------------------------------------------
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
              #"prev_taxnode_id",
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
  
  #
  # add new MSL to [taxonomy_toc]
  #
  cat("insert into [taxonomy_toc] ([tree_id],[msl_release_num],[comments]) ",
      "values (", paste(
       sort(levels(as.factor(newMSL$tree_id)),decreasing = T)[1],
        sql_msl_num,
        "NULL",
      sep=","),
      ")\n",
      file=sqlout
  )
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
  # ------------------------------------------------------------------------------
  #
  # SQL to update out_ in prevMSL
  #
  # column names and data codings in newMSL should match the taxonomy_node table. 
  #
  # convert newMSL to a data.frame, so we can convert factor columns to character, 
  # then generate SQL
  #
  # ------------------------------------------------------------------------------
  
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

