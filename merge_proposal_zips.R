
params = list(
  # inputs
  prev_msl=37
  ,next_msl=38
  ,msl_name="2022"
  ,msl_notes="EC 54, 2022 **PROVISIONAL** (MSL #38)"
  ,taxnode_delta=100000
  ,db_rank_fname="./current_msl/taxonomy_level.txt"
  ,db_molecule_fname="./current_msl/taxonomy_molecule.txt"
  ,prev_taxa_fname="./current_msl/taxonomy_node_export.txt"
  ,cv_xlsx="./TP_Template_Excel_module_2022_v2.xlsx"
  ,cv_sheet="Menu Items (Do not change)"
  ,VMR_filename="./current_msl/VMR_21-221122_MSL37.xlsx"
  ,templateURL="https://ictv.global/taxonomy/templates"
  ,proposals_dir="./proposals2"
  ,out_dir="./results2"
  # output files
  ,dest_msl=38
  ,merged="load_next_msl.txt"
  ,status="merged_status.txt"
  # debug output: 0=none, 1=some, 2=details
  ,verbose=1
)


# QC MSL38: https://uab-lefkowitz.atlassian.net/browse/IVK-123
# Merge++:  https://uab-lefkowitz.atlassian.net/browse/IVK-22

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


#
# WARNING: we use data.TABLE instead of data.FRAME
#
# this allows modification in place of a (data) passed
# to a subroutine (pass-by-reference feature)
library(data.table)

library(yaml)
library(tidyverse)
library(readxl)
library(writexl) # another option library(openxlsx)
#library(gtools) # for mixedsort/mixedorder
library(DescTools) # AscToChar

library(knitr) #kable
# debug - echo everything
knitr::opts_chunk$set(echo = TRUE)

#### load cache (xxx/.Rdata) ####

# rm("xlxsList","changeList") # to re-run xlsx file loading & QC, delete these caches
# rm("changeList") # to re-run only QC, delete this cache
cacheFilename=paste0(params$proposals_dir,"/.RData")
if(FALSE && file.exists(cacheFilename)) {
  cat("Loading image ", cacheFilename, "...\n")
  load(file=cacheFilename)
  params$out_dir=str_replace(params$proposals_dir,"proposals","results")
  rm("changeList") # to re-run only QC, delete this cache
  cat("RM(changeList) # re-run Qn")
}
cat("PROPOSALS_DIR:",params$proposals_dir,"\n")
cat("OUT_DIR      :",params$out_dir,"\n")

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
  # curMSL = curMSL; prev_msl=params$prev_msl; dest_msl=params$next_msl; taxnode_delta=params$taxnode_delta

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
  newMSL[,(fkTaxnodeIdCols) := lapply(.SD, function(id) id+params$taxnode_delta),.SDcols=fkTaxnodeIdCols]
  
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
# ```
# # load previous MSLs
# ```{r load prev MSL, echo=F}

#
# this uses the DATABASE schema for the taxonomy_node table (ie, naming convention)
#

taxonomy_node_seq = "select  * from taxonomy_node_export"
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
taxonomyDt=fread(file=params$prev_taxa_fname,
                 header=FALSE,col.names=names(taxonomy_node_names),colClasses=as.character(taxonomy_node_names),
                 stringsAsFactors=FALSE,na.strings=c("","NULL"))
cat("Previous taxa:",dim(taxonomyDt), " from ",params$prev_taxa_fname,"\n")


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

rankCV = read.delim(file=params$db_rank_fname,header=TRUE, stringsAsFactors=TRUE,na.strings=c("NULL"))
rownames(rankCV) = rankCV$id
cat("RankCV: ", dim(rankCV), " from ",params$db_rank_fname,"\n")
dbCvList[["rank"]] = rankCV
dbCvMapList[["rank"]] = rankCV$id
names(dbCvMapList[["rank"]]) = rankCV$name


moleculeCV = read.delim(file=params$db_molecule_fname,header=TRUE, stringsAsFactors=TRUE,na.strings=c("NULL"))
rownames(moleculeCV) = moleculeCV$id
cat("MoleculeCV: ", dim(moleculeCV), " from ",params$db_molecule_fname,"\n")
dbCvList[["molecule"]] = moleculeCV
dbCvMapList[["molecule"]] = moleculeCV$id
names(dbCvMapList[["molecule"]]) = moleculeCV$abbrev

# ```
# # Load CVs for Proposal QC
# ```{r load cvs}

#
# this uses the PROPOSAL.XLSX schema (naming convention)
#

proposalCV = suppressMessages(data.frame(read_excel(params$cv_xlsx,sheet = params$cv_sheet,col_names = FALSE)))
#cvDf = data.f rame(trib[,])  # remove "select one" line

cvList=list()
for(cv_col in 1:ncol(proposalCV)) {
  cv_name = proposalCV[1,cv_col]
  cv = proposalCV[,cv_col][-1]
  cvList[[cv_name]]=c(cv[!is.na(cv)],NA)
}

# map to actual input xlsx column names
cvNameMap = c(
  "Genome coverage"=    "genomeCoverage",
  "Genome composition"= "molecule",
  "Host/Source"=        "hostSource",
  "Change"=             "change",
  "Rank"=               "rank"  
)
names(cvList)=cvNameMap[names(cvList)]

# remove ("Please select",NA) from "change" & "rank" CVs - that is a required field
for( cv in c("change","rank") ) {
  isRemoveTerm = cvList[[cv]] %in% c("Please select",NA) 
  cvList[[cv]] = cvList[[cv]][!isRemoveTerm]
}

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
# load hostSource CV from VMR (most up-to-date)
#
vmrDf = data.frame(
    # suppress errors about column names
    # https://github.com/tidyverse/readxl/issues/580#issuecomment-519804058
    suppressMessages(
      read_excel(
        path    = params$VMR_filename,
        sheet   = "Terms",
        trim_ws = TRUE,
        na      = "Please select",
        skip    = 0,
        range   = cell_cols("A:AO"),
        col_names = TRUE
      )
    )
  )
cvList[["hostSource"]] = union(
    cvList[["hostSource"]],
    vmrDf[,"Host.Source"]
)
cat("VMR(hostSource): ", length(vmrDf[,"Host.Source"]), " from ",params$VMR_filename," [Terms]\n")

# ```
# 
# # scan for proposal .xlsx files
# ```{r scan for proposal matching xlsx and docx}
# 
# error reporting data frame
# 
allErrorDf = data.table(
  "folder" = character(),
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
  "notes" = character()
)
loadErrorDf = allErrorDf %>% filter(FALSE)
dir.create(params$out_dir,recursive=T,showWarnings = F)

#
#### SECTION scan DOCX #### 
#
# break filename into proposal code, filenaem and folder name
#
proposals = data.frame(docxpath=list.files(path=params$proposals_dir,pattern="20[0-9][0-9].[0-9A-Z]+.*.v*.*.docx", recursive=T, full.names=TRUE) )
proposals$path     = gsub(     "^(.*/)(20[0-9]+.[0-9A-Z]+).[^/]+.docx$","\\1",proposals$docxpath)
                    # look for folder names "Words.. (?) proposals" - there may be folders above and below that we don't care about
proposals$folder   = gsub(".*/([^/]+\\([A-Z]\\)[^/]+).*/(20[0-9]+.[0-9A-Z]+).[^/]+.docx$","\\1",proposals$docxpath)
proposals$code     = gsub(".*/([^/]+)/(20[0-9]+.[0-9A-Z]+).[^/]+.docx$","\\2",proposals$docxpath)
proposals$basename = gsub(".*/([^/]+)/(20[0-9]+.[0-9A-Z]+.[^/]+).docx$","\\2",proposals$docxpath)
proposals$docx     = paste0(proposals$basename,".docx")
# remove all *.Ud.* files
proposals = proposals[grep(proposals$docx, pattern=".*\\.Ud\\..*", invert=T),]

# check for duplicate Proposal IDs
dups = duplicated(proposals$code)
allDups =proposals$code %in% proposals$code[dups]
if(sum(dups) > 0) {
  errorDf = proposals[allDups, c("folder", "code", "docx")]
  errorDf$level = "ERROR"
  errorDf$error = "DUPCODE.DOCX"
  errorDf$message = "duplicate proposal ID"
  loadErrorDf = rbindlist(list(loadErrorDf, errorDf),fill=TRUE)
  # output
  kable(errorDf,caption = paste0("QC01: ERRORS: dupliate docx proposal IDs"))
  write_xlsx(x=errorDf,path=file.path(params$out_dir,"QC01.docx_duplicate_ids.xlsx"))
  write_xlsx(x=loadErrorDf,path=file.path(params$out_dir,"QC.summary.xlsx"))

}
rownames(proposals)=proposals$code

#
#### SECTION: scan XLSX #### 
#
xlsxs = data.frame(xlsxpath=list.files(path=params$proposals_dir,pattern="^20[0-9][0-9].[^.]+.*.v*.*.xlsx", recursive=T, full.names=TRUE) )
xlsxs$path     = gsub(     "^(.*/)(20[0-9]+.[0-9A-Z]+).[^/]+.xlsx$","\\1",xlsxs$xlsxpath)
# look for folder names "Words.. (?) proposals" - there may be folders above and below that we don't care about
xlsxs$folder   = gsub(".*/([^/]+\\([A-Z]\\)[^/]+).*/(20[0-9]+.[0-9A-Z]+).[^/]+.xlsx$","\\1",xlsxs$xlsxpath)
xlsxs$basename = gsub(".*/([^/]+)/(20[0-9]+.[0-9A-Z]+.[^/]+).xlsx$","\\2",xlsxs$xlsxpath)
xlsxs$code     = gsub(".*/([^/]+)/(20[0-9]+.[0-9A-Z]+).[^/]+.xlsx$","\\2",xlsxs$xlsxpath)
xlsxs$xlsx     = paste0(xlsxs$basename,".xlsx")
# remove all *.Ud.* files
xlsxs = xlsxs[grep(xlsxs$xlsx, pattern=".*\\.Ud\\..*", invert=T),]

# ignore "Suppl" files 
sups = grep(xlsxs$xlsx,pattern="_Suppl." )
xlsxs=xlsxs[-(sups),]

# QC for duplicate codes
dups = duplicated(xlsxs$code)
allDups =xlsxs$code %in% xlsxs$code[dups]
if( sum(dups) > 0 ) {
  # error details
  errorDf = xlsxs[allDups, c("folder", "code", "xlsx")]
  errorDf$level = "ERROR"
  errorDf$error = "DUPCODE.XLSX"
  for( code in xlsxs$code[dups] ) {
    # build list of all xlsxs using duplicate codes
    errorDf[errorDf$code==code,"message"] = 
      paste0("duplicate proposal ID: ",
             paste(xlsxs$xlsx[xlsxs$code==code],collapse=","))
    
  }
  # append to global list
  loadErrorDf = rbindlist(list(loadErrorDf, errorDf),fill=TRUE)
  # output
  kable(errorDf,caption = paste0("QC01: ERRORS: dupliate docx proposal IDs"))
  write_xlsx(x=errorDf,path=file.path(params$out_dir,"QC01.docx_duplicate_ids.xlsx"))
  write_xlsx(x=loadErrorDf,path=file.path(params$out_dir,"QC.summary.xlsx"))

  kable(caption=paste0("ERROR: XLSX dupliate proposal IDs"),
        x=proposals[proposals$code %in% proposals$code[dups],])
}
rownames(xlsxs) = xlsxs$basename

#
# merge XLSX list into DOCX list, verify
#
proposals$xlsx = xlsxs[proposals$basename,"xlsx"]
proposals$xlsxpath = xlsxs[proposals$basename,"xlsxpath"]
missing= is.na(proposals$xlsx)
if( sum(missing) > 0 ) {
  errorDf= proposals[missing,c("folder","code","docx")]
  # suggest possible matches based on ID
  errorDf$xlsx= NA 
  errorDf$row = NA
  errorDf$level= "ERROR"
  errorDf$error = "XLSX.MISSING"
  errorDf$message = "DOCX has no matching XLSX"
  errorDf$notes= "Suggestions: contact corresponding author"
  for(row in rownames(errorDf) ) {
    guesses = sum(xlsxs$code==errorDf[row,"code"])
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
      errorDf[row,"notes"] = paste("SUGGESTIONS:",paste(paste0(xlsxs[xlsxs$xlsxID==errorDf[row,"docxID"],"basename"],".xslx"),collapse=", "))
    } else {
      # 
      # can't find the xlsx
      #
      # let the error pass through
      proposals[]
    }
  }
  # append to global list
  loadErrorDf = rbindlist(list(loadErrorDf, errorDf),fill=TRUE)
  # QC output
  kable(errorDf,caption = paste0("QC02: ERRORS: DOCX without matching XLSX"))
  write_xlsx(x=errorDf,path=file.path(params$out_dir,"QC02.docx_without_matching_xlsx.xlsx"))
  write_xlsx(x=errorDf,path=file.path(params$out_dir,"QC.summary.xlsx"))

  }
# ```
# 
# # Load summary
# 
# ```{r load_summary,echo=FALSE}
# kable(caption=paste0("SUMMARY: Proposal .xlsx file found in ",params$proposals_dir,"/"),
#       x=data.frame(nProposals=summary(as.factor(proposals$folder))),)
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

# dput(unname(as.vector(df[2,])))
xlsx_row3_v1=c(
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
xlsx_row3_v2=c(
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
xlsx_col_info = data.frame(row.names=xlsx_change_colnames,
                           pattern = rep(NA_character_,length(xlsx_change_colnames)),
                           pat_warn = rep(NA_character_,length(xlsx_change_colnames))
                           )
# most taxa can only be alphanumeric plus a hyphen
xlsx_col_info[c(xlsx_change_colnames[c(xlsx_change_srcCols,xlsx_change_destCols)],"rank"),"pattern"] = "([^[:alnum:]-]+)"
xlsx_col_info[c(xlsx_change_colnames[c(xlsx_change_srcCols,xlsx_change_destCols)],"rank"),"pat_warn"] = "non-(AlphaNumeric or hyphen)"
xlsx_col_info[c(xlsx_change_colnames[c(xlsx_change_srcCols,xlsx_change_destCols)],"rank"),"pat_code"] = "XLSX.NON_ALPHA-NUMERIC"
xlsx_col_info[c(xlsx_change_colnames[c(xlsx_change_srcCols,xlsx_change_destCols)],"rank"),"pat_replace"] = ""
# allow internal spaces on species and columns with lists/isolate names
xlsx_col_info[c("srcSpecies","species"),"pattern"] = "([^[:alnum:] -]+)"
xlsx_col_info[c("srcSpecies","species"),"pat_warn"] = "non-(AlphaNumeric,hyphen,space)"
xlsx_col_info[c("srcSpecies","species"),"pat_code"] = "XLSX.NON_ALPHA-NUMERIC-SPACE"
# allow internal spaces on species and columns with lists/isolate names
xlsx_col_info[c("change"),"pattern"] = "([^[:alnum:] ;-]+)"
xlsx_col_info[c("change"),"pat_warn"] = "non-(AlphaNumeric,hyphen,space,semicolon)"
xlsx_col_info[c("change"),"pat_code"] = "XLSX.NON_ALPHA-NUMERIC-SPACE-SEMI"
# semi-colon separated lists
xlsx_col_info[c( "Abbrev"),"pattern"] = "([^[:alnum:];_/ . -]+)"
xlsx_col_info[c( "Abbrev"),"pat_warn"] = "Should be semicolon-separated list"
xlsx_col_info[c( "Abbrev"),"pat_code"] = "XLSX.SEMICOLON-SEP-LIST"
xlsx_col_info[c( "exemplarAccession"),"pattern"] = "([^[:alnum:]:;_/ .-]+)"
xlsx_col_info[c( "exemplarAccession"),"pat_warn"] = "Should be semicolon-separated list of [name:]accession"
xlsx_col_info[c( "exemplarAccession"),"pat_code"] = "XLSX.SEMICOLON-SEP-LIST-NAMED"
# exemplarIsolate and exemplarName are anything printable. 
xlsx_col_info[c("exemplarName","exemplarIsolate"),"pattern"] = "([^[:print:]]+)"
xlsx_col_info[c("exemplarName","exemplarIsolate"),"pat_warn"] = "unprintable/non-ASCII"
xlsx_col_info[c("exemplarName","exemplarIsolate"),"pat_code"] = "XLSX.NON_PRINTABLE"

# others just must be printable
otherCols = is.na(xlsx_col_info$pattern)
xlsx_col_info[otherCols,"pattern"] = "([^[:print:]]+)"
xlsx_col_info[otherCols,"pat_warn"] = "unprintable/non-ASCII"
xlsx_col_info[otherCols,"pat_code"] = "XLSX.NON_PRINTABLE"

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
    return(NA)
  }
}
getParentTaxon = function(realmSpecies) {
  #
  # get the NAME of the penultimate low ranked taxon
  #
  nonNAs = !is.na(realmSpecies[realmSpeciesCols])
  if( sum(nonNAs) > 1 ) {
    return(realmSpecies[sort(which(nonNAs),decreasing=T)[2]])
  } else {
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
    return(NA)
  }
}
getParentLineage = function(realmSpecies) {
  nonNAs = !is.na(realmSpecies[realmSpeciesCols])
  if( sum(nonNAs) > 0 ) {
    rank=getTaxonRank(realmSpecies)
    return(getLineage(realmSpecies,masks=c(rank)))
  } else {
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
    return(NA)
  }
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
  # read XLSX file, no column names
  df = data.frame(
    # suppress errors about column names
    # https://github.com/tidyverse/readxl/issues/580#issuecomment-519804058
    suppressMessages(
      read_excel(
        proposals[code,"xlsxpath"],
        trim_ws = TRUE,
        na = "Please select",
        skip = 2,
        range = cell_cols("A:AO"),
        col_names = FALSE
      )
    )
  )
  # return value
  return(df)
}


#
# QC the .xlsx data at the sheet level
#
# returns: errorDf
addError=function(errorDf,code,row,change,rank,taxon,levelStr,errorCode,errorStr,notes) {
  nextErrorDf = data.frame(
    folder = proposals[code,]$folder,
    code = code,
    row = row,
    change = change, 
    rank=rank,
    taxon = taxon,
    docx = proposals[code,]$docx,
    xlsx = proposals[code,]$xlsx,
    level = levelStr,
    error = errorCode,
    message = errorStr, 
    notes = notes)
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
                                paste(paste0("column ",toupper(c(letters,paste0("a",letters)))[which(!row2v1match)],"='",proposalDf[2,which(!row2v1match)],"' instead of '",xlsx_row2_v1[which(!row2v1match)],"'"),collapse="; ")
      )
    } 
    if( row2v1matchCt < row2v2matchCt ) { 
      templateLine2Error= paste("similar to v2, but with ", 
                                paste(paste0("column ",toupper(c(letters,paste0("a",letters)))[which(!row2v2match)],"='",proposalDf[2,which(!row2v2match)],"' instead of '",xlsx_row2_v2[which(!row2v2match)],"'"),collapse="; ")
      )
    } 
  }
  if(params$verbose>1) { cat("     ", code, " XLSX.Row 2: ",templateLine2Version," ",templateLine2Error,"\n")}
  
  # check row 3 to validate version of templates
  row3v1match = proposalDf[3, seq(from = 1, to = min(length(xlsx_row3_v1), ncol(proposalDf)))] == xlsx_row3_v1
  row3v1matchCt = sum(row3v1match, na.rm=T)
  row3v2match = proposalDf[3, seq(from = 1, to = min(length(xlsx_row3_v2), ncol(proposalDf)))] == xlsx_row3_v2
  row3v2matchCt = sum(row3v2match, na.rm=T)
  templateLine3Error = ""
  if (sum(row3v1match,na.rm=T) == sum(!is.na(xlsx_row3_v1))) {
    templateLine3Version = "v1"
  } else if (sum(row3v2match, na.rm=T) == sum(!is.na(xlsx_row3_v2))) {
    templateLine3Version = "v2"
  } else {
    templateLine3Version = "unrecognized"
    # report the unexpected values
    if( row3v1matchCt > row3v2matchCt ) { 
      templateLine3Error= paste("similar to v1, but with ", 
                                paste(paste0("column ",toupper(c(letters,paste0("a",letters)))[which(!row3v1match)],"='",proposalDf[3,which(!row3v1match)],"' instead of '",xlsx_row3_v1[which(!row3v1match)],"'"),collapse="; ")
      )
    } 
    if( row3v1matchCt < row3v2matchCt ) { 
      templateLine3Error= paste("similar to v2, but with ", 
                                paste(paste0("column ",toupper(c(letters,paste0("a",letters)))[which(!row3v2match)],"='",proposalDf[3,which(!row3v2match)],"' instead of '",xlsx_row3_v2[which(!row3v2match)],"'"),collapse="; ")
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
  if(params$verbose){cat("     ", code, " XLSX template ",templateVersion,"\n")}
  
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
      errorDf=addError(errorDf,code,codeRow,NA,NA,NA,"INFO","XLSX.CODE_MISS", "XLSX code missing", 
                       paste0("XLSX cell ",codeCell,
                              " is ", "'",codeValue, 
                              "'; replace with the actual code: '",code,"'")
      )
    } else {
      errorDf=addError(errorDf,code,codeRow,NA,NA,NA,"WARNING", "XLSX.CODE_BAD","XLSX code wrong", 
                       paste0("XLSX cell ",codeCell,
                              " does not match proposal code from filename: ", 
                              "'",codeValue, "' should be '", code,"' ")
      )
    }
  }
  
  #
  # extract & standardize
  # cols&rows with change data 
  #
  changeDf = data.frame()

  # map columns
  firstDataRow=4
  if(templateVersion=="v1") { changeDf = proposalDf[firstDataRow:nrow(proposalDf),xlsx_v1_change_cols] }
  if(templateVersion=="v2") { changeDf = proposalDf[firstDataRow:nrow(proposalDf),xlsx_v2_change_cols] }
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
  # extract src/dest taxon names for error reporting
  #
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
  # QC spaces, quotes, etc
  # 
  # TODO: remove leading/trailing spaces, quotes
  for( col in rownames(xlsx_col_info) ) {
    # col = "srcSpecies" # debug
    #  col = "hostSource" # debug

    #
    # non-breaking space
    # AscToChar(202)
    pattern="([ ]+)"; pat_warn="non-breaking space character"; pat_replace=" "
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(length(qc.matches)>0) { 
      if(params$verbose) { cat("WARNING:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") }
       errorDf=addError(errorDf,code,rownames(changeDf)[qc.matches],
                       changeDf$change[qc.matches],changeDf$rank[qc.matches],changeDf$.changeTaxon[qc.matches],
                       "WARNING","XLSX.NB_SPACE", paste("XLSX has",pat_warn),
                       paste0(paste(col,gsub(pattern,"[\\1]",changeDf[qc.matches,col]),sep=":")," (replacing with '",pat_replace,"')")
      )
      # backup original values
      changeDf[qc.matches,paste0(col,"_orig")] == changeDf[qc.matches,col]
      # remove non-ascii chars
      changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    }     

    #
    # newline
    #
    pattern="([\r\n]+)"; pat_warn="newline character(s)";pat_replace=";"
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(length(qc.matches)>0) { 
      if(params$verbose) { cat("WARNING:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") }
      errorDf=addError(errorDf,code,rownames(changeDf)[qc.matches],
                       changeDf$change[qc.matches],changeDf$rank[qc.matches],changeDf$.changeTaxon[qc.matches],
                       "WARNING","XLSX.NEWLINE", paste("XLSX has",pat_warn),
                       paste0(paste(col,gsub(pattern,"[\\1]",changeDf[qc.matches,col]),sep=":")," (replacing with '",pat_replace,"')")
      )
      # backup original values
      changeDf[qc.matches,paste0(col,"_orig")] == changeDf[qc.matches,col]
      # remove non-ascii chars
      changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    }     
    
    #
    # long-dashes
    # "–" and "—"
    pattern=paste0("([",AscToChar(207),AscToChar(208),"]+)"); pat_warn="long-dash[en/em dashes]";pat_replace="-"
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(length(qc.matches)>0) { 
      if(params$verbose) { cat("WARNING:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") }
      errorDf=addError(errorDf,code,rownames(changeDf)[qc.matches],
                       changeDf$change[qc.matches],changeDf$rank[qc.matches],changeDf$.changeTaxon[qc.matches],
                       "WARNING","XLSX.NEWLINE", paste("XLSX has",pat_warn),
                       paste0(paste(col,gsub(pattern,"[\\1]",changeDf[qc.matches,col]),sep=":")," (replacing with '",pat_replace,"')")
      )
      # backup original values
      changeDf[qc.matches,paste0(col,"_orig")] == changeDf[qc.matches,col]
      # remove non-ascii chars
      changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    }     
 
    #
    # curvy quotes
    # AscToChar(210)AscToChar(211)
    pattern=paste0("([“”]+)"); pat_warn="curvy quotes";pat_replace='"'
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(length(qc.matches)>0) { 
      if(params$verbose) { cat("WARNING:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") }
      errorDf=addError(errorDf,code,rownames(changeDf)[qc.matches],
                       changeDf$change[qc.matches],changeDf$rank[qc.matches],changeDf$.changeTaxon[qc.matches],
                       "WARNING","XLSX.CURVY_DOUBLE_QUOTES", paste("XLSX has",pat_warn),
                       paste0(paste(col,gsub(pattern,"[\\1]",changeDf[qc.matches,col]),sep=":")," (replacing with '",pat_replace,"')")
      )
      # backup original values
      changeDf[qc.matches,paste0(col,"_orig")] == changeDf[qc.matches,col]
      # remove non-ascii chars
      changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    }     
    #
    # curvy single-quotes
    # AscToChar(212)AscToChar(213)
    pattern=paste0("([‘’]+)"); pat_warn="curvy single quotes";pat_replace="'"
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(length(qc.matches)>0) { 
      if(params$verbose) { cat("WARNING:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") }
      errorDf=addError(errorDf,code,rownames(changeDf)[qc.matches],
                       changeDf$change[qc.matches],changeDf$rank[qc.matches],changeDf$.changeTaxon[qc.matches],
                       "WARNING","XLSX.CURVY_SINGLE_QUOTES", paste("XLSX has",pat_warn),
                       paste0(paste(col,gsub(pattern,"[\\1]",changeDf[qc.matches,col]),sep=":")," (replacing with '",pat_replace,"')")
      )
      # backup original values
      changeDf[qc.matches,paste0(col,"_orig")] == changeDf[qc.matches,col]
      # remove non-ascii chars
      changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    }     
    
    #
    # leading white space
    #
    pattern="^([ \t]+)"; pat_warn="leading whitespace"; pat_replace=""
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(length(qc.matches)>0) { 
      if(params$verbose) { cat("WARNING:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") }
      errorDf=addError(errorDf,code,rownames(changeDf)[qc.matches],
                       changeDf$change[qc.matches],changeDf$rank[qc.matches],changeDf$.changeTaxon[qc.matches],
                       "WARNING","XLSX.LEAD_SPACE", paste("XLSX has",pat_warn),
                       paste0(paste(col,gsub(pattern,"[\\1]",changeDf[qc.matches,col]),sep=":")," (replacing with '",pat_replace,"')")
      )
      # backup original values
      changeDf[qc.matches,paste0(col,"_orig")] == changeDf[qc.matches,col]
      # remove non-ascii chars
      changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    } 
    #
    # trailing white space
    #
    pattern="([ \t]+)$"; pat_warn="trailing whitespace"; pat_replace=""
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(length(qc.matches)>0) { 
      if(params$verbose) { cat("WARNING:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") }
      errorDf=addError(errorDf,code,rownames(changeDf)[qc.matches],
                       changeDf$change[qc.matches],changeDf$rank[qc.matches], changeDf$.changeTaxon[qc.matches],
                       "WARNING","XLSX.TRAIL_SPACE", paste("XLSX has",pat_warn),
                       paste0(paste(col,gsub(pattern,"[\\1]",changeDf[qc.matches,col]),sep=":")," (replacing with '",pat_replace,"')")
      )
      # backup original values
      changeDf[qc.matches,paste0(col,"_orig")] == changeDf[qc.matches,col]
      # remove non-ascii chars
      changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
    } 
    #
    # quotes
    #
    if( col != "comments" ) {
      pattern='(["]+)'; pat_warn="quote"; pat_replace=""
      qc.matches =grep(changeDf[,col],pattern=pattern)
      if(length(qc.matches)>0) { 
        if(params$verbose) { cat("WARNING:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") }
        errorDf=addError(errorDf,code,rownames(changeDf)[qc.matches],
                         changeDf$change[qc.matches],changeDf$rank[qc.matches],changeDf$.changeTaxon[qc.matches],
                         "WARNING","XLSX.QUOTES", paste("XLSX has",pat_warn),
                         paste0(paste(col,gsub(pattern,"[\\1]",changeDf[qc.matches,col]),sep=":")," (replacing with '",pat_replace,"')")
        )
        # backup original values
       changeDf[qc.matches,paste0(col,"_orig")] == changeDf[qc.matches,col]
        # remove non-ascii chars
        changeDf[,col] = gsub(pattern,pat_replace,changeDf[,col])
      } 
    }
    # 
    # more rigorous for some columns
    #
    pattern=xlsx_col_info[col,"pattern"]; pat_warn=xlsx_col_info[col,"pat_warn"]
    qc.matches =grep(changeDf[,col],pattern=pattern)
    if(length(qc.matches)>0) { 
      if(params$verbose) { cat("WARNING:",code,"has",length(qc.matches),"cells with",pat_warn,"in column",col,"\n") }
      errorDf=addError(errorDf,code,rownames(changeDf)[qc.matches],
                       changeDf$change[qc.matches],changeDf$rank[qc.matches], changeDf$.changeTaxon[qc.matches],
                       "WARNING",xlsx_col_info[col,"pat_code"], paste("XLSX has",pat_warn),
                       paste(col,gsub(pattern,"[\\1]",changeDf[qc.matches,col]),sep=":"))
      # backup original values
      changeDf[qc.matches,paste0(col,"_orig")] == changeDf[qc.matches,col]
      # remove non-ascii chars
      changeDf[,col] = gsub(pattern,"",changeDf[,col])
    }
  }  # for col
  
  #
  # QC controlled vocabularies
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
  # check that all the terms are legal
  #
  # remove lines containing invalid terms
  #
  for(cv in names(cvList)) {
    # cv="molecule" # debug
    
    # compare to CV terms, removing all spaces and capitalization
    isTermPerfect = changeDf[,cv] %in% cvList[[cv]]
    isTermClose = 
      tolower(gsub(pattern=" ",replacement="",changeDf[,cv])) %in% 
      tolower(gsub(pattern=" ",replacement="",cvList[[cv]]))

    # 
    # reject unknown terms
    #
    badTerms=changeDf[!isTermClose,]
    if(nrow(badTerms)>0) {
      # complain
      if(params$verbose) { cat("ERROR:",code,"illegal term in CV column '",cv,"':",
                              "[",paste(paste0(rownames(badTerms),":",badTerms[,cv]), collapse=","),"] on rows (",paste(rownames(badTerms),collapse=","),")\n") }
       errorDf=addError(errorDf,code,rownames(badTerms),
                       badTerms$change,badTerms$rank, badTerms$.changeTaxon,
                       "ERROR","XLSX.INVALID_TERM",paste("XLSX incorrect term in column",cv),
                     paste("XLSX incorrect value [",badTerms[,cv],"]. Valid terms: [",paste0(cvList[[cv]],collapse=","),"]")) 
      # for optional CVs, NA the term, for obligatory ones, flag the line
      if(cv %in% c("genomeCoverage","molecule","hostSource") ) {
          changeDf[rownames(badTerms),cv] = NA
      } else {
         changeDf[rownames(badTerms),".noErrors"] = FALSE
      }
      changeDf[rownames(badTerms),".errors"] = paste("CV '",cv,"' incorrect value [",badTerms[,cv],"]. Valid terms: [",paste0(cvList[[cv]],collapse=","),"]")
      #browser()
    } # if badTerms
    #
    # whine about terms with spacing/caps problems
    #
    flawedTerms=changeDf[isTermClose & !isTermPerfect,] #cv, drop=FALSE]
    if(nrow(flawedTerms)>0) {
      # complain
      if(params$verbose) { cat("WARNING:",code,"typo in CV column '",cv,"':",
                               "[",paste(paste0(rownames(flawedTerms),":",flawedTerms[,cv]), collapse=","),"] on rows (",paste(rownames(flawedTerms),collapse=","),")\n") }
      errorDf=addError(errorDf,code,rownames(flawedTerms),
                       flawedTerms$change, flawedTerms$rank, flawedTerms$.changeTaxon,
                       "WARNING","XLSX.TYPO_TERM",paste("XLSX term with typo (space,caps) in column",cv),
                       paste("XLSX incorrect value [",flawedTerms[,cv],"]. Valid terms: [",paste0(cvList[[cv]],collapse=","),"]")) 
    } # if flawedTerms
  } # for cvList
  
  # return warnings, if any
  return(list(errorDf=errorDf,changeDf=changeDf))
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
      # try to load proposal
      #
      if(is.null(xlsxList[[code]])){
        # load raw xlsx file
        xlsxList[[code]] = load_proposal(code)
        cat("# LOADED: ",code,"\n")
      } else {
        cat("# FROM CACHE: ",code,"\n")
      }
      proposalDf = xlsxList[[code]]
      #
      # qc
      #
      results = qc_proposal(code,proposalDf)
      errorDf = results[["errorDf"]]
      cat("# QCed: ",code," with ",nrow(errorDf)," errors/warnings\n")
      if(nrow(errorDf)>0){.GlobalEnv$loadErrorDf=rbindlist(list(.GlobalEnv$loadErrorDf,errorDf),fill=TRUE)}
      changeDf = results[["changeDf"]]
      if(!is.null(changeDf)) {
        changeList[[code]]=changeDf
      }
    } # LOAD if exists(xlsx)
  } # LOAD if not in changeList
}
# load summary
cat("changeList: ",paste(names(changeList)),"\n")

# error summary
write_xlsx(x=loadErrorDf,path=file.path(params$out_dir,"QC.summary.xlsx"))
kable(allErrorDf,caption = paste0("QC: Summary of ERRORs and WARNINGs"))
 
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
    .GlobalEnv$newMSL[kid_rows,"lineage"] = paste(parent_lineage,.GlobalEnv$newMSL[kid_rows,"name"],sep=";")

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
#  "Exemplar \r\nvirus name"           = "exemplarName",                        "exemplarName"              NA                          "exemplar_name"
# "Virus name abbreviation"            = "Abbrev",                              "Abbrev"                    "abbrev_csv"                "abbrev_csv"
# "Exemplar\r\nisolate designation"    = "exemplarIsolate",                     "exemplarIsolate"           "isolate_csv"               "isolate_csv"
# "Genome coverage"                    = "genomeCoverage",                      "isComplete"                NA                          "genome_coverage"
# "Genome composition"                 = "molecule",                            "molecule"                  "molecule_id"               "molecule_id"
# "Host/Source"                        = "hostSource",                          "hostSource"                NA                          "host_source"
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
      ,"comments" = "comments"
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
  
  # action mappings - externalize or move up
  actionCV = c(
    "create" = "new"
    ,"create new" = "new"
    ,"rename" = "rename"
    ,"abolish" = "abolish"
    ,"move" = "move"
    ,"move; rename" = "move"
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
      errorDf = with(changeDf[linenum,],addError(errorDf,code,linenum,change,rank,.changeTaxon,"
                                                 ERROR","ACTION.UNK","Unknown action/change",action))
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

    # debug
    if(change$.destTaxon %in% c('Innmovirus hailarense','Innmovirus') ) {
      print(paste0(code,":",linenum," ", actionClean, "(",srcLineage,", ",destLineage,")"))
      print(paste("srcNewMatches=",sum(srcNewMatches),", destCurMatches=",sum(destCurMatches),",destOldMatches=",sum(destOldMatches),",destNewMatches=",sum(destNewMatches)))
      #browser()
    }
    
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
    if(!is.na(destTaxonName) && !(actionClean %in% c("move")) )  {
      # check if exists in the current MSL
      if(sum(destCurMatches)>0) {
        matchDf = .GlobalEnv$curMSL[destCurMatches,c("msl_release_num","lineage")]
        matchDf = .GlobalEnv$curMSL[order(matchDf$msl_release_num,decreasing = T)[1],]
        matchLineage = paste0( "MSL", matchDf$msl_release_num,":",matchDf$lineage)
        errorDf=addError(errorDf,code,linenum,change$change,change$rank,change$.changeTaxon,
                         "ERROR", "DEST.IN_CUR", 
                         paste0("Change=",toupper(action),", but 'proposed taxonomy' already exists"), 
                         paste0(", proposedTaxonomy=", destTaxonName, 
                                "; existingTaxon=", matchLineage)
        )
        next;
      }
      # check if exists in a historical MLS (not the most recent)
      if(sum(destOldMatches)>0) {
        matchDf = .GlobalEnv$oldMSLs[destOldMatches,c("msl_release_num","lineage")]
        matchDf = .GlobalEnv$oldMSLs[order(matchDf$msl_release_num,decreasing = T)[1],]
        matchLineage = paste0( "MSL", matchDf$msl_release_num,":",matchDf$lineage)
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "ERROR", "DEST.IN_OLD", 
                         paste0("Change=",toupper(action),",  but 'proposed taxonomy' existed historically"), 
                         paste0("proposedTaxonomy=", destTaxonName, 
                                "; existingTaxon=", matchLineage)
        )
        next;
      }
      # check if already created in newMSL
      if(sum(destNewMatches)>0) {
        matchDf = .GlobalEnv$newMSL[destNewMatches,c("msl_release_num","lineage")]
        matchDf = .GlobalEnv$newMSL[order(matchDf$msl_release_num,decreasing = T)[1],]
        matchLineage = paste0( "MSL", matchDf$msl_release_num,":",matchDf$lineage)
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "ERROR", "DEST.IN_NEW", 
                         paste0("Change=",toupper(action),", but 'proposed taxonomy' already created in new MSL"), 
                         paste0("proposedTaxonomy=", destTaxonName, 
                                "; existingTaxon=", matchLineage,
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
    #### CREATE ####
    #
    # TODO: also handle SPLIT?
    #
    # Note: doesn't fix left_idx/right_idx, etc
    #  -------------------------------------------------------------------------
    if(actionClean %in% c("new") ) {

      # check if srcTaxon was specified in xlsx (shouldn't be)
      if(!is.na(srcTaxonName)) {
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                         "WARNING", "CREATE.W_SRC", "Change=CREATE, but 'current taxonomy' columns are not empty", 
                         paste0("currentTaxonomy=", change$.srcLineage,
                                ", prposedTaxonomy=",change$.destLineage)
        )
       }
     
      #
      # check if same accession number already exists
      #
      isDupAccession = (.GlobalEnv$newMSL$genbank_accession_csv == change$exemplarAccession)
      if(sum(isDupAccession, na.rm=TRUE)>0) {
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                         "ERROR", "CREATE.DUP_ACC", "Change=CREATE, a species with this accession number already exists", 
                         paste0("accession=", change$exemplarAccession, ", existingSpecies=",.GlobalEnv$newMSL[isDupAccession,]$lineage)
        )
        ## QQQ what proposal created this? (from this round? historically? Need a function: last_modified())
        next;
        
      }
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
         errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                          "ERROR", "CREATE.PARENT_NO_EXIST", "Change=CREATE, but parent rank taxon does not exist", 
                          paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage)
         )
         next;
      } else if(sum(parentDestNewMatches)>1) {
         # this should never happen
         errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                          "ERROR", "CREATE.PARENT_MANY", "Change=CREATE, multiple taxa exist with parent name", 
                          paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage)
         )
         next;
      } else {
        # WARN: check if taxon rank matches change rank
        if ( tolower(destTaxonRank) != tolower(change$rank) ) {
          errorDf=addError(errorDf,code,linenum,change$change,change$rank,change$.destTaxon,
                           "WARNING", "CREATE.RANK_MISMATCH", 
                           "Change=CREATE, but proposed taxon rank does not match [rank] column", 
                           paste0("ankColumn=", change$rank,
                                  ", proposedTaxonRank=",destTaxonRank,
                                  ", proposedTaxonomy=", destLineage)
          )
        }

        # 
        # create new taxon
        #
        
        
        # get parent
        newTaxon = .GlobalEnv$newMSL[parentDestNewMatches,]
        
        # WARN if PARENT_LINEAGE is not expected, AND USE PARENT LINEAGE
        if( !is.na(destParentLineage) && (newTaxon$lineage != destParentLineage) ) {
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.destTaxon,
                           "WARNING", "CREATE.PARENT_LINEAGE", 
                           "Change=CREATE, proposed parent taxon exists, but not with expected lineage, using observed lineage", 
                           paste0("proposedParentLineage=", destParentLineage,
                                  ", observedParentLineage=",newTaxon$lineage,
                                  ", otherProposals=",newTaxon$prev_proposals)
          )
          destLineage=destParentLineage
        }
        
        # add new info - primary columns
        newTaxon[1,"in_change"]   = "new"
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
        newTaxon[1,xlsx2dbMap["molecule"]] = dbCvMapList[["molecule"]][change[1,"molecule"] ]
      
        # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
        newTaxon[1,xlsx2dbMap["hostSource"]] = change[1,"hostSource"] 

        # comments
        newTaxon[1,xlsx2dbMap["comments"]]= change[1,"comments"]

        ## for species only
        if( destTaxonRank == "species" ) {
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
        #browser()
        .GlobalEnv$newMSL <- rbindlist(list(.GlobalEnv$newMSL,newTaxon),fill=TRUE)
      } # create new taxon
    } else  if(actionClean %in% c("rename") ) {
    
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
                         paste0("currentTaxonomy=", srcLineage, ", destTaxonomy=", destLineage)
        )
      }
      
      
      #
      # find the target to rename in NEW
      #
      srcNewTarget=(.GlobalEnv$newMSL$name==as.character(srcTaxonName))
      if(params$verbose) {print(paste0("RENAME: ",code," line ",linenum," '",destTaxonName, "' findTarget(",srcTaxonName,")=",sum(srcNewTarget)))}
    
      if(sum(srcNewTarget)==0) {
         errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                          "ERROR", "RENAME.NO_EXIST", "Change=RENAME, but taxon does not exist", 
                          paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName))
         next;
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
        
        # change the name and update lineage
        .GlobalEnv$newMSL[srcNewTarget,"name"] = destTaxonName
        .GlobalEnv$newMSL[srcNewTarget,"lineage"] = paste(newMSL[srcNewParent,"lineage"],destTaxonName,sep=";") # should we do this? 
        
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
        ## ZZZ add comments to out_notes?
      } # target found
      # RENAME
    } else if(actionClean %in% c("abolish") ) { 
      #  -------------------------------------------------------------------------
      #
      # ABOLISH 
      #
      #  -------------------------------------------------------------------------
      # check if srcTaxon was specified in xlsx (required)
      if(is.na(srcTaxonName)) {
        errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                         "ERROR", "ABOLISH.WO_SRC", "Change=ABOLISH, but 'current taxonomy' columns are empty", 
                         paste0("currentTaxonomy=", srcLineage, ", destTaxonomy=", destLineage)
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
         errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                          "ERROR", "ABOLISH.NO_EXIST", "Change=ABOLISH, but taxon does not exist", 
                          paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName))
         next;
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
       
    } else if(actionClean %in% c("move") ) { 
      #  -------------------------------------------------------------------------
      #
      # MOVE (w/ or w/o RENAME) 
      #
      #  -------------------------------------------------------------------------
      # check if srcTaxon was specified in xlsx (required)
      if(is.na(srcTaxonName)) {
        # they didn't specify the source, but since it's a MOVE
        # we should be able to find it
        guessSource = (.GlobalEnv$newMSL$name==as.character(destTaxonName))
        if( sum(guessSource,na.rm=TRUE) == 1) {
          # found it, so WARNING
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "WARNING", "MOVE.NO_SRC", "Change=MOVE, but 'current taxonomy' columns are empty", 
                           paste0("destTaxonomy=", destLineage,
                                  "; we guess you meant=", .GlobalEnv$newMSL$lineage[guessSource])
          )
          srcTaxonName=.GlobalEnv$newMSL$name[guessSource]
          srcTaxonRank=.GlobalEnv$newMSL$rank[guessSource]
          srcLineage  =.GlobalEnv$newMSL$lineage[guessSource]
        } else {
            # didn't find it, or found multiple, so ERROR
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "WARNING", "MOVE.NO_SRC_GUESS", "Change=MOVE, but 'current taxonomy' columns are empty", 
                           paste0("destTaxonomy=", destLineage,
                                  "; no existing taxon '",destTaxonName,"' found")
          )
          next
        }
      }
      
      #
      # validate that this move does NOT change rank (not a promote/demote)
      #
      if( srcTaxonRank != rankClean || destTaxonRank != rankClean ) {
        # src/dest ranks mis-match (promote/demote)
        if( srcTaxonRank != destTaxonRank ) {
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "ERROR", "MOVE.DIFF_RANK", 
                           "Change=MOVE, but 'current taxonomy' and 'proposed taxonomy' ranks differ; use promote or demote, not move, to change rank", 
                           paste0("srcTaxonomy=[", srcTaxonRank,"]",srcLineage,
                                  ", destTaxonomy=[", destTaxonRank,"]", destLineage)
          )
          next;
        } else {
          # rank column doesn't match src & dest ranks (but they match each other)
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "WARNING", "MOVE.RANK_WRONG", 
                           "Change=MOVE, but 'rank' column does not match 'current taxonomy' or 'proposed taxonomy'; use promote or demote, not move, to change rank", 
                           paste0("rank=", rankClean, 
                                  ", srcTaxonomy=[", srcTaxonRank,"]", srcLineage,
                                  ", destTaxonomy=[", destTaxonRank,"]", destLineage)
          )
        }
        
      }
      #
      # find the SRC target to MOVE in NEW and CUR
      #
      srcNewTarget =(.GlobalEnv$newMSL$name==as.character(srcTaxonName))
      srcPrevTarget=(.GlobalEnv$curMSL$taxnode_id==.GlobalEnv$newMSL[srcNewTarget,]$prev_taxnode_id)

      if(params$verbose) {print(paste0("MOVE: ",code," line ",linenum," '",srcTaxonName, "' findTarget(",srcTaxonName,")=",sum(srcNewTarget),"/",sum(srcPrevTarget)))}
    
      if(sum(srcNewTarget,na.rm=TRUE)==0) {
         errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                          "ERROR", "MOVE.NO_EXIST", "Change=MOVE, but taxon does not exist", 
                         paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName))
         next;
      } else if(sum(srcNewTarget,na.rm=TRUE)>1) {
         errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                          "ERROR", "MOVE.MANY", "Change=MOVE, multiple taxa exist with name", 
                         paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName
                                , " matches ", sum(srcNewTarget) ))
         next;
      }
      #
      # check if same accession number already exists
      #
      isDupAccession = (.GlobalEnv$newMSL$genbank_accession_csv == change$exemplarAccession & !srcNewTarget)
      if(sum(isDupAccession, na.rm=TRUE)>0) {
        errorDf=addError(errorDf,code,linenum,change$change,change$rank,change$.changeTaxon,
                         "ERROR", "MOVE.DUP_ACC", "Change=CREATE, a species with this accession number already exists", 
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
         errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                          "ERROR", "MOVE.PARENT_NO_EXIST", "Change=MOVE, but parent rank taxon does not exist", 
                         paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage)
         )
         next;
      } else if(sum(parentDestNewMatches,na.rm=TRUE)>1) {
         # this should never happen
         errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                          "ERROR", "MOVE.PARENT_MANY", "Change=MOVE, multiple taxa exist with parent name", 
                         paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage)
         )
         next;
      } else {
        # WARN: check if DEST rank matches change rank
        if ( tolower(destTaxonRank) != tolower(change$rank) ) {
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "WARNING", "MOVE.RANK_MISMATCH_DEST", 
                           "Change=MOVE, but proposed taxon rank does not match [rank] column", 
                           paste0("rankColumn=", change$rank,
                                  ", proposedTaxonRank=",destTaxonRank,
                                  ", proposedTaxonomy=", destLineage)
          )
        }
        # WARN: check if SRC rank matches change rank
        if ( tolower(srcTaxonRank) != tolower(change$rank) ) {
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon, 
                           "WARNING", "MOVE.RANK_MISMATCH_SRC", 
                           "Change=MOVE, but 'current taxonomy' taxon rank does not match [rank] column", 
                           paste0("rankColumn=", change$rank,
                                  ", currentTaxonRank=",srcTaxonRank,
                                  ", currentTaxonomy=", srcLineage)
          )
        }

        # 
        # MOVE the taxon
        #
        
        
        # get parent
        parentTaxon = .GlobalEnv$newMSL[parentDestNewMatches,]
        
        # WARN if PARENT_LINEAGE is not expected
        if( !is.na(destParentLineage) && (parentTaxon$lineage != destParentLineage) ) {
          errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                           "WARNING", "MOVE.PARENT_LINEAGE", 
                           "Change=MOVE, proposed parent taxon exists, but not with expected lineage", 
                           paste0("proposedParentLineage=", destParentLineage,
                                  ", observedParentLineage=",parentTaxon$lineage,
                                  ", otherProposals=",parentTaxon$prev_proposals)
          )
        }
        
        # add new info - primary columns
        .GlobalEnv$newMSL[srcNewTarget,"name"]       = destTaxonName
        # this part shouldn't be a change unless this is promote/demote....
        .GlobalEnv$newMSL[srcNewTarget,"level_id"]   = rankCV$id[rankCV$name==destTaxonRank]
        .GlobalEnv$newMSL[srcNewTarget,"rank"]       = destTaxonRank
        .GlobalEnv$newMSL[srcNewTarget,"parent_id"]  = parentTaxon[1,"taxnode_id"]


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
        
        # info-only columns - wont be saved to DB
        .GlobalEnv$newMSL[srcNewTarget,"lineage"]    = destLineage
 
        # QQQQQ RECURSE TO SET LINEAGE OF KIDS
        taxnode_id = .GlobalEnv$newMSL$taxnode_id[srcNewTarget]
        #browser()
        update_lineage(taxnode_id,destLineage)
      }
      #
      # update curMSL [out_*]
      #
      .GlobalEnv$curMSL[srcPrevTarget,"out_updated"] = TRUE  # admin; mark this to save to db
      .GlobalEnv$curMSL[srcPrevTarget,"out_change"] = "move"
      .GlobalEnv$curMSL[srcPrevTarget,"out_filename"] = proposalZip
      .GlobalEnv$curMSL[srcPrevTarget,"out_target"] = destLineage
      .GlobalEnv$curMSL[srcPrevTarget,"out_notes"] = paste0("linenum=",linenum) # add comments?
       
    } else {
      #  -------------------------------------------------------------------------
      #
      # ACTION not implemented.
      #
      #  -------------------------------------------------------------------------
      errorDf=addError(errorDf,code,linenum, change$change,change$rank,change$.changeTaxon,
                       "ERROR", "CHANGE.UNIMP", 
                       paste0("Change=",toupper(action)," is NOT (yet) implemented"), 
                         paste0("taxon=", srcTaxonName, ", lineage=",srcLineage,", proposedTaxon=", destTaxonName)
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
  # did proposal create empty (non-species) taxa?
  #
  kidCounts = apply(.GlobalEnv$newMSL[,"taxnode_id"],1,
                    function(parent_id,MSL) { sum(MSL$parent_id==parent_id)},
                    MSL=.GlobalEnv$newMSL)
  # scan for taxa with no kids, not species, not yet reported.
  emptyTaxa = .GlobalEnv$newMSL[kidCounts == 0,]%>% filter(level_id != 600) %>% filter(is.na(.emptyReported))
  if( nrow(emptyTaxa) > 0 ) {
    errorDf=addError(errorDf,code,NA, NA,NA,NA,
                     "ERROR", "PROPOSAL.EMPTY_TAXA", 
                     paste0("Proposal created empty (non-species) taxa"), 
                     paste0("the ",emptyTaxa$rank," '",emptyTaxa$name,"' is empty - it does not contain ny lower taxons")
    )
    # mark empty taxa so we don't re-report them
    .GlobalEnv$newMSL[.GlobalEnv$newMSL$name %in% emptyTaxa$name,".emptyReported"] = xlsxs[code,"xlsx"]
    #browser()
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
# extract current taxonomy
#
.GlobalEnv$oldMSLs = subset(taxonomyDt, msl_release_num < params$prev_msl)
.GlobalEnv$curMSL = subset(taxonomyDt, msl_release_num==params$prev_msl)
# add accounting columns
.GlobalEnv$curMSL[,"out_updated"] = FALSE

#
# copy prev MSL to make new MSL,
# to which we will try and apply these edits
#
.GlobalEnv$newMSL=createNewMSL(.GlobalEnv$curMSL,params$prev_msl, params$next_msl, params$taxnode_delta)

#
# DEBUG ONLY - delete all load/QC errors in allErrorDf
.GlobalEnv$allErrorDf = loadErrorDf %>% filter(TRUE)
# debug 
#cat("allErrorDf:",tracemem(allErrorDf),"\n");
summary(as.factor(errorDf$code)); summary(as.factor(allErrorDf$code));
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
    cat("# START PROC: ",code," with ", nrow(changeDf), " changes\n")
    
    # debug if code of interest
    #if( code %in% c("2022.043B.fixed") ) { cat("DEBUG ",code,"\n"); if(interactive()) {browser() }}
    
    
    # Apply the Change
    results=apply_changes(code,proposals[code,"basename"],changeDf)
    
    # Internal sanity check
    if(!(20070000 %in% .GlobalEnv$newMSL$ictv_id)) { 
      cat("!!!!!!!!!!!!!!!! LOST ROOT NODE in ",code," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      if(interactive()) {browser()}
    }
    errorDf = results[["errorDf"]]
    cat("# PROCed: ",code," with ",nrow(errorDf)," errors/warnings\n")
    if(nrow(errorDf)>0){.GlobalEnv$allErrorDf=rbindlist(list(.GlobalEnv$allErrorDf,errorDf),fill=TRUE)}
    
    processed=processed+1
  } # if have proposal XLSX
} # proposal code
cat("# DONE. Found",length(rownames(proposals)),"proposals; Processed",processed,", skipped",skipped,"\n")

#### SECTION write error list #####
# write current error list
# TODO: add line breaks, bold
# openxlsx
# r2excel: (install fails?) http://www.sthda.com/english/wiki/r2excel-read-write-and-format-easily-excel-files-using-r-software#install-and-load-r2excel-package
# xlsx (requires java)
# XLConnect (requires rJava)

prettyErrorDf = data.frame(allErrorDf %>% filter(FALSE))
prettyRow = 0
prevCode = ""
errorSortCols = allErrorDf[,c("code","row")]
errorSortCols$rown = as.integer(errorSortCols$row)
errorsSorted = do.call(order,errorSortCols[,c("code","rown")])

for(i in seq(1,nrow(allErrorDf)) ) {
  row=errorsSorted[i]
  # add blank line and header when document changes
  if(allErrorDf[row,"code"]!= prevCode) { 
    prevCode = allErrorDf[row,"code"]

    prettyErrorDf[prettyRow,c("folder")] = c(allErrorDf[row,"folder"])
    prettyRow=prettyRow+1
    
    prettyErrorDf[prettyRow,c("folder","code","xlsx")] = c(
      allErrorDf[row,"folder"],
      allErrorDf[row,"code"],
      ifelse(is.na(allErrorDf[row,"docx"]),
             allErrorDf[row,"xlsx"],
             allErrorDf[row,"docx"]
      )
    )
    prettyRow=prettyRow+1
    }
  
  # copy other lines as-is
  prettyErrorDf[prettyRow,]=allErrorDf[row,]
  prettyRow=prettyRow+1
}
# this should get chunked into files/worksheets by "folder"
for(committee in levels(as.factor(prettyErrorDf$folder)) ) {
  folderFilename = str_replace_all(str_replace(committee,pattern="(.*) \\(([A-Z])\\) .*","\\2 \\1")," ","_")
  filename = file.path(params$out_dir,
                       paste0("QC.pretty_summary.",folderFilename,".xlsx")
  )
  prettyCols = grep(names(prettyErrorDf),pattern="(code|docx|folder)",invert=T,value=T)
  
  select = (prettyErrorDf$folder == committee)
  write_xlsx( x=prettyErrorDf[select,prettyCols],path=filename)
  cat("Wrote: ", filename, " (",nrow(prettyErrorDf[select,]),"rows )\n")
}
filename = file.path(params$out_dir,"QC.pretty_summary.all.xlsx")
prettyCols = grep(names(prettyErrorDf),pattern="(code|docx)",invert=T,value=T)
write_xlsx( x=prettyErrorDf[,prettyCols],path=filename)
cat("Wrote: ", filename, " (",nrow(prettyErrorDf),"rows)\n")

write_xlsx( x=allErrorDf,                path=file.path(params$out_dir,"QC.summary.xlsx"))
write_delim(x=allErrorDf,                file=file.path(params$out_dir,"QC.summary.tsv"), delim="\t")
# ```
# 
# # known problems
# 
# ```{r final_stats}

dim(curMSL)
dim(newMSL)
dim(allErrorDf)

print("What is that ?space?: it's not a space, nor a tab:")
newMSL[grep(newMSL$lineage, pattern="Orthornavirae;"),c("taxnode_id","parent_id","lineage","name","cleaned_name","in_filename","in_notes")]
#                                                                                                         lineage         name cleaned_name
# 1: Riboviria;Orthornavirae ;Kitrinoviricota;Flasuviricetes;Amarillovirales;Flaviviridae;Pestivirus;Pestivirus L Pestivirus L   Pestivirus
# 2: Riboviria;Orthornavirae ;Kitrinoviricota;Flasuviricetes;Amarillovirales;Flaviviridae;Pestivirus;Pestivirus M Pestivirus M   Pestivirus
# 3: Riboviria;Orthornavirae ;Kitrinoviricota;Flasuviricetes;Amarillovirales;Flaviviridae;Pestivirus;Pestivirus N Pestivirus N   Pestivirus
# 4: Riboviria;Orthornavirae ;Kitrinoviricota;Flasuviricetes;Amarillovirales;Flaviviridae;Pestivirus;Pestivirus O Pestivirus O   Pestivirus
# 5: Riboviria;Orthornavirae ;Kitrinoviricota;Flasuviricetes;Amarillovirales;Flaviviridae;Pestivirus;Pestivirus P Pestivirus P   Pestivirus
# 6: Riboviria;Orthornavirae ;Kitrinoviricota;Flasuviricetes;Amarillovirales;Flaviviridae;Pestivirus;Pestivirus Q Pestivirus Q   Pestivirus
# 7: Riboviria;Orthornavirae ;Kitrinoviricota;Flasuviricetes;Amarillovirales;Flaviviridae;Pestivirus;Pestivirus R Pestivirus R   Pestivirus
# 8: Riboviria;Orthornavirae ;Kitrinoviricota;Flasuviricetes;Amarillovirales;Flaviviridae;Pestivirus;Pestivirus S Pestivirus S   Pestivirus

newMSL[taxnode_id==202203151, c("taxnode_id","parent_id","lineage","name","cleaned_name","in_filename","in_notes")]
#   taxnode_id parent_id                                                                                        lineage       name cleaned_name
#1:  202203151 202203068 Riboviria;Orthornavirae;Kitrinoviricota;Flasuviricetes;Amarillovirales;Flaviviridae;Pestivirus Pestivirus   Pestivirus
  
# fix it
#(gsub(pattern="[^a-zA-Z]",replacement="",x$kingdom))

#
# but need to scan and report the problems before fixing them!!!!

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
#   3. single-quote values
#
# ------------------------------------------------------------------------------
# cat(paste0('"',paste(names(newMSL),collapse='",\n"'),'"'))
colList = c("taxnode_id",
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
            "out_notes"
            #"lineage", # computed by trigger in db
            #"cleaned_name", # computed by trigger in db
            #"rank", # should be in level_id
            # "molecule", # should be in molecule_id
            # program admin columns - not in db
            #"out_updated",
            #"prev_taxnode_id",
            #"prev_proposals",
            # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
            #"host_source",
            # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
            #"exemplar_name",
            # NOTE: column does not (yet) exist in [taxonomy_node], only in [load_next_msl##]
            #"genome_coverage",
            #"comments" # doesn't exist in db, should be notes???
            #"lineage"  # computed by trigger in db
            )
newSqlFilename = "results/msl_load.sql"
sqlout=file(newSqlFilename,"wt")
cat("Writing ",newSqlFilename,"\n")

#
# output start transaction
# 
cat("begin transaction\n", file=sqlout)
cat("-- rollback transaction\n", file=sqlout)

#
# add new MSL to [taxonomy_toc]
#
cat("insert into [taxonomy_toc] ([tree_id],[msl_release_num],[comments]) ",
    "values (", paste(
     sort(levels(as.factor(newMSL$tree_id)),decreasing = T)[1],
      params$next_msl,
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
for( col in colList)  {
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
batchSize = 200 # rows
for(  row in order(newMSL$level_id) )  {
  #row=head(order(newMSL$level_id),n=1) # debug
  
  # insert several rows per batch insert  statement
  # much faster - fewer trigger calls
  # but makes localizing errors harder
  if( rowCount %% batchSize == 1 ) {
    cat(paste0("insert into [taxonomy_node] ",
               "([",
               paste0(colList,collapse="],["),
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
               ifelse(is.na(newMslStr[row, colList]), "NULL",
                      paste0("'",
                             # escape apostrophies as double-appostrophies (MSSQL)
                             sub(
                               "'", "''", newMslStr[row, colList]
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
where msl_release_num = ",params$next_msl,"
group by level_id
order by level_id
"), file=sqlout)

# QC query
cat(paste("-- QC Query
select 
     level_id,in_change, ct=count(*)
from taxonomy_node
where msl_release_num = ",params$next_msl,"
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
cat("begin transaction\n", file=sqlout)
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
                               sub("'","''", curMslStr[row,curMslColList])
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
exec rebuild_delta_nodes NULL

-- rebuild merge/split (seconds) \n
EXECUTE [dbo].[rebuild_node_merge_split] 

", file=sqlout)
cat(paste("
-- NOW check if all newMSL have delta in, and prevMSL have delta out
select 'prevMSL w/o delta to new', count(*) from taxonomy_node
where msl_release_num = ",params$next_msl-1,"
and taxnode_id not in (select prev_taxid from taxonomy_node_delta)
union all
select  'newMSL w/o delta to prev',  count(*) from taxonomy_node
where msl_release_num = ",params$next_msl,"
and taxnode_id not in (select new_taxid from taxonomy_node_delta)
"),
file=sqlout)

# QC query
cat(paste("-- QC Query
select 
     level_id,out_change, ct=count(*)
from taxonomy_node
where msl_release_num = ",params$next_msl-1,"
group by level_id,out_change
order by level_id,out_change

"), file=sqlout)
#
# close/flush file
#
close(sqlout)
cat("WROTE   ", newSqlFilename, "\n")

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

rdataFilename = paste0(params$proposals_dir,"/.RData")
save.image(file=rdataFilename)
cat("WROTE",rdataFilename,"\n")

