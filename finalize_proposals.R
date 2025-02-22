#
# prepare proposal files for posting
#
# 1. clean up names (remove .v#[.fix])
# 2. copy to final directory structure w/ filenames
# 3. zip up
# 4. copy zips to final merged dir structure. 
# -------------------------------------------------------
# Naming convention
# -------------------------------------------------------
# The naming convention (unofficial) is to use a period to 
# separate different “fields” in the file name, with the fields being 
#    Year.ProposalNumberSC.status.version.TaxonName_SummaryChanges
# where
#   N: New
#   U: Unaccepted (but may be reconsidered)            
#   A: Accepted
#   R: Ratified
# -------------------------------------------------------

suppressPackageStartupMessages(library(dplyr))

#### Parse Args ####
suppressPackageStartupMessages(library("optparse"))
option_list <- list( 
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output"),
  make_option(c("-m","--mode"), default="scan", dest="mode",
              help=paste("build mode: scan, build, rezip", 
                         "\n\t\t\tscan: just look for files",
                         "\n\t\t\tbuild: copy and rename files, then zip",
                         "\n\t\t\trezip: just re-zip already renamed files")),
  # inputs
  make_option(c("-s","--src"), default="pending/all_proposals", dest="src_dir",
              help="source proposals dir"),
  # outputs
  make_option(c("-o","--out"), default="pending/proposalsFinal", dest="dest_dir",
              help="output (renamed) proposals dir"),
  make_option(c("-z","--zip"), default="pending/proposalsFinalZips", dest="zips_dir",
              help="output (renamed) proposals dir"),
  make_option(c("-t","--tmp"), default="pending/proposalsFinalZips/tmp", dest="ztmp_dir",
              help="temp dir for staging files to zip")
  
)
params <- parse_args(OptionParser(option_list=option_list))

# needs to be externalized to a file, so can be shared with merge_prosal_zips.R
sc2destFolder = c(
  "S"="Animal +ssRNA (S) proposals",
  "D"="Animal DNA viruses and Retroviruses (D) proposals",
  "M"="Animal dsRNA and -ssRNA (M) proposals",
  "A"="Archaeal viruses (A) proposals",
  "B"="Bacterial viruses (B) proposals",
  "F"="Fungal and protist virus (F) proposals",
  "G"="General (G) proposals",
  "P"="Plant virus (P) proposals"
)

status2text = c(
  "N"="New",
  "U"="Unaccepted (but may be reconsidered)" ,          
  "Ud"="Unaccepted (diferred)" ,          
  "A"="Accepted",
  "R"="Ratified"
)

#
# scan for proposal codes
#
proposals = data.frame(path=list.files(path=params$src_dir,pattern="2[0-9][0-9][0-9]\\.[0-9]{3}[A-Z][X]{0,1}\\..*\\.(doc|docx|xls|ppt|xlsx|pptx|pdf|png|zip)$", recursive=T, full.names=TRUE) )
cat(paste0("SCANNED ",params$src_dir, "/ FOUND ",dim(proposals)[1]," proposal documents\n"))
proposals$filename = gsub("^.*/","",proposals$path)

# filter editor temp files
proposals = proposals[grep("^~",proposals$filename, invert=T),]

# extract code
proposals$code = gsub("^([0-9]+\\.[0-9]+[A-Z]X*)\\.([NUAR]).*$","\\1",proposals$filename)
proposals$status = gsub("^([0-9]+\\.[0-9]+[A-Z]X*)\\.([NUAR]).*$","\\2",proposals$filename)
proposals$sc =   gsub("^2[0-9][0-9][0-9]\\.[0-9]{3}([A-Z])[X]{0,1}$","\\1",proposals$code)

# get list of uniq codes
allCodes = unique(proposals[,c("code","sc")])
allCodes$sc = gsub("^.*([A-Z])X*$","\\1",allCodes$code)
rownames(allCodes)=allCodes$code

# QC
badSC = ! proposals$sc %in% names(sc2destFolder)
if( sum(badSC) > 0) {
  cat(paste0("### ERROR: ",sum(badSC)," codes include terminal letters that aren't Study Section abbreviations:","\n"))
  cat(paste0(proposals[badSC,c("sc","code","filename","path")],"\n"))
  return(1)     
}
badStatus = ! proposals$status %in% names(status2text)
if( sum(badStatus) > 0) {
  cat(paste0("### ERROR: ",sum(badStatus)," filenames include invalid status lettters (valid status:", paste(paste0(names(status2text),"='",status2text,"'"),collapse=","),"):","\n"))
  cat(paste0(proposals[badStatus,c("status","code","filename","path")],"\n"))
  return(1)     
}
# filter out "Unaccepted*"
proposals = proposals %>% filter(!(status %in% c("U","Ud")))
# generate "cleaned" filenames
# remove .A. 
# remove .v#
# remove .fix
#proposals$cleanFilename = gsub("^([0-9]+\\.[0-9]+[A-Z]X*\\.[A-Z]+)\\.v[0-9]+(\\.fix)*(\\..*)$","\\1\\3",proposals$filename)
proposals$finalFilename = gsub("^([0-9]+\\.[0-9]+[A-Z]X*)\\.[A-Z]+(\\.v[0-9]+)*(\\.fix)*(\\..*)$","\\1\\4",proposals$filename)

#
# END SCAN
# 
if( params$mode == "scan" ) {
  print("#SCAN complete.")
  return(0)
}

#
# remove and re-create target directory structure
#
# remove dirs
if( params$mode != "rezip" ) {
  cat("CMD: rm -rf", params$dest_dir, "\n")
  system(paste0("rm -rf '",params$dest_dir,"'"), intern=F,ignore.stdout=F, ignore.stderr=F,wait=T)
}
cat("CMD: rm -rf", params$zips_dir, "\n")
system(paste0("rm -rf '",params$zips_dir,"'"), intern=F,ignore.stdout=F, ignore.stderr=F,wait=T)
# re-create dirs
for(dirPath in c(params$dest_dir, params$zips_dir, params$ztmp_dir, paste0(params$dest_dir,"/",sc2destFolder))) {
  if (!dir.exists(dirPath)){
    dir.create(dirPath)
    if(params$verbose) {print(paste0("mkdir ", dirPath))}
  }
}

#
# for each CODE
#
# 1. copy files into SC folder, and clean names
# 2. create zip file in download folder
#
proposals$cleanPath = paste0(params$dest_dir,"/",sc2destFolder[proposals$sc],"/",proposals$finalFilename)
proposals$finalPath = paste0(params$ztmp_dir, "/",proposals$finalFilename)
for(curCode in allCodes$code) {
  cat(paste0("### code=[[",curCode,"]] ###\n"))

  # find fileset for proposal (could include Supp files)
  codeFiles = proposals %>% filter(code == curCode)
  cat(paste0("file count=",nrow(codeFiles),"\n"))
  
  for(fileRow in rownames(codeFiles)) {
    # copy to clean folders area
    if(params$mode != "rezip") {
      cmdStr = paste0("       cp -a '",codeFiles[fileRow,]$path, "' '", codeFiles[fileRow,]$cleanPath,"'")
      if(params$verbose) {cat(paste0("   CMD: ", cmdStr, "\n"))}
      system(cmdStr)
    }

    #  copy to tmp area with finalFilename, instead of cleanFilename
    cmdStr = paste0("       cp -a '",codeFiles[fileRow,]$path, "' '", codeFiles[fileRow,]$finalPath,"'")
    if(params$verbose) {cat(paste0("   CMD: ", cmdStr,"\n"))}
    system(cmdStr)
  }
  
  # create zip files
  # -j : filename only, no path stored
  # -df : [MacOS] store only data-fork of file, for export to other FS.
  primaryFile = codeFiles[grep(pattern="Suppl", codeFiles$filename, invert = T),]
  primaryFile = primaryFile[grep(pattern=".docx$",primaryFile$filename),]
  zipFilename = paste0(params$zips_dir,"/",sub(".docx$",".zip",primaryFile$finalFilename))
  cmdStr = paste0("zip -j '",zipFilename,"' '",
                  paste0(codeFiles$finalPath,collapse="' '"),
                  "'")
  if(params$verbose) {cat(paste0(" CMD: ", cmdStr, "\n"))}
  system(cmdStr)
}

#
# clean up tmp dir
#
cat(paste0("CMD: rm -rf ", params$ztmp_dir),"\n")
system(paste0("rm -rf '",params$ztmp_dir,"'"), intern=F,ignore.stdout=F, ignore.stderr=F,wait=T)
