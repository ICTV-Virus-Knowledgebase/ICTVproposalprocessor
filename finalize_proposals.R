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

library(dplyr)

MSL_dir = "MSL39v6"
# if true, rebuild just zips
# if false, rebuild from Pending_Proposals/
preserve_final = F


params = list(
  # inputs
  proposals_dir=file.path(MSL_dir,"Pending_Proposals"),
  # outputs
  dest_dir     =file.path(MSL_dir,"proposalsFinal"),
  download_dir =file.path(MSL_dir,"proposalsFinalZips"),
  # temp files
  tmp_dir      =file.path(MSL_dir,"proposalsFinalZips/tmp")
  
)
 
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
# remove and re-create target directory structure
#
# remove dirs
if( !preserve_final ) {
  cat("CMD: rm -rf", params$dest_dir, "\n")
  system(paste0("rm -rf '",params$dest_dir,"'"), intern=F,ignore.stdout=F, ignore.stderr=F,wait=T)
}
cat("CMD: rm -rf", params$download_dir, "\n")
system(paste0("rm -rf '",params$download_dir,"'"), intern=F,ignore.stdout=F, ignore.stderr=F,wait=T)
# re-create dirs
for(dirPath in c(params$dest_dir, params$download_dir, params$tmp_dir, paste0(params$dest_dir,"/",sc2destFolder))) {
  if (!dir.exists(dirPath)){
      dir.create(dirPath)
      print(paste0("mkdir ", dirPath))
  }
}
#
# scan for proposal codes
#
proposals = data.frame(path=list.files(path=params$proposals_dir,pattern="2[0-9][0-9][0-9]\\.[0-9]{3}[A-Z][X]{0,1}\\..*\\.(doc|docx|xls|ppt|xlsx|pptx|pdf|png|zip)$", recursive=T, full.names=TRUE) )
print(paste0("SCANNED ",params$proposals_dir, "/ FOUND ",dim(proposals)[1]," proposal documents"))
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
  print(paste0("### ERROR: ",sum(badSC)," codes include terminal letters that aren't Study Section abbreviations:"))
  print(proposals[badSC,c("sc","code","filename","path")])
  return(1)     
}
badStatus = ! proposals$status %in% names(status2text)
if( sum(badStatus) > 0) {
  print(paste0("### ERROR: ",sum(badStatus)," filenames include invalid status lettters (valid status:", paste(paste0(names(status2text),"='",status2text,"'"),collapse=","),"):"))
  print(proposals[badStatus,c("status","code","filename","path")])
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
# for each CODE
#
# 1. copy files into SC folder, and clean names
# 2. create zip file in download folder
#
proposals$cleanPath = paste0(params$dest_dir,"/",sc2destFolder[proposals$sc],"/",proposals$finalFilename)
proposals$finalPath = paste0(params$tmp_dir, "/",proposals$finalFilename)
for(curCode in allCodes$code) {
  print(paste0("### code=[[",curCode,"]] ###"))

  # find fileset for proposal (could include Supp files)
  codeFiles = proposals %>% filter(code == curCode)
  print(paste0("file count=",nrow(codeFiles)))
  
  for(fileRow in rownames(codeFiles)) {
    # copy to clean folders area
    if(!preserve_final) {
      cmdStr = paste0("       cp -a '",codeFiles[fileRow,]$path, "' '", codeFiles[fileRow,]$cleanPath,"'")
      print(paste0("   CMD: ", cmdStr))
      system(cmdStr)
    }

    #  copy to tmp area with finalFilename, instead of cleanFilename
    cmdStr = paste0("       cp -a '",codeFiles[fileRow,]$path, "' '", codeFiles[fileRow,]$finalPath,"'")
    print(paste0("   CMD: ", cmdStr))
    system(cmdStr)
  }
  
  # create zip files
  # -j : filename only, no path stored
  # -df : [MacOS] store only data-fork of file, for export to other FS.
  primaryFile = codeFiles[grep(pattern="Suppl", codeFiles$filename, invert = T),]
  primaryFile = primaryFile[grep(pattern=".docx$",primaryFile$filename),]
  zipFilename = paste0(params$download_dir,"/",sub(".docx$",".zip",primaryFile$finalFilename))
  cmdStr = paste0("zip -j '",zipFilename,"' '",
                  paste0(codeFiles$finalPath,collapse="' '"),
                  "'")
  print(paste0("   CMD: ", cmdStr))
  system(cmdStr)
}

#
# clean up tmp dir
#
print(paste0("CMD: rm -rf ", params$tmp_dir))
system(paste0("rm -rf '",params$tmp_dir,"'"), intern=F,ignore.stdout=F, ignore.stderr=F,wait=T)
