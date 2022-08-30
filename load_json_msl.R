#
# JSON: read current_taxonomy
# 

#
# RJSON - produces nested list
# 
library("rjson")
nonSpecies = fromJSON(file="current_msl/nonSpeciesTaxa_083022.json")
names(nonSpecies)
# [1] "has_species"    "name"           "parentDistance" "rankIndex"      "rankName"      
# [6] "taxNodeID"      "children" 


#
# JSONlite 
# 
library("jsonlite")
library("readr")
fname="current_msl/nonSpeciesTaxa_083022.json"
nonSpecies = jsonlite::fromJSON(read_file(fname))
names(nonSpecies)
#
