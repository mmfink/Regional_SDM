# File: 1_pointsInPolys_cleanBkgPts.r
# Purpose: 
# 1. Check input presence reach dataset for missing columns or data
# 2. Removing any reaches from the background dataset that are adjacent to
#  the input presence reach dataset.

library(RSQLite)
library(dplyr) # adding dplyr here as it may be causing downstream effects due to a bad interaction with 'sf'
library(sf)
library(stringr)


####
# Assumptions
# - the csv is named with the species code that is used in the lookup table (e.g. glypmuhl.shp)
# - There is lookup data in the sqlite database to link to other element information (full name, common name, etc.)
# - the csv has at least these fields EO_ID_ST, SNAME, SCOMNAME, COMID, OBSDATE, group_id, huc12

####
#### load input reaches ----
###
## two lines need your attention. The one directly below (loc_scripts)
## and about line 38 where you choose which file to use


# set the working directory to the location of the csv of species by reaches
setwd(loc_model)
dir.create(paste0(model_species,"/inputs/presence"), recursive = T, showWarnings = F)
dir.create(paste0(model_species,"/inputs/model_input"), showWarnings = F)
setwd(paste0(loc_model,"/",model_species,"/inputs"))
# changing to this WD temporarily allows for presence file to be either in presence folder or specified with full path name

# load data, QC ----
presReaches <- read.csv(nm_presFile)

#check for proper column names. If no error from next code block, then good to go
#presPolys$RA <- presPolys$SFRACalc
shpColNms <- names(presPolys)
desiredCols <- c("UID", "GROUP_ID", "SPECIES_CD", "RA", "OBSDATE")
# check if all required names are in file
if("FALSE" %in% c(desiredCols %in% shpColNms)) {
  stop(paste0("Column(s) are missing or incorrectly named: ", paste(desiredCols[!desiredCols %in% shpColNms], collapse = ", ")))
} else {
  print("Required columns are present")
}

# check if all columns have complete data
if(any(is.na(presPolys[,c("UID", "GROUP_ID", "SPECIES_CD", "RA")]))) {
  stop("The columns 'UID','GROUP_ID','SPECIES_CD', and 'RA' (SFRACalc) cannot have NA values.")
}

# check if file already exists; if it does, stop and print error
#if (!file.copy(nm_presFile, paste0(baseName, ".csv"))) {
#  stop("A file already exists with that name: '", 
#       paste0(getwd(), "/", baseName, ".csv"), "'. Rename or delete it to continue.")
#}
# set wd to inputs again
# setwd(paste0(loc_model,"/",model_species,"/inputs"))
# file in place, read it
#presReaches <- read.csv(paste0("presence/", baseName, ".csv"))

# arrange, pare down columns
presReaches <- presReaches[,desiredCols]

# set date/year column to [nearest] year, rounding when day is given
presReaches$OBSDATE <- as.character(presReaches$OBSDATE)
presReaches$date <- NA
for (d in 1:length(presReaches$OBSDATE)) {
  dt <- NA
  do <- presReaches$OBSDATE[d]
  if (grepl("^[0-9]{4}.{0,3}$", do)) {
    # year only formats
    dt <- as.numeric(substring(do,1,4))
  } else {
    if (grepl("^[0-9]{4}[-|/][0-9]{1,2}[-|/][0-9]{1,2}", do)) {
      # ymd formats
      try(dt <- as.Date(do), silent = TRUE)
    } else if (grepl("^[0-9]{1,2}[-|/][0-9]{1,2}[-|/][0-9]{4}", do)) {
      # mdy formats
      try(dt <- as.Date(do, format = "%m/%d/%Y"), silent = TRUE)
    }
    # if still no match, or if failed
    if (is.na(dt)) {
      if (grepl("[0-9]{4}", do)) {
        # use first 4-digit sequence as year
        dt <- regmatches(do, regexpr("[0-9]{4}", do))
        if (as.integer(dt) < 1900 | as.integer(dt) > format(Sys.Date(), "%Y")) {
          # years before 1900 and after current year get discarded
          dt <- Sys.Date()
        } else {
          dt <- as.Date(paste0(dt,"-01-01"))
        }
      }
      # put additional date formats here
    }
    # give up and assign current date
    if (is.na(dt)) {
      dt <- Sys.Date()
    }
    dt <- round(as.numeric(format(dt, "%Y")) + (as.numeric(format(dt,"%j"))/365.25))
  }
  presReaches$date[d] <- dt
}
desiredCols <- c(desiredCols, "date")

# just in case convert column names to lowercase
names(presReaches) <- tolower(names(presReaches))

###
# remove reaches from background dataset that have presence of the target species in the reach

# read in the shapefile, get the attribute data
dbEV <- dbConnect(SQLite(),dbname=nm_bkg[1])
SQLQuery <- paste0("SELECT * FROM ",nm_bkg[2]," WHERE COMID IN ('", paste(presReaches$comid, collapse = "','"),"')") 
shapef <- dbGetQuery(dbEV, SQLQuery)
SQLQuery <- paste0("SELECT proj4string p FROM lkpCRS WHERE table_name = '", nm_bkg[2], "';") 
proj4 <- dbGetQuery(dbEV, SQLQuery)$p
# shapef <- st_read(nm_allflowlines)
names(shapef) <- tolower(names(shapef))
shapef <- st_sf(shapef[c("comid", "huc12","wacomid")], geometry=st_as_sfc(shapef$wkt), crs=proj4)
# testcatchments <- shapef@data
shapef$huc12 <- str_pad(shapef$huc12, 12, pad=0)

# get huc12s, geom
pres.geom <- merge(shapef, presReaches, by = "comid")

# define project background
# subset background reaches by HUC2 to prevent predictions into basics where the species is not known to occur
presHUCs <- pres.geom$huc12

# test at what level HUCS are the same, and choose that level to run the predictions at.  
# For example, if all know occurences are within the same HUC6, then the study area will be clipped to that HUC6. 
# If they are not same at any level, then the model will be run at the full extant of the predictor layer.  
# THis is used to define the project background below.
if (is.null(huc_level)) {
  if(length(unique(substr(presHUCs,1,8)))==1){
    huc_level <- 8
  } else if(length(unique(substr(presHUCs,1,6)))==1){
    huc_level <- 6  
  } else if(length(unique(substr(presHUCs,1,4)))==1){
    huc_level <- 4  
  } else if(length(unique(substr(presHUCs,1,2)))==1){
    huc_level <- 2 
  } else {
    huc_level <- 4 # changed from 2 to try to narrow up the prediction area
  }
  fn_args$huc_level <- huc_level
  save(fn_args, file = paste0(loc_model, "/" , model_species, "/runSDM_paths.Rdata"))
}
message("Using huc_level of ", huc_level , "...")

# create background geom based on HUCsubset
dbEV <- dbConnect(SQLite(),dbname=nm_bkg[1])
if (huc_level != 0) {
  # subset to huc if requested
  HUCsubset <- unique(substr(presHUCs, 1, huc_level)) # subset to number of huc digits
  SQLQuery <- paste0("SELECT * FROM ",nm_bkg[2]," WHERE substr(huc12, 1, ", huc_level, ") IN ('", paste(HUCsubset, collapse = "','"),"');") 
} else {
  # otherwise take all reaches
  SQLQuery <- paste0("SELECT * FROM ",nm_bkg[2],";") 
}
shapef <- dbGetQuery(dbEV, SQLQuery)
names(shapef) <- tolower(names(shapef))
SQLQuery <- paste0("SELECT proj4string p FROM lkpCRS WHERE table_name = '", nm_bkg[2], "';") 
proj4 <- dbGetQuery(dbEV, SQLQuery)$p
shapef <- st_sf(shapef[c("comid", "huc12")], geometry = st_as_sfc(shapef$wkt), crs = proj4)


#EObyRA <- unique(shp_expl[,c("expl_id", "group_id","ra")])
shp_expl$minSamps[shp_expl$ra == "very high"] <- 5
shp_expl$minSamps[shp_expl$ra == "high"] <- 4
shp_expl$minSamps[shp_expl$ra == "medium"] <- 3
shp_expl$minSamps[shp_expl$ra == "low"] <- 2
shp_expl$minSamps[shp_expl$ra == "very low"] <- 1

shp_expl$finalSampNum <- ifelse(shp_expl$PolySampNum < shp_expl$minSamps, 
                                shp_expl$minSamps, 
                                shp_expl$PolySampNum)

ranPts <- st_sample(shp_expl, size = shp_expl$finalSampNum * 2)
ranPts.sf <- st_sf(ranPts)
names(ranPts.sf) <- "geometry"
st_geometry(ranPts.sf) <- "geometry"

ranPts.joined <- st_join(ranPts.sf, shp_expl)

# check for polys that didn't get any points
polysWithNoPoints <- shp_expl[!shp_expl$expl_id %in% ranPts.joined$expl_id,]
if(nrow(polysWithNoPoints) > 0){
  stop("One or more polygons didn't get any points placed in them.")
############### Aquatics ##################
# find presence and presence-adjacent reaches by intersection
#bkgd.int <- st_intersects(st_zm(shapef), st_zm(pres.geom) , sparse = F)
#bkgd.geom <- shapef[!apply(bkgd.int, 1, FUN = any),]
#if (length(bkgd.geom$geometry) > 3000) { 
#  bkgd.geom <- bkgd.geom[sort(sample(as.numeric(row.names(bkgd.geom)), size = 3000, replace = F)),]
}

# write species reach data
#st_write(pres.geom,paste("presence/", baseName,"_prepped.shp",sep=""))
# wtite background reach data
#st_write(bkgd.geom, paste("model_input/", baseName,"_bkgd_clean.shp",sep=""))
############### Aquatics ##################
#  remove extras using straight table work

# this randomly assigns digits to each point by group (stratum) then next row only takes 
# members in group that are less than target number of points
rndid <- with(ranPts.joined, ave(expl_id, stratum, FUN=function(x) {sample.int(length(x))}))
ranPts.joined2 <- ranPts.joined[rndid <= ranPts.joined$finalSampNum,]

# get actual finalSampNum
#ranPts.joined2 <- ranPts.joined[0,]
# #### this is slow! ####
# for (ex in 1:length(shp_expl$geometry)) {
#   s1 <- shp_expl[ex,]
#   samps <- row.names(ranPts.joined[ranPts.joined$expl_id==s1$expl_id,])
#   if (length(samps) > s1$finalSampNum) samps <- sample(samps, size = s1$finalSampNum) # samples to remove
#   ranPts.joined2 <- rbind(ranPts.joined2, ranPts.joined[samps,])
# }
ranPts.joined <- ranPts.joined2
#rm(ex, s1, samps, ranPts.joined2)
rm(rndid, ranPts.joined2)


#check for cases where sample smaller than requested
# how many points actually generated?
ptCount <- table(ranPts.joined$expl_id)
targCount <- shp_expl[order(shp_expl$expl_id),"finalSampNum"]
overUnderSampled <- ptCount - targCount$finalSampNum

#positive vals are oversamples, negative vals are undersamples
print(table(overUnderSampled))
# If you get large negative values then there are many undersamples and 
# exploration might be worthwhile

names(ranPts.joined) <- tolower(names(ranPts.joined))

colsToKeep <- c("stratum", tolower(desiredCols))
ranPts.joined <- ranPts.joined[,colsToKeep]

# name of random points output shapefile
nm.RanPtFile <- paste(loc_model, model_species,"inputs/presence",paste(baseName, "_RanPts.shp", sep = ""), sep = "/")
# write it out
st_write(ranPts.joined, nm.RanPtFile, driver="ESRI Shapefile", delete_layer = TRUE)

###
### remove Coincident Background points ----
###

# get range info from the DB (as a list of HUCs)
db <- dbConnect(SQLite(),dbname=nm_db_file)
SQLquery <- paste0("SELECT huc10_id from lkpRange
                   inner join lkpSpecies on lkpRange.EGT_ID = lkpSpecies.EGT_ID
                   where lkpSpecies.sp_code = '", model_species, "';")
hucList <- dbGetQuery(db, statement = SQLquery)$huc10_id
dbDisconnect(db)
rm(db)

op <- options()
options(useFancyQuotes = FALSE) #need straight quotes for query
# get the background data from the DB
db <- dbConnect(SQLite(), nm_bkgPts[1])
qry <- paste0("SELECT * from ", nm_bkgPts[2], " where substr(huc12,1,10) IN (", paste(sQuote(hucList), collapse = ", ", sep = "")," );")
bkgd <- dbGetQuery(db, qry)
tcrs <- dbGetQuery(db, paste0("SELECT proj4string p from lkpCRS where table_name = '", nm_bkgPts[2], "';"))$p
samps <- st_sf(bkgd, geometry = st_as_sfc(bkgd$wkt, crs = tcrs))
options(op)
rm(op)

# reduce the number of bkg points if huge
# use the greater of 20 * pres points or 50,000
bkgTarg <- max(nrow(ranPts.joined) * 20, 50000)
if(nrow(samps) > bkgTarg){
  samps <- samps[sample(nrow(samps), bkgTarg),]
}

# find coincident points ----
polybuff <- st_transform(shp_expl, st_crs(samps))
polybuff <- st_buffer(polybuff, dist = 30)

coincidentPts <- unlist(st_contains(polybuff, samps, sparse = TRUE))

# remove them (if any)
if (length(coincidentPts) > 0) backgSubset <- samps[-coincidentPts,] else backgSubset <- samps

#write it up and do join in sqlite (faster than downloading entire att set)
st_geometry(backgSubset) <- NULL
tmpTableName <- paste0(nm_bkgPts[2], "_", baseName)
dbWriteTable(db, tmpTableName, backgSubset, overwrite = TRUE)

# do the join, get all the data back down
qry <- paste0("SELECT * from ", tmpTableName, " INNER JOIN ", nm_bkgPts[2], "_att on ",
              tmpTableName,".fid = ", nm_bkgPts[2], "_att.fid;")
bgSubsAtt <- dbGetQuery(db, qry)
# delete the table on the db
dbRemoveTable(db, tmpTableName)
dbDisconnect(db)

# # remove NAs #### convert this to removing columns instead?
# bgSubsAtt$bulkDens <- as.numeric(bgSubsAtt$bulkDens)
# bgSubsAtt$clay <- as.numeric(bgSubsAtt$clay) 
# bgSubsAtt$soil_pH <- as.numeric(bgSubsAtt$soil_pH) 
# bgSubsAtt$flowacc <- as.numeric(bgSubsAtt$flowacc) 
# 
# bgSubsAtt <- bgSubsAtt[complete.cases(bgSubsAtt),]
# nrow(bgSubsAtt)

dbName <- paste0(loc_model, "/", model_species, "/inputs/model_input/", baseName, "_att.sqlite")
db <- dbConnect(SQLite(), dbName)
dbWriteTable(db, paste0(nm_bkgPts[2], "_clean"), bgSubsAtt, overwrite = TRUE)

dbDisconnect(db)

rm(db, do, dt, qry, tmpTableName)
