library(HMDHFDplus)
library(curl)
library(XLConnect)
library(tibble)
library(magrittr)
library(dplyr)


config <- yaml::read_yaml("hmd_credentials.yaml")

if (config$HMD_username=="YOUR_USERNAME" | config$HMD_password=="YOUR_PASSWORD"){
  cat("To download the required data you must provide your Human Mortality Database \n
       username and password in as fields in the hmd_credentials.yaml file \n")
  cat("If you do not have a username you can set one up at http://www.mortality.org/")
}


## download HMD data --------------------------------------------------------------------
# UK data (not GBR as code suggest) - single year of age and single year exposure data
exp_hmd <- readHMDweb(CNTRY = "GBR_NP", item = "Exposures_1x1", fixup = TRUE,
                      username = config$HMD_username,
                      password = config$HMD_password)

# UK death data single year of age and single year exposure data
deaths_hmd  <- readHMDweb(CNTRY = "GBR_NP", item = "Deaths_1x1", fixup = TRUE,
                          username = config$HMD_username,
                          password = config$HMD_password)

dir.create(file.path("data", "hmd"), recursive=T)
saveRDS(exp_hmd, "data/hmd/exposures_hmd.rds")
saveRDS(deaths_hmd, "data/hmd/deaths_hmd.rds")


# download ons files ---------------------------------------------------------------

url <- paste0("https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/",
              "populationandmigration/populationprojections/datasets/",
              "z1zippedpopulationprojectionsdatafilesuk/2014based/",
              "rft-table-z1-open-data---uk.zip")
dl_folder <- file.path("data", "ons")
dir.create(dl_folder, recursive = T)
curl_download(url, file.path(dl_folder, "NPP14.zip"))
unzip(zipfile = file.path(dl_folder, "NPP14.zip"), exdir=dl_folder)
# Each of the xml files unzipped here contain one forecast scenario. 
# The files with names containing plp, ppp and php are the low, principal, and 
# high scenarios respectively.
# The sheet "mortality assumptions" contains the required qx forecasts. 
# Because of difficulties in parsing the xml file into a R dataframe, the 
# relevant sheets had to be saved out manually to an csv file through excel
# This is the only manual step required, and the outputed csv files are included
# with this code. 
# The 5 lines of code above are not therefore required, but are included to document
# the source of the data.


## life expectancy -------------------------------------------------------------
# this section of the code is reliant on the urls and formatting of the below ONS 
# data files remaining the same.

# life expectancy forecast, 2014 based, principal projections
url <- paste0("https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/",
              "birthsdeathsandmarriages/lifeexpectancies/datasets/",
              "expectationoflifeprincipalprojectionunitedkingdom/2014based/",
              "wukprincipal14.xls")


curl_download(url, file.path(dl_folder, "wukprincipal14.xls"))

# life expectancy forecast, 2014 based, low e0 projections
url <- paste0("https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/",
              "birthsdeathsandmarriages/lifeexpectancies/datasets/",
              "expectationoflifelowlifeexpectancyvariantunitedkingdom/2014based/",
              "wuklle14.xls")
curl_download(url,file.path(dl_folder, "wuklle14.xls"))

# life expectancy forecast, 2014 based, high e0 projections
url <- paste0("https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/",
              "birthsdeathsandmarriages/lifeexpectancies/datasets/",
              "expectationoflifehighlifeexpectancyvariantunitedkingdom/2014based/",
              "wukhle14.xls")

curl_download(url,file.path(dl_folder,"wukhle14.xls"))

# read in data ons data and format --------------------------------------------------
mp_le <- readWorksheetFromFile("data/ons/wukprincipal14.xls",
                               sheet=3, startRow=9)
mp_le_df <- tibble(Year=as.numeric(mp_le[1,2:dim(mp_le)[2]]),
                   e0=as.numeric(mp_le[3,2:dim(mp_le)[2]]))

fp_le <- readWorksheetFromFile("data/ons/wukprincipal14.xls",
                               sheet=4, startRow=8)

fp_le_df <- tibble(Year=as.numeric(fp_le[2,2:dim(fp_le)[2]]),
                   e0=as.numeric(fp_le[4,2:dim(fp_le)[2]]))


mh_le <- readWorksheetFromFile("data/ons/wukhle14.xls",
                               sheet=3, startRow=8)
mh_le_df <- tibble(Year=as.numeric(mh_le[2,2:dim(mh_le)[2]]),
                   e0=as.numeric(mh_le[4,2:dim(mh_le)[2]]))

ml_le <- readWorksheetFromFile("data/ons/wuklle14.xls",
                               sheet=3, startRow=8)
ml_le_df <- tibble(Year=as.numeric(ml_le[2,2:dim(ml_le)[2]]),
                   e0=as.numeric(ml_le[4,2:dim(ml_le)[2]]))

me0_df <- rbind(ml_le_df %>% mutate(`ONS Variant` = "Low"),
                mp_le_df %>% mutate(`ONS Variant` = "Principal"),
                mh_le_df %>% mutate(`ONS Variant` = "High")) %>%
  mutate(Sex="Male")


fh_le <- readWorksheetFromFile("data/ons/wukhle14.xls",
                               sheet=4, startRow=8)
fh_le_df <- tibble(Year=as.numeric(fh_le[2,2:dim(fh_le)[2]]),
                   e0=as.numeric(fh_le[4,2:dim(fh_le)[2]]))

fl_le <- readWorksheetFromFile("data/ons/wuklle14.xls",
                               sheet=4, startRow=8)
fl_le_df <- tibble(Year=as.numeric(fl_le[2,2:dim(fl_le)[2]]),
                   e0=as.numeric(fl_le[4,2:dim(fl_le)[2]]))

fe0_df <- rbind(fl_le_df %>% mutate(`ONS Variant` = "Low"),
                fp_le_df %>% mutate(`ONS Variant` = "Principal"),
                fh_le_df %>% mutate(`ONS Variant` = "High")) %>%
  mutate(Sex="Female")

# combine
e0_df <- rbind(fe0_df, me0_df)

# save out data --------
saveRDS(e0_df, file = "data/ons/e0_df.rds")
