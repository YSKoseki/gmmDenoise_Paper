# 09_SampleData_DRA005106.R
# Last updated on 2022.4.16 by YK
# An R script to make a sample data table from SRA BioSample xml document
# # R 4.1.2
# 
# To get SRA BioSample xml document:
# 1. Go to the page:
#    https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJDB5158
# 2. Click the 'BioSample' link on the right side of the page, which opens the page: 
#    https://www.ncbi.nlm.nih.gov/biosample?LinkName=bioproject_biosample_all&from_uid=388425
# 3. Click the 'Send to:▼' and choose 'File' in the new window.
# 4. Select 'Full XML (text)' as the format and press 'Create File'
# Reference: https://www.biostars.org/p/279582/

# Packages required
library(xml2); packageVersion("xml2") # 1.3.3
library(tidyverse); packageVersion("tidyverse") # 1.3.1

# Input and output
## Input: metadata.csv obtained by the shell script "01_SRAfastq_DRA005106.sh"
path_input_csv <- "./01-SRAfastq_DRA005106/metadata.csv"
## Input: biosample_result.xml downloaded at https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJDB5158
path_input_xml <- "./01-SRAfastq_DRA005106/biosample_result.xml"
## Output
path_output <- "./09-SampleDat_DRA005106"

# Making sample data table
## Extract river names, latitudes, and longitudes from the BioSample data
biosample <- read_xml(path_input_xml)
BioSample <- xml_find_all(biosample, "//BioSample") %>% xml_attr("accession")
location <- xml_find_all(biosample, "///Attributes") %>% xml_text() %>% 
  str_split("Japan:", simplify=TRUE) %>% 
  `[`(, 2) %>% 
  str_split("stream water", simplify=TRUE) %>% 
  `[`(, 1)
river <- str_replace(location, "River", "River;") %>% 
  str_split(";", simplify=TRUE) %>% 
  `[`(, 1) %>% 
  if_else(.=="", NA_character_, .) %>% 
  str_replace_all("^ ", "") # fixes some data having a blank before character string 
latlon <- str_split(location, "River", simplify=TRUE) %>%
  `[`(, 2) %>%
  str_split(" ", simplify=TRUE) %>% 
  `[`(, c(1, 3))
##  Tabulate the data
sample_tib <- bind_cols(biosample=BioSample,
                        river=river,
                        lat=as.numeric(latlon[, 1]), lon=as.numeric(latlon[, 2]))
## Add a column of sample (run) IDs to the table 
metadat_tib <- read_csv(path_input_csv)
sample_tib <- mutate(metadat_tib, biosample=BioSample, sample=Run) %>% 
  select(biosample, sample) %>% 
  left_join(sample_tib, by="biosample") %>% 
  select(-biosample) %>% 
  arrange(sample)

# Write csv file
dir.create(path_output, recursive=TRUE)
write.csv(sample_tib, paste0(path_output, "/02-sample_tib.csv"))

# Save data
## Save R objects
path_saved_object <- paste0(path_output, "/01-Saved_object")
dir.create(path_saved_object, recursive=TRUE)
saveRDS(sample_tib, paste0(path_saved_object, "/sample_tib", ".obj"))
## Save workspace and session info
save.image(paste0(path_output, "/sampledat.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_output, "/sampledat.info"))
