####  Moisture response curve experiment - script - Juan F. Dueñas juanfduenas@proton.me

####### Load Packages we will work with. If you have not installed them yet, do it.
pkgs <- c("tidyverse", "mgcv", "gratia", "DHARMa", "ggpubr")

vapply(pkgs, FUN = library, FUN.VALUE = logical(1L), logical.return = TRUE, character.only = TRUE)

#~~~~ create path object ~ 
path <- "/home/juanfduenas/Documents/Cliwac" # Linux desktop

# Figure 1 ------ Control soil

mwd <- read_delim(paste(path,"/SR/mwd.csv", sep = ""), delim = ",") %>% # read .csv file - 
       filter(treatment%in%"Ctrl")%>%
       select(potID, moistureL, whc, MWD)%>%
       rename(dl=moistureL)

wsa <- read_delim(paste(path,"/SR/wsa.csv", sep = ""), delim = ",") %>% # read .csv file - 
       filter(treatment%in%"ctrl")%>%
       select(potID, WSA)

imp <- read_delim(paste(path,"/SR/imp.csv", sep = ""), delim = ",")%>%
       filter(treatment%in%"ctrl")%>%
       select(potID, wdpt)
