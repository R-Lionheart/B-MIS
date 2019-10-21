library(data.table)
library(ggplot2)
library(parallel)
library(stringr)
library(tidyverse)
options(scipen=999)

## This BMIS is for Wei's Eddy Transect data.

# Things to return 
# IS_inspectPlot (plot to make sure there aren't any internal standards we should kick out)
# QuickReport (% that picked BMIS, with cut off values)
# ISTest_plot (plot to evaluate if you cut off is appropriate)
# BMIS_normalizedData (tibble with the info you actually want!)


# Imports -----------------------------------------------------------------

Wei.SampKey_all <- read.csv("data/Sample.Key.EddyTransect.csv") 
Wei.Internal.Standards <- read.csv("data/Ingalls_Lab_Standards.csv") %>%
  filter(Column == "HILIC") %>%
  filter(z == 1)

# Positive data only. Formerly known as xcms.dat_pos
Wei.transect.pos <- read.csv("data/Wei_Transect_QC.csv", header = TRUE) %>% 
  slice(-1:-6) %>%
  select(-c(Description, Value)) 


# Change class + adjust data. Set cutoff values -----------------------------------------------------------------
Wei.transect.pos <- Wei.transect.pos %>%
  filter(!str_detect(ReplicateName, "Blk")) %>%
  filter(!str_detect(ReplicateName, "Std")) %>%
  mutate(ReplicateName = as.character(ReplicateName)) %>%
  mutate(Metabolite.name = as.character(Metabolite.name)) %>%
  mutate(RTValue = as.numeric(RTValue)) %>%
  mutate(AreaValue = as.numeric(AreaValue)) %>%
  mutate(SNValue = as.numeric(SNValue)) %>% # ORIGINAL STOPS HERE- BELOW IS REMOVING THE INGALLS_ PREFIX FROM METABOLITE NAME
  mutate(Metabolite.name = ifelse(str_detect(Metabolite.name, "Ingalls_"), sapply(strsplit(Metabolite.name, "_"), `[`, 2), Metabolite.name))


cut.off <- 0.4 # 40% decrease in RSD of pooled injections, aka improvement cutoff
cut.off2 <- 0.1 # RSD minimum

# Match transect data with Internal Standards list -----------------------------------------------------------------
Wei.transect.withIS <- Wei.transect.pos %>%
  filter(Metabolite.name %in% Wei.Internal.Standards$Compound.Name) 

Wei.transect.NoIS <- Wei.transect.pos %>%
  filter(!Metabolite.name %in% Wei.Internal.Standards$Compound.Name) 


# Read in Internal Standard data -----------------------------------------------------------------
# If injection volume is known, add in here.
Wei.IS.data <- Wei.transect.withIS %>%
  # select(ReplicateName, Metabolite.name, AreaValue) %>% # Original, non-QC'd AreaValue 
  select(ReplicateName, Metabolite.name, Area.with.QC) %>%
  mutate(MassFeature = Metabolite.name) %>%
  select(-Metabolite.name) 

Wei.IS.data$ReplicateName <- gsub("^.{0,1}", "", Wei.IS.data$ReplicateName)
  

Wei.SampKey <- Wei.SampKey_all %>%
  filter(Sample.Name %in% Wei.IS.data$ReplicateName) %>% # Drops standards from SampKey_all
  select(Sample.Name, Bio.Normalization) %>%
  # filter(!is.na(Bio.Normalization)) %>% # Unnecessary for transect dataset
  mutate(MassFeature = "Inj_vol",
         Area.with.QC = Bio.Normalization,
         ReplicateName = Sample.Name) %>%
  select(ReplicateName, Area.with.QC, MassFeature)

Wei.IS.data <- rbind(Wei.IS.data, Wei.SampKey) 

# Extraction replication of Internal Standards -----------------------------------------------------------------
IS_inspectPlot <- ggplot(Wei.IS.data, aes(x = ReplicateName, y = AreaValue)) + 
  geom_bar(stat = "identity") + 
  facet_wrap( ~MassFeature, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 5), 
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10))+
  ggtitle("IS Raw Areas")
#print(IS_inspectPlot)


# Edit data so names match-----------------------------------------------------------------
Wei.IS.data <- Wei.IS.data %>% 
  mutate(ReplicateName = ReplicateName %>%
           str_replace("-",".")) 

Wei.transect.long  <- Wei.transect.NoIS %>%
  rename(MassFeature = Metabolite.name) %>%
  select(ReplicateName, MassFeature, Area.with.QC)

Wei.transect.long$ReplicateName <- gsub("^.{0,1}", "", Wei.transect.long$ReplicateName)


# Caluclate mean values for each IS----------------------------------------------------------------
Wei.IS.means <- Wei.IS.data %>% 
  filter(!grepl("_Blk_", ReplicateName)) %>%
  mutate(MassFeature = as.factor(MassFeature)) %>%
  group_by(MassFeature) %>%
  summarise(Average.Area = mean(as.numeric(Area.with.QC), na.rm = TRUE)) %>%
  mutate(MassFeature = as.character(MassFeature))


# Normalize to each internal Standard----------------------------------------------------------------
Wei.binded <- rbind(Wei.IS.data, Wei.transect.long) %>%
  arrange(MassFeature)

Split_Dat <- list()

for (i in 1:length(unique(Wei.IS.data$MassFeature))) {
  Split_Dat[[i]] <- Wei.binded %>% 
    mutate(MIS = unique(Wei.IS.data$MassFeature)[i]) %>%
    left_join(Wei.IS.data %>% 
                rename(MIS = MassFeature, IS_Area = AreaValue) %>% 
                select(MIS, ReplicateName, IS_Area), by = c("ReplicateName", "MIS")) %>%
    left_join(Wei.IS.means %>% 
                rename(MIS = MassFeature), by = "MIS") %>%
    mutate(Adjusted_Area = AreaValue/IS_Area*Average.Area)
}


Wei.area.norm <- do.call(rbind, Split_Dat) %>% 
  select(-IS_Area, -Average.Area) 


# Standardize name structure to: Date_type_ID_replicate_anythingextraOK) ----------------------------------------------------------------
Wei.mydata_new <- Wei.area.norm %>% 
  separate(ReplicateName, c("runDate", "type", "SampID","replicate"), "_") %>%
  mutate(Run.Cmpd = paste(Wei.area.norm$ReplicateName, Wei.area.norm$MassFeature))


# Find the B-MIS for each MassFeature----------------------------------------------------------------

# Look only at the Pooled samples, to get a lowest RSD of the pooled possible (RSD_ofPoo), 
# then choose which IS reduces the RSD the most (Poo.Picked.IS) 
Wei.poodat <- Wei.mydata_new %>%
  filter(type == "Poo") %>%
  group_by(SampID, MassFeature, MIS) %>%
  summarise(RSD_ofPoo_IND = sd(Adjusted_Area, na.rm = TRUE) / mean(Adjusted_Area, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo_IND = ifelse(RSD_ofPoo_IND == "NaN", NA, RSD_ofPoo_IND)) %>%
  group_by(MassFeature, MIS) %>%
  summarise(RSD_ofPoo =  mean(RSD_ofPoo_IND, na.rm = TRUE))

Wei.poodat <- Wei.poodat %>% 
  left_join(Wei.poodat %>% group_by(MassFeature) %>%
              summarise(Poo.Picked.IS = unique(MIS)[which.min(RSD_ofPoo)] [1]))


# Get the original RSD, calculate RSD change, decide if MIS is acceptable----------------------------------------------------------------
Wei.poodat <- left_join(Wei.poodat, Wei.poodat %>%
                      filter(MIS == "Inj_vol" ) %>%
                      mutate(Orig_RSD = RSD_ofPoo) %>%
                      select(-RSD_ofPoo, -MIS)) %>%
  mutate(del_RSD = (Orig_RSD - RSD_ofPoo)) %>%
  mutate(percentChange = del_RSD/Orig_RSD) %>%
  mutate(accept_MIS = (percentChange > cut.off & Orig_RSD > cut.off2)) 


# Change the BMIS to "Inj_vol" if the BMIS is not an acceptable----------------------------------------------------------------

# Adds a column that has the BMIS, not just Poo.picked.IS
# Changes the FinalBMIS to inject_volume if its no good

Wei.fixedpoodat <- Wei.poodat %>%
  #filter(MIS == "Poo.Picked.IS") %>% # original from krh
  #filter(MIS == Poo.Picked.IS) %>%
  mutate(FinalBMIS = ifelse(accept_MIS == "FALSE", "Inj_vol", Poo.Picked.IS)) %>%
  mutate(FinalRSD = RSD_ofPoo) 

Wei.newpoodat <- Wei.poodat %>% 
  left_join(Wei.fixedpoodat %>% select(MassFeature, FinalBMIS)) %>%
  filter(MIS == FinalBMIS) %>%
  mutate(FinalRSD = RSD_ofPoo)

Try <- Wei.newpoodat %>% 
  filter(FinalBMIS != "Inj_vol")

QuickReport <- print(paste("% of MFs that picked a BMIS", 
                           length(Try$MassFeature) / length(Wei.newpoodat$MassFeature), 
                           "RSD improvement cutoff", cut.off,
                           "RSD minimum cutoff", cut.off2,
                           sep = " "))


# Evaluate the results of your BMIS cutoff----------------------------------------------------------------
IS_toISdat <- Wei.mydata_new %>%
  filter(MassFeature %in% Wei.IS.data$MassFeature) %>%
  select(MassFeature, MIS, Adjusted_Area, type) %>%
  filter(type == "Smp") %>%
  group_by(MassFeature, MIS) %>%
  summarise(RSD_ofSmp = sd(Adjusted_Area)/mean(Adjusted_Area)) %>%
  left_join(Wei.poodat %>% select(MassFeature, MIS, RSD_ofPoo, accept_MIS))

injectONlY_toPlot <- IS_toISdat %>%
  filter(MIS == "Inj_vol") 


ISTest_plot <- ggplot() +
  geom_point(dat = IS_toISdat, shape = 21, color = "black", size = 2,aes(x = RSD_ofPoo, y = RSD_ofSmp, fill = accept_MIS)) + 
  scale_fill_manual(values=c("white","dark gray")) +
  geom_point(dat = injectONlY_toPlot, aes(x = RSD_ofPoo, y = RSD_ofSmp), size = 3) +
  facet_wrap(~ MassFeature)
#print(ISTest_plot)


# Return data that is normalized via BMIS----------------------------------------------------------------


## original
Wei.BMIS_normalizedData <- Wei.newpoodat %>% select(MassFeature, FinalBMIS, Orig_RSD, FinalRSD) %>%
  left_join(Wei.mydata_new %>% rename(FinalBMIS = MIS)) %>%
  unique() 
  #filter(!MassFeature %in% Wei.IS.data$MassFeature)
##


write.csv(Wei.BMIS_normalizedData, file = "~/Downloads/Wei_Transect_BMISd.csv")




