library(data.table)
library(ggplot2)
library(parallel)
library(stringr)
library(tidyverse)
options(scipen=999)

## This BMIS is for Wei's eddy center data, with fixed valine!
## The file comes from Skyline.

# Imports -----------------------------------------------------------------

Wei.eddycenter.SampKey_all <- read.csv("data/Wei_EddyCenter_HILICPosNeg_QC_FixValine.csv", skip = 1) %>%
  select(Replicate.Name) %>%
  rename(ReplicateName = Replicate.Name) %>%
  mutate(Bio.Normalization = ifelse(str_detect(ReplicateName, "Half"), 0.5, 1.0)) %>%
  mutate(ReplicateName = substring(ReplicateName, 8)) %>%
  unique()
         
         
Wei.Internal.Standards <- read.csv("data/Ingalls_Lab_Standards.csv") %>%
  filter(Column == "HILIC") %>%
  filter(Compound.Type == "Internal Standard") 

trimws(Wei.Internal.Standards$Compound.Name, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")


# HILICPos and HILICNeg data
Wei.eddycenter <- read.csv("data/Wei_EddyCenter_HILICPosNeg_QC_FixValine.csv", skip = 1) %>%
  mutate(Replicate.Name = substring(Replicate.Name, 8)) %>%
  rename(ReplicateName = Replicate.Name) %>%
  mutate(Metabolite.name = Mass.Feature) %>%
  select(-Mass.Feature)

# Change class + adjust data. Set cutoff values -----------------------------------------------------------------
Wei.eddycenter <- Wei.eddycenter %>%
  filter(!str_detect(ReplicateName, "Blk")) %>%
  filter(!str_detect(ReplicateName, "Std")) %>%
  mutate(ReplicateName = as.character(ReplicateName)) %>%
  mutate(Metabolite.name = as.character(Metabolite.name)) %>%
  mutate(Retention.Time = as.numeric(Retention.Time)) %>%
  mutate(Area.with.QC = as.numeric(Area.with.QC)) %>% 
  select(Metabolite.name, ReplicateName, Area, Area.with.QC, Column) %>%
  mutate(Metabolite.name = ifelse(str_detect(Metabolite.name, "Ingalls_"), sapply(strsplit(Metabolite.name, "_"), `[`, 2), Metabolite.name))

cut.off <- 0.3 # % decrease in RSD of pooled injections, aka improvement cutoff
cut.off2 <- 0.1 # RSD minimum

# Match eddycenter data with Internal Standards list -----------------------------------------------------------------
Wei.eddycenter.withIS <- Wei.eddycenter %>%
  filter(Metabolite.name %in% Wei.Internal.Standards$Compound.Name) 

Wei.eddycenter.NoIS <- Wei.eddycenter %>%
  filter(!Metabolite.name %in% Wei.Internal.Standards$Compound.Name) 


# Read in Internal Standard data -----------------------------------------------------------------
# If injection volume is known, add in here.
Wei.eddycenter.IS.data <- Wei.eddycenter.withIS %>%
  select(ReplicateName, Metabolite.name, Area.with.QC) %>%
  mutate(MassFeature = Metabolite.name) %>%
  select(-Metabolite.name) 

Wei.eddycenter.SampKey <- Wei.eddycenter.SampKey_all %>%
  filter(ReplicateName %in% Wei.eddycenter.IS.data$ReplicateName) %>% 
  select(ReplicateName, Bio.Normalization) %>%
  mutate(MassFeature = "Inj_vol",
         Area.with.QC = Bio.Normalization) %>%
  select(ReplicateName, Area.with.QC, MassFeature)

Wei.eddycenter.IS.data <- rbind(Wei.eddycenter.IS.data, Wei.eddycenter.SampKey) 


# Extraction replication of Internal Standards -----------------------------------------------------------------
IS_inspectPlot <- ggplot(Wei.eddycenter.IS.data, aes(x = ReplicateName, y = Area.with.QC)) + 
  geom_bar(stat = "identity") + 
  facet_wrap( ~MassFeature, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 5), 
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10))+
  ggtitle("IS Raw Areas")
print(IS_inspectPlot)


# Edit data so names match-----------------------------------------------------------------
Wei.eddycenter.IS.data <- Wei.eddycenter.IS.data %>% 
  mutate(ReplicateName = ReplicateName %>%
           str_replace("-",".")) %>%
  arrange(ReplicateName)

Wei.eddycenter.long  <- Wei.eddycenter.NoIS %>%
  rename(MassFeature = Metabolite.name) %>%
  select(ReplicateName, MassFeature, Area.with.QC) %>%
  arrange(ReplicateName)


test_IS.data <- unique(Wei.eddycenter.IS.data$ReplicateName)
test_long.data <- unique(Wei.eddycenter.long$ReplicateName)
all.equal(test_IS.data, test_long.data)


# Caluclate mean values for each IS----------------------------------------------------------------
Wei.eddycenter.IS.means <- Wei.eddycenter.IS.data %>% 
  filter(!grepl("_Blk_", ReplicateName)) %>%
  mutate(MassFeature = as.factor(MassFeature)) %>%
  group_by(MassFeature) %>%
  summarise(Average.Area = mean(as.numeric(Area.with.QC), na.rm = TRUE)) %>%
  mutate(MassFeature = as.character(MassFeature))

Wei.eddycenter.IS.means[is.na(Wei.eddycenter.IS.means)] <- NA


# Normalize to each internal Standard----------------------------------------------------------------
Wei.eddycenter.binded <- rbind(Wei.eddycenter.IS.data, Wei.eddycenter.long) %>%
  arrange(MassFeature)

Split_Dat <- list()

for (i in 1:length(unique(Wei.eddycenter.IS.data$MassFeature))) {
  Split_Dat[[i]] <- Wei.eddycenter.binded %>% 
    mutate(MIS = unique(Wei.eddycenter.IS.data$MassFeature)[i]) %>%
    left_join(Wei.eddycenter.IS.data %>% 
                rename(MIS = MassFeature, IS_Area = Area.with.QC) %>% 
                select(MIS, ReplicateName, IS_Area), by = c("ReplicateName", "MIS")) %>%
    left_join(Wei.eddycenter.IS.means %>% 
                rename(MIS = MassFeature), by = "MIS") %>%
    mutate(Adjusted_Area = Area.with.QC/IS_Area*Average.Area)
}


Wei.eddycenter.area.norm <- do.call(rbind, Split_Dat) %>% 
  select(-IS_Area, -Average.Area) 


# Standardize name structure to: Date_type_ID_replicate_anythingextraOK) ----------------------------------------------------------------
Wei.eddycenter.mydata_new <- Wei.eddycenter.area.norm %>% 
  separate(ReplicateName, c("type", "SampID", "replicate"), "_") %>%
  mutate(Run.Cmpd = paste(Wei.eddycenter.area.norm$ReplicateName, Wei.eddycenter.area.norm$MassFeature))


# Find the B-MIS for each MassFeature----------------------------------------------------------------

# Look only at the Pooled samples, to get a lowest RSD of the pooled possible (RSD_ofPoo), 
# then choose which IS reduces the RSD the most (Poo.Picked.IS) 
Wei.eddycenter.poodat <- Wei.eddycenter.mydata_new %>%
  filter(type == "Poo") %>%
  group_by(SampID, MassFeature, MIS) %>%
  summarise(RSD_ofPoo_IND = sd(Adjusted_Area, na.rm = TRUE) / mean(Adjusted_Area, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo_IND = ifelse(RSD_ofPoo_IND == "NaN", NA, RSD_ofPoo_IND)) %>% 
  group_by(MassFeature, MIS) %>%
  summarise(RSD_ofPoo =  mean(RSD_ofPoo_IND, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo = ifelse(RSD_ofPoo == "NaN", NA, RSD_ofPoo)) # New addition to transform NaNs to NAs


Wei.eddycenter.poodat <- Wei.eddycenter.poodat %>% 
  left_join(Wei.eddycenter.poodat %>% group_by(MassFeature) %>%
              summarise(Poo.Picked.IS = unique(MIS)[which.min(RSD_ofPoo)] [1]))


# Get the original RSD, calculate RSD change, decide if MIS is acceptable----------------------------------------------------------------
Wei.eddycenter.poodat <- left_join(Wei.eddycenter.poodat, Wei.eddycenter.poodat %>%
                                   filter(MIS == "Inj_vol" ) %>%
                                   mutate(Orig_RSD = RSD_ofPoo) %>%
                                   select(-RSD_ofPoo, -MIS)) %>%
  mutate(del_RSD = (Orig_RSD - RSD_ofPoo)) %>%
  mutate(percentChange = del_RSD/Orig_RSD) %>%
  mutate(accept_MIS = (percentChange > cut.off & Orig_RSD > cut.off2)) 


# Change the BMIS to "Inj_vol" if the BMIS is not an acceptable----------------------------------------------------------------

# Adds a column that has the BMIS, not just Poo.picked.IS
# Changes the FinalBMIS to inject_volume if its no good

Wei.eddycenter.fixedpoodat <- Wei.eddycenter.poodat %>%
  filter(MIS == Poo.Picked.IS) %>% # If not commented out, this line filters out all but three compounds.
  mutate(FinalBMIS = ifelse(accept_MIS == "FALSE", "Inj_vol", Poo.Picked.IS)) %>%
  mutate(FinalRSD = RSD_ofPoo) 

Wei.newpoodat <- Wei.eddycenter.poodat %>% 
  left_join(Wei.eddycenter.fixedpoodat %>% select(MassFeature, FinalBMIS)) %>%
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
IS_toISdat <- Wei.eddycenter.mydata_new %>%
  filter(MassFeature %in% Wei.eddycenter.IS.data$MassFeature) %>%
  select(MassFeature, MIS, Adjusted_Area, type) %>%
  filter(type == "Smp") %>%
  group_by(MassFeature, MIS) %>%
  summarise(RSD_ofSmp = sd(Adjusted_Area)/mean(Adjusted_Area)) %>%
  left_join(Wei.eddycenter.poodat %>% select(MassFeature, MIS, RSD_ofPoo, accept_MIS))

injectONlY_toPlot <- IS_toISdat %>%
  filter(MIS == "Inj_vol") 


ISTest_plot <- ggplot() +
  geom_point(dat = IS_toISdat, shape = 21, color = "black", size = 2,aes(x = RSD_ofPoo, y = RSD_ofSmp, fill = accept_MIS)) + 
  scale_fill_manual(values=c("white","dark gray")) +
  geom_point(dat = injectONlY_toPlot, aes(x = RSD_ofPoo, y = RSD_ofSmp), size = 3) +
  facet_wrap(~ MassFeature)
print(ISTest_plot)


# Return data that is normalized via BMIS----------------------------------------------------------------


## original
Wei.eddycenter.BMIS_normalizedData <- Wei.newpoodat %>% select(MassFeature, FinalBMIS, Orig_RSD, FinalRSD) %>%
  left_join(Wei.eddycenter.mydata_new, by = "MassFeature") %>%
  filter(MIS == FinalBMIS) %>%
  unique() 

write.csv(Wei.eddycenter.BMIS_normalizedData, file = "~/Downloads/Wei_Eddycenter_BMISd_withQC_fixValine.csv")




