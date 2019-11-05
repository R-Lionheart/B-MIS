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

Wei.transect.SampKey_all <- read.csv("data/Sample.Key.EddyTransect.csv") %>%
  mutate(Sample.Name = Sample.Name %>%
           str_replace("-","."))
Wei.transect.SampKey_all[Wei.transect.SampKey_all == "180821_Poo_MesoScopeQC_1a"] <- "180821_Poo_MesoScopeQC_1"

Wei.Internal.Standards <- read.csv("data/Ingalls_Lab_Standards.csv") %>%
  filter(Column == "HILIC") %>%
  filter(Compound.Type == "Internal Standard")

trimws(Wei.Internal.Standards$Compound.Name, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")


# HILICPos and HILICNeg data
Wei.transect <- read.csv("data/Wei_Transect_HILICPosNeg_QC.csv", header = TRUE) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value))


# Change class + adjust data. Set cutoff values -----------------------------------------------------------------
Wei.transect <- Wei.transect %>%
  filter(!str_detect(ReplicateName, "Blk")) %>%
  filter(!str_detect(ReplicateName, "Std")) %>%
  mutate(ReplicateName = as.character(ReplicateName)) %>%
  mutate(Metabolite.name = as.character(Metabolite.name)) %>%
  mutate(RTValue = as.numeric(RTValue)) %>%
  mutate(AreaValue = as.numeric(AreaValue)) %>%
  mutate(SNValue = as.numeric(SNValue)) %>%
  mutate(Metabolite.name = ifelse(str_detect(Metabolite.name, "Ingalls_"), sapply(strsplit(Metabolite.name, "_"), `[`, 2), Metabolite.name))

cut.off <- 0.3 # 30% decrease in RSD of pooled injections, aka improvement cutoff
cut.off2 <- 0.1 # RSD minimum

# Match transect data with Internal Standards list -----------------------------------------------------------------
Wei.transect.withIS <- Wei.transect %>%
  filter(Metabolite.name %in% Wei.Internal.Standards$Compound.Name)

Wei.transect.NoIS <- Wei.transect %>%
  filter(!Metabolite.name %in% Wei.Internal.Standards$Compound.Name)


# Read in Internal Standard data -----------------------------------------------------------------
# If injection volume is known, add in here.
Wei.transect.IS.data <- Wei.transect.withIS %>%
  select(ReplicateName, Metabolite.name, Area.with.QC) %>%
  mutate(MassFeature = Metabolite.name) %>%
  select(-Metabolite.name) %>%
  filter(!MassFeature == "Guanosine Monophosphate, 15N5")

Wei.transect.IS.data$ReplicateName <- gsub("^.{0,1}", "", Wei.transect.IS.data$ReplicateName)

Wei.transect.SampKey <- Wei.transect.SampKey_all %>%
  filter(Sample.Name %in% Wei.transect.IS.data$ReplicateName) %>% # Drops standards from SampKey_all
  select(Sample.Name, Bio.Normalization) %>%
  mutate(MassFeature = "Inj_vol",
         Area.with.QC = Bio.Normalization,
         ReplicateName = Sample.Name) %>%
  select(ReplicateName, Area.with.QC, MassFeature)

Wei.transect.IS.data <- rbind(Wei.transect.IS.data, Wei.transect.SampKey)

IS_Issues <- Wei.transect.IS.data[is.na(Wei.transect.IS.data$Area.with.QC),]


# Extraction replication of Internal Standards -----------------------------------------------------------------
IS_inspectPlot <- ggplot(Wei.transect.IS.data, aes(x = ReplicateName, y = Area.with.QC)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~MassFeature, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 5),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10))+
  ggtitle("IS Raw Areas")
print(IS_inspectPlot)


# Edit data so names match-----------------------------------------------------------------
Wei.transect.IS.data <- Wei.transect.IS.data %>%
  mutate(ReplicateName = ReplicateName %>%
           str_replace("-",".")) %>%
  arrange(ReplicateName)

Wei.transect.long  <- Wei.transect.NoIS %>%
  rename(MassFeature = Metabolite.name) %>%
  select(ReplicateName, MassFeature, Area.with.QC) %>%
  arrange(ReplicateName)
Wei.transect.long$ReplicateName <- gsub("^.{0,1}", "", Wei.transect.long$ReplicateName)

# Test that names are equal
test_IS.data <- as.data.frame(unique(Wei.transect.IS.data$ReplicateName))
test_IS.data <- test_IS.data %>% mutate(ReplicateName = as.character(unique(Wei.transect.IS.data$ReplicateName)))
test_long.data <- as.data.frame(unique(Wei.transect.long$ReplicateName))
test_long.data <- test_long.data %>% mutate(ReplicateName = as.character(unique(Wei.transect.long$ReplicateName)))

identical(test_IS.data$ReplicateName, test_long.data$ReplicateName)

# Caluclate mean values for each IS----------------------------------------------------------------
Wei.transect.IS.means <- Wei.transect.IS.data %>%
  filter(!grepl("_Blk_", ReplicateName)) %>%
  mutate(MassFeature = as.factor(MassFeature)) %>%
  group_by(MassFeature) %>%
  summarise(Average.Area = mean(as.numeric(Area.with.QC), na.rm = TRUE)) %>%
  mutate(MassFeature = as.character(MassFeature))

Wei.transect.IS.means[is.na(Wei.transect.IS.means)] <- NA



# Normalize to each internal Standard----------------------------------------------------------------
Wei.transect.binded <- rbind(Wei.transect.IS.data, Wei.transect.long) %>%
  arrange(MassFeature)

Split_Dat <- list()

for (i in 1:length(unique(Wei.transect.IS.data$MassFeature))) {
  Split_Dat[[i]] <- Wei.transect.binded %>%
    mutate(MIS = unique(Wei.transect.IS.data$MassFeature)[i]) %>%
    left_join(Wei.transect.IS.data %>%
                rename(MIS = MassFeature, IS_Area = Area.with.QC) %>%
                select(MIS, ReplicateName, IS_Area), by = c("ReplicateName", "MIS")) %>%
    left_join(Wei.transect.IS.means %>%
                rename(MIS = MassFeature), by = "MIS") %>%
    mutate(Adjusted_Area = Area.with.QC/IS_Area*Average.Area)
}


Wei.transect.area.norm <- do.call(rbind, Split_Dat) %>%
  select(-IS_Area, -Average.Area)
test_wei.is.data<- Wei.transect.IS.data[rowSums(is.na(Wei.transect.IS.data)) > 0,]

# Standardize name structure to: Date_type_ID_replicate_anythingextraOK) ----------------------------------------------------------------
Wei.transect.mydata_new <- Wei.transect.area.norm %>%
  separate(ReplicateName, c("runDate", "type", "SampID", "replicate"), "_") %>%
  mutate(Run.Cmpd = paste(Wei.transect.area.norm$ReplicateName, Wei.transect.area.norm$MassFeature))


# Find the B-MIS for each MassFeature----------------------------------------------------------------

# Look only at the Pooled samples, to get a lowest RSD of the pooled possible (RSD_ofPoo),
# then choose which IS reduces the RSD the most (Poo.Picked.IS)
Wei.transect.poodat <- Wei.transect.mydata_new %>%
  filter(type == "Poo") %>%
  group_by(SampID, MassFeature, MIS) %>%
  summarise(RSD_ofPoo_IND = sd(Adjusted_Area, na.rm = TRUE) / mean(Adjusted_Area, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo_IND = ifelse(RSD_ofPoo_IND == "NaN", NA, RSD_ofPoo_IND)) %>%
  group_by(MassFeature, MIS) %>%
  summarise(RSD_ofPoo =  mean(RSD_ofPoo_IND, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo = ifelse(RSD_ofPoo == "NaN", NA, RSD_ofPoo)) # New addition to transform NaNs to NAs


Wei.transect.poodat <- Wei.transect.poodat %>%
  left_join(Wei.transect.poodat %>% group_by(MassFeature) %>%
              summarise(Poo.Picked.IS = unique(MIS)[which.min(RSD_ofPoo)] [1]))


# Get the original RSD, calculate RSD change, decide if MIS is acceptable----------------------------------------------------------------
Wei.transect.poodat <- left_join(Wei.transect.poodat, Wei.transect.poodat %>%
                                   filter(MIS == "Inj_vol" ) %>%
                                   mutate(Orig_RSD = RSD_ofPoo) %>%
                                   select(-RSD_ofPoo, -MIS)) %>%
  mutate(del_RSD = (Orig_RSD - RSD_ofPoo)) %>%
  mutate(percentChange = del_RSD/Orig_RSD) %>%
  mutate(accept_MIS = (percentChange > cut.off & Orig_RSD > cut.off2))


# Change the BMIS to "Inj_vol" if the BMIS is not an acceptable----------------------------------------------------------------

# Adds a column that has the BMIS, not just Poo.picked.IS
# Changes the FinalBMIS to inject_volume if its no good

Wei.transect.fixedpoodat <- Wei.transect.poodat %>%
  filter(MIS == Poo.Picked.IS) %>% # If not commented out, this line filters out all but three compounds.
  mutate(FinalBMIS = ifelse(accept_MIS == "FALSE", "Inj_vol", Poo.Picked.IS)) %>%
  mutate(FinalRSD = RSD_ofPoo)

Wei.newpoodat <- Wei.transect.poodat %>%
  left_join(Wei.transect.fixedpoodat %>% select(MassFeature, FinalBMIS)) %>%
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
IS_toISdat <- Wei.transect.mydata_new %>%
  filter(MassFeature %in% Wei.transect.IS.data$MassFeature) %>%
  select(MassFeature, MIS, Adjusted_Area, type) %>%
  filter(type == "Smp") %>%
  group_by(MassFeature, MIS) %>%
  summarise(RSD_ofSmp = sd(Adjusted_Area, na.rm = TRUE)/mean(Adjusted_Area, na.rm = TRUE)) %>%
  left_join(Wei.transect.poodat %>% select(MassFeature, MIS, RSD_ofPoo, accept_MIS))

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
Wei.transect.BMIS_normalizedData <- Wei.newpoodat %>% select(MassFeature, FinalBMIS, Orig_RSD, FinalRSD) %>%
  left_join(Wei.transect.mydata_new, by = "MassFeature") %>%
  filter(MIS == FinalBMIS) %>%
  unique()

#write.csv(Wei.transect.BMIS_normalizedData, file = "~/Downloads/Wei_Transect_BMISd_withQCNov4.csv")


