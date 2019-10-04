
#Define all your inputs here
# samp.key.file <- "Sample_key.csv"
# is.names.file <- "InternalStandardNames.csv"
# xcms.dat.pos.file <- "HILICPos_IntegrationsBigPeaksWTargeted2.csv"
# xcms.dat.neg.file <- "HILICNeg_IntegrationsBigPeaksWTargeted.csv"
# cut.off <- 0.0 
# cut.off2 <- 0.00 

#Things to return 
#IS_inspectPlot (plot to make sure there aren't any internal standards we should kick out)
#QuickReport (% that picked BMIS, with cut off values)
#ISTest_plot (plot to evaluate if you cut off is appropriate)
#BMIS_normalizedData (tibble with the info you actually want!)


BMIS <- function(samp.key.file = "MetaData/Sample_key.csv",
                 is.names.file = "MetaData/InternalStandardNames.csv",
                 xcms.dat.pos.file = "RawOutput/HILICPos_IntegrationsBigPeaksWTargeted2.csv",
                 xcms.dat.neg.file = "RawOutput/HILICNeg_IntegrationsBigPeaksWTargeted.csv",
                 cut.off = 0.0,
                 cut.off2 = 0.00 ) {
  
  library(ggplot2)
  library(tidyverse)
  
#Import data - set filenames within this chunk for xcms output, sample key, and ISdata
SampKey_all <- read_csv(samp.key.file) 
IS_names <- read_csv(is.names.file)
xcms.dat_pos <- read_csv(xcms.dat.pos.file) %>% 
  mutate(Column = "HILICPos") %>% select(-Protein)
xcms.dat_neg <- read_csv(xcms.dat.neg.file) %>% 
  mutate(Column = "HILICNeg") 
xcms.dat <- rbind(xcms.dat_pos, xcms.dat_neg) %>%
    filter(!str_detect(`Replicate Name`, "Blk")) %>%
    filter(!str_detect(`Replicate Name`, "Std")) %>%
    mutate(Area = as.numeric(Area))
ISdatfull <- xcms.dat %>%
  filter(`Precursor Ion Name` %in% IS_names$Internal_Standards)
xcms.dat <- xcms.dat %>%
  filter(!`Precursor Ion Name` %in% IS_names$Internal_Standards)


#Read in Internal Standard data, add in injec_volume data from Sample Key
IS.dat <- ISdatfull %>%
  select(`Replicate Name`, `Precursor Ion Name`, Area) %>%
  mutate(MassFeature = `Precursor Ion Name`) %>%
  select(-`Precursor Ion Name`)

SampKey <- SampKey_all %>%
  filter(Sample.Name %in% IS.dat$`Replicate Name`) %>%
  select(Sample.Name, Injec_vol) %>%
  filter(!is.na(Injec_vol))%>%
  mutate(MassFeature = "Inj_vol",
         Area = Injec_vol,
         `Replicate Name` = Sample.Name) %>%
  select(`Replicate Name`, Area, MassFeature)

IS.dat <- rbind(IS.dat, SampKey) %>% mutate(Column = "HILICPos")


#Look at extraction replication of the Internal Standards----
IS_inspectPlot <- ggplot(IS.dat, aes(x=`Replicate Name`, y=Area)) + 
  geom_bar(stat="identity") + 
  facet_wrap( ~MassFeature, scales="free_y")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 5), 
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10))+
  ggtitle("IS Raw Areas")

#Edit data so names match
IS.dat <- IS.dat %>% mutate(Replicate.Name = `Replicate Name` %>%
                              str_replace("-","."))  %>%
  select(Area, Replicate.Name, MassFeature, Column)
xcms.long <- xcms.dat %>%
  rename(Replicate.Name = `Replicate Name`,
         MassFeature = `Precursor Ion Name`) %>%
  select(Replicate.Name, MassFeature,  Column, Area)
xcms.long <- xcms.long %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("170410_Poo_April11AqExtractsFull_",
                       "170410_Poo_April11AqExtracts_Full") %>%
           str_replace("170410_Poo_April11AqExtractsHalf_",
                       "170410_Poo_April11AqExtracts_Half")) 

IS.dat <- IS.dat %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("\\.","-") %>%
           str_replace("170410_Poo_April11AqExtractsFull_",
                       "170410_Poo_April11AqExtracts_Full") %>%
           str_replace("170410_Poo_April11AqExtractsHalf_",
                       "170410_Poo_April11AqExtracts_Half")) 

#Calculate mean values for each IS----
IS.means <- IS.dat %>% filter(!grepl("_Blk_", Replicate.Name)) %>%
  mutate(MassFeature = as.factor(MassFeature))%>%
  group_by(MassFeature) %>%
  summarise(ave = mean(as.numeric(Area)))



#Normalize to each internal Standard----
binded <- rbind(IS.dat, xcms.long)
Split_Dat <- list()
for (i in 1:length(unique(IS.dat$MassFeature))){
  Split_Dat[[i]] <- binded %>% mutate(MIS = unique(IS.dat$MassFeature)[i]) %>%
    left_join(IS.dat %>% 
                rename(MIS = MassFeature, IS_Area = Area) %>% 
                select(MIS, Replicate.Name, IS_Area)) %>%
    left_join(IS.means %>% 
                rename(MIS = MassFeature)) %>%
    mutate(Adjusted_Area = Area/IS_Area*ave)
}
area.norm <- do.call(rbind, Split_Dat) %>% select(-IS_Area, -ave)
  
  
#Break Up the Names (Name structure must be:  Date_type_ID_replicate_anythingextraOK)----
mydata_new <- area.norm %>% separate(Replicate.Name, 
                                     c("runDate",
                                       "type","SampID","replicate"),"_") %>%
  mutate(Run.Cmpd = paste(area.norm$Replicate.Name,area.norm$MassFeature))
  
  
  #Find the B-MIS for each MassFeature----
  #Look only the Pooled samples, to get a lowest RSD of the pooled possible (RSD_ofPoo), 
  #then choose which IS reduces the RSD the most (Poo.Picked.IS) 
poodat <- mydata_new %>%
  filter(type == "Poo")%>%
  group_by(SampID, MassFeature, MIS) %>%
  summarise(RSD_ofPoo_IND = sd(Adjusted_Area, 
                               na.rm = TRUE)/mean(Adjusted_Area, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo_IND = ifelse(RSD_ofPoo_IND == "NaN", NA, RSD_ofPoo_IND)) %>%
  group_by(MassFeature, MIS) %>%
  summarise(RSD_ofPoo =  mean(RSD_ofPoo_IND, na.rm = TRUE))

poodat <- poodat %>% left_join(poodat %>%
                                 group_by(MassFeature) %>%
                                 summarise(Poo.Picked.IS = 
                                             unique(MIS)[which.min(RSD_ofPoo)][1]))

#Get the starting point of the RSD (Orig_RSD), calculate the change in the RSD, say if the MIS is acceptable----
poodat <- left_join(poodat, poodat %>%
                      filter(MIS == "Inj_vol" ) %>%
                      mutate(Orig_RSD = RSD_ofPoo) %>%
                      select(-RSD_ofPoo, -MIS)) %>%
  mutate(del_RSD = (Orig_RSD - RSD_ofPoo)) %>%
  mutate(percentChange = del_RSD/Orig_RSD) %>%
  mutate(accept_MIS = (percentChange > cut.off & Orig_RSD > cut.off2)) 
  
  #Change the BMIS to "Inj_vol" if the BMIS is not an acceptable -----
  #Adds a column that has the BMIS, not just Poo.picked.IS
  #Changes the finalBMIS to inject_volume if its no good
fixedpoodat <- poodat %>%
  filter(MIS == Poo.Picked.IS)%>%
  mutate(FinalBMIS = ifelse((accept_MIS == "FALSE"), "Inj_vol", Poo.Picked.IS), 
         FinalRSD = RSD_ofPoo) 
newpoodat <- poodat %>% left_join(fixedpoodat %>% select(MassFeature, FinalBMIS)) %>%
  filter(MIS == FinalBMIS) %>%
  mutate(FinalRSD = RSD_ofPoo)
Try <- newpoodat %>% filter(FinalBMIS != "Inj_vol")
QuickReport <- paste("% of MFs that picked a BMIS", 
                       length(Try$MassFeature) / length(newpoodat$MassFeature), 
                       "RSD improvement cutoff", cut.off,
                       "RSD minimum cutoff", cut.off2,
                       sep = " ")
  
#Evaluate the results of your BMIS cutoff-----
  IS_toISdat <- mydata_new %>%
    filter(MassFeature %in% IS.dat$MassFeature) %>%
    select(MassFeature, MIS, Adjusted_Area, type) %>%
    filter(type == "Smp") %>%
    group_by(MassFeature, MIS) %>%
    summarise(RSD_ofSmp = sd(Adjusted_Area)/mean(Adjusted_Area)) %>%
    left_join(poodat %>% select(MassFeature, MIS, RSD_ofPoo, accept_MIS))
  
injectONlY_toPlot <- IS_toISdat %>%
    filter(MIS == "Inj_vol" ) 
  
  
ISTest_plot <- ggplot()+
    geom_point(dat = IS_toISdat, shape = 21, color = "black", size = 2,aes(x = RSD_ofPoo, y = RSD_ofSmp, fill = accept_MIS))+ 
    scale_fill_manual(values=c("white","dark gray"))+
    geom_point(dat = injectONlY_toPlot, aes(x = RSD_ofPoo, y = RSD_ofSmp), size = 3) +
    facet_wrap(~ MassFeature)
  
  #Get all the data back - and keep only the MF-MIS match set for the BMIS----
  #Add a column to the longdat that has important information from the FullDat_fixed, 
  #then only return data that is normalized via B-MIS normalization
BMIS_normalizedData <- newpoodat %>% select(MassFeature, FinalBMIS, Orig_RSD, FinalRSD) %>%
    left_join(mydata_new %>% rename(FinalBMIS = MIS)) %>% unique() %>%
    filter(!MassFeature %in% IS.dat$MassFeature)
  
  
  
BMISlist <- list(IS_inspectPlot, QuickReport, ISTest_plot, BMIS_normalizedData)
  
  return(invisible(BMISlist))
  
}
  
  
  