
# 
# system("xcode-select --install")
# if(!requireNamespace("Rcpp", quietly = TRUE))
#   install.packages("Rcpp")
# install.packages("/Applications/Wind API.app/Contents/R/WindR_1.0.8.tar.gz", repos=NULL)


rm(list = ls(all = TRUE))


# "rfPermute",
libraries = c("purrr","stringr","haven","sandwich","lmtest","nnet","rms","grf","varImp", "sandwich", "lmtest", "Hmisc", "ggplot2",
              "readxl", "plm", "psych", "sjmisc","RODBC","stringr",
              "iml","data.table", "igraph","timeDate", "stringr", "graphics","magick", "scales", "tidyr", "zoo", "foreach", "doParallel",
              "xgboost","shapr","randomForest", "rpart", "quantreg", "readxl","dplyr", "xlsx", "psych","qgraph", "gganimate","av",
              "gifski", "strex","matrixStats","tools","Hmisc","vars","aTSA","quantreg","rapport","sjmisc","haven","foreign")

lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)


wdir = "/Users/ruting/Documents/macbook/PcBack/23.Topic_wkh_cooperation/FRM"
setwd(wdir)
dir.create("Input/StockPrice")


Stockinfor <- read_excel(paste0(wdir,'/Input/StockInfor_2025.xlsx')) %>% data.frame()
Stockinfor = Stockinfor[c(1:(nrow(Stockinfor)-2)),]
colnames(Stockinfor) = c('ID','Shortname_CN','ShortName_EN','ListDate','ListLocation','Name_EN',
                         'Ownership','Prefeture','Prefeture_Code','Province','City','NEIC_B','NEIC_M','NEIC_S',
                         'NEIC_B_Code','NEIC_M_Code','NEIC_S_Code')

# 
# library(WindR)
# ret = w.start()
# print(ret)
# ret = w.isconnected()
# print(ret)

# read emission 
CarbonEmission <- read_excel(paste0(wdir,'/Input/CarbonData/上市公司碳排放/CNE_CemissDerive.xlsx')) %>% data.frame()

CarbonEmission = CarbonEmission[-c(1,2),]
CarbonEmission_keep = CarbonEmission[,c('EndDate','Symbol','CEmission','CemissionIntensity','NEIndustryCode','NEIndustryName')]
CarbonEmission_keep[,c('CEmission','CemissionIntensity')] = lapply(CarbonEmission_keep[,c('CEmission','CemissionIntensity')], as.numeric)
CarbonEmission_keep$Year = as.numeric(substr(CarbonEmission_keep$EndDate,1,4))
Carbon_Rank = data.frame(Symbol = unique(CarbonEmission_keep$Symbol), ID = NA, CarbonIntensity_Mean = NA, 
                         CarbonEmi_Mean = NA, Start = NA, End = NA)

for (iR in c(1: length(Carbon_Rank$ID))){
  iCode = Carbon_Rank$Symbol[iR]
  Temp = CarbonEmission_keep[CarbonEmission_keep$Symbol == iCode, c('EndDate','Symbol','Year','CEmission','CemissionIntensity')]
  if (as.numeric(substr(iCode,1,1) == 6)){
    ID = paste0(iCode,'.SH')
  }else if (substr(iCode,1,1) == '0' | substr(iCode,1,1) == '3'){
    ID = paste0(iCode,'.SZ')
  }else if (substr(iCode,1,1) == '8'){
    ID = paste0(iCode,'.BJ')
  }
  Carbon_Rank$ID[iR] = ID
  Carbon_Rank[iR, c('CarbonIntensity_Mean','CarbonEmi_Mean')] = colMeans(Temp[,c('CemissionIntensity','CEmission')],na.rm = TRUE)
  Carbon_Rank[iR,c('Start','End')] = c(min(Temp$Year),max(Temp$Year))
}

Carbon_Rank = Carbon_Rank[Carbon_Rank$Start<= 2010 & Carbon_Rank$End >= 2021 & (!is.na(Carbon_Rank$CarbonIntensity_Mean)),]
Unique_Indu = unique(CarbonEmission_keep[,c('Symbol','NEIndustryCode','NEIndustryName')]) %>%data.frame()

duplicated_symbols <- duplicated(Unique_Indu$Symbol)

Unique_Indu<- Unique_Indu[!duplicated_symbols, ]

# Top_Carbon Low_Carbon
Carbon_Rank = merge(Carbon_Rank,Unique_Indu, all.x = TRUE)

saveRDS(Carbon_Rank, file = "Output/Carbon_Rank.rds") 


quantiles <- quantile(Carbon_Rank$CarbonIntensity_Mean, probs = c(0.3, 0.7))

intensity_groups <- list(
  Low_Intensity = Carbon_Rank$ID[Carbon_Rank$CarbonIntensity_Mean < quantiles[1]],
  Medium_Intensity = Carbon_Rank$ID[Carbon_Rank$CarbonIntensity_Mean >= quantiles[1] 
                                    & Carbon_Rank$CarbonIntensity_Mean <= quantiles[2]],
  High_Intensity = Carbon_Rank$ID[Carbon_Rank$CarbonIntensity_Mean > quantiles[2]]
)

intensity_groups_Indus <- list(
  Low_Intensity = Carbon_Rank[Carbon_Rank$CarbonIntensity_Mean < quantiles[1], c('NEIndustryCode','NEIndustryName')],
  Medium_Intensity = Carbon_Rank[Carbon_Rank$CarbonIntensity_Mean >= quantiles[1] 
                                 & Carbon_Rank$CarbonIntensity_Mean <= quantiles[2],c('NEIndustryCode','NEIndustryName')],
  High_Intensity = Carbon_Rank[Carbon_Rank$CarbonIntensity_Mean > quantiles[2],c('NEIndustryCode','NEIndustryName')]
)


High_Intensity_Indus <- intensity_groups_Indus[[3]] %>%
  group_by(NEIndustryCode, NEIndustryName) %>%
  summarise(num = n(), .groups = 'drop') %>%
  arrange(desc(num))

Low_Intensity_Indus <- intensity_groups_Indus[[1]] %>%
  group_by(NEIndustryCode, NEIndustryName) %>%
  summarise(num = n(), .groups = 'drop') %>%
  arrange(desc(num))

Median_Intensity_Indus <- intensity_groups_Indus[[2]] %>%
  group_by(NEIndustryCode, NEIndustryName) %>%
  summarise(num = n(), .groups = 'drop') %>%
  arrange(desc(num))

write.csv(High_Intensity_Indus, 
          file = paste0('Output/HighCarbonIndustry.csv'), 
          row.names = FALSE,
          fileEncoding = "GB18030")

write.csv(Low_Intensity_Indus, file = paste0('Output/LowCarbonIndustry.csv'),
          row.names = FALSE,
          fileEncoding = "GB18030")

write.csv(Median_Intensity_Indus, file = paste0('Output/MedianCarbonIndustry.csv'), 
          row.names = FALSE,
          fileEncoding = "GB18030")


# Low_CarbonIntensi_Indus =Carbon_Rank[Carbon_Rank$ID %in% intensity_groups$Low_Intensity,]
# High_CarbonIntensi_Indus = Carbon_Rank[Carbon_Rank$ID %in% intensity_groups$High_Intensity,]

#test WSD function

# Type = c("High_CarbonIntensity","Low_CarbonIntensity")
Type = c("Medium_Intensity")
for (iType in Type){
  if (iType == "High_CarbonIntensity"){
    Stocklist = intensity_groups$High_Intensity
    
  }else if(iType == "Low_CarbonIntensity"){
    Stocklist = intensity_groups$Low_Intensity
  }else{
    Stocklist = intensity_groups$Medium_Intensity
  }
  dir.create(paste0("Input/StockPrice_",iType))
  write.csv(Stocklist, file = paste0('Input/Stocklist_',iType,'.csv'), row.names = FALSE)
  # for (iCode in Stocklist){
  #   w_wsd_data<-w.wsd(iCode, "open,high,low,close,volume,amt,ev", "2005-01-01", "2025-03-05", "unit=1;Currency=CNY;PriceAdj=F")
  #   w_wsd_data = w_wsd_data$Data
  #   if (nrow(w_wsd_data) > 0){
  #     write.csv(w_wsd_data, file = paste0('Input/StockPrice_',iType,'/',iCode,'.csv'), row.names = FALSE)
  #     print(paste0(which(Stocklist == iCode),'/', length(Stocklist)))
  #   }
  #   
  # }
}

# download from wind by wind platform




