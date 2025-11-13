# 0211: overlapping problem 
# 1013: double sort modification
# Positive Jump of FRM > premium return, VaR of premium?
  
  
rm(list = ls(all = TRUE))

libraries = c("dynlm","readxl","SHAPforxgboost","shapr","xgboost","readxl","iml","ggplot2", "data.table", "igraph","timeDate", "stringr", "graphics","magick", "scales", "tidyr", "zoo","xts", "foreach", "doParallel",
              "xgboost","shapr","randomForest", "rpart", "quantreg", "readxl","dplyr", "xlsx", "psych","qgraph", "gganimate","av",
              "gifski", "strex","matrixStats","tools","Hmisc","vars","aTSA","quantreg","rapport","sjmisc","haven","foreign","e1071")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})

lapply(libraries, library, quietly = TRUE, character.only = TRUE)

library(dplyr)
library(haven)
wdir = "/Users/ruting/Documents/macbook/PcBack/23.Topic_wkh_cooperation/FRM"
save =  paste0(wdir, "/Output/EfficiencyTest")

setwd(wdir)
date_start_source = 20140704
date_end_source =  20250127
save_date = 20250318

# gen volatility of VAR and \varepsilon

channel = c("HighCarbonIntens", "LowCarbonIntens")
for (iCha in channel){
  
  output_path = paste0("Output/", iCha, "/FitQr")
  
  FitQr = readRDS(paste0(output_path,"/FitQr_",iCha,".rds"))
  N_h = length(FitQr)
  #Transform the list of lambdas into a wide dataset
  stock_names = vector()
  for (t in 1:N_h) stock_names = c(stock_names, attributes(FitQr[[t]])$dimnames[[2]])
  
  stock_names = unique(stock_names)
  N_names = length(stock_names)
  FitQr_wide = matrix(NA, N_h, N_names+1)
  
  for (k in 1:N_names) 
    for (t in 1:N_h) 
      if (stock_names[k] %in% attributes(FitQr[[t]])$dimnames[[2]]) 
        FitQr_wide[t, k+1] = FitQr[[t]][, stock_names[k]]
  
  FitQr_wide = round(FitQr_wide, digits = 6)
  FitQr_wide[, 1] = names(FitQr)
  colnames(FitQr_wide) = c("date", stock_names)
  
  # 先把矩阵转为 data.frame，确保列名存在
  FitQr_wide <- as.data.frame(FitQr_wide)
  
  # 提取数值部分（排除第一列日期）
  num_part <- FitQr_wide[, -1]
  
  # 统计每列中 NA 的比例
  zero_na_ratio <- colMeans( is.na(num_part))
  
  # 找出 0 比例 <= 0.2 的列
  keep_cols <- names(zero_na_ratio[zero_na_ratio < 0.1])
  FitQr_wide = FitQr_wide[,c("date",keep_cols)]
  
  # volatility
  all_std = FitQr_wide
  all_std[,-1] = NA
  for (i in c(63:nrow(all_std))){
    all_std[i,-1] =  apply(FitQr_wide[(i-62):i,-1], 2, sd,na.rm = TRUE)
  }

  write.csv(all_std, paste0(output_path,"/Volatility_VaR5.csv"), quote = FALSE)
}


# generate Sharpe ratio volatility 

# channel = c("HighCarbonIntens", "LowCarbonIntens")
channel = c("LowCarbonIntens")
for (iCha in channel){
  input_path = paste0("Input/", iCha, "/", date_start_source, "-", date_end_source)
  output_path = paste0("Output/", iCha, "/Lambda")
  
  # data
  fun_loc = function(x){
    loc = which(f$Stkcd == x)
    return(loc)
  }
  
  # volatility 63 days
  mktcap = read.csv(file = paste0(input_path, "/", iCha, "_Mktcap_", 
                                  date_end_source, ".csv"), header = TRUE) %>% data.frame()
  stock_prices = read.csv(file = paste0(input_path, "/", iCha, "_Price_", 
                                        date_end_source, ".csv"), header = TRUE)%>% data.frame()
  macro = read.csv(file = paste0(input_path, "/", iCha, "_Macro_", 
                                 date_end_source, ".csv"), header = TRUE)%>% data.frame()
  
  # crix = read.csv(file = paste0(input_path, "/vcrix.csv"), header = TRUE)%>% data.frame()
  # crix = crix[,-c(1,2)]
  
  # calculate weighted index
  all_prices = stock_prices
  all_prices[, -1] = sapply(all_prices[, -1], as.numeric)
  # write.csv(all_prices, paste0(wdir, "/MethodAdd/Price_",save_date,".csv"), quote = FALSE)
  
  all_return = all_prices[-1, ]
  all_return[,-1] = 0
  
  all_return[, -1] = 
    as.data.frame(diff(log(as.matrix(all_prices[, -1])))) 
  all_return[is.na(all_return)] = 0
  # write.csv(all_return, paste0(wdir, "/MethodAdd/Return_",save_date,".csv"), quote = FALSE)
  
  # volatility
  all_std = all_return
  all_std[,-1] = NA
  for (i in c(63:nrow(all_std))){
    all_std[i,-1] =  apply(all_return[(i-62):i,-1], 2, sd)
  }
  
  # Sharpe
  all_mean_63 = all_return
  all_mean_63[,-1] = NA
  for (i in c(63:nrow(all_mean_63))){
    all_mean_63[i,-1] =  apply(all_return[(i-62):i,-1], 2, mean)
  }
  
  all_Sharpe_63 = all_mean_63
  all_Sharpe_63[,-1] = all_mean_63[,-1] / all_std[,-1]
  all_Sharpe_63 = all_Sharpe_63[-c(1:62),]
  
  # Sharpe volatility
  all_Sharpe_vol_63 = all_Sharpe_63
  for (i in c(63:nrow(all_Sharpe_vol_63))){
    all_Sharpe_vol_63[i,-1] =  apply(all_Sharpe_63[(i-62):i,-1], 2, sd)
  }
  
  all_Sharpe_vol_63 = all_Sharpe_vol_63[-c(1:62),]
  
  all_std = all_std[all_std$Date>= all_Sharpe_vol_63$Date[1], ]
  all_Sharpe_63 = all_Sharpe_63[all_Sharpe_63$Date>= all_Sharpe_vol_63$Date[1], ]
  
  write.csv(all_std, paste0(save, "/",channel,"Volatility.csv"), quote = FALSE)
  write.csv(all_Sharpe_63, paste0(save, "/",channel,"/Sharpe.csv"), quote = FALSE)
  write.csv(all_Sharpe_vol_63, paste0(save, "/",channel,"/SharpeVolatility.csv"), quote = FALSE)
  
}

# benchmark volatility
all_Sharpe_vol_High = read.csv(file = paste0(save, "/HighCarbonIntens/SharpeVolatility.csv"), header = TRUE)%>% data.frame()
all_Sharpe_vol_High = all_Sharpe_vol_High[,-1]
all_Sharpe_vol_Low = read.csv(file = paste0(save, "/LowCarbonIntens/SharpeVolatility.csv"), header = TRUE)%>% data.frame()
all_Sharpe_vol_Low = all_Sharpe_vol_Low[,-1]


all_Vol_High = read.csv(file = paste0(save, "/HighCarbonIntens/Volatility.csv"), header = TRUE)%>% data.frame()
all_Vol_High = all_Vol_High[,-1]
all_Vol_Low = read.csv(file = paste0(save, "/LowCarbonIntens/Volatility.csv"), header = TRUE)%>% data.frame()
all_Vol_Low = all_Vol_Low[,-1]

mean_Sharpe_vol_Low = data.frame(Date = all_Sharpe_vol_Low$Date, mean_Sharpe_vol_Low = apply(all_Sharpe_vol_Low[,-1], 1, mean,na.rm = TRUE))
mean_vol_Low = data.frame(Date = all_Vol_Low$Date, mean_vola_Low = apply(all_Vol_Low[,-1], 1, mean, na.rm = TRUE))


Sample = colnames(all_Sharpe_vol_High)[-1]
Sample = substr(Sample,2,10)

# Financial reports
dir_FR = "Input/FinancialReport"

cIndi = c("Asset","Equity","Leverage","ROE","CapitalExpenditure")
Data_IF = list()
for(iIF in cIndi){
  IF = read_xlsx(paste0(dir_FR, "/",iIF,".xlsx"))%>% data.frame()
  IF = IF[c(1:(nrow(IF)-2)),-2]
  colnames(IF) = c("ID",c(2000:2024))
  Data_IF[[iIF]] <- IF
  
}

# add Market Value 
input_path = paste0("Input/HighCarbonIntens/", date_start_source, "-", date_end_source)
mktcap = read.csv(file = paste0(input_path, "/HighCarbonIntens_Mktcap_", 
                                date_end_source, ".csv"), header = TRUE) %>% data.frame()
unique_Year = unique(as.numeric(substr(mktcap$Date,1,4)))

mktcap_Year = mktcap
mktcap_Year = mktcap_Year[1:length(unique_Year),]
mktcap_Year$Date = unique_Year
mktcap_Year[,-1] = NA
for (iY in c(1:nrow(mktcap_Year))){
  Temp = mktcap[as.numeric(substr(mktcap$Date,1,4)) == unique_Year[iY],]
  mktcap_Year[iY,-1] = Temp[nrow(Temp),-1]
}
colnames(mktcap_Year)[1] = "Year"

for (i in c(2:ncol(all_Sharpe_vol_High))){
  Temp = cbind(all_Sharpe_vol_High[,c(1,i)],all_Vol_High[,i])
  colnames(Temp)[-1] = c("Sharpe_vol_High","Vol_High")
  Temp$Year = as.numeric(substr(Temp$Date,1,4))
  iID = Sample[i-1]
  Result_List <- list()
  for(iIF in cIndi){
    Temp_IF <- Data_IF[[iIF]]  # 获取当前指标的数据
    Temp_IF <- Temp_IF[Temp_IF$ID == iID, ]  # 筛选 ID
    Temp_IF <- data.frame(Year = 2000:2024, Value = as.numeric(Temp_IF[-1]))  # 创建 DataFrame
    colnames(Temp_IF)[2] <- iIF  # 重新命名列
    Temp_IF[,3] = NA
    Temp_IF[c(2:nrow(Temp_IF)),3] = Temp_IF[1:(nrow(Temp_IF)-1),2]
    colnames(Temp_IF)[3] = paste0(iIF,"_LastY")
    # 存入 list
    Result_List[[iIF]] <- Temp_IF
  }
  
  # 合并所有数据为一个 data.frame
  Final_Data <- Reduce(function(x, y) merge(x, y, by = "Year", all = TRUE), Result_List)
  
  # add MV
  Temp_MV = mktcap_Year[,c(1,i)]
  colnames(Temp_MV)[2] = "MV"
  Temp_MV$MV = Temp_MV$MV * 10^8
  Temp_MV[,3] = NA
  Temp_MV[c(2:nrow(Temp_MV)),3] = Temp_MV[1:(nrow(Temp_MV)-1),2]
  colnames(Temp_MV)[3] = "MV_LastY"
  
  Temp = merge(Temp, Final_Data, all.x = TRUE)
  Temp = merge(Temp, Temp_MV, all.x = TRUE)

  Temp$ID = iID
  
  if (i == 2){
    Data_Reg = Temp
  }else{
    Data_Reg = rbind(Data_Reg,Temp) 
  }
}

Data_Reg = merge(Data_Reg, mean_Sharpe_vol_Low, all.x = TRUE)
Data_Reg = merge(Data_Reg, mean_vol_Low, all.x = TRUE)

Data_Reg$Date = as.Date(Data_Reg$Date)

# add FRM 1% 5% 10%
channel = c("HighCarbonIntens", "LowCarbonIntens")
for (iCha in channel){
  output_path = paste0("Output/", iCha)
  FRM_5 = read.csv(file = paste0(output_path,"/Lambda/FRM_",iCha,"_index.csv"), header = TRUE) %>% data.frame()
  colnames(FRM_5)[2] = paste0("FRM_5_",iCha)
  FRM_1 = read.csv(file = paste0(output_path,"/Sensitivity/tau=1/s=63/Lambda/FRM_",iCha,"_index.csv"), header = TRUE) %>% data.frame()
  colnames(FRM_1)[2] = paste0("FRM_1_",iCha)
  FRM_10 = read.csv(file = paste0(output_path,"/Sensitivity/tau=10/s=63/Lambda/FRM_",iCha,"_index.csv"), header = TRUE) %>% data.frame()
  colnames(FRM_10)[2] = paste0("FRM_10_",iCha)
  
  FRM_out = cbind(FRM_5,FRM_1[,2],FRM_10[,2])
  colnames(FRM_out) = c("Date",paste0("FRM_5_",iCha),paste0("FRM_1_",iCha),paste0("FRM_10_",iCha))
  if (iCha == channel[1]){
    FRM_Reg = FRM_out
  }else{
    FRM_Reg = cbind(FRM_Reg, FRM_out[,-1])
  }
}  

FRM_Reg$FRM_5_High_Low = FRM_Reg$FRM_5_HighCarbonIntens / FRM_Reg$FRM_5_LowCarbonIntens
FRM_Reg$FRM_10_High_Low = FRM_Reg$FRM_10_HighCarbonIntens / FRM_Reg$FRM_10_LowCarbonIntens
FRM_Reg$FRM_1_High_Low = FRM_Reg$FRM_1_HighCarbonIntens / FRM_Reg$FRM_1_LowCarbonIntens
FRM_Reg$Date = as.Date(FRM_Reg$Date)
saveRDS(FRM_Reg, file = paste0(save, "/FRM_Reg.rds")) 

Data_Reg = merge(Data_Reg, FRM_Reg, all.x = TRUE)

write_dta(Data_Reg, paste0(save,"/InsampleBeta_CarbonRisk.dta"))
saveRDS(Data_Reg, file = paste0(save, "/InsampleBeta_CarbonRisk.rds")) 

# volatility of VaR
all_Vol_High = read.csv(file = paste0("Output/HighCarbonIntens/FitQr/Volatility_VaR5.csv"), header = TRUE)%>% data.frame()
all_Vol_High = all_Vol_High[,-1]
colnames(all_Vol_High)[1] = "Date"

all_Vol_Low = read.csv(file = paste0("Output/LowCarbonIntens/FitQr/Volatility_VaR5.csv"), header = TRUE)%>% data.frame()
all_Vol_Low = all_Vol_Low[,-1]
colnames(all_Vol_Low)[1] = "Date"

mean_vol_High = data.frame(Date = all_Vol_High$Date, mean_vola_High= apply(all_Vol_High[,-1], 1, mean, na.rm = TRUE))
mean_vol_Low = data.frame(Date = all_Vol_Low$Date, mean_vola_Low = apply(all_Vol_Low[,-1], 1, mean, na.rm = TRUE))

Sample = colnames(all_Vol_High)[-1]
Sample = substr(Sample,2,10)

# Financial reports
dir_FR = "Input/FinancialReport"
library(readxl)
cIndi = c("Asset","Equity","Leverage","ROE","CapitalExpenditure")
Data_IF = list()
for(iIF in cIndi){
  IF = read_xlsx(paste0(dir_FR, "/",iIF,".xlsx"))%>% data.frame()
  IF = IF[c(1:(nrow(IF)-2)),-2]
  colnames(IF) = c("ID",c(2000:2024))
  Data_IF[[iIF]] <- IF
  
}

# add Market Value 
input_path = paste0("Input/HighCarbonIntens/", date_start_source, "-", date_end_source)
mktcap = read.csv(file = paste0(input_path, "/HighCarbonIntens_Mktcap_", 
                                date_end_source, ".csv"), header = TRUE) %>% data.frame()
unique_Year = unique(as.numeric(substr(mktcap$Date,1,4)))

mktcap_Year = mktcap
mktcap_Year = mktcap_Year[1:length(unique_Year),]
mktcap_Year$Date = unique_Year
mktcap_Year[,-1] = NA
for (iY in c(1:nrow(mktcap_Year))){
  Temp = mktcap[as.numeric(substr(mktcap$Date,1,4)) == unique_Year[iY],]
  mktcap_Year[iY,-1] = Temp[nrow(Temp),-1]
}
colnames(mktcap_Year)[1] = "Year"

for (i in c(2:ncol(all_Vol_High))){
  Temp = data.frame(Date = all_Vol_High[, 1], Vol_High_VaR = all_Vol_High[, i])
  Temp$Year = as.numeric(substr(Temp$Date,1,4))
  iID = Sample[i-1]
  Result_List <- list()
  for(iIF in cIndi){
    Temp_IF <- Data_IF[[iIF]]  # 获取当前指标的数据
    Temp_IF <- Temp_IF[Temp_IF$ID == iID, ]  # 筛选 ID
    Temp_IF <- data.frame(Year = 2000:2024, Value = as.numeric(Temp_IF[-1]))  # 创建 DataFrame
    colnames(Temp_IF)[2] <- iIF  # 重新命名列
    Temp_IF[,3] = NA
    Temp_IF[c(2:nrow(Temp_IF)),3] = Temp_IF[1:(nrow(Temp_IF)-1),2]
    colnames(Temp_IF)[3] = paste0(iIF,"_LastY")
    # 存入 list
    Result_List[[iIF]] <- Temp_IF
  }
  
  # 合并所有数据为一个 data.frame
  Final_Data <- Reduce(function(x, y) merge(x, y, by = "Year", all = TRUE), Result_List)
  
  # add MV
  Temp_MV = mktcap_Year[,c(1,i)]
  colnames(Temp_MV)[2] = "MV"
  Temp_MV$MV = Temp_MV$MV * 10^8
  Temp_MV[,3] = NA
  Temp_MV[c(2:nrow(Temp_MV)),3] = Temp_MV[1:(nrow(Temp_MV)-1),2]
  colnames(Temp_MV)[3] = "MV_LastY"
  
  Temp = merge(Temp, Final_Data, all.x = TRUE)
  Temp = merge(Temp, Temp_MV, all.x = TRUE)
  
  Temp$ID = iID
  
  if (i == 2){
    Data_Reg = Temp
  }else{
    Data_Reg = rbind(Data_Reg,Temp) 
  }
}

Data_Reg = merge(Data_Reg, mean_vol_Low, all.x = TRUE)
Data_Reg$Date = as.Date(Data_Reg$Date)

# add FRM 1% 5% 10%
channel = c("HighCarbonIntens", "LowCarbonIntens")
for (iCha in channel){
  output_path = paste0("Output/", iCha)
  FRM_5 = read.csv(file = paste0(output_path,"/Lambda/FRM_",iCha,"_index.csv"), header = TRUE) %>% data.frame()
  colnames(FRM_5)[2] = paste0("FRM_5_",iCha)
  FRM_1 = read.csv(file = paste0(output_path,"/Sensitivity/tau=1/s=63/Lambda/FRM_",iCha,"_index.csv"), header = TRUE) %>% data.frame()
  colnames(FRM_1)[2] = paste0("FRM_1_",iCha)
  FRM_10 = read.csv(file = paste0(output_path,"/Sensitivity/tau=10/s=63/Lambda/FRM_",iCha,"_index.csv"), header = TRUE) %>% data.frame()
  colnames(FRM_10)[2] = paste0("FRM_10_",iCha)
  
  FRM_out = cbind(FRM_5,FRM_1[,2],FRM_10[,2])
  colnames(FRM_out) = c("Date",paste0("FRM_5_",iCha),paste0("FRM_1_",iCha),paste0("FRM_10_",iCha))
  if (iCha == channel[1]){
    FRM_Reg = FRM_out
  }else{
    FRM_Reg = cbind(FRM_Reg, FRM_out[,-1])
  }
}  

FRM_Reg$FRM_5_High_Low = FRM_Reg$FRM_5_HighCarbonIntens / FRM_Reg$FRM_5_LowCarbonIntens
FRM_Reg$FRM_10_High_Low = FRM_Reg$FRM_10_HighCarbonIntens / FRM_Reg$FRM_10_LowCarbonIntens
FRM_Reg$FRM_1_High_Low = FRM_Reg$FRM_1_HighCarbonIntens / FRM_Reg$FRM_1_LowCarbonIntens
FRM_Reg$Date = as.Date(FRM_Reg$Date)

Data_Reg = merge(Data_Reg, FRM_Reg, all.x = TRUE)

write_dta(Data_Reg, paste0(save,"/InsampleBeta_CarbonRisk_VaR.dta"))
saveRDS(Data_Reg, file = paste0(save, "/InsampleBeta_CarbonRisk_VaR.rds")) 

# aggregate
Data_Reg = readRDS(paste0(save, "/InsampleBeta_CarbonRisk.rds"))
Unique_FRM = unique(Data_Reg[,c("Date","FRM_5_High_Low","FRM_10_High_Low","FRM_1_High_Low")])
Unique_FRM$Date = as.character(Unique_FRM$Date )
macro = read.csv(file = paste0('Input/HighCarbonIntens/20140704-20250127/HighCarbonIntens_Macro_20250127.csv'), header = TRUE) %>% data.frame()

merged_df <- list(
  mean_vol_Low,
  mean_vol_High,
  Unique_FRM,
  macro
) %>%
  reduce(left_join, by = "Date")

merged_df$Date <- as.Date(merged_df$Date)
merged_df <- merged_df %>%
  arrange(Date)
write_dta(merged_df, paste0(save,"/Insample_VolatilityVaR_CarbonRisk.dta"))
saveRDS(merged_df, file = paste0(save, "/Insample_VolatilityVaR_CarbonRisk.rds")) 



# volatility, sharpe vola,FRM, macro
Data_Reg = readRDS(paste0(save, "/InsampleBeta_CarbonRisk.rds"))
Unique_FRM = unique(Data_Reg[,c("Date","FRM_5_High_Low","FRM_10_High_Low","FRM_1_High_Low")])
Unique_FRM$Date = as.character(Unique_FRM$Date )

macro = read.csv(file = paste0('Input/HighCarbonIntens/20140704-20250127/HighCarbonIntens_Macro_20250127.csv'), header = TRUE) %>% data.frame()
all_Sharpe_vol_High = read.csv(file = paste0(save, "/HighCarbonIntens/SharpeVolatility.csv"), header = TRUE)%>% data.frame()
all_Sharpe_vol_High = all_Sharpe_vol_High[,-1]
all_Sharpe_vol_Low = read.csv(file = paste0(save, "/LowCarbonIntens/SharpeVolatility.csv"), header = TRUE)%>% data.frame()
all_Sharpe_vol_Low = all_Sharpe_vol_Low[,-1]

all_Vol_High = read.csv(file = paste0(save, "/HighCarbonIntens/Volatility.csv"), header = TRUE)%>% data.frame()
all_Vol_High = all_Vol_High[,-1]
all_Vol_Low = read.csv(file = paste0(save, "/LowCarbonIntens/Volatility.csv"), header = TRUE)%>% data.frame()
all_Vol_Low = all_Vol_Low[,-1]

mean_Sharpe_vol_Low = data.frame(Date = all_Sharpe_vol_Low$Date, mean_Sharpe_vol_Low = apply(all_Sharpe_vol_Low[,-1], 1, mean,na.rm = TRUE))
mean_vol_Low = data.frame(Date = all_Vol_Low$Date, mean_vola_Low = apply(all_Vol_Low[,-1], 1, mean, na.rm = TRUE))

mean_Sharpe_vol_High = data.frame(Date = all_Sharpe_vol_High$Date, mean_Sharpe_vol_High = apply(all_Sharpe_vol_High[,-1], 1, mean,na.rm = TRUE))
mean_vol_High = data.frame(Date = all_Vol_High$Date, mean_vola_High = apply(all_Vol_High[,-1], 1, mean, na.rm = TRUE))

library(dplyr)
library(purrr)

merged_df <- list(
  mean_Sharpe_vol_Low,
  mean_vol_Low,
  mean_Sharpe_vol_High,
  mean_vol_High,
  Unique_FRM,
  macro
) %>%
  reduce(left_join, by = "Date")

merged_df$Date <- as.Date(merged_df$Date)
merged_df <- merged_df %>%
  arrange(Date)
write_dta(merged_df, paste0(save,"/Insample_Volatility_CarbonRisk.dta"))
saveRDS(merged_df, file = paste0(save, "/Insample_Volatility_CarbonRisk.rds")) 

#  Update 1013
#  From the Large / Small group divide into H M L carbon 
# sort portfolio 
# sort by Market value
# High_Intensity_Indus = read.csv(file = paste0('Output/HighCarbonIndustry.csv'), header = TRUE) %>% data.frame()
# Low_Intensity_Indus = read.csv(file = paste0('Output/LowCarbonIndustry.csv'), header = TRUE) %>% data.frame()
# Median_Intensity_Indus = read.csv(file = paste0('Output/MedianCarbonIndustry.csv'), header = TRUE) %>% data.frame()

# read 
channel = c("HighCarbonIntens", "LowCarbonIntens","MedCarbonIntens")
IDlist = list()
for (iCha in channel){
  input_path = paste0("Input/", iCha, "/", date_start_source, "-", date_end_source)
  
  mktcap = read.csv(file = paste0(input_path, "/", iCha, "_Mktcap_", 
                                  date_end_source, ".csv"), header = TRUE) %>% data.frame()
  stock_prices = read.csv(file = paste0(input_path, "/", iCha, "_Price_", 
                                        date_end_source, ".csv"), header = TRUE)%>% data.frame()
  
  all_prices = stock_prices
  all_prices[, -1] = sapply(all_prices[, -1], as.numeric)
  # write.csv(all_prices, paste0(wdir, "/MethodAdd/Price_",save_date,".csv"), quote = FALSE)
  
  all_return = all_prices[-1, ]
  all_return[,-1] = 0
  
  all_return[, -1] = 
    as.data.frame(diff(log(as.matrix(all_prices[, -1])))) 
  all_return[is.na(all_return)] = 0
  
  Mean_mktcap = colMeans(mktcap[,-1],na.rm=TRUE)
  IDlist[[iCha]] = colnames(all_prices)[-1]
  if(iCha == channel[1]){
    all_return_Total = all_return
    Mean_mktcap_Total = Mean_mktcap
  }else{
    all_return_Total = cbind(all_return_Total,all_return[,-1])
    Mean_mktcap_Total = c(Mean_mktcap_Total, Mean_mktcap)
  }
}
colnames(all_return_Total)[-1] = gsub("X", "", colnames(all_return_Total)[-1])
names(Mean_mktcap_Total)= gsub("X", "",names(Mean_mktcap_Total))

Mean_mktcap_Total = sort(Mean_mktcap_Total)

Small = names(Mean_mktcap_Total)[1:floor(length(Mean_mktcap_Total)/2)]
Large = names(Mean_mktcap_Total)[(floor(length(Mean_mktcap_Total)/2)+1):length(Mean_mktcap_Total)] 

# each group divide into high, median,low 
# combine with carbon emission data
Carbon_Rank = readRDS("Output/Carbon_Rank.rds")
FRM_Reg = readRDS(paste0(save,"/FRM_Reg.rds"))

double_sort <- function(Type, Carbon_Rank){
  Out = data.frame(ID = Type)
  Out$ID <- gsub("X", "", Out$ID)
  Out = merge(Out, Carbon_Rank, all.x = TRUE)
  
  quantiles <- quantile(Out$CarbonIntensity_Mean, probs = c(0.3, 0.7))
  
  intensity_groups <- list(
    Low_Intensity = Out$ID[Out$CarbonIntensity_Mean < quantiles[1]],
    Medium_Intensity = Out$ID[Out$CarbonIntensity_Mean >= quantiles[1] 
                                       & Out$CarbonIntensity_Mean <= quantiles[2]],
    High_Intensity = Out$ID[Out$CarbonIntensity_Mean > quantiles[2]]
  )
  
  Out$CarbonType = 0
  Out$CarbonType[Out$ID %in% intensity_groups$High_Intensity]= 1
  Out$CarbonType[Out$ID %in% intensity_groups$Medium_Intensity]= 2
  Out$CarbonType[Out$ID %in% intensity_groups$Low_Intensity]= 3
  
  Out = Out[,c("ID","CarbonType")]
  return(Out)
}

Small_Carbon = double_sort(Small, Carbon_Rank)
Large_Carbon = double_sort(Large, Carbon_Rank)

Small_High = Small_Carbon$ID[Small_Carbon$CarbonType == 1]
Large_High = Large_Carbon$ID[Large_Carbon$CarbonType == 1]

Small_Low = Small_Carbon$ID[Small_Carbon$CarbonType == 3]
Large_Low = Large_Carbon$ID[Large_Carbon$CarbonType == 3]

# 存储不同类型的加权收益
Port_results = list()

for (iType in c("Large_Low","Small_Low","Large_High","Small_High")) {
  
  stock_list = get(iType)  # 确保变量存在
  
  # 选取对应股票的收益
  Port_return = all_return_Total[, c("Date", stock_list)]
  
  # 选取对应股票的市值
  # 提示缺失股票
  missing_stocks <- setdiff(stock_list, names(Mean_mktcap_Total))
  if(length(missing_stocks) > 0){
    warning("以下股票缺失市值，将被忽略: ", paste(missing_stocks, collapse=", "))
  }
  
  Port_MV = Mean_mktcap_Total[names(Mean_mktcap_Total) %in% stock_list]
  Port_MV = Port_MV[stock_list]  # 确保顺序匹配
  Port_MV = Port_MV / sum(Port_MV)  # 归一化权重
  
  # 计算加权收益
  Weighted_return = as.matrix(Port_return[, -1]) %*% Port_MV  
  Weighted_return = data.frame(Date = Port_return$Date, Return = Weighted_return)
  
  # 存入列表
  Port_results[[iType]] = Weighted_return
}

CarbonPremia = 0.5*(Port_results[["Large_Low"]]["Return"]+Port_results[["Small_Low"]]["Return"])- 0.5*(Port_results[["Large_High"]]["Return"]+Port_results[["Small_High"]]["Return"])
CarbonPremia = data.frame(Date = Port_results[["Large_Low"]]["Date"], Premia = CarbonPremia)
CarbonPremia$Date = as.Date(CarbonPremia$Date)

FRM_Reg$Date = as.Date(FRM_Reg$Date)
CarbonPremia = merge(CarbonPremia, FRM_Reg, all.x = TRUE)
CarbonPremia = CarbonPremia[(!is.na(CarbonPremia$FRM_5_HighCarbonIntens))&(!is.infinite(CarbonPremia$Return)),]

macro = read.csv(file = paste0('Input/HighCarbonIntens/20140704-20250127/HighCarbonIntens_Macro_20250127.csv'), header = TRUE) %>% data.frame()
macro$Date = as.Date(macro$Date)
CarbonPremia = merge(CarbonPremia, macro, all.x = TRUE)

Events <- read_excel("Input/Event/Important_Carbon_Events.xlsx") %>% data.frame
Date_Events = as.Date(Events$Date)

CarbonPremia$Event = 0
CarbonPremia$Event[CarbonPremia$Date %in% Date_Events] = 1

for (d in Date_Events[Date_Events %in% CarbonPremia$Event[CarbonPremia$Event == 0]]) {
  # 找到 CarbonPremia 中比事件日期大的最小日期
  next_date <- min(CarbonPremia$Date[CarbonPremia$Date > d])
  if (!is.infinite(next_date)) {
    CarbonPremia$Event[CarbonPremia$Date == next_date] <- 1
  }
}

# 假设数据是时间序列 ts 对象
R_C_ts <- ts(CarbonPremia$Return)  # R_t^C
Mkt_ts <- ts(CarbonPremia$MKreturn)   # 市场因子
Event_ts <- ts(CarbonPremia$Event)    # 碳因子

p <- 4  # 假设滞后阶数为 3

# 构建滞后变量
df <- CarbonPremia %>%
  arrange(Date) %>%
  mutate(
    R_C_lag1 = lag(Return, 1),
    Mkt_lag1 = lag(MKreturn, 1), Mkt_lag2 = lag(MKreturn, 2), Mkt_lag3 = lag(MKreturn, 3), Mkt_lag4 = lag(MKreturn, 4),
    Event_lag1 = lag(Event, 1), Event_lag2 = lag(Event, 2), Event_lag3 = lag(Event, 3), Event_lag4 = lag(Event, 4)
  )

# 回归
model <- lm(Return ~ R_C_lag1 + Mkt_lag1 + Mkt_lag2 + Mkt_lag3 + Mkt_lag4 +
              Event_lag1 + Event_lag2 + Event_lag3 + Event_lag4,
            data = df, na.action = na.exclude)

CarbonPremia$CResidual <- residuals(model)

# CoVaR of FRM on Premium?
# Functions
calcRMA <- function(serie,tau, h = 63){
  n <- length(serie)
  rma <- rep(NA, n)
  first_nonNA <- which(!is.na(serie))[1]
  for(i in (first_nonNA + h - 1):n){
    window <- serie[(i - h + 1):i]
    rma[i] <- mean(window^2, na.rm = TRUE)^ 0.5 * qnorm(tau)+mean(window,na.rm = TRUE)
  }
  return(rma)
}


emaWeight = function(h, gamma){
  weight = gamma^(h: 1) * (1 - gamma)
  return(weight)
}

calcEMA <- function(serie, tau, gamma = 0.94, h = 63){

  n <- length(serie)
  ema <- rep(NA, n)
  first_nonNA <- which(!is.na(serie))[1]
  weight <- emaWeight(h, gamma)
  for(i in (first_nonNA + h - 1):n){
    window <- serie[(i - h + 1):i]
    ema[i] <- sum(window^2 * weight, na.rm = TRUE)^ 0.5 * qnorm(tau)+mean(window)
  }
  return(ema)
}


# calculate VaR
vari = "Return"
CarbonPremia$VaRRMAY_95 = calcRMA(CarbonPremia[, vari], 0.95) 
CarbonPremia$VaREMAY_95 = calcEMA(CarbonPremia[, vari], 0.95) 

CarbonPremia$VaRRMAY_99 = calcRMA(CarbonPremia[, vari], 0.99) 
CarbonPremia$VaREMAY_99 = calcEMA(CarbonPremia[, vari], 0.99) 

CarbonPremia$VaRRMAY_90 = calcRMA(CarbonPremia[, vari], 0.90) 
CarbonPremia$VaREMAY_90 = calcEMA(CarbonPremia[, vari], 0.90) 

CarbonPremia$VaRRMAY_5 = calcRMA(CarbonPremia[, vari], 0.05) 
CarbonPremia$VaREMAY_5 = calcEMA(CarbonPremia[, vari], 0.05) 

CarbonPremia$VaRRMAY_1 = calcRMA(CarbonPremia[, vari], 0.01) 
CarbonPremia$VaREMAY_1 = calcEMA(CarbonPremia[, vari], 0.01) 

CarbonPremia$VaRRMAY_10 = calcRMA(CarbonPremia[, vari], 0.1) 
CarbonPremia$VaREMAY_10 = calcEMA(CarbonPremia[, vari], 0.1) 

# Quantile R

h <- 63
n <- nrow(CarbonPremia)
first_nonNA <- which(!is.na(CarbonPremia$Return))[1]

# 初始化存储列
taus <- c(0.01, 0.05, 0.10, 0.90, 0.95, 0.99)
for (tau in taus) {
  CarbonPremia[[paste0("QrVaR_", tau * 100)]] <- NA
}

# 滚动分位数回归
for (i in (first_nonNA + h - 1):n) {
  window <- CarbonPremia[(i - h + 1):i, ]
  
  for (tau in taus) {
    qr_model <- rq(
      Return ~ Change_TY3M + Slope + TED + RealEstate_excess + 
        MKreturn + MKvol + CarbonVol_Shenzhen + CarbonVol_Guangdong + CarbonVol_Hubei,
      tau = tau, data = window
    )
    
    # 使用窗口最后一行数据预测下一期VaR
    CarbonPremia[i, paste0("QrVaR_", tau * 100)] <- 
      predict(qr_model, newdata = CarbonPremia[i, ])
  }
}

write_dta(CarbonPremia, paste0(save,"/InsampleBeta_Premia_1014_daily.dta"))

# Plot Premia  vs FRM rolling window = 63 or loess 
nLag = 63
# CarbonPremia_Month = rollapply(CarbonPremia[,-1], width = nLag, FUN = mean, align = "right", fill = NA)%>% data.frame()
# CarbonPremia_Month$Date = CarbonPremia[,1]
CarbonPremia_Month = CarbonPremia
CarbonPremia_Month[,"Return"] = rollapply(CarbonPremia[,"Return"], width = nLag, FUN = mean, align = "right", fill = NA)%>% data.frame()

CarbonPremia_Month = CarbonPremia_Month[-c(1:(nLag-1)),]
CarbonPremia_Month$Return = CarbonPremia_Month$Return*100
CarbonPremia_Month$Year = format(CarbonPremia_Month$Date, "%Y")

# # add macro
# macro = read.csv(file = paste0('Input/HighCarbonIntens/20140704-20250127/HighCarbonIntens_Macro_20250127.csv'), header = TRUE) %>% data.frame()
# macro$Date = as.Date(macro$Date)
# CarbonPremia_Month = merge(CarbonPremia_Month, macro, all.x = TRUE)

write_dta(CarbonPremia_Month, paste0(save,"/InsampleBeta_Premia_1014_month.dta"))
# Plot Premia  vs FRM rolling window = 63 or loess 
nLag = 63
# CarbonPremia_Month = rollapply(CarbonPremia[,-1], width = nLag, FUN = mean, align = "right", fill = NA)%>% data.frame()
# CarbonPremia_Month$Date = CarbonPremia[,1]
CarbonPremia_Month = CarbonPremia
CarbonPremia_Month[,"Return"] = rollapply(CarbonPremia[,"Return"], width = nLag, FUN = mean, align = "right", fill = NA)%>% data.frame()

CarbonPremia_Month = CarbonPremia_Month[-c(1:(nLag-1)),]
CarbonPremia_Month$Return = CarbonPremia_Month$Return*100
CarbonPremia_Month$Year = format(CarbonPremia_Month$Date, "%Y")

# # add macro
# macro = read.csv(file = paste0('Input/HighCarbonIntens/20140704-20250127/HighCarbonIntens_Macro_20250127.csv'), header = TRUE) %>% data.frame()
# macro$Date = as.Date(macro$Date)
# CarbonPremia_Month = merge(CarbonPremia_Month, macro, all.x = TRUE)

write_dta(CarbonPremia_Month, paste0(save,"/InsampleBeta_Premia_1014_month.dta"))

# plot
Events <- read_excel("Input/Event/Important_Carbon_Events.xlsx") %>% data.frame
Date_Events = as.Date(Events$Date)


png(paste0("Output/Compare_FRM_CarbonPremia_1013.png"), width = 900, height = 600, bg = "transparent")

# geom_vline(
#   xintercept = Date_Events, 
#   linetype = "dashed", 
#   color = "blue", 
#   size = 0.5  # Vertical dashed lines for event dates
# )
ggplot(CarbonPremia_Month, aes(x = Date)) +
  # Primary axis (CarbonRisk or other metrics)
  geom_line(aes(y = FRM_1_High_Low, color = "CarbonRisk_tau_1"), linewidth = 1.2) +
  geom_line(aes(y = FRM_5_High_Low, color = "CarbonRisk_tau_5"), linewidth = 1.2) +
  geom_line(aes(y = FRM_10_High_Low, color = "CarbonRisk_tau_10"), linewidth = 1.2) +
  # Add dashed horizontal lines at y = 0 and y = 1
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 1) + 
  # Secondary axis (Return mapped to Carbon Premia)
  geom_line(aes(y = Return, color = "Carbon Premia"), linewidth = 1.2) + 
  # geom_vline(xintercept = Date_Events, linetype = "dashed", color = "blue", size = 0.5)+
  # Set x-axis breaks
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")+
   # Primary y-axis for CarbonRisk and secondary y-axis for Return (Carbon Premia)
  scale_y_continuous(name = "Carbon Risk", 
                     sec.axis = sec_axis(trans = ~ ./100, name = "Carbon Premia")) +  # No transformation, just use original scale for Carbon Premia
  # Custom theme
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent"),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )

dev.off()



Fun_ROOS = function (Index, i_type, select, horizon){
  R_oos_out = matrix( nrow = length(horizon), ncol = length(select)+1)
  
  R_oos_out[,1] = horizon
  colnames(R_oos_out) = c('lagDay', paste0('Roos_', select))
  nstart = 1
  
  for (iLag in horizon){

    for (j in select){
      temp_reg = cbind(Index[(iLag+1):nrow(Index), i_type], Index[1:(nrow(Index)-iLag), j])
      colnames(temp_reg) = c('y','x')
      
      temp_reg = data.frame(temp_reg)
      temp_reg = temp_reg[!is.na(temp_reg[,'x']),]
      
      for (i in c(63:(nrow(temp_reg)-1))){
        
        fit = lm(y ~ x, temp_reg[(i-62):(i-1), ])
        
        temp_reg[i, paste0('p_vol_',j)] = fit$coefficients[1] + fit$coefficients[2]*temp_reg$x[i]
        temp_reg[i, paste0('diff_s_',j)] = (temp_reg$y[i]-temp_reg[i,paste0('p_vol_',j)])^2
        # temp_reg[i, paste0('diff_m_',j)] = (temp_reg$y[i]-mean(temp_reg$y))^2
        temp_reg[i, paste0('diff_m_',j)] = (temp_reg$y[i]-mean(temp_reg$y[63:nrow(temp_reg)]))^2
      }
      
      R_oos_out[nstart, paste0('Roos_', j)] = 1 - sum(temp_reg[,paste0('diff_s_',j)], na.rm = TRUE)/sum(temp_reg[,paste0('diff_m_',j)], na.rm = TRUE)
      
    }
    
    nstart = nstart + 1
    print(paste0('finish ',iLag, ' Type ', i_type))
  }
  write.csv(R_oos_out, paste0("Output/Prediction_Continuous",i_type,"_ROO.csv"), quote = FALSE)
  Mean_Roos = rbind(colMeans(R_oos_out[,-1]),colSds(R_oos_out[,-1]))
  write.csv(Mean_Roos, paste0("Output/Mean_Continuous",i_type,"_vol_ROO.csv"), quote = FALSE)
  
}

horizon = c(5,21,63,84,105)
select = c("FRM_1_High_Low","FRM_5_High_Low","FRM_10_High_Low")
Fun_ROOS(CarbonPremia_Month, 'Return', select, horizon)


#  combine with macro
macro = read.csv(file = paste0('Input/HighCarbonIntens/20140704-20250127/HighCarbonIntens_Macro_20250127.csv'), header = TRUE) %>% data.frame()
macro$Date = as.Date(macro$Date)
Shap_CarbonPremia = merge(CarbonPremia, macro, all.x = TRUE)

# Shapley for carbon risk
Fun_Shapley_ML <-function(iLag, Index,i_type){
  
  
  
  temp_reg = cbind(Testdata[(iLag+1):nrow(Index), i_type], Testdata[1:(nrow(Index)-iLag), select])
  colnames(temp_reg)[1] = i_type
  temp_reg = temp_reg[complete.cases(temp_reg),]
  trainval = as.matrix(temp_reg[,-1])
  
  param_list <- list(objective = "reg:squarederror",  # For regression
                     eta = 0.02,
                     max_depth = 100,
                     gamma = 0.01,
                     subsample = 0.8,
                     colsample_bytree = 0.86)
  
  for (i in c(1:500)){
    mod <- xgboost::xgboost(data = trainval, 
                            label = as.matrix(temp_reg[, i_type]), 
                            params = param_list, nrounds = 200,
                            verbose = FALSE,
                            early_stopping_rounds = 8)
    shap_long <- shap.prep(xgb_model = mod, X_train = trainval)
    Shapley_Mean_Temp = unique(shap_long[,c('variable','mean_value')])
    colnames(Shapley_Mean_Temp)[2] = paste0('mean_value_',i)
    if (i == 1){
      Shapley_Mean = Shapley_Mean_Temp
    }else{
      Shapley_Mean = merge(Shapley_Mean, Shapley_Mean_Temp)
    }
  }
  Shapley_Mean$mean_value = rowMeans(Shapley_Mean[,-1])
  Shapley_Mean$variable = as.character(Shapley_Mean$variable)
  Shapley_Mean$variable[Shapley_Mean$variable == 'FRM_5P'] = 'FRM 5%'
  Shapley_Mean$variable[Shapley_Mean$variable == 'FRM_25P'] = 'FRM 25%'
  Shapley_Mean$variable[Shapley_Mean$variable == 'FRM_50P'] = 'FRM 50%'
  Shapley_Mean$variable[Shapley_Mean$variable == 'TCI'] = 'Total Connectedness'
  Shapley_Mean$variable[Shapley_Mean$variable == 'DAG'] = 'Bayesian Graphical VAR'
  Shapley_Mean$variable[Shapley_Mean$variable == 'PCA_1'] = 'Principal Components'
  Shapley_Mean$variable[Shapley_Mean$variable == 'GDC'] = 'Granger Causality'
  
  Shapley_Mean$Shapley_Mean = round(Shapley_Mean$mean_value,5)
  
  p <- ggplot(Shapley_Mean,aes(x= reorder(variable, Shapley_Mean),y=Shapley_Mean, fill = variable))+geom_bar(stat="identity", width=0.5)+
    coord_flip()+
    geom_text(aes(label=Shapley_Mean), vjust=-2, size=5)+
    labs(y = "Shapley Value", size = 20)+
    labs(x = "Risk Measures", size = 20)+
    theme(axis.text.x = element_text(size =15, angle = 30, hjust = 1), 
          axis.text.y = element_text(size = 15), 
          axis.title =  element_text(size=,face = "bold"),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 35, face = "bold"),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          axis.line = element_line(colour = "black"),legend.position="none")
  
  # png(paste0(save, "/Shapley/SHAP_",i_type,"_lag_",iLag,".png"), width = 900, height = 600, bg = "transparent")
  # 
  # print(p)
  # dev.off()
  
  ggsave(paste0(save, "/Shapley/SHAP_",i_type,"_lag_",iLag,".eps"), plot = p, width = 12, height = 8, bg = "transparent",dpi = 800)
  
  # write.csv(Shapley_Mean, paste0(save, "/Shapley/ShapleyCalData_",i_type,"_lag_",iLag,"_",save_date,".csv"), quote = FALSE)
  
}

Fun_Shapley_fixed_ML <-function(iLag, Index,i_type, select){
  
  temp_reg = cbind(Index[(iLag+1):nrow(Index), i_type], Index[1:(nrow(Index)-iLag), select])
  colnames(temp_reg)[1] = i_type
  temp_reg = temp_reg[complete.cases(temp_reg),]
  trainval = as.matrix(temp_reg[,-1])
  
  
  param_list <- list(objective = "reg:squarederror",  # For regression
                     eta = 0.02,
                     max_depth = 100,
                     gamma = 0.01,
                     subsample = 0.8,
                     colsample_bytree = 0.86)
  
  for (i in c(1:500)){
    mod <- xgboost::xgboost(data = trainval, 
                            label = as.matrix(temp_reg[, i_type]), 
                            params = param_list, nrounds = 200,
                            verbose = FALSE,
                            early_stopping_rounds = 8)
    shap_long <- shap.prep(xgb_model = mod, X_train = trainval)
    Shapley_Mean_Temp = unique(shap_long[,c('variable','mean_value')])
    colnames(Shapley_Mean_Temp)[2] = paste0('mean_value_',i)
    if (i == 1){
      Shapley_Mean = Shapley_Mean_Temp
    }else{
      Shapley_Mean = merge(Shapley_Mean, Shapley_Mean_Temp)
    }
  }
  
  Shapley_Mean$mean_value = rowMeans(Shapley_Mean[,-1])
  Shapley_Mean$variable = as.character(Shapley_Mean$variable)
 
  Shapley_Mean$Shapley_Mean = round(Shapley_Mean$mean_value,5)
  
  p = ggplot(Shapley_Mean,aes(x= reorder(variable, Shapley_Mean),y=Shapley_Mean, fill = variable))+geom_bar(stat="identity", width=0.5)+
    coord_flip()+
    geom_text(aes(label=Shapley_Mean), vjust=-2, size=5)+
    labs(y = "Shapley Value", size = 20)+
    labs(x = "Risk Measures", size = 20)+
    theme(axis.text.x = element_text(size =15, angle = 30, hjust = 1),
          axis.text.y = element_text(size = 15),
          axis.title =  element_text(size=,face = "bold"),
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 35, face = "bold"),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          axis.line = element_line(colour = "black"),legend.position="none")
  p
  
  ggsave(paste0(save, "/Shapley/SHAP_",i_type,"_lag_",iLag,".eps"), plot = p, width = 12, height = 8, bg = "transparent",dpi = 800)
  ggsave(paste0(save, "/Shapley/SHAP_",i_type,"_lag_",iLag,".png"), plot = p, width = 12, height = 8, bg = "transparent",dpi = 800)
  
  write.csv(Shapley_Mean, paste0(save, "/Shapley/ShapleyCalData_",i_type,"_lag_",iLag,".csv"), quote = FALSE)
  
}

select = c("Change_TY3M","Slope", "TED", "RealEstate_excess","MKreturn","MKvol",
           "CarbonVol_Shenzhen","CarbonVol_Guangdong","CarbonVol_Hubei")

Fun_Shapley_fixed_ML(iLag = 0, Index = Shap_CarbonPremia,
                     i_type = "FRM_5_High_Low", select)







CarbonPremia_Month$FRM_5_High_Low_lag21 = lag(CarbonPremia_Month$FRM_5_High_Low,21)

CarbonPremia_Month$FRM_1_High_Low_lag21 = lag(CarbonPremia_Month$FRM_1_High_Low,21)
CarbonPremia_Month$FRM_10_High_Low_lag21 = lag(CarbonPremia_Month$FRM_10_High_Low,21)

cor_test_FRM_1 = cor.test(CarbonPremia_Month$Return, CarbonPremia_Month$FRM_1_High_Low_lag21, 
                          method = "pearson", use = "complete.obs")

cor_test_FRM_5 = cor.test(CarbonPremia_Month$Return, CarbonPremia_Month$FRM_5_High_Low_lag21, 
                          method = "pearson", use = "complete.obs")

cor_test_FRM_10 = cor.test(CarbonPremia_Month$Return, CarbonPremia_Month$FRM_10_High_Low_lag21, 
                           method = "pearson", use = "complete.obs")

# Print the results
cor_test_FRM_1
cor_test_FRM_5
cor_test_FRM_10


model <- lm(Return ~ FRM_5_High_Low_lag21, data = CarbonPremia_Month)
summary(model)

Y = ts(CarbonPremia_Month$Return)
X = ts(CarbonPremia_Month$FRM_1_High_Low)
data = cbind(X,Y)
data = data.frame(data)

lag_selection <- VARselect(data, lag.max = 20, type = "const")
optimal_lag_aic <- which.min(lag_selection$criteria["AIC(n)", ])
optimal_lag_aic
# Test Granger causality
# granger_test <- grangertest(Y~X, order = optimal_lag_aic)
granger_test <- grangertest(Y~X, order = 5)

print(granger_test)

data$Y_lag1 <- lag(data$Y, 1)
data$Y_lag2 <- lag(data$Y, 2)
data$Y_lag3 <- lag(data$Y, 3)
data$Y_lag4 <- lag(data$Y, 4)
data$Y_lag5 <- lag(data$Y, 5)

data$X_lag1 <- lag(data$X, 1)
data$X_lag2 <- lag(data$X, 2)
data$X_lag2 <- lag(data$X, 2)
data$X_lag3 <- lag(data$X, 3)
data$X_lag4 <- lag(data$X, 4)
data$X_lag5 <- lag(data$X, 5)

model <- lm(Y ~ Y_lag1 + Y_lag2 + Y_lag3+ Y_lag4+ Y_lag5 + X_lag1 + X_lag2 + X_lag3 + X_lag4 + X_lag5, data = data)

# View the summary to get the coefficients
summary(model)

plot(Y)
plot(X)

# 结果：carbon_premia_results[[“Low_Large”]]、carbon_premia_results[[“Low_Small”]] 等分别存储不同类型的加权收益

mktcap[, -1] = sapply(mktcap[, -1], as.numeric)
mktcap_weight = mktcap
mktcap_weight[,-1] = 0
mktcap_weight[,-1] = mktcap[,-1]/rowSums(mktcap[, -1],na.rm = TRUE)
mktcap_weight = mktcap_weight[-1,]

all_return_temp = all_return
all_return_temp[,-1]= all_return_temp[,-1]*mktcap_weight[,-1]
all_return_temp[is.na(all_return_temp)] = 0
Fin_Index = data.frame(date = all_return_temp$date, Return = rowSums(all_return_temp[, -1]))


Fun_Vola = function (data){

  vola = NA
  for (i in c(63 : length(data)) ){
    vola = c(vola, sd(data[(i-62) : i]))
  }
  vola = vola[-1]
  
  return(vola)
}

# Fin_Index$MKVola = apply(Fin_Index$Return, width = 63, FUN = sd, fill = NA, align = "right") 
Fin_Index$MKVola = NA
Fin_Index$MKVola[63:nrow(Fin_Index)] <- Fun_Vola(Fin_Index$Return)
# Fin_Index = merge(Fin_Index,crix, all.x = TRUE)


Test_Fin_Index = Fin_Index
Test_Fin_Index$MKVola_L_10 = NA
Test_Fin_Index$MKVola_L_25 = NA
Test_Fin_Index$MKVola_L_110 = NA


Test_Fin_Index[11:nrow(Test_Fin_Index), 'MKVola_L_10'] = Test_Fin_Index[1:(nrow(Test_Fin_Index)-10), 'MKVola']
Test_Fin_Index[26:nrow(Test_Fin_Index), 'MKVola_L_25'] = Test_Fin_Index[1:(nrow(Test_Fin_Index)-25), 'MKVola']
Test_Fin_Index[64:nrow(Test_Fin_Index), 'MKVola_L_63'] = Test_Fin_Index[1:(nrow(Test_Fin_Index)-63), 'MKVola']
Test_Fin_Index[111:nrow(Test_Fin_Index), 'MKVola_L_110'] = Test_Fin_Index[1:(nrow(Test_Fin_Index)-110), 'MKVola']


# read FRM
FRM_5P  = data.frame(read.csv(paste0(wdir,'/Output/Crypto/Lambda/FRM_Crypto_index.csv')))
colnames(FRM_5P)[2] = 'FRM_5P'

FRM_25P  = data.frame(read.csv(paste0(wdir,'/Output/Crypto/Lambda/Quantiles/q25_lambda.csv')))
colnames(FRM_25P)[2] = 'FRM_25P'

FRM_50P  = data.frame(read.csv(paste0(wdir,'/Output/Crypto/Lambda/Quantiles/q50_lambda.csv')))
colnames(FRM_50P)[2] = 'FRM_50P'


# read BGVAR
BGVAR  = data.frame(read_excel(paste0(wdir,'/MethodAdd/BGVAR_20230514.xlsx')))
colnames(BGVAR)[1] = colnames(FRM_5P)[1]

# read GDC
GDC  = data.frame(read_excel(paste0(wdir,'/MethodAdd/rolling_GDC_total_20230514.xlsx')))
colnames(GDC)[1] = colnames(FRM_5P)[1]
GDC[,1] = as.character(GDC[,1])
colnames(GDC)[2] = 'GDC'
# GDC$date = paste0(substr(GDC$date,1,4), substr(GDC$date,6,7), substr(GDC$date,9,10))
# GDC$date = as.numeric(GDC$date)

# read PCA
PCA  = data.frame(read_excel(paste0(wdir,'/MethodAdd/PCA_20230514.xlsx')))
PCA = PCA[,-1]
PCA[,1] = as.character(PCA[,1])
colnames(PCA)[1] = colnames(FRM_5P)[1]
# PCA$date = paste0(substr(PCA$date,1,4), substr(PCA$date,6,7), substr(PCA$date,9,10))
# PCA$date = as.numeric(PCA$date)

# read DY
DY  = data.frame(read_excel(paste0(wdir,'/MethodAdd/DY_TCI.xlsx')))
DY = DY[,-1]
colnames(DY)[1] = colnames(FRM_5P)[1]

# read MV
MVWealth  = read.csv(paste0(wdir,'/Output/Crypto/Add_PortfolioTest/PortWealth.csv'), header = TRUE)
MVWealth = MVWealth[,c('date','Return_MV')]
MVWealth$date = as.character(as.Date(as.character(MVWealth[,1]),  "%Y%m%d"))

Index = merge(FRM_5P, FRM_25P, by = intersect(names(FRM_5P), names(FRM_25P)),all = FALSE,sort = TRUE)
Index = merge(Index, FRM_50P, by = intersect(names(Index), names(FRM_50P)),all = FALSE,sort = TRUE)
Index = merge(Index, BGVAR, by = intersect(names(Index), names(BGVAR)),all = FALSE,sort = TRUE)
Index = merge(Index, GDC, by = intersect(names(Index), names(GDC)),all = FALSE,sort = TRUE)
Index = merge(Index, PCA, by = intersect(names(Index), names(PCA)),all = FALSE,sort = TRUE)
Index = merge(Index, DY, by = intersect(names(Index), names(DY)),all = FALSE,sort = TRUE)
Index = merge(Index, Fin_Index, by = intersect(names(Index), names(Fin_Index)),all = FALSE,sort = TRUE)

Index_MVWealth = merge(MVWealth, Index, by = intersect(names(MVWealth), names(Index)),all.x = TRUE,sort = TRUE)
Index_MVWealth = Index_MVWealth[63:nrow(MVWealth),]
# x_type = c('MKVola', 'Return')
select = c('FRM_5P','FRM_25P','FRM_50P','DAG','GDC','PCA_1','TCI')

# continuous: from 1day to 110day

Fun_Vola_dynamic = function (data,iLag){
  
  vola = NA
  for (i in c(iLag : length(data)) ){
    vola = c(vola, sd(data[(i-iLag+1) : i]))
  }
  vola = vola[-1]
  
  return(vola)
}




Fun_ROOS_fixed = function (Index, i_type, select, horizon){
  R_oos_out = matrix( nrow = length(horizon), ncol = length(select)+1)
  
  R_oos_out[,1] = horizon
  colnames(R_oos_out) = c('lagDay', paste0('Roos_', select))
  nstart = 1
  
  for (iLag in horizon){
    
    if (i_type == 'Return'){
      Testdata = Index
    }else if(i_type == 'MKVola'){
      Testdata = Index
      Testdata$MKVola = NA
      Testdata$MKVola[63:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, 63)
      # Testdata$MKVola[iLag:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, 21)
    }else{
      Testdata = Index
      Testdata$MVVola = NA
      Testdata$MVVola[63:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return_MV, 63)
    }
    
    for (j in select){
      temp_reg = cbind(Testdata[(iLag+1):nrow(Index), i_type], Testdata[1:(nrow(Index)-iLag), j])
      colnames(temp_reg) = c('y','x')
      
      temp_reg = data.frame(temp_reg)
      temp_reg = temp_reg[!is.na(temp_reg[,'x']),]
      
      for (i in c(63:(nrow(temp_reg)-1))){
        
        fit = lm(y ~ x, temp_reg[(i-62):(i-1), ])
        
        temp_reg[i, paste0('p_vol_',j)] = fit$coefficients[1] + fit$coefficients[2]*temp_reg$x[i]
        temp_reg[i, paste0('diff_s_',j)] = (temp_reg$y[i]-temp_reg[i,paste0('p_vol_',j)])^2
        # temp_reg[i, paste0('diff_m_',j)] = (temp_reg$y[i]-mean(temp_reg$y))^2
        temp_reg[i, paste0('diff_m_',j)] = (temp_reg$y[i]-mean(temp_reg$y[63:nrow(temp_reg)]))^2
      }
      
      R_oos_out[nstart, paste0('Roos_', j)] = 1 - sum(temp_reg[,paste0('diff_s_',j)], na.rm = TRUE)/sum(temp_reg[,paste0('diff_m_',j)], na.rm = TRUE)
      
    }
    
    nstart = nstart + 1
    print(paste0('finish ',iLag, ' Type ', i_type))
  }
  
  if ('MKVola_L_10' %in% select){
    save_csv =  paste0(save, "/Self_Prediction_Continuous",i_type,"_Fixed63_ROO_",save_date,".csv")
    save_csv_mean =  paste0(save, "/Self_Mean_Continuous",i_type,"_Fixed63_ROO_",save_date,".csv")
  }else{
    save_csv =  paste0(save, "/Prediction_Continuous",i_type,"_Fixed63_ROO_",save_date,".csv")
    save_csv_mean =  paste0(save, "/Mean_Continuous",i_type,"_Fixed63_ROO_",save_date,".csv")
  }
  write.csv(R_oos_out, save_csv, quote = FALSE)
  Mean_Roos = rbind(colMeans(R_oos_out[,-1]),colSds(R_oos_out[,-1]))
  write.csv(Mean_Roos, save_csv_mean, quote = FALSE)
  
}


Fun_R2 = function (Index, i_type, select, horizon){
  R2 = matrix( nrow = length(horizon), ncol = length(select)+1)
  
  R2[,1] = horizon
  colnames(R2) = c('lagDay', paste0('R2_', select))
  nstart = 1
  
  for (iLag in horizon){
    
    if (i_type == 'Return'){
      Testdata = Index
    }else if(i_type == 'MKVola'){
      Testdata = Index
      Testdata$MKVola = NA
      Testdata$MKVola[iLag:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, iLag)
      # Testdata$MKVola[iLag:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, 21)
    }else{
      Testdata = Index
      Testdata$MVVola = NA
      Testdata$MVVola[iLag:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return_MV, iLag)
    }
    
    
    for (j in select){
      temp_reg = cbind(Testdata[(iLag+1):nrow(Index), i_type], Testdata[1:(nrow(Index)-iLag), j])
      colnames(temp_reg) = c('y','x')
      temp_reg = data.frame(temp_reg)
      temp_reg = temp_reg[!is.na(temp_reg[,'x']),]
      
      for (i in c(63:(nrow(temp_reg)-1))){
        
        fit = lm(y ~ x, temp_reg[(i-61):i, ])
        
        temp_reg[i, paste0('p_vol_',j)] = fit$coefficients[1] + fit$coefficients[2]*temp_reg$x[i]
        temp_reg[i, paste0('diff_s_',j)] = (temp_reg$y[i]-temp_reg[i,paste0('p_vol_',j)])^2
        # temp_reg[i, paste0('diff_m_',j)] = (temp_reg$y[i]-mean(temp_reg$y))^2
        temp_reg[i, paste0('diff_m_',j)] = (temp_reg$y[i]-mean(temp_reg$y[63:nrow(temp_reg)]))^2
      }
      
      R2[nstart, paste0('R2_', j)] = 1 - sum(temp_reg[,paste0('diff_s_',j)], na.rm = TRUE)/sum(temp_reg[,paste0('diff_m_',j)], na.rm = TRUE)
      
    }
    
    nstart = nstart + 1
    print(paste0('finish ',iLag, ' Type ', i_type))
  }
  write.csv(R2, paste0(save, "/Prediction_Continuous",i_type,"_R2_",save_date,".csv"), quote = FALSE)
  Mean_R2 = rbind(colMeans(R2[,-1]),colSds(R2[,-1]))
  write.csv(Mean_R2, paste0(save, "/Mean_Continuous",i_type,"_R2_",save_date,".csv"), quote = FALSE)
  
}


Fun_R2_fixed = function (Index, i_type, select, horizon){
  R2 = matrix( nrow = length(horizon), ncol = length(select)+1)
  
  R2[,1] = horizon
  colnames(R2) = c('lagDay', paste0('R2_', select))
  nstart = 1
  
  for (iLag in horizon){
    
    if (i_type == 'Return'){
      Testdata = Index
    }else if(i_type == 'MKVola'){
      Testdata = Index
      Testdata$MKVola = NA
      Testdata$MKVola[63:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, 63)
      # Testdata$MKVola[iLag:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, 21)
    }else{
      Testdata = Index
      Testdata$MVVola = NA
      Testdata$MVVola[63:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return_MV, 63)
    }
    
    
    for (j in select){
      temp_reg = cbind(Testdata[(iLag+1):nrow(Index), i_type], Testdata[1:(nrow(Index)-iLag), j])
      colnames(temp_reg) = c('y','x')
      temp_reg = data.frame(temp_reg)
      temp_reg = temp_reg[!is.na(temp_reg[,'x']),]
      
      for (i in c(63:(nrow(temp_reg)-1))){
        
        fit = lm(y ~ x, temp_reg[(i-61):i, ])
        
        temp_reg[i, paste0('p_vol_',j)] = fit$coefficients[1] + fit$coefficients[2]*temp_reg$x[i]
        temp_reg[i, paste0('diff_s_',j)] = (temp_reg$y[i]-temp_reg[i,paste0('p_vol_',j)])^2
        # temp_reg[i, paste0('diff_m_',j)] = (temp_reg$y[i]-mean(temp_reg$y))^2
        temp_reg[i, paste0('diff_m_',j)] = (temp_reg$y[i]-mean(temp_reg$y[63:nrow(temp_reg)]))^2
      }
      
      R2[nstart, paste0('R2_', j)] = 1 - sum(temp_reg[,paste0('diff_s_',j)], na.rm = TRUE)/sum(temp_reg[,paste0('diff_m_',j)], na.rm = TRUE)
      
    }
    
    nstart = nstart + 1
    print(paste0('finish ',iLag, ' Type ', i_type))
  }
  
  if ('MKVola_L_10' %in% select){
    save_csv =  paste0(save, "/Self_Prediction_Continuous",i_type,"_Fixed63_R2_",save_date,".csv")
    save_csv_mean =  paste0(save, "/Self_Mean_Continuous",i_type,"_Fixed63_R2_",save_date,".csv")
  }else{
    save_csv =  paste0(save, "/Prediction_Continuous",i_type,"_Fixed63_R2_",save_date,".csv")
    save_csv_mean =  paste0(save, "/Mean_Continuous",i_type,"_Fixed63_R2_",save_date,".csv")
  }
  write.csv(R2, save_csv, quote = FALSE)
  Mean_R2 = rbind(colMeans(R2[,-1]),colSds(R2[,-1]))
  write.csv(Mean_R2, save_csv_mean, quote = FALSE)
}


# horizon = c(10:110)
horizon = c(5:110)
Fun_ROOS(Index, 'MKVola', select, horizon)
Fun_R2(Index, 'MKVola', select, horizon)

Fun_ROOS_fixed(Index, 'MKVola', select, horizon)
Fun_R2_fixed(Index, 'MKVola', select, horizon)

Fun_ROOS(Index = Index_MVWealth, 'MVVola', select, horizon)
Fun_ROOS_fixed(Index = Index_MVWealth, 'MVVola', select, horizon)

Fun_ROOS(Index, 'Return', select, horizon = c(10,25,63,110))
Fun_R2(Index, 'Return', select, horizon = c(10,25,63,110))
Fun_ROOS_fixed(Index, 'Return', select, horizon = c(10,25,63,110))
Fun_R2_fixed(Index, 'Return', select, horizon = c(10,25,63,110))

Fun_ROOS_fixed(Index = Test_Fin_Index, 'MKVola', select = c('MKVola_L_10','MKVola_L_25', 'MKVola_L_63', 'MKVola_L_110'), horizon = c(0,1))
Fun_R2_fixed(Index = Test_Fin_Index, 'MKVola', select = c('MKVola_L_10','MKVola_L_25','MKVola_L_63','MKVola_L_110'), horizon = c(0,1))


# Output 
# Portfolio Prediction

Fun_ROOS_Performance <- function(i_type, fix){
  if (fix == FALSE){
    R_oos_out  = data.frame(read.csv(paste0(save, "/Prediction_Continuous",i_type,"_vol_ROO_",save_date,".csv")))
  }else{
    R_oos_out  = data.frame(read.csv(paste0(save, "/Prediction_Continuous",i_type,"_Fixed63_ROO_",save_date,".csv")))
  }
 
  
  # within 10 days 1 month 5weeks 3 months 5months
  Mean_horizon = c(10, 25, 63, 110)
  Roo_mean =  matrix(nrow = length(Mean_horizon), ncol = ncol(R_oos_out)-1)
  Roo_mean[,1] = Mean_horizon
  colnames(Roo_mean) = c('Horizon',colnames(R_oos_out)[3:ncol(R_oos_out)]) 
  # for (i in c(1:nrow(Roo_mean))){
  #   Roo_mean[i,-1] = colMeans(R_oos_out[2:which(R_oos_out$lagDay == Mean_horizon[i]),3:ncol(R_oos_out)])
  #   
  # }
  for (i in c(1:nrow(Roo_mean))){
    Roo_mean[i,-1] = as.matrix(R_oos_out[R_oos_out$lagDay == Mean_horizon[i],3:ncol(R_oos_out)])

  }
  
  if (fix == FALSE){
    write.csv(Roo_mean, paste0(save, "/R2/Table_Output_",i_type,"_vol_ROO_",save_date,".csv"), quote = FALSE)
  }else{
    write.csv(Roo_mean, paste0(save, "/R2/Table_Output_",i_type,"_Fixed63_ROO_",save_date,".csv"), quote = FALSE)
  }
  
  
}

Fun_R2_Performance <- function(i_type, fix){
  
  if (fix == FALSE){
    R_oos_out  =  data.frame(read.csv(paste0(save, "/Prediction_Continuous",i_type,"_R2_",save_date,".csv")))
  }else{
    R_oos_out  = data.frame(read.csv(paste0(save, "/Prediction_Continuous",i_type,"_Fixed63_R2_",save_date,".csv")))
  }

  # within 10 days 1 month 5weeks 3 months 5months
  Mean_horizon = c(10, 25, 63, 110)
  Roo_mean =  matrix(nrow = length(Mean_horizon), ncol = ncol(R_oos_out)-1)
  Roo_mean[,1] = Mean_horizon
  colnames(Roo_mean) = c('Horizon',colnames(R_oos_out)[3:ncol(R_oos_out)]) 
  # for (i in c(1:nrow(Roo_mean))){
  #   Roo_mean[i,-1] = colMeans(R_oos_out[2:which(R_oos_out$lagDay == Mean_horizon[i]),3:ncol(R_oos_out)])
  #   
  # }
  for (i in c(1:nrow(Roo_mean))){
    Roo_mean[i,-1] = as.matrix(R_oos_out[R_oos_out$lagDay == Mean_horizon[i],3:ncol(R_oos_out)])
    
  }
  
  if (fix == FALSE){
    write.csv(Roo_mean, paste0(save, "/R2/Table_Output_",i_type,"_R2_",save_date,".csv"), quote = FALSE)
  }else{
    write.csv(Roo_mean, paste0(save, "/R2/Table_Output_",i_type,"_Fixed63_R2_",save_date,".csv"), quote = FALSE)
  }
}

Fun_ROOS_Performance('MVVola',fix = FALSE)
Fun_ROOS_Performance('MVVola',fix = TRUE)

Fun_ROOS_Performance('MKVola',fix = FALSE)
Fun_ROOS_Performance('MKVola',fix = TRUE)

Fun_R2_Performance('MKVola',fix = TRUE)
Fun_ROOS_Performance('Return',fix = TRUE)
Fun_R2_Performance('Return',fix = TRUE)


# ML prediction SVR, rf 
# SHAP value 
# prediction error: MAE RMSE
# https://github.com/liuyanguu/SHAPforxgboost

Fun_ROOS_ML <-function(iLag, Index,i_type){
  R_oos = data.frame(matrix(0,1,4))
  colnames(R_oos) = c('Risk','R_oos_xgboost','MAE_xgboost','RMSE_xgboost')

  if (i_type == 'Return'){
    Testdata = Index
  }else if(i_type == 'MKVola'){
    Testdata = Index
    Testdata$MKVola = NA
    Testdata$MKVola[iLag:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, iLag)
    # Testdata$MKVola[iLag:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, 21)
  }else{
    Testdata = Index
    Testdata$MVVola = NA
    Testdata$MVVola[iLag:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return_MV, iLag)
  }
  
  
    temp_reg = cbind(Testdata[(iLag+1):nrow(Index), i_type], Testdata[1:(nrow(Index)-iLag), select])
    temp_reg = temp_reg[complete.cases(temp_reg),]
    
    colnames(temp_reg)[1:8] = c('y','x1','x2','x3','x4','x5','x6','x7')

    for (i in c(63:(nrow(temp_reg)-1))){
      
      trainval = temp_reg[(i-62):(i-1),c('x1','x2','x3','x4','x5','x6','x7')]
      trainval = sapply(trainval, as.numeric)
      trainval = as.matrix(trainval) 
      
      xgboost_model <- xgboost(data = trainval,
                               label = as.matrix(temp_reg[(i-63+1):(i-1), 'y']),
                               nround = 20,
                               verbose = FALSE)
      prediction_xgboost <- predict(xgboost_model, newdata = as.matrix(temp_reg[i, c('x1','x2','x3','x4','x5','x6','x7')]))

      temp_reg[i, 'PreXgboost'] = prediction_xgboost
      temp_reg[i, 'diffXgboost'] = temp_reg$y[i]-temp_reg[i,'PreXgboost']
      temp_reg[i, 'diffXgboost_s'] = (temp_reg[i, 'diffXgboost'])^2
      
      temp_reg[i, 'diff_m'] = (temp_reg$y[i]-mean(temp_reg$y[63:nrow(temp_reg)]))^2
    }
    
    R_oos[2] = 1 - sum(temp_reg[,'diffXgboost_s'], na.rm = TRUE)/sum(temp_reg[,'diff_m'], na.rm = TRUE)
    R_oos[3] = mean(abs(temp_reg[,'diffXgboost']),na.rm = TRUE)
    R_oos[4] = sqrt(mean(temp_reg[,'diffXgboost_s'],na.rm = TRUE))

    write.csv(R_oos, paste0(save, "/Shapley/MLPrediction_",i_type,"_lag_",iLag,"_vol_ROO_",save_date,".csv"), quote = FALSE)

}

Fun_ROOS_fixed_ML <-function(iLag, Index,i_type){
  R_oos = data.frame(matrix(0,1,4))
  colnames(R_oos) = c('Risk','R_oos_xgboost','MAE_xgboost','RMSE_xgboost')
  
  if (i_type == 'Return'){
    Testdata = Index
  }else if(i_type == 'MKVola'){
    Testdata = Index
    Testdata$MKVola = NA
    Testdata$MKVola[63:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, 63)
    # Testdata$MKVola[iLag:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, 21)
  }else{
    Testdata = Index
    Testdata$MVVola = NA
    Testdata$MVVola[63:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return_MV, 63)
  }
  
  
  temp_reg = cbind(Testdata[(iLag+1):nrow(Index), i_type], Testdata[1:(nrow(Index)-iLag), select])
  temp_reg = temp_reg[complete.cases(temp_reg),]
  
  colnames(temp_reg)[1:8] = c('y','x1','x2','x3','x4','x5','x6','x7')
  
  
  
  for (i in c(63:(nrow(temp_reg)-1))){
    
    trainval = temp_reg[(i-62):(i-1),c('x1','x2','x3','x4','x5','x6','x7')]
    trainval = sapply(trainval, as.numeric)
    trainval = as.matrix(trainval) 
    
    xgboost_model <- xgboost(data = trainval,
                             label = as.matrix(temp_reg[(i-63+1):(i-1), 'y']),
                             nround = 20,
                             verbose = FALSE)
    prediction_xgboost <- predict(xgboost_model, newdata = as.matrix(temp_reg[i, c('x1','x2','x3','x4','x5','x6','x7')]))
    
    temp_reg[i, 'PreXgboost'] = prediction_xgboost
    temp_reg[i, 'diffXgboost'] = temp_reg$y[i]-temp_reg[i,'PreXgboost']
    temp_reg[i, 'diffXgboost_s'] = (temp_reg[i, 'diffXgboost'])^2
    
    temp_reg[i, 'diff_m'] = (temp_reg$y[i]-mean(temp_reg$y[63:nrow(temp_reg)]))^2
  }
  
  R_oos[2] = 1 - sum(temp_reg[,'diffXgboost_s'], na.rm = TRUE)/sum(temp_reg[,'diff_m'], na.rm = TRUE)
  R_oos[3] = mean(abs(temp_reg[,'diffXgboost']),na.rm = TRUE)
  R_oos[4] = sqrt(mean(temp_reg[,'diffXgboost_s'],na.rm = TRUE))
  
  write.csv(R_oos, paste0(save, "/Shapley/MLPrediction_",i_type,"_lag_",iLag,"_Fixed63_ROO_",save_date,".csv"), quote = FALSE)
  
}

Fun_Shapley_ML <-function(iLag, Index,i_type){
  
  if (i_type == 'Return'){
    Testdata = Index
  }else if(i_type == 'MKVola'){
    Testdata = Index
    Testdata$MKVola = NA
    Testdata$MKVola[iLag:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, iLag)
    # Testdata$MKVola[iLag:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, 21)
  }else{
    Testdata = Index
    Testdata$MVVola = NA
    Testdata$MVVola[iLag:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return_MV, iLag)
  }
  
  
  temp_reg = cbind(Testdata[(iLag+1):nrow(Index), i_type], Testdata[1:(nrow(Index)-iLag), select])
  colnames(temp_reg)[1] = i_type
  temp_reg = temp_reg[complete.cases(temp_reg),]
  trainval = as.matrix(temp_reg[,-1])
  
  
  param_list <- list(objective = "reg:squarederror",  # For regression
                     eta = 0.02,
                     max_depth = 100,
                     gamma = 0.01,
                     subsample = 0.8,
                     colsample_bytree = 0.86)
 
  for (i in c(1:500)){
    mod <- xgboost::xgboost(data = trainval, 
                            label = as.matrix(temp_reg[, i_type]), 
                            params = param_list, nrounds = 200,
                            verbose = FALSE,
                            early_stopping_rounds = 8)
    shap_long <- shap.prep(xgb_model = mod, X_train = trainval)
    Shapley_Mean_Temp = unique(shap_long[,c('variable','mean_value')])
    colnames(Shapley_Mean_Temp)[2] = paste0('mean_value_',i)
    if (i == 1){
      Shapley_Mean = Shapley_Mean_Temp
    }else{
      Shapley_Mean = merge(Shapley_Mean, Shapley_Mean_Temp)
    }
  }
  Shapley_Mean$mean_value = rowMeans(Shapley_Mean[,-1])
  Shapley_Mean$variable = as.character(Shapley_Mean$variable)
  Shapley_Mean$variable[Shapley_Mean$variable == 'FRM_5P'] = 'FRM 5%'
  Shapley_Mean$variable[Shapley_Mean$variable == 'FRM_25P'] = 'FRM 25%'
  Shapley_Mean$variable[Shapley_Mean$variable == 'FRM_50P'] = 'FRM 50%'
  Shapley_Mean$variable[Shapley_Mean$variable == 'TCI'] = 'Total Connectedness'
  Shapley_Mean$variable[Shapley_Mean$variable == 'DAG'] = 'Bayesian Graphical VAR'
  Shapley_Mean$variable[Shapley_Mean$variable == 'PCA_1'] = 'Principal Components'
  Shapley_Mean$variable[Shapley_Mean$variable == 'GDC'] = 'Granger Causality'
  
  Shapley_Mean$Shapley_Mean = round(Shapley_Mean$mean_value,5)
  
  p <- ggplot(Shapley_Mean,aes(x= reorder(variable, Shapley_Mean),y=Shapley_Mean, fill = variable))+geom_bar(stat="identity", width=0.5)+
    coord_flip()+
    geom_text(aes(label=Shapley_Mean), vjust=-2, size=5)+
    labs(y = "Shapley Value", size = 20)+
    labs(x = "Risk Measures", size = 20)+
    theme(axis.text.x = element_text(size =15, angle = 30, hjust = 1), 
          axis.text.y = element_text(size = 15), 
          axis.title =  element_text(size=,face = "bold"),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 35, face = "bold"),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          axis.line = element_line(colour = "black"),legend.position="none")
  
  # png(paste0(save, "/Shapley/SHAP_",i_type,"_lag_",iLag,".png"), width = 900, height = 600, bg = "transparent")
  # 
  # print(p)
  # dev.off()
  
  ggsave(paste0(save, "/Shapley/SHAP_",i_type,"_lag_",iLag,".eps"), plot = p, width = 12, height = 8, bg = "transparent",dpi = 800)
  
  # write.csv(Shapley_Mean, paste0(save, "/Shapley/ShapleyCalData_",i_type,"_lag_",iLag,"_",save_date,".csv"), quote = FALSE)
  
}

Fun_Shapley_fixed_ML <-function(iLag, Index,i_type){
  
  if (i_type == 'Return'){
    Testdata = Index
  }else if(i_type == 'MKVola'){
    Testdata = Index
    Testdata$MKVola = NA
    Testdata$MKVola[63:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, 63)
    # Testdata$MKVola[iLag:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return, 21)
  }else{
    Testdata = Index
    Testdata$MVVola = NA
    Testdata$MVVola[63:nrow(Testdata)] <- Fun_Vola_dynamic(Testdata$Return_MV, 63)
  }
  
  
  temp_reg = cbind(Testdata[(iLag+1):nrow(Index), i_type], Testdata[1:(nrow(Index)-iLag), select])
  colnames(temp_reg)[1] = i_type
  temp_reg = temp_reg[complete.cases(temp_reg),]
  trainval = as.matrix(temp_reg[,-1])
  
  
  param_list <- list(objective = "reg:squarederror",  # For regression
                     eta = 0.02,
                     max_depth = 100,
                     gamma = 0.01,
                     subsample = 0.8,
                     colsample_bytree = 0.86)
  
  for (i in c(1:500)){
    mod <- xgboost::xgboost(data = trainval, 
                            label = as.matrix(temp_reg[, i_type]), 
                            params = param_list, nrounds = 200,
                            verbose = FALSE,
                            early_stopping_rounds = 8)
    shap_long <- shap.prep(xgb_model = mod, X_train = trainval)
    Shapley_Mean_Temp = unique(shap_long[,c('variable','mean_value')])
    colnames(Shapley_Mean_Temp)[2] = paste0('mean_value_',i)
    if (i == 1){
      Shapley_Mean = Shapley_Mean_Temp
    }else{
      Shapley_Mean = merge(Shapley_Mean, Shapley_Mean_Temp)
    }
  }
  Shapley_Mean$mean_value = rowMeans(Shapley_Mean[,-1])
  Shapley_Mean$variable = as.character(Shapley_Mean$variable)
  Shapley_Mean$variable[Shapley_Mean$variable == 'FRM_5P'] = 'FRM 5%'
  Shapley_Mean$variable[Shapley_Mean$variable == 'FRM_25P'] = 'FRM 25%'
  Shapley_Mean$variable[Shapley_Mean$variable == 'FRM_50P'] = 'FRM 50%'
  Shapley_Mean$variable[Shapley_Mean$variable == 'TCI'] = 'Total Connectedness'
  Shapley_Mean$variable[Shapley_Mean$variable == 'DAG'] = 'Bayesian Graphical VAR'
  Shapley_Mean$variable[Shapley_Mean$variable == 'PCA_1'] = 'Principal Components'
  Shapley_Mean$variable[Shapley_Mean$variable == 'GDC'] = 'Granger Causality'
  
  Shapley_Mean$Shapley_Mean = round(Shapley_Mean$mean_value,5)
  
  # png(paste0(save, "/Shapley/SHAP_fixed63_",i_type,"_lag_",iLag,".png"), width = 900, height = 600, bg = "transparent")
  # 
  # print(ggplot(Shapley_Mean,aes(x= reorder(variable, Shapley_Mean),y=Shapley_Mean, fill = variable))+geom_bar(stat="identity", width=0.5)+
  #         coord_flip()+
  #         geom_text(aes(label=Shapley_Mean), vjust=-2, size=5)+
  #         labs(y = "Shapley Value", size = 20)+
  #         labs(x = "Risk Measures", size = 20)+
  #         theme(axis.text.x = element_text(size =15, angle = 30, hjust = 1), 
  #               axis.text.y = element_text(size = 15), 
  #               axis.title =  element_text(size=,face = "bold"),
  #               panel.grid.major =element_blank(), 
  #               panel.grid.minor = element_blank(),
  #               plot.title = element_text(hjust = 0.5,size = 35, face = "bold"),
  #               panel.background = element_rect(fill = "transparent",colour = NA),
  #               plot.background = element_rect(fill = "transparent",colour = NA),
  #               axis.line = element_line(colour = "black"),legend.position="none")
  # )
  # dev.off()
  
  # Shapley_Mean <- read.csv(paste0(save, "/Shapley/ShapleyCalData_fixed63_",i_type,"_lag_",iLag,"_",save_date,".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  # 
  p = ggplot(Shapley_Mean,aes(x= reorder(variable, Shapley_Mean),y=Shapley_Mean, fill = variable))+geom_bar(stat="identity", width=0.5)+
          coord_flip()+
          geom_text(aes(label=Shapley_Mean), vjust=-2, size=5)+
          labs(y = "Shapley Value", size = 20)+
          labs(x = "Risk Measures", size = 20)+
          theme(axis.text.x = element_text(size =15, angle = 30, hjust = 1),
                axis.text.y = element_text(size = 15),
                axis.title =  element_text(size=,face = "bold"),
                panel.grid.major =element_blank(),
                panel.grid.minor = element_blank(),
                plot.title = element_text(hjust = 0.5,size = 35, face = "bold"),
                panel.background = element_rect(fill = "transparent",colour = NA),
                plot.background = element_rect(fill = "transparent",colour = NA),
                axis.line = element_line(colour = "black"),legend.position="none")
  
  ggsave(paste0(save, "/Shapley/SHAP_",i_type,"_fixed63_lag_",iLag,".eps"), plot = p, width = 12, height = 8, bg = "transparent",dpi = 800)
  
  # write.csv(Shapley_Mean, paste0(save, "/Shapley/ShapleyCalData_fixed63_",i_type,"_lag_",iLag,"_",save_date,".csv"), quote = FALSE)
  
}

Fun_ROOS_fixed_ML(iLag = 10, Index, 'MKVola')
Fun_Shapley_fixed_ML(iLag = 10, Index, 'MKVola')

Fun_ROOS_ML(iLag = 10, Index, 'MKVola')

Fun_Shapley_ML(iLag = 10, Index, 'MKVola')
Fun_Shapley_ML(iLag = 25, Index, 'MKVola')
Fun_Shapley_ML(iLag = 63, Index, 'MKVola')
Fun_Shapley_ML(iLag = 110, Index, 'MKVola')

Fun_Shapley_ML(iLag = 10, Index_MVWealth, 'MVVola')
Fun_Shapley_ML(iLag = 25, Index_MVWealth, 'MVVola')
Fun_Shapley_ML(iLag = 63, Index_MVWealth, 'MVVola')
Fun_Shapley_ML(iLag = 110, Index_MVWealth, 'MVVola')

  
# Generate beta test part 
Index$MKVola_10[10:nrow(Index)] <- Fun_Vola_dynamic(Index$Return, iLag = 10)
Index$MKVola_25[25:nrow(Index)] <- Fun_Vola_dynamic(Index$Return, iLag = 25)
Index$MKVola_63[63:nrow(Index)] <- Fun_Vola_dynamic(Index$Return, iLag = 63)
Index$MKVola_110[110:nrow(Index)] <- Fun_Vola_dynamic(Index$Return, iLag = 110)

fun_lag <- function(x,lag){
  return(paste0(x,'_L',lag))
} 

pre_col_L10 = fun_lag(select, 10)
pre_col_L25 = fun_lag(select, 25)
pre_col_L63 = fun_lag(select, 63)
pre_col_L110 = fun_lag(select, 110)

pre_col = c(pre_col_L10, pre_col_L25, pre_col_L110)
Index[, pre_col] = NA

Index[11:nrow(Index), pre_col_L10] = Index[1:(nrow(Index)-10), select]
Index[26:nrow(Index), pre_col_L25] = Index[1:(nrow(Index)-25), select]
Index[64:nrow(Index), pre_col_L63] = Index[1:(nrow(Index)-63), select]
Index[111:nrow(Index), pre_col_L110] = Index[1:(nrow(Index)-110), select]


write_dta(Index, paste0(save, "/InsampleBetaTest_",save_date,".dta"))


# plot R2
save_csv =  paste0(save, "/Back_otherPrediction/Prediction_ContinuousMKVola_Fixed63_ROO_",save_date,".csv")
R_oos_out = data.frame(read.csv(save_csv))
# png(paste0(save,  "/R2/Roos_Dynamic_MKVola.png"), width = 900, height = 600, bg = "transparent")
# Plot https://www.cnblogs.com/biostat-yu/p/13839621.html

# R_oos_out = data.frame(R_oos_out)
p=ggplot(R_oos_out, aes(x=lagDay)) +
  geom_line(aes(y = Roos_FRM_5P), color="#009933", linewidth = 0.75) +
  geom_line(aes(y = Roos_FRM_25P), color="#996600", linewidth = 0.75) + 
  geom_line(aes(y = Roos_FRM_50P), color="#33CCCC", linewidth = 0.75) +
  geom_line(aes(y = Roos_DAG), color="#FF6666", linewidth = 0.75) +
  geom_line(aes(y = Roos_GDC), color="#666666", linewidth = 0.75) +
  geom_line(aes(y = Roos_PCA_1), color="#9966CC", linewidth = 0.75) +
  geom_line(aes(y = Roos_TCI), color="magenta", linewidth = 0.75) +
  scale_x_continuous(limits = c(5, 110), breaks = seq(10,110,by=10))+
  scale_y_continuous(limits = c(0.5, 0.8), breaks = seq(-0.4,0.8,by=0.1))+
  labs(y = "Out-of-Sample R2", size = 20)+
  labs(x = "Prediction Horizon (day)", size = 20)+
  theme( panel.grid=element_blank(),
         legend.position = "none",
         axis.title =  element_text(size=,face = "bold"),
         panel.background = element_rect(fill = "transparent",colour = NA),
         plot.background = element_rect(fill = "transparent",colour = NA),
         legend.box.background = element_rect(fill = "transparent"),
         axis.line = element_line(colour = "black"),
         axis.text.x = element_text(size=14),
         axis.text.y = element_text(size=14 ) )
# print(p)
# dev.off()
ggsave(paste0(save,  "/R2/Roos_Dynamic_MKVola.eps"), plot = p, width = 12, height = 8, bg = "transparent", dpi = 800)

# Plot beta's trend
Beta = read_excel( paste0(save, "/InsampleBeta/InsampleBeta_dynamic.xls")) %>% data.frame()
colnames(Beta) = c('Year','Beta','lb','ub')
Beta = Beta[Beta$Year<= 2021, ]
# png(paste0(save, "/InsampleBeta/InsamplePredictFRM_5P.png"), 
#     width = 900, height = 600, bg = "transparent")

p <- ggplot(Beta, aes(x = Year)) +
  geom_segment(aes(xend = Year, y = lb, yend = ub), linewidth = 0.75) +
  geom_point(aes(y = Beta), size = 1.5) +
  geom_line(aes(y = Beta), linewidth = 0.75) +
  scale_x_continuous(breaks = unique(Beta$Year), labels = unique(Beta$Year), expand = c(0.1, 0.11)) +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),axis.line = element_line(colour = "black"),axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))+
  labs(x = "Year",
       y = "Coefficient") 
# dev.off()

ggsave(paste0(save, "/InsampleBeta/InsamplePredictFRM_5P.eps"), plot = p, width = 12, height = 8, bg = "transparent", dpi = 800)



#  predict vcrix
x_type = c('vcrix')

Index_crix = merge(crix, Index, by = intersect(names(crix), names(Index)),all.x = TRUE,sort = TRUE)
for (i_type in x_type){
  for (j in pre_col){
    temp_reg = Index_crix[,c(i_type,j)]
    temp_reg = temp_reg[!is.na(temp_reg[,j]),]
    colnames(temp_reg) = c('y','x')
    for (i in c(63:(nrow(temp_reg)-1))){
      
      fit = lm(y ~ x, temp_reg[(i-62):(i-1), ])
      
      temp_reg[i, paste0('p_vol_',j)] = fit$coefficients[1] + fit$coefficients[2]*temp_reg$x[i]
      temp_reg[i, paste0('diff_s_',j)] = (temp_reg$y[i]-temp_reg[i,paste0('p_vol_',j)])^2
      # temp_reg[i, paste0('diff_m_',j)] = (temp_reg$y[i]-mean(temp_reg$y))^2
      temp_reg[i, paste0('diff_m_',j)] = (temp_reg$y[i]-mean(temp_reg$y[63:nrow(temp_reg)]))^2
    }
    
    R_oos[1,j] = 1 - sum(temp_reg[,paste0('diff_s_',j)], na.rm = TRUE)/sum(temp_reg[,paste0('diff_m_',j)], na.rm = TRUE)
    R_oos[2,j] =sum(temp_reg[,paste0('diff_s_',j)], na.rm = TRUE)
    R_oos[3,j] =sum(temp_reg[,paste0('diff_m_',j)], na.rm = TRUE)
  }
  write.csv(R_oos, paste0(wdir, "/Output/Crypto/Prediction_",i_type,"_vol_ROO_",save_date,".csv"), quote = FALSE)
  
  
  R_2 = data.frame(matrix(0,3,length(pre_col)))
  colnames(R_2) = pre_col
  for (j in pre_col){
    temp_reg = Index_crix[,c(i_type,j)]
    temp_reg = temp_reg[!is.na(temp_reg[,j]),]
    colnames(temp_reg) = c('y','x')
    for (i in c(63:(nrow(temp_reg)-1))){
      
      fit = lm(y ~ x, temp_reg[(i-61):i, ])
      
      temp_reg[i, paste0('p_vol_',j)] = fit$coefficients[1] + fit$coefficients[2]*temp_reg$x[i]
      temp_reg[i, paste0('diff_s_',j)] = (temp_reg$y[i]-temp_reg[i,paste0('p_vol_',j)])^2
      temp_reg[i, paste0('diff_m_',j)] = (temp_reg$y[i]-mean(temp_reg$y[63:nrow(temp_reg)]))^2
      # temp_reg[i, paste0('diff_m_',j)] = (temp_reg$y[i]-mean(temp_reg$y[i:nrow(temp_reg)]))^2
    }
    
    R_2[1,j] = 1 - sum(temp_reg[,paste0('diff_s_',j)], na.rm = TRUE)/sum(temp_reg[,paste0('diff_m_',j)], na.rm = TRUE)
    R_2[2,j] =sum(temp_reg[,paste0('diff_s_',j)], na.rm = TRUE)
    R_2[3,j] =sum(temp_reg[,paste0('diff_m_',j)], na.rm = TRUE)
  }
  write.csv(R_2, paste0(wdir, "/Output/Crypto/Prediction_",i_type,"_vol_R2_",save_date,".csv"), quote = FALSE)
}


