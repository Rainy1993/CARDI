rm(list = ls(all = TRUE))

# wdir = "//Users/ruting/Documents/Github/FRM/FRM_Quantlet/FRM_All"
wdir = "/Users/ruting/Documents/macbook/PcBack/23.Topic_wkh_cooperation/FRM/"

#Check if package is installed, if not: install, either way: load 
libraries = c("iml","ggplot2", "data.table", "igraph","timeDate", "stringr", "graphics","magick", "scales", "tidyr", "zoo", "foreach", "doParallel",
              "xgboost","shapr","randomForest", "rpart", "quantreg", "readxl","dplyr", "xlsx", "psych","qgraph", "gganimate","av",
              "gifski", "strex","matrixStats","tools","qgraph","haven")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

save_date = 20251015
save =  paste0(wdir, "/Output/EfficiencyTest/", save_date)

channel_High = "HighCarbonIntens"
channel_Low = "LowCarbonIntens"

CarbonPremia_Month = read_dta("Output/EfficiencyTest/InsampleBeta_Premia_1014_month.dta")
CarbonPremia_Month = data.frame(CarbonPremia_Month)
CarbonPremia_Month$Date = as.character(CarbonPremia_Month$Date)

cVaRlist = c("QrVaR_1", "QrVaR_5", "QrVaR_10")
CarbonPremia_Month[,cVaRlist] = CarbonPremia_Month[,cVaRlist]*100

CarbonPremia_Month = CarbonPremia_Month[,c("Date",c("Return", cVaRlist))]

s = 63

M_macro = 9
setwd(wdir)

Degree_Cal <- function(input_path){
  
  #Create a list of files in the folder and extract dates from the names
  file_list = list.files(path = paste0(input_path, "/Adj_Matrices"))
  file_list = file_list[file_list!="Fixed"]
  dates = as.character(str_first_number(file_list), format = "%Y%m%d")
  dates = as.Date(dates, format = "%Y%m%d")
  N = length(file_list)
  
  #Create a list of all network graphs
  allgraphs = lapply(1:N, function(i) {
    data = read.csv(paste0(input_path, "/Adj_Matrices/", file_list[i]), row.names = 1)
    M_stock = ncol(data)-M_macro
    adj_matrix = data.matrix(data[1:M_stock, 1:M_stock])
    q = qgraph(adj_matrix, layout = "circle", details = TRUE, 
               vsize = c(5,15), DoNotPlot = TRUE)
    return(q)
  })
  
  allcentralities = centrality(allgraphs)
  
  compute_eigencentrality_list <- function(graph_list, use_weights = FALSE) {
    N <- length(graph_list)
    eigencentrality_list <- vector("list", N)
    log_msgs <- character(0)
    
    for (i in seq_len(N)) {
      g <- tryCatch(as.igraph(graph_list[[i]]), error = function(e) NULL)
      
      if (is.null(g)) {
        log_msgs <- c(log_msgs, paste("❌ Graph", i, ": Conversion to igraph failed."))
        eigencentrality_list[[i]] <- NA
        next
      }
      
      w <- if (use_weights) E(g)$weight else NA
      
      ec <- tryCatch({
        eigen_centrality(g, weights = w)$vector
      }, error = function(e) {
        log_msgs <- c(log_msgs, paste("❌ Graph", i, ": Eigen centrality failed -", e$message))
        rep(NA, vcount(g))
      })
      
      if (length(ec) != vcount(g)) {
        log_msgs <- c(log_msgs, paste("⚠️ Graph", i, ": Centrality length mismatch. Expected", vcount(g), "but got", length(ec)))
        ec <- rep(NA, vcount(g))
      }
      
      eigencentrality_list[[i]] <- ec
    }
    
    
    return(eigencentrality_list)
  }
  
  eigencentrality <- compute_eigencentrality_list(allgraphs, use_weights = TRUE)
  
  
  outdegree_avg = sapply(1:N, function(i) mean(allcentralities[[i]]$OutDegree))
  indegree_avg = sapply(1:N, function(i) mean(allcentralities[[i]]$InDegree))
  closeness_avg = sapply(1:N, function(i) mean(allcentralities[[i]]$Closeness))
  betweenness_avg = sapply(1:N, function(i) mean(allcentralities[[i]]$Betweenness))
  eigenvector_avg = sapply(1:N, function(i) mean(eigencentrality[[i]]))
  
  # Degree_Index = cbind(FRM_Reg, outdegree_avg, indegree_avg, closeness_avg, eigenvector_avg)
  Degree_Index = cbind(FRM_Reg[,1],outdegree_avg, indegree_avg, closeness_avg, eigenvector_avg)
  
  return(Degree_Index)
}

tau_list = c(0.05,0.01,0.1)
for (tau in tau_list){
  
  if (tau == 0.05) input_path_High = paste0("Output/", channel_High) else 
    input_path_High = paste0("Output/", channel_High, "/Sensitivity/tau=", 100*tau, "/s=", s)
  output_path_High = input_path_High
  
  if (tau == 0.05) input_path_Low = paste0("Output/", channel_Low) else 
    input_path_Low = paste0("Output/", channel_Low, "/Sensitivity/tau=", 100*tau, "/s=", s)
  output_path_Low = input_path_Low
  
  #Set ggplot2 theme
  theme_set(theme_classic())
  
  #List of centrality types and their numbers
  centralitylist = list("OutDegree" = 1, "InDegree" = 2, "Closeness" = 3, 
                        "Betweenness" = 4, "InInfluence" = 5, "OutInfluence" = 6)
  
  #Read historical carbon indicator 
  # FRM_index = read.csv(paste0(input_path, "/Lambda/FRM_", channel, "_index.csv"))
  
  
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
  # FRM_Reg$Date = as.Date(FRM_Reg$Date)
  ## Calculate centralities
  Degree_Index_High = Degree_Cal(input_path_High)
  Degree_Index_Low = Degree_Cal(input_path_Low)
  
  colnames(Degree_Index_High)[1] = "Date"
  colnames(Degree_Index_High)[-1] = paste0(colnames(Degree_Index_High)[-1],"_High")
  
  colnames(Degree_Index_Low)[1] = "Date"
  colnames(Degree_Index_Low)[-1] = paste0(colnames(Degree_Index_Low)[-1],"_Low")
  
  Degree_Index = cbind(Degree_Index_High, Degree_Index_Low[,-1])
  Degree_Index = data.frame(Degree_Index)
  
  Degree_Index = left_join(Degree_Index,FRM_Reg[,c("Date","FRM_5_High_Low","FRM_1_High_Low","FRM_10_High_Low")],by = "Date" )
  
  Degree_Index[,-1] = lapply(Degree_Index[,-1],as.numeric)
  Degree_Index$Degree_High_Low = Degree_Index$indegree_avg_High / Degree_Index$indegree_avg_Low
  Degree_Index$Eigen_High_Low = Degree_Index$eigenvector_avg_High / Degree_Index$eigenvector_avg_Low
  
  Degree_Index = left_join(Degree_Index,CarbonPremia_Month, by = "Date")
  
  # corelationship
  lag = c(0, 63, 84, 105)
  
  for (iLag in lag){
    
    if (tau == 0.05){
      Lag_Index = Degree_Index[c(1:(nrow(Degree_Index)-iLag)),c('FRM_5_High_Low','Degree_High_Low','Eigen_High_Low')]
      
    }
    
    if (tau == 0.01){
      Lag_Index = Degree_Index[c(1:(nrow(Degree_Index)-iLag)),c('FRM_1_High_Low','Degree_High_Low','Eigen_High_Low')]
      
    }
    
    if (tau == 0.1){
      Lag_Index = Degree_Index[c(1:(nrow(Degree_Index)-iLag)),c('FRM_10_High_Low','Degree_High_Low','Eigen_High_Low')]
      
    }
   
    cormatrix = cbind(Degree_Index[(iLag+1):nrow(Degree_Index),'Return'],Lag_Index)
    cormatrix <- cormatrix[, sapply(cormatrix, is.numeric)]
    
    # 删除含 NA / NaN / Inf 的行（安全处理）
    cormatrix_clean <- cormatrix[apply(cormatrix, 1, function(row) all(is.finite(row))), ]
    colnames(cormatrix_clean)[1] = 'LowCarbonPremium'
    colnames(cormatrix_clean)[2] = 'Cardi'
    
    corr = cor(cormatrix_clean,  method = "pearson")
    corr
    
    
    p_value_matrix <- matrix(NA, nrow = ncol(corr), ncol = ncol(corr))
    
    # Calculate p-values for each pair of variables
    for (i in 1:(ncol(cormatrix_clean) - 1)) {
      for (j in (i + 1):ncol(cormatrix_clean)) {
        cor_test_result <- cor.test(cormatrix_clean[, i], cormatrix_clean[, j], method = "pearson")
        p_value_matrix[i, j] <- p_value_matrix[j, i] <- cor_test_result$p.value
      }
    }
    p_value_matrix
    
    write.csv(corr,paste0(wdir, "Output/Centrality/Centralilty_Corr_",iLag,"_tau",100*tau,".csv"),quote = FALSE)
    write.csv(p_value_matrix,paste0(wdir, "Output/Centrality/Centralilty_CorrP_",iLag,"_tau",100*tau,".csv"),quote = FALSE)
    
  }
  
  # add macro
  macro = read.csv(file = paste0('Input/HighCarbonIntens/20140704-20250127/HighCarbonIntens_Macro_20250127.csv'), header = TRUE) %>% data.frame()
  Degree_Index = left_join(Degree_Index,macro,by = "Date")
  
  write_dta(Degree_Index, paste0(save,"/Insample_Centrality_tau",100*tau,".dta"))
  
}



# backup: VAR test
library(vars)
iLag = 63
Lag_Index = Degree_Index[c(1:(nrow(Degree_Index)-iLag)),c('FRM_5_High_Low','Degree_High_Low','Eigen_High_Low')]

cormatrix = cbind(Degree_Index[(iLag+1):nrow(Degree_Index),'Return'],Lag_Index[,c('FRM_5_High_Low','Eigen_High_Low','Degree_High_Low')])
cormatrix <- cormatrix[, sapply(cormatrix, is.numeric)]

# 删除含 NA / NaN / Inf 的行（安全处理）
cormatrix_clean <- cormatrix[apply(cormatrix, 1, function(row) all(is.finite(row))), ]
# cormatrix_clean = cormatrix_clean[cormatrix_clean$Close_High_Low!= 0,]
colnames(cormatrix_clean)[1] = 'LowCarbonPremium'
colnames(cormatrix_clean)[2] = 'Cardi'


# Step 2: 转换为时间序列对象（假设季度数据）
ts_data <- ts(cormatrix_clean, frequency = 4)

# Step 3: 滞后阶选择（建议清洗后再做）
library(vars)
lagselect <- VARselect(ts_data, lag.max = 5, type = "const")
print(lagselect$selection)

# Step 4: 构建 VAR 模型（假设选择了 p = 2）
var_model <- VAR(ts_data, p = 5, type = "const")
summary(var_model)

# Sigma <- cov(ts_data)
# kappa(Sigma)  # 条件数，越大说明越不稳定（>10^10就需要警惕）

# 因果检验
causality(var_model, cause = "Degree_High_Low")
causality(var_model, cause = "Cardi")

# IRF
plot(irf(var_model, impulse = "Degree_High_Low", response = c("Cardi", "LowCarbonPremium"), n.ahead = 5, boot = TRUE))

# FEVD
plot(fevd(var_model, n.ahead = 10))


cent_plot = function(cent_type, lambda) {
  # cent_string = deparse(substitute(cent_type))
  cent_string = cent_type
  lambda = Degree_Index$FRM_5_High_Low
  
  png(paste0(wdir, "Output/Centrality/FRM_", cent_string,".png"), 
      width = 900, height = 600, bg = "transparent")
  par(mar = c(5, 4, 4, 4) + 0.3)
  # plot(lambda[-outliers], type = "l", col = "blue", xlab = "", 
  #      ylab = "FRM index", xaxt = "n", lwd = 2)
  plot(lambda, type = "l", col = "blue", xlab = "", 
       ylab = "FRM based Carbon risk", xaxt = "n", lwd = 2)
  par(new = TRUE)
  # plot(cent_type[-outliers], type = "l", col = "red", axes = FALSE, 
  #      xlab = "", ylab = "", xaxt = "n")
  plot(Degree_Index[,cent_type], type = "l", col = "red", axes = FALSE, 
       xlab = "", ylab = "", xaxt = "n")
  # axis(side = 4, at = pretty(range(cent_type[-outliers])))
  axis(side = 4, at = pretty(range(Degree_Index[,cent_type])))
  
  # ll = which(FRM_index$date[-outliers] %in% plot_labels)
  # ll = which(FRM_index$date %in% plot_labels)
  ll = seq(120, nrow(Degree_Index), 120)
  axis(1, at = ll, labels = Degree_Index$Date[ll])
  mtext(paste0(gsub("_.*", "", cent_string) %>% toTitleCase, " centrality"), 
        side = 4, line = 3)
  dev.off()
}

# cent_plot('EDC_simp')
# cent_plot('EDH')
# cent_plot('BEKK_A')
# cent_plot('BEKK_G')
cent_plot('outdegree_avg')
cent_plot('indegree_avg')
cent_plot('betweenness_avg')
cent_plot('closeness_avg')
cent_plot('eigenvector_avg')
# cent_plot(degree_avg, FRM_index$frm)

# corr 
# cormatrix = cbind(FRM_index$frm,outdegree_avg, closeness_avg, eigenvector_avg,betweenness_avg)

cormatrix = Degree_Index[,c('FRM_5_High_Low','indegree_avg','eigenvector_avg')]
cormatrix = cormatrix[!is.na(cormatrix$eigenvector_avg),]
colnames(cormatrix)[1] = 'FRM based carbon risk'
corr = cor(cormatrix,  method = "pearson")
corr

p_value_matrix <- matrix(NA, nrow = ncol(cormatrix), ncol = ncol(cormatrix))

# Calculate p-values for each pair of variables
for (i in 1:(ncol(cormatrix) - 1)) {
  for (j in (i + 1):ncol(cormatrix)) {
    cor_test_result <- cor.test(cormatrix[, i], cormatrix[, j], method = "pearson")
    p_value_matrix[i, j] <- p_value_matrix[j, i] <- cor_test_result$p.value
  }
}

write.csv(corr,paste0(wdir, "Output/Centrality/Centralilty_FRM_Corr.csv"),quote = FALSE)
write.csv(p_value_matrix,paste0(wdir, "Output/Centrality/Centralilty_FRM_CorrP.csv"),quote = FALSE)

# compare with carbon premium



#Restructure list into individual node centralities

## Plot FRM index vs all centralities

# cent_plot = function(cent_type, lambda) {
#   cent_string = deparse(substitute(cent_type))
#   png(paste0(output_path, "/Centrality/FRM_", cent_string,".png"), 
#       width = 900, height = 600, bg = "transparent")
#   par(mar = c(5, 4, 4, 4) + 0.3)
#   # plot(lambda[-outliers], type = "l", col = "blue", xlab = "", 
#   #      ylab = "FRM index", xaxt = "n", lwd = 2)
#   plot(lambda, type = "l", col = "blue", xlab = "", 
#        ylab = "FRM index", xaxt = "n", lwd = 2)
#   par(new = TRUE)
#   # plot(cent_type[-outliers], type = "l", col = "red", axes = FALSE, 
#   #      xlab = "", ylab = "", xaxt = "n")
#   plot(cent_type, type = "l", col = "red", axes = FALSE, 
#        xlab = "", ylab = "", xaxt = "n")
#   # axis(side = 4, at = pretty(range(cent_type[-outliers])))
#   axis(side = 4, at = pretty(range(cent_type)))
#   
#   # ll = which(FRM_index$date[-outliers] %in% plot_labels)
#   # ll = which(FRM_index$date %in% plot_labels)
#   ll = seq(120, nrow(FRM_index), 120)
#   axis(1, at = ll, labels = FRM_index$date[ll])
#   mtext(paste0(gsub("_.*", "", cent_string) %>% toTitleCase, " centrality"), 
#         side = 4, line = 3)
#   dev.off()
# }
# 
# cent_plot(outdegree_avg, FRM_index$frm)
# cent_plot(indegree_avg, FRM_index$frm)
# cent_plot(betweenness_avg, FRM_index$frm)
# cent_plot(closeness_avg, FRM_index$frm)
# cent_plot(eigenvector_avg, FRM_index$frm)
# # cent_plot(degree_avg, FRM_index$frm)
# 
# # corr 
# cormatrix = cbind(FRM_index$frm,outdegree_avg, closeness_avg, eigenvector_avg,betweenness_avg)
# colnames(cormatrix)[1] = 'FRM'
# corr = cor(cormatrix,  method = "pearson")
# corr
# 
# p_value_matrix <- matrix(NA, nrow = ncol(cormatrix), ncol = ncol(cormatrix))
# 
# # Calculate p-values for each pair of variables
# for (i in 1:(ncol(cormatrix) - 1)) {
#   for (j in (i + 1):ncol(cormatrix)) {
#     cor_test_result <- cor.test(cormatrix[, i], cormatrix[, j], method = "pearson")
#     p_value_matrix[i, j] <- p_value_matrix[j, i] <- cor_test_result$p.value
#   }
# }
# 
# write.csv(corr,paste0(output_path,"/Centrality/Centralilty_FRM_Corr.csv"),quote = FALSE)
# write.csv(p_value_matrix,paste0(output_path,"/Centrality/Centralilty_FRM_CorrP.csv"),quote = FALSE)
# 

# compare with VCRIX and CRIX
VIX_file = paste0("Input/", channel, "/VCRIX")
CRIX = read.csv(paste0(VIX_file, "/CRIX.csv")) 
loc = which(CRIX$price>200000) 
CRIX$price[loc] = CRIX$price[loc] /1000
CRIX$CRIX_return = c(0,diff(log(CRIX$price)))
CRIX$vol_CRIX = NA
for (i in c(63:nrow(CRIX))){
  CRIX$vol_CRIX[i] = (sd(CRIX$CRIX_return[(i-62):i]))^2
}

FRM_index = read.csv(paste0(input_path, "/Lambda/FRM_", channel, "_index.csv"))
FRM_index = merge(FRM_index, CRIX, all.x = TRUE, all.y = FALSE)
FRM_index = FRM_index[!is.na(FRM_index$vol_CRIX),]
FRM_index = FRM_index[FRM_index$date>'2018-11-01',]

png(paste0(output_path, "/FRM_CRIX_Volatility.png"), 
    width = 900, height = 600, bg = "transparent")
par(mar = c(5, 4, 4, 4) + 0.3)
plot(FRM_index$frm, type = "l", col = "blue", xlab = "", 
     ylab = "FRM index", xaxt = "n", lwd = 2)
par(new = TRUE)
plot(FRM_index$vol_CRIX, type = "l", col = "red", axes = FALSE, 
     xlab = "", ylab = "", xaxt = "n")
axis(side = 4, at = pretty(range(FRM_index$vol_CRIX)))

ll = seq(120, nrow(FRM_index), 120)
axis(1, at = ll, labels = FRM_index$date[ll])
mtext("CRIX log return rolling variance", side = 4, line = 3)
dev.off()


# corr 
cormatrix = FRM_index[,c('frm','vol_CRIX')]
colnames(cormatrix)[1] = 'FRM'
corr = cor(cormatrix,  method = "pearson")
corr

p_value_matrix <- matrix(NA, nrow = ncol(cormatrix), ncol = ncol(cormatrix))

# Calculate p-values for each pair of variables
for (i in 1:(ncol(cormatrix) - 1)) {
  for (j in (i + 1):ncol(cormatrix)) {
    cor_test_result <- cor.test(cormatrix[, i], cormatrix[, j], method = "pearson")
    p_value_matrix[i, j] <- p_value_matrix[j, i] <- cor_test_result$p.value
  }
}
write.csv(corr,paste0(output_path,"/FRM_CRIX_Volatility_Corr.csv"),quote = FALSE)
write.csv(p_value_matrix,paste0(output_path,"/FRM_CRIX_Volatility_CorrP.csv"),quote = FALSE)


