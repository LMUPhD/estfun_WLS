library(openxlsx)

setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS/res")



for(n in c("500","1000","2000")){
  for(s in c("1","2","4","6")){
    
    load(paste0("univariate/only betas/240210_",n,"_",s,"_univ_betas.RData"))
    
    wb=loadWorkbook("univariate/only betas/univariate_beta.xlsx")
    sh=addWorksheet(wb,paste0("n=",n,"_sc=",s))
    writeData(wb, sh, resu, startCol = 1, startRow = 1, rowNames = T)
    writeData(wb, sh, resu_t, startCol = 1, startRow = 10, rowNames = T)
    
    
    saveWorkbook(wb, "univariate/only betas/univariate_beta.xlsx", overwrite = TRUE)
    
  }
}


