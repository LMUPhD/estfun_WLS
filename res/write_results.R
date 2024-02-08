library(openxlsx)

setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS/res")


for(n in c("3000")){
  for(s in c("1","2","4","6")){
    
    load(paste0("complexWLS/univariate/231206_",n,"_",s,"_result_univ_complexWLS.RData"))
    
    wb=loadWorkbook("complexWLS/univariate_complex.xlsx")
    sh=addWorksheet(wb,paste0("n=",n,"_sc=",s))
    writeData(wb, sh, resu, startCol = 1, startRow = 1, rowNames = T)
    writeData(wb, sh, resu_t, startCol = 1, startRow = 10, rowNames = T)
    
    
    saveWorkbook(wb, "complexWLS/univariate_complex.xlsx", overwrite = TRUE)
    
  }
}


for(n in c("3000")){
  for(s in c("1","2","4","6")){
    
    load(paste0("complexWLS/multivariate/231201_",s,"_",n,"_result_multi_complexWLS.RData"))
    
    wb=loadWorkbook("complexWLS/multivariate_complex.xlsx")
    sh=addWorksheet(wb,paste0("n=",n,"_sc=",s))
    writeData(wb, sh, resu, startCol = 1, startRow = 1, rowNames = T)
    writeData(wb, sh, resu_t, startCol = 1, startRow = 10, rowNames = T)
    
    
    saveWorkbook(wb, "complexWLS/multivariate_complex.xlsx", overwrite = TRUE)
    
  }
}


