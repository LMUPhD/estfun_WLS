library(openxlsx)

setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS/res")



for(n in c("500","1000")){
  for(s in c("1","2","4")){
    
    load(paste0("0724_yves_complex/240715_",s,"_",n,"_multi_none_complex.RData"))
    
    wb=loadWorkbook("0724_yves_complex/multivariate_GEE.xlsx")
    sh=addWorksheet(wb,paste0("n=",n,"_sc=",s))
    writeData(wb, sh, resu, startCol = 1, startRow = 1, rowNames = T)
    writeData(wb, sh, resu_t, startCol = 1, startRow = 10, rowNames = T)
    
    
    saveWorkbook(wb, "0724_yves_complex/multivariate_GEE.xlsx", overwrite = TRUE)
    
  }
}


