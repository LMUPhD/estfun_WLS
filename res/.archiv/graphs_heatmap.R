

library(ggplot2)
library(ComplexHeatmap)
library(circlize)
setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS/res")

################################################################################
######################## Multivariate Plot #####################################
################################################################################

cond = c("mirt","WLS")
stat = c("maxBB","ordmax","ordL2BB","catL2BB")
for(c in cond){
  for(t in stat){
    assign(paste0("heatmap_",c,"_",t),matrix(NA,nrow=4,ncol=3))
  }
}


anzahl = c("500","1000","2000")
schwellen = c("1","2","4","6")
for(n in 1:3){
  for(s in 1:4){

    load(paste0("simpleWLS/multivariate/231125_",schwellen[s],"_",anzahl[n],"_result_multi_simpleWLS.RData"))
    
    
    for(c in 1:2){
      for(t in 1:4){
        mat = get(paste0("heatmap_",cond[c],"_",stat[t]))
        mat[s,n] <- resu[t,c]
        if(n==1 & s==1){
          colnames(mat) <- c("n=500","n=1000","n=2000")
          rownames(mat) <- c("k=1","k=2","k=4","k=6")
        }
        assign(paste0("heatmap_",cond[c],"_",stat[t]),mat)
      }
    }
  }
}



######## Draw Plot
mat11 = heatmap_mirt_maxBB
mat12 = heatmap_mirt_ordmax
mat13 = heatmap_mirt_ordL2BB
mat14 = heatmap_mirt_catL2BB
mat21 = heatmap_WLS_maxBB
mat22 = heatmap_WLS_ordmax
mat23 = heatmap_WLS_ordL2BB
mat24 = heatmap_WLS_catL2BB

#color scheme
col_fun = colorRamp2(c(0, 0.05, 0.1), c("yellow", "white", "darkred"))


#save as pdf
pdf(file="multivariate_plot.pdf",width=10,height=10)


#draw grid
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 3,widths=c(1,1,0.2))))

#annotations
ha = HeatmapAnnotation(foo = anno_block(labels="mirt",gp = gpar(fill = "lightgrey")))
hl_maxbb = rowAnnotation(foo = anno_block(labels="maxBB",gp = gpar(fill = "lightgrey")))
hl_ordmax = rowAnnotation(foo = anno_block(labels="ordmax",gp = gpar(fill = "lightgrey")))
hl_ordL2BB = rowAnnotation(foo = anno_block(labels="ordL2BB",gp = gpar(fill = "lightgrey")))
hl_catl2BB = rowAnnotation(foo = anno_block(labels="catL2BB",gp = gpar(fill = "lightgrey")))

#first columns heatmaps
h11 = Heatmap(mat11, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2),top_annotation=ha, left_annotation = hl_maxbb,cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat11[i, j],2), x, y, gp = gpar(fontsize = 10))})
h12 = Heatmap(mat12, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2), left_annotation = hl_ordmax,                 cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat12[i, j],2), x, y, gp = gpar(fontsize = 10))})
h13 = Heatmap(mat13, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2), left_annotation = hl_ordL2BB,                cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat13[i, j],2), x, y, gp = gpar(fontsize = 10))})
h14 = Heatmap(mat14, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2), left_annotation = hl_catl2BB,                cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat14[i, j],2), x, y, gp = gpar(fontsize = 10))})
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(h11%v%h12%v%h13%v%h14, newpage = FALSE)
upViewport()

#second column heatmaps
ha = HeatmapAnnotation(foo = anno_block(labels="WLS",gp = gpar(fill = "lightgrey")))
h21 = Heatmap(mat21, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2),top_annotation=ha,cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat21[i, j],2), x, y, gp = gpar(fontsize = 10))})
h22 = Heatmap(mat22, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2)                  ,cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat22[i, j],2), x, y, gp = gpar(fontsize = 10))})
h23 = Heatmap(mat23, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2)                  ,cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat23[i, j],2), x, y, gp = gpar(fontsize = 10))})
h24 = Heatmap(mat24, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2)                  ,cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat24[i, j],2), x, y, gp = gpar(fontsize = 10))})
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(h21%v%h22%v%h23%v%h24, newpage = FALSE)
upViewport()

#legend
lgd = Legend(at = c(0, 0.05, 0.1), col_fun = col_fun, legend_height = unit(5, "cm"),title = "alpha\nerror")
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
grid.draw(lgd)
upViewport()

##
upViewport()
dev.off()




################################################################################
######################### Univariate Plot ######################################
################################################################################


cond = c("mirt","WLS")
stat = c("maxBB","meanL2BB","supLM","ordL2BB","ordmax","catL2BB")
for(c in cond){
  for(t in stat){
    assign(paste0("heatmap_",c,"_",t),matrix(NA,nrow=4,ncol=3))
  }
}


anzahl = c("500","1000","2000")
schwellen = c("1","2","4","6")
for(n in 1:3){
  for(s in 1:4){
    
    load(paste0("simpleWLS/univariate/231124_",anzahl[n],"_",schwellen[s],"_result_univ_simpleWLS.RData"))
    
    
    for(c in 1:2){
      for(t in 1:6){
        mat = get(paste0("heatmap_",cond[c],"_",stat[t]))
        mat[s,n] <- resu[t,c]
        if(n==1 & s==1){
          colnames(mat) <- c("n=500","n=1000","n=2000")
          rownames(mat) <- c("k=1","k=2","k=4","k=6")
        }
        assign(paste0("heatmap_",cond[c],"_",stat[t]),mat)
      }
    }
  }
}



######## Draw Plot
mat11 = heatmap_mirt_supLM
mat12 = heatmap_mirt_maxBB
mat13 = heatmap_mirt_meanL2BB;mat13=mat13[1:3,]
mat14 = heatmap_mirt_ordmax
mat15 = heatmap_mirt_ordL2BB
mat16 = heatmap_mirt_catL2BB

mat21 = heatmap_WLS_supLM
mat22 = heatmap_WLS_maxBB
mat23 = heatmap_WLS_meanL2BB;mat23=mat23[1:3,]
mat24 = heatmap_WLS_ordmax
mat25 = heatmap_WLS_ordL2BB
mat26 = heatmap_WLS_catL2BB

#color scheme
col_fun = colorRamp2(c(0, 0.05, 0.1), c("yellow", "white", "darkred"))


#save as pdf
pdf(file="univariate_plot.pdf",width=10,height=10)


#draw grid
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 3,widths=c(1,1,0.2))))

#annotations first col
ha = HeatmapAnnotation(foo = anno_block(labels="mirt",gp = gpar(fill = "lightgrey")))
hl_suplm = rowAnnotation(foo = anno_block(labels="suplm",gp = gpar(fill = "lightgrey")))
hl_maxbb = rowAnnotation(foo = anno_block(labels="maxBB",gp = gpar(fill = "lightgrey")))
hl_meanl2bb = rowAnnotation(foo = anno_block(labels="meanL2BB",gp = gpar(fill = "lightgrey")))
hl_ordmax = rowAnnotation(foo = anno_block(labels="ordmax",gp = gpar(fill = "lightgrey")))
hl_ordL2BB = rowAnnotation(foo = anno_block(labels="ordL2BB",gp = gpar(fill = "lightgrey")))
hl_catl2BB = rowAnnotation(foo = anno_block(labels="catL2BB",gp = gpar(fill = "lightgrey")))

#first columns heatmaps
h11 = Heatmap(mat11, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2),top_annotation=ha, left_annotation = hl_suplm,cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat11[i, j],2), x, y, gp = gpar(fontsize = 10))})
h12 = Heatmap(mat12, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2), left_annotation = hl_maxbb,                  cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat12[i, j],2), x, y, gp = gpar(fontsize = 10))})
h13 = Heatmap(mat13, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2), left_annotation = hl_meanl2bb,               cell_fun = function(j, i, x, y, width, height, fill) { if(mat13[i, j] < 0.5){grid.text( round(mat13[i, j],2), x, y, gp = gpar(fontsize = 10))}})
h14 = Heatmap(mat14, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2), left_annotation = hl_ordmax,                 cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat14[i, j],2), x, y, gp = gpar(fontsize = 10))})
h15 = Heatmap(mat15, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2), left_annotation = hl_ordL2BB,                cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat15[i, j],2), x, y, gp = gpar(fontsize = 10))})
h16 = Heatmap(mat16, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2), left_annotation = hl_catl2BB,                cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat16[i, j],2), x, y, gp = gpar(fontsize = 10))})

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(h11%v%h12%v%h13%v%h14%v%h15%v%h16, newpage = FALSE)
upViewport()

#annotations second col
ha = HeatmapAnnotation(foo = anno_block(labels="WLS",gp = gpar(fill = "lightgrey")))

#second column heatmaps
h21 = Heatmap(mat21, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2),top_annotation=ha, cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat21[i, j],2), x, y, gp = gpar(fontsize = 10))})
h22 = Heatmap(mat22, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2),                   cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat22[i, j],2), x, y, gp = gpar(fontsize = 10))})
h23 = Heatmap(mat23, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2),                   cell_fun = function(j, i, x, y, width, height, fill) { if(mat23[i, j] < 0.5){grid.text(round(mat23[i, j],2), x, y, gp = gpar(fontsize = 10))}})
h24 = Heatmap(mat24, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2),                   cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat24[i, j],2), x, y, gp = gpar(fontsize = 10))})
h25 = Heatmap(mat25, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2),                   cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat25[i, j],2), x, y, gp = gpar(fontsize = 10))})
h26 = Heatmap(mat26, border=T, col = col_fun, show_heatmap_legend = FALSE,cluster_rows = FALSE,show_column_dend = FALSE,rect_gp = gpar(col = "white", lwd = 2),                   cell_fun = function(j, i, x, y, width, height, fill) { grid.text(round(mat26[i, j],2), x, y, gp = gpar(fontsize = 10))})
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(h21%v%h22%v%h23%v%h24%v%h25%v%h26, newpage = FALSE)
upViewport()

#legend
lgd = Legend(at = c(0, 0.05, 0.1), col_fun = col_fun, legend_height = unit(5, "cm"),title = "alpha\nerror")
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
grid.draw(lgd)
upViewport()

##
upViewport()
dev.off()


