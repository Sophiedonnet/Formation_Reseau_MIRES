rm(list=ls())
library(GREMLIN)
library(blockmodels)
library(igraph)
library(alluvial)
library(ggplot2)
source('function_for_blockmodels.R')

setwd("C:/Users/Donnet/Dropbox/WORK_DROPBOX/RECHERCHE_en_cours/CONFERENCE_EXPOSES/EXPOSES/2019/2019-05-GDR-Resodiv/codes")
setwd("~/Dropbox/WORK_DROPBOX/RECHERCHE_en_cours/CONFERENCE_EXPOSES/EXPOSES/2019/2019-05-GDR-Resodiv/codes")
file_data <- "ssna_adj_30n_all_ind_ssloop.csv"


############ DATA VANUATU RELATIONS
data_relation <-  read.csv(file  = file_data)
data_relation <- data_relation[, -1]
data_relation <- as.matrix(data_relation[1:30,1:30])
data_relation_bin <- (data_relation>0)* matrix(1,30,30)

rownames(data_relation_bin) <- 1:30
colnames(data_relation_bin) <- 1:30


g <- graph_from_adjacency_matrix(data_relation_bin,mode = "directed", diag=FALSE )
plot.igraph(g, edge.label.cex=0.7,arrow.size  = 0.01, arrow.width = 0.5)

plotMatrix(Mat = data_relation_bin,rowFG = 'foyer', colFG  = 'foyer')



#--------- SBM 

sbm.relations <- BM_bernoulli("SBM",data_relation_bin)
sbm.relations$estimate()
Q = which.max(sbm.relations$ICL)
paramEstimSBM <- extractParamBM(sbm.relations,Q)

plotMatrix(data_relation_bin,'Foyer','Foyer', fileNameSave = NULL, clustering = list(row = paramEstimSBM$Z))


G <- graph_from_adjacency_matrix(paramEstimSBM$alpha, mode = c("directed"), weighted = TRUE, diag = TRUE)
G %>% set_edge_attr("curved", value=100)
plot.igraph(G,vertex.size = paramEstimSBM$pi*100,vertex.color = 1:Q)

plotNetBM(sbm.relations,Q)

##########################################" 
#--------------------- LBM 
inventory_species <- read.csv("~/Dropbox/WORK_DROPBOX/RECHERCHE_en_cours/CONFERENCE_EXPOSES/EXPOSES/2019/2019-05-GDR-Resodiv/codes/inventory_species.csv")
inventory <- inventory_species[, -1]
inventory <- as.matrix(inventory[1:30,1:37])
inventory_bin <- (inventory > 0)* matrix(1,30,37)

rownames(inventory_bin) <- 1:30
colnames(inventory_bin) <- 1:37


g <- graph_from_incidence_matrix(inventory_bin)
V(g)$type
plot(g, layout = layout_as_bipartite, vertex.color=c("green","cyan")[V(g)$type+1])
plotMatrix(inventory_bin,'Foyer','Espèce', fileNameSave = NULL)


#--------- SBM 

lbm <- BM_bernoulli("LBM",inventory_bin)
lbm$estimate()
Q = which.max(lbm$ICL)
paramEstimLBM <- extractParamBM(lbm,Q)

plotMatrix(inventory_bin,'Foyer','Espece', fileNameSave = NULL, clustering = list(row = paramEstimLBM$ZRow, col = paramEstimLBM$ZCol))

G <- graph_from_incidence_matrix(paramEstimLBM$alpha, weighted = TRUE)
plot.igraph(G,vertex.size = c(paramEstimLBM$piRow, paramEstimLBM$piCol)*100,vertex.color = c(rep('green',3),rep("cyan",2)), layout  = layout_as_bipartite)







