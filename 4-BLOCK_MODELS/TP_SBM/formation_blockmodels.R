## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)
#setwd("D:/WORK_ALL/RECHERCHE/FORMATION_RECHERCHE/FORMATION_RESEAUX/2019_06_MIRES/Formation_Reseau_MIRES/4-BLOCK_MODELS/TP_SBM")


## ----import dataset,  echo = TRUE----------------------------------------
# Relation between individuals 
adjacency_individuals <- as.matrix(read.csv("ssna_adj_30n_all_ind_ssloop.csv", row.names=1))
adjacency_individuals_binary <- 1 * (adjacency_individuals > 0)
dim(adjacency_individuals)
isSymmetric(adjacency_individuals_binary)
diag(adjacency_individuals_binary)

# Inventory of species 

inventory_species <- read.csv("inventory_species.csv", row.names=1)
incidence_species <- as.matrix(inventory_species)
incidence_species_binary <- 1 * (incidence_species > 0)
dim(incidence_species)



## ----load pack,  echo = TRUE,message=FALSE-------------------------------
library(blockmodels)
library(igraph)
library(alluvial)
#library(ggplot2)


## ----source function,  echo = TRUE, message=FALSE------------------------
source('function_for_blockmodels.R')


## ----adjbin plot, echo=TRUE----------------------------------------------
plotMatrix(Mat = adjacency_individuals_binary,rowFG = 'individuals', colFG  = 'individuals')



## ----SBM, echo=TRUE, eval = TRUE-----------------------------------------
sbm.bern <- BM_bernoulli("SBM",adjacency_individuals_binary)
## ----estimate SBM, echo=TRUE, eval = FALSE-------------------------------
sbm.bern$estimate()

## ----save SBM, echo=FALSE, eval = FALSE----------------------------------
## save(sbm.bern,file = 'res_sbm_bern.Rdata')
## ----load SBM, echo=FALSE, eval = TRUE-----------------------------------
##load(file = 'res_sbm_bern.Rdata')

## ----select SBM, echo=TRUE, eval = TRUE----------------------------------
bestmodel = which.max(sbm.bern$ICL)
bestmodel


## ----extract param SBM, echo=TRUE, eval = TRUE---------------------------
paramEstimSBM.bern <- extractParamBM(sbm.bern,bestmodel)
paramEstimSBM.bern$pi
paramEstimSBM.bern$alpha
paramEstimSBM.bern$Z


## ----plot org adj bin,  echo=TRUE, eval = TRUE---------------------------
plotMatrix(adjacency_individuals_binary,'individuals','individuals', fileNameSave = NULL, clustering = list(row = paramEstimSBM.bern$Z))


## ----plot BM network adj,  echo=TRUE, eval = FALSE-----------------------
G <- graph_from_adjacency_matrix(paramEstimSBM.bern$alpha, mode = c("directed"), weighted = TRUE, diag = TRUE)
plot.igraph(G,vertex.size = paramEstimSBM.bern$pi*100,edge.curved = 0.2,edge.width = abs(E(G)$weight)*10,vertex.color = 1:paramEstimSBM.bern$Q, layout = layout_nicely)


## ----list names blocks,  echo=TRUE, eval = TRUE--------------------------
lapply(1:paramEstimSBM.bern$Q,function(q){rownames(adjacency_individuals_binary)[paramEstimSBM.bern$Z == q]})



########################### Passage au Poisson = graphe valué

## ----adj plot, echo=TRUE-------------------------------------------------
plotMatrix(Mat = adjacency_individuals,rowFG = 'individual', colFG  = 'individual')




## ----SBM Poisson, echo=TRUE, eval = TRUE---------------------------------
sbm.poiss <- BM_poisson("SBM",adjacency_individuals)


## ----estimate SBM poisson, echo=TRUE, eval = FALSE-----------------------
sbm.poiss$estimate()
## save(sbm.poiss,file = 'res_sbm_poiss.Rdata')


## ----sload SBM poisson,echo=FALSE, eval =TRUE----------------------------
##load('res_sbm_poiss.Rdata')


## ----select SBM poisson, echo=TRUE, eval = TRUE--------------------------
Q = which.max(sbm.poiss$ICL)
Q


## ----extract param SBM poisson, echo=TRUE, eval = TRUE-------------------
paramEstimSBM.pois <- extractParamBM(sbm.poiss,Q)
paramEstimSBM.pois$pi
paramEstimSBM.pois$lambda
paramEstimSBM.pois$Z


## ----plot org  ,  echo=TRUE, eval = TRUE---------------------------------
plotMatrix(adjacency_individuals,'individuals','individuals', fileNameSave = NULL, clustering = list(row = paramEstimSBM.pois$Z))


## ----plot BM network  poisson, echo=TRUE, eval = FALSE-------------------
G <- graph_from_adjacency_matrix(paramEstimSBM.pois$lambda, mode = c("directed"), weighted = TRUE, diag = TRUE)
plot.igraph(G,vertex.size=paramEstimSBM.pois$pi*100,edge.curved = 0.1,edge.width=abs(E(G)$weight)*0.5,vertex.color=1:paramEstimSBM.pois$Q, layout=layout_nicely)


## ----list names blocks Poisson,  echo=TRUE, eval = TRUE------------------
lapply(1:paramEstimSBM.pois$Q,function(q){rownames(adjacency_individuals_binary)[paramEstimSBM.pois$Z == q]})


## ----alluvial, echo=TRUE,eval=TRUE---------------------------------------
A <- as.data.frame(table(paramEstimSBM.bern$Z,paramEstimSBM.pois$Z))
colnames(A)=c('SBM Bern',"SBM Poisson","Freq")
w   <- which(A$Freq!=0)
A <- A[w,]
alluvial(A[,c(1,2)],freq=A$Freq)



######################## 
## ----incidence plot, echo=TRUE-------------------------------------------
plotMatrix(Mat = incidence_species_binary,rowFG = 'individuals', colFG  = 'species')


## ----LBM,echo=TRUE,eval=FALSE--------------------------------------------
lbm_bern <- BM_bernoulli("LBM",incidence_species_binary)
lbm_bern$estimate()

## ----LBM save,echo=FALSE,eval=FALSE--------------------------------------
##save(lbm_bern,file = 'res_lbm_bern.Rdata')


## ----LBM load,echo=FALSE,eval=TRUE---------------------------------------
##load(file = 'res_lbm_bern.Rdata')


## ----select lBM  , echo=TRUE, eval = TRUE--------------------------------
BestModel = which.max(lbm_bern$ICL) 
BestModel 
paramEstimLBM_bern <- extractParamBM(lbm_bern,BestModel )
paramEstimLBM_bern$Q



## ----extract param LBM, echo=TRUE, eval = TRUE---------------------------
paramEstimLBM_bern$piRow
paramEstimLBM_bern$piCol
paramEstimLBM_bern$alpha
paramEstimLBM_bern$ZRow
paramEstimLBM_bern$ZCol


## ----plot org  incidence,  echo=TRUE, eval = TRUE------------------------
plotMatrix(incidence_species_binary,'individuals','species', fileNameSave = NULL, clustering = list(row = paramEstimLBM_bern$ZRow,col = paramEstimLBM_bern$ZCol))

## ----plot BM network incidence,  echo=TRUE, eval = TRUE------------------
G <- graph_from_incidence_matrix(paramEstimLBM_bern$alpha, weighted = TRUE)
plot(G,vertex.size=c(paramEstimLBM_bern$piRow*100, paramEstimLBM_bern$piCol*100), vertex.shape=c("circle", "square")[V(G)$type +1], edge.width=abs(E(G)$weight*2),vertex.color=1:sum(paramEstimLBM_bern$Q), layout=layout.bipartite)



###############################" 
## ----incidence plot , echo=TRUE------------------------------------------
plotMatrix(Mat = incidence_species,rowFG = 'individuals', colFG  = 'species')


## ----LBM pois,echo=TRUE,eval=FALSE---------------------------------------
lbm_pois <- BM_poisson("LBM",incidence_species)
lbm_pois$estimate()

## ----LBM pois save,echo=FALSE,eval=FALSE---------------------------------
## save(lbm_pois,file = 'res_lbm_pois.Rdata')


## ----LBM pois load,echo=FALSE,eval=TRUE----------------------------------
##load(file = 'res_lbm_pois.Rdata')


## ----select lBM pois, echo=TRUE, eval = TRUE-----------------------------
BestModel = which.max(lbm_pois$ICL) 
BestModel 
paramEstimLBM_pois <- extractParamBM(lbm_pois,BestModel )
paramEstimLBM_pois$Q



## ----extract param LBM pois, echo=TRUE, eval = TRUE----------------------
paramEstimLBM_pois$piRow
paramEstimLBM_pois$piCol
paramEstimLBM_pois$alpha
paramEstimLBM_pois$ZRow
paramEstimLBM_pois$ZCol


## ----plot org  lbm pois,  echo=TRUE, eval = TRUE-------------------------
plotMatrix(incidence_species,'individuals','species', fileNameSave = NULL, clustering = list(row = paramEstimLBM_pois$ZRow,col = paramEstimLBM_pois$ZCol))

## ----plot BM network pois,  echo=TRUE, eval = TRUE-----------------------
G <- graph_from_incidence_matrix(paramEstimLBM_pois$alpha, weighted = TRUE)
plot(G,vertex.size = c(paramEstimLBM_pois$piRow*100, paramEstimLBM_pois$piCol * 100), vertex.shape = c("circle", "square")[V(G)$type + 1], edge.width=abs(E(G)$weight * 2),vertex.color = 1:sum(paramEstimLBM_pois$Q), layout = layout.bipartite)

