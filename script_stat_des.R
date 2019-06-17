### Manipulation de réseaux setwd("C:/Massol/Enseignement/Formation réseaux/Formation Resodiv Juin 2019/Données Sophie")

library(igraph)

data.network<-read.table("ssna_adj_30n_all_ind_ssloop.csv",sep=",",header=T,row.names=1)
network<-graph_from_adjacency_matrix(as.matrix(data.network), mode = "directed",weighted=T)

plot(network,layout = layout_with_kk,vertex.label=substr(V(network)$name,5,7))

data.network.bin<-data.network
data.network.bin[data.network.bin>0]<-1
network.bin<-graph_from_adjacency_matrix(as.matrix(data.network.bin), mode = "directed")
network.bin.undirected<-graph_from_adjacency_matrix(as.matrix(data.network.bin), mode = "undirected")

gorder(network.bin)
gsize(network.bin)
edge_density(network.bin)
gsize(network.bin)/(gorder(network.bin)*(gorder(network.bin)-1))
reciprocity(network.bin)
transitivity(network,type="barrat")

### Degrés, centralités

degree(network.bin)
degree(network.bin,mode="in")
degree(network.bin,mode="out")
degree_distribution(network.bin)
degree_distribution(network.bin,mode="in")
degree_distribution(network.bin,mode="out")

plot(degree(network.bin,mode="in")~degree(network.bin,mode="out"))

pseudolog<-function(x) ifelse(x==0,0,log(x))
notzero<-function(x) ifelse(x==0,1,x)

dfun<-function(mat) {
	n<-dim(mat)[1]
	a<-sum(mat)
	ai<-notzero(apply(mat,1,sum))
	aj<-notzero(apply(mat,2,sum))
	ai_mat<-matrix(rep(ai,n),nrow=n,ncol=n,byrow=F)
	aj_mat<-matrix(rep(aj,n),nrow=n,ncol=n,byrow=T)
	elem<-(mat/ai_mat)*sapply(a*mat/(ai_mat*aj_mat),pseudolog)
	d<-apply(elem,1,sum)
	dmax<-pseudolog(a/ai)
	d/dmax
}

dprime.out<-dfun(data.network)
dprime.in<-dfun(t(data.network))

plot(dprime.in~dprime.out,xlim=c(0,1),ylim=c(0,1))

ev.cent<-eigen_centrality(network.bin.undirected)
clo.cent<-closeness(network.bin.undirected,normalized=T)
be.cent<-betweenness(network.bin.undirected,normalized=T)

ev.vect<-ev.cent$vector
plot(ev.vect~degree(network.bin))

kb.cent<-sapply((1:100)/20,function(x) alpha_centrality(network.bin, alpha = x, exo = 1))
alpha<-(1:100)/20

regressioncoef<-function(var) {
gls(var ~ ev.vect)$coef[2]
}
alpha.kb.ev.coef<-sapply(1:100,function(y) regressioncoef(kb.cent[,y]))

plot(alpha.kb.ev.coef~alpha)

kleinberg.auth<-authority_score(network.bin)
kleinberg.hub<-hub_score(network.bin)

plot(kleinberg.auth$vector~degree(network.bin,mode="in"))

### Modularité
EB.mod<-cluster_edge_betweenness(network.bin.undirected)
LE.mod<-cluster_leading_eigen(network.bin.undirected)
ML.mod<-cluster_louvain(network.bin.undirected)

plot(LE.mod,network.bin,layout = layout_with_mds,vertex.label=substr(V(network)$name,5,7))

IM.mod<-cluster_infomap(network.bin)
plot(IM.mod,network.bin,layout = layout_with_mds,vertex.label=substr(V(network)$name,5,7))


### Randomisation et tests
sample.config.directed<-lapply(1:100,function(x) sample_degseq(degree(network.bin,mode="out"), degree(network.bin,mode="in"), method = "simple.no.multiple"))

sample.config.undirected<-lapply(1:100,function(x) sample_degseq(degree(network.bin.undirected), method = "vl"))

plot(sample.config.undirected[[5]],layout = layout_with_kk,vertex.label=NULL)

LE.mod$mod
LE.mods<-sapply(1:100,function(x) cluster_leading_eigen(sample.config.undirected[[x]])$mod)
LE.mod.distrib<-ecdf(LE.mods)
1-LE.mod.distrib(LE.mod$mod)

plot(density(LE.mods))
abline(v=LE.mod$mod,col="red")

data.bip<-read.table("inventory_species.csv",sep=",",header=T,row.names=1)
bip.bin<-graph_from_incidence_matrix(as.matrix(data.bip),weighted=NULL)
bip.mod<-cluster_leading_eigen(bip.bin)$mod

library(vegan)
sample.bip.config<-simulate(nullmodel(data.bip,"curveball"),nsim=1000)
dim(sample.bip.config)
sample.mods<-sapply(1:1000,function(x) cluster_leading_eigen(graph_from_incidence_matrix(as.matrix(sample.bip.config[,,x])))$mod)
sample.mod.distrib<-ecdf(sample.mods)
1-sample.mod.distrib(bip.mod)
