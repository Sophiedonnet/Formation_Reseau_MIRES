## Fonctions d'échantillonnages de réseau



##construit l'ego-network d'un graphe associe a un sous-ensemble de noeud
#' sample.snowball
#'
#' @param graph un objet igraph
#' @param node un vecteur de noeud exploré
#' @param mode "in", "out" ou "all" type de voisinage choisi
#'
#' @return un objet igraph correspondant au graph échantillonné
#

build_egonetwork=  function(graph,node,mode="all"){
  neighbors.sampled = neighbors.build(graph,node,mode)
  graph.sampled <- induced_subgraph(graph,vids=unique(c(node, neighbors.sampled)))
  V(graph.sampled)$explored = rep(F,length(V(graph.sampled)))
  V(graph.sampled)$explored[V(graph.sampled)[node$name]] <-T
  V(graph.sampled)$color= ifelse(V(graph.sampled)$explored==T,"red","blue")
  x= V(graph.sampled)
  y= V(graph.sampled)[V(graph.sampled)$explored==F]
  if (mode=="in") {graph.sampled=delete_edges(graph.sampled,E(graph.sampled)[x %->% y])}
  if (mode=="out") {graph.sampled=delete_edges(graph.sampled,E(graph.sampled)[y %->% x])}
  if (mode=="all") {graph.sampled=delete_edges(graph.sampled,E(graph.sampled)[y %--% y])}
  return(graph.sampled)
}


##simulation d'un échantillonnage d'ego-network
#' sample.egonetwork
#'
#' @param graph un objet igraph
#' @param budget un entier correspondant au nombre de noeud échantillonné
#' @param mode "in", "out" ou "all" type de voisinage choisi
#'
#' @return un objet igraph correspondant au graph échantillonné
#

sample.egonetwork = function(graph, budget, mode="all"){
  nodes.sampled = sample(V(graph), budget)
  network.sampled = build_egonetwork(graph,nodes.sampled,mode)
  return(network.sampled)
}

##  Simulation d'un échantillonnage snowball
#' sample.snowball
#'
#' @param graph un objet igraph
#' @param depth un entier correspondant au nombre d'étape de l'échantillonage snowball
#' @param budget un vecteur d'entier de longueur supérieur ou égale à depth
#' @param mode "in", "out" ou "all" type de voisinage choisi
#'
#' @return un objet igraph correspondant au graph échantillonné
#
# Détails: à chaque étape, les noeuds sont échantillonné uniformément parmi-ceux qui sont 
# découvert mais non explorés. S'il n'y a plus assez de noeuds découverts, de nouveaux noeuds sont tirés 
# uniformément parmi ceux non découverts

sample.snowball= function(graph, depth, budget, mode="all" ){
  nodes.sampled = sample(V(graph), budget[1])
  for (i in (2:depth)){
    neighb= neighbors.build(graph,nodes.sampled,mode)
    nodes.sampled.snow=   sample(neighb , min(budget[i],length(neighb)))
    if (budget[i]>length(neighb)) {
      nodes.sampled.snow= c(nodes.sampled.snow, sample(difference(V(graph),nodes.sampled.snow) , budget[i]- length(neighb))) 
      }
    nodes.sampled= unique(c(nodes.sampled,nodes.sampled.snow))
  }
  network.sampled = build_egonetwork(graph,nodes.sampled,mode)
 return(network.sampled)
}

##construction du voisinage d'ordre 1 d'un vecteur de noeud
#' neighbors.build
#'
#' @param graph un objet igraph
#' @param node un vecteur de noeud de graph
#' @param mode "in", "out" ou "all" type de voisinage choisi
#'
#' @return le vecteur des voisins d'au moins un noeud de node qui n'appartiennent pas a node

neighbors.build = function(graph,node,mode="all"){
  list.neighbors = sapply(node, function(x){ as_ids(neighbors(graph,x,mode))})
  vect.neighbors= V(graph)[unique(unlist(list.neighbors,use.names=F))]
  return(difference(unique(c(vect.neighbors,node)),node))
}



## Calcul du sous-graphe induit par les noeud exploré
#' induced.sample
#'
#' @param graph un objet igraph dont les noeuds contienne l'attribut explored
#'
#' @return un objet igraph

induced.sample = function(graph){
  network.induced = induced_subgraph(graph,V(graph)[V(graph)$explored==T])  
  return(network.induced)
}


