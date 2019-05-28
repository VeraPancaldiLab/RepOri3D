
###function to calculate ChAS of many features

calc_assort=function(G, data){
names=colnames(data)
ass=list()
G_epi=list()
for (i in c(1:ncol(data))){
#print(i)
#print(names[i])
G=set.vertex.attribute(G, names[i],value= data[V(G)$name, names[i]])
attsel=which( names(vertex.attributes(G))==names[i])
G_epi[[names[i]]]=delete.vertices(G, V(G)[is.na(vertex.attributes(G)[[attsel]]) ])
ass[[names[i]]]=assortativity(G_epi[[names[[i]]]], types1=vertex.attributes(G_epi[[names[i]]])[[attsel]], directed=F)
#print (ass[[names[i]]])
}

return(ass)
}
