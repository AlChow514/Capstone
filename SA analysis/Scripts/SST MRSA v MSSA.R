# ---
# title: "SST MRSA v MSSA"
# author: "Casey Cazer"
# Last Updated: "May 4, 2020"
# ---


#identify the edges in SST contributed by MRSA vs MSSA isolates
#if I was to first divide SST isolates into MRSA and MSSA and then run entire pipeline, the graphs may look very differnt
##because changes in number of isolates, prevalences, etc will led to different eCSR p-values, cLift, etc

#I want to take the edges or sets found in SST and then calculate new cLift values from SST-MRSA vs SST-MSSA

#first divide SST db and transactions
SST.MRSA_db <- subset(SA.SST_db, Oxacillin==TRUE)
SST.MSSA_db <- subset(SA.SST_db, Oxacillin==FALSE)
SST.MRSA <- subset(SA.SST, items %in% "Oxacillin")
SST.MSSA <- subset(SA.SST, !items %in% "Oxacillin")

#get SST itemsets excluding those with cLift bootstrap crossing 1
SST.setList <- LIST(items(SA.SST_bestsets[!(SA.SST_bestsets_boot$cLiftBoot0.025<=1 & SA.SST_bestsets_boot$cLiftBoot0.975>=1)]), decode=FALSE)
SST.sets <- SA.SST_bestsets[!(SA.SST_bestsets_boot$cLiftBoot0.025<=1 & SA.SST_bestsets_boot$cLiftBoot0.975>=1)]


#calculate cLift from MRSA and MSSA subsets
#can't use cLift function because it reuses transaction set in interestMeasure
SST.MRSA.cLift <- as.data.frame(matrix(nrow=length(SST.setList), ncol=2))
SST.MSSA.cLift <- as.data.frame(matrix(nrow=length(SST.setList), ncol=2))

    #copied procedure from Mining Functions 1-29-2020; changed reuse=TRUE to reuse=FALSE
  AMsInSet <- lapply(SST.setList, function(i) itemLabels(SST.MRSA)[i])
  OnlyAMsInSet <- lapply(AMsInSet, function(i) subset(SST.MRSA_db, select=i) %>% na.omit()) #select only AM columns in the set and isolates that had been tested against all of them (no NA)
  setCount <- interestMeasure(SST.sets, "support", SST.MRSA, reuse=FALSE) * length(SST.MRSA) #number of isolates with the set (all AM in the set must have been tested in these isolates)
  numerator <- setCount/ldply(OnlyAMsInSet, nrow) #P(set | all AM in set are tested)
  denominator <- ldply(OnlyAMsInSet, function(i) prod(sapply(i, mean))) #find the support of each individual AM in the set (mean applied to logical vector) and multiply them together
  SST.MRSA.cLift[,2] <- numerator/denominator
  SST.MRSA.cLift[,1] <- ldply(AMsInSet, function(i) paste(unlist(i), collapse=","))
  colnames(SST.MRSA.cLift) <- c("items", "cLift")

  
  AMsInSet <- lapply(SST.setList, function(i) itemLabels(SST.MSSA)[i])
  OnlyAMsInSet <- lapply(AMsInSet, function(i) subset(SST.MSSA_db, select=i) %>% na.omit()) #select only AM columns in the set and isolates that had been tested against all of them (no NA)
  setCount <- interestMeasure(SST.sets, "support", SST.MSSA, reuse=FALSE) * length(SST.MSSA) #number of isolates with the set (all AM in the set must have been tested in these isolates)
  numerator <- setCount/ldply(OnlyAMsInSet, nrow) #P(set | all AM in set are tested)
  denominator <- ldply(OnlyAMsInSet, function(i) prod(sapply(i, mean))) #find the support of each individual AM in the set (mean applied to logical vector) and multiply them together
  SST.MSSA.cLift[,2] <- numerator/denominator
  SST.MSSA.cLift[,1] <- ldply(AMsInSet, function(i) paste(unlist(i), collapse=","))
  colnames(SST.MSSA.cLift) <- c("items", "cLift")
  

  rm(AMsInSet, OnlyAMsInSet, setCount, numerator, denominator)

#store the cLift values in the keep_edges dataframe for plotting
  ##note that other QM are not replaced and are therefore not valid for the MRSA/MSSA subsets
  SST.edges <- filter(keep_edges, Category=="SA.SST") #select the egde list for SA.SST
  SST.index <- match(SST.MRSA.cLift$items, SST.edges$items) #match the itemsets from SA.edges to SST.MRSA.cLift (same as SST.MSSA.cLift)
  SST.MRSA.edges <- SST.edges #duplicate to different variable
  SST.MRSA.edges$cLift <- replace(SST.MRSA.edges$cLift, SST.index, SST.MRSA.cLift$cLift) #replace cLift in SST.MRSA.edges with cLift from SST.MRSA.cLift, matching from index (visually confirmed indexing is correct)
  SST.MSSA.edges <- SST.edges
  SST.MSSA.edges$cLift <- replace(SST.MSSA.edges$cLift, SST.index, SST.MSSA.cLift$cLift)
  
  #drop rows where itemset was not found (all were found in MRSA, but not all in MSSA)
  ##although this is redundant with processing in plot_sets_weighted
  SST.MRSA.edges <- filter(SST.MRSA.edges, !is.na(cLift))
  SST.MSSA.edges <- filter(SST.MSSA.edges, !is.na(cLift))
  
#plotting, copied from AM circle plots 1-30-20
  QM.names <- c("cLift")
  QM.index <- match(QM.names, colnames(SST.MRSA.edges))
SST.MRSA_edges <- plot_sets_weighted(SST.MRSA.edges, "items", "Category", QM.index, QM.names, AM_class)
SST.MSSA_edges <- plot_sets_weighted(SST.MSSA.edges, "items", "Category", QM.index, QM.names, AM_class)

SST.MRSA_edges$cLiftBin <- cut(SST.MRSA_edges$cLift, breaks=c(0,0.25,0.5,0.75,1,2,5,10,Inf))
SST.MSSA_edges$cLiftBin <- cut(SST.MSSA_edges$cLift, breaks=c(0,0.25,0.5,0.75,1,2,5,10,Inf))


#copied from AM circle plots 1-30-20
#MRSA subset
  graph <- graph_from_data_frame(d=as.matrix(SST.MRSA_edges[,3:4]), vertices=vertices$AM, directed=FALSE)
  
  E(graph)$width <- SST.MRSA_edges[,"n"]+1 #start weight at 2
  E(graph)$width <- replace(E(graph)$width, E(graph)$width>10, 10) #cap weight at 10
  E(graph)$QM <- as.character(SST.MRSA_edges[,"cLiftBin"])
  E(graph)$color <- edge_colors[match(E(graph)$QM, names(edge_colors))]
  
  #vertex is pie chart of R, S, NA
  dat <- SST.MRSA_db
  dat <- dat[,colnames(dat) %in% vertices$AM] #drop cols not in vertices (because all susceptible)
  SST.MRSA.vertices <- vertices
  SST.MRSA.vertices$pieNA <- sapply(dat, function(x) sum(is.na(x))/length(x))
  SST.MRSA.vertices$pieR <- sapply(dat, function(x) sum(x, na.rm=TRUE)/length(x))
  SST.MRSA.vertices$pieS <- 1-SST.MRSA.vertices$pieNA - SST.MRSA.vertices$pieR
  
  #make vector = length of number of vertices plus S, NA. 
  SST.MRSA.vertices$pieVal <- mapply(append, lapply(1:nrow(SST.MRSA.vertices), function(x) matrix(0, nrow=1, ncol=nrow(SST.MRSA.vertices))), 
                                 mapply(append, as.list(SST.MRSA.vertices$pieS), as.list(SST.MRSA.vertices$pieNA), SIMPLIFY=FALSE), SIMPLIFY=FALSE)
  
  #then fill in appropriate R value for corresponding color values
  for (j in 1:nrow(SST.MRSA.vertices)){
    SST.MRSA.vertices$pieVal[[j]][j] <- SST.MRSA.vertices$pieR[j]
  }
  names(SST.MRSA.vertices$pieVal) <- SST.MRSA.vertices$AM
  
  #assign to graph
  V(graph)$label <- as.character(SST.MRSA.vertices$label)
  V(graph)$label.color <- "black"
  V(graph)$label.cex <- 1.25
  V(graph)$frame.color <- SST.MRSA.vertices$color
  V(graph)$shape <- "pie"
  V(graph)$pie <- SST.MRSA.vertices$pieVal
  V(graph)$pie.color <- longcolors
  
  #title
  graph_attr(graph, "main") <- "SSSI MRSA"
  
  #layout
  graph$layout <- good.coords
  
  #save plots
  label <- "g_AM_SST.MRSA"
  assign(label, graph)
  
  #save as png to maintain proportions across graphs
  filename <- paste('./figures/FigS4_B_SST.MRSA_AMsets_cLift.png', sep="")
  Cairo(file=filename, type="png", units="in", width=7, height=8, pointsize=8, res=300)
  plot.new()
  par(cex.main=4)
  plot(graph)
  mtext(LETTERS[2], side=3, adj=0, cex=4)
  dev.off()

#MSSA subset
  graph <- graph_from_data_frame(d=as.matrix(SST.MSSA_edges[,3:4]), vertices=vertices$AM, directed=FALSE)
  
  E(graph)$width <- SST.MSSA_edges[,"n"]+1 #start weight at 2
  E(graph)$width <- replace(E(graph)$width, E(graph)$width>10, 10) #cap weight at 10
  E(graph)$QM <- as.character(SST.MSSA_edges[,"cLiftBin"])
  E(graph)$color <- edge_colors[match(E(graph)$QM, names(edge_colors))]
  
  #vertex is pie chart of R, S, NA
  dat <- SST.MSSA_db
  dat <- dat[,colnames(dat) %in% vertices$AM] #drop cols not in vertices (because all susceptible)
  SST.MSSA.vertices <- vertices
  SST.MSSA.vertices$pieNA <- sapply(dat, function(x) sum(is.na(x))/length(x))
  SST.MSSA.vertices$pieR <- sapply(dat, function(x) sum(x, na.rm=TRUE)/length(x))
  SST.MSSA.vertices$pieS <- 1-SST.MSSA.vertices$pieNA - SST.MSSA.vertices$pieR
  
  #make vector = length of number of vertices plus S, NA. 
  SST.MSSA.vertices$pieVal <- mapply(append, lapply(1:nrow(SST.MSSA.vertices), function(x) matrix(0, nrow=1, ncol=nrow(SST.MSSA.vertices))), 
                                     mapply(append, as.list(SST.MSSA.vertices$pieS), as.list(SST.MSSA.vertices$pieNA), SIMPLIFY=FALSE), SIMPLIFY=FALSE)
  
  #then fill in appropriate R value for corresponding color values
  for (j in 1:nrow(SST.MSSA.vertices)){
    SST.MSSA.vertices$pieVal[[j]][j] <- SST.MSSA.vertices$pieR[j]
  }
  names(SST.MSSA.vertices$pieVal) <- SST.MSSA.vertices$AM
  
  #assign to graph
  V(graph)$label <- as.character(SST.MSSA.vertices$label)
  V(graph)$label.color <- "black"
  V(graph)$label.cex <- 1.25
  V(graph)$frame.color <- SST.MSSA.vertices$color
  V(graph)$shape <- "pie"
  V(graph)$pie <- SST.MSSA.vertices$pieVal
  V(graph)$pie.color <- longcolors
  
  #title
  graph_attr(graph, "main") <- "SSSI MSSA"
  
  #layout
  graph$layout <- good.coords
  
  #save plots
  label <- "g_AM_SST.MSSA"
  assign(label, graph)
  
  #save as png to maintain proportions across graphs
  filename <- paste('./figures/FigS4_A_SST.MSSA_AMsets_cLift.png', sep="")
  Cairo(file=filename, type="png", units="in", width=7, height=8, pointsize=8, res=300)
  plot.new()
  par(cex.main=4)
  plot(graph)
  mtext(LETTERS[1], side=3, adj=0, cex=4)
  dev.off()
  
#plot together
filename <- "figures/FigS4_SST MRSA v MSSA AM sets.png"

Cairo(file=filename, type="png", units="in", width=6, height=3, pointsize=8, res=300)
plot.new()
grid.arrange(rasterGrob(readPNG("figures/FigS4_A_SST.MSSA_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS4_B_SST.MRSA_AMsets_cLift.png")),
             ncol=2, nrow=1)
dev.off()