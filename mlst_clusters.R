setwd("~/Documents/Postdoc/assigning/revisions/")
require(igraph)
require(mclust)

get_dists = function(dataframe) {
  dists = matrix(data=NA, nrow=nrow(dataframe), ncol=nrow(dataframe),
                 dimnames=list(rownames(dataframe), rownames(dataframe)))
  for (i in 1:nrow(dataframe)) {
   for (j in 1:nrow(dataframe)) {
      dists[i,j] = sum(dataframe[i,] == dataframe[j,])
   }
  }
  dists = ncol(dataframe) - dists
  diag(listeria_mlst_dists) = NA # remove self edges
  return(dists)
}

eval_mlst = function(mlst_dists, changes, outfile, pp_clusters, change=TRUE)
{
  strain.graph = graph_from_adjacency_matrix(mlst_dists <= changes, mode="undirected")
  strain.components<-components(strain.graph)
  if (change == TRUE)
  {
    taxon = gsub(pattern = '_', replacement = "#", x = names(strain.components$membership))
    taxon = sub(pattern = '#', replacement = "_", x = taxon)
    taxon = paste0(taxon, '.contigs_velvet.fa')
  }
  else
  {
    taxon = names(strain.components$membership)
  }
  mlst_diff = data.frame(Taxon=taxon,
                         cluster_mlst=strain.components$membership,
                         stringsAsFactors = F)
  write.table(mlst_diff,file=outfile,quote = F,col.names = T,row.names = F,sep=",")
  
  mergey = merge(pp_clusters, mlst_diff, by = "Taxon")
  print(strain.components$no)
  print(adjustedRandIndex(mergey$Cluster, mergey$cluster_mlst))
}

# Listeria MLST
listeria_poppunk = read.csv("~/Documents/Postdoc/assigning/revisions/listeria_db_revisions/listeria_db_revisions_clusters.csv", 
                            stringsAsFactors=FALSE)
#listeria_poppunk$Taxon = gsub(pattern = '.contigs.velvet.fa', replacement = "", x = listeria_poppunk$Taxon)

listeria_mlst <- read.delim("~/Documents/Postdoc/assigning/revisions/listeria_st/listeria_mlst.txt", 
                            stringsAsFactors=FALSE)
rownames(listeria_mlst) = listeria_mlst$Sample
listeria_mlst = listeria_mlst[,2:8]
listeria_mlst_dists = get_dists(listeria_mlst)

# same (ST)
eval_mlst(listeria_mlst_dists, 0, "listeria_st/mlst_0_clusters.csv", listeria_poppunk)
# 1 diff
eval_mlst(listeria_mlst_dists, 1, "listeria_st/mlst_1_clusters.csv", listeria_poppunk)
# 3 diffs
eval_mlst(listeria_mlst_dists, 3, "listeria_st/mlst_3_clusters.csv", listeria_poppunk)

# Listeria cgMLST
listeria_cgmlst <- read.delim("~/Documents/Postdoc/assigning/revisions/listeria_st/listeria_cgmlst.txt", stringsAsFactors=FALSE)
rownames(listeria_cgmlst) = gsub(pattern = '.contigs.velvet.fa', replacement = "", x = listeria_cgmlst$FILE)
rownames(listeria_cgmlst) = gsub(pattern = '#', replacement = "_", x = rownames(listeria_cgmlst))
listeria_cgmlst = listeria_cgmlst[,-1]
listeria_cgmlst_dists = get_dists(listeria_cgmlst)

# same (ST)
eval_mlst(listeria_cgmlst_dists, 0, "listeria_st/cgmlst_0_clusters.csv", listeria_poppunk)
# 1 diff
eval_mlst(listeria_cgmlst_dists, 1, "listeria_st/cgmlst_1_clusters.csv", listeria_poppunk)
# 3 diffs
eval_mlst(listeria_cgmlst_dists, 3, "listeria_st/cgmlst_3_clusters.csv", listeria_poppunk)
# 10 diffs
eval_mlst(listeria_cgmlst_dists, 10, "listeria_st/cgmlst_10_clusters.csv", listeria_poppunk)
# 30 diffs
eval_mlst(listeria_cgmlst_dists, 30, "listeria_st/cgmlst_30_clusters.csv", listeria_poppunk)

# Listeria MLST
listeria_poppunk = read.csv("~/Documents/Postdoc/assigning/revisions/listeria_db_revisions/listeria_db_revisions_clusters.csv", 
                            stringsAsFactors=FALSE)
#listeria_poppunk$Taxon = gsub(pattern = '.contigs.velvet.fa', replacement = "", x = listeria_poppunk$Taxon)

listeria_mlst <- read.delim("~/Documents/Postdoc/assigning/revisions/listeria_st/listeria_mlst.txt", 
                            stringsAsFactors=FALSE)
rownames(listeria_mlst) = listeria_mlst$Sample
listeria_mlst = listeria_mlst[,2:8]
listeria_mlst_dists = get_dists(listeria_mlst)

# same (ST)
eval_mlst(listeria_mlst_dists, 0, "listeria_st/mlst_0_clusters.csv", listeria_poppunk)
# 1 diff
eval_mlst(listeria_mlst_dists, 1, "listeria_st/mlst_1_clusters.csv", listeria_poppunk)
# 3 diffs
eval_mlst(listeria_mlst_dists, 3, "listeria_st/mlst_3_clusters.csv", listeria_poppunk)

# E coli MLST
ecoli_poppunk = read.csv("~/Documents/Postdoc/assigning/other_species/e_coli/ecoli_fit_final/ecoli_fit_final_clusters.csv", 
                            stringsAsFactors=FALSE)
#ecoli_poppunk$Taxon = gsub(pattern = 'assemblyfind_E_coli.lane.list/', replacement = "", x = ecoli_poppunk$Taxon)
#ecoli_poppunk$Taxon = gsub(pattern = '.contigs_velvet.fa', replacement = "", x = ecoli_poppunk$Taxon)
#ecoli_poppunk$Taxon = gsub(pattern = '#', replacement = "_", x = ecoli_poppunk$Taxon)

  ecoli_mlst.dists = read.csv("~/Documents/Postdoc/assigning/revisions/ecoli_st/ecoli_mlst.dists.csv", 
                                header=FALSE, stringsAsFactors=FALSE)
  sample_names = paste0("assemblyfind_E_coli.lane.list/", ecoli_mlst.dists[1,-1], ".contigs_velvet.fa")
  ecoli_mlst.dists = ecoli_mlst.dists[-1,-1]
  ecoli_mlst = data.matrix(ecoli_mlst.dists)
  colnames(ecoli_mlst) = sample_names
  rownames(ecoli_mlst) = sample_names
  
  # same (ST)
  eval_mlst(ecoli_mlst, 0, "ecoli_st/mlst_0_clusters.csv", ecoli_poppunk, change=FALSE)
  # 1 diff
  eval_mlst(ecoli_mlst, 1, "ecoli_st/mlst_1_clusters.csv", ecoli_poppunk, change=FALSE)
  # 3 diffs
  eval_mlst(ecoli_mlst, 3, "ecoli_st/mlst_3_clusters.csv", ecoli_poppunk, change=FALSE)

# E coli cgMLST
ecoli_cgmlst.dists = read.csv("~/Documents/Postdoc/assigning/revisions/ecoli_st/ecoli_cgmlst.dists.csv", 
         header=FALSE, stringsAsFactors=FALSE)
sample_names = paste0("assemblyfind_E_coli.lane.list/", ecoli_cgmlst.dists[1,-1])
#sample_names = gsub(pattern = '.contigs.velvet.fa', replacement = "", x = ecoli_cgmlst.dists[1,])
#sample_names = gsub(pattern = '#', replacement = "_", x = sample_names)
ecoli_cgmlst.dists = ecoli_cgmlst.dists[-1,-1]
ecoli_cgmlst = data.matrix(ecoli_cgmlst.dists)
colnames(ecoli_cgmlst) = sample_names
rownames(ecoli_cgmlst) = sample_names

# same (ST)
eval_mlst(ecoli_cgmlst, 0, "ecoli_st/cgmlst_0_clusters.csv", ecoli_poppunk, change=FALSE)
# 1 diff
eval_mlst(ecoli_cgmlst, 1, "ecoli_st/cgmlst_1_clusters.csv", ecoli_poppunk, change=FALSE)
# 3 diffs
eval_mlst(ecoli_cgmlst, 3, "ecoli_st/cgmlst_3_clusters.csv", ecoli_poppunk, change=FALSE)
# 10 diffs
eval_mlst(ecoli_cgmlst, 10, "ecoli_st/cgmlst_10_clusters.csv", ecoli_poppunk, change=FALSE)
# 30 diffs
eval_mlst(ecoli_cgmlst, 30, "ecoli_st/cgmlst_30_clusters.csv", ecoli_poppunk, change=FALSE)
# 100 diffs
eval_mlst(ecoli_cgmlst, 100, "ecoli_st/cgmlst_100_clusters.csv", ecoli_poppunk, change=FALSE)
# 100 diffs
eval_mlst(ecoli_cgmlst, 300, "ecoli_st/cgmlst_300_clusters.csv", ecoli_poppunk, change=FALSE)
