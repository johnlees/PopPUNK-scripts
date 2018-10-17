require(ape)
require(ggplot2)
setwd("~/Documents/Postdoc/assigning/revisions/")

# (modified) Gerry's code for polyphyly
n_non_mono_split <- function(clusters, phy){
  clusters <- split(names(clusters), clusters)
  clusters <- clusters[lapply(clusters, length) > 1]
  pairs <- do.call(rbind, lapply(clusters, function(tips) {n_non_mono_split_c(phy, tips) }))
  pairs_filter = pairs[pairs[,2] > 1,]
  print(sum(pairs_filter[,1])/sum(pairs_filter[,2]))
  return(sum(pairs_filter[,1] > 0))
}

n_non_mono_split_c <- function (phy, tips, reroot = !is.rooted(phy), plot = FALSE, 
                                ...) 
{
  if (!inherits(phy, "phylo")) 
    stop("object 'phy' is not of class 'phylo'")
  n <- length(phy$tip.label)
  if (length(tips) %in% c(1L, n)) 
    return(TRUE)
  ROOT <- n + 1L
  if (is.numeric(tips)) {
    if (any(tips > n)) 
      stop("incorrect tip#: should not be greater than the number of tips")
    tips <- sort(as.integer(tips))
  }
  if (is.character(tips)) 
    tips <- which(phy$tip.label %in% tips)
  if (reroot) {
    outgrp <- phy$tip.label[-tips][1]
    phy <- root(phy, outgroup = outgrp, resolve.root = TRUE)
    rerooted <- TRUE
  }
  else rerooted <- FALSE
  phy <- reorder(phy)
  seq.nod <- .Call(seq_root2tip, phy$edge, n, phy$Nnode)
  sn <- seq.nod[tips]
  newroot <- ROOT
  i <- 2
  repeat {
    x <- unique(unlist(lapply(sn, "[", i)))
    if (length(x) != 1) 
      break
    newroot <- x
    i <- i + 1
  }
  desc <- which(unlist(lapply(seq.nod, function(x) any(x %in% 
                                                         newroot))))
  if (plot) {
    zoom(phy, tips, subtree = FALSE, ...)
    if (rerooted) 
      mtext("Input tree arbitrarily rerooted", side = 1, 
            cex = 0.9)
  }
  
  initial.n.pairs <- choose(length(tips),2)
  
  jumps <- which(diff(tips) != 1)
  sub.sizes <- c(jumps, length(tips)) - c(0, jumps)
  kept.pairs <- sum(unlist(lapply(sub.sizes, choose, 2)))
  missing.pairs <- choose(length(tips), 2) - kept.pairs
  
  # return(missing.pairs)
  return(c(missing.pairs, initial.n.pairs))
  # return(length(jumps))
}

do_compare = function(treefile, pp_clusters, baps_clusters, baps_level = 1)
{
  tr = read.tree(treefile)
  tr$tip.label = gsub(pattern = '#', replacement = "_", x = tr$tip.label)
  
  poppunk = read.csv(pp_clusters, stringsAsFactors=FALSE)
  poppunk$Taxon = gsub(pattern = '.contigs.velvet.fa', replacement = "", x = poppunk$Taxon)
  poppunk$Taxon = gsub(pattern = '#', replacement = "_", x = poppunk$Taxon)
  poppunk$Taxon = sub(pattern = ".*/", replacement = '', x=poppunk$Taxon, perl = TRUE)
  poppunk$Taxon = sub(pattern = ".scaffold.fasta", replacement = '', x=poppunk$Taxon)
  poppunk$Taxon = sub(pattern = "fa", replacement = '', x=poppunk$Taxon)
  
  singleton_clusters = names(which(table(poppunk$Cluster) == 1))
  pp_cluster_assignments = poppunk[!poppunk$Cluster %in% singleton_clusters, 'Cluster']
  names(pp_cluster_assignments) = poppunk[!poppunk$Cluster %in% singleton_clusters, 'Taxon']
  pp_tr = drop.tip(tr,  poppunk[poppunk$Cluster %in% singleton_clusters, 'Taxon'])
  
  print(paste("PopPUNK", length(which(table(poppunk$Cluster) > 1))))
  print(n_non_mono_split(pp_cluster_assignments, pp_tr))
  
  hierbaps_partition <- read.csv(baps_clusters, stringsAsFactors=FALSE)
  baps_cluster_assignments = hierbaps_partition[,baps_level+1]
  names(baps_cluster_assignments) = gsub(pattern = '#', replacement = "_", x = hierbaps_partition$Isolate)
  
  print(paste("hierBAPS", length(which(table(baps_cluster_assignments) > 1))))
  print(n_non_mono_split(baps_cluster_assignments, tr))
}

do_compare("../other_species/staph/core_gene_aln.tree.midpoint.nwk", 
           "~/Documents/Postdoc/assigning/other_species/staph/staph_fit_final/staph_fit_final_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/staph/hierbaps_partition.csv",
           1)

do_compare("../other_species/e_coli/core_gene_aln.midpoint_root.nwk", 
           "~/Documents/Postdoc/assigning/other_species/e_coli/e_coli_db_K3/e_coli_db_K3_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/e_coli/hierbaps_partition.csv",
           1)

do_compare("../other_species/salmonella/microreact_tree.nwk", 
           "~/Documents/Postdoc/assigning/other_species/salmonella/salmonella_db/salmonella_db_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/salmonella/roary_results/hierbaps_partition.csv",
           1)

do_compare("../other_species/listeria/core_gene_tree.midpoint.nwk", 
           "~/Documents/Postdoc/assigning/revisions/listeria_db_revisions/listeria_db_revisions_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/listeria/roary_results/hierbaps_partition.csv",
           2)

do_compare("../other_species/h_flu/core_gene_aln.treefile.midpoint.nwk", 
           "~/Documents/Postdoc/assigning/other_species/h_flu/haem_fit_final/haem_fit_final_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/h_flu/hierbaps_partition.csv",
           1)

do_compare("../other_species/mening/core_gene_aln.treefile.midpoint.nwk", 
           "~/Documents/Postdoc/assigning/other_species/mening/mening_fit_final/mening_fit_final_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/mening/hierbaps_partition.csv",
           1)

do_compare("../other_species/gono_new/core_aln.nwk", 
           "~/Documents/Postdoc/assigning/other_species/gono_new/gono_refine_new/gono_refine_new_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/gono_new/roary_results/hierbaps_partition.csv",
           1)

do_compare("../other_species/s_pyogenes/microreact_ml_tree.nwk", 
           "~/Documents/Postdoc/assigning/other_species/s_pyogenes/pyogenes_refine_final/pyogenes_refine_final_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/s_pyogenes/old_roary/hierbaps_partition.csv",
           1)

do_compare("../SPARC_new/sparc_microreact.nwk", 
           "../SPARC_new/SPARC_refine/SPARC_refine_clusters.csv",
           "../SPARC_new/roary_results/hierbaps_partition.csv",
           1)

do_compare("../other_species/tb/roary_results/tb_tree.midpoint.nwk", 
           "tb_db_revisions/tb_db_revisions_clusters.csv",
           "../other_species/tb/roary_results/hierbaps_partition.csv",
           1)

avg_distances = function(dist_mat, clusters, cluster_type)
{
  #dist_mat = cophenetic.phylo(phy) * length
  # set diag to NA
  diag(dist_mat) = NA
  within_tot = 0
  within_sum = 0
  # iterate over clusters
  within = rep(NA, length(dist_mat))
  within_idx = 1
  for (clust in factor(clusters))
  {
    # sum distances within cluster
    inclust = names(clusters)[clusters == as.numeric(clust)]
    if (length(inclust > 1))
    {
      cls_idx = which(colnames(dist_mat) %in% inclust)
      within_tot = within_tot + sum(dist_mat[cls_idx, cls_idx], na.rm = TRUE)
      within_sum = within_sum + length(cls_idx)
      within_vals = as.numeric(dist_mat[cls_idx, cls_idx])
      # set distances within cluster to NA
      dist_mat[cls_idx, cls_idx] = NA
      
      within[within_idx:(within_idx+length(within_vals)-1)] = within_vals
      within_idx = within_idx + length(within_vals)
    }
  }
  within = within[!is.na(within)]
  
  # average sum
  print(within_tot/within_sum)
  # average remaining non-NA distances
  between = as.numeric(dist_mat[!is.na(dist_mat)])
  print(mean(dist_mat, na.rm = TRUE))
  return(data.frame(Clusters=cluster_type, 
                    Comparison=c(rep("Within", length(within)), rep("Between", length(between))),
                    SNP_distance=c(within, between)))
}

do_distances = function(distfile, pp_clusters, baps_clusters, baps_level, species)
{
  #tr = read.tree(treefile)
  #tr$tip.label = gsub(pattern = '#', replacement = "_", x = tr$tip.label)
  
  #dists = read.table(distfile, sep = ",", header=TRUE)
  dists = read.csv(distfile, header=FALSE, stringsAsFactors=FALSE)
  colnames(dists) = gsub(pattern = '#', replacement = "_", x = dists[1,])
  colnames(dists)[1] = gsub(pattern = '_ ', replacement = "", x = colnames(dists)[1])
  dists = dists[-1,]
  dist_mat = data.matrix(dists)
  rownames(dist_mat) = colnames(dists)
  
  poppunk = read.csv(pp_clusters, stringsAsFactors=FALSE)
  poppunk$Taxon = gsub(pattern = '.contigs.velvet.fa', replacement = "", x = poppunk$Taxon)
  poppunk$Taxon = gsub(pattern = '#', replacement = "_", x = poppunk$Taxon)
  poppunk$Taxon = sub(pattern = ".*/", replacement = '', x=poppunk$Taxon, perl = TRUE)
  poppunk$Taxon = sub(pattern = ".scaffold.fasta", replacement = '', x=poppunk$Taxon)
  poppunk$Taxon = sub(pattern = ".fa", replacement = '', x=poppunk$Taxon)
  pp_cluster_assignments = poppunk$Cluster
  names(pp_cluster_assignments) = poppunk$Taxon
  
  df1 = avg_distances(dist_mat, pp_cluster_assignments, "PopPUNK")

  hierbaps_partition <- read.csv(baps_clusters, stringsAsFactors=FALSE)
  baps_cluster_assignments = hierbaps_partition[,baps_level+1]
  names(baps_cluster_assignments) = gsub(pattern = '#', replacement = "_", x = hierbaps_partition$Isolate)
  
  df2 = avg_distances(dist_mat, baps_cluster_assignments, "RhierBAPS")
  
  ggplot(rbind(df1, df2)) + 
    geom_violin(mapping=aes(x=Clusters, y=SNP_distance, colour=Comparison), scale="width") +
    ggtitle(species) + 
    ylab("SNP distance") + 
    scale_colour_manual(values=c("#999999", "#E69F00")) +
    theme_bw(base_size = 16) +
    theme(plot.title=element_text(face="italic"), legend.position="none")
  ggsave(paste0(species, "_dists.pdf"), device="pdf")
}

# avg distances
#gono_aln_length = 35897
#ecoli_aln_length = 241750
#salmonella_aln_length = 471808
#tb_aln_length = 11871
#pyogenes_aln_length = 80227
#staph_aln_length = 849438 # or 37308 if snps...
#listeria_aln_length = 2035979 # or 143851 if snps...
#hflu_aln_length = 1150000 # or 104261 if snps...
#mening_aln_length = 943285 # or 80593 if snps...
#pneumo_aln_length = 22464

do_distances("staph_dists.csv", 
           "~/Documents/Postdoc/assigning/other_species/staph/staph_fit_final/staph_fit_final_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/staph/hierbaps_partition.csv",
           1,
           "Staphylococcus aureus")

do_distances("ecoli_dists.csv", 
           "~/Documents/Postdoc/assigning/other_species/e_coli/e_coli_db_K3/e_coli_db_K3_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/e_coli/hierbaps_partition.csv",
           1,
           "Escherichia coli")

do_distances("salmonella_dists.csv", 
           "~/Documents/Postdoc/assigning/other_species/salmonella/salmonella_db/salmonella_db_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/salmonella/roary_results/hierbaps_partition.csv",
           1,
           "Salmonella enterica")

do_distances("listeria_dists.csv", 
           "~/Documents/Postdoc/assigning/revisions/listeria_db_revisions/listeria_db_revisions_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/listeria/roary_results/hierbaps_partition.csv",
           2,
           "Listeria monocytogenes")

do_distances("hflu_dists.csv", 
           "~/Documents/Postdoc/assigning/other_species/h_flu/haem_fit_final/haem_fit_final_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/h_flu/hierbaps_partition.csv",
           1,
           "Haemophilus influenzae")

do_distances("mening_dists.csv", 
           "~/Documents/Postdoc/assigning/other_species/mening/mening_fit_final/mening_fit_final_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/mening/hierbaps_partition.csv",
           1,
           "Neisseria meningitidis")

do_distances("gono_dists.csv", 
           "~/Documents/Postdoc/assigning/other_species/gono_new/gono_refine_new/gono_refine_new_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/gono_new/roary_results/hierbaps_partition.csv",
           1,
           "Neisseria gonorrhoeae")

do_distances("pyogenes_dists.csv", 
           "~/Documents/Postdoc/assigning/other_species/s_pyogenes/pyogenes_refine_final/pyogenes_refine_final_clusters.csv",
           "~/Documents/Postdoc/assigning/other_species/s_pyogenes/old_roary/hierbaps_partition.csv",
           1,
           "Streptococcus pyogenes")

do_distances("pneumo_dists.csv", 
           "../SPARC_new/SPARC_refine/SPARC_refine_clusters.csv",
           "../SPARC_new/roary_results/hierbaps_partition.csv",
           1,
           "Streptococcus pneumoniae")

do_distances("tb_dists.csv", 
           "tb_db_revisions/tb_db_revisions_clusters.csv",
           "../other_species/tb/roary_results/hierbaps_partition.csv",
           1,
           "Mycobacterium tuberculosis")
