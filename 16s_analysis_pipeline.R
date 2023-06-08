# edited by Dr. Bin Zhu at Virginia Commonwealth University
##### install packages #####
ip <- as.data.frame(installed.packages())
ip <- ip$Package

if (sum(ip == "rstudioapi") == 0) {
  install.packages("rstudioapi")
}

if (sum(ip == "vegan") == 0) {
  install.packages("vegan")
}

if (sum(ip == "Rtsne") == 0) {
  install.packages("Rtsne")
}

if (sum(ip == "qvalue") == 0) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("qvalue")
}

if (sum(ip == "limma") == 0) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("limma")
}

if (sum(ip == "GUniFrac") == 0) {
  BiocManager::install("GUniFrac")
}

if (sum(ip == "tidyr") == 0) {
  install.packages("tidyr")
}

if (sum(ip == "compositions") == 0) {
  BiocManager::install("compositions")
}

if (sum(ip == "reshape2") == 0) {
  install.packages("reshape2")
}

if (sum(ip == "cp4p") == 0) {
  BiocManager::install("cp4p")
}

if (sum(ip == "ggplot2") == 0) {
  install.packages("ggplot2")
}

if (sum(ip == "Hmisc") == 0) {
  BiocManager::install("Hmisc")

}

if (sum(ip == "Matrix") == 0) {
  install.packages("Matrix")
}

if (sum(ip == "ggpubr") == 0) {
  install.packages("ggpubr")
}

if (sum(ip == "corrplot") == 0) {
  install.packages("corrplot")
}

if (sum(ip == "pheatmap") == 0) {
  BiocManager::install("pheatmap")
}

if (sum(ip == "RColorBrewer") == 0) {
  install.packages("RColorBrewer")
}

if (sum(ip == "igraph") == 0) {
  install.packages("igraph")
}

if (sum(ip == "ALDEx2") == 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("ALDEx2")
}

##### input functions #####
# functions
{
  ######## library ########
  library(vegan) # for diversity
  library(stringr) 
  library(ggplot2)
  library(Rtsne) # for t-SNE
  library(ALDEx2) # for diffential abundance
  library(GUniFrac) # for 'Rarefy'
  library(tidyr) # for 'gather'
  library(compositions) # for 'clr'
  library(reshape2) # for 'melt' in bar plot
  library(cp4p) # for fdr
  library(parallel)
  library(ggpubr)
  # for network
  library(Hmisc)
  library(Matrix)  
  library(corrplot)
  library(pheatmap)
  library(RColorBrewer)
  
  
  
  ######### prepare reads table and metadata #########
  prepare_reads_table = function(reads_table, metadata, total_reads_threshold = 5000, species_threshold = 0.01, mc.cores = 1) {
    
    if (ncol(metadata) == 1) {
      metadata$No_use = NA
    }
    
    # total reads threshold
    keep <- colSums(reads_table) >= total_reads_threshold 
    reads_table <- reads_table[,keep]
    metadata <- metadata[keep,]
    
    # get abundance table
    reads_table_abundance = sweep(reads_table,2,colSums(reads_table),"/")
    
    # species threshold
    keep <- rep(T, nrow(reads_table_abundance))
    
    trials = c(1: nrow(reads_table_abundance))
    func_1 = function(trial) {
      c = sum(reads_table_abundance[trial,] >= species_threshold) / ncol(reads_table_abundance) >= 0.05        # input
      d = sum(reads_table_abundance[trial,] >= species_threshold/10) / ncol(reads_table_abundance) >= 0.15        # input
      keep = c|d
      return(keep)
    }
    keep = mclapply(trials, func_1, mc.cores = mc.cores)
    keep = unlist(keep)
    
    reads_table <- reads_table[keep,]
    
    # total reads threshold
    keep <- colSums(reads_table) >= total_reads_threshold 
    reads_table <- reads_table[,keep]
    metadata <- metadata[keep,]
    metadata$No_use = NULL
    
    return(c(list(reads_table = reads_table),list(metadata = metadata)))
  }
  
  
  ######### alpha diversity ### alpha_diversity ### input samples in cols ###########
  alpha_diversity = function(reads_table, metadata = NA, factor_name = NA, paired = F, order = NA, rarefy_to = NA) {
  
  if (is.na(metadata)[1]) {
    print('no metadata')
    return(NA)
  }
  
  if (is.na(factor_name)) {
    factor_name = ''
  }
  
  # rarefy to normalize data
  reads_table = t(reads_table)
  
  if (is.na(rarefy_to)) {
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  } else {
    reads_table = Rarefy(reads_table, depth = rarefy_to)
  }
  
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  # calculate diversity
  alpha.shannon_diversity <- data.frame(diversity(reads_table))
  alpha.evenness <- alpha.shannon_diversity/log(specnumber(reads_table))
  alpha.ovserved_OTU <- data.frame(colSums(t(reads_table) != 0))
  
  alpha = as.data.frame(matrix(data = NA,ncol=3,nrow = nrow(reads_table)))
  colnames(alpha) = c('alpha.shannon','alpha.evenness','alpha.ovserved_OTU')
  
  alpha$alpha.shannon <- alpha.shannon_diversity$diversity.reads_table.
  alpha$alpha.evenness <- alpha.evenness$diversity.reads_table.
  alpha$alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table.....0.
  
  metadata = cbind(metadata, alpha)
  
  colnames(metadata)[1] = 'factor'
  
  metadata = metadata[order(metadata$factor),]
  
  if (is.na(order)[1] ) {
    metadata$factor <- factor(metadata$factor , levels = unique(metadata$factor))
  } else {
    metadata$factor <- factor(metadata$factor , levels = order)
  }
  
  alpha.shannon = ggplot(metadata, aes(x=factor, y=alpha.shannon)) + geom_violin(trim=T, aes(fill=factor))+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1)+ theme_bw()+
    labs(x = NULL, y = "Shannon index", fill=factor_name)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  
  alpha.evenness = ggplot(metadata, aes(x=factor, y=alpha.evenness)) + geom_violin(trim=T, aes(fill=factor))+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1)+theme_bw()+
    labs(x = NULL, y = "Evenness", fill=factor_name)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  alpha.ovserved_OTU = ggplot(metadata, aes(x=factor, y=alpha.ovserved_OTU)) +geom_violin(trim=T, aes(fill=factor))+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1)+theme_bw()+
    labs(x = NULL, y = "Number of observed taxa", fill=factor_name)+
    geom_jitter(size = 0.1)+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  # calculate significance
  factor_levels = unique(metadata$factor)
  n = length(factor_levels)
  
  Shannon_sig = as.data.frame(matrix(data = NA, nrow =n, ncol = n))
  colnames(Shannon_sig) = factor_levels
  row.names(Shannon_sig) = factor_levels
  
  Evenness_sig = as.data.frame(matrix(data = NA, ncol=n, nrow = n))
  colnames(Evenness_sig) = factor_levels
  row.names(Evenness_sig) = factor_levels
  
  OTU_sig = as.data.frame(matrix(data = NA, ncol=n, nrow = n))
  colnames(OTU_sig) = factor_levels
  row.names(OTU_sig) = factor_levels
  
  for (a in 1:(n-1)) {
    for (b in (a+1) : n) {
      factor_level1 <- subset(metadata,  factor == factor_levels[a],
                              drop = TRUE)
      factor_level2 <- subset(metadata,  factor == factor_levels[b],
                              drop = TRUE)
      
      Shannon_sig[a,b] <- wilcox.test(factor_level1$alpha.shannon, 
                                      factor_level2$alpha.shannon, paired = paired)$p.value
      Evenness_sig[a,b] <- wilcox.test(factor_level1$alpha.evenness, 
                                       factor_level2$alpha.evenness, paired = paired)$p.value
      OTU_sig[a,b] <- wilcox.test(factor_level1$alpha.ovserved_OTU, 
                                  factor_level2$alpha.ovserved_OTU, paired = paired)$p.value
      
    }
  }
  output = c(list(alpha = alpha), list(shannon = alpha.shannon), 
             list(evenness =alpha.evenness) , list(ovserved_OTU =alpha.ovserved_OTU),
             list(sig_Shannon = Shannon_sig),list(sig_Evenness = Evenness_sig),
             list(sig_OTU = OTU_sig))
  
  return(output)
}

  
  
  
  
  
  
  
  
  
  
  ######### beta diversity ### beta_diversity ### input samples in cols; metadata and factor_name are needed; order of factors could be set; ref_group is for setting the reference of bc distance in different groups; can skip from NMDS; output bc distance, within sample distance, distance amoung groups and NMDS #####
  beta_diversity = function(reads_table, metadata = NA, factor_name = NA, order = NA, NMDS_skip = T, 
                            ref_group = NA, rarefy_to = NA, pheatmap_fontsize = 5,treeheight = 50, 
                            pheatmap_y = T, mc.cores = parameter_24) {

    # pretreat data
    {
      if (is.na(metadata)[1]) {
        print('no metadata')
        return(NA)
      }
      
      if (is.na(factor_name)[1]) {
        print('no factor name')
        return(NA)
      }
      
      # rarefy to normalize data
      reads_table = as.data.frame(t(reads_table))
      
      if (is.na(rarefy_to)) {
        reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
      } else {
        reads_table = Rarefy(reads_table, depth = rarefy_to)
      }
      
      reads_table <- reads_table$otu.tab.rff
      reads_table <- as.data.frame(reads_table)
      
      metadata=as.matrix(metadata)
      
      # Bray_Curtis
      Bray_Curtis <- as.matrix(vegdist(reads_table, method = "bray", binary=FALSE))
      Bray_Curtis <- as.data.frame(Bray_Curtis)
      
      Bray_Curtis_2 = Bray_Curtis
      Bray_Curtis_2[row(Bray_Curtis_2) <= col(Bray_Curtis_2)] =NA
      
      # within sample distance
      group_dis = gather(Bray_Curtis_2)
      group_dis$key2 = rep(row.names(Bray_Curtis_2),ncol(Bray_Curtis_2))
      
      Source = matrix(data = NA, ncol = length(metadata), nrow = length(metadata))
      
      for (a in 1:length(metadata)) {
        Source[a,] = metadata
      }
      Source = gather(as.data.frame(Source))
      group_dis$Source = Source$value
      group_dis$Target = rep(metadata,length(metadata))
      
      group_dis = group_dis[!is.na(group_dis$value),]
      group_dis$Source <- as.factor(group_dis$Source)
    }
    
    
    # within-group distance
    {
      keep = group_dis$Source == group_dis$Target
      within_dis = group_dis[keep,]
      keep = within_dis$key != within_dis$key2
      within_dis = within_dis[keep,]
      #  within_dis$value = as.numeric(as.character(within_dis$value))
      
      if (!is.na(order)[1]) {
        within_dis$Source = as.factor(within_dis$Source)
        within_dis$Source = factor(within_dis$Source, levels= order)
      }
      
      within_dis_p = ggplot(within_dis, aes(x=Source, y=value)) + geom_violin(trim=T, aes(fill = Source))+
        geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
        geom_jitter(size = 0.1)+ theme_bw()+
        labs(x = NULL, y = "Within sample distance", fill=factor_name)+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
      
      within_dis_p_2 = ggplot(within_dis, aes(x=Source, y=value)) + geom_violin(trim=T, aes(fill = Source))+
        geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) + theme_bw()+
        labs(x = NULL, y = "Within sample distance", fill=factor_name)+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

      # significance among within sample distance
      
      # Multiple Response Permutation Procedure (MRPP) provides a test of whether there is a significant
      # difference between two or more groups of sampling units. This difference may be one of location
      # (differences in mean) or one of spread (differences in within-group distance; cf. Warton et al. 2012).
      
      group_level = unique(metadata[,1])
      n = length(group_level)
      
      within_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
      colnames(within_dis_sig) = group_level
      row.names(within_dis_sig) = group_level
      for (a in 1:(n-1)) {
        for (b in (a+1): n) {
          keep = metadata[,1] == group_level[a] | metadata[,1] == group_level[b]
          metadata_2 = metadata[keep,1]
          reads_table_2 = reads_table[keep,]
          data = mrpp(reads_table_2, as.matrix(metadata_2), permutations = 999, distance = "bray",
                      weight.type = 1, strata = NULL, parallel = getOption("mc.cores"))
          within_dis_sig[a,b] <- data$Pvalue
          
        }
      }
    }
    
    # distance among groups
    group_level = unique(group_dis$Source)
    n = length(group_level)
    
    {
      keep = group_dis$Source == ref_group
      group_dis_2 = group_dis[keep,]
      
      if (!is.na(order)[1]) {
        group_dis_2$Target = as.factor(group_dis_2$Target)
        group_dis_2$Target = factor(group_dis_2$Target, levels= order)
      }
      
      group_dis_2_p = ggplot(group_dis_2, aes(x=Target, y=value)) +
        geom_boxplot(aes(fill = Target),outlier.shape=NA) +
        labs(x = NULL, y = paste0("Distance to ",ref_group))+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7))
      
      # significance of group sample distance, adonis test
      group_level = unique(metadata)
      n = length(group_level)
      
      group_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
      colnames(group_dis_sig) = group_level
      row.names(group_dis_sig) = group_level
      
      for (a in 1:(n-1)) {
        for (b in (a+1):n) {
          keep = metadata == group_level[a] | metadata == group_level[b]
          metadata_2 = as.character(metadata[keep])
          reads_table_2 = reads_table[keep,]
          
          metadata_2 = as.data.frame(metadata_2)
          pvalue <- adonis2(reads_table_2 ~ ., data = metadata_2, method = "bray", parallel = mc.cores)[1,5] 
          group_dis_sig[a,b] = pvalue
          group_dis_sig[b,a] = pvalue
        }
      }
      
    } 
    
    # Running Nonmetric Multidimensional Scaling (NMDS) Ordination
    if (NMDS_skip == T) {
      # output
      output = c(Bray_Curtis = list(Bray_Curtis), group_dis_2_p = list(group_dis_2_p),
                 group_dis_sig = list(group_dis_sig), within_dis_p = list(within_dis_p), 
                 within_dis_sig = list(within_dis_sig))
      return(output)
      
    } else {
      
      colnames(metadata)[1] = 'factor'
      
      NMDS <-
        metaMDS(Bray_Curtis,
                distance = "bray",
                k = 2,
                maxit = 999, 
                trymax = 20,
                wascores = TRUE)
      
      mds_data <- as.data.frame(NMDS$points)
      mds_data$factor <- metadata
      
      if (!is.na(order)[1]) {
        mds_data$factor = as.factor(mds_data$factor)
        mds_data$factor = factor(mds_data$factor, levels= order)
      }
      
      colnames(mds_data)[1]='NMDS1'
      colnames(mds_data)[2]='NMDS2'
      
      NMDS = ggplot(mds_data, aes(x = NMDS1, y = NMDS2, color = factor)) +
        geom_point(size = 0.3)+
        scale_colour_discrete(factor_name)+
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7))
      
      NMDS_2 = ggplot(mds_data, aes(x =NMDS1, y = NMDS2, color = factor)) +
        geom_point(size = 0.3)+
        scale_colour_discrete(factor_name)+
        stat_ellipse(type = "t")+
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7))
      
      
      # output
      output = c(Bray_Curtis = list(Bray_Curtis), group_dis_2_p = list(group_dis_2_p), within_dis_p = list(within_dis_p), within_dis_p_2 = list(within_dis_p_2),
                 within_dis_sig = list(within_dis_sig), group_dis_sig = list(group_dis_sig), NMDS =list(NMDS), NMDS_2 =list(NMDS_2))
      return(output)
    }
    
  }
  
  
  
  
  
  
  ######### vagitype #######################
  vagitype <- function(reads_table, th = 0.3) {
    # reads_table=reads_table2
    
    reads_table_abundance <- matrix(data =0, ncol = ncol(reads_table),nrow = nrow(reads_table))
    
    for (a in 1:ncol(reads_table)) {
      reads_table_abundance[,a] <- reads_table[,a] / colSums(reads_table)[a]
    }
    reads_table_abundance <- as.data.frame(reads_table_abundance)
    row.names(reads_table_abundance) = row.names(reads_table)
    colnames(reads_table_abundance) = colnames(reads_table)
    
    
    vagitype_1 = matrix(data = NA, ncol =1, nrow = ncol(reads_table_abundance))
    colnames(vagitype_1) = 'Vagitype'
    
    for (a in 1:ncol(reads_table_abundance)) {
      if (max(reads_table_abundance[,a]) < th) {
        vagitype_1[a,1]= 'No type'
      } else {
        n= which(reads_table_abundance[,a] == max(reads_table_abundance[,a]), arr.ind=TRUE)
        vagitype_1[a,1]= row.names(reads_table_abundance)[n][1]
      }
    }
    
    return(vagitype_1)
  }
  
  
  
  
  
  ######### heat map ########
  heatmap_plot <- function(reads_table, metadata, type_th = 0.1, vagitype_num = 5, abundant_taxa_number = 20, width=parameter_15, height=parameter_16) { 
    # reads_table = reads_table_all
    # metadata = metadata_all
    reads_table_abundance <-  sweep(reads_table,2,colSums(reads_table),"/")
    
    # assign vagitype
    reads_table = as.data.frame(t(reads_table_abundance))
    reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
    mytypes <- apply(reads_table,1,which.max)   # find the most abundant taxa
    maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes), ncol=2)]  # the abundance of the most abundant taxa
    mytypes <- colnames(reads_table)[mytypes]   # find the name of the most abundant taxa
    mytypes[maxprop < type_th] <- "No Type"
    uniqueTypes <- unique(mytypes)
    top_species = as.data.frame(sort(colSums(reads_table),decreasing=T))
    keep = mytypes %in% row.names(top_species)[1:vagitype_num]
    mytypes[!keep] = 'Others'
    
    metadata$Vagitype = mytypes
    row.names(metadata) = row.names(reads_table)
    #colnames(metadata)[1] = metadata_name
    
    top_abundant_taxa = rowSums(reads_table_abundance)
    top_abundant_taxa = as.data.frame(sort(top_abundant_taxa, decreasing = T))
    top_abundant_taxa = row.names(top_abundant_taxa)[1:abundant_taxa_number]
    reads_table_abundance = reads_table_abundance[row.names(reads_table_abundance) %in% top_abundant_taxa, ]
    
    save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
      stopifnot(!missing(x))
      stopifnot(!missing(filename))
      pdf(filename, width=width, height=height)
      grid::grid.newpage()
      grid::grid.draw(x$gtable)
      dev.off()
    }
    
    heatmap_plot = pheatmap(as.matrix(reads_table_abundance),cluster_rows = T, cluster_cols =T, clustering_distance_cols = "euclidean", clustering_method = "complete",
             annotation_col = metadata ,show_colnames = F,treeheight_col = 50)
    save_pheatmap_pdf(heatmap_plot, "Heatmap_clustering.pdf", width=parameter_15, height=parameter_16)
    
    
  }
  
  ######### network rcorr ############# sample in cols #######
  newwork_rcorr <- function(reads_table, correlation = "spearman",  pvalue = 0.05, cor_parameter= 0, style = 1, bar_max = 2, bar_min = -2, pheatmap_fontsize = 5, treeheight = 50, alpha = 0.05) {
    #  reads_table = reads_table_2
    
    reads_table = as.data.frame(t(reads_table))
    reads_table = reads_table + 0.5
    reads_table <- clr(reads_table)      ### CLR normalization in rows
    
    # pvalue is the pvalue threshold for cooccurence
    # cor_parameter is the cooccurence value threshold
    # style is the style of the heatmap
    
    #reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
    #reads_table <- reads_table$otu.tab.rff
    #reads_table <- as.data.frame(reads_table)
    
    reads_table = as.matrix((reads_table))
    
    otu.cor <- rcorr(reads_table, type=correlation)
    
    otu.pval <- forceSymmetric(otu.cor$P)
    otu.pval <- otu.pval@x
    otu.pval <- adjust.p(otu.pval, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
    otu.pval <- otu.pval$adjp
    otu.pval <- otu.pval$adjusted.p
    
    p.yes <- otu.pval< pvalue  
    
    r.val = otu.cor$r # select all the correlation values 
    p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion 
    
    p.yes.r <- abs(p.yes.r) > cor_parameter # output is logical vector
    p.yes.rr <- p.yes.r*r.val # use logical vector for subscripting.
    
    p.yes.rr[is.na(p.yes.rr)] = 0
    
    keep = abs(colSums(p.yes.rr)) > 0
    p.yes.rr = p.yes.rr[,keep]
    p.yes.rr = p.yes.rr[keep,]
    
    # gephi output
    gephi_p.yes.rr = as.matrix(p.yes.rr)
    for (a in 1:nrow(p.yes.rr)) {
      for (b in 1:ncol(p.yes.rr)) {
        if (a >=b) {
          gephi_p.yes.rr[a,b] = NA
        }
      }
    }
    gephi_p.yes.rr <- as.data.frame(gephi_p.yes.rr)
    gephi_p.yes.rr = gather(gephi_p.yes.rr)
    gephi_p.yes.rr$Taxa = rep(row.names(p.yes.rr),ncol(p.yes.rr))
    gephi_p.yes.rr = gephi_p.yes.rr[!is.na(gephi_p.yes.rr[,2]),]
    gephi_p.yes.rr$value = as.numeric(as.character(gephi_p.yes.rr$value))
    gephi_p.yes.rr = gephi_p.yes.rr[gephi_p.yes.rr$value != 0,] 
    gephi_p.yes.rr = gephi_p.yes.rr[gephi_p.yes.rr[,2] >= cor_parameter | gephi_p.yes.rr[,2] <= -cor_parameter,]
    colnames(gephi_p.yes.rr) = c('Source','Weight','Target')
    # heatmap
    library(igraph)
    
    if (style == 1) {
      if (bar_max == 2) {
        bar_max = max(p.yes.rr,na.rm = T)
        bar_min = min(p.yes.rr,na.rm = T)
        paletteLength <- 50
        myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
        myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                      seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
        
        myBreaks <- unique(myBreaks)
        
        if (sum(myBreaks < 0) == 0) {
          myColor <- colorRampPalette(c("white", "red"))(paletteLength)
          p <- pheatmap(p.yes.rr, color=myColor, fontsize = pheatmap_fontsize, 
                        treeheight_row = treeheight, treeheight_col = treeheight)
        } else if (sum(myBreaks > 0) == 0) {
          myColor <- colorRampPalette(c("blue", "white"))(paletteLength)
          p <- pheatmap(p.yes.rr, color=myColor, fontsize = pheatmap_fontsize, 
                        treeheight_row = treeheight, treeheight_col = treeheight)
        } else {
          p <- pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                        treeheight_row = treeheight, treeheight_col = treeheight)
        }
      } else {
        paletteLength <- 50
        myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
        myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                      seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
        
        if (sum(myBreaks < 0) == 0) {
          myColor <- colorRampPalette(c("white", "red"))(paletteLength)
          p <- pheatmap(p.yes.rr, color=myColor, fontsize = pheatmap_fontsize, 
                        treeheight_row = treeheight, treeheight_col = treeheight)
        } else if (sum(myBreaks > 0) == 0) {
          myColor <- colorRampPalette(c("blue", "white"))(paletteLength)
          p <- pheatmap(p.yes.rr, color=myColor, fontsize = pheatmap_fontsize, 
                        treeheight_row = treeheight, treeheight_col = treeheight)
        } else {
          p <- pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                        treeheight_row = treeheight, treeheight_col = treeheight)
        }
      }
      
    } else {
      p.yes.rr <- as.matrix(p.yes.rr)
      p = corrplot(p.yes.rr, type="upper", order="hclust", insig = "blank",tl.col = "black",na.label = "o")
    }
    
    detach("package:igraph", unload = TRUE, force = T)
    return(c(list(p = p), list(gephi_input = gephi_p.yes.rr), list(cor_matrix = p.yes.rr)))
  }
  
  
  
  
  
  
  ######### diffential abundance ### dif_abundance; sample in column; paired text or not; change order of species; fold change threshold; if style = 1, give dot plot, else give box plot ####################
  dif_abundance <- function(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, paired_test = F, order_reverse = F, style = 1, order = NA) {
    #   style =1
    #   reads_table= reads_table
    #   metadata= metadata$Flag
    #   paired_test = F
    #   order_reverse = F
    #   style =2 
    #   order = c('TB','PTB')
    #   fold_change_th = 1
    #   pvalue_th = 0.05
    
    conds <- metadata
    
    x <- aldex.clr(reads_table, conds, mc.samples=128, denom="all", verbose=F)
    
    # paired Wilcoxon Rank Sum test and Welch's t-test
    x.tt <- aldex.ttest(x, paired.test= paired_test)
    
    x.effect <- aldex.effect(x)
    
    x.all <- data.frame(cbind(x.tt,x.effect))
    x.all$`-Log10(adj-pvalue)` <- -log10(x.all$wi.eBH)
    x.all$Taxa = row.names(x.all)
    
    if (max(abs(x.all$diff.btw)) < fold_change_th | min(x.tt$wi.eBH) > 0.05) {
      print('No taxon has significant abundance change')
      return(list(data = x.all))
      
    }
    
    if (order_reverse == T) {
      x.all$diff.btw = -x.all$diff.btw
    }
    
    # draw figure
    das <- x.all[(x.all$`-Log10(adj-pvalue)` >= -log10(pvalue_th) & (x.all$diff.btw >=fold_change_th | x.all$diff.btw <=-fold_change_th)),]
    
    if (nrow(das)==0) {
      print('No taxon has significant abundance change')
      return(list(data = x.all))
    }
    
    das$Species <- row.names(das)
    das <- das[order(das$diff.btw),] 
    
    metadata = as.factor(metadata)
    lev = levels(metadata)
    
    if (order_reverse == T) {
      das$Color <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[2]), paste0("Enriched in ",lev[1]))  # above / below avg flag
    } else {
      das$Color <- ifelse(das$diff.btw < 0, paste0("Enriched in ",lev[1]), paste0("Enriched in ",lev[2]))  # above / below avg flag
      
    }
    
    das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.
    
    das$diff.btw[das$diff.btw == Inf] = 10
    das$diff.btw[das$diff.btw == -Inf] = -10
    das$diff.btw[das$diff.btw <= -10] = -10
    das$diff.btw[das$diff.btw >= 10] = 10
    
    if (!is.na(order)[1]) {
      if (sum(sort(unique(metadata)) == order) == 2) {
        das$diff.btw = -das$diff.btw
      }
    }
    
    if (style == 1) {
      theme_set(theme_bw())  
      
      p <- ggplot(das, aes(Species, diff.btw)) + 
        geom_point(aes(col=Color, size=`-Log10(adj-pvalue)`)) + 
        coord_flip() +          # convert x y axis
        labs(x = 'Taxa', y = "Median difference in clr values")+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7))  
      
      
    } else {
      taxa_list = row.names(das)
      
      keep = which(row.names(reads_table) %in% taxa_list)
      
      reads_table_abundance = sweep(reads_table,2,colSums(reads_table),"/")
      reads_table_2 <- reads_table_abundance[keep,]
      reads_table_2 <- as.data.frame(t(reads_table_2))
      
      #   colnames(reads_table_2)=taxa_list
      reads_table_3 = gather(reads_table_2)
      reads_table_3$Type = rep(conds, length(taxa_list))
      colnames(reads_table_3) = c('Taxa', 'Abundance','Type')
      
      reads_table_3$Taxa = str_replace_all(reads_table_3$Taxa,'.*__','')
      
      if (!is.na(order)[1]) {
        reads_table_3$Type = factor(reads_table_3$Type, levels = order)
      }
      
      p <- ggplot(reads_table_3, aes(x=Taxa, y=Abundance,fill=Type)) +
        geom_boxplot(outlier.shape = NA)+
        ylab("Abundance (%)")+
        theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1))
    }
    
    return(c(list(p = p), list(data = x.all)))
    
  }
  
  
  ##### One_factor_script #####
  One_factor_script = function(reads_table, metadata, ref_group, order, VMB = parameter_28) {
    factor_1 = colnames(metadata)[1]
    # reads table filter
    data = prepare_reads_table(reads_table, metadata, total_reads_threshold = parameter_21, species_threshold = parameter_22) 
    reads_table = data$reads_table
    metadata = data$metadata
    
    # alpha diversity
    {
      alpha = alpha_diversity(reads_table, metadata = metadata, factor_name = factor_1, paired = F,order = order)
      alpha$shannon 
#      ggsave('alpha_shannon.pdf',width=parameter_1, height=parameter_2)
#      alpha$evenness  
#      ggsave('alpha_evenness.pdf',width=parameter_1, height=parameter_2)
#      alpha$ovserved_OTU 
#      ggsave('alpha_ovserved_OTU.pdf',width=parameter_1, height=parameter_2)
      
      ggarrange(alpha$shannon, alpha$evenness, alpha$ovserved_OTU,
                          labels = c(" ", " ", " "),
                          ncol = 3, nrow = 1)
      ggsave('alpha.pdf',width=parameter_1*3, height=parameter_2)
      
      
      write.csv(alpha$sig_Shannon,'alpha_shannon.csv', row.names = T, quote = F)
      write.csv(alpha$sig_Evenness,'alpha_evenness.csv', row.names = T, quote = F)
      write.csv(alpha$sig_OTU,'alpha_ovserved_OTU.csv', row.names = T, quote = F)
    }
    
    # beta diversity
    {
      beta = beta_diversity(reads_table, metadata = as.character(metadata[,1]), factor_name = factor_1, 
                            order = order, NMDS_skip = F, ref_group = ref_group, 
                            rarefy_to = NA, pheatmap_fontsize = 5, treeheight = 5, 
                            pheatmap_y = T)
      
      beta$NMDS_2
      ggsave('NMDS.pdf',width=parameter_3, height=parameter_4)
      write.csv(beta$group_dis_sig, 'NMDS_Adonis_test.csv')
      
      pdf('Within_group_distance.pdf',width=parameter_1, height=parameter_2) 
      print(beta$within_dis_p)
      dev.off()
      
      pdf('Within_group_distance_2.pdf',width=parameter_1, height=parameter_2) 
      print(beta$within_dis_p_2)
      dev.off()
      write.csv(beta$within_dis_sig, 'Within_group_distance_MRPP_test.csv')
    }

    metadata_level = length(unique(metadata[,1]))
    metadata_level_2 = unique(metadata[,1])
    
    # composition bar plot
    {
      type_th = parameter_14; taxa_num = parameter_13; VMB = parameter_28
      reads_table_2 = reads_table; metadata_2 = metadata
        if (VMB == F) {
          colnames(metadata_2)[1] = 'factor_1'
          
          # convert reads to relative abundance
          {
            reads_table_2_abundance <- sweep(reads_table_2,2,colSums(reads_table_2),"/")
            mypropdata = as.data.frame(t(reads_table_2_abundance))
            mypropdata <- mypropdata[,apply(mypropdata,2,sum) > 0]
          }
          
          # assign vagitype
          {
            mytypes <- apply(mypropdata,1,which.max)
            maxprop <- mypropdata[matrix(c(1:nrow(mypropdata),mytypes), ncol=2)]
            mytypes <- colnames(mypropdata)[mytypes]
            mytypes[maxprop < type_th] <- "No Type"
            
            top_mytypes = as.data.frame(table(mytypes))
            top_mytypes = top_mytypes[order(top_mytypes$Freq, decreasing = T),]
            keep = mytypes %in% top_mytypes$mytypes[1:taxa_num]
            mytypes[!keep] = 'Others'
          }
          
          # set type order
          {
            uniqueTypes <- unique(mytypes)
            
            myTypeOrder <- uniqueTypes
            if (length(grep("No Type", myTypeOrder)) > 0) {
              myTypeOrder <- c(myTypeOrder[-grep("No Type", myTypeOrder)], "No Type")
            }
            if (length(grep("Others", myTypeOrder)) > 0) {
              myTypeOrder <- c(myTypeOrder[-grep("Others", myTypeOrder)], "Others")
            }
            myTypeOrder <- data.frame(typeOrder=c(1:length(myTypeOrder)),row.names=myTypeOrder)
            
            mypropdata <- 100*mypropdata
          }
          
          # order samples by type then proportion
          {
            mypropdata <- mypropdata[order(myTypeOrder[mytypes,],-maxprop),]  # myTypeOrder[mytypes,]: find the order number of mytypes in the myTypeOrder list    # order(a,b): order a then order b
            barplottypes <- mytypes[order(myTypeOrder[mytypes,],-maxprop)]
            
            x = metadata_2; x$X =NA
            x = x[match(row.names(mypropdata), row.names(x)),]
            x1 = order(x$factor_1)
            
            mypropdata = mypropdata[x1,]
            x = x[x1,]
            
            row.names(mypropdata) = paste(x$factor_1, row.names(mypropdata), sep='_')
            barplottypes = barplottypes[x1]
          }

          # set color
          {
            mycolors <- rep(NA, ncol(mypropdata))
            keep = as.data.frame(colSums(mypropdata))
            keep$x = NA
            keep = keep[order(keep$`colSums(mypropdata)`, decreasing = T),]
            taxa_on_legend = row.names(keep)[1:taxa_num]
            keep = colnames(mypropdata) %in% taxa_on_legend
            mycolors[!keep] = "gray"
            
            n = sum(is.na(mycolors))
            
            library(RColorBrewer)
            qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
            col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

            mycolors_2 = col_vector[1:n]
            mycolors[is.na(mycolors)] = mycolors_2
          } 
          
          # ggplot
          {
            library(ggplot2)
            library(reshape2) # for 'melt' in bar plot
            #p
            mycolors <- factor(mycolors)
            mycolordf <- data.frame(mycolors, Taxa=colnames(mypropdata))
            sampleorder <- factor(rownames(mypropdata))
            
            ggplotdata <- melt(as.matrix(mypropdata),  
                               id.vars=names(mypropdata), 
                               varnames=c("SampleID", "Taxa"), 
                               value.name="ATprop")
            ggplotdata$SampleID <- factor(ggplotdata$SampleID, levels=sampleorder)
            ggplotdata <- merge(ggplotdata, mycolordf, by="Taxa")   # add color information by searching taxa
            ggplotdata <- ggplotdata[ggplotdata$ATprop != 0,]
            
            p1 = ggplot(ggplotdata, aes(SampleID, ATprop, fill=mycolors, group=ATprop)) + 
              geom_bar(stat="identity", position="stack", width=1) +
              scale_fill_manual(values = levels(mycolors)) +
              labs(x="Sample", y="Relative Abundance") + 
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                    axis.ticks=element_blank())  
          } 

          # legend figure
          {
            mycolordf_2 = mycolordf[mycolordf$Taxa %in% taxa_on_legend,]
            mycolordf_2 = mycolordf_2[!duplicated(mycolordf_2$mycolors),]
            mycolordf_2$Taxa[mycolordf_2$mycolors == 'gray'] = 'Others'
            mycolordf_2$value =1
            
            taxa_on_legend = data.frame(taxa_on_legend = taxa_on_legend, order = c(1:length(taxa_on_legend)))
            
            mycolordf_2$order = NA
            for (a in 1:nrow(mycolordf_2)) {
              if (mycolordf_2$Taxa[a] == 'Others') {
                mycolordf_2$order[a] = nrow(taxa_on_legend) +1
              } else {
                mycolordf_2$order[a] = taxa_on_legend$order[which(taxa_on_legend$taxa_on_legend == mycolordf_2$Taxa[a])] 
              }
            }
            mycolordf_2 = mycolordf_2[order(mycolordf_2$order),]
            
            mycolordf_2$Taxa = factor(mycolordf_2$Taxa, levels = mycolordf_2$Taxa)
            mycolordf_2$mycolors = as.character(mycolordf_2$mycolors)
            mycolordf_2[nrow(mycolordf_2) +1,] = NA
            mycolordf_2$mycolors[nrow(mycolordf_2)] = 'gray'
            
            mycolordf_2$Taxa = as.character(mycolordf_2$Taxa)
            mycolordf_2$Taxa[nrow(mycolordf_2)] = 'Others'
            mycolordf_2$Taxa = factor(mycolordf_2$Taxa, levels = mycolordf_2$Taxa)
            mycolordf_2$value[nrow(mycolordf_2)] = 1
            mycolordf_2$order[nrow(mycolordf_2)] = nrow(mycolordf_2)
            
            p2 = ggplot(mycolordf_2, aes(x = Taxa, y = value, fill = Taxa)) +
              geom_bar(stat = "identity")+
              scale_fill_manual(values= mycolordf_2$mycolors)
          }

          pdf(paste0('bar_plot_1.pdf') ,width=parameter_5, height=parameter_6) 
          print(p1)
          dev.off()
          
          pdf(paste0('bar_plot_2.pdf') ,width=parameter_5, height=8) 
          print(p2)
          dev.off()
        } else {
          colnames(metadata_2)[1] = 'factor_1'
          
          # convert reads to relative abundance
          {
            reads_table_2_abundance <- sweep(reads_table_2,2,colSums(reads_table_2),"/")
            mypropdata = as.data.frame(t(reads_table_2_abundance))
            mypropdata <- mypropdata[,apply(mypropdata,2,sum) > 0]
          }
          
          # assign vagitype
          {
            mytypes <- apply(mypropdata,1,which.max)
            maxprop <- mypropdata[matrix(c(1:nrow(mypropdata),mytypes), ncol=2)]
            mytypes <- colnames(mypropdata)[mytypes]
            mytypes[maxprop < type_th] <- "No Type"
            
            top_mytypes = as.data.frame(table(mytypes))
            top_mytypes = top_mytypes[order(top_mytypes$Freq, decreasing = T),]
            keep = mytypes %in% top_mytypes$mytypes[1:taxa_num]
            mytypes[!keep] = 'Others'
          }
          
          # set type order
          {
            uniqueTypes <- unique(mytypes)
            lactoTypes <- c(
              grep("crispatus", uniqueTypes, value=TRUE),
              grep("jensenii", uniqueTypes, value=TRUE), 
              grep("gasseri", uniqueTypes, value=TRUE), 
              grep("iners", uniqueTypes, value=TRUE)
            ) 
            
            nonLactoTypes <- c(
              grep("Gardnerella", uniqueTypes, value=TRUE), 
              grep("BVAB", uniqueTypes, value=TRUE), 
              grep("Sneathia", uniqueTypes, value=TRUE)
            )
            
            myTypeOrder <- c(lactoTypes, nonLactoTypes)
            myTypeOrder <- c(myTypeOrder, setdiff(uniqueTypes,myTypeOrder))
            if (length(grep("No Type", myTypeOrder)) > 0) {
              myTypeOrder <- c(myTypeOrder[-grep("No Type", myTypeOrder)], "No Type")
            }
            if (length(grep("Others", myTypeOrder)) > 0) {
              myTypeOrder <- c(myTypeOrder[-grep("Others", myTypeOrder)], "Others")
            }
            myTypeOrder <- data.frame(typeOrder=c(1:length(myTypeOrder)),row.names=myTypeOrder)
            
            mypropdata <- 100*mypropdata
          }
          
          # order samples by type then proportion
          {
            mypropdata <- mypropdata[order(myTypeOrder[mytypes,],-maxprop),]  # myTypeOrder[mytypes,]: find the order number of mytypes in the myTypeOrder list    # order(a,b): order a then order b
            barplottypes <- mytypes[order(myTypeOrder[mytypes,],-maxprop)]
            
            x = metadata_2; x$X =NA
            x = x[match(row.names(mypropdata), row.names(x)),]
            x1 = order(x$factor_1)
            
            mypropdata = mypropdata[x1,]
            x = x[x1,]
            
            row.names(mypropdata) = paste(x$factor_1, row.names(mypropdata), sep='_')
            mypropdata = mypropdata[,]
            
            barplottypes = barplottypes[x1]
          }
          
          
          # set color
          {
            
            mycolors <- rep("gray", ncol(mypropdata))
            # for genus-level
            mycolors[grep("Lactobacillus",   colnames(mypropdata))]      <- "gray"
            mycolors[grep("Sneathia",   colnames(mypropdata))]           <- "#9467bd"
            mycolors[grep("Gardnerella",     colnames(mypropdata))]      <- "#d62728"
            mycolors[grep("Lachnospiraceae", colnames(mypropdata))]      <- "#ff7f0e"
            mycolors[grep("BVAB", colnames(mypropdata))]                 <- "#ff7f0e"
            mycolors[grep("Prevotella",        colnames(mypropdata))]    <- "#1f77b4"
            mycolors[grep("Atopobium",        colnames(mypropdata))]     <- "#c49c94"
            mycolors[grep("Fannyhessea",        colnames(mypropdata))]     <- "#c49c94"
            mycolors[grep("Streptococcus",        colnames(mypropdata))]     <- "#9edae5"
            mycolors[grep("Dialister",        colnames(mypropdata))]     <- "#dbdb8d"
            mycolors[grep("Megasphaera",        colnames(mypropdata))]     <- "#f7b6d2"
            mycolors[grep("Coriobacteriaceae",        colnames(mypropdata))]     <- "#bdc2cc"
            mycolors[grep("Ureaplasma",        colnames(mypropdata))]     <- "#2ca02c"
            mycolors[grep("Mobiluncus",        colnames(mypropdata))]     <- "#c4c392"
            mycolors[grep("TM7",        colnames(mypropdata))]     <- "#e377c2"
            mycolors[grep("Fusobacterium",        colnames(mypropdata))]     <- "#c8d045"
            mycolors[grep("Aerococcus",        colnames(mypropdata))]     <- "#a5acaf"
            mycolors[grep("Coriobacteriales",        colnames(mypropdata))]     <- "#698838"
            mycolors[grep("Tissierellia",        colnames(mypropdata))]     <- "#ce4946"
            mycolors[grep("Veillonellaceae",        colnames(mypropdata))]     <- "#a37b5d"
            mycolors[grep("Gemella",        colnames(mypropdata))]     <- "#d59693"
            mycolors[grep("Winkia",        colnames(mypropdata))]     <- "#c84b9d"
            mycolors[grep("Peptoniphilus",        colnames(mypropdata))]     <- "#d06465"
            
            mycolors[grep("crispatus",   colnames(mypropdata))]           <- "#fff5aa"
            mycolors[grep("iners",   colnames(mypropdata))]               <- "#aec7e8"
            mycolors[grep("gasseri",   colnames(mypropdata))]             <- "#c5b0d5"
            mycolors[grep("delbrueckii",   colnames(mypropdata))]         <- "black"
            mycolors[grep("jensenii",   colnames(mypropdata))]            <- "#ffbb78"
            mycolors[grep("sanguinegens",        colnames(mypropdata))] <- "#ff9896"
            mycolors[grep("Streptococcus",     colnames(mypropdata))]     <- "orange"
            mycolors[grep("Prevotella amnii",     colnames(mypropdata))]     <- "#c7c7c7"
            mycolors[grep("Prevotella_amnii",     colnames(mypropdata))]     <- "#c7c7c7"
            mycolors[grep("Prevotella bivia",     colnames(mypropdata))]     <- "#17becf"
            mycolors[grep("Prevotella_bivia",     colnames(mypropdata))]     <- "#17becf"
            
            
            keep = as.data.frame(colSums(mypropdata))
            keep$x = NA
            keep = keep[order(keep$`colSums(mypropdata)`, decreasing = T),]
            taxa_on_legend = row.names(keep)[1:taxa_num]
            keep = colnames(mypropdata) %in% taxa_on_legend
            mycolors[!keep] = "gray"
          } 
          
          # ggplot
          {
            library(ggplot2)
            library(reshape2) # for 'melt' in bar plot
            #p
            mycolors <- factor(mycolors)
            mycolordf <- data.frame(mycolors, Taxa=colnames(mypropdata))
            sampleorder <- factor(rownames(mypropdata))
            
            ggplotdata <- melt(as.matrix(mypropdata),  
                               id.vars=names(mypropdata), 
                               varnames=c("SampleID", "Taxa"), 
                               value.name="ATprop")
            ggplotdata$SampleID <- factor(ggplotdata$SampleID, levels=sampleorder)
            ggplotdata <- merge(ggplotdata, mycolordf, by="Taxa")   # add color information by searching taxa
            ggplotdata <- ggplotdata[ggplotdata$ATprop != 0,]
            
            p1 = ggplot(ggplotdata, aes(SampleID, ATprop, fill=mycolors, group=ATprop)) + 
              geom_bar(stat="identity", position="stack", width=1) +
              scale_fill_manual(values = levels(mycolors)) +
              labs(x="Sample", y="Relative Abundance") + 
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                    axis.ticks=element_blank())  
          } 
          
          # legend figure
          {
            mycolordf_2 = mycolordf[mycolordf$Taxa %in% taxa_on_legend,]
            mycolordf_2 = mycolordf_2[!duplicated(mycolordf_2$mycolors),]
            mycolordf_2$Taxa[mycolordf_2$mycolors == 'gray'] = 'Others'
            mycolordf_2$value =1
            
            taxa_on_legend = data.frame(taxa_on_legend = taxa_on_legend, order = c(1:length(taxa_on_legend)))
            n = which(str_detect(taxa_on_legend$taxa_on_legend,"crispatus" ))
            if (length(n)!=0){taxa_on_legend$order[n]=-100}
            n = which(str_detect(taxa_on_legend$taxa_on_legend,"jensenii" ))
            if (length(n)!=0){taxa_on_legend$order[n]=-99}
            n = which(str_detect(taxa_on_legend$taxa_on_legend,"gasseri" ))
            if (length(n)!=0){taxa_on_legend$order[n]=-98}
            n = which(str_detect(taxa_on_legend$taxa_on_legend,"iners" ))
            if (length(n)!=0){taxa_on_legend$order[n]=-97}
            
            n = which(str_detect(taxa_on_legend$taxa_on_legend,"Gardnerella" ))
            if (length(n)!=0){taxa_on_legend$order[n]=-10}
            n = which(str_detect(taxa_on_legend$taxa_on_legend,"BVAB" ))
            if (length(n)!=0){taxa_on_legend$order[n]=-9}
            n = which(str_detect(taxa_on_legend$taxa_on_legend,"Sneathia" ))
            if (length(n)!=0){taxa_on_legend$order[n]=-8}
            
            
            mycolordf_2$order = NA
            for (a in 1:nrow(mycolordf_2)) {
              if (mycolordf_2$Taxa[a] == 'Others') {
                mycolordf_2$order[a] = nrow(taxa_on_legend) +1
              } else {
                mycolordf_2$order[a] = taxa_on_legend$order[which(taxa_on_legend$taxa_on_legend == mycolordf_2$Taxa[a])] 
              }
            }
            mycolordf_2 = mycolordf_2[order(mycolordf_2$order),]
            
            mycolordf_2$Taxa = factor(mycolordf_2$Taxa, levels = mycolordf_2$Taxa)
            mycolordf_2$mycolors = as.character(mycolordf_2$mycolors)
            
            p2 = ggplot(mycolordf_2, aes(x = Taxa, y = value, fill = Taxa)) +
              geom_bar(stat = "identity")+
              scale_fill_manual(values= mycolordf_2$mycolors)
          }
          
          output = c(p1 = list(p1), p2 = list(p2))
          
          pdf(paste0('bar_plot_1.pdf') ,width=parameter_5, height=parameter_6) 
          print(output$p1)
          dev.off()
          
          pdf(paste0('bar_plot_2.pdf') ,width=parameter_5, height=8) 
          print(output$p2)
          dev.off()
        }
        
      }

    # abundance difference
    if (metadata_level == 2) {
      abundance_result = dif_abundance(reads_table,as.character(metadata[,1]), 
                                       pvalue_th = parameter_25, fold_change_th = parameter_26, paired_test = F, 
                                       order_reverse = F, style = 1, 
                                       order = order)
      write.csv(abundance_result$data,'abundance_difference.csv', row.names = T, quote = F)
      
      if (length(abundance_result) == 2) {
        abundance_result$p
        ggsave('abundance_difference.pdf',width=parameter_7, height=parameter_8)

      }
    } else {
      for (x in 1:(metadata_level-1)) {
        factor_levels = unique(metadata[,1])
        factor_levels = as.character(factor_levels[-which(factor_levels == ref_group)])
        
        keep = metadata[,1] %in% c(ref_group, factor_levels[x])
        metadata_2 = metadata[keep,]
        reads_table_2 = reads_table[,keep]
        abundance_result = dif_abundance(reads_table_2, as.character(metadata_2), 
                                         pvalue_th = parameter_25, fold_change_th = parameter_26, paired_test = F, 
                                         order_reverse = F, style = 1, 
                                         order = order[c(1,x+1)])
        
        write.csv(abundance_result$data,paste0('abundance_difference_',factor_levels[x],'.csv'), row.names = T, quote = F)
        if (length(abundance_result) == 2) {
          abundance_result$p
          ggsave(paste0('abundance_difference_',factor_levels[x],'.pdf'),width=parameter_7, height=parameter_8)
          
        }
      }
      
    }
    
    # heatmap
    {
      heatmap_plot(reads_table, metadata, type_th = 0.1, vagitype_num = 5, abundant_taxa_number = 20, width=18, height=8)
    }
    
    
    # network
    {
      network_result = newwork_rcorr(reads_table, correlation = parameter_20, pvalue = 0.05, 
                                     cor_parameter= 0, style = 1, bar_max = 2, bar_min = -2, 
                                     pheatmap_fontsize = parameter_18, treeheight = parameter_19, alpha = 0.05)
      save_pheatmap_pdf <- function(x, filename, width=8, height=8) {
        stopifnot(!missing(x))
        stopifnot(!missing(filename))
        pdf(filename, width=width, height=height)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
      }
      save_pheatmap_pdf(network_result$p, "network.pdf",width=parameter_9, height=parameter_10)
      
      # network of abundant taxa
      data = prepare_reads_table(reads_table, metadata, total_reads_threshold = parameter_21, species_threshold = parameter_23, mc.cores = parameter_24) 
      reads_table_2 = data$reads_table
      
      network_result = newwork_rcorr(reads_table_2, correlation = parameter_20, pvalue = 0.05,
                                     cor_parameter= 0, style = 1, bar_max = 2, bar_min = -2, 
                                     pheatmap_fontsize = parameter_18, treeheight = parameter_19, alpha = 0.05)
      save_pheatmap_pdf(network_result$p, "network_abundant_taxa.pdf",width=parameter_11, height=parameter_12)
      
    }

  }
}

##### input data #####
{
  library(rstudioapi)
  path_input = dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(path_input)
  
  reads_table_all = read.csv('reads_table.csv', row.names = 1, header = T)  ### input reads table file name
  metadata_all = read.csv('metadata.csv', row.names = 1, header = T)  ### input reads table file name
  
  if (sum(colnames(reads_table_all) %in% row.names(metadata_all)) == 0) {
    reads_table_all = as.data.frame(t(reads_table_all))
  }
  
  if (sum(colnames(reads_table_all) %in% row.names(metadata_all)) == 0) {
    print('Sample_ID in the metadata is not the same as that in the reads table')
  }
  
  keep = colnames(reads_table_all) %in% row.names(metadata_all)
  reads_table_all = reads_table_all[,keep]
  reads_table_all = reads_table_all[row.names(metadata_all)]
  
  {
    parameter_list = read.csv('parameter_list.csv', row.names = 1, header = T)
    
    parameter_20 = as.character(parameter_list$value[row.names(parameter_list) == 'parameter_20']);
    parameter_28 = as.character(parameter_list$value[row.names(parameter_list) == 'parameter_28']);
    if (parameter_28 == "Yes") {
      parameter_28 = T
    } else {parameter_28 = F}
    
    parameter_27 = as.character(parameter_list$value[row.names(parameter_list) == 'parameter_27']);
    if (!is.na(parameter_27)[1]) {
      parameter_27 = strsplit(parameter_27,',|;| ')
      parameter_27 = parameter_27[[1]]
    }
    
    parameter_list$value = as.numeric(as.character(parameter_list$value))
    
    parameter_1 = parameter_list$value[row.names(parameter_list) == 'parameter_1'];
    parameter_2 = parameter_list$value[row.names(parameter_list) == 'parameter_2'];
    parameter_3 = parameter_list$value[row.names(parameter_list) == 'parameter_3'];
    parameter_4 = parameter_list$value[row.names(parameter_list) == 'parameter_4'];
    parameter_5 = parameter_list$value[row.names(parameter_list) == 'parameter_5'];
    parameter_6 = parameter_list$value[row.names(parameter_list) == 'parameter_6'];
    parameter_7 = parameter_list$value[row.names(parameter_list) == 'parameter_7'];
    parameter_8 = parameter_list$value[row.names(parameter_list) == 'parameter_8'];
    parameter_9 = parameter_list$value[row.names(parameter_list) == 'parameter_9'];
    parameter_10 = parameter_list$value[row.names(parameter_list) == 'parameter_10'];
    parameter_11 = parameter_list$value[row.names(parameter_list) == 'parameter_11'];
    parameter_12 = parameter_list$value[row.names(parameter_list) == 'parameter_12'];
    parameter_13 = parameter_list$value[row.names(parameter_list) == 'parameter_13'];
    parameter_14 = parameter_list$value[row.names(parameter_list) == 'parameter_14'];
    parameter_15 = parameter_list$value[row.names(parameter_list) == 'parameter_15'];
    parameter_16 = parameter_list$value[row.names(parameter_list) == 'parameter_16'];
    parameter_17 = parameter_list$value[row.names(parameter_list) == 'parameter_17'];
    parameter_18 = parameter_list$value[row.names(parameter_list) == 'parameter_18'];
    parameter_19 = parameter_list$value[row.names(parameter_list) == 'parameter_19'];
    parameter_21 = parameter_list$value[row.names(parameter_list) == 'parameter_21'];
    parameter_22 = parameter_list$value[row.names(parameter_list) == 'parameter_22'];
    parameter_23 = parameter_list$value[row.names(parameter_list) == 'parameter_23'];
    parameter_24 = parameter_list$value[row.names(parameter_list) == 'parameter_24'];
    parameter_25 = parameter_list$value[row.names(parameter_list) == 'parameter_25'];
    parameter_26 = parameter_list$value[row.names(parameter_list) == 'parameter_26'];
    
    
    
  }
}

##### analysize results #####
setwd(path_input)
path_output = path_input
file_list = list.files(path = ".")
path_output = paste0(path_input,'/results')
path_output_1 = path_output
n=1
name_1 = 'results'
while (sum(!(name_1 %in% file_list)) == 0) {
  path_output_1 = paste0(path_output,n)
  name_1 = c(name_1,paste0('results',n))
  n = n+1
}
dir.create(path_output_1)
setwd(path_output_1)


#ref_group <- readline(prompt= paste0("Enter a control group (", Factor_levels, '):   '))

Factor_levels = as.character(unique(metadata_all[,1]))
Factor_levels = Factor_levels[order(Factor_levels)]

if (is.na(parameter_27)[1] | sum(Factor_levels %in% parameter_27) != length(Factor_levels)) {
  order = Factor_levels
} else {
  order = parameter_27
}

ref_group = order[1]

reads_table = reads_table_all
metadata = metadata_all
One_factor_script(reads_table, metadata, ref_group, order, VMB = parameter_28)










