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
    trials = c(1: ncol(reads_table))
    func_1 = function(trial) {
      reads_table_abundance <- reads_table[,trial] / colSums(reads_table)[trial]
      return(reads_table_abundance)
    }
    reads_table_abundance = mclapply(trials, func_1, mc.cores = mc.cores)
    reads_table_abundance = unlist(reads_table_abundance)
    
    reads_table_abundance_2 = as.data.frame(matrix(data =0, ncol = ncol(reads_table), nrow = nrow(reads_table)))
    row.names(reads_table_abundance_2) <- row.names(reads_table)
    colnames(reads_table_abundance_2) <- colnames(reads_table)
    
    for (a in 1: ncol(reads_table)) {
      reads_table_abundance_2[,a] = reads_table_abundance[((a-1)*nrow(reads_table)+1):(a*nrow(reads_table))]
    }
    reads_table_abundance = reads_table_abundance_2
    
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
      print('no factor name')
      return(NA)
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
    
    
    alpha.shannon = ggplot(metadata, aes(x=factor, y=alpha.shannon)) +
      geom_boxplot(aes(fill = factor),outlier.shape=NA) +
      labs(x = NULL, y = "Shannon index", fill=factor_name)+ 
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            axis.text.x=element_blank(),
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))
    
    alpha.evenness = ggplot(metadata, aes(x=factor, y=alpha.evenness)) +
      geom_boxplot(aes(fill = factor),outlier.shape=NA) +
      labs(x = NULL, y = "Evenness", fill=factor_name)+ 
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            axis.text.x=element_blank(),
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))
    
    alpha.ovserved_OTU = ggplot(metadata, aes(x=factor, y=alpha.ovserved_OTU)) +
      geom_boxplot(aes(fill = factor),outlier.shape=NA) +
      labs(x = NULL, y = "Observed OTU", fill=factor_name)+ 
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            axis.text.x=element_blank(),
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))
    # geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5)+
    
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
  beta_diversity = function(reads_table, metadata = NA, factor_name = NA, order = NA, NMDS_skip = T, ref_group = NA, rarefy_to = NA, pheatmap_fontsize = 5,treeheight = 50, pheatmap_y = T) {
    #reads_table; metadata = as.character(metadata[,1]); factor_name = factor_1; 
    #order = order; NMDS_skip = F; ref_group = ref_group; 
    #rarefy_to = NA; pheatmap_fontsize = 5; treeheight = 5; 
    #pheatmap_y = T
    
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
          pvalue <- adonis2(reads_table_2 ~ ., data = metadata_2, method = "bray")[1,5] 
          group_dis_sig[a,b] = pvalue
          group_dis_sig[b,a] = pvalue
        }
      }
      
    } 
    
    # Running Nonmetric Multidimensional Scaling (NMDS) Ordination
    if (NMDS_skip == T) {
      # output
      output = c(Bray_Curtis = list(Bray_Curtis), group_dis_2_p = list(group_dis_2_p),
                 group_dis_sig = list(group_dis_sig))
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
      output = c(Bray_Curtis = list(Bray_Curtis), group_dis_2_p = list(group_dis_2_p),
                 group_dis_sig = list(group_dis_sig), NMDS =list(NMDS), NMDS_2 =list(NMDS_2))
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
  
  
  
  
  
  ######### bar plot ########## samples on columns ###########
  barplot <- function(reads_table, type_th = 0.1, taxa_num = 9, VMB = T) { # pass proportion data, samples on rows, taxa on columns.

    if (VMB == F) {
      reads_table_abundance <- matrix(data =0, ncol = ncol(reads_table),nrow = nrow(reads_table))
      
      for (a in 1:ncol(reads_table)) {
        reads_table_abundance[,a] <- reads_table[,a] / colSums(reads_table)[a]
        
      }
      row.names(reads_table_abundance) = row.names(reads_table)
      colnames(reads_table_abundance) = colnames(reads_table) 
      
      reads_table = as.data.frame(t(reads_table_abundance))
      
      reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
      
      # assign vagitype
      mytypes <- apply(reads_table,1,which.max)   # find the most abundant taxa
      maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes), ncol=2)]  # the abundance of the most abundant taxa
      mytypes <- colnames(reads_table)[mytypes]   # find the name of the most abundant taxa
      mytypes[maxprop < type_th] <- "No Type"
      uniqueTypes <- unique(mytypes)
      
      top_species = as.data.frame(sort(colSums(reads_table),decreasing=T))
      
      keep = mytypes %in% row.names(top_species)[1:taxa_num]
      
      mytypes[!keep] = 'Others'
      
      if (length(uniqueTypes) < (taxa_num+1)) {
        myTypeOrder <- c( row.names(top_species)[1:(length(uniqueTypes)-1)], "No Type")
      } else {
        myTypeOrder <- c( row.names(top_species)[1:taxa_num], "No Type","Others")
      }
      myTypeOrder <- data.frame(typeOrder=c(1:length(myTypeOrder)),row.names=myTypeOrder)
      
      reads_table <- 100*reads_table
      
      # order data by type then proportion
      reads_table <- reads_table[order(myTypeOrder[mytypes,],-maxprop),]  # reorder the samples by vagitype and the abundance of dominant taxa
      barplottypes <- mytypes[order(myTypeOrder[mytypes,],-maxprop)]   # get the list of vagitypes of samples
      
      # color
      mycolors <- rep("gray", ncol(reads_table))
      
      color_code_list = c("#3264B8","#FC9800","#C80B0B","#FDFFBA","#55ff42",
                          "#F7A0A0","#D7A7F1","#337F08","#00FFE8","#9C58FE","#A68300",'gray')
      
      color_code_list = color_code_list[1:(taxa_num+1)]
      
      for (a in 1:(length(color_code_list)-1)) {
        mycolors[which(colnames(reads_table) == row.names(top_species)[a])] <- color_code_list[a]
      }
      
      mycolors <- factor(mycolors)
      mycolordf <- data.frame(mycolors, Taxa=colnames(reads_table))
      
      mycolordf2 = mycolordf[mycolordf$mycolors != 'gray',]
      mycolordf2 = mycolordf2[1:taxa_num,]
      mycolordf2[(taxa_num+1),1] = "gray"
      mycolordf2[(taxa_num+1),2] = "Others"
      mycolors_label = mycolordf2$Taxa[order(mycolordf2$mycolors)]
      
      sampleorder <- factor(rownames(reads_table))
      
      ggplotdata <- melt(as.matrix(reads_table),  
                         id.vars=names(reads_table), 
                         varnames=c("SampleID", "Taxa"), 
                         value.name="ATprop")
      ggplotdata$SampleID <- factor(ggplotdata$SampleID, levels=sampleorder)
      ggplotdata <- merge(ggplotdata, mycolordf, by="Taxa")
      ggplotdata <- ggplotdata[ggplotdata$ATprop != 0,]
      
      p <- ggplot(ggplotdata, aes(SampleID, ATprop, fill=mycolors, group=ATprop)) + 
        geom_bar(stat="identity", position="stack", width=1) +
        scale_fill_manual(values = levels(mycolors),labels = mycolors_label) +
        labs(x="Sample", y="Relative Abundance") + 
        theme(axis.text.x=element_blank(),
              axis.ticks=element_blank()) 
      
      
      #,labels = mycolors_label
      output = c(p = list(p), color = list(mycolordf))
      return(output)
    } else {
      # convert reads to relative abundance
      {
        reads_table_abundance <- matrix(data =0, ncol = ncol(reads_table),nrow = nrow(reads_table))
        for (a in 1:ncol(reads_table)) {
          reads_table_abundance[,a] <- reads_table[,a] / colSums(reads_table)[a]
          
        }
        row.names(reads_table_abundance) = row.names(reads_table)
        colnames(reads_table_abundance) = colnames(reads_table) 
        mypropdata = as.data.frame(t(reads_table_abundance))
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
      mypropdata <- mypropdata[order(myTypeOrder[mytypes,],-maxprop),]  # myTypeOrder[mytypes,]: find the order number of mytypes in the myTypeOrder list    # order(a,b): order a then order b
      barplottypes <- mytypes[order(myTypeOrder[mytypes,],-maxprop)]
      
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
        mycolors[grep("Streptococcus",        colnames(mypropdata))]     <- "#9edae5"
        mycolors[grep("Dialister",        colnames(mypropdata))]     <- "#dbdb8d"
        mycolors[grep("Megasphaera",        colnames(mypropdata))]     <- "#f7b6d2"
        mycolors[grep("Coriobacteriaceae",        colnames(mypropdata))]     <- "#bdc2cc"
        mycolors[grep("Ureaplasma",        colnames(mypropdata))]     <- "#2ca02c"
        mycolors[grep("Mobiluncus",        colnames(mypropdata))]     <- "#c4c392"
        mycolors[grep("TM7",        colnames(mypropdata))]     <- "#e377c2"
        mycolors[grep("Fusobacterium",        colnames(mypropdata))]     <- "#c8d045"
        
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
        
        keep = colnames(mypropdata) %in% row.names(myTypeOrder)
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
          theme(axis.text.x=element_blank(), 
                axis.ticks=element_blank())  
        
      } 
      
      # legend figure
      {
        mycolordf_2 = mycolordf[mycolordf$Taxa %in% row.names(myTypeOrder),]
        mycolordf_2 = mycolordf_2[!duplicated(mycolordf_2$mycolors),]
        mycolordf_2$Taxa[mycolordf_2$mycolors == 'gray'] = 'Others'
        mycolordf_2$value =1
        
        mycolordf_2$order = myTypeOrder[mycolordf_2$Taxa,]
        mycolordf_2 = mycolordf_2[order(mycolordf_2$order),]
        
        mycolordf_2$Taxa = factor(mycolordf_2$Taxa, levels = mycolordf_2$Taxa)
        mycolordf_2$mycolors = as.character(mycolordf_2$mycolors)
        
        p2 = ggplot(mycolordf_2, aes(x = Taxa, y = value, fill = Taxa)) +
          geom_bar(stat = "identity")+
          scale_fill_manual(values= mycolordf_2$mycolors)
      }
      
      output = c(p1 = list(p1), p2 = list(p2))
      return(output)
    }
    
  }
  
  
  
  
  ######### heat map ########
  heatmap_plot <- function(reads_table, metadata, type_th = 0.1, vagitype_num = 5, abundant_taxa_number = 20, width=18, height=8) { # pass proportion data, samples on columns, taxa on rows.
    # reads_table = reads_table_all
    # metadata = metadata_all
    reads_table_abundance <- matrix(data =0, ncol = ncol(reads_table),nrow = nrow(reads_table))
    
    for (a in 1:ncol(reads_table)) {
      reads_table_abundance[,a] <- reads_table[,a] / colSums(reads_table)[a]
      
    }
    row.names(reads_table_abundance) = row.names(reads_table)
    colnames(reads_table_abundance) = colnames(reads_table) 
    
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
    save_pheatmap_pdf(heatmap_plot, "Heatmap_clustering.pdf", width=width, height=height)
    
    
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
      
      reads_table_abundance = get_abundance_table(reads_table, mc.cores = 8)
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
    data = prepare_reads_table(reads_table, metadata, total_reads_threshold = parameter_21, species_threshold = parameter_22, mc.cores = parameter_24) 
    reads_table = data$reads_table
    metadata = data$metadata
    
    # alpha diversity
    alpha = alpha_diversity(reads_table, metadata = as.character(metadata[,1]), factor_name = factor_1, paired = F,order = order)
    alpha$shannon + theme_bw() 
    ggsave('alpha_shannon.pdf',width=parameter_1, height=parameter_2)
    alpha$evenness + theme_bw() 
    ggsave('alpha_evenness.pdf',width=parameter_1, height=parameter_2)
    alpha$ovserved_OTU + theme_bw() 
    ggsave('alpha_ovserved_OTU.pdf',width=parameter_1, height=parameter_2)
    write.csv(alpha$sig_Shannon,'alpha_shannon.csv', row.names = T, quote = F)
    write.csv(alpha$sig_Evenness,'alpha_evenness.csv', row.names = T, quote = F)
    write.csv(alpha$sig_OTU,'alpha_ovserved_OTU.csv', row.names = T, quote = F)
    
    # beta diversity
    beta = beta_diversity(reads_table, metadata = as.character(metadata[,1]), factor_name = factor_1, 
                          order = order, NMDS_skip = F, ref_group = ref_group, 
                          rarefy_to = NA, pheatmap_fontsize = 5, treeheight = 5, 
                          pheatmap_y = T)
    
    beta$NMDS_2
    ggsave('NMDS_2.pdf',width=parameter_3, height=parameter_4)
    
    beta$group_dis_sig
    write.csv(beta$group_dis_sig, 'Adonis_test.csv')

    metadata_level = length(unique(metadata[,1]))
    metadata_level_2 = unique(metadata[,1])
    
    # composition bar plot
    if (VMB == T) {
      for (a in 1: metadata_level) {
        keep = metadata[,1] == metadata_level_2[a]
        reads_table_2 = reads_table[,keep]
        
        barplot_result = barplot(reads_table_2, type_th = 0.3, taxa_num = 9) 

        pdf(paste0('bar_plot_', metadata_level_2[a], '_1.pdf') ,width=parameter_5, height=parameter_6) 
        print(barplot_result$p1)
        dev.off()
        
        pdf(paste0('bar_plot_', metadata_level_2[a], '_2.pdf') ,width=parameter_5, height=parameter_6) 
        print(barplot_result$p2)
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
          ggsave('abundance_difference.pdf',width=parameter_7, height=parameter_8)
          
        }
      }
      
    }
    
    # network
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
      parameter_27 = strsplit(parameter_27,',')
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

if (ncol(metadata_all) == 1) {
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
  
  
  Factor_levels = as.character(unique(metadata_all[,1]))
  Factor_levels = paste(Factor_levels, collapse = '/')
  ref_group <- readline(prompt= paste0("Enter a control group (", Factor_levels, '):   '))
  
  Factor_levels = as.character(unique(metadata_all[,1]))
  Factor_levels = Factor_levels[order(Factor_levels)]
  
  if (is.na(parameter_27)[1] | sum(Factor_levels %in% parameter_27) != length(Factor_levels)) {
    order = c(ref_group, Factor_levels[-which(Factor_levels == ref_group)])
  } else {
    order = parameter_27
  }
  
  # heatmap
  heatmap_plot(reads_table_all, metadata_all,type_th = parameter_14, vagitype_num = parameter_17, 
               width=parameter_15, height=parameter_16, abundant_taxa_number = parameter_13) 
  
  One_factor_script(reads_table_all, metadata_all, ref_group, order)
} else {
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
  
  factor_number = ncol(metadata_all)
  
  # heatmap
  heatmap_plot(reads_table_all, metadata_all,type_th = parameter_14, vagitype_num = parameter_17,
               width=parameter_15, height=parameter_16, abundant_taxa_number = parameter_13) 
  
  for (factor_analysis in seq(factor_number)) {
    metadata_all_2 = metadata_all
    metadata_all_2[,1] = metadata_all_2[,factor_analysis]
    colnames(metadata_all_2)[1] = colnames(metadata_all_2)[factor_analysis]
    metadata_all_2[,c(2:ncol(metadata_all_2))] = NULL

    file_list = list.files(path = ".")
    path_output_2 = paste0(path_output_1,'/',colnames(metadata_all_2)[1])
    if (sum(colnames(metadata_all_2)[1] %in% file_list) == 0) {
      dir.create(path_output_2)
      setwd(path_output_2)
    }
    
    Factor_levels = as.character(unique(metadata_all_2[,1]))
    Factor_levels = paste(Factor_levels, collapse = '/')
    ref_group <- readline(prompt= paste0("Enter a control group (", Factor_levels, '):   '))
    
    Factor_levels = as.character(unique(metadata_all_2[,1]))
    Factor_levels = Factor_levels[order(Factor_levels)]
    
    if (is.na(parameter_27)[1] | sum(Factor_levels %in% parameter_27) != length(Factor_levels)) {
      order = c(ref_group, Factor_levels[-which(Factor_levels == ref_group)])
    } else {
      order = parameter_27
    }
    
    One_factor_script(reads_table_all, metadata_all_2, ref_group, order)
  }
}








