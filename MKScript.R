##Script for running McDonald-Kreitman Test on Melospiza data

#poorly written by David JX Tan
#based on math stolen from:
# https://www.picb.ac.cn/evolgen/softwares/download/NeutralityTest/MKtest.htm
# https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test

library(seqinr)
library(dplyr)
library(gtools)
library(stringr)
codontable <- read.csv("C:/Users/david/Dropbox/melospiza-genes/Genes_Renamed/codontable.csv",header=T)

#this is the workhorse function that calculates the number of synonymous and nonsynonymous sites for each pairwise codon comparison
codon_syn_nonsyn <- function(codon1,codon2) {
  s=0
  n=0
  counter=0
  nmut <- as.numeric(adist(codon1,codon2))
  #if there is just one substitution, the math is pretty simple
  if (nmut == 1) {
    #check to see if the substitution is synonymous or non-synonymous
    if (filter(codontable,Codon==codon1)$AA==filter(codontable,Codon==codon2)$AA) {
      s = 1
    } else {
      n = 1
    }
  } else if (nmut == 2) { #when there are 2 substitutions in a locus, there are 2 possible pathways to consider
    s2=0
    n2=0
    transformation <- unlist(strsplit(attr(adist(codon1,codon2,counts=T,costs=list(ins=10000,del=10000,sub=0)),"trafos"),""))
    subst_pos <- str_which(transformation,"S")
    perm_matrix <- gtools::permutations(n=length(subst_pos),r=length(subst_pos),subst_pos)
    for (j in 1:2) { #iterate over each pathway
      codoni1 <- codon1
      substr(codoni1,perm_matrix[j,1],perm_matrix[j,1]) <- substr(codon2,perm_matrix[j,1],perm_matrix[j,1]) #execute the first substitution and check for synonymy
      if (codoni1 == "taa" | codoni1 == "tag" | codoni1 == "tga") next
      if (filter(codontable,Codon==codon1)$AA==filter(codontable,Codon==codoni1)$AA) {
        s2 = s2+1
      } else {
        n2 = n2+1
      }
      if (filter(codontable,Codon==codoni1)$AA==filter(codontable,Codon==codon2)$AA) {
        s2 = s2+1
      } else {
        n2 = n2+1
      }
      counter = counter +1
    }
    s = s2/counter #since there are 2 pathways, we have do divide the number of s and n sites by 2
    n = n2/counter
  } else if (nmut == 3) { #This one is a huge PITA
    s3 = 0
    n3 = 0
    subst_pos <- c(1,2,3)
    perm_matrix <- gtools::permutations(n=length(subst_pos),r=length(subst_pos),subst_pos)
    for (k in 1:6) {
      codoni1 <- codon1
      codoni2 <- codon2
      substr(codoni1,perm_matrix[k,1],perm_matrix[k,1]) <- substr(codon2,perm_matrix[k,1],perm_matrix[k,1])
      if (codoni1 == "taa" | codoni1 == "tag" | codoni1 == "tga") next
      if (filter(codontable,Codon==codon1)$AA==filter(codontable,Codon==codoni1)$AA) {
        s3 = s3+1
      } else {
        n3 = n3+1
      }
      codoni2 <- codoni1
      substr(codoni2,perm_matrix[k,2],perm_matrix[k,2]) <- substr(codon2,perm_matrix[k,2],perm_matrix[k,2])
      if (codoni2 == "taa" | codoni2 == "tag" | codoni2 == "tga") next
      if (filter(codontable,Codon==codoni2)$AA==filter(codontable,Codon==codoni1)$AA) {
        s3 = s3+1
      } else {
        n3 = n3+1
      }
      if (filter(codontable,Codon==codoni2)$AA==filter(codontable,Codon==codon2)$AA) {
        s3 = s3+1
      } else {
        n3 = n3+1
      }
      counter = counter +1
    }
    s = s3/counter
    n = n3/counter
  }
  return(c(s,n))
}

#this function runs the MK test on each gene 
calc_melospiza_MK = function(fastafile) {
  #read fasta file for specific gene
  gene <- seqinr::read.fasta(fastafile)
  
  #calculate length of the alignment
  numcodons <- length(gene$maxima)/3
  
  #let's torture some sequence data into codons
  maxima_UAM_codons <- as.vector(sapply(split(gene$maxima_UAM31500,ceiling(seq_along(gene$maxima_UAM31500)/3)),paste0,collapse=""))
  maxima_codons <- as.vector(sapply(split(gene$maxima,ceiling(seq_along(gene$maxima)/3)),paste0,collapse=""))
  georgiana_codons <- as.vector(sapply(split(gene$georgiana,ceiling(seq_along(gene$georgiana)/3)),paste0,collapse=""))
  #convert codon positions into a dataframe
  gene_codons <- as.data.frame(cbind(maxima_UAM_codons,maxima_codons,georgiana_codons))
  #get rid of sites with Ns or gaps
  gene_codons <- filter(gene_codons, !grepl('N|-',maxima_codons) & !grepl('N|-',maxima_UAM_codons) & !grepl('N|-',georgiana_codons))
  #remove STOP codons
  gene_codons <- filter(gene_codons, !grepl("taa|tag|tga",maxima_codons) & !grepl("taa|tag|tga",maxima_UAM_codons) & !grepl("taa|tag|tga",georgiana_codons))
  #remove incomplete codons
  gene_codons <- filter(gene_codons, !nchar(maxima_codons) < 3 | !nchar(maxima_UAM_codons) < 3 | !nchar(georgiana_codons) < 3)
  
  
  #filter for polymorphic sites within maxima
  p_sites <- filter(gene_codons,maxima_UAM_codons!=maxima_codons)
  #filter for divergent sites between maxima and georgiana
  d_sites <- d_sites <- filter(gene_codons,maxima_UAM_codons!=georgiana_codons & maxima_codons!=georgiana_codons)
  
  if (nrow(d_sites)==0 | nrow(p_sites)==0) {
    return(c(NA,NA,NA,NA,NA,NA,NA))
  } else {
    #set up variables 
    Ds=0 #divergent synonymous substitutions
    Dn=0 #divergent nonsynonymous substitutions
    Ps=0 #polymorphic synonymous substitutions
    Pn=0 #polymorphic nonsynonymous substitutions
    
    #here we calculate the Ps and Pn values
    for (i in 1:nrow(p_sites)) {
      ns <- codon_syn_nonsyn(p_sites[i,1],p_sites[i,2])
      nmut <- adist(p_sites[i,1],p_sites[i,2])
      syn = as.numeric(ns[1] * nmut / (ns[1] + ns[2]))
      nsyn = as.numeric(ns[2] * nmut / (ns[1] + ns[2]))
      Ps = Ps + syn
      Pn = Pn + nsyn
    }
    
    for (x in 1:nrow(d_sites)) {
      #count number of fixed sites
      comp13 <- unlist(strsplit(attr(adist(d_sites[x,1],d_sites[x,3],counts = T,costs = list(ins=10000,del=10000,sub=0)),"trafos"),""))
      comp23 <- unlist(strsplit(attr(adist(d_sites[x,2],d_sites[x,3],counts = T,costs = list(ins=10000,del=10000,sub=0)),"trafos"),""))
      site1 <- floor(sum(stringr::str_count(c(comp13[1],comp23[1]),"S"))/2)
      site2 <- floor(sum(stringr::str_count(c(comp13[2],comp23[2]),"S"))/2)
      site3 <- floor(sum(stringr::str_count(c(comp13[3],comp23[3]),"S"))/2)
      nfixed <- sum(site1,site2,site3)
      
      #if both maxima are identical, only 1 pairwise comparison is needed
      if (d_sites[x,1]==d_sites[x,2]) {
        ns <- codon_syn_nonsyn(d_sites[x,1],d_sites[x,3])
        syn = as.numeric(ns[1] * nfixed / (ns[1] + ns[2]))
        nsyn = as.numeric(ns[2] * nfixed / (ns[1] + ns[2]))
        Ds = Ds + syn
        Dn = Dn + nsyn
      } else if (d_sites[x,1]!=d_sites[x,2]) { #if both maxima codons are different, then we need to calculate every pairwise comparison between maxima and georgiana
        ns_1 <- codon_syn_nonsyn(d_sites[x,1],d_sites[x,3])
        ns_2 <- codon_syn_nonsyn(d_sites[x,2],d_sites[x,3])
        total_sites = ns_1[1] + ns_1[2] + ns_2[1] + ns_2[2]
        syn = as.numeric((ns_1[1] + ns_2[1]) * nfixed / total_sites)
        nsyn = as.numeric((ns_1[2] + ns_2[2]) * nfixed / total_sites)
        Ds = Ds + syn
        Dn = Dn + nsyn
      }
    }
    #calculate metrics 
    neutrality = (Pn/Ps)/(Dn/Ds)
    alpha = 1 - (Ds*Pn)/(Dn*Ps)
    fisher = fisher.test(matrix(c(Dn,Pn,Ds,Ps),ncol=2))
    output <- data.frame(Ds,Dn,Ps,Pn,neutrality,alpha,fisher$p.value)
    colnames(output) <- c("Ds","Dn","Ps","Pn","neutrality","alpha","fisher.pval")
    return(output)
  }
}


genelist <- c("AFAP1","AGRP","ALX1","APOD","ASIP","BCO1","BCO2","CCDC149","CD36","CLOCK","CREB1","CRY1","CRY2","DLK1","GGA1","HMGA2","HPS5","KCTD21","LAP3","LCORL","LGI2","MC1R","MITF","MMP17","MYOF","NPAS3","PCSK2","PER2","PER3","QDPR","SLC2A1","SLIT2","STARD3","TGFB2","TYR","WAPL","WNK2","YPEL1")
mktest_out <- data.frame()
for (i in genelist) {
  mk_results <- calc_melospiza_MK(paste0("C:/Users/david/Dropbox/melospiza-genes/cds_aligned/",i,"_aligned.fasta"))
  mktest_out <- rbind(mktest_out,c(i,as.numeric(mk_results)))
}
colnames(mktest_out) <- c("Gene","Ds","Dn","Ps","Pn","neutrality","alpha","fisher.pval")
write.csv(mktest_out,"C:/Users/david/Dropbox/melospiza-genes/cds_aligned/mktest_results.csv")
