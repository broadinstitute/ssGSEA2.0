
# $Id: ssGSEA.R 66 2015-09-24 20:46:31Z manidr $
print ("$Id: ssGSEA.R 66 2015-09-24 20:46:31Z manidr $")


##
## Single sample GSEA
## Original code written by Pablo Tamayo. Adapted with additional modifications
## Referece:
#    1. Abazeed, M. E., Adams, D. J., Hurov, K. E., Tamayo, P., Creighton, C. J., Sonkin, D., et al. (2013). 
#       Integrative Radiogenomic Profiling of Squamous Cell Lung Cancer. Cancer Research, 73(20), 6289–6298. 
#       http://doi.org/10.1158/0008-5472.CAN-13-1616
#    2. Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005). 
#       Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. 
#       Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545–15550.
#


library(MASS)
library(verification)



Rutil.path <- switch (Sys.info()[['sysname']],
                      Windows = {'//flynn-cifs/prot_proteomics/LabMembers/manidr/R-utilities'},
                      Darwin = {'/Volumes/prot_proteomics/LabMembers/manidr/R-utilities'},
                      Linux = {'/prot/proteomics/LabMembers/manidr/R-utilities'})
# default genesets (currently MSigDB C2)
gsea.dbs <- c (file.path (Rutil.path, 'GSEA', 'c2.all.v4.0.symbols.gmt'))



ssGSEA <- function ( 
  input.ds,                      # input data file in gct format, first column (Name) must contain gene symbols
  output.prefix,                 # prefix used for output tables
  gene.set.databases=gsea.dbs,   # list of genesets (in gmt format) to evaluate enrichment on
  gene.set.selection  = "ALL",   # "ALL" or list with names of gene sets
  sample.norm.type    = "rank",  # sample normalization: "rank", "log", "log.rank" or "none"
  weight              = 0,       # when weight==0, all genes have the same weight; if weight>0, 
                                 # actual values matter, and can change the resulting score
  statistic           = "area.under.RES",
                                 # alternatives: "Kolmogorov-Smirnov", "area.under.RES"
  output.score.type   = "NES",   # "ES" or "NES"
  nperm               = 100,     # number of random permutations for NES case
  combine.mode        = "combine.off",  
                                 # "combine.off" do not combine *_UP and *_DN versions in a single score.
                                 # "combine.replace" combine *_UP and *_DN versions in a single score.
                                 # "combine.add" combine *_UP and *_DN versions in a single score and 
                                 #    add it but keeping the individual *_UP and *_DN versions.
  min.overlap         = 10,
  correl.type         = "rank",  # correlation type: "rank", "z.score", "symm.rank"
  fdr.pvalue          = TRUE,    # output adjusted (FDR) p-values
  global.fdr          = FALSE    # if TRUE calculate global FDR; else calculate FDR sample-by-sample
) { 
  ## single sample GSEA
  # for results similar to the Java version, use: weight=0;
  # when weight==0, sample.norm.type and correl.type do not matter;
  # when weight > 0, the combination of sample.norm.type and correl.type
  #  dictate how the gene expression values in input.ds are transformed
  #  to obtain the score -- use this setting with care (the transformations
  #  can skew scores towards +ve or -ve values)
  #  . sample.norm.type=='none' uses actual expression values; combined
  #    with correl.type=='rank', genes are weighted by actual values
  #  . sample.norm.type=='rank' weights genes proportional to rank
  #  . sample.norm.type=='log' can be used for log-transforming input data
  #  . correl.type=='z.score' standardizes the (normalized) input values 
  #    before using them to calculate scores.
  
  # Load libraries
  library(gtools)
  library(verification)
  library(RColorBrewer)
  
  # Read input dataset  
  dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read dataset (GCT format)
  m <- data.matrix(dataset$ds)
  gene.names <- dataset$row.names
  gene.descs <- dataset$descs
  sample.names <- dataset$names
  Ns <- length(m[1,])
  Ng <- length(m[,1])
  temp <- strsplit(input.ds, split="/")  # Extract input file name
  s <- length(temp[[1]])
  input.file.name <- temp[[1]][s]
  temp <- strsplit(input.file.name, split=".gct")
  input.file.prefix <-  temp[[1]][1]
  
  # Sample normalization
  if (sample.norm.type == "rank") {
    # column rank normalization 
    for (j in 1:Ns) {  
      m[,j] <- rank(m[,j], ties.method = "average")
    }
    m <- 10000*m/Ng
  } else if (sample.norm.type == "log.rank") {
    # column rank normalization 
    for (j in 1:Ns) {  
      m[,j] <- rank(m[,j], ties.method = "average")
    }
    m <- log(10000*m/Ng + exp(1))
  } else if (sample.norm.type == "log") {
    m[m < 1] <- 1
    m <- log(m + exp(1))
  } else if (sample.norm.type == "none") {
    # keep original value -- do not transform
  }
  
  # Read gene set databases
  max.G <- 0
  max.N <- 0
  for (gsdb in gene.set.databases) {
    GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
    max.G <- max(max.G, max(GSDB$size.G))
    max.N <- max.N +  GSDB$N.gs
  }
  N.gs <- 0
  gs <- matrix("null", nrow=max.N, ncol=max.G)
  gs.names <- vector(length=max.N, mode="character")
  gs.descs <- vector(length=max.N, mode="character")
  size.G <- vector(length=max.N, mode="numeric")
  start <- 1
  for (gsdb in gene.set.databases) {
    GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
    N.gs <- GSDB$N.gs 
    gs.names[start:(start + N.gs - 1)] <- GSDB$gs.names
    gs.descs[start:(start + N.gs - 1)] <- GSDB$gs.desc
    size.G[start:(start + N.gs - 1)] <- GSDB$size.G
    gs[start:(start + N.gs - 1), 1:max(GSDB$size.G)] <- GSDB$gs[1:N.gs, 1:max(GSDB$size.G)]
    start <- start + N.gs
  }
  N.gs <- max.N
  
  # Select desired gene sets
  if (gene.set.selection[1] == "ALL") {
    gene.set.selection <- unique(gs.names)
  } 
  
  locs <- match(gene.set.selection, gs.names)
  # print (rbind(gene.set.selection, locs))
  N.gs <- sum(!is.na(locs))
  if(N.gs > 1) { 
    gs <- gs[locs,]
  } else { 
    gs <- t(as.matrix(gs[locs,]))   # Force vector to matrix if only one gene set specified
  }
  gs.names <- gs.names[locs]
  gs.descs <- gs.descs[locs]
  size.G <- size.G[locs]
  
  # Check for redundant gene sets 
  tab <- as.data.frame(table(gs.names))
  ind <- order(tab[, "Freq"], decreasing=T)
  tab <- tab[ind,]
  print(tab[1:10,])
  print(paste("Total gene sets:", length(gs.names)))
  print(paste("Unique gene sets:", length(unique(gs.names))))
  
  # Loop over gene sets
  score.matrix <- matrix(0, nrow=N.gs, ncol=Ns)
  pval.matrix <- matrix(0, nrow=N.gs, ncol=Ns)
  for (gs.i in 1:N.gs) {
    gene.set <- gs[gs.i, 1:size.G[gs.i]]
    gene.overlap <- intersect(gene.set, gene.names)
    print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
    if (length(gene.overlap) < min.overlap) { 
      score.matrix[gs.i, ] <- rep(NA, Ns)
      pval.matrix[gs.i, ] <- rep(NA, Ns)
      next
    } else {
      gene.set.locs <- match(gene.overlap, gene.set)
      gene.names.locs <- match(gene.overlap, gene.names)
      msig <- m[gene.names.locs,]
      msig.names <- gene.names[gene.names.locs]
      if (output.score.type == "ES") {
        OPAM <- project.geneset (data.array = m, gene.names = gene.names, n.cols = Ns, 
                                   n.rows = Ng, weight = weight, statistic = statistic,
                                   gene.set = gene.overlap, nperm = 1, correl.type = correl.type)
        score.matrix[gs.i,] <- OPAM$ES.vector
      } else if (output.score.type == "NES") {
        OPAM <- project.geneset (data.array = m, gene.names = gene.names, n.cols = Ns, 
                                   n.rows = Ng, weight = weight, statistic = statistic,
                                   gene.set = gene.overlap, nperm = nperm, correl.type = correl.type)
        score.matrix[gs.i,] <- OPAM$NES.vector
      }
      pval.matrix[gs.i,] <- OPAM$p.val.vector
    }
  }
  
  
  locs <- !is.na(score.matrix[,1])
  print(paste("N.gs before overlap prunning:", N.gs))
  N.gs <- sum(locs)
  print(paste("N.gs after overlap prunning:", N.gs))
  score.matrix <- score.matrix[locs,]
  pval.matrix <- pval.matrix[locs,]
  gs.names <- gs.names[locs]
  gs.descs <- gs.descs[locs]
  
  initial.up.entries <- 0
  final.up.entries <- 0
  initial.dn.entries <- 0
  final.dn.entries <- 0
  combined.entries <- 0
  other.entries <- 0
  
  if (combine.mode == "combine.off") {
    score.matrix.2 <- score.matrix
    pval.matrix.2 <- pval.matrix
    gs.names.2 <- gs.names
    gs.descs.2 <- gs.descs
  } else if ((combine.mode == "combine.replace") || (combine.mode == "combine.add")) {
    fisher.pval <- function(p) {
      Xsq <- -2*sum(log(p))
      p.val <- pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
      return(p.val)
    }
    
    score.matrix.2 <- NULL
    pval.matrix.2 <- NULL
    gs.names.2 <- NULL
    gs.descs.2 <- NULL
    
    add.entry.2 <- function (s, p, n, d) {
      score.matrix.2 <<- rbind (score.matrix.2, s)
      pval.matrix.2 <<- rbind (pval.matrix.2, p)
      gs.names.2 <<- c (gs.names.2, n)
      gs.descs.2 <<- c (gs.descs.2, d)
    }
    
    k <- 1
    for (i in 1:N.gs) {
      temp <- strsplit(gs.names[i], split="_") 
      body <- paste(temp[[1]][seq(1, length(temp[[1]]) -1)], collapse="_")
      suffix <- tail(temp[[1]], 1)
      print(paste("i:", i, "gene set:", gs.names[i], "body:", body, "suffix:", suffix))
      if (suffix == "UP") {  # This is an "UP" gene set
        initial.up.entries <- initial.up.entries + 1
        target <- paste(body, "DN", sep="_")
        loc <- match(target, gs.names)            
        if (!is.na(loc)) {   # found corresponding "DN" gene set: create combined entry
          score <- score.matrix[i,] - score.matrix[loc,]
          pval <- sapply (1:Ns, function (k) fisher.pval (c (pval.matrix[i,k], pval.matrix[loc,k]))) # combine UP and DN p-values
          add.entry.2 (score, pval, body, paste(gs.descs[i], "combined UP & DN"))
          combined.entries <- combined.entries + 1
          if (combine.mode == "combine.add") {  # also add the "UP entry
            add.entry.2 (score.matrix[i,], pval.matrix[i,], gs.names[i], gs.descs[i])
            final.up.entries <- final.up.entries + 1
          }
        } else {   # did not find corresponding "DN" gene set: create "UP" entry
          add.entry.2 (score.matrix[i,], pval.matrix[i,], gs.names[i], gs.descs[i])
          final.up.entries <- final.up.entries + 1
        }
      } else if (suffix == "DN") { # This is a "DN" gene set
        initial.dn.entries <- initial.dn.entries + 1
        target <- paste(body, "UP", sep="_")
        loc <- match(target, gs.names)            
        if (is.na(loc)) { # did not find corresponding "UP" gene set: create "DN" entry
          add.entry.2 (score.matrix[i,], pval.matrix[i,], gs.names[i], gs.descs[i])
          final.dn.entries <- final.dn.entries + 1
        } else { # it found corresponding "UP" gene set
          if (combine.mode == "combine.add") { # create "DN" entry
            add.entry.2 (score.matrix[i,], pval.matrix[i,], gs.names[i], gs.descs[i])
            final.dn.entries <- final.dn.entries + 1
          }
        }
      } else { # This is neither "UP nor "DN" gene set: create individual entry
        add.entry.2 (score.matrix[i,], pval.matrix[i,], gs.names[i], gs.descs[i])
        other.entries <- other.entries + 1
      }
    } # end for loop over gene sets
    print(paste("initial.up.entries:", initial.up.entries))
    print(paste("final.up.entries:", final.up.entries))
    print(paste("initial.dn.entries:", initial.dn.entries))
    print(paste("final.dn.entries:", final.dn.entries))
    print(paste("other.entries:", other.entries))
    print(paste("combined.entries:", combined.entries))
    
    print(paste("total entries:", length(score.matrix.2[,1])))
  }            
  
  # Make sure there are no duplicated gene names after adding entries
  unique.gene.sets <- unique(gs.names.2)
  locs <- match(unique.gene.sets, gs.names.2)
  score.matrix.2 <- score.matrix.2[locs,]
  gs.names.2 <- gs.names.2[locs]
  gs.descs.2 <- gs.descs.2[locs]
  
  # Final count 
  tab <- as.data.frame(table(gs.names.2))
  ind <- order(tab[, "Freq"], decreasing=T)
  tab <- tab[ind,]
  print(tab[1:20,])
  print(paste("Total gene sets:", length(gs.names.2)))
  print(paste("Unique gene sets:", length(unique(gs.names.2))))
  
  V.GCT <- data.frame(score.matrix.2)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- gs.names.2
  write.gct.ssgsea(gct.data.frame=V.GCT, descs=gs.descs.2, filename=paste (output.prefix, '.gct', sep=''))  

  P.GCT <- data.frame(pval.matrix.2)
  names(P.GCT) <- sample.names
  row.names(P.GCT) <- gs.names.2
  write.gct.ssgsea(gct.data.frame=P.GCT, descs=gs.descs.2, filename=paste (output.prefix, '-pvalues.gct', sep=''))  
  
  if (fdr.pvalue) {
    F.GCT <- P.GCT
    if (global.fdr) 
      F.GCT <- matrix ( p.adjust(unlist (F.GCT), method='fdr'), 
                        ncol=ncol(F.GCT))
    else 
      for (i in 1:ncol(F.GCT)) F.GCT[,i] <- p.adjust (F.GCT[,i], method='fdr')
    write.gct.ssgsea(gct.data.frame=F.GCT, descs=gs.descs.2, filename=paste (output.prefix, '-fdr-pvalues.gct', sep=''))      
  }
} 




project.geneset <- function (
  data.array,
  gene.names,
  n.cols,
  n.rows,
  weight = 0,
  statistic = "Kolmogorov-Smirnov",  # alternatives: "Kolmogorov-Smirnov", "area.under.RES"
  gene.set,
  nperm = 200,
  correl.type  = "rank"              # "rank", "z.score", "symm.rank"
) { 
  ## fast implementation of GSEA
  
  ES.vector <- vector(length=n.cols)
  NES.vector <- vector(length=n.cols)
  p.val.vector <- vector(length=n.cols)
  correl.vector <- vector(length=n.rows, mode="numeric")
  
  # Compute ES score for signatures in each sample
  phi <- array(0, c(n.cols, nperm))
  for (sample.index in 1:n.cols) {
    gene.list <- order(data.array[, sample.index], decreasing=T)                
    gene.set2 <- match(gene.set, gene.names)
    
    gsea.score <- function (ordered.gene.list) {
      # function to calculate GSEA enrichment score
      if (weight == 0) {
        correl.vector <- rep(1, n.rows)
      } else if (weight > 0) {
        if (correl.type == "rank") {
          correl.vector <- data.array[ordered.gene.list, sample.index]
        } else if (correl.type == "symm.rank") {
          correl.vector <- data.array[ordered.gene.list, sample.index]
          correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows/2)], 
                                  correl.vector,
                                  correl.vector + correl.vector - correl.vector[ceiling(n.rows/2)]) 
        } else if (correl.type == "z.score") {
          x <- data.array[ordered.gene.list, sample.index]
          correl.vector <- (x - mean(x))/sd(x)
        }
      }
      tag.indicator <- sign(match(ordered.gene.list, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
      no.tag.indicator <- 1 - tag.indicator 
      N <- length(ordered.gene.list) 
      Nh <- length(gene.set2) 
      Nm <-  N - Nh 
      orig.correl.vector <- correl.vector
      if (weight == 0) correl.vector <- rep(1, N)   # unweighted case
      ind = which(tag.indicator==1)
      correl.vector <- abs(correl.vector[ind])^weight
      
      
      sum.correl = sum(correl.vector)
      up = correl.vector/sum.correl     # "up" represents the peaks in the mountain plot
      gaps = (c(ind-1, N) - c(0, ind))  # gaps between ranked pathway genes
      down = gaps/Nm
      
      RES = cumsum(c(up,up[Nh])-down)
      valleys = RES[1:Nh]-up
      
      max.ES = max(RES)
      min.ES = min(valleys)
      
      if( statistic == "Kolmogorov-Smirnov" ){
        if( max.ES > -min.ES ){
          ES <- signif(max.ES, digits=5)
          arg.ES <- which.max(RES)
        } else{
          ES <- signif(min.ES, digits=5)
          arg.ES <- which.min(RES)
        }
      }
      
      if( statistic == "area.under.RES"){
        if( max.ES > -min.ES ){
          arg.ES <- which.max(RES)
        } else{
          arg.ES <- which.min(RES)
        }
        gaps = gaps+1
        RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)
        ES = sum(RES)
      }
      gsea.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator)
      return (gsea.results)
    }
    
    GSEA.results <- gsea.score (gene.list)
    ES.vector[sample.index] <- GSEA.results$ES
    
    if (nperm == 0) {
      NES.vector[sample.index] <- ES.vector[sample.index]
      p.val.vector[sample.index] <- 1
    } else {
      for (r in 1:nperm) {
        reshuffled.gene.labels <- sample(1:n.rows)
        
        GSEA.results = gsea.score (reshuffled.gene.labels)
        phi[sample.index, r] <- GSEA.results$ES
      }
      
      if (ES.vector[sample.index] >= 0) {
        pos.phi <- phi[sample.index, phi[sample.index, ] >= 0]
        if (length(pos.phi) == 0) pos.phi <- 0.5
        pos.m <- mean(pos.phi)
        NES.vector[sample.index] <- ES.vector[sample.index]/pos.m
        s <- sum(pos.phi >= ES.vector[sample.index])/length(pos.phi)
        p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
      } else {
        neg.phi <-  phi[sample.index, phi[sample.index, ] < 0]
        if (length(neg.phi) == 0) neg.phi <- 0.5 
        neg.m <- mean(neg.phi)
        NES.vector[sample.index] <- ES.vector[sample.index]/abs(neg.m)
        s <- sum(neg.phi <= ES.vector[sample.index])/length(neg.phi)
        p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
      }
    }
  }
  
  return(list(ES.vector = ES.vector, NES.vector =  NES.vector, p.val.vector = p.val.vector))
}



##
## Support Functions
##

MSIG.Gct2Frame <- function(filename = "NULL") { 
  # Reads a gene expression dataset in GCT format and converts it into an R data frame  
  ds <- read.delim (filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T)
  descs <- ds[,1]
  ds <- ds[-1]
  row.names <- row.names(ds)
  names <- names(ds)
  return (list (ds = ds, row.names = row.names, descs = descs, names = names))
}


Read.GeneSets.db <- function (gs.db, thres.min = 2, thres.max = 2000, gene.names = NULL) {
  ## read gmt files
  temp <- readLines(gs.db)
  max.Ng <- length(temp)
  temp.size.G <- vector(length = max.Ng, mode = "numeric") 
  for (i in 1:max.Ng) {
    temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
  }
  max.size.G <- max(temp.size.G)      
  gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
  temp.names <- vector(length = max.Ng, mode = "character")
  temp.desc <- vector(length = max.Ng, mode = "character")
  gs.count <- 1
  for (i in 1:max.Ng) {
    gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
    gene.set.name <- gs.line[1] 
    gene.set.desc <- gs.line[2] 
    gene.set.tags <- vector(length = gene.set.size, mode = "character")
    for (j in 1:gene.set.size) {
      gene.set.tags[j] <- gs.line[j + 2]
    }
    if (is.null(gene.names)) {
      existing.set <- rep(TRUE, length(gene.set.tags))
    } else {
      existing.set <- is.element(gene.set.tags, gene.names)
    }
    set.size <- length(existing.set[existing.set == T])
    if ((set.size < thres.min) || (set.size > thres.max)) next
    temp.size.G[gs.count] <- set.size
    gs[gs.count,] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
    temp.names[gs.count] <- gene.set.name
    temp.desc[gs.count] <- gene.set.desc
    gs.count <- gs.count + 1
  }
  Ng <- gs.count - 1
  gs.names <- vector(length = Ng, mode = "character")
  gs.desc <- vector(length = Ng, mode = "character")
  size.G <- vector(length = Ng, mode = "numeric") 
  
  gs.names <- temp.names[1:Ng]
  gs.desc <- temp.desc[1:Ng]
  size.G <- temp.size.G[1:Ng]
  
  return(list(N.gs = Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc, 
              size.G = size.G, max.N.gs = max.Ng))
}


write.gct.ssgsea <- function(gct.data.frame, descs = "", filename) 
{
  f <- file(filename, "w")
  cat("#1.2", "\n", file = f, append = TRUE, sep = "")
  cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
  cat("Name", "\t", file = f, append = TRUE, sep = "")
  cat("Description", file = f, append = TRUE, sep = "")
  
  names <- names(gct.data.frame)
  cat("\t", names[1], file = f, append = TRUE, sep = "")
  
  if (length(names) > 1) {
    for (j in 2:length(names)) {
      cat("\t", names[j], file = f, append = TRUE, sep = "")
    }
  }
  cat("\n", file = f, append = TRUE, sep = "\t")
  
  oldWarn <- options(warn = -1)
  m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] +  2)
  m[, 1] <- row.names(gct.data.frame)
  if (length(descs) > 1) {
    m[, 2] <- descs
  } else {
    m[, 2] <- row.names(gct.data.frame)
  }
  index <- 3
  for (i in 1:dim(gct.data.frame)[2]) {
    m[, index] <- gct.data.frame[, i]
    index <- index + 1
  }
  write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
  close(f)
  options(warn = 0)
}




##
## Usage Example:
#  Input dataset: example-data.gct
#  Outputs: example-pathways.gct and example-pvalues.gct
#  Genesets: MSigDB C2 (curated pathways from KEGG, BioCarta and Reactome)
#
## R Code:
#    > source ('ssGSEA.R')
#    > ssGSEA ('example-data.gct', 'example')
#



