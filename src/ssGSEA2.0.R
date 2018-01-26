## #######################################################################################################
## 20161013 modified by Karsten Krug
## adapt the ssGSEA code to:
##   1) work with site specific signature sets
##   2) take directionality of regulation into account
##   3) multi-threaded using 'doParallel'
##   4) handle missing values


## Single sample GSEA
## Original code written by Pablo Tamayo. Adapted with additional modifications
## Referece:
##    1. Abazeed, M. E., Adams, D. J., Hurov, K. E., Tamayo, P., Creighton, C. J., Sonkin, D., et al. (2013).
##       Integrative Radiogenomic Profiling of Squamous Cell Lung Cancer. Cancer Research, 73(20), 6289–6298.
##       http://doi.org/10.1158/0008-5472.CAN-13-1616
##    2. Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005).
##       Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.
##       Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545–15550.

#source('gct-io.r')
require(pacman)
p_load_current_gh('cmap/cmapR')
require(cmapR)

ssGSEA2 <- function (
                     input.ds,                      ## input data file in gct format, first column (Name) must contain gene symbols
                     output.prefix,                 ## prefix used for output tables
                     gene.set.databases=gsea.dbs,   ## list of genesets (in gmt format) to evaluate enrichment on
                     sample.norm.type=c("rank", "log", "log.rank", "none"),  ## sample normalization
                     weight= 0,                     ## when weight==0, all genes have the same weight; if weight>0,
                                                    ## actual values matter, and can change the resulting score
                     statistic           = c("area.under.RES", "Kolmogorov-Smirnov"), ## test statistic
                     output.score.type   = c("NES", "ES"),
                     nperm               = 1000,    ## number of random permutations for NES case
                     combine.mode        = c("combine.off", "combine.replace", "combine.add"),
                     ## "combine.off" do not combine *_UP and *_DN versions in a single score.
                     ## "combine.replace" combine *_UP and *_DN versions in a single score.
                     ## "combine.add" combine *_UP and *_DN versions in a single score and
                     ##    add it but keeping the individual *_UP and *_DN versions.
                     min.overlap         = 10,
                     correl.type         = c("rank", "z.score", "symm.rank"),  ## correlation type: "rank", "z.score", "symm.rank"
                     fdr.pvalue          = TRUE,    ## output adjusted (FDR) p-values
                     global.fdr          = FALSE,   ## if TRUE calculate global FDR; else calculate FDR sample-by-sample

                     par=F,
                     spare.cores=1,
                     log.file='run.log') {
    ## #######################################################################
    ## single sample GSEA
    ## for results similar to the Java version, use: weight=0;
    ## when weight==0, sample.norm.type and correl.type do not matter;
    ## when weight > 0, the combination of sample.norm.type and correl.type
    ##  dictate how the gene expression values in input.ds are transformed
    ##  to obtain the score -- use this setting with care (the transformations
    ##  can skew scores towards +ve or -ve values)
    ##  . sample.norm.type=='none' uses actual expression values; combined
    ##    with correl.type=='rank', genes are weighted by actual values
    ##  . sample.norm.type=='rank' weights genes proportional to rank
    ##  . sample.norm.type=='log' can be used for log-transforming input data
    ##  . correl.type=='z.score' standardizes the (normalized) input values
    ##    before using them to calculate scores.

    ## ###################################################
    ## initialize log-file
    cat('##', format(Sys.time()), '\n', file=log.file)

    ## ###################################################
    ## arguments
    sample.norm.type <- match.arg( sample.norm.type)
    statistic <- match.arg(statistic)
    output.score.type <- match.arg(output.score.type)
    combine.mode <- match.arg(combine.mode)
    correl.type <- match.arg(correl.type)

    ## ###################################################
    ## Load libraries
    p_load(gtools)
    p_load(verification)
    p_load(RColorBrewer)

    ## ####################################################
    ##            import dataset
    ## ####################################################
    #dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read dataset (GCT format)
    #m <- data.matrix(dataset$ds)
    #gene.names <- dataset$row.names
    #gene.descs <- dataset$descs
    #sample.names <- dataset$names
    dataset <- try(parse.gctx(input.ds))
    
    m <- dataset@mat
    gene.names <- dataset@rid
    gene.descs <- dataset@rdesc
    sample.names <- dataset@cid
    sample.descs <- dataset@cdesc
    ##View(sample.descs)
    # remove id column. will be repeated otherwise.
    if('id' %in% colnames(sample.descs))
      sample.descs <- sample.descs[, -which(colnames(sample.descs) == 'id')]
    gct.version <- dataset@version
    gct.src <- dataset@src

    ## number of sample columns
    Ns <- ncol(m)
    ## number of genes/ptm sites
    Ng <- nrow(m)

    ## Extract input file name
    input.file.prefix <-  sub('.*/(.*)\\.gct$', '\\1', input.ds)


    ## ###################################################
    ##         import gene set databases
    ## ###################################################
    GSDB <- vector('list', length(gene.set.databases))
    names(GSDB) <- gene.set.databases

    for (gsdb in gene.set.databases)
        GSDB[[gsdb]] <- Read.GeneSets.db2(gsdb, thres.min = min.overlap, thres.max = 2000)

    for(i in 1:length(GSDB)){
        if(i == 1){
            gs <- GSDB[[i]]$gs
            N.gs <- GSDB[[i]]$N.gs
            gs.names <- GSDB[[i]]$gs.names
            gs.descs <- GSDB[[i]]$gs.desc
            size.G <-  GSDB[[i]]$size.G
        } else {
            gs <- append(gs, GSDB[[i]]$gs)
            gs.names <- append(gs.names, GSDB[[i]]$gs.names)
            gs.descs <- append(gs.descs, GSDB[[i]]$gs.desc)
            size.G <- append(size.G, GSDB[[i]]$size.G)
            N.gs <- N.gs + GSDB[[i]]$N.gs
        }
    }


    ## ####################################################
    ##            Sample normalization
    ## -take care of missing values already here
    ## ####################################################
    if (sample.norm.type == "rank") {
        m <- apply(m, 2, function(x){
            x.rank=rep(NA, length(x))
            keep.idx=which(!is.na(x))
            x.rank[keep.idx ]=rank(x[keep.idx], ties.method="average")
            return(x.rank)
        })
        m <- 10000*m/Ng

    } else if (sample.norm.type == "log.rank") {
        m <- apply(m, 2, function(x){
            x.rank=rep(NA, length(x))
            keep.idx=which(!is.na(x))
            x.rank[keep.idx ]=rank(x[keep.idx], ties.method="average")
            return(x.rank)
        })
        m <- log(10000*m/Ng + exp(1))

    } else if (sample.norm.type == "log") {
        m[m < 1] <- 1
        m <- log(m + exp(1))

    } ##else if (sample.norm.type == "none") {
        ## keep original value -- do not transform
    ##}
    tt <- Sys.time()

    ## ####################################################
    ## calculate overlap between rownames of the dataset and
    ## the each gene set
    ## - in case of PTM signatures extract
    ##   the directionality information
    ## ####################################################
    if(length(grep(';u|;d', gs[[1]])) > 0 ){  ## PTMsigDB
        gs <- lapply(gs, function(x) x[ sub(';u$|;d$','', x) %in% gene.names ])
        gs.direction <- lapply(gs, function(x) sub('^.*;(d|u)$', '\\1', x))
        gs <- lapply(gs, function(x) sub(';u$|;d$','', x))
    } else { ## MSigDB
        gs.direction <- NULL
        gene.set.direction <- NULL
        gs <- lapply(gs, function(x) intersect(x, gene.names))
    }

    ## ###########################################
    ## calculate the overlap
    size.ol.G <- sapply(gs, length)  ## original: require at least 'min.overlap' unique GENE SET members

    cat('MSigDB import: ',  file=log.file, append=T)
    cat(Sys.time()-tt, '\n',  file=log.file, append=T)

    ## ###########################################
    ## remove gene sets with unsufficient overlap
    ## index of gene sets to test
    keep.idx <- which(size.ol.G >= min.overlap)

    ## ###########################################
    ## stop if no overlap has been found
    if(length(keep.idx) == 0)
        stop('No overlap to any gene sets found in your data set!\nPossible reasons: 1) organism other than human; 2) wrong format of gene names/site ids!\n')

    ## #####################################################
    ## update all data
    gs.names <- gs.names[keep.idx]
    gs <- gs[keep.idx]
    gs.descs <- gs.descs[keep.idx]
    size.G <- size.G[keep.idx]
    size.ol.G <- size.ol.G[keep.idx]

    if(!is.null(gs.direction))
        gs.direction <- gs.direction[keep.idx]

    ## final number of gene sets to test
    N.gs <- length(keep.idx)

    cat('Testing', length(keep.idx), 'gene sets:\n',  file=log.file, append=T)

    ## #########################################
    ## check for redundant signature sets
    gene.set.selection <- unique(gs.names)
    locs <- match(gene.set.selection, gs.names)
    if(length(locs) < N.gs){

        gs <- gs[locs]
        gs.names <- gs.names[locs]
        gs.descs <- gs.descs[locs]
        size.G <- size.G[locs]
    }

    tt <- Sys.time()

    ## ####################################################
    ##
    ##   function executed per gene set
    ##
    ## ####################################################
    project.geneset <- function (data.array,
                             gene.names,
                             n.cols,
                             n.rows,
                             weight = 0,
                             statistic = "Kolmogorov-Smirnov",   ## alternatives: "Kolmogorov-Smirnov", "area.under.RES"
                             gene.set,
                             nperm = 200,
                             correl.type  = "rank",              ## "rank", "z.score", "symm.rank"
                             gene.set.direction=NULL,             ## direction of regulation; has to be in same order than 'gene.set'
                             min.overlap
                             ) {



        ## #############################################################################
        ##
        ##           function to calculate GSEA enrichment score
        ## - apply correlation scheme and weighting
        ## - calculate ES
        ## ############################################################################
            gsea.score <- function (ordered.gene.list, gene.set2, weight, n.rows, correl.type, gene.set.direction, data.expr) {

                ##################################################################
                ## function to calculate ES score
                score <- function(max.ES, min.ES, RES, gaps, valleys, statistic){
                    ## KM
                    if( statistic == "Kolmogorov-Smirnov" ){
                        if( max.ES > -min.ES ){
                            ES <- signif(max.ES, digits=5)
                            arg.ES <- which.max(RES)
                        } else{
                            ES <- signif(min.ES, digits=5)
                            arg.ES <- which.min(RES)
                        }
                    }
                    ## AUC
                    if( statistic == "area.under.RES"){
                        if( max.ES > -min.ES ){
                            arg.ES <- which.max(RES)
                        } else{
                            arg.ES <- which.min(RES)
                        }
                        gaps = gaps+1
                        RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES) - c(valleys,0) ) * (gaps)
                        ES = sum(RES)
                    }
                    return(list(RES=RES, ES=ES, arg.ES=arg.ES))
                } ## end function score


                ## #######################################
                ## weighting
                ## #######################################
                if (weight == 0) {

                    correl.vector <- rep(1, n.rows)

                } else if (weight > 0) {
                    ## if weighting is used (weight > 0), bring
                    ## 'correl.vector' into the same order
                    ## as the ordered gene list
                    if (correl.type == "rank") {
                        ##correl.vector <- data.array[ordered.gene.list, sample.index]
                        correl.vector <- data.expr[ordered.gene.list]

                    } else if (correl.type == "symm.rank") {
                        ##correl.vector <- data.array[ordered.gene.list, sample.index]
                        correl.vector <- data.expr[ordered.gene.list]

                        correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows/2)],
                                                correl.vector,
                                                correl.vector + correl.vector - correl.vector[ceiling(n.rows/2)])
                    } else if (correl.type == "z.score") {
                        ##x <- data.array[ordered.gene.list, sample.index]
                        x <- data.expr[ordered.gene.list]
                        correl.vector <- (x - mean(x))/sd(x)
                    }
                }

                ## length of gene list - equals number of rows in input matrix
                N = length(ordered.gene.list)

                ## #####################################
                ## directionality of the gene set
                if(!is.null(gene.set.direction)){

                    ########################################
                    ## up-regulated part
                    ########################################
                    u.idx <- which(gene.set.direction=='u')

                    ## number of 'u' features
                    Nh.u <- length(u.idx)
                    Nm.u <-  N - Nh.u

                    if(length(u.idx) > 1){

                        ## locations of 'up' features
                        tag.u <- sign( match(ordered.gene.list, gene.set2[ u.idx ], nomatch=0) )
                        ind.u = which(tag.u == 1)


                        ## extract and apply weighting
                        correl.vector.u <- correl.vector[ind.u]
                        correl.vector.u <- abs(correl.vector.u)^weight           ## weighting

                        sum.correl.u <- sum(correl.vector.u)

                        up.u <- correl.vector.u/sum.correl.u         ## steps up in th random walk
                        gaps.u <- (c(ind.u-1, N) - c(0, ind.u))      ## gaps between hits
                        down.u <- gaps.u/Nm.u                        ## steps down in the random walk

                        RES.u <- cumsum(c(up.u,up.u[Nh.u])-down.u)
                        valleys.u = RES.u[1:Nh.u]-up.u

                        max.ES.u = max(RES.u)
                        min.ES.u = min(valleys.u)

                        ## calculate final score
                        score.res <- score(max.ES.u, min.ES.u, RES.u[1:Nh.u], gaps.u, valleys.u, statistic)
                        ES.u <- score.res$ES
                        arg.ES.u <- score.res$arg.ES
                        RES.u <- score.res$RES

                    } else {
                        correl.vector.u <- rep(0, N)
                        ES.u=0
                        RES.u=0
                        ind.u=NULL
                        arg.ES.u=NA
                        up.u=0
                        down.u=0
                    }

                    ## ######################################
                    ## down-regulated part
                    ## ######################################
                    d.idx <- which(gene.set.direction=='d')
                    Nh.d <- length(d.idx)
                    Nm.d <-  N - Nh.d

                    if(length(d.idx) > 1){

                        ## locations of 'd' features
                        tag.d <- sign( match(ordered.gene.list, gene.set2[ d.idx ], nomatch=0) )
                        ind.d = which(tag.d == 1)

                        ## extract and apply weighting
                        correl.vector.d <- correl.vector[ind.d]
                        correl.vector.d <- abs(correl.vector.d)^weight           ## weighting

                        sum.correl.d <- sum(correl.vector.d)

                        up.d <- correl.vector.d/sum.correl.d
                        gaps.d <- (c(ind.d-1, N) - c(0, ind.d))
                        down.d <- gaps.d/Nm.d

                        RES.d <- cumsum(c(up.d,up.d[Nh.d])-down.d)
                        valleys.d = RES.d[1:Nh.d]-up.d

                        max.ES.d = max(RES.d)
                        min.ES.d = min(valleys.d)

                        ## calculate final score
                        score.res <- score(max.ES.d, min.ES.d, RES.d[1:Nh.d], gaps.d, valleys.d, statistic)
                        ES.d <- score.res$ES
                        arg.ES.d <- score.res$arg.ES
                        RES.d <- score.res$RES


                    } else {
                        correl.vector.d <- rep(0, N)
                        ES.d=0
                        RES.d=0
                        ind.d=NULL
                        arg.ES.d=NA
                        up.d=0
                        down.d=0
                    }
                    ## ###########################
                    ## combine the results
                    ES <- ES.u - ES.d
                    RES <- list(u=RES.u, d=RES.d)
                    arg.ES <- c(arg.ES.u, arg.ES.d)
                    ##tag.indicator <- rep(0, N)
                    correl.vector = list(u=correl.vector.u, d=correl.vector.d)
                    ##if(!is.null(ind.u))tag.indicator[ind.u] <- 1
                    ##if(!is.null(ind.d))tag.indicator[ind.d] <- -1
                    ind <- list(u=ind.u, d=ind.d)
                    step.up <- list(u=up.u, d=up.d )
                    ##                   step.down <- list(u=down.u, d=down.d )
                    step.down <- list(u=1/Nm.u, d=1/Nm.d)
                    gsea.results = list(ES = ES, ES.all = list(u=ES.u, d=ES.d), arg.ES = arg.ES, RES = RES, indicator = ind, correl.vector = correl.vector, step.up=step.up, step.down=step.down)

                    ## ##############################################################
                    ##
                    ##      original ssGSEA code without directionality
                    ##
                    ## ##############################################################

                } else { ## end  if(!is.null(gene.set.direction))

                    Nh <- length(gene.set2)
                    Nm <-  N - Nh

                    ## #####################################
                    ## match gene set to data
                    tag.indicator <- sign(match(ordered.gene.list, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag)
                    ## positions of gene set in ordered gene list
                    ind = which(tag.indicator==1)
                    ## 'correl.vector' is now the size of 'gene.set2'
                    correl.vector <- abs(correl.vector[ind])^weight
                    ## sum of weights
                    sum.correl = sum(correl.vector)

                    #########################################
                    ## determine peaks and valleys
                    ## divide correl vector by sum of weights
                    up = correl.vector/sum.correl     # "up" represents the peaks in the mountain plot
                    gaps = (c(ind-1, N) - c(0, ind))  # gaps between ranked pathway genes
                    down = gaps/Nm

                    RES = cumsum(c(up,up[Nh])-down)
                    valleys = RES[1:Nh]-up

                    max.ES = max(RES)
                    min.ES = min(valleys)

                    ## calculate final score
                    score.res <- score(max.ES, min.ES, RES[1:Nh], gaps, valleys, statistic)

                    ES <- score.res$ES
                    arg.ES <- score.res$arg.ES
                    RES <- score.res$RES

                    gsea.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = ind, correl.vector = correl.vector, step.up=up, step.down=1/Nm)
                } ## end else

                return (gsea.results)
            }
        ## end function 'gsea.score'
        ## #######################################################################################################################################################
        ##debug(gsea.score)

        ## fast implementation of GSEA)
        ES.vector <- NES.vector <- p.val.vector <- rep(NA, n.cols)
        correl.vector <- rep(NA, n.rows)

        ## Compute ES score for signatures in each sample
        phi <- array(NA, c(n.cols, nperm))

        ## list to store random walk accross samples
        random.walk <- vector('list', n.cols)

        ## locations of gene set in input data (before ranking/ordering)
        ## 'gene.names' is in the same order as the input data
        gene.names.all <- gene.names
        gene.set.all <- gene.set
        gene.set.direction.all <- gene.set.direction
        n.rows.all <- n.rows

        ## loop over columns/samples
        for (sample.index in 1:n.cols) {

            ## ######################################################
            ##         handle missing values appropriatly
            ## ######################################################

            ## reset values
            gene.names <- gene.names.all
            gene.set <- gene.set.all
            gene.set.direction <- gene.set.direction.all
            n.rows <- n.rows.all

            ## extract expression values of current sample
            data.expr <- data.array[ , sample.index]

            ## index of missing values
            na.idx <- which( is.na(data.expr) | is.infinite(data.expr) )

            ## if there are missing values ...
            if( length(na.idx) > 0){

                ## update input data + gene names
                ## those are in the same order
                data.expr <- data.expr[-na.idx]
                gene.names <- gene.names[-na.idx]

                ## update numbers
                n.rows=length(data.expr)

                ## extract gene set members present in data
                gene.set.no.na <- which(gene.set %in% gene.names)

                ## update gene sets + gene set directions
                ## those are in the same order
                gene.set <- gene.set[ gene.set.no.na ]
                gene.set.direction <- gene.set.direction[ gene.set.no.na ]

                ##cat('gene set overlap:', length(gene.set), '\n', file=log.file, append=T)
                cat('gene set overlap:', sum(gene.names %in% gene.set), '\n', file=log.file, append=T)
            }

            ##if(length(gene.set) < min.overlap){
            if(sum(gene.names %in% gene.set) < min.overlap){
                random.walk[[sample.index]] <- NA

            } else {

                ## locations of gene set in input data (before ranking/ordering)
                ## 'gene.names' is in the same order as the input data
                ## gene.set2 <- match(gene.set, gene.names)
                ## gene.set2 <- which( gene.names %in% gene.set ) ## to work with redundant gene names, e.g. phospho-data
                ## this takes care about redundant gene lists, the approach above does not return the locations of 'gene.set' in 'gene.names' but the indices of 'gene.names' present in 'gene.set'
                gene.set2 <- as.numeric( unlist(sapply( gene.set, function(x) which(gene.names == x) )))

                ## order of ranks, list is now ordered, elements are locations of the ranks in
                ## original data,
                gene.list <- order(data.expr, decreasing=T)


                ## ##############################################
                ## calculate enrichment score
                GSEA.results <- gsea.score (gene.list, gene.set2, weight, n.rows, correl.type, gene.set.direction, data.expr)
                ES.vector[sample.index] <- GSEA.results$ES

                ## store the gsea results to
                random.walk[[sample.index]] <- GSEA.results

                ## ##############################################
                ## no permutations: - ES and NES are the same
                ##                  - all p-values are 1
                if (nperm == 0) {
                    NES.vector[sample.index] <- ES.vector[sample.index]
                    p.val.vector[sample.index] <- 1

                    ## ##############################################
                    ## do permutations
                } else {

                    ES.tmp = sapply(1:nperm,  function(x) gsea.score(sample(1:n.rows), gene.set2, weight, n.rows, correl.type, gene.set.direction, data.expr)$ES)
                    phi[sample.index, ] <- unlist(ES.tmp)

                    ## #######################################################
                    ## calculate NES and p-values
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
        }
        return(list(ES.vector = ES.vector, NES.vector =  NES.vector, p.val.vector = p.val.vector, random.walk = random.walk))
    } ## end function 'project.geneset'
    ## ######################################################################################################



    ## #########################################################
    ##
    ##                 multicore version
    ##
    if(par){
        p_load(doParallel)
        p_load(foreach)

        ## register cores
        cl <- makeCluster(detectCores() - spare.cores)
        registerDoParallel(cl)

        ######################
        ## parallel loop
        tmp <-  foreach(gs.i = 1:N.gs) %dopar% {

            ##cat( (gs.i / N.gs)*100, '%\n')
            cat(names(gs)[gs.i],  '\n' , file=log.file, append=T)
            gene.overlap <- gs[[gs.i]]

            if(!is.null(gs.direction))
                gene.set.direction <- gs.direction[[gs.i]]

            if (output.score.type == "ES") {

                OPAM <- project.geneset (data.array = m, gene.names = gene.names, n.cols = Ns, n.rows= Ng, weight = weight, statistic = statistic, gene.set = gene.overlap, nperm = 0, correl.type = correl.type, gene.set.direction = gene.set.direction, min.overlap = min.overlap)

            } else if (output.score.type == "NES") {
                OPAM <- project.geneset (data.array = m, gene.names = gene.names, n.cols = Ns, n.rows= Ng, weight = weight, statistic = statistic, gene.set = gene.overlap, nperm = nperm, correl.type = correl.type, gene.set.direction = gene.set.direction, min.overlap = min.overlap)
        }
        OPAM
        }

        on.exit(stopCluster(cl))

        names(tmp) <- gs.names

        ###################################
        ## serial loop
    } else { ## end if 'par'

        tt <- Sys.time()
        tmp <- lapply(1:N.gs, function(gs.i){

            cat( gs.names[gs.i], '  ' )
            cat( (gs.i / N.gs)*100, '%\n')

            cat(names(gs)[gs.i],  '\n' , file=log.file, append=T)

            gene.overlap <- gs[[gs.i]]

            if(!is.null(gs.direction))
                gene.set.direction <- gs.direction[[gs.i]]

            if (output.score.type == "ES") {
                OPAM <- project.geneset (data.array = m, gene.names = gene.names, n.cols = Ns, n.rows= Ng, weight = weight, statistic = statistic, gene.set = gene.overlap, nperm = 1, correl.type = correl.type, gene.set.direction = gene.set.direction, min.overlap = min.overlap)

            } else if (output.score.type == "NES") {
                OPAM <- project.geneset (data.array = m, gene.names = gene.names, n.cols = Ns, n.rows= Ng, weight = weight, statistic = statistic, gene.set = gene.overlap, nperm = nperm, correl.type = correl.type, gene.set.direction = gene.set.direction, min.overlap = min.overlap)

            }
        OPAM
        })
        names(tmp) <- gs.names
    } ## end else


    ## #######################################
    ## extract scores and pvalues
    ##  and generate matrices
    tmp.pval <- lapply(tmp, function(x)x$p.val.vector)
    pval.matrix <- matrix(unlist(tmp.pval), byrow=T, nrow=N.gs)
    if (output.score.type == "ES"){
        tmp.es <- lapply(tmp, function(x)x$ES.vector)
        score.matrix <- matrix(unlist(tmp.es), byrow=T, nrow=N.gs)
    }
    if (output.score.type == "NES"){
        tmp.nes <- lapply(tmp, function(x)x$NES.vector)
        score.matrix <- matrix(unlist(tmp.nes), byrow=T, nrow=N.gs)
    }
    

    random.walk <- lapply(tmp, function(x) x$random.walk)

    cat('main loop: ')
    cat(Sys.time()-tt, '\n')

    initial.up.entries <- 0
    final.up.entries <- 0
    initial.dn.entries <- 0
    final.dn.entries <- 0
    combined.entries <- 0
    other.entries <- 0

    ## #####################################
    ## don't combine
    if (combine.mode == "combine.off") {

        score.matrix.2 <- score.matrix
        pval.matrix.2 <- pval.matrix
        gs.names.2 <- gs.names
        gs.descs.2 <- gs.descs

        ## ####################################
        ## combine replace
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
    } ## end combine results


    ## ####################################################################
    ## Make sure there are no duplicated gene set names after adding entries
    unique.gene.sets <- unique(gs.names.2)
    locs <- match(unique.gene.sets, gs.names.2)

    score.matrix.2 <- data.frame( score.matrix.2[locs, ], stringsAsFactors=F )#[locs, ]
    pval.matrix.2 <- data.frame( pval.matrix.2[locs, ], stringsAsFactors=F )#[locs, ]
    
    ## to make it work with single vector inputs
    #score.matrix.2 <- data.frame(score.matrix.2, stringsAsFactors=F)
    #pval.matrix.2 <- data.frame(pval.matrix.2, stringsAsFactors=F)
    
    gs.names.2 <- gs.names.2[locs]
    gs.descs.2 <- gs.descs.2[locs]

    ## #################################################
    ##  remove emtpy rows (gene set did not achieve sufficient
    ##  overlap in any sample column)
    ## #################################################
    locs <- which( unlist( apply(score.matrix.2, 1, function(x) sum(is.na(x))/length(x)) ) < 1 )

    score.matrix.2 <- score.matrix.2[locs, ]
    pval.matrix.2 <- pval.matrix.2[locs, ]

    gs.names.2 <- gs.names.2[locs]
    gs.descs.2 <- gs.descs.2[locs]

    # ###############################################
    # prepare for new gct export (R CMAP functions)
    gs.descs.2 <- data.frame(Description=gs.descs.2, stringsAsFactors = F)
    #if(gct.version == '#1.2'){ 
    #  gs.descs.2 <- data.frame(id=gs.names.2, Description=gs.descs.2, stringsAsFactors = F)
    #} else {
    #  gs.descs.2 <- data.frame(gs.descs.2, stringsAsFactors = F)
    #  colnames(gs.descs.2) <- 
    #  }
    ## #################################################
    ## Final count
    print(paste("Total gene sets:", length(gs.names.2)))
    print(paste("Unique gene sets:", length(unique(gs.names.2))))

    #################################################
    ## score matrix
    V.GCT <- new('GCT')
    V.GCT@mat <- data.matrix(score.matrix.2)
    V.GCT@rid <- gs.names.2
    V.GCT@cid <- sample.names
    V.GCT@rdesc <- gs.descs.2
    #View(sample.descs)
    if(nrow(sample.descs) > 0)
      V.GCT@cdesc <- data.frame(sample.descs, stringsAsFactors = F)
    #else  
    #  V.GCT@cdesc <- NULL
    V.GCT@src <- gct.src
    V.GCT@version <- gct.version
    write.gct(V.GCT, paste (output.prefix, '_', paste(dim(V.GCT@mat), collapse ='x'),'.gct', sep=''))
    
    #V.GCT <- data.frame(score.matrix.2)
    #colnames(V.GCT) <- sample.names
    #rownames(V.GCT) <- gs.names.2
    
    #write.gct(gct.data.frame=V.GCT, descs=gs.descs.2, filename=paste (output.prefix, '.gct', sep=''))
    ##View(V.GCT)

    #################################################
    ## p-value matrix
    #P.GCT <- data.frame(pval.matrix.2)
    #colnames(P.GCT) <- sample.names
    #rownames(P.GCT) <- gs.names.2
    #write.gct(gct.data.frame=P.GCT, descs=gs.descs.2, filename=paste (output.prefix, '-pvalues.gct', sep=''))
    P.GCT <- new('GCT')
    P.GCT@mat <- data.matrix(pval.matrix.2)
    P.GCT@rid <- gs.names.2
    P.GCT@cid <- sample.names
    P.GCT@rdesc <- gs.descs.2
    if(nrow(sample.descs) > 0)
      P.GCT@cdesc <- data.frame(sample.descs, stringsAsFactors = F)
    P.GCT@src <- gct.src
    P.GCT@version <- gct.version
    write.gct(P.GCT, paste (output.prefix, '-pvalues' ,'_', paste(dim(P.GCT@mat), collapse ='x'),'.gct', sep=''))
    
    
    ##################################################
    ## p-value correction
    if (fdr.pvalue) {
        #F.GCT <- P.GCT
        F.GCT.mat <- P.GCT@mat
        if (global.fdr)
            F.GCT.mat <- matrix ( p.adjust(unlist (F.GCT.mat), method='fdr'),
                             ncol=ncol(F.GCT.mat))
        else
            for (i in 1:ncol(F.GCT.mat))
                F.GCT.mat[,i] <- p.adjust (F.GCT.mat[,i], method='fdr')

        ## export
        #write.gct(gct.data.frame=F.GCT, descs=gs.descs.2, filename=paste (output.prefix, '-fdr-pvalues.gct', sep=''))
        F.GCT <- new('GCT')
        F.GCT@mat <- F.GCT.mat
        F.GCT@rid <- gs.names.2
        F.GCT@cid <- sample.names
        F.GCT@rdesc <- gs.descs.2
        if(nrow(sample.descs) > 0)
            F.GCT@cdesc <- data.frame(sample.descs, stringsAsFactors = F)
        F.GCT@src <- gct.src
        F.GCT@version <- gct.version
        write.gct(F.GCT, paste (output.prefix, '-fdr-pvalues' ,'_', paste(dim(F.GCT@mat), collapse ='x'),'.gct', sep=''))
            
    }

    return(random.walk)
}



##########################################################################################
## Support Functions
##
##
##########################################################################################


#############################################################
##
## Reads a gene expression dataset in GCT 1.2 format and converts
## it into an R data frame
## - modified by kk
#############################################################
# MSIG.Gct2Frame <- function(filename = "NULL") {
# 
#   
#     ##KK 20161208
#     ##ds <- read.delim (filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T)
#     ds <- read.delim (filename, header=T, sep="\t", skip=2, row.names=NULL, blank.lines.skip=T, comment.char="", as.is=T)
#     ##ds <- read.delim(filename, skip=2, row.names=NULL, blank.lines.skip=T, stringsAsFactors = F)
#   
#     names <- colnames(ds)[3:ncol(ds)]
#     
#     descs <- ds[,2]
#     ##row.names <- row.names(ds)
#     row.names <- ds[,1]
# 
#     ds <- data.frame(ds[,-c(1,2)], stringsAsFactors=F)
# 
#   
#   return (list (ds = ds, row.names = row.names, descs = descs, names = names))
# }

#############################################################
##
##  import gene sets into R workspace
##
## 20161013 modified by kk
#############################################################
Read.GeneSets.db2 <- function (gs.db, thres.min = 2, thres.max = 2000) {
    ## read gmt files
    temp <- readLines(gs.db)
    temp <- strsplit(temp, '\t')
    temp.size.G <- sapply(temp, function(x) length(x)-2)

    ## filter gene sets according to size
    rm.idx <- which(temp.size.G < thres.min | temp.size.G > thres.max)
    if(length(rm.idx) > 0){
        temp <- temp[-rm.idx]
        temp.size.G <- temp.size.G[-rm.idx]
    }

    max.Ng <- length(temp)         ## number of signature sets
    temp.size.G <- sapply(temp, function(x) length(x)-2)
    max.size.G <- max(temp.size.G) ## maximal size

    gs <- lapply(temp, function(x)x[3:length(x)])
    gs.names <- sapply(temp, function(x)x[1])
    gs.desc <- sapply(temp, function(x)x[2])
    size.G <- temp.size.G
    names(gs) <- names(gs.names) <- names(gs.desc) <- names(size.G) <- gs.names


  return(list(N.gs = max.Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc,
              size.G = size.G, max.N.gs = max.Ng))
}
######################################################################
##
##           export dataframe to gct 1.2 format
## 20170302 modified by kk
######################################################################
# write.gct <- function(gct.data.frame, names = NULL, descs = NULL, filename)
# {
#     gct <- c()
# 
#     ## #################################
#     ##       add names and desc
#     if(is.null(names))
#         names <- rownames(gct.data.frame)
# 
#     df <- data.frame(Name=names, Description=descs, gct.data.frame)
# 
#     ## #################################
#     ##       header
#     gct[1] <- "#1.2"
#     gct[2] <- paste(dim(gct.data.frame), collapse='\t')
#     gct[3] <- paste( colnames(df) , collapse='\t')
# 
#     ## ##################################
#     ##      data frame
#     gct <- c(gct, unlist(apply(df, 1, paste, collapse='\t')))
# 
#     ## ###################################
#     ##        export
#     writeLines(gct, con=filename)
# 
# }

#################################################
##   Given a string and a number of characters
##   the function chops the string to the
##   specified number of characters and adds
##   '...' to the end.
## parameter
##   string     - character
##   nChar      - numeric
## value
##   string of 'nChar' characters followed
##     by '...'
##################################################
chopString <- function(string, nChar=10, add.dots=T)
{
    string.trim <- strtrim(string, nChar)

    if(add.dots)
        string.trim[ which(nchar(string) > nChar) ] <-  paste(string.trim[which(nchar(string) > nChar) ], '...')
    if(!add.dots)
        string.trim[ which(nchar(string) > nChar) ] <-  paste(string.trim[which(nchar(string) > nChar) ])

    return(string.trim)
}

##########################################################################################################
## translate a color name into rgb space
##
## changelog:  20100929 implementation
##########################################################################################################
my.col2rgb <- function(color, alpha=80, maxColorValue=255){

    out <- vector( "character", length(color) )

    for(col in 1:length(color)){

        col.rgb <- col2rgb(color[col])

        out[col] <- rgb(col.rgb[1], col.rgb[2], col.rgb[3], alpha=alpha, maxColorValue=maxColorValue)

    }
    return(out)
}
