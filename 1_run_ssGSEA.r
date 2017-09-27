################################################################################################################
## Filename: 1_run_ssGSEA.r
## Created: September 09, 2017
## Author(s): Karsten Krug
##
## Purpose: 
##      - Wrapper to ssGSEA script to perform single sample Gene Set Enrichment analysis.
##          
## Instructions:  
##      - Just source the script into an R-session:
##          - RStudio: open file and press 'Source' in the upper right part of the edityor window 
##          - R-GUI: drap and drop this file into an R-GUI windoe
##      - In order to specify your input files and databases the script will invoke two 
##        Windows file dialogs.
##      - The first dialog let's you choose a folder containing input files in GTC v1.2 format. 
##        The script will loop over all gct files in this directory and run ssGSEA on each file 
##        separately.
##      - The second dialog window let's the user choose a gene set database such as MSigDB. 
##        Some default database can be found in the 'db' subfolder. 
##
################################################################################################################
rm(list=ls())
script.dir <- dirname(sys.frame(1)$ofile) ## get folder the script

## ##########################################################
##  define parameters below:
## ##########################################################

## ##############################
## sGSEA / PSEA parameters
sample.norm.type    = "rank"
weight              = 0.75
statistic           = "area.under.RES"    ##"Kolmogorov-Smirnov"
output.score.type   = "NES"               ## 'ES' or 'NES'
nperm               = 1e3                 ## No of permutations
min.overlap         = 10                  ## gene centric 10
correl.type         = "z.score"
par                 = T
spare.cores         = 1

## #####################################################################
##   end paramaters
## - in a perfect world users don't have to worry about the stuff below...
## #####################################################################
## directory with gct files
##gct.dir <- "C:/Users/Karsten/Dropbox/Manuscripts/20170630_PTM-GSEA_manuscript/CCLE/"

gct.dir.ok=F
while(!gct.dir.ok){
  gct.dir <- choose.dir(default='.', caption = 'Choose directory containing gct files. As of now only gst v1.2 files are supported.')
  if(length(grep('\\.gct', dir(gct.dir))) > 0)
    gct.dir.ok=T
  }

## directory to write output
out.dir <- gct.dir

## MSigDB
##gsea.db = "h.all.v6.0.symbols.gmt" ## gene centric
gene.set.databases = choose.files(default = paste( script.dir, 'db/c2.cp.v6.0.symbols.gmt', sep='/' ), caption='Choose gene set database in gmt format. See Broad\'s MSigDB website for details.')

##gsea.db = 'C:/Users/Karsten/Dropbox/Manuscripts/20170630_PTM-GSEA_manuscript/db/ptm.sig.db.all.uniprot.human.v1.3.gmt' ## site centric
##gsea.db <- 'H:/Projects/PTM_GSEA/20170920_PTMsigDB_v1.5/ptm.sig.db.all.uniprot.human.v1.5.gmt'



## ######################################################################
##                          START
## ######################################################################
##gene.set.databases  = ifelse(length(grep('/',gsea.db)) == 0, paste(script.dir,'db', gsea.db, sep='/'), gsea.db)

source(paste(script.dir, 'src/ssGSEA_PSEA.R', sep='/'))

## #############################################
## prepare output folder
setwd(out.dir )

date.str <- paste(sub(' .*', '', Sys.time()), sep='_')
dir.create(date.str)
setwd(date.str)


## #############################################
## import signature database
signat.all <- readLines(gene.set.databases)
signat.all <- strsplit(signat.all, '\t')
names(signat.all) <- sapply(signat.all, function(x)x[1])
signat.all <- lapply(signat.all, function(x) x[-c(1,2)])



## save parameters used for ssGSEA
param.str = c(
    paste('##', Sys.time()),
    paste('gct.directory:', gct.dir, sep='\t'),
    paste('output.directory:', out.dir, sep='\t'),
    paste('gene.set.database:',gene.set.databases, sep='\t'),
    paste('sample.norm.type:', sample.norm.type, sep='\t'),
    paste('weight:', weight, sep='\t'),
    paste('statistic:', statistic, sep='\t'),
    paste('output.score.type', output.score.type, sep='\t'),
    paste('nperm:', nperm, sep='\t'),
    paste('min.overlap:', min.overlap, sep='\t'),
    paste('correl.type:', correl.type, sep='\t'),
    paste('run.parallel:', par, sep='\t')
   )
writeLines(param.str, con='parameters.txt')


## identify all gct files
gct.files <- dir(gct.dir, pattern='\\.gct$', full.names=T)
##names(gct.files) <- sub('_gene_centric.*','', sub('\\.gct$', '', sub('.*/','', gct.files)))
names(gct.files) <- paste(  sub('\\.gct$', '', sub('.*/','', gct.files)), 'ssGSEA', sep='_' )


## #####################################
## loop over gct files and run ssGSEA
for(i in names(gct.files)){


    ## create sub folders if more than one gct file was found
    if(length(gct.files) > 1){
        subdir=sub(',|\\.|:|;|/', '_', i)
        dir.create(subdir)
        setwd(subdir)
    }

    ## ########################################
    ## ssGSEA

    ## input data set
    input.ds <- gct.files[i]

    cat('Running ssSGEA on:', sub('.*/', '', input.ds), '\n\n')

    ## run ssGSEA
    gsea.res <- ssGSEA2(input.ds, gene.set.databases=gene.set.databases, sample.norm.type=sample.norm.type, weight=weight,statistic=statistic, output.score.type = output.score.type, nperm  = nperm, min.overlap  = min.overlap, correl.type = correl.type, output.prefix = paste(i), par=par, spare.cores=spare.cores )

    ## save object
    save(gsea.res, file=paste(i, '.RData', sep=''))

    ## #########################################################
    ##              rank plots
    ## #########################################################
    
    ## input dataset
    input = read.delim(input.ds, stringsAsFactors=F, skip=2, row.names=NULL)

    ## gene/site ids
    gn.input <- input[, 1]

    ## expression data only
    input <- input[, -c(1,2)]

    ## sample names
    all.samp <- names(input)

    ## import enrichment scores and p-values
    gsea.score <- read.delim(paste( i, '.gct', sep=''), stringsAsFactors=F, skip=2, row.names='Name')
    gsea.pval <-  read.delim(paste( i, '-pvalues.gct', sep=''), stringsAsFactors=F, skip=2, row.names='Name')

    ## gene set names
    all.gs <- rownames(gsea.score)

    ## keep only scored signatures
    signat <- signat.all[all.gs]

    ## create sub-folder
    dir.create('rank-plots')
    
    ## loop over gene sets
    for(gs in 1:length(all.gs)){

        gs.name <- all.gs[gs]

        pdf(paste('rank-plots/', chopString( gsub('\\:|\\/\\\t', ' ', gs.name), nChar=20, add.dots=F) ,'_2.pdf', sep=''), 9.5, 9.5)
        par(mfrow=c(3, 3))
        for(samp in 1:length(all.samp)){

            ## extract results
            samp.name <- all.samp[samp]

            ## gsea results
            score <- gsea.score[gs.name, samp.name]
            pval <- gsea.pval[gs.name, samp.name]

            ## extract data
            data.expr <- input[, samp.name ]

            valid.idx <- which( !(is.na( data.expr ) | is.infinite(data.expr)) )

            data.expr <- data.expr[ valid.idx ]
            gn <- gn.input[ valid.idx ]

            ## order
            ord.idx <- order(data.expr, decreasing=T)

            ##gn <- row.names(input)[ord.idx]
            gn <- gn[ ord.idx ]
            data.expr <- data.expr[ ord.idx ]


            plot( data.expr, pch=20, col='darkgrey', lwd=4, type='l', xlab='Rank', ylab='Expression', main=paste(gs.name, samp.name, sep='\n'), ylim=range(data.expr), yaxs='i')
						abline(h=0, lty='dashed', lwd=2, col='grey70')

            ## #########################################################
            ##  ptm signatures?
            if(length(grep(';u$|;d$', signat[[gs.name]], value=T)) > 0){

                ## locations
                gsea.tmp.u <- sub(';u$','',grep(';u$', signat[[gs.name]], value=T))
                loc.u <- na.omit(match(gsea.tmp.u, gn))

                gsea.tmp.d <- sub(';d$','',grep(';d$',  signat[[gs.name]], value=T))
                loc.d <- na.omit(match(gsea.tmp.d, gn))

                if(!is.null(loc.u)){

                    rug(loc.u, col='darkred', side=3, lwd=3, ticksize=0.02)
                    points(loc.u, data.expr[loc.u], col=my.col2rgb('darkred',  150), pch=16, cex=2)
                }
                if(!is.null(loc.d)){
                    rug(loc.d, col='darkblue', side=1, lwd=3, ticksize=0.02)
                    points(loc.d, data.expr[loc.d], col=my.col2rgb('darkblue',  150), pch=16, cex=2)

                }
                ## some info
                legend('bottom', legend=c(paste('No. down-regulated in signature:', length(grep(';d$', signat[[gs.name]]))),
                                          paste('No. found in data set:', length(loc.d))
                                          ), inset=.05, bty='n', text.col='darkblue')

                legend('top', legend=c(paste('No. up-regulated in signature:', length(grep(';u$', signat[[gs.name]]))),
                                       paste('No. found in data set:', length(loc.u))
                                       ), inset=.05, bty='n', text.col='darkred')
            } else {## end if signature

                ## ####################################################
                ## regular gene set
                loc <- which(gn %in% signat[[gs.name]])
                rug(loc, col=my.col2rgb('darkred',  50), side=3, lwd=2, ticksize=0.02)
								points(loc, data.expr[loc], col=my.col2rgb('darkred',  150), pch=16, cex=2)

                ## box plot
                loc.quart <- quantile(loc)
                rug(loc.quart, col='darkblue', side=3, lwd=2, ticksize=0.03)
                rect( loc.quart[2], max(data.expr)-0.04*max(data.expr-min(data.expr)), loc.quart[4], max(data.expr), border='darkblue', lwd=2, col=NA )
                rect( loc.quart[1], max(data.expr)-0.02*max(data.expr-min(data.expr)), loc.quart[2], max(data.expr)-0.02*max(data.expr-min(data.expr)), border='darkblue', lwd=2 )
                rect( loc.quart[4], max(data.expr)-0.02*max(data.expr-min(data.expr)), loc.quart[5], max(data.expr)-0.02*max(data.expr-min(data.expr)), border='darkblue', lwd=2 )

                ## some info
                legend('bottom', legend=c(paste('No. in signature:', length( signat[[gs.name]])),
                                          paste('No. found in data set (non-redund.):', sum(signat[[gs.name]] %in% gn)),
                                          paste('No. found in data set (redundant):', length(loc))
                                       ), inset=.05, bty='n', text.col='darkred')

            }
            legend('right', legend=paste('NES=', round(score, 3), ' (p=', round(pval, 5), ')', sep=''), bty='n', inset=.2, cex=1.5)


            ##ssGSEA.plot( data.expr, gs.name, samp.name, signat, score, pval )

	}
	par(mfrow=c(1, 1))
	dev.off()
    } ## end loop over gene sets

    if(length(gct.files) > 1)
        setwd('..')
}






