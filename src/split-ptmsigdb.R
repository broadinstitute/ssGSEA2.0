#create separate databases for PTMSigDB based on the main 4 categories: PATH, PERT, KINASE, and DISEASE
library(cmapR)
library(glue)

#database directory
db.dir <- "../db/ptmsigdb/v2.0.0/all"

#get all files
gmt.files <- list.files(path=db.dir,pattern=".gmt")

#for each file, split and save into correct directory
sub.classes <- c("PATH","PERT","KINASE","DISEASE")
for(f in gmt.files){
  gmt <- parse_gmt(file.path(db.dir,f))
  for(s in sub.classes){
    sub.gmt <- gmt[grepl(s,names(gmt))]
    fn <- gsub(".all",glue(".{s}"),f)
    sub.dir <- gsub("all",s,db.dir)
    write_gmt(sub.gmt,file.path(sub.dir,fn))
  }
}