
## From https://github.com/cmap/l1ktools/tree/master/R/cmap/io.R
## For reading/writing gct v1.2, v1.3 and gctx (binary HDF5) files
## With modifications ...


## install and load libraries automatically
if (!require("pacman")) install.packages ("pacman")
pacman::p_load (rhdf5) 

# load dependencies
# library(rhdf5)

## NOTE: If using functions from io.r, SOURCE THIS FILE FIRST

########################################
### GCT class and method definitions ###
########################################

setClass("GCT",
         representation(
           mat = "matrix",
           rid = "vector",
           cid = "vector",
           rdesc = "data.frame",
           cdesc = "data.frame",
           version = "character",
           src = "character"
         )
)


#### define some helper methods for parsing gctx files ###

# helper function to set all the row and column annotations to the correct data type
fix.datatypes <- function(meta) {
  # turn all warnings to errors so we can use the try statement to grab strings
  options(warn = 2)
  for (field.name in names(meta)) {
    # get the field (convert to string to deal with factors)
    field = as.character (meta[[field.name]])
    # check if it's numeric
    try({field = as.numeric(field)}, silent = TRUE)
    if (is.numeric(field)) {
      #check if it's an integer
      int.field = NULL
      try({int.field = as.integer(field)}, silent = TRUE)
      if ( ! is.null(int.field) && identical(int.field, field) )
        field = int.field
    }
    # insert back into the annotations
    meta[[field.name]] = field
  }
  options(warn = 0)
  return(meta)
}


# helper function for parsing row or column metadata
read.gctx.meta <- function(gctx_path, dimension="row", ids=NULL, set_annot_rownames=T) {
  if (!(dimension %in% c("row", "col"))) {
    stop("dimension can be either row or col")
  }
  if (dimension == "row") {
    name <- "0/META/ROW"
  } else {
    name <- "0/META/COL"
  }
  raw_annots <- h5read(gctx_path, name=name) # returns a list
  fields <- names(raw_annots)
  # define an empty data frame of the correct dimensions
  annots <-  data.frame(matrix(nrow=length(raw_annots[[fields[1]]]), ncol=length(fields)))
  names(annots) <-  fields
  # loop through each field and fill the annots data.frame
  for (i in 1:length(fields)) {
    field <- fields[i]
    # remove any trailing spaces
    annots[,i] <- gsub("\\s*$", "", raw_annots[[field]], perl=T)
  } 
  annots <- fix.datatypes(annots)
  # subset to the provided set of ids, if given
  if (is.null(ids)) {
    ids <- as.character(annots$id)
  } else {
    ids <- ids
  }
  # make sure annots row ordering matches that of ids
  annots <- subset_to_ids(annots, ids)
  # use the id field to set the rownames
  if (set_annot_rownames) {
    rownames(annots) <- as.character(annots$id)
  }
  return(annots)
}


# helper function for reading gctx row/col ids
read.gctx.ids <- function(gctx_path, dimension="row") {
  if (!(dimension %in% c("row", "col"))) {
    stop("dimension can be either row or col")
  }
  if (dimension == "row") {
    name <- "0/META/ROW/id"
  } else {
    name <- "0/META/COL/id"
  }
  # remove any spaces
  ids <- gsub("\\s*$", "", h5read(gctx_path, name=name), perl=T)
  return(ids)
}


subset_to_ids <- function(df, ids) {
  # helper function to do a robust df subset
  check_colnames("id", df)
  newdf <- data.frame(df[match(ids, df$id), ])
  names(newdf) <- names(df)
  return(newdf)
}


check_colnames <- function(test_names, df, throw_error=T) {
  # check whether test_names are valid names in df
  # throw error if specified
  diffs <- setdiff(test_names, names(df))
  if (length(diffs) > 0) {
    if (throw_error) {
      stop(paste("the following column names are not found in", deparse(substitute(df)), ":",
                 paste(diffs, collapse=" "), "\n"))
    } else {
      return(F)
    }
  } else {
    return(T)
  }
}

# define the initialization method for the class
setMethod("initialize",
          signature = "GCT",
          definition = function(.Object, src = NULL, rid = NULL, cid = NULL, set_annot_rownames = T) {
            # create empty object if src==NULL -- use to create a new gct file
            if (is.null (src)) {
              return (.Object)
            }
            # check to make sure it's either .gct or .gctx
            if (! (grepl(".gct$", src) || grepl(".gctx$", src) ))
              stop("Either a .gct or .gctx file must be given")
            if (grepl(".gct$", src)) {
              if ( ! is.null(rid) || !is.null(cid) )
                stop("rid and cid values may only be given for .gctx files, not .gct files")
              # parse the .gct
              .Object@src = src
              # get the .gct version by reading first line
              .Object@version = scan(src, what = "", nlines = 1, sep = "\t", quiet = TRUE)[1]
              # get matrix dimensions by reading second line
              dimensions = scan(src, what = double(0), nlines = 1, skip = 1, sep = "\t", quiet = TRUE)
              nrmat = dimensions[1]
              ncmat = dimensions[2]
              if (length(dimensions)==4) {
                # a #1.3 file
                message("parsing as GCT v1.3")
                nrhd <- dimensions[3]
                nchd <- dimensions[4]
              } else {
                # a #1.2 file
                message("parsing as GCT v1.2")
                nrhd <- 1
                nchd <- 0
              }
              message(paste(src, nrmat, "rows,", ncmat, "cols,", nrhd, "row descriptors,", nchd, "col descriptors"))
              # read in header line
              header = scan(src, what = "", nlines = 1, skip = 2, sep = "\t", quote = NULL, quiet = TRUE)
              # construct row header and column id's from the header line
              if ( nrhd > 0 ) {
                rhd <- header[2:(nrhd+1)]
                cid <- header[-(nrhd+1):-1]
                col_offset <- 1
              }
              else {
                if ("Description" %in% header) {
                  # check for presence of description column in v1.2 files
                  col_offset <- 2
                }
                rhd = NULL
                cid = header[(1+col_offset):length(header)]
              }
              # read in the next set of headers (column annotations) and shape into a matrix
              if ( nchd > 0 ) {
                header = scan(src, what = "", nlines = nchd, skip = 3, sep = "\t", 
                              quote = NULL, quiet = TRUE)		
                header = matrix(header, nrow = nchd, 
                                ncol = ncmat + nrhd + 1, byrow = TRUE)
                # extract the column header and column descriptions
                chd = header[,1]
                cdesc = header[,-(nrhd+1):-1]
                # need to transpose in the case where there's only one column annotation
                if ( nchd == 1 )
                  cdesc = t(cdesc)
              }
              else {
                chd = NULL
                cdesc = data.frame()
              }
              # read in the data matrix and row descriptions, shape into a matrix
              mat = scan(src, what = "", nlines = nrmat, 
                         skip = 3 + nchd, sep = "\t", quote = NULL, quiet = TRUE)
              mat = matrix(mat, nrow = nrmat, ncol = ncmat + nrhd + col_offset, 
                           byrow = TRUE)
              # message(paste(dim(mat), collapse="\t"))
              # Extract the row id's row descriptions, and the data matrix
              rid = mat[,1]
              if ( nrhd > 0 ) {
                # need as.matrix for the case where there's only one row annotation
                rdesc = as.matrix(mat[,2:(nrhd + 1)])
                mat = matrix(as.numeric(mat[,-(nrhd + 1):-1]),
                             nrow = nrmat, ncol = ncmat)
              }
              else {
                rdesc = data.frame()
                mat = matrix(as.numeric(mat[, (1+col_offset):ncol(mat)]), nrow = nrmat, ncol = ncmat)
              }
              # assign names to the data matrix and the row and column descriptions
              # message(paste(dim(mat), collapse="\t"))
              dimnames(mat) = list(rid, cid)
              if ( nrhd > 0 ) {
                dimnames(rdesc) = list(rid,rhd)
                rdesc = as.data.frame(rdesc, stringsAsFactors = FALSE)
              }
              if ( nchd > 0 ) {
                cdesc = t(cdesc)
                dimnames(cdesc) = list(cid,chd)
                cdesc = as.data.frame(cdesc, stringsAsFactors = FALSE)
              }
              # assign to the GCT slots
              .Object@mat = mat
              .Object@rid = rownames(mat)
              .Object@cid = colnames(mat)
              .Object@rdesc = fix.datatypes(rdesc)
              .Object@cdesc = fix.datatypes(cdesc)
              # add id columns to rdesc and cdesc
              .Object@rdesc$id <- rownames(.Object@rdesc)
              .Object@cdesc$id <- rownames(.Object@cdesc)                  
              return(.Object)
            }
            else { 
              # parse the .gctx
              message(paste("reading", src))
              .Object@src = src
              # if the rid's or column id's are .grp files, read them in
              if ( length(rid) == 1 && grepl(".grp$", rid) )
                rid <- parse.grp(rid)
              if ( length(cid) == 1 && grepl(".grp$", cid) )
                cid <- parse.grp(cid)
              # get the row and column ids
              all_rid <- read.gctx.ids(src, dimension="row")
              all_cid <- read.gctx.ids(src, dimension="col")
              # if rid or cid specified, read only those rows/columns
              # if already numeric, use as is
              # else convert to numeric indices
              if (!is.null(rid)) {
                if (is.numeric(rid)) {
                  ridx <- rid
                } else {
                  ridx <- match(rid, all_rid)
                }
              } else {
                ridx <- seq_along(all_rid)
              }
              if (!is.null(cid)) {
                if (is.numeric(cid)) {
                  cidx <- cid
                } else {
                  cidx <- match(cid, all_cid)
                }
              } else {
                cidx <- seq_along(all_cid)
              }
              # subset the character ids to the ones we want
              rid_keep <- all_rid[ridx]
              cid_keep <- all_cid[cidx]
              # read the data matrix
              .Object@mat <- h5read(src, name="0/DATA/0/matrix", index=list(ridx, cidx))
              # set the row and column ids
              .Object@rid <- rid_keep
              .Object@cid <- cid_keep
              colnames(.Object@mat) <- all_cid[cidx]
              rownames(.Object@mat) <- all_rid[ridx]
              # get the meta data
              .Object@rdesc <- read.gctx.meta(src, dimension="row", ids=rid_keep, set_annot_rownames=set_annot_rownames)
              .Object@cdesc <- read.gctx.meta(src, dimension="col", ids=cid_keep, set_annot_rownames=set_annot_rownames)
              # close any open handles and return the object
              H5close()
              return(.Object)
            }
          }
)


# function to parse a GCT(X)
# just instantiates a new GCT object
parse.gctx <- function(fname, rid = NULL, cid = NULL, set_annot_rownames = T) {
  ds <- new("GCT", src = fname, rid = rid, cid = cid, set_annot_rownames = set_annot_rownames)
  return(ds)
}


append.dim <- function(ofile, mat, extension="gct") {
  nc <- ncol(mat)
  nr <- nrow(mat)
  outFile <- basename(ofile)
  filename <- strsplit(outFile,'.',fixed=T)[[1]][1]
  ofile <- file.path(dirname(ofile),
                     sprintf('%s_n%dx%d.%s',filename,
                             nc, nr, extension))
  return(ofile)
}


# subset a gct object (sample subset)
subset.gct <- function (gct, index) {
  # returns a gct object that contains only the subset specified by index
  # index must be a boolean vector with length equal to # columns of the data matrix
  gct@mat <- gct@mat [,index]
  gct@cid <- gct@cid [index]
  if (nrow (gct@cdesc) > 0) gct@cdesc <- gct@cdesc [index,]
  
  return (gct)
}


add.cols.gct <- function (gct, dx, dx.annot=NULL) {
  # returns a gct object with the provided data frame dx added to the original 
  # appended to the (end of the ) data matrix
  # if dx.annot is provided it is added as annotation for the samples in dx
  # (else all annotation columns are set to NA)
  gct@mat <- cbind (gct@mat, as.matrix (dx))
  gct@cid <- c (gct@cid, colnames (dx))
  if (nrow (gct@cdesc) > 0) {
    if (is.null (dx.annot)) {
      nc <- ncol(gct@cdesc)
      dx.annot <- data.frame (matrix (rep (NA, nc), ncol=nc))
    }
    nr <- nrow (gct@cdesc)
    gct@cdesc <- rbind (gct@cdesc, dx.annot)
    gct@cdesc [seq (from=nr+1, length.out=ncol(dx)), 'id'] <- colnames (dx)
    rownames (gct@cdesc) <- gct@cdesc [,'id']
  }
  
  return (gct)
}

# subset a gct object (gene/row subset)
row.subset.gct <- function (gct, index) {
  # returns a gct object that contains only the row subset specified by index
  # index must be a boolean vector with length equal to # rows of the data matrix
  gct@mat <- gct@mat [index,]
  gct@rid <- gct@rid [index]
  if (nrow (gct@rdesc) > 0) gct@rdesc <- gct@rdesc [index,]
  
  return (gct)
}


# rearrange a gct object columns
rearrange.gct <- function (gct, index, new.cid=NULL) {
  # returns a gct object that contains the samples (cols) rearranged as in index
  # index must be an integer vector with items indicating column numbers
  gct@mat <- gct@mat [,index]
  gct@cid <- gct@cid [index]
  if (nrow (gct@cdesc) > 0) gct@cdesc <- gct@cdesc [index,]
  if (!is.null (new.cid)) {
    new.cid <- as.character (new.cid)
    colnames (gct@mat) <- new.cid
    gct@cid <- new.cid
    if (nrow (gct@cdesc) > 0) {
      gct@cdesc[,'id'] <- new.cid
      rownames (gct.cdesc) <- new.cid
    }
  }
  
  return (gct)
}


# write a gct file to disk
write.gct <- function(ds, ofile, precision=5, appenddim=F, ver=NULL) {
  # gct must contain the following fields
  #          mat: Numeric data matrix [RxC]
  #          rid: Cell array of row ids
  #          rdesc: Cell array of row annotations
  #          cid: Cell array of column ids
  #          cdesc: Cell array of column annotations
  #          version: GCT version string
  #          src: Source filename
  # version is decided based on ds@version unless explicitly specified by ver
  
  
  # append the dimensions of the data set, if desired
  if (appenddim) ofile <- append.dim(ofile, ds@mat, extension="gct")
  
  # detect version (unless specified)
  if (is.null (ver)) ver <- ifelse (ds@version == "#1.3", 3, 2)
  
  precision = floor(precision)
  cat(sprintf('Saving file to %s\n',ofile))
  nr <- nrow(ds@mat)
  nc <- ncol(ds@mat)
  cat(sprintf('Dimensions of matrix: [%dx%d]\n',nr,nc))
  cat(sprintf('Setting precision to %d\n',precision))
  
  # open file      
  if (ver==3) {
    nrdesc = dim(ds@rdesc)[2]
    ncdesc = dim(ds@cdesc)[2]
    colkeys = colnames(ds@cdesc)
    # append header
    cat(sprintf('#1.%d\n%d\t%d\t%d\t%d', ver, nr, nc, nrdesc, ncdesc),
        file=ofile,sep='\n')      
    # line 3: sample row desc keys and sample names
    cat(paste(c('id',colnames(ds@rdesc),ds@cid),collapse='\t'),
        file=ofile,sep='\n',append=T)
    # line 4 + ncdesc: sample desc
    filler = 'na'
    for (ii in 1:ncdesc) {
      if (is.numeric(ds@cdesc[,ii])) {
        cat(paste(c(colkeys[ii],rep(filler,nrdesc),
                    round(ds@cdesc[,ii],precision)),
                  collapse='\t'),
            file=ofile,sep='\n',append=T)  
      } else {
        cat(paste(c(colkeys[ii],rep(filler,nrdesc),
                    ds@cdesc[,ii]),
                  collapse='\t'),
            file=ofile,sep='\n',append=T)
      }
    }
    
    for (ii in 1:nr) {    
      # print rows
      cat(paste(c(ds@rid[ii],
                  ds@rdesc[ii,],
                  round(ds@mat[ii,],precision)),collapse='\t'),
          sep='\n',file=ofile,append=T)
    }
  } else {
    # assume ver 1.2 and below, ignore descriptors
    # append header
    cat(sprintf('#1.%d\n%d\t%d', ver, nr, nc),
        file=ofile,sep='\n')      
    # line 3: sample row desc keys and sample names
    cat(paste(c('id','Description',ds@cid),collapse='\t'),
        file=ofile,sep='\n',append=T)
    
    for (ii in 1:nr) {    
      # print rows
      cat(paste(c(ds@rid[ii],
                  ds@rdesc[ii, 'Description'],
                  round(ds@mat[ii,],precision)),collapse='\t'),
          sep='\n',file=ofile,append=T)
    }
  }
  
  cat(sprintf('Saved.\n'))  
}


# write a GCTX object
write.gctx <- function(ds, ofile, appenddim=T, compression_level=6, matrix_only=F) {
  if (appenddim) ofile <- append.dim(ofile, ds@mat, extension="gctx")
  # check if the file already exists
  if (file.exists(ofile)) {
    message(paste(ofile, "exists, removing"))
    file.remove(ofile)
  }
  message(paste("writing", ofile))
  
  # start the file object
  h5createFile(ofile)
  
  # create all the necessary groups
  h5createGroup(ofile, "0")
  h5createGroup(ofile, "0/DATA")
  h5createGroup(ofile, "0/DATA/0")
  h5createGroup(ofile, "0/META")
  h5createGroup(ofile, "0/META/COL")
  h5createGroup(ofile, "0/META/ROW")
  
  # H5Gcreate(fid, "0")
  # H5Gcreate(fid, "0/DATA")
  # H5Gcreate(fid, "0/DATA/0")
  # H5Gcreate(fid, "0/META")
  # H5Gcreate(fid, "0/META/COL")
  # H5Gcreate(fid, "0/META/ROW")
  
  # create and write matrix data, using
  # chunking if dimensions exceed 1000
  # assume values are 32 bit (4 bytes each), so we can fit 1024 / 4 = 256 values in 1 KB (1024 bytes)
  row_chunk_size <- min(nrow(ds@mat), 1000)
  # column chunk, such that row * col <= 1024
  # should play with these values
  col_chunk_size <- min(floor(1024 / row_chunk_size), ncol(ds@mat))
  chunking <- c(row_chunk_size, col_chunk_size) 
  message(paste(c("chunk sizes:", chunking), collapse="\t"))
  h5createDataset(ofile, "0/DATA/0/matrix", dim(ds@mat), chunk=chunking, level=compression_level)
  h5write(ds@mat, ofile, "0/DATA/0/matrix")
  
  # write annotations
  h5write(ds@rid, ofile, "0/META/ROW/id")
  h5write(ds@cid, ofile, "0/META/COL/id")
  
  if (!matrix_only) {
    write.meta(ofile, ds@cdesc, dimension="column")
    write.meta(ofile, ds@rdesc, dimension="row")
  }
  
  # close any open handles
  H5close()
  
  # add the version annotation and close
  fid <- H5Fopen(ofile)
  h5writeAttribute("GCTX1.0", fid, "version")
  H5close()
  
}


# helper function to write a data.frame of meta data to gctx object
# makes an HDF5 entry for each column of the data.frame
write.meta <- function(ofile, df, dimension="row") {
  path <- if ((dimension=="row")) "0/META/ROW/" else "0/META/COL/"
  # loop through all columns
  fields <- names(df)
  if (length(fields) > 0) {
    for (i in 1:length(fields)) {
      field <- fields[i]
      v <- df[, i]
      # convert factors to character
      if(class(v) == "factor" || class(v) == "AsIs") {
        v <- as.character(v)
      }
      h5write(v, ofile, paste(path, field, sep=""))
    }
  }
}


###########################################
### functions for other CMap file types ###
###########################################

### function to read a .grp file and return a vector ###
parse.grp <- function(fname) {
  grp <- scan(fname, what = "", quote = NULL, quiet = TRUE)
  return(grp)
}


### function to write a .grp file
write.grp <- function(vals, fname) {
  if (is.list(vals)) vals <- unlist(vals)
  if (!is.vector(vals)) vals <- as.vector(vals)
  write(vals, fname, ncolumns=1)
}


### function to read a .gmx file and return a list ###
parse.gmx <- function(fname) {
  tmp <- read.table(fname, sep = "\t", 
                    header = TRUE, stringsAsFactors = FALSE)
  # preallocate a list for the gmx
  L <- list()
  # loop over the first row of the .gmx
  for ( n in names(tmp) ) {
    # get all the values; remove empties at the end
    values <- tmp[[n]][-1]
    remove.idx <- values == ""
    values <- values[!remove.idx]
    # put in a list
    L[[n]] <- list(head = n,
                   desc = tmp[[n]][1], 
                   len = length(values), 
                   entry = values)
  }
  return(L)
}


### function to read a .gmt file and return a a list ###
parse.gmt <- function(fname) {
  gmt.lines <- scan(fname, what = "", sep = "\n",
                    quote = NULL, quiet = TRUE)
  tmp <- lapply(gmt.lines, function(x) unlist(strsplit(x, "\t")))
  mk.gmt.entry <- function(x) {
    L <- list()
    L[["head"]] <- x[1]
    L[["desc"]] <- x[2]
    l.entry <- x[-c(1:2)]
    idx <- l.entry != ""
    L[["entry"]] <- l.entry[idx]
    L[["len"]] <- length(L[["entry"]])
    return(L)
  }
  L <- lapply(tmp, function(x) mk.gmt.entry(x))
  names(L) <- unlist(lapply(L, function(x) x$head))
  return(L)
}


### function for writing nested list objects as gmt files
write.gmt <- function(lst, fname) {
  # assumes that each element of the list will have the fields
  # head, desc, entry
  if (file.exists(fname)) {
    message(paste(fname, "exists, deleting..."))
    file.remove(fname)
  }
  for (i in 1:length(lst)) {
    el <- lst[[i]]
    ncolumns <- 2 + length(el$entry)
    write(c(el$head, el$desc, el$entry), file=fname, sep="\t", append=T, ncolumns=ncolumns)
  }
}


########################################
### Other Misc. utility functions ######
########################################


### function to write tab-delimited text files ###
write.tbl <- function(tbl, ofile, col.names = TRUE, row.names = FALSE) {
  write.table(tbl, file = ofile, sep = "\t", quote = FALSE, 
              col.names = col.names, row.names = row.names)
}

# for backwards compatibility
mktbl <- write.tbl