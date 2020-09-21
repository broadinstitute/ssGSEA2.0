# This is a copy of the file 'io.R' downloaded from https://github.com/cmap/cmapR GitHub repository on Oct 1, 2019.

# BSD 3-Clause License
# 
# Copyright (c) 2017, Connectivity Map (CMap) at the Broad Institute, Inc. 
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   
#   * Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# * Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#          SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


########################################
### GCT class and method definitions ###
########################################

#' An S4 class to represent a GCT object
#' 
#' @slot mat a numeric matrix
#' @slot rid a character vector of row ids
#' @slot cid a character vector of column ids
#' @slot rdesc a \code{data.frame} of row descriptors
#' @slot rdesc a \code{data.frame} of column descriptors
#' @slot src a character indicating the source (usually file path) of the data
#' 
#' @description The GCT class serves to represent annotated
#'   matrices. The \code{mat} slot contains said data and the
#'   \code{rdesc} and \code{cdesc} slots contain data frames with
#'   annotations about the rows and columns, respectively
#'   
#' @seealso \code{\link{parse.gctx}}, \code{\link{write.gctx}}, \code{\link{read.gctx.meta}}, \code{\link{read.gctx.ids}}
#' @seealso \link{http://clue.io/help} for more information on the GCT format
methods::setClass("GCT",
         methods::representation(
             mat = "matrix",
             rid = "character",
             cid = "character",
             rdesc = "data.frame",
             cdesc = "data.frame",
             version = "character",
             src = "character"
         )
)


# set up methods for checking GCT validity
methods::setValidity("GCT",
  function(object) {
    # check whether dimensions of various
    # slots are in sync
    nrows <- nrow(object@mat)
    ncols <- ncol(object@mat)
    if (nrows != length(object@rid)) {
      return("rid must be the same length as number of matrix rows")
    }
    if (ncols != length(object@cid)) {
      return("cid must be the same length as number of matrix columns")
    }
    if (length(object@cid) > length(unique(object@cid))) {
      return("cid must be unique")
    }
    if (length(object@rid) > length(unique(object@rid))) {
      return("rid must be unique")
    }
    if (nrow(object@cdesc) != ncols & nrow(object@cdesc) != 0) {
      return("cdesc must either have 0 rows or the same number of rows as matrix has columns")
    }
    if (nrow(object@rdesc) != nrows & nrow(object@rdesc) != 0) {
      return("rdesc must either have 0 rows or the same number of rows as matrix has rows")
    }
    else {
      return(T)
    }
  }
)

suppressMessages({
  # set method for displaying a GCT object
  # just use the 'str' function to show its structure
  setMethod("show", methods::signature("GCT"), function(object) {
    utils::str(object)
  })

  # dim, nrow and ncol to display the # of rows and columns
  # for a GCT object's matrix
  setMethod("ncol", methods::signature("GCT"), function(x) {
    ncol(x@mat)
  })
  setMethod("nrow", methods::signature("GCT"), function(x) {
    nrow(x@mat)
  })
  setMethod("dim", methods::signature("GCT"), function(x) {
    dim(x@mat)
  })
  setMethod("range", methods::signature("GCT"), function(x, na.rm=F, finite=F) {
    range(x@mat, na.rm=na.rm, finite=finite)
  })
  setMethod("max", methods::signature("GCT"), function(x, na.rm=F) {
    max(x@mat, na.rm=na.rm)
  })
  setMethod("min", methods::signature("GCT"), function(x, na.rm=F) {
    min(x@mat, na.rm=na.rm)
  })
  setMethod("diag", methods::signature("GCT"), function(x) {
    diag(x@mat)
  })
})


#### define some helper methods for parsing gctx files ###

#' Adjust the data types for columns of a meta data frame
#' 
#' @description GCT(X) parsing initially returns data frames
#'   of row and column descriptors where all columns are of
#'   type character. This is inconvenient for analysis, so
#'   the goal of this function is to try and guess the
#'   appropriate data type for each column.
#'   
#' @param meta a data.frame
#' 
#' @details This is a low-level helper function
#'   which most users will not need to access directly
#' 
#' @return meta the same data frame with (potentially) adjusted
#'   column types.
#'   
#' @examples 
#' # meta data table with all character types
#' str(cdesc_char)
#' fixed <- cmapR:::fix.datatypes(cdesc_char)
#' # note how some column classes have changed
#' str(fixed)
#' 
#' @family GCTX parsing functions
#' @keywords internal
fix.datatypes <- function(meta) {
    for (field.name in names(meta)) {
        # get the field values
        field <- meta[[field.name]]
        # check if it's numeric. data may come in as a string
        # but actually contains numeric values. if so, as.numeric
        # will not result in a vector of NA values
        field.as.numeric <- suppressWarnings(as.numeric(field))
        if (!any(is.na(field.as.numeric))) {
          field <- field.as.numeric
        }
        if (is.numeric(field)) {
            # check if it's an integer. data may be floats but
            # if we coerce to an integer and the difference from
            # original values is zero, that means data are actually
            # integers. integer conversion will return NA if there
            # are any issues.
            field.as.integer <- suppressWarnings(as.integer(field))
            if (!any(is.na(field.as.integer))) {
              # integer conversion was fine, lets see if the
              # values are altered
              diffs <- field - field.as.integer
              if (all(diffs == 0)) {
                # converting to integer didn't change value,
                # set field to integer values
                field <- field.as.integer
              }
            }
        }
        # insert back into the annotations
        meta[[field.name]] <- field
    }
    return(meta)
}


#' Parse row or column metadata from GCTX files
#' 
#' @param gctx_path the path to the GCTX file
#' @param dimension which metadata to read (row or column)
#' @param ids a character vector of a subset of row/column ids
#'   for which to read the metadata
#'  
#' @return a \code{data.frame} of metadata
#' 
#' @examples
#' gct_file <- system.file("extdata", "modzs_n272x978.gctx", package="cmapR") 
#' # row meta
#' row_meta <- read.gctx.meta(gct_file)
#' str(row_meta)
#' # column meta
#' col_meta <- read.gctx.meta(gct_file, dimension="column")
#' str(col_meta)
#' # now for only the first 10 ids
#' col_meta_first10 <- read.gctx.meta(gct_file, dimension="column", ids=col_meta$id[1:10])
#' str(col_meta_first10)
#' 
#' @family GCTX parsing functions
#' @export
read.gctx.meta <- function(gctx_path, dimension="row", ids=NULL) {
  if (!file.exists(gctx_path)) {
    stop(paste(gctx_path, "does not exist"))
  }
  if (dimension=="column") dimension <- "col"
  if (!(dimension %in% c("row", "col"))) {
    stop("dimension can be either row or col")
  }
  if (dimension == "row") {
    name <- "0/META/ROW"
  } else {
    name <- "0/META/COL"
  }
  raw_annots <- rhdf5::h5read(gctx_path, name=name) # returns a list
  fields <- names(raw_annots)
  # define an empty data frame of the correct dimensions
  annots <-  data.frame(matrix(nrow=length(raw_annots[[fields[1]]]), ncol=length(fields)))
  names(annots) <-  fields
  # loop through each field and fill the annots data.frame
  for (i in 1:length(fields)) {
    field <- fields[i]
    # remove any trailing spaces
    # and cast as vector
    annots[,i] <- as.vector(gsub("\\s*$", "", raw_annots[[field]], perl=T))
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
  annots$id <- as.character(annots$id)
  return(annots)
}


#' Read GCTX row or column ids
#' 
#' @param gctx_path path to the GCTX file
#' @param dimension which ids to read (row or column)
#' 
#' @return a character vector of row or column ids from the provided file
#' 
#' @examples 
#' gct_file <- system.file("extdata", "modzs_n272x978.gctx", package="cmapR")
#' # row ids
#' rid <- read.gctx.ids(gct_file)
#' head(rid)
#' # column ids
#' cid <- read.gctx.ids(gct_file, dimension="column")
#' head(cid)
#' 
#' @family GCTX parsing functions
#' @export
read.gctx.ids <- function(gctx_path, dimension="row") {
  if (!file.exists(gctx_path)) {
    stop(paste(gctx_path, "does not exist"))
  }
  if (dimension=="column") dimension <- "col"
  if (!(dimension %in% c("row", "col"))) {
    stop("dimension can be either row or col")
  }
  if (dimension == "row") {
    name <- "0/META/ROW/id"
  } else {
    name <- "0/META/COL/id"
  }
  # remove any spaces
  ids <- gsub("\\s*$", "", rhdf5::h5read(gctx_path, name=name), perl=T)
  # cast as character
  ids <- as.character(ids)
  return(ids)
}

#' Return a subset of requested GCTX row/colum ids
#' out of the universe of all ids
#' 
#' @details This is a low-level helper function
#'   which most users will not need to access directly
#' 
#' @param ids vector of requested ids. If \code{NULL}, no
#'   subsetting is performed
#' @param all_ids vector of universe of ids
#' @param type flag indicating the type of ids being processed
#' 
#' @return a list with the following elements
#'  \code{ids}: a character vector of the processed ids
#'  \code{idx}: an integer list of their corresponding indices in \code{all_ids}
#' 
#' @examples 
#' gct_file <- system.file("extdata", "modzs_n272x978.gctx", package="cmapR")
#' ids <- read.gctx.ids(gct_file)
#' processed_ids <- cmapR:::process_ids(ids[1:10], ids)
#' str(processed_ids)
#' 
#' @family GCTX parsing functions
#' @keywords internal
process_ids <- function(ids, all_ids, type="rid") {
  if (!is.null(ids)) {
    if (is.numeric(ids)) {
      # is it numeric?
      idx <- ids
      is_invalid_idx <- (idx > length(all_ids)) | (idx <= 0)
      invalid_idx <- idx[is_invalid_idx]
      if (all(is_invalid_idx)) {
        stop(paste("none of the requested", type, "indices were found in the dataset"))
      }
      if (any(is_invalid_idx)) {
        # requested indices are outside of the possible range
        warning(paste("the following ", type, " were are outside possible range and will be ignored:\n",
                      paste(invalid_idx, collapse="\n"), sep="")) 
      }
      idx <- idx[!is_invalid_idx]
    } else {
      # assume its a character
      idx <- match(ids, all_ids)
      if (all(is.na(idx))) {
        stop(paste("none of the requested", type, "were found in the dataset"))
      }
      if (any(is.na(idx))) {
        ids_not_found <- ids[is.na(idx)]
        warning(paste("the following ", type, " were not found and will be ignored:\n",
                      paste(ids_not_found, collapse="\n"), sep=""))
      }
      idx <- idx[!is.na(idx)]
    }
  } else {
    # ids were null, just return an index vector
    # allong all_ids
    idx <- seq_along(all_ids)
  }
  # subset the character ids to the ones we want
  id_keep <- as.character(all_ids[idx])
  return(list(idx=idx, ids=id_keep))
}


# define the initialization method for the GCT class
methods::setMethod("initialize",
          signature = "GCT",
          definition = function(.Object, mat=NULL, rdesc=NULL, cdesc=NULL, src=NULL, rid=NULL, cid=NULL,
                                matrix_only=FALSE) {
            # if we were supplied a matrix and annotations, use them
            if (!is.null(mat)) {
              .Object@mat <- mat
              # if given rid and cid, use those as well
              if (!is.null(rid)) {
                .Object@rid <- rid
              } else {
                .Object@rid <- rownames(mat)
              }
              if (!is.null(cid)) {
                .Object@cid <- cid
              } else {
                .Object@cid <- colnames(mat)
              }
            }
            if (!is.null(rdesc)) {
              .Object@rdesc <- rdesc
            }
            if (!is.null(cdesc)) {
              .Object@cdesc <- cdesc
            } else if (!is.null(src)) {
              # we were not given a matrix, were we given a src file?
              # check to make sure it's either .gct or .gctx
              if (! (grepl(".gct$", src) || grepl(".gctx$", src) ))
                stop("Either a .gct or .gctx file must be given")
              if (grepl(".gct$", src)) {
                if ( ! is.null(rid) || !is.null(cid) )
                  warning(paste("rid and cid values may only be given for .gctx files, not .gct files\n",
                                "ignoring"))
                # parse the .gct
                .Object@src <- src
                # get the .gct version by reading first line
                .Object@version <- scan(src, what = "", nlines = 1, sep = "\t", quiet = TRUE)[1]
                # get matrix dimensions by reading second line
                dimensions <- scan(src, what = double(0), nlines = 1, skip = 1, sep = "\t", quiet = TRUE)
                nrmat <- dimensions[1]
                ncmat <- dimensions[2]
                if (length(dimensions)==4) {
                  # a #1.3 file
                  message("parsing as GCT v1.3")
                  nrhd <- dimensions[3]
                  nchd <- dimensions[4]
                } else {
                  # a #1.2 file
                  message("parsing as GCT v1.2")
                  nrhd <- 0
                  nchd <- 0
                }
                message(paste(src, nrmat, "rows,", ncmat, "cols,", nrhd, "row descriptors,", nchd, "col descriptors"))
                # read in header line
                header <- scan(src, what = "", nlines = 1, skip = 2, sep = "\t", quote = NULL, quiet = TRUE)
                # construct row header and column id's from the header line
                if ( nrhd > 0 ) {
                  rhd <- header[2:(nrhd+1)]
                  cid <- header[-(nrhd+1):-1]
                  col_offset <- 1
                }
                else {
                  if (any(grepl("description", header, ignore.case=T))) {
                    # check for presence of description column in v1.2 files
                    col_offset <- 2
                  } else {
                    col_offset <- 1
                  }
                  rhd <- NULL
                  cid <- header[(1+col_offset):length(header)]
                }
                # read in the next set of headers (column annotations) and shape into a matrix
                if ( nchd > 0 ) {
                  header <- scan(src, what = "", nlines = nchd, skip = 3, sep = "\t", 
                                quote = NULL, quiet = TRUE)		
                  header <- matrix(header, nrow = nchd, 
                                  ncol = ncmat + nrhd + 1, byrow = TRUE)
                  # extract the column header and column descriptions
                  chd <- header[,1]
                  cdesc <- header[,-(nrhd+1):-1]
                  # need to transpose in the case where there's only one column annotation
                  if ( nchd == 1 )
                    cdesc <- t(cdesc)
                }
                else {
                  chd = NULL
                  cdesc <- data.frame(id=cid)
                }
                # read in the data matrix and row descriptions, shape into a matrix
                mat <- scan(src, what = "", nlines = nrmat, 
                           skip = 3 + nchd, sep = "\t", quote = NULL, quiet = TRUE)
                mat <- matrix(mat, nrow = nrmat, ncol = ncmat + nrhd + col_offset, 
                             byrow = TRUE)
                # message(paste(dim(mat), collapse="\t"))
                # Extract the row id's row descriptions, and the data matrix
                rid <- mat[,1]
                if ( nrhd > 0 ) {
                  # need as.matrix for the case where there's only one row annotation
                  rdesc <- as.matrix(mat[,2:(nrhd + 1)])
                  mat <- matrix(as.numeric(mat[,-(nrhd + 1):-1]),
                               nrow = nrmat, ncol = ncmat)
                }
                else {
                  rdesc <- data.frame(id=rid)
                  mat <- matrix(as.numeric(mat[, (1+col_offset):ncol(mat)]), nrow = nrmat, ncol = ncmat)
                }
                # assign names to the data matrix and the row and column descriptions
                # message(paste(dim(mat), collapse="\t"))
                dimnames(mat) <- list(rid, cid)
                if ( nrhd > 0 ) {
                  dimnames(rdesc) <- list(rid, rhd)
                  rdesc <- as.data.frame(rdesc, stringsAsFactors = FALSE)
                }
                if ( nchd > 0 ) {
                  cdesc <- t(cdesc)
                  dimnames(cdesc) <- list(cid,chd)
                  cdesc <- as.data.frame(cdesc, stringsAsFactors = FALSE)
                }
                # assign to the GCT slots
                .Object@mat <- mat
                .Object@rid <- rownames(mat)
                .Object@cid <- colnames(mat)
                if (!matrix_only) {
                  # return annotations as well as matrix
                  .Object@rdesc <- fix.datatypes(rdesc)
                  .Object@cdesc <- fix.datatypes(cdesc)
                  # add id columns to rdesc and cdesc
                  .Object@rdesc$id <- rownames(.Object@rdesc)
                  .Object@cdesc$id <- rownames(.Object@cdesc)
                }
              }
              else { 
                # parse the .gctx
                message(paste("reading", src))
                .Object@src <- src
                # if the rid's or column id's are .grp files, read them in
                if ( length(rid) == 1 && grepl(".grp$", rid) )
                  rid <- parse.grp(rid)
                if ( length(cid) == 1 && grepl(".grp$", cid) )
                  cid <- parse.grp(cid)
                # get all the row and column ids
                all_rid <- read.gctx.ids(src, dimension="row")
                all_cid <- read.gctx.ids(src, dimension="col")
                # if rid or cid specified, read only those rows/columns
                # if already numeric, use as is
                # else convert to numeric indices
                processed_rids <- process_ids(rid, all_rid, type="rid")
                processed_cids <- process_ids(cid, all_cid, type="cid")
                # read the data matrix
                .Object@mat <- rhdf5::h5read(src, name="0/DATA/0/matrix",
                                      index=list(processed_rids$idx, processed_cids$idx))
                # set the row and column ids, casting as characters
                .Object@rid <- processed_rids$ids
                .Object@cid <- processed_cids$ids
                rownames(.Object@mat) <- processed_rids$ids
                colnames(.Object@mat) <- processed_cids$ids
                # get the meta data
                if (!matrix_only) {
                  .Object@rdesc <- read.gctx.meta(src, dimension="row", ids=processed_rids$ids)
                  .Object@cdesc <- read.gctx.meta(src, dimension="col", ids=processed_cids$ids)
                }
                else {
                  .Object@rdesc <- data.frame(id=.Object@rid, stringsAsFactors = F)
                  .Object@cdesc <- data.frame(id=.Object@cid, stringsAsFactors = F)
                }
                # close any open handles and return the object
                if(utils::packageVersion('rhdf5') < "2.23.0") {
                    rhdf5::H5close()
                } else {
                  rhdf5::h5closeAll()
                }
                message("done")
              }
            }
            # finally, make sure object is valid before returning
            ok <- methods::validObject(.Object)
            return(.Object)
          }
)


#' Parse a GCTX file into the workspace as a GCT object
#' 
#' @param fname path to the GCTX file on disk
#' @param rid either a vector of character or integer
#'   row indices or a path to a grp file containing character
#'   row indices. Only these indicies will be parsed from the
#'   file.
#' @param cid either a vector of character or integer
#'   column indices or a path to a grp file containing character
#'   column indices. Only these indicies will be parsed from the
#'   file.
#' @param matrix_only boolean indicating whether to parse only
#'   the matrix (ignoring row and column annotations)
#'
#' @details \code{parse.gctx} also supports parsing of plain text
#'   GCT files, so this function can be used as a general GCT parser.
#' 
#' @examples 
#' gct_file <- system.file("extdata", "modzs_n272x978.gctx", package="cmapR")
#' (ds <- parse.gctx(gct_file))
#' 
#' # matrix only
#' (ds <- parse.gctx(gct_file, matrix_only=TRUE))
#' 
#' # only the first 10 rows and columns
#' (ds <- parse.gctx(gct_file, rid=1:10, cid=1:10))
#' 
#' @family GCTX parsing functions
#' @export
parse.gctx <- function(fname, rid=NULL, cid=NULL, matrix_only=FALSE) {
    ds <- methods::new("GCT",
              src = fname,
              rid = rid,
              cid = cid,
              matrix_only = matrix_only)
    return(ds)
}

#' Append matrix dimensions to filename
#' 
#' @param ofile the file name
#' @param mat the matrix
#' @param extension the file extension
#' 
#' @return a character string of the filename with
#'   matrix dimensions appended
#' 
#' @details This is a helper function that most users
#'   will not use directly
#'   
#' @examples 
#' (filename <- cmapR:::append.dim("my.gctx.filename", matrix(nrow=10, ncol=15)))
#'   
#' @keywords internal
#' @family GCTX parsing functions
append.dim <- function(ofile, mat, extension="gct") {
  nc <- ncol(mat)
  nr <- nrow(mat)
  filename <- basename(ofile)
  if (grepl("n[0-9]+x[0-9]+\\.gct", filename)) {
    # already has a dimensions token, ignore
    filename <- sub("_n[0-9]+x[0-9]+\\.gct.*", "", filename)
  } 
  filename <- file.path(dirname(ofile),
                    sprintf('%s_n%dx%d.%s',filename,
                            nc, nr, extension))
  return(filename)
}


#' Write a GCT object to disk in GCT format
#' 
#' @param ds the GCT object
#' @param ofile the desired output filename
#' @param precision the numeric precision at which to
#'   save the matrix. See \code{details}.
#' @param appenddim boolean indicating whether to append
#'   matrix dimensions to filename
#' @param ver the GCT version to write. See \code{details}.
#' 
#' @details Since GCT is text format, the higher \code{precision}
#'   you choose, the larger the file size.
#'   \code{ver} is assumed to be 3, aka GCT version 1.3, which supports
#'   embedded row and column metadata in the GCT file. Any other value
#'   passed to \code{ver} will result in a GCT version 1.2 file which
#'   contains only the matrix data and no annotations.
#'
#' @return NULL
#' 
#' @examples 
#' \dontrun{
#' write.gct(ds, "dataset", precision=2)
#' }
#' @family GCTX parsing functions
#' @export
write.gct <- function(ds, ofile, precision=4, appenddim=T, ver=3) {
  if (!class(ds)=="GCT") {
    stop("ds must be a GCT object")
  }
  # make sure it's valid
  ok <- methods::validObject(ds)
  # append the dimensions of the data set, if desired
  if (appenddim) ofile <- append.dim(ofile, ds@mat, extension="gct")
  
  precision = floor(precision)
  cat(sprintf('Saving file to %s\n',ofile))
  nr <- nrow(ds@mat)
  nc <- ncol(ds@mat)
  cat(sprintf('Dimensions of matrix: [%dx%d]\n',nr,nc))
  cat(sprintf('Setting precision to %d\n',precision))
  
  # open file and write   
  if (ver==3) {
    # remove the 'id' columns
    ds@cdesc$id <- NULL
    ds@rdesc$id <- NULL
    # get the counts of meta data fields
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
    if (ncdesc > 0) {
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
                  ds@rdesc[ii, 2],
                  round(ds@mat[ii,],precision)),collapse='\t'),
          sep='\n',file=ofile,append=T)
    }
  }

  cat(sprintf('Saved.\n'))  
}


#' Write a GCT object to disk in GCTX format
#' 
#' @param ds a GCT object
#' @param ofile the desired file path for writing
#' @param appenddim boolean indicating whether the
#'   resulting filename will have dimensions appended
#'   (e.g. my_file_n384x978.gctx)
#' @param compression_level integer between 1-9 indicating
#'   how much to compress data before writing. Higher values
#'   result in smaller files but slower read times.
#' @param matrix_only boolean indicating whether to write
#'   only the matrix data (and skip row, column annotations)
#' @param max_chunk_kb for chunking, the maximum number of KB
#'   a given chunk will occupy
#' 
#' @examples 
#' \dontrun{
#' # assume ds is a GCT object
#' write.gctx(ds, "my/desired/outpath/and/filename")
#' }
#' @family GCTX parsing functions
#' @export
write.gctx <- function(ds, ofile, appenddim=T, compression_level=0, matrix_only=F,
                       max_chunk_kb=1024) {
  if (!class(ds)=="GCT") {
    stop("ds must be a GCT object")
  }
  # make sure it's valid
  ok <- methods::validObject(ds)
  # add dimensions to filename if desired
  if (appenddim) ofile <- append.dim(ofile, ds@mat, extension="gctx")
  # check if the file already exists
  if (file.exists(ofile)) {
    message(paste(ofile, "exists, removing"))
    file.remove(ofile)
  }
  message(paste("writing", ofile))
  
  # start the file object
  rhdf5::h5createFile(ofile)
  
  # create all the necessary groups
  rhdf5::h5createGroup(ofile, "0")
  rhdf5::h5createGroup(ofile, "0/DATA")
  rhdf5::h5createGroup(ofile, "0/DATA/0")
  rhdf5::h5createGroup(ofile, "0/META")
  rhdf5::h5createGroup(ofile, "0/META/COL")
  rhdf5::h5createGroup(ofile, "0/META/ROW")
  
  # create and write matrix data, using chunking
  bits_per_element <- switch(storage.mode(ds@mat),
                             "double" = 64,
                             "integer" = 32)
  elem_per_kb <- max_chunk_kb * 8 / bits_per_element
  # assume matrix is of dimensions row_dim x col_dim 
  row_dim <- nrow(ds)
  col_dim <- ncol(ds)
  row_chunk_size <- min(row_dim, 1000)
  # column chunk, such that row * col <= max_chunk_kb
  col_chunk_size <- min(((max_chunk_kb * elem_per_kb) %/% row_chunk_size), col_dim)
  chunking <- c(row_chunk_size, col_chunk_size) 
  message(paste(c("chunk sizes:", chunking), collapse="\t"))
  rhdf5::h5createDataset(ofile, "0/DATA/0/matrix", dim(ds@mat), chunk=chunking, level=compression_level)
  rhdf5::h5write.default(ds@mat, ofile, "0/DATA/0/matrix")
  
  # write annotations
  rhdf5::h5write.default(as.character(ds@rid), ofile, "0/META/ROW/id")
  rhdf5::h5write.default(as.character(ds@cid), ofile, "0/META/COL/id")
  
  if (!matrix_only) {
    write.gctx.meta(ofile, ds@cdesc, dimension="column")
    write.gctx.meta(ofile, ds@rdesc, dimension="row")
  }

  # close any open handles
  if(utils::packageVersion('rhdf5') < "2.23.0") {
    rhdf5::H5close()
  } else {
    rhdf5::h5closeAll()
  }

  # add the version annotation and close
  fid <- rhdf5::H5Fopen(ofile)
  rhdf5::h5writeAttribute.character("GCTX1.0", fid, "version")
  rhdf5::H5Fclose(fid)

}

#' Update the matrix of an existing GCTX file
#' 
#' @param x an array of data
#' @param ofile the filename of the GCTX to update
#' @param rid integer indices or character ids of the rows
#'   to update
#' @param cid integer indices or character ids of the columns
#'   to update
#' 
#' @details Overwrite the rows and columns of \code{ofile} 
#' as indicated by \code{rid} and \code{cid} respectively.
#' \code{rid} and \code{cid} can either be integer indices
#' or character ids corresponding to the row and column ids
#' in \code{ofile}.
#' 
#' @examples
#' \dontrun{
#' m <- matrix(rnorm(20), nrow=10)
#' # update by integer indices
#' update.gctx(m, ofile="my.gctx", rid=1:10, cid=1:2)
#' # update by character ids
#' row_ids <- letters[1:10]
#' col_ids <- LETTERS[1:2]
#' update.gctx(m, ofile="my.gctx", rid=row_ids, cid=col_ids)
#' }
#' @export
update.gctx <- function(x, ofile, rid=NULL, cid=NULL) {
  # x must be numeric
  stopifnot(is.numeric(x))
  # must give us at least one of rid or cid
  if (is.null(rid) && is.null(cid)) {
    stop("one of rid or cid must not be NULL")
  }
  # make sure the dimensions of x agree with the 
  # number of rid and cid supplied
  if (is.matrix(x)) {
    stopifnot(all(dim(x) == c(length(rid), length(cid))))
  } else {
    # x is a vector, so we must be updating in one dimension
    if(!is.null(rid) & !is.null(cid)) {
      stop(paste("x is a vector so you can only update in one dimension",
                 "(only one of rid or cid can be non-NULL)", sep="\n"))
    }
    if (is.null(rid)) {
      stopifnot(length(cid) == length(x))
    }
    if (is.null(cid)) {
      stopifnot(length(rid) == length(x))
    }
  }
  # get the dimensions of the file in question
  info <- rhdf5::h5dump(ofile, load=F)
  dims <- as.numeric(unlist(strsplit(info[["0"]][["DATA"]][["0"]][["matrix"]][["dim"]], " x ")))
  # get the full space of row and column ids
  all_rid <- cmapR::read.gctx.ids(ofile, dim="row")
  all_cid <- cmapR::read.gctx.ids(ofile, dim="col")
  # helper functions to validate integer and character ids
  validate_integer_ids <- function(ids, maxdim, which_dim) {
    stopifnot(all(ids > 0))
    out_of_range <- setdiff(ids, seq_len(maxdim))
    if (length(out_of_range) > 0) {
      stop(paste("the following", which_dim, "indices are out of range\n",
                 paste(out_of_range, collapse="\n")))
    }
  }
  validate_character_ids <- function(ids, all_ids, which_dim) {
    out_of_range <- setdiff(ids, all_ids)
    if (length(out_of_range) > 0) {
      stop(paste("the following", which_dim, "ids do not exist in the dataset\n",
                 paste(out_of_range, collapse="\n")))
    }
  }
  # given integer ids
  if (is.integer(rid)) {
    validate_integer_ids(rid, dims[1], "row")
    ridx <- rid
  }
  if (is.integer(cid)) {
    validate_integer_ids(cid, dims[2], "column")
    cidx <- cid
  }
  # given character ids
  if (is.character(rid)) {
    validate_character_ids(rid, all_rid, "row")
    ridx <- match(rid, all_rid)
  }
  if (is.character(cid)) {
    validate_character_ids(cid, all_cid, "column")
    cidx <- match(cid, all_cid)
  }
  # make the updates to the specified rows/columns
  rhdf5::h5write.default(x, ofile, "0/DATA/0/matrix", index=list(ridx, cidx))
  # close any open handles
  if(utils::packageVersion('rhdf5') < "2.23.0") {
    rhdf5::H5close()
  } else {
    rhdf5::h5closeAll()
  }
}

#' Write a \code{data.frame} of meta data to GCTX file
#' 
#' @param ofile the desired file path for writing
#' @param df the \code{data.frame} of annotations
#' @param dimension the dimension to annotate
#'   (row or column)
#' 
#' @examples 
#' \dontrun{
#' # assume ds is a GCT object
#' cmapR:::write.gctx.meta("/my/file/path", cdesc_char, dimension="col")
#' }
#' @family GCTX parsing functions
#' @keywords internal
write.gctx.meta <- function(ofile, df, dimension="row") {
  path <- if ((dimension=="row")) "0/META/ROW/" else "0/META/COL/"
  # loop through all columns
  fields <- names(df)
  if (length(fields) > 0) {
    for (i in 1:length(fields)) {
      field <- fields[i]
      # if this is the id field, skip b/c that field is special
      # and is written as part of write.gctx
      if (field == "id") next
      v <- df[, i]
      # convert factors to character
      if(class(v) == "factor" || class(v) == "AsIs") {
        v <- as.character(v)
      }
      rhdf5::h5write.default(v, ofile, paste(path, field, sep=""))
    }
  }
}


###########################################
### functions for other CMap file types ###
###########################################

#' Read a GRP file and return a vector of its contents
#' @param fname the file path to be parsed
#' @return a vector of the contents of \code{fname}
#' @examples 
#' grp_path <- system.file("extdata", "lm_epsilon_n978.grp", package="cmapR")
#' values <- parse.grp(grp_path)
#' str(values)
#' @family CMap parsing functions
#' @seealso \link{http://clue.io/help} for details on the GRP file format
#' @export
parse.grp <- function(fname) {
  grp <- scan(fname, what = "", quote = NULL, quiet = TRUE, sep="\n")
  return(grp)
}


#' Write a vector to a GRP file
#' 
#' @param vals the vector of values to be written
#' @param fname the desired file name
#' 
#' @examples 
#' \dontrun{
#' write.grp(letters, "letter.grp")
#' }
#' 
#' @family CMap parsing functions
#' @seealso \link{http://clue.io/help} for details on the GRP file format
#' @export
write.grp <- function(vals, fname) {
  if (is.list(vals)) vals <- unlist(vals)
  if (!is.vector(vals)) vals <- as.vector(vals)
  write(vals, fname, ncolumns=1)
}


#' Read a GMX file and return a list
#' 
#' @param fname the file path to be parsed
#' 
#' @return a list of the contents of \code{fname}. See details.
#' 
#' @details \code{parse.gmx} returns a nested list object. The top
#'   level contains one list per column in \code{fname}. Each of 
#'   these is itself a list with the following fields:
#'   - \code{head}: the name of the data (column in \code{fname})
#'   - \code{desc}: description of the corresponding data
#'   - \code{len}: the number of data items
#'   - \code{entry}: a vector of the data items
#' 
#' @examples 
#' gmx_path <- system.file("extdata", "lm_probes.gmx", package="cmapR")
#' gmx <- parse.gmx(gmx_path)
#' str(gmx)
#' 
#' @family CMap parsing functions
#' @seealso \link{http://clue.io/help} for details on the GMX file format
#' @export
parse.gmx <- function(fname) {
    tmp <- utils::read.table(fname, sep = "\t", 
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


#' Read a GMT file and return a list
#' @param fname the file path to be parsed
#' 
#' @return a list of the contents of \code{fname}. See details.
#' 
#' @details \code{parse.gmt} returns a nested list object. The top
#'   level contains one list per row in \code{fname}. Each of 
#'   these is itself a list with the following fields:
#'   - \code{head}: the name of the data (row in \code{fname})
#'   - \code{desc}: description of the corresponding data
#'   - \code{len}: the number of data items
#'   - \code{entry}: a vector of the data items
#' 
#' @examples 
#' gmt_path <- system.file("extdata", "query_up.gmt", package="cmapR")
#' gmt <- parse.gmt(gmt_path)
#' str(gmt)
#' 
#' @family CMap parsing functions
#' @seealso \link{http://clue.io/help} for details on the GMT file format
#' @export
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


#' Write a nested list to a GMT file
#' 
#' @param lst the nested list to write. See \code{details}.
#' @param fname the desired file name
#' 
#' @details \code{lst} needs to be a nested list where each 
#'   sub-list is itself a list with the following fields:
#'   - \code{head}: the name of the data
#'   - \code{desc}: description of the corresponding data
#'   - \code{len}: the number of data items
#'   - \code{entry}: a vector of the data items
#' 
#' @examples 
#' \dontrun{
#' write.gmt(gene_set, "gene_set.gmt")
#' }
#' 
#' @family CMap parsing functions
#' @seealso \link{http://clue.io/help} for details on the GMT file format
#' @export
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

#' Read an LXB file and return a matrix
#' 
#' @param lxb_path the path to the lxb file
#' @param columns which columns in the lxb file to retain
#' @param newnames what to name these columns in the returned matrix
#' 
#' @examples 
#' lxb_path <- system.file("extdata", "example.lxb", package="cmapR")
#' lxb_data <- lxb2mat(lxb_path)
#' str(lxb_data)
#' 
#' @importFrom prada readFCS exprs
#' 
#' @family CMap parsing functions
#' @export
lxb2mat <- function(lxb_path, columns=c("RID", "RP1"),
                    newnames=c("barcode_id", "FI")) {
  message(paste("reading", lxb_path))
  # suppressing warning about signed integers since
  # lxb data will be unsigned
  lxb <- suppressWarnings(prada::readFCS(lxb_path))
  m <- prada::exprs(lxb)[, columns]
  keep_idx <- m[, 1] != 0
  m <- m[keep_idx, ]
  colnames(m) <- newnames
  return(m)
}

########################################
### Other Misc. utility functions ######
########################################


#' Write a \code{data.frame} to a tab-delimited text file
#' 
#' @param tbl the \code{data.frame} to be written
#' @param ofile the desired file name
#' @param ... additional arguments passed on to \code{write.table}
#' 
#' @details This method simply calls \code{write.table} with some
#'   preset arguments that generate a unquoated, tab-delimited file
#'   without row names.
#' 
#' @examples 
#' \dontrun{
#' write.tbl(cdesc_char, "col_meta.txt")
#' }
#' 
#' @seealso \code{\link{write.table}}
#' @export
write.tbl <- function(tbl, ofile, ...) {
    utils::write.table(tbl, file = ofile, sep="\t", quote=F,
      col.names=T, row.names=F, ...)
}
