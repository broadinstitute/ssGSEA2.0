# This is a copy of the file 'utils.R' downloaded from https://github.com/cmap/cmapR GitHub repository on Oct 1, 2019.

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


#' Transform a GCT object in to a long form \code{\link{data.table}} (aka 'melt')
#' 
#' @description Utilizes the \code{\link{data.table::melt}} function to transform the
#'   matrix into long form. Optionally can include the row and column
#'   annotations in the transformed \code{\link{data.table}}.
#'   
#' @param g the GCT object
#' @param keep_rdesc boolean indicating whether to keep the row
#'   descriptors in the final result
#' @param keep_cdesc boolean indicating whether to keep the column
#'   descriptors in the final result
#' @param remove_symmetries boolean indicating whether to remove
#'   the lower triangle of the matrix (only applies if \code{g@mat} is symmetric)
#' @param suffixes the character suffixes to be applied if there are
#'   collisions between the names of the row and column descriptors
#' @param ... further arguments passed along to \code{data.table::merge}
#'   
#' @return a \code{\link{data.table}} object with the row and column ids and the matrix
#'   values and (optinally) the row and column descriptors
#'   
#' @examples 
#' # simple melt, keeping both row and column meta
#' head(melt.gct(ds))
#' 
#' # update row/colum suffixes to indicate rows are genes, columns experiments
#' head(melt.gct(ds, suffixes = c("_gene", "_experiment")))
#' 
#' # ignore row/column meta
#' head(melt.gct(ds, keep_rdesc = FALSE, keep_cdesc = FALSE))
#' 
#' @family GCT utilities
#' @export
setGeneric("melt.gct", function(g, suffixes=NULL, remove_symmetries=F,
                                keep_rdesc=T, keep_cdesc=T, ...) {
  standardGeneric("melt.gct")
})
setMethod("melt.gct", signature("GCT"),
          function(g, suffixes, remove_symmetries=F, keep_rdesc=T, keep_cdesc=T, ...) {
          # melt a gct object's matrix into a data.frame and merge row and column
          # annotations back in, using the provided suffixes
          # assumes rdesc and cdesc data.frames both have an 'id' field.
          # merges row and/or column annotations into the melted matrix as indicated by
          # keep_rdesc and keep_cdesc, respectively.
          # if remove_symmetries, will check whether matrix is symmetric
          # and return only values corresponding to the upper triangle
          # check whether rdesc/cdesc are empty
          # if so, fill with id column
          if (nrow(g@rdesc) == 0) g@rdesc <- data.frame(id=g@rid)
          if (nrow(g@cdesc) == 0) g@cdesc <- data.frame(id=g@cid)
          # first, check if matrix is symmetric
          # if it is, use only the upper triangle
          message("melting GCT object...")
          mat <- g@mat
          if (remove_symmetries & isSymmetric(mat)) {
            mat[upper.tri(mat, diag=F)] <- NA
          }
          mat <- data.table::data.table(mat)
          mat$rid <- g@rid
          d <- data.table::melt(mat, id.vars="rid")
          data.table::setattr(d, "names", c("id.x", "id.y", "value"))
          d$id.x <- as.character(d$id.x)
          d$id.y <- as.character(d$id.y)
          # standard data.frame subset here to comply with testthat
          d <- subset(d, !is.na(value))
          if (keep_rdesc && keep_cdesc) {
            # merge back in both row and column descriptors
            data.table::setattr(d, "names", c("id", "id.y", "value"))
            d <- merge(d, data.table::data.table(g@rdesc), by="id", all.x=T, ...)
            data.table::setnames(d, "id", "id.x")
            data.table::setnames(d, "id.y", "id")
            d <- merge(d, data.table::data.table(g@cdesc), by="id", all.x=T, ...)
            data.table::setnames(d, "id", "id.y")
          } else if (keep_rdesc) {
            # keep only row descriptors
            rdesc <- data.table::data.table(g@rdesc)
            data.table::setnames(rdesc, "id", "id.x")
            d <- merge(d, rdesc, by="id.x", all.x=T, ...)
          } else if (keep_cdesc) {
            # keep only column descriptors
            cdesc <- data.table::data.table(g@cdesc)
            data.table::setnames(cdesc, "id", "id.y")
            d <- merge(d, cdesc, by="id.y", all.x=T, ...)
          }
          # use suffixes if provided
          if (!is.null(suffixes) & length(suffixes) == 2) {
            newnames <- gsub("\\.x", suffixes[1], names(d))
            newnames <- gsub("\\.y", suffixes[2], newnames)
            data.table::setattr(d, "names", newnames)
          }
          message("done")
          return(d)
})


#' Check if x is a whole number
#'
#' @param x number to test
#' @param tol the allowed tolerance
#' @return boolean indicating whether x is tol away from a whole number value
#' @examples
#' is.wholenumber(1)
#' is.wholenumber(0.5)
#' @export
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {
  return(abs(x - round(x)) < tol)
}

#' Check whether \code{test_names} are columns in the \code{\link{data.frame}} df
#' @param test_names a vector of column names to test
#' @param df the \code{\link{data.frame}} to test against
#' @param throw_error boolean indicating whether to throw an error if
#'   any \code{test_names} are not found in \code{df}
#' @return boolean indicating whether or not all \code{test_names} are
#'   columns of \code{df}
#' @examples 
#' check_colnames(c("pert_id", "pert_iname"), cdesc_char)            # TRUE
#' check_colnames(c("pert_id", "foobar"), cdesc_char, throw_error=FALSE) # FALSE, suppress error
#' @export
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

#' Do a robust \code{\link{data.frame}} subset to a set of ids
#' @param df \code{\link{data.frame}} to subset
#' @param ids the ids to subset to
#' @return a subset version of \code{df}
#' @keywords internal
subset_to_ids <- function(df, ids) {
  # helper function to do a robust df subset
  check_colnames("id", df)
  newdf <- data.frame(df[match(ids, df$id), ])
  names(newdf) <- names(df)
  return(newdf)
}


#' Subset a gct object using the provided row and column ids
#'
#' @param g a gct object
#' @param rid a vector of character ids or integer indices for ROWS
#' @param cid a vector of character ids or integer indices for COLUMNS
#' @examples
#' # first 10 rows and columns by index
#' (a <- subset.gct(ds, rid=1:10, cid=1:10))
#' 
#' # first 10 rows and columns using character ids
#' (b <- subset.gct(ds, rid=ds@rid[1:10], cid=ds@cid[1:10]))
#' 
#' identical(a, b) # TRUE
#' 
#' @family GCT utilities
#' @export
setGeneric("subset.gct", function(g, rid=NULL, cid=NULL) {
  standardGeneric("subset.gct")
})
setMethod("subset.gct", signature("GCT"),
          function(g, rid, cid) {
          # ids can either be a vector of character strings corresponding
          # to row / column ids in the gct object, or integer vectors
          # corresponding to row / column indices
          if (is.null(rid)) rid <- g@rid
          if (is.null(cid)) cid <- g@cid
          # see whether we were given characters or integers
          # and handle accordingly
          process_ids <- function(ids, ref_ids, param) {
            # simple helper function to handle id/idx conversion
            # for character or integer ids
            if (is.character(ids)) {
              idx <- match(ids, ref_ids)
            } else if (all(is.wholenumber(ids))) {
              idx <- ids
            } else {
              stop(paste(param, "must be character or integer"))
            }
            idx <- idx[!is.na(idx)]
            ids <- ref_ids[idx]
            return(list(ids=ids, idx=idx))
          }
          processed_rid <- process_ids(rid, g@rid, "rid")
          processed_cid <- process_ids(cid, g@cid, "cid")
          rid <- processed_rid$ids
          ridx <- processed_rid$idx
          cid <- processed_cid$ids
          cidx <- processed_cid$idx
          sdrow <- setdiff(rid, g@rid)
          sdcol <- setdiff(cid, g@cid)
          if (length(sdrow) > 0) warning("the following rids were not found:\n", paste(sdrow, collapse="\n"))
          if (length(sdcol) > 0) warning("the following cids were not found:\n", paste(sdcol, collapse="\n"))
          newg <- g
          # make sure ordering is right
          rid <- g@rid[ridx]
          cid <- g@cid[cidx]
          newg@mat <- matrix(g@mat[ridx, cidx], nrow=length(rid), ncol=length(cid))
          colnames(newg@mat) <- cid
          rownames(newg@mat) <- rid
          # cdesc <- data.frame(g@cdesc)
          # rdesc <- data.frame(g@rdesc)
          # make sure annotations row ordering matches
          # matrix, rid, and cid
          newg@cdesc <- subset_to_ids(g@cdesc, cid)
          newg@rdesc <- subset_to_ids(g@rdesc, rid)
          newg@rid <- rid
          newg@cid <- cid
          if (any(dim(newg@mat) == 0)) {
            warning("one or more returned dimension is length 0
                    check that at least some of the provided rid and/or
                    cid values have matches in the GCT object supplied")
          }
          return(newg)
})

#' Merge two GCT objects together
#'
#' @param g1 the first GCT object
#' @param g2 the second GCT object
#' @param dimension the dimension on which to merge (row or column)
#' @param matrix_only boolean idicating whether to keep only the
#'   data matrices from \code{g1} and \code{g2} and ignore their
#'   row and column meta data
#' @examples
#' # take the first 10 and last 10 rows of an object
#' # and merge them back together
#' (a <- subset.gct(ds, rid=1:10))
#' (b <- subset.gct(ds, rid=969:978))
#' (merged <- merge.gct(a, b, dimension="row"))
#' 
#' @family GCT utilities
#' @export
setGeneric("merge.gct", function(g1, g2, dimension="row", matrix_only=F) {
  standardGeneric("merge.gct")
})
setMethod("merge.gct", signature("GCT", "GCT"),
          function(g1, g2, dimension, matrix_only) {
          # given two gcts objects g1 and g2, merge them
          # on the specified dimension
          if (dimension == "column") dimension <- "col"
          if (dimension == "row") {
            message("appending rows...")
            newg <- g1
            # we're just appending rows so don't need to do anything
            # special with the rid or rdesc. just cat them
            newg@rid <- c(g1@rid, g2@rid)
            newg@rdesc <- data.frame(rbind(data.table::data.table(g1@rdesc), data.table::data.table(g2@rdesc), fill=T))
            # need figure out the index for how to sort the columns of
            # g2@mat so that they are in sync with g1@mat
            idx <- match(g1@cid, g2@cid)
            newg@mat <- rbind(g1@mat, g2@mat[, idx])
            if (!matrix_only) {
              # apply the same sort order to the rows of g2@cdesc so that
              # it's in sync with the final merged matrix
              # figure out which fields are common and keep from the first gct
              cmn_names <- intersect(names(g1@cdesc), names(g2@cdesc))
              newg@cdesc <- cbind(g1@cdesc, g2@cdesc[idx, !(names(g2@cdesc) %in% cmn_names)])
            } else {
              newg@cdesc <- data.frame()
            }
          }
          else if (dimension == "col") {
            message("appending columns...")
            newg <- g1
            # we're just appending columns so don't need to do anything
            # special with cid or cdesc. just cat them
            newg@cid <- c(g1@cid, g2@cid)
            newg@cdesc <- data.frame(rbind(data.table::data.table(g1@cdesc), data.table::data.table(g2@cdesc), fill=T))
            # need figure out the index for how to sort the rows of
            # g2@mat so that they are in sync with g1@mat
            idx <- match(g1@rid, g2@rid)
            newg@mat <- cbind(g1@mat, g2@mat[idx, ])
            if (!matrix_only) {
              # apply the same sort order to the rows of g2@rdesc so that
              # it's in sync with the final merged matrix
              # figure out which fields are common and keep from the first gct
              cmn_names <- intersect(names(g1@rdesc), names(g2@rdesc))
              newg@rdesc <- cbind(g1@rdesc, g2@rdesc[idx, !(names(g2@rdesc) %in% cmn_names)])
            } else {
              newg@rdesc <- data.frame()
            }
          } else {
            stop("dimension must be either row or col")
          }
          return(newg)
})


#' Merge two \code{\link{data.frame}}s, but where there are common fields
#' those in \code{x} are retained and those in \code{y} are dropped.
#' 
#' @param x the \code{\link{data.frame}} whose columns take precedence
#' @param y another \code{\link{data.frame}}
#' @param by a vector of column names to merge on
#' @param allow.cartesian boolean indicating whether it's ok
#'   for repeated values in either table to merge with each other
#'   over and over again.
#' @param as_data_frame boolean indicating whether to ensure
#'   the returned object is a \code{\link{data.frame}} instead of a \code{\link{data.table}}.
#'   This ensures compatibility with GCT object conventions,
#'   that is, the \code{\link{rdesc}} and \code{\link{cdesc}} slots must be strictly
#'   \code{\link{data.frame}} objects.
#'   
#' @return a \code{\link{data.frame}} or \code{\link{data.table}} object
#' 
#' @examples 
#' (x <- data.table(foo=letters[1:10], bar=1:10))
#' (y <- data.table(foo=letters[1:10], bar=11:20, baz=LETTERS[1:10]))
#' # the 'bar' column from y will be dropped on merge
#' cmapR:::merge_with_precedence(x, y, by="foo")
#'
#' @keywords internal
#' @seealso data.table::merge
merge_with_precedence <- function(x, y, by, allow.cartesian=T,
                                  as_data_frame = T) {
  trash <- check_colnames(by, x)
  trash <- check_colnames(by, y)
  # cast as data.tables
  x <- data.table::data.table(x)
  y <- data.table::data.table(y)
  # get rid of row names
  data.table::setattr(x, "rownames", NULL)
  data.table::setattr(y, "rownames", NULL)
  common_cols <- intersect(names(x), names(y))
  y_keepcols <- unique(c(by, setdiff(names(y), common_cols)))
  y <- y[, y_keepcols, with=F]
  # if not all ids match, issue a warning
  if (!all(x[[by]] %in% y[[by]])) {
    warning("not all rows of x had a match in y. some columns may contain NA")
  }
  # merge keeping all the values in x, making sure that the
  # resulting data.table is sorted in the same order as the 
  # original object x
  merged <- merge(x, y, by=by, allow.cartesian=allow.cartesian, all.x=T)
  if (as_data_frame) {
    # cast back to a data.frame if requested
    merged <- data.frame(merged)
  }
  return(merged)
}


#' Add annotations to a GCT object
#' 
#' @description Given a GCT object and either a \code{\link{data.frame}} or
#' a path to an annotation table, apply the annotations to the
#' gct using the given \code{keyfield}.
#' 
#' @param g a GCT object
#' @param annot a \code{\link{data.frame}} or path to text table of annotations
#' @param dimension either 'row' or 'column' indicating which dimension
#'   of \code{g} to annotate
#' @param keyfield the character name of the column in \code{annot} that 
#'   matches the row or column identifiers in \code{g}
#'   
#' @return a GCT object with annotations applied to the specified
#'   dimension
#'   
#' @examples 
#' \dontrun{
#'  g <- parse.gctx('/path/to/gct/file')
#'  g <- annotate.gct(g, '/path/to/annot')
#' }
#' 
#' @family GCT utilities
#' @export
setGeneric("annotate.gct", function(g, annot, dimension="row", keyfield="id") {
  standardGeneric("annotate.gct")
})
setMethod("annotate.gct", signature("GCT"),
          function(g, annot, dimension, keyfield) {
          if (!(any(class(annot) == "data.frame"))) {
            # given a file path, try to read it in
            annot <- fread(annot)
          } else {
            # convert to data.table
            annot <- data.table::data.table(annot)
          }
          # convert the keyfield column to id for merging
          # assumes the gct object has an id field in its existing annotations
          if (!(keyfield %in% names(annot))) {
            stop(paste("column", keyfield, "not found in annotations"))
          } 
          # rename the column to id so we can do the merge
          annot$id <- annot[[keyfield]]
          if (dimension == "column") dimension <- "col"
          if (dimension == "row") {
            orig_id <- g@rdesc$id
            merged <- merge_with_precedence(g@rdesc, annot, by="id", allow.cartesian=T,
                                            as_data_frame=T)
            idx <- match(orig_id, merged$id)
            merged <- merged[idx, ]
            g@rdesc <- merged
          } else if (dimension == "col") {
            orig_id <- g@cdesc$id
            merged <- merge_with_precedence(g@cdesc, annot, by="id", allow.cartesian=T,
                                            as_data_frame=T)
            idx <- match(orig_id, merged$id)
            merged <- merged[idx, ]
            g@cdesc <- merged
          } else {
            stop("dimension must be either row or column")
          }
          return(g)
})


#' Transpose a GCT object
#' 
#' @param g the \code{GCT} object
#' 
#' @return a modified verion of the input \code{GCT} object
#'   where the matrix has been transposed and the row and column
#'   ids and annotations have been swapped.
#'   
#' @examples 
#' transpose.gct(ds)
#' 
#' @family GCT utilties
#' @export
setGeneric("transpose.gct", function(g) {
  standardGeneric("transpose.gct")
})
setMethod("transpose.gct", signature("GCT"), function(g) {
  # transpose matrix
  g@mat <- t(g@mat)
  # create new data
  rid.new <- g@cid
  cid.new <- g@rid
  rdesc.new <- g@cdesc
  cdesc.new <- g@rdesc
  # overwrite g
  g@rid <- rid.new
  g@cid <- cid.new
  g@rdesc <- rdesc.new
  g@cdesc <- cdesc.new
  return(g)
})


#' Convert a GCT object's matrix to ranks
#' 
#' @param g the \code{GCT} object to rank
#' @param dim the dimension along which to rank
#'   (row or column)
#' @param decreasing boolean indicating whether
#'    higher values should get lower ranks
#' 
#' @return a modified version of \code{g}, with the
#'   values in the matrix converted to ranks
#'   
#' @examples 
#' (ranked <- rank.gct(ds, dim="column"))
#' # scatter rank vs. score for a few columns
#' plot(ds@mat[, 1:3], ranked@mat[, 1:3],
#'   xlab="score", ylab="rank")
#' 
#' @family GCT utilities
#' @export
setGeneric("rank.gct", function(g, dim="col", decreasing=T) {
  standardGeneric("rank.gct")
})
setMethod("rank.gct", signature("GCT"), function(g, dim, decreasing=T) {
  # check to make sure dim is allowed
  if (dim=="column") dim <- "col"
  if (!(dim %in% c("row","col"))){
    stop('Dim must be one of row, col')
  }
  # rank along the specified axis. transpose if ranking rows so that the data 
  # comes back in the correct format
  if (decreasing) {
    g@mat <- -1*g@mat
  }
  if (dim == 'row'){
    g@mat <- t(apply(g@mat, 1, function(x) rank(x)))
  } else {
    g@mat <- (apply(g@mat, 2, function(x) rank(x)))
  }
  # done
  return(g)
})

# TODO: update to catch the error thrown in example
#' Check for duplicates in a vector
#' @param x the vector
#' @param name the name of the object to print
#'   in an error message if duplicates are found
#' @return NULL
#' @examples 
#' check_dups(c("a", "b", "c", "a", "d"))
#' @export
check_dups <- function(x, name="") {
  if (anyDuplicated(x)) {
    stop(paste(name, "has duplicated values:",
               paste(x[duplicated(x)], collaps="\n"),
               sep="\n"))
  }
}

#' Pad a matrix with additional rows/columns of NA values
#' 
#' @param m a matrix with unique row and column names
#' @param row_universe a vector with the universe of possible
#'   row names
#' @param col_universe a vector with the universe of possible
#'   column names
#' 
#' @return a matrix
#' 
#' @examples 
#' m <- matrix(rnorm(10), nrow=2)
#' rownames(m) <- c("A", "B")
#' colnames(m) <- letters[1:5]
#' na_pad_matrix(m, row_universe=LETTERS, col_universe=letters)
#' 
#' @export
na_pad_matrix <- function(m, row_universe=NULL, col_universe=NULL) {
  # make sure row/col names are assigned and unique
  if (is.null(colnames(m)) || is.null(rownames(m))) {
    stop("m must have unique row and column names assigned")
  }
  check_dups(rownames(m), name="m rownames")
  check_dups(colnames(m), name="m colnames")
  # get original row and col names
  orig_rows <- rownames(m)
  orig_cols <- colnames(m)
  # figure out which new rows to add
  if (is.null(row_universe)) {
    rows_to_add <- c()
  } else {
    rows_to_add <- setdiff(row_universe, orig_rows)
  }
  # figure out which new columns to add
  if (is.null(col_universe)) {
    cols_to_add <- c()
  } else {
    cols_to_add <- setdiff(col_universe, orig_cols)
  }
  # add new rows
  new_rows <- matrix(NA, ncol=ncol(m), nrow=length(rows_to_add))
  m <- rbind(m, new_rows)
  rownames(m) <- c(orig_rows, rows_to_add)
  # add new columns
  new_cols <- matrix(NA, ncol=length(cols_to_add), nrow=nrow(m))
  m <- cbind(m, new_cols)
  colnames(m) <- c(orig_cols, cols_to_add)
  return(m)
}

#' Align the rows and columns of two (or more) matrices
#' 
#' @param m1 a matrix with unique row and column names
#' @param m2 a matrix with unique row and column names
#' @param ... additional matrices with unique row and
#'   column names
#' @param L a list of matrix objects. If this is given,
#'   m1, m2, and ... are ignored
#' @param na.pad boolean indicating whether to pad the
#'   combined matrix with NAs for rows/columns that are
#'   not shared by m1 and m2.
#' @param as.3D boolean indicating whether to return the
#'   result as a 3D array. If FALSE, will return a list.
#'   
#' @return an object containing the aligned matrices. Will
#'   either be a list or a 3D array
#' 
#' @examples 
#' # construct some example matrices
#' m1 <- matrix(rnorm(20), nrow=4)
#' rownames(m1) <- letters[1:4]
#' colnames(m1) <- LETTERS[1:5]
#' m2 <- matrix(rnorm(20), nrow=5)
#' rownames(m2) <- letters[1:5]
#' colnames(m2) <- LETTERS[1:4]
#' m1
#' m2
#' 
#' # align them, padding with NA and returning a 3D array
#' align_matrices(m1, m2)
#' 
#' # align them, not padding and retuning a list
#' align_matrices(m1, m2, na.pad=FALSE, as.3D=FALSE)
#' 
#' @export
align_matrices <- function(m1, m2, ..., L=NULL, na.pad=TRUE, as.3D=TRUE) {
  # get the additional matrices if given
  if (!is.null(L)) {
    if (is.list(L) && all(unlist(lapply(L, is.matrix)))) {
      matrices <- L
    } else {
      stop("L must be a list of matrices")
    }
  } else {
    matrices <- list(m1, m2, ...) 
  }
  n_matrices <- length(matrices)
  # make sure row/col names are assigned and unique
  lapply(1:n_matrices, function(i) {
    if (is.null(colnames(matrices[[i]])) ||
        is.null(rownames(matrices[[i]]))) {
      stop(paste("matrix", i, "must have unique row and column names"))
    }
    check_dups(rownames(matrices[[i]]), name=paste("matrix", i, "rownames"))
    check_dups(colnames(matrices[[i]]), name=paste("matrix", i, "colnames"))
  })
  # figure out the common rows and columns
  common_rows <- sort(Reduce(intersect, lapply(matrices, rownames)))
  common_cols <- sort(Reduce(intersect, lapply(matrices, colnames)))
  # if we're not NA padding, this is all we need to do
  if (!na.pad) {
    matrices <- lapply(matrices, function(m) {
      m[common_rows, common_cols]
    })
  } else {
    # transform the matrices so that they contain the
    # union of rows/cols padded with NA where needed
    row_universe <- sort(Reduce(union, lapply(matrices, rownames)))
    col_universe <- sort(Reduce(union, lapply(matrices, colnames)))
    matrices <- lapply(matrices, function(m) {
      padded <- na_pad_matrix(m, row_universe=row_universe,
                    col_universe=col_universe)
      # rearrange the rows and columns so they're in a consistent order
      # for each matrix
      padded[row_universe, col_universe]
    })
  }
  # if we're not converting to 3D array, return a list
  if (!as.3D) {
    return(matrices)
  } else {
    # initialize an empty 3D array
    arr3d <- array(NA, dim=c(length(row_universe),
                             length(col_universe),
                             length(matrices)),
                   # set the dimnames using the first matrix
                   # b/c we assume they're the same for all
                   # matrices
                   dimnames=list(rownames(matrices[[1]]),
                              colnames(matrices[[1]]),
                              names(matrices)))
    # and fill with the aligned matrices
    for (i in 1:length(matrices)) {
      arr3d[, , i] <- matrices[[i]]
    }
    return(arr3d)
  }
}

#' Exract elements from a GCT matrix
#' 
#' @param g the GCT object
#' @param row_field the column name in rdesc to search on
#' @param col_field the column name in cdesc to search on
#' @param rdesc a \code{data.frame} of row annotations
#' @param cdesc a \code{data.frame} of column annotations
#' @param row_keyfield the column name of \code{rdesc} to use
#'    for annotating the rows of \code{g}
#' @param col_keyfield the column name of \code{cdesc} to use
#'    for annotating the rows of \code{g}
#' 
#' @description extract the elements from a \code{GCT} object
#'   where the values of \code{row_field} and \code{col_field}
#'   are the same. A concrete example is if \code{g} represents
#'   a matrix of signatures of genetic perturbations, and you wan
#'   to extract all the values of the targeted genes.
#'   
#' @return a list of the following elements
#' \describe{
#'   \item{mask}{a logical matrix of the same dimensions as
#'         \code{ds@mat} indicating which matrix elements have
#'         been extracted}
#'  \item{idx}{an array index into \code{ds@mat}
#'         representing which elements have been extracted}
#'  \item{vals}{a vector of the extracted values}
#'  }
#' 
#' @examples
#' # get the values for all targeted genes from a 
#' # dataset of knockdown experiments 
#' res <- extract.gct(kd_gct, row_field="pr_gene_symbol", col_field="pert_mfc_desc")
#' str(res)
#' stats::quantile(res$vals)
#' 
#' @export
extract.gct <- function(g, row_field, col_field,
                        rdesc=NULL, cdesc=NULL,
                        row_keyfield="id", col_keyfield="id") {
  # annotate the gct object if external annotations have been provided
  if (!is.null(rdesc)) {
    g <- annotate.gct(g, rdesc, dim="row", keyfield=row_keyfield)
    rdesc <- g@rdesc
  }
  if (!is.null(cdesc)) {
    g <- annotate.gct(g, cdesc, dim="col", keyfield=col_keyfield)
    cdesc <- g@cdesc
  }
  rdesc <- data.table::data.table(g@rdesc)
  cdesc <- data.table::data.table(g@cdesc)
  # what are the common values
  common_vals <- intersect(rdesc[[row_field]], cdesc[[col_field]])
  mask <- matrix(F, nrow=nrow(g), ncol=ncol(g))
  for (v in common_vals) {
    ridx <- which(rdesc[[row_field]] == v)
    cidx <- which(cdesc[[col_field]] == v)
    mask[ridx, cidx] <- T
  }
  idx <- which(mask, arr.ind=T)
  vals <- g@mat[mask]
  # data.frame containing the extracted values
  # alongside their row and column annotations
  df = cbind(
    {
      x <- rdesc[idx[, 1], ]
      data.table::setattr(x, "names", paste("row", names(x), sep="_"))
      x
    },
    {
      y <- cdesc[idx[, 2], ]
      data.table::setattr(y, "names", paste("col", names(y), sep="_"))
      y
    })
  df$value <- vals
  return(list(
    mask = mask,
    idx = idx,
    vals = vals,
    df = df 
  ))
}

