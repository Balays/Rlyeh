library(data.table)

merge_list <- function(dtl, by, all = TRUE,
                       suffixes = NULL) {
  stopifnot(length(dtl)  > 1,
            is.character(by))
  
  ## put a sensible name on every element – used to build suffixes
  if (is.null(names(dtl))) names(dtl) <- paste0("dt", seq_along(dtl))
  
  ## define one suffix per table if user hasn’t supplied any
  if (is.null(suffixes))
    suffixes <- paste0(".", names(dtl))
  
  if (length(suffixes) != length(dtl))
    stop("`suffixes` must have the same length as `dtl`")
  
  ## 1. rename the non-key columns up-front -------------------------------
  dtl <- Map(function(DT, sfx) {
    DT <- copy(DT)                         # never modify originals
    nn  <- setdiff(names(DT), by)
    setnames(DT, nn, paste0(nn, sfx))
    DT
  }, dtl, suffixes)
  
  ## 2. pairwise merges through Reduce ------------------------------------
  Reduce(function(x, y)
    merge(x, y, by = by, all = all),   # key columns are identical
    dtl)
}

## ------------------ Example ---------------------------------------------
DT1 <- data.table(id = 1:3, value = 10:12)
DT2 <- data.table(id = 2:4, value = 20:22)
DT3 <- data.table(id = 1:4, value = 30:33)

out <- merge_list(list(A = DT1, B = DT2, C = DT3),
                  by = "id")

out
#    id value.A value.B value.C
# 1:  1      10      NA      30
# 2:  2      11      20      31
# 3:  3      12      21      32
# 4:  4      NA      22      33
