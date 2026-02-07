write_named_object_safe <- function(
    x,
    name,
    dir = ".",
    format = c("tsv", "rds"),
    compress_tsv = FALSE,
    verbose = TRUE,
    ...
) {
  format <- match.arg(format)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  if (format == "rds") {
    fn <- file.path(dir, paste0(name, ".rds"))
    saveRDS(x, fn)
  } else {
    ext <- if (compress_tsv) ".tsv.gz" else ".tsv"
    fn  <- file.path(dir, paste0(name, ext))
    data.table::fwrite(
      x,
      file = fn,
      sep = "\t",
      compress = if (compress_tsv) "gzip" else "none",
      ...
    )
  }

  if (verbose) message("[ok] wrote: ", fn)
  invisible(fn)
}


dump_objects_if_exist <- function(
    obj_names,
    dir = ".",
    format = c("tsv", "rds"),
    compress_tsv = FALSE,
    envir = parent.frame(),
    verbose = TRUE,
    ...
) {
  format <- match.arg(format)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  for (nm in obj_names) {
    if (!exists(nm, envir = envir, inherits = TRUE)) {
      if (verbose) message("[skip] missing: ", nm)
      next
    }

    x <- get(nm, envir = envir, inherits = TRUE)

    if (is.null(x)) {
      if (verbose) message("[skip] NULL: ", nm)
      next
    }

    write_named_object_safe(
      x,
      name = nm,
      dir = dir,
      format = format,
      compress_tsv = compress_tsv,
      verbose = verbose,
      ...
    )
  }

  invisible(TRUE)
}

