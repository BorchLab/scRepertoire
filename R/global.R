.onLoad <- function (libname, pkgname) {
  utils::globalVariables(c(
    ".", "Abundance", "CTaa", "CTnt", "Clones", "Freq", "Sample",
    "barcode", "cdr3", "cdr3_nt", "cdr3_nt1", "cdr3_nt2",
    "ci_lower", "ci_upper", "clone", "cloneSize",
    "clonalFrequency", "clonalProportion",
    "cluster", "contigType", "group", "include", "ncells",
    "position", "productive", "prop", "scaled", "size",
    "total", "value", "values"
  ))
  invisible ()
}
