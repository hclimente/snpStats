.onLoad <- function(libname, package) {
  library.dynam("snpAssoc", package)
  methods:::bind_activation(TRUE)
}

.Last.lib <- function(libname, package) {
  methods:::bind_activation(FALSE)
}
