.onLoad <- function(libname, pkgname) {
  methods:::bind_activation(TRUE)
}

.Last.lib <- function(libpath) {
  methods:::bind_activation(FALSE)
}
