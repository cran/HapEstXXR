################################################################################
# library(HapEstXXR)
# Created: November 30, 2011
#
################################################################################

# mit NAMESPACE
.onLoad <- function (libname, pkgname)
{
  library.dynam("HapEstXXR", pkgname, libname)

}
# .onGenerics <- TRUE
.onUnLoad <- function (libpath) { library.dynam.unload ("HapEstXXR",libpath) }

# ohne NAMESPACE
#.First.lib <- function(lib, pkg) {
#
#  library.dynam("HapEst", pkg, lib)
#  require("MASS")
#  require("haplo.stats")
#  require("gap")
#}
