# zzz.r HapEstXXR
# mit NAMESPACE
.onLoad <- function (libname, pkgname)
{
  print("R/HapEstXXR starts.")
  require(survival)
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
