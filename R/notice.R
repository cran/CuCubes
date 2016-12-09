.onAttach <- function(libname, pkgname) {
  if (!is.loaded('r_cucubes')) {
    packageStartupMessage("You are running basic CuCubes library. CUDA-accelerated CuCubes library can be obtained at https://featureselector.uco.uwb.edu.pl/pub/cucubes/")
  }
}
