BI-MAP is distributed as LGPL-licensed source code and have to be compiled at your system to be used.

BI-MAP was developed and tested to run on Linux operating system, the other OSes were not tested, although it should be possible to compile BI-MAP on the other platforms.

Please contact the authors if you have problems installing the software.

## Software Requirements ##

To compile and run BI-MAP the following packages (including their _-devel_ versions) need to be installed:
  * **[GCC](http://gcc.gnu.org/) 4.4+** C++ compiler
  * **[cmake](http://cmake.org) 2.6+** build system
  * **[boost](http://www.boost.org) 1.50+** set of C++ libraries (for parallel version of BI-MAP you would need to `boost::mpi` library)
  * **[GSL](http://www.gnu.org/s/gsl/) 1.10+** Various mathematical and statistical functions
  * **[R](http://r-project.org) 2.13+** (optional) for results analysis and visualization
  * **[OpenMPI](http://www.open-mpi.org) 1.5+** (optional) for parallel execution of BI-MAP (other implementations of MPI should also work)
  * **[GoogleTest](http://code.google.com/p/googletest/) 1.5+** for automatic unit testing

The following **R** packages are required to use BI-MAP with R:
  * **[Rcpp](http://dirk.eddelbuettel.com/code/rcpp.html) 0.9+** for loading BI-MAP results into R session (optional)
  * ...

## Compilation ##

**Step 1** Unpack the source archive
```
> tar -xzf bi-map-<version>.tar.gz
```

**Step 2** Create a build folder
```
> cd bi-map-<version>
bi-map> mkdir build
bi-map> cd build
```

**Step 3** Configure a build system environment
```
build> cmake
```
Make sure all the required packages are installed. Please check the output of cmake to find out what packages are missing in your system. Sometimes cmake is not able to automatically find the installed package and its location should be manually specified. Please consult [cmake documentation](http://www.cmake.org/cmake/help/v2.8.8/cmake.html#command:find_package) for the configuration options for a specific package. Upon successful configuration the /UNIX Makefile/ is generated.

**Step 4** Build the package.
```
build> make
```
Upon successful execution the `build/src` folder should contain `BIMAP-sampler` executable.