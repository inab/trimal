Basic Installation
==================

The simplest way to compile this package is:

  1. 'cd' to the repository folder.

  2. Type 'cmake .' to create a custom makefile. The CMake build will autodetect the best CPU features available. They can be individually disabled to build a version of the binary without SIMD (for testing):

    cmake . -DDISABLE_SSE2=1 -DDISABLE_AVX2=1 -DDISABLE_NEON=1 will build without SIMD
    cmake . -DDISABLE_AVX2=1 -DDISABLE_NEON=1 will build with SSE2 only
    cmake . -DDISABLE_SSE2=1 -DDISABLE_NEON=1 will build with AVX2 only
    cmake . -DDISABLE_SSE2=1 -DDISABLE_AVX2=1 will build with NEON only
    cmake . will build with AVX2, SSE2 and NEON (and effectively use AVX2)
  
  3. Type 'make' to compile.

  3. Optionally, run trimAl/readAl with the examples into the 'dataset' 
     directory to check the correct installation.

By default, 'make' compiles the source code of trimAl and readAl in the
bin directory. After that, you can either add to PATH the bin directory
or move these files to '/usr/local/bin' or to '/usr/bin' using root privileges.

For more information, please, refer to the documentation:
https://trimal.readthedocs.io