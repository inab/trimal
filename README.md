Basic Installation
==================

The simplest way to compile this package is:

  1. 'cd' to the repository folder.
  2. Run this to clone the `cpu_features` submodule:

    git submodule update --init

  3. Type 'cmake .' to create a custom makefile. The CMake build will autodetect the best CPU features available. They can be individually disabled to build a version of the binary without SIMD (for testing):

      3.1.  Build with AVX2, SSE2 and NEON (and effectively use AVX2 on x86 and NEON in arm)

            cmake .
     
      3.2.  Build with SSE2 only
     
            cmake . -DDISABLE_AVX2=1 -DDISABLE_NEON=1
    
      3.3.  Build with AVX2 only

            cmake . -DDISABLE_SSE2=1 -DDISABLE_NEON=1

      3.4. Build with NEON only
          
            cmake . -DDISABLE_SSE2=1 -DDISABLE_AVX2=1

      3.5. Build without SIMD
    
            cmake . -DDISABLE_SSE2=1 -DDISABLE_AVX2=1 -DDISABLE_NEON=1
      
  4. Type 'make' to compile.

  5. Optionally, run trimAl/readAl with the examples into the 'dataset' 
     directory to check the correct installation.

By default, 'make' compiles the source code of trimAl and readAl in the
bin directory. After that, you can either add to PATH the bin directory
or move these files to '/usr/local/bin' or to '/usr/bin' using root privileges.

For more information, please, refer to the documentation:
https://trimal.readthedocs.io
