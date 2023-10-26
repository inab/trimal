Installation
*************

Version 1.4
============
The simplest way to compile this package is:

  1. Go to the directory containing the package's source code ('source').
     ::

     $ cd trimal/source

  2. Compile the package.
     ::

     $ make

  3. Optionally, run trimAl/readAl with the examples into the 'dataset' 
     directory to check the correct installation. It should return the original MSA.
     ::

     $ cd ..
     $ source/trimal -in dataset/example.004.AA.fasta


By default, 'make' compiles the source code of trimAl and readAl in the
current directory. After that, you can either add to PATH the current
directory or move these files to '/usr/local/bin' or to '/usr/bin' using
root privileges.


Version 2.0
============
The simplest way to compile this package is:

  1. Go to the repository folder.
     ::

     $ cd trimal

  2. Create a custom makefile. The CMake build will autodetect the best CPU features available. They can be individually disabled:
     
     2.1. Build with AVX2, SSE2 and NEON (and effectively use AVX2)

     ::

     $ cmake . 

     2.2. Build with SSE2 only

     ::

     $ cmake . -DDISABLE_AVX2=1 -DDISABLE_NEON=1

     2.3. Build with AVX2 only

     ::

     $ cmake . -DDISABLE_SSE2=1 -DDISABLE_NEON=1

     2.4. Build with NEON only

     ::

     $ cmake . -DDISABLE_SSE2=1 -DDISABLE_AVX2=1

     2.5. Build without SIMD

     ::

     $ cmake . -DDISABLE_SSE2=1 -DDISABLE_AVX2=1 -DDISABLE_NEON=1
  
  3. Compile the package.
     ::

     $ make

  4. Optionally, run trimAl/readAl with the examples into the 'dataset' 
     directory to check the correct installation. It should return the original MSA.
     ::

     $ bin/trimal -in dataset/example.004.AA.fasta

By default, 'make' compiles the source code of trimAl and readAl in the
bin directory. After that, you can either add to PATH the bin directory
or move these files to '/usr/local/bin' or to '/usr/bin' using root privileges.