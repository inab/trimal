Installation
*************

.. warning::
    There may be some dissimilarity in the final MSA after trimming using different versions,
    due to bug fixes incorporated during the development.


Version 1.4.1 (latest)
========================
The binaries can be downloaded for `Linux <https://github.com/inab/trimal/releases/download/v1.4.1/trimAl_Linux_x86-64.zip>`_
and `Windows <https://github.com/inab/trimal/releases/download/v1.4.1/trimAl_Windows_x86-64.zip>`_. You may also download
the source code from `Github repository <https://github.com/inab/trimal/releases/tag/v1.4.1>`_ and then compile it yourself.
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

Development branch
======================================
This is a `development branch <https://github.com/inab/trimal>`_ which includes mainly bug fixes and some new features. However,
no new features have been added for a long time since there is another branch which includes much more new features
as well as performance improvements. 

You may clone this development branch using git::

  $ git clone https://github.com/inab/trimal.git

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


Version 2.0.0 (Release candidate)
=================================
This branch includes much more new features as well as performance improvements. The source code is available at
the `Github repository 2.0_RC branch <https://github.com/inab/trimal/tree/2.0_RC>`_.
You may clone it using git and then checkout the 2.0_RC branch::

  $ git clone https://github.com/inab/trimal.git
  
::
  
  $ cd trimal

::

  $ git checkout 2.0_RC

The simplest way to compile this package is:

  1. Go to the repository folder.
     ::

     $ cd trimal

  2. Run this to clone the `cpu_features` submodule:
     ::

     $ git submodule update --init

  3. Create a custom makefile. The CMake build will autodetect the best CPU features available. They can be individually disabled:
     
     3.1. Build with AVX2, SSE2 and NEON (and effectively use AVX2 on x86 and NEON in arm)

     ::

     $ cmake . 

     3.2. Build with SSE2 only

     ::

     $ cmake . -DDISABLE_AVX2=1 -DDISABLE_NEON=1

     3.3. Build with AVX2 only

     ::

     $ cmake . -DDISABLE_SSE2=1 -DDISABLE_NEON=1

     3.4. Build with NEON only

     ::

     $ cmake . -DDISABLE_SSE2=1 -DDISABLE_AVX2=1

     3.5. Build without SIMD

     ::

     $ cmake . -DDISABLE_SSE2=1 -DDISABLE_AVX2=1 -DDISABLE_NEON=1
  
  4. Compile the package.
     ::

     $ make

  5. Optionally, run trimAl/readAl with the examples into the 'dataset' 
     directory to check the correct installation. It should return the original MSA.
     ::

     $ bin/trimal -in dataset/example.004.AA.fasta

By default, 'make' compiles the source code of trimAl and readAl in the
bin directory. After that, you can either add to PATH the bin directory
or move these files to '/usr/local/bin' or to '/usr/bin' using root privileges.



Version 1.2
============
The source code is available `here <_static/trimal.v1.2rev59.tar.gz>`_.

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
     $ source/trimal -in dataset/example1.phy


By default, 'make' compiles the source code of trimAl and readAl in the
current directory. After that, you can either add to PATH the current
directory or move these files to '/usr/local/bin' or to '/usr/bin' using
root privileges.
