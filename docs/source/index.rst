
.. raw:: html

   <a href="https://github.com/inab/trimal">
      <i class="fa fa-github" style="text-decoration: none;"></i>
   </a>
   
   

Welcome to trimAl's documentation!

This is currently under development. For now you can visit `old documentation <http://inab.github.io/trimal/index.html>`_ and our `Github repository <https://github.com/inab/trimal>`_.


==================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Introduction
==================
This is trimAl's information page. You can also find information related to readAl, a MSA format conversor.
trimAl is a tool for the automated removal of spurious sequences or poorly aligned regions from a multiple sequence alignment.

trimAl can consider several parameters, alone or in multiple combinations, in order to select the most-reliable positions in the alignment.
These include the proportion of sequences with a gap, the level of residue similarity and, if several alignments for the same set of sequences are provided, the consistency level of columns among alignments.
Moreover, trimAl allows to manually select a set of columns and sequences to be removed from the alignment.

trimAl implements a series of automated algorithms that trim the alignment searching for optimum thresholds based on inherent characteristics of the input alignment, to be used so that the signal-to-noise ratio after alignment trimming phase is increased.

Among trimAl's additional features, trimAl allows getting the complementary alignment (columns that were trimmed), to compute statistics from the alignment, to select the output file format , to get a summary of trimAl's trimming in HTML and SVG formats, and many other options. 

Installation
==================


Version 1.4
---------------------
The simplest way to compile this package is:

  1. 'cd' to the directory containing the package's source code ('source').

  2. Type 'make' to compile the package.

  3. Optionally, run trimAl/readAl with the examples into the 'dataset' 
     directory to check the correct installation.

By default, 'make' compiles the source code of trimAl and readAl in the
current directory. After that, you can either add to PATH the current
directory or move these files to '/usr/local/bin' or to '/usr/bin' using
root privileges.


Version 2.0
---------------------
The simplest way to compile this package is:

  1. 'cd' to the repository folder.

  2. Type 'cmake .' to create a custom makefile. The CMake build will autodetect the best CPU features available. They can be individually disabled to build a version of the binary without SIMD (for testing)::

      $ cmake . -DDISABLE_SSE2=1 -DDISABLE_AVX2=1 -DDISABLE_NEON=1 will build without SIMD
      $ cmake . -DDISABLE_AVX2=1 -DDISABLE_NEON=1 will build with SSE2 only
      $ cmake . -DDISABLE_SSE2=1 -DDISABLE_NEON=1 will build with AVX2 only
      $ cmake . -DDISABLE_SSE2=1 -DDISABLE_AVX2=1 will build with NEON only
      $ cmake . will build with AVX2, SSE2 and NEON (and effectively use AVX2)
  
  3. Type 'make' to compile.

  4. Optionally, run trimAl/readAl with the examples into the 'dataset' 
     directory to check the correct installation.

By default, 'make' compiles the source code of trimAl and readAl in the
bin directory. After that, you can either add to PATH the bin directory
or move these files to '/usr/local/bin' or to '/usr/bin' using root privileges.





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
