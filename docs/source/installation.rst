Installation
*************


Version 1.5.0
========================
The binaries can be downloaded for `Linux <https://github.com/inab/trimal/releases/download/v1.5.0/trimAl_Linux_x86-64.zip>`_,
`MacOS <https://github.com/inab/trimal/releases/download/v1.5.0/trimAl_MacOS_x86-64.zip>`_ and `Windows <https://github.com/inab/trimal/releases/download/v1.5.0/trimAl_Windows_x86-64.zip>`_. You may also download
the source code from `Github repository <https://github.com/inab/trimal/releases/tag/v1.5.0>`_ and then compile it yourself.
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
