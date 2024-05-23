
.. raw:: html

   <a href="https://github.com/inab/trimal">
      <i class="fa fa-github" style="text-decoration: none;"></i>
   </a>
   
.. toctree::
   :hidden:
   :maxdepth: 3
   
   installation
   usage
   scores
   algorithms
   benchmarking


Welcome to trimAl's documentation!

This is currently under development. For now you can visit `old documentation <http://inab.github.io/trimal/index.html>`_ and our `Github repository <https://github.com/inab/trimal>`_.

trimAl
==================
This is trimAl's information page. You can also find information related to readAl, a MSA format conversor.
trimAl is a tool for the automated removal of spurious sequences or poorly aligned regions from a multiple sequence alignment.

trimAl can consider several parameters, alone or in multiple combinations, in order to select the most-reliable positions in the alignment.
These include the proportion of sequences with a gap, the level of residue similarity and, if several alignments for the same set of sequences are provided, the consistency level of columns among alignments.
Moreover, trimAl allows to manually select a set of columns and sequences to be removed from the alignment.

trimAl implements a series of automated algorithms that trim the alignment searching for optimum thresholds based on inherent characteristics of the input alignment, to be used so that the signal-to-noise ratio after alignment trimming phase is increased.

Among trimAl's additional features, trimAl allows getting the complementary alignment (columns that were trimmed), to compute statistics from the alignment, to select the output file format , to get a summary of trimAl's trimming in HTML and SVG formats, and many other options. 



Publications
============
- `trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses <https://doi.org/10.1093/bioinformatics/btp348>`_
- `trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses (pdf) <https://academic.oup.com/bioinformatics/article-pdf/25/15/1972/48994574/bioinformatics_25_15_1972.pdf>`_
- `Supplementary material (pdf) <_static/supplementary_material.pdf>`_

Citing this tool
==================
::

   Capella-Gutiérrez, S., Silla-Martínez, J. M., & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics (Oxford, England), 25(15), 1972–1973. https://doi.org/10.1093/bioinformatics/btp348


License
========
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, the last available version.

Contact
========
If you have any doubt or problem feel free to open a `Github issue <https://github.com/inab/trimal/issues>`_ in the repository.

- :email:`Nicolás Díaz Roussel <nicolas.diazroussel@bsc.es>`
- :email:`Salvador Capella Gutierrez <salvador.capella@bsc.es>`


Development team
================
- Nicolás Díaz Roussel
- Salvador Capella Gutiérrez
- Toni Gabaldón
- Víctor Fernández Rodríguez (former)
- Jose M. Silla Martínez (former)

