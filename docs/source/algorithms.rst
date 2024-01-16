Trimming algorithms
***********************

Manual methods
========================

Custom columns
------------------------
This algorithm eliminates a specified set of columns defined by the user. The set of columns to be removed should be provided as individual column numbers separated by commas, and/or as blocks of consecutive columns indicated by the first and last column numbers separated by a hyphen. In the following example::

-selectcols { n,l,m-k }

where n and l are interpreted single column numbers, while m-k is a range of columns (from column m to column k, both included) to be deleted. Note that column numbering starts from 0. For instance, the command::

-selectcols { 2,7,20-25,80-100 }

will remove columns 2 and 7, along with two blocks of columns ranging from column 20 to 25 and 80 to 100, respectively.


Threshold-based trimming
------------------------
The user can choose to remove all columns that do not meet a specified threshold or a combination of thresholds. The gap threshold (*-gt*) and similarity threshold (*-st*) represent the minimum values of the respective scores explained above and can be used individually or in combination. Like the scores they refer to, both thresholds range from 0 to 1.

trimAl provides two shortcuts to commonly used thresholds: *-nogaps* (equivalent to *-gt 1*), which deletes all columns with at least one gap, and *-noallgaps*, which removes columns composed solely of gaps.

In addition, the user can set a conservation threshold (*-cons*), indicating the minimum percentage of columns from the input alignment that should be retained in the trimmed alignment. This threshold is defined between 0 and 100 and takes precedence over all other thresholds. If any other threshold would result in a trimmed alignment with fewer columns than specified by the conservation threshold, trimAl adds more columns to meet the conservation threshold. These columns are added based on their scores, with a preference for columns with higher scores. In the case of equal scores, columns adjacent to already selected column-blocks and closer to the center of the alignment are added first, 
prioritizing the extension of longer and central blocks.


When provided with a set of multiple sequence alignments, trimAl calculates a consistency score for each alignment in the set. Subsequently, the alignment with the highest score is selected. The chosen alignment can undergo various trimming methods, one of which involves removing columns that exhibit lower consistency across the other alignments. To achieve this, the user can utilize the *-ct* parameter to define the minimum values for the consistency score, within the range of 0 to 1. Any columns not meeting this specified value will be removed. Alternatively, the conservation score, as explained previously, can also be employed here. Moreover, it can be used in conjunction with gap and/or similarity methods.


Automated methods
========================

Gappyout method
------------------------
This method relies on the gap distribution within the multiple sequence alignment (MSA). This method relies on the gap distribution within the Multiple Sequence Alignment (MSA). Initially, the method calculates gap scores for all columns and arranges them based on this score, generating a plot depicting potential gap score thresholds versus the percentage of the alignment below each threshold (see :numref:`gappyout-figure`). In the subsequent step, for every set of three consecutive points on this plot, trimAl computes the slopes between the first and third point, represented by blue lines. Following a comparison of all slopes, trimAl identifies the point with the maximal variation between consecutive slopes, indicated by a vertical red line in :numref:`gappyout-figure`.


.. _gappyout-figure:
.. figure:: _static/gappyout_plot.png
    :name: gappyout-plot
    :width: 500px
    :align: center
    :alt: Gappyout plot

    Example of an internal trimAl plot showing possible gap score thresholds (y axis) versus percentages of alignment length below that threshold (x axis). Thin blue lines indicate slopes computed by the program. The vertical red line indicates the cut-off point
    selected by the gappyout algorithm.


After determining a gap score cut-off point, trimAl removes all columns that do not meet this specified value (see :numref:`example-figure`). In practical terms, this method effectively identifies the bimodal distribution of gap scores (columns rich in gaps and columns with fewer gaps) within an alignment. Subsequently, it eliminates the mode associated with a higher concentration of gaps. Our benchmarks indicate that this method efficiently eliminates a significant portion of poorly aligned regions.


.. _example-figure:
.. figure:: _static/gappyout_example.png
    :name: gappyout-example
    :width: 500px
    :align: center
    :alt: Gappyout example MSA

    An example of an alignment trimmed with the gappyout method. Conserved (grey) and trimmed (white) columns are indicated. This figure has been generated with trimAl -htmlout option.

