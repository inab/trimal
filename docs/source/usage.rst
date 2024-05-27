Usage
***********************

Here we show all the options that trimAl offers to the user. At the bottom there are some :ref:`examples <some-examples>` of how to use the tool.

Basic usage
=================
    ::
    
    $ trimal -in <inputfile> -out <outputfile> -<trimming_method>


Help Options
=================
.. option:: -h

    Print this information and show some examples.

.. option:: --version

    Print the trimAl version.


Input-Output Options
====================
.. option:: -in <inputfile>

        Input file in several formats.
        Available input formats: clustal, fasta, nexus, phylip32, phylip40, pir

.. option:: -out <outputfile>
        
        Output alignment in the same input format (default stdout).

.. option:: -nbrf

        Output file in NBRF/PIR format

.. option:: -mega

        Output file in MEGA format

.. option:: -nexus

        Output file in NEXUS format

.. option:: -clustal

        Output file in CLUSTAL format

.. option:: -fasta

        Output file in FASTA format

.. option:: -fasta_m10

        Output file in FASTA format.
        Sequences name length up to 10 characters.
        
.. option:: -phylip

        Output file in PHYLIP/PHYLIP4 format.

.. option:: -phylip_m10

        Output file in PHYLIP/PHYLIP4 format.
        Sequences name length up to 10 characters.

.. option:: -phylip_paml

        Output file in PHYLIP format compatible with PAML.

.. option:: -phylip_paml_m10

        Output file in PHYLIP format compatible with PAML.
        Sequences name length up to 10 characters.

.. option:: -phylip3.2

        Output file in PHYLIP3.2 format.

.. option:: -phylip3.2_m10

        Output file in PHYLIP3.2 format.
        Sequences name length up to 10 characters.

Report Output
====================
.. option:: -htmlout <outputfile>

        Get a summary of trimal's work in an HTML file.

.. option:: -colnumbering

        Get the relationship between the columns in the old and new alignment.

Compare Set Options
====================
.. option:: -compareset <inputfile>

        Input list of paths for the files containing the alignments to compare.

.. option:: -forceselect <inputfile>

        Force selection of the given input file in the files comparison method.

Backtranslation Options
=========================
.. option:: -backtrans <inputfile>

        Use a Coding Sequences file to get a backtranslation for a given AA alignment.

.. option:: -ignorestopcodon

        Ignore stop codons in the input coding sequences.
        
.. option:: -splitbystopcodon

        Split input coding sequences up to first stop codon appearance.

Trimming Parameters
=======================

.. option:: --alternative_matrix <name>

        Select an alternative similarity matrix already loaded. Only available 'degenerated_nt_identity'.

.. option:: -matrix <inputfile>

        Input file for user-defined similarity matrix (default is Blosum62).

.. option:: -block <n>

        Minimum column block size to be kept in the trimmed alignment.
        Available with manual and automatic (gappyout) methods.

.. option:: -keepheader

        Keep original sequence header including non-alphanumeric characters.
        Only available for input FASTA format files.

.. option:: -keepseqs

        Keep sequences even if they are composed only by gaps.

.. option:: -complementary

        Get the complementary alignment in residues.
        Reverses the effect of residue trimming:
        All residues that were to be removed are kept and vice versa.

.. option:: -terminalonly

        Only columns out of internal boundaries
        (first and last column without gaps) are
        candidates to be trimmed depending on the applied method.

Trimming Methods
==================

Manual Selection
------------------

.. option:: -selectcols { n,l,m-k }

        Selection of columns to be removed from the alignment.
        Range: [0 - (Number of Columns - 1)]. (see User Guide).

.. option:: -selectseqs { n,l,m-k }

        Selection of sequences to be removed from the alignment.
        Range: [0 - (Number of Sequences - 1)]. (see User Guide).

Manual Trimming - Thresholds
-----------------------------

.. option:: -gt -gapthreshold <n>

        1 - (fraction of gaps in the column).
        Range: [0 - 1]
        Not compatible with -gat.

.. option:: -st -simthreshold <n>

        Minimum average similarity required.
        Range: [0 - 1]

.. option:: -ct -conthreshold <n>

        Minimum consistency value required.
        Range: [0 - 1]

.. option:: -cons <n>

        Minimum percentage of positions
        in the original alignment to conserve.
        Range: [0 - 100]

.. option:: -clusters <n>

        Get the most Nth representatives sequences from a given alignment.
        Range: [1 - (Number of sequences)]

.. option:: -maxidentity <n>

        Get the representatives sequences for a given identity threshold.
        Range: [0 - 1].


Overlap Trimming
------------------

    Overlap is defined as having a gap in both positions,
    an indetermination in both positions, or a residue in both positions.
    It's main purpose is to remove sequences which share only a reduced region,
    whereas the other regions are not shared with the rest of sequences
    in the alignment and filled with gaps.
    Both overlap thresholds (-resoverlap and -seqoverlap) must be provided jointly.

    Ex: Sp8 may be removed from the alignment depending on the thresholds.

    Sp8    =====GLG===========TKSD---NNNNNNNNNNNNNNNNWV=================

    Sp17   --FAYTAPDLLL-IGFLLKTV-ATFG=================DTWFQLWQGLDLNKMPVF

    Sp10   ======DPAVL--FVIMLGTI-TKFS=================SEWFFAWLGLEINMMVII
    
    Sp26   AAAAAAAAALLTYLGLFLGTDYENFA=================AAAANAWLGLEINMMAQI

.. option:: -resoverlap <n>

        Minimum overlap of a positions with other positions in the column
        to be considered a "good position".
        Range: [0 - 1]. (see User Guide).

.. option:: -seqoverlap <n>

        Minimum percentage of "good positions" that a sequence must have
        in order to be conserved.
        Range: [0 - 100](see User Guide).

.. option:: -nogaps

        Remove all positions with gaps in the alignment.

.. option:: -noallgaps

        Remove columns composed only by gaps.

Automated
------------

.. option:: -gappyout

        Use automated selection on "gappyout" mode.
        This method only uses information based on gaps' distribution.

.. option:: -strict

        Use automated selection on "strict" mode.

.. option:: -strictplus

        Use automated selection on "strictplus" mode.
        Optimized for Neighbour Joining phylogenetic tree reconstruction.

.. option:: -automated1

        Use a heuristic selection of the automatic method
        based on similarity statistics. (see User Guide).
        Optimized for Maximum Likelihood phylogenetic tree reconstruction.


Window frame
==================

Window frame size, score of position i is the average of the window (i - n) to (i + n).
Only compatible with manual methods.

.. option:: -w <n>

        General window frame size, applied to all stats.
            Not compatible with specific sizes.

.. option:: -gw <n>

        Window frame size applied to Gaps.

.. option:: -sw <n>

        Window frame size applied to Similarity.

.. option:: -cw <n>

        Window frame size applied to Consistency.

Statistics Output
==================

Statistics to be calculated and outputted by trimAl

.. option:: -sgc

        Print gap scores for each column in the input alignment.

.. option:: -sgt

        Print accumulated gap scores for the input alignment.

.. option:: -ssc

        Print similarity scores for each column in the input alignment.

.. option:: -sst

        Print accumulated similarity scores for the input alignment.

.. option:: -sfc

        Print sum-of-pairs scores for each column from the selected alignment.

.. option:: -sft

        Print accumulated sum-of-pairs scores for the selected alignment.

.. option:: -sident
    
        Print identity scores for all sequences in the input alignment.
        (see User Guide).

.. option:: -soverlap

        Print overlap scores matrix for all sequences in the input alignment.
        (see User Guide).


.. _some-examples:
Some Examples
======================

1. Removes all positions in the alignment with gaps in 10% or more of
   the sequences, unless this leaves less than 60% of original alignment.
   In such case, print the 60% best (with less gaps) positions.
   ::

   $ trimal -in <inputfile> -out <outputfile> -gt 0.9 -cons 60
        
2. As above but, the gap score is averaged over a window starting
   3 positions before and ending 3 positions after each column.
   ::
   
   $ trimal -in <inputfile> -out <outputfile> -gt 0.9 -cons 60 -w 3
        
3. Use an automatic method to decide optimal thresholds, based in the gap scores
   from input alignment. (see User Guide for details).
   ::
   
   $ trimal -in <inputfile> -out <outputfile> -gappyout
        
4. Use automatic methods to decide optimal thresholds, based on the combination
   of gap and similarity scores. (see User Guide for details).
   ::
   
   $ trimal -in <inputfile> -out <outputfile> -strictplus
        
5. Use an heuristic to decide the optimal method for trimming the alignment.
   (see User Guide for details).
   ::
   
   $ trimal -in <inputfile> -out <outputfile> -automated1
        
6. Use residues and sequences overlap thresholds to delete some sequences from the
   alignment. (see User Guide for details).
   ::
   
   $ trimal -in <inputfile> -out <outputfile> -resoverlap 0.8 -seqoverlap 75
        
7. Selection of columns to be deleted from the alignment. The selection can
   be a column number or a column number interval. Start from 0
   ::
   
   $ trimal -in <inputfile> -out <outputfile> -selectcols { 0,2,3,10,45-60,68,70-78 }
        
8. Get the complementary alignment from the alignment previously trimmed.
   ::

   $ trimal -in <inputfile> -out <outputfile> -selectcols { 0,2,3,10,45-60,68,70-78 } -complementary

9. Selection of sequences to be deleted from the alignment. Start from 0
   ::

   $ trimal -in <inputfile> -out <outputfile> -selectseqs { 2,4,8-12 }

10. Select the 5 most representative sequences from the alignment
    ::
        
    $ trimal -in <inputfile> -out <outputfile> -selectseqs { 2,4,8-12 }
