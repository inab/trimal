Scores
***********************

To implement various trimming algorithms and heuristics (covered in the next section), trimAl utilizes several scores.

Gap Score
========================
The gap score for a column of size :math:`n` is calculated as the fraction of positions in the column without a gap. The formula is as follows:

    .. math::

        gapScore = \frac{\sum_{i=1}^{n} res_i \neq gap}{n}


.. _identity_score:
Identity Score
========================
The identity score for any pair of sequences in the MSA is the ratio of identical residues by the sequence length, skipping positions where both residues are gaps or indeterminations.

    .. math::

        seqIdentity[seq_i][seq_j] = \frac{\sum_{k=1}^{n} {res_{ik}} == {res_{jk}}}{\sum_{k=1}^{n} {res_{ik}} \neq (gap\ or\ indetermination)\ OR\ res_{jk} \neq (gap\ or\ indetermination)}

    Here:

    - :math:`seq_i` and :math:`seq_j` are the sequences to be compared.
    - :math:`n` is the length of the MSA.
    - :math:`res_{ik}` and :math:`res_{jk}` are the residues in position :math:`k` of sequences :math:`seq_i` and :math:`seq_j`.


Residue Similarity Score
========================

The residue similarity score involves Mean Distance (MD) scores, inspired by Thompson et al. (2001). This score is calculated for each column as follows:

1. **Define a substitution matrix:** By default, trimAl employs the BLOSUM62 matrix for amino acid sequences and the identity matrix for nucleotide residues. However, users can input alternative matrices. We denote it as :math:`subsMat`, where :math:`subsMat[res_k][res_i]` represents the substitution score between residues at positions :math:`k` and :math:`i`.

2. **Calculate the distance matrix:** For each pair of residues :math:`res_i`, :math:`res_j`, of the alphabet (i.e., amino acids or nucleotides), the Euclidean distance is computed. We denote it as :math:`distMat`, where :math:`distMat[res_i][res_j]` represents the Euclidean distance between residues :math:`res_i` and :math:`res_j`. The Euclidean distance measures dissimilarity between residues based on their substitution scores.

    .. math::

        distMat[res_i][res_j] = \sqrt{\sum_{k=1}^{n} (subsMat[res_k][res_i] - subsMat[res_k][res_j])^2}

    Here:

    - :math:`res_i`, :math:`res_j`, and :math:`res_k` are residues at positions :math:`i`, :math:`j`, and :math:`k` of the alphabet (e.g., alanine, valine, and lysine).
    - :math:`n` is the number of residues of the alphabet.

3. **Calculate sequence identity:** As explained in :ref:`Identity Score <identity_score>`.

4. **Calculate Q value:** For each column of the MSA, :math:`Q` is defined as a weighted measure of the distances between all the residues of a column, taking into account some identity factor (:math:`1 - seqIdentity`) between the sequences of the compared residues.

    .. math::

        Q_{k} = \frac{\sum_{i=1}^{n} \sum_{j=1}^{n} distMat[res_{ik}][res_{jk}] \times identityFactor[seq_i][seq_j]}{\sum_{i=1}^{n} \sum_{j=1}^{n} identityFactor[seq_i][seq_j]}

    Here:

    - :math:`n` is the number of sequences in the MSA.
    - :math:`distMat[res_{ik}][res_{jk}]` is the distance between residues in position :math:`k` of sequences :math:`seq_i` and :math:`seq_j`.
    - :math:`identityFactor[seq_i][seq_j] = 1 - seqIdentity[seq_i][seq_j]`.

5. **Calculate MD for one specific column:** Finally, the :math:`MD` for one specific column is the exponential function of - :math:`Q`.

    .. math::

        MD_{k} = e^{-Q_k}

If the gap score for the column is equal to or less than 0.2, the :math:`MD` score is set to zero. This adjustment penalizes columns with numerous gaps, preventing those with few residues from receiving inflated scores. This penalty was introduced to address potential issues when trimming alignments based solely on similarity information.


Consistency Score
========================
To calculate the column consistency score, it requires multiple alignments for the same set of sequences in the same order. Once you have selected a reference alignment, each aligned residue pair in the reference alignment is compared with the corresponding pairs in the other alignments. For every occurrence of a matched aligned residue pair in the other alignments, a score of 1 is added to the cumulative score. The final column consistency score is obtained by dividing the cumulative score by the total number of alignments considered and the total number of pairs in the reference alignment. This ensures that the final score ranges from 0 (indicating no aligned pair found in the other alignments) to 1 (implying all alignments share the same pairs, signifying full consistency).

    .. math::

        consistencyScore = \frac{\sum_{i=1}^{n}\sum_{j=1}^{m} matchedPairs[MSA_i][pos_j]}{n \cdot m}

    Here:

    - :math:`n` is the number of alignments considered.
    - :math:`m` is the number of pairs in the reference alignment.
    - :math:`matchedPairs[MSA_i][pos_j]` is the the count of matched pairs between the reference alignment and the :math:`i`-th alignment at position :math:`j`.


Overlap Score
========================
The overlap score of a sequence is defined as the ratio of valid residues to the total number of residues in the sequence. The validity of a residue is determined through the following process:

    1. Count the number of sequences in the multiple sequence alignment (MSA) that have the same residue as the reference sequence at a specific position. Include cases where both the reference residue and the considered residue are neither indeterminate nor gaps.
    2. If the number of occurrences is greater than or equal to the residue overlap threshold, the residue is considered valid.


    .. math::

        {overlapScore_i} = \frac{validResidues_i}{n}

    Here:

    - :math:`i` is index of the sequence.
    - :math:`n` is the lenght of the MSA.
    - :math:`validResidues_i` is the number of valid residues in sequence :math:`i`.
    
        


Window size
========================
This value determines the extent of columns on each side of a specified position that trimAl must take into account when calculating gap, similarity, or consistency scores for that position. If a window size is specified, trimAl calculates the average value across all the considered columns.
