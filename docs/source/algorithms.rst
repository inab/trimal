Algorithms
***********************

To implement various trimming algorithms and heuristics (covered in the next section), trimAl utilizes several scores.

Gap Score
========================
The gap score for a column of size :math:`n` is calculated as the fraction of positions in the column without a gap. The formula is as follows:

    .. math::

        \text{gapScore} = \frac{\sum_{\text{i}=1}^{n} res_i \neq gap}{n}


.. _identity_score:

Identity Score
========================
The identity score for any pair of sequences in the MSA is the ratio of identical residues by the sequence length, skipping positions where both residues are gaps or unknown.

    .. math::

        \text{seqIdentity}[seq_i][seq_j] = \frac{\sum_{\text{k}=1}^{n} {res_{ik}} == {res_{jk}}}{\sum_{\text{k}=1}^{n} {res_{ik}} \neq (gap\ or\ unknown)\ OR\ res_{jk} \neq (gap\ or\ unknown)}

    Here:

    - :math:`seq_i` and :math:`seq_j` are the sequences to be compared.
    - :math:`n` is the lenght of the MSA.
    - :math:`res_{ik}` and :math:`res_{jk}` are the residues in position :math:`k` of sequences :math:`seq_i` and :math:`seq_j`.


Residue Similarity Score
========================

The residue similarity score involves Mean Distance (MD) scores, inspired by Thompson et al. (2001). This score is calculated for each column as follows:

1. **Define a substitution matrix:** By default, trimAl employs the BLOSUM62 matrix for amino acid sequences and the identity matrix for nucleotide residues. However, users can input alternative matrices. We denote it as :math:`subsMat`, where :math:`subsMat[res_k][res_i]` represents the substitution score between residues at positions :math:`k` and :math:`i`.

2. **Calculate the distance matrix:** For each pair of residues :math:`res_i`, :math:`res_j`, of the alphabet (i.e., amino acids or nucleotides), the Euclidean distance is computed. We denote it as :math:`distMat`, where :math:`distMat[res_i][res_j]` represents the Euclidean distance between residues :math:`res_i` and :math:`res_j`. The Euclidean distance measures dissimilarity between residues based on their substitution scores.

    .. math::

        \text{distMat}[res_i][res_j] = \sqrt{\sum_{\text{k}=1}^{n} (\text{subsMat}[res_k][res_i] - \text{subsMat}[res_k][res_j])^2}

    Here:

    - :math:`res_i`, :math:`res_j`, and :math:`res_k` are residues at positions :math:`i`, :math:`j`, and :math:`k` of the alphabet (e.g., alanine, valine, and lysine).
    - :math:`n` is the number of residues of the alphabet.

3. **Calculate sequence identity:** As explained in :ref:`Identity Score <identity_score>`.

4. **Calculate Q value:** For each column of the MSA, :math:`Q` is defined as a weighted measure of the distances between all the residues of a column, taking into account some identity factor (:math:`1 - seqIdentity`) between the sequences of the compared residues.

    .. math::

        Q_{\text{k}} = \frac{\sum_{i=1}^{n} \sum_{j=1}^{n} \text{distMat}[res_{ik}][res_{jk}] \times \text{identityFactor}[seq_i][seq_j]}{\sum_{i=1}^{n} \sum_{j=1}^{n} \text{identityFactor}[seq_i][seq_j]}

    Here:

    - :math:`n` is the number of sequences in the MSA.
    - :math:`distMat[res_{ik}][res_{jk}]` is the distance between residues in position :math:`k` of sequences :math:`seq_i` and :math:`seq_j`.
    - :math:`identityFactor[seq_i][seq_j] = 1 - seqIdentity[seq_i][seq_j]`.

4. **Calculate MD for one specific column:** Finally, the :math:`MD` for one specific column is the exponential function of - :math:`Q`.

    .. math::

        MD_{\text{k}} = e^{-Q_k}

In trimAl v1.2, if the gap score for the column is equal to or less than 0.2, the :math:`MD` score is set to zero. This adjustment penalizes columns with numerous gaps, preventing those with few residues from receiving inflated scores. This penalty was introduced to address potential issues when trimming alignments based solely on similarity information.
