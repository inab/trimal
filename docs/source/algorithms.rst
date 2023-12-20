Algorithms
***********************

Residue Similarity Score
========================

The residue similarity score involves Mean Distance (MD) scores, inspired by Thompson et al. (2001). This score is calculated for each column as follows:

1. **Define a substitution matrix:** By default, trimAl employs the BLOSUM62 matrix for amino acid sequences and the identity matrix for nucleotide residues. However, users can input alternative matrices. We denote it as :code:`subsMat`, where :code:`subsMat[res_k][res_i]` represents the substitution score between residues at positions :sub:`k` and :sub:`i`.

2. **Calculate the distance matrix:** For each pair of residues :sub:`res_i`, :sub:`res_j`, of the alphabet (i.e., amino acids or nucleotides), the Euclidean distance is computed. We denote it as :code:`distMat`, where :code:`distMat[:sub:`res_i`][:sub:`res_j`]` represents the Euclidean distance between residues :sub:`res_i` and :sub:`res_j`. The Euclidean distance measures dissimilarity between residues based on their substitution scores.

    .. math::

        \text{distMat}[:\text{sub:`res_i`}][:\text{sub:`res_j`}] = \sqrt{\sum_{\text{sub:`res_k`}=1}^{n} (\text{subsMat}[:\text{sub:`res_k`}][:\text{sub:`res_i`}] - \text{subsMat}[:\text{sub:`res_k`}][:\text{sub:`res_j`}])^2}

    Here:

    - :sub:`res_i`, :sub:`res_j`, and :sub:`res_k` are residues at positions :sub:`i`, :sub:`j`, and :sub:`k` of the alphabet (e.g., alanine, valine, and lysine).
    - :math:`n` is the number of residues of the alphabet.

3. **Calculate Q value for each column:** Then, this distance matrix is used to calculate the :math:`Q` value for each column of the multiple sequence alignment. :math:`Q` is a weighted measure of the distances between all the residues of a column, taking into account some identity factor between the sequences of the compared residues.

    .. math::

        Q_{\text{column}} = \frac{\sum_{i=1}^{n} \sum_{j=1}^{n} \text{distMat}[res_i][res_j] \times \text{identityFactor}[seq_i][seq_j]}{\sum_{i=1}^{n} \sum_{j=1}^{n} \text{identityFactor}[seq_i][seq_j]}

    Here:

    - :math:`n` is the number of sequences in the MSA.
    - \text{distMat}[res_i][res_j] is the distance between residues of sequences :sub:`i` and :sub:`j`.
    - \text{identityFactor}[seq_i][seq_j] = 100 - PCID_{ij}, where PCID_{ij} is the percentage of sequence identity between sequences :sub:`s_i` and :sub:`s_j`.

4. **Calculate MD for one specific column:** Finally, the :math:`MD` for one specific column is the exponential function of - :math:`Q`.

    .. math::

        MD_{\text{column}} = e^{-Q}

In trimAl v1.2, if the gap score for the column is equal to or less than 0.2, the :math:`MD` score is set to zero. This adjustment penalizes columns with numerous gaps, preventing those with few residues from receiving inflated scores. This penalty was introduced to address potential issues when trimming alignments based solely on similarity information.
