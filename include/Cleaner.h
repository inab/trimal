/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2009-2019
        Fernandez-Rodriguez V.  (victor.fernandez@bsc.es)
        Capella-Gutierrez S.    (salvador.capella@bsc.es)
        Gabaldon, T.            (tgabaldon@crg.es)

    This file is part of trimAl/readAl.

    trimAl/readAl are free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl/readAl are distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl/readAl. If not, see <http://www.gnu.org/licenses/>.

***************************************************************************** */

#ifndef TRIMAL_CLEANER_H
#define TRIMAL_CLEANER_H

#include <functional>
#include <algorithm>
#include <queue>

// Forward declaration
class Alignment;

 /**
 \brief Submodule contained by \link Alignment Alignment \endlink.
 \par "Complementary parameter"
  Cleaning methods may contain a boolean flag argument called<b> 'complementary' </b>.\n
  If set to<i> true</i>, the function will return a new alignment that contains
  the sequences and residues that method would<i> reject</i>.\n
  Otherwise it returns the cleaned version of
  the original alignment this object belongs to.
  \n\n
  The use of this flag - complementary - is <b>not recommended.</b>\n
  Complementarity is now performed when ALL methods have been applied.
 */
class Cleaner {

public:

    /**
     \brief Flag for trimming only on the terminal positions.
     */
    bool terminalGapOnly;
     /**
     \brief Flag for keeping sequences even when they are composed only by gaps.
     */
    bool keepSequences;
     /**
      \brief Block size to use on the cleaning methods.
      */
    int blockSize;
     /**
     \brief left boundary of the alignment
     */
    int left_boundary;
     /**
      \brief right boundary of the alignment
      */
    int right_boundary;
     /**
     \brief Method that selects the best
      cleaning method based on statistics of the alignment.
     */
    int selectMethod();

     /**
     \brief Method to clean an alignment. \n
     It removes sequences that<b> overpass or equals </b>a certain threshold. \n
     The function detects if it would clean too many sequences, 
      and relaxes the threshold until we have enough sequences 
      to achieve the given percentage desired to keep. \n
     To achieve this, the method starts in the middle of the sequence and, 
      alternating sides, adds column one by one.\n
     \param cut
      Cut value to use.
      If a column has a value lower or equal to the cut value, it is removed.
     \param baseLine
      Percent of sequences to keep
     \param gInCol
      Vector that contains the values that will be tested for each column.
     \param complementary
      Whether or not to return the complementary version of the trimmed alignment.
     \return Pointer to the cleaned alignment.
     */
    Alignment *cleanByCutValueOverpass(double cut, float baseLine,
                                          const int *gInCol, bool complementary);

     /**
     \brief Method to clean an alignment. \n
     It removes sequences that<b> fall behind </b>a certain threshold. \n
     The function detects if it would clean too many sequences,
      and relaxes the threshold until we have enough sequences
      to achieve the given percentage desired to keep. \n
     To achieve this, the method starts in the middle of the sequence and,
      alternating sides, adds column one by one.
     \param cut
      Gap cut value to use.
     \param baseLine
      Percent of sequences to keep
     \param ValueVect
      Vector that contains the gaps present on each column.
     \param complementary
      Whether or not to return the complementary version of the trimmed alignment.
     \return Pointer to the cleaned alignment.
     */
    Alignment *cleanByCutValueFallBehind(float cut, float baseLine,
                                            const float *ValueVect,
                                            bool complementary);

     /**
    \brief Method to clean an alignment. \n
    It removes sequences that<b> overpass or equals </b> the gap threshold
      but <b> fall behind </b> the similarity threshold. \n
    The function detects if it would clean too many sequences,
      and relaxes both thresholds until we have enough sequences
      to achieve the given percentage desired to keep. \n
    To achieve this, the method starts in the middle of the sequence and,
      alternating sides, adds column one by one.
    \param cutGaps
      Gap cut value to use.
    \param baseLine
      Percent of sequences to keep
    \param gInCol
      Vector that contains the gaps present on each column.
    \param MDK_Win
      Vector that contains similarity value for each column.
    \param cutCons
      Minimum similarity value to keep a column in the alignment.
    \param complementary
      Whether or not to return the complementary version of the trimmed alignment.
    \return Pointer to the cleaned alignment.
    */
    Alignment *cleanByCutValueOverpassOrEquals(
            double cutGaps, const int *gInCol, float baseLine,
            float cutCons, const float *MDK_Win, bool complementary);

     /**
     \brief Method to clean an alignment. It carries out strict and strictplus.\n
     It removes sequences that<b> overpass </b>the gap threshold but
      <b> fall behind </b>the similarity threshold. \n
     The function recovers those columns that, by themselves would be rejected,
      but it's neighbours (3 of 4) don't.\n
     Column blocks that don't have a minimum size set by the method itself,
      will be removed too.
     \param gapCut
      Gap cut value to use.
     \param gInCol
      Vector that contains the gaps present on each column.
     \param simCut
      Minimum similarity value to keep a column.
     \param MDK_W
      Vector that contains the similarity of each column.
     \param complementary
      Whether or not to return the complementary version of the trimmed alignment.
     \param variable
      Whether to use a variable block length.
      If false, block will be size 5.
      Else, it will use 1% of the alignment length, with a minimum of 3 and maximum of 12.
      This value will be overwritten if blockSize (of this object) is bigger than 0.
     \return Pointer to the cleaned alignment.
     */

    Alignment *cleanStrict(int gapCut, const int *gInCol, float simCut,
                              const float *MDK_W, bool complementary, bool variable);

     /**
      \brief Method to trim an alignment based on a minimum sequence overlap threshold.\n
       The method selects a combination of parameters
       to maximize the final number of columns in the new alignment.\n
      \param minimumOverlap
       Min overlap to keep a sequence.
      \param overlapSeq
       Vector containing the overlap for each column.
      \param complementary
       Whether or not to return the complementary version of the trimmed alignment.
      \return Pointer to the cleaned alignment.
    */
    Alignment *cleanOverlapSeq(float minimumOverlap, float *overlapSeq,
                                  bool complementary);

     /**
      \brief Method to trim an alignment based on the gap distribution values.
       Column blocks that don't have a minimum size set by the method itself,
       will be removed too.
      \param baseLine
       Minimum percentage of columns to conserve in the new alignment.
      \param gapsPct
       Maximum percentage of gaps per column.
      \param complementary
       Whether or not to return the complementary version of the trimmed alignment.
      \return Pointer to the cleaned alignment.
    */
    Alignment *cleanGaps(float baseLine, float gapsPct,
                            bool complementary);

     /**
      \brief Method to trim an alignment based on the similarity distribution values.
      \param baseLine
       Minimum percentage of columns to conserve in the new alignment.
      \param conservationPct
       Minimum value of similarity per column to keep the column.
      \param complementary
       Whether or not to return the complementary version of the trimmed alignment.
      \return Pointer to the cleaned alignment.
    */
    Alignment *cleanConservation(float baseLine, float conservationPct,
                                    bool complementary);

     /**
      \brief Method to trim an alignment based on the
       similarity and gaps distribution values.\n
      \param baseLine
       Minimum percentage of columns to conserve in the new alignment.
      \param GapsPct
       Maximum percentage of gaps per column.
      \param conservationPct
       Minimum value of similarity per column to keep the column.
      \param complementary
       Whether or not to return the complementary version of the trimmed alignment.
      \par "Take in mind"
       If baseLine is too strict, the other two will be relaxed
       to obtain the minimum percentage desired.
      \return Pointer to the cleaned alignment.
*/
    Alignment *clean(float baseLine, float GapsPct,
                        float conservationPct, bool complementary);

     /**
      \brief Method to trim an alignment based on consistency values 
       obtained from a dataset of alignments.\n
       The function computes the optimal parameter combination values
       to trim an alignment based on the consistency value 
       from the comparison among a dataset of alignments with the same sequences.
      \param cutpoint 
       Hint of Gap cut point.
       May be used if it's lower than the minimum percentage threshold.
      \param baseLine 
       Minimum percentage of columns to conserve in the new alignment.
      \param vectValues
       Vector with alignment consistency values
      \param complementary
       Whether or not to return the complementary version of the trimmed alignment.
      \return Pointer to the cleaned alignment.
*/
    Alignment *cleanCompareFile(float cutpoint, float baseLine, float *vectValues, bool complementary);

     /**
      \brief Method to compute the overlap values.\n
        See Alignment::overlaps
      \param overlap
       Overlap threshold.
      \param[out] spuriousVector
       Pointer to the spuriousVector to fill.
      \return   <b> True </b> if the calculation went ok.\n
                <b> False </b> otherwise.\n
      This should happen only if you pass a null pointer instead of a spuriousVector.
*/
    bool calculateSpuriousVector(float overlap, float *spuriousVector);

     /**
      \brief Method to remove sequences missaligned
       with the rest of sequences in the alignment.\n
       For each residue in the sequence, it tests it's similarity.
       If the similarity of that residue is higher than overlapColumn value,
       it counts as a hit for the sequence.\n
      After calculating the number of hits for the sequence,
      it removes the sequence if it has a proportion hits/residues
      lower tan minimumOverlap.
      \param overlapColumn
      Minimum similarity value that a residue needs to be considered a hit.
      \param minimumOverlap
      Minimum proportion of hits that a sequence needs to be kept in the new alignment.
      \param complementary
      Whether or not to return the complementary version of the trimmed alignment.
      \return Pointer to the cleaned alignment.
      \
      */
    Alignment *cleanSpuriousSeq(float overlapColumn, float minimumOverlap, bool complementary);

     /**
      \brief Method that carries the gappyout approach.\n
      This methods calculates the slope in gaps distribution
       on the original alignment.\n
      Then, it compares groups of three consecutive residues,
       searching for the group with the most abrupt change in slope.\n
      When found, the first residue is taken as the cutpoint for the sequences.
      \param complementary
       Whether or not to return the complementary version of the trimmed alignment.
      */
    Alignment *clean2ndSlope(bool complementary);

     /**
     \brief Method to clean an alignment. It carries out strict and strictplus.\n
     The method:
       - Computes gaps values and gap cut point of the alignment.
       - Computes similarity values and similarity cut point of the alignment.
       - Calls the cleanStrict method with these values and returns its output.
     \param complementary
      Whether or not to return the complementary version of the trimmed alignment.
     \param variable Whether to use a variable block length.
     If false, block will be size 5.
     Else, it will use 1% of the alignment length, with a minimum of 3 and maximum of 12.
     This value will be overwritten if blockSize (of this object) is bigger than 0.
     \return Pointer to the cleaned alignment returned by cleanStrict
      using the parameters calculated in this function.
     */
    Alignment *cleanCombMethods(bool complementary, bool variable);

     /**
      \brief Method to remove columns composed only by gaps\n
      This method is specially useful when we remove missaligned sequences
       from a given alignment.
      \param complementary
       Whether or not to return the complementary version of the trimmed alignment.
       \n Although this method contains a complementary flag,
       setting this up would return an alignment full of gaps-only columns.
      */
    Alignment *cleanNoAllGaps(bool complementary);

     /**
      \brief Method to remove columns, expressed as ranges.
      \param columns
       Vector containing the columns to remove.
      \param init
       Where does the vector start. Set to 1 if the vector contains its size as first element.
      \param size
       Size of the columns vector.
      \param complementary
       Whether or not to return the complementary version of the trimmed alignment.
      */
    Alignment *removeColumns(int *columns, int init, int size, bool complementary);

     /**
      \brief Method to remove sequences, expressed as ranges.
      \param seqs
       Vector containing the sequences to remove.
      \param init
       Where does the vector start.
       Set to 1 if the vector contains its size as first element.
      \param size
       Size of the columns vector.
      \param complementary
       Whether or not to return the complementary version of the trimmed alignment.
      */
    Alignment *removeSequences(int *seqs, int init, int size, bool complementary);

     /**
      \brief Method to select the most representative sequence (the longest one)
       for each cluster from the input alignment to generate a new alignment.
      \param identityThreshold
      Threshold used to assign sequences to clusters.\n
       If identity between representative sequence of the cluster and
       sequence to assign is superior to this threshold and
       no other cluster has a better identity with its representative,
       it will be assigned to that cluster.
      */
    Alignment *getClustering(float identityThreshold);

     /**
      \brief Method that calculates the optimal cut point for a given clusters number.\n
       The idea is to obtain a cutpoint that can be used to obtain
       as sequences as clusterNumber
      \param clusterNumber
       Number of representative sequences to obtain.
      \return CutPoint to obtain clusterNumber sequences.
      */
    float getCutPointClusters(int clusterNumber);

     /**
      \brief Method to remove blocks of columns with no rejected residue
        smaller than a given size.
      \param blockSize Minimum size to remove a block from the alignment
      \param original Alignment to apply the removal.
       Minimum size a block has to be to be kept.
      */
    void removeSmallerBlocks(int blockSize, Alignment &original);


     /**
      \brief Method to detect right and left borders.
       Borders are the first column found with no gaps.\n
      Everything between the borders are kept in the trimmed alignments.
      \return   <b>True</b> if everything went ok.
      \n        <b> False </b>if it was not possible to calculate the gap stats.
      */
    bool removeOnlyTerminal();

     /**
     \brief Method that identifies and removes columns and sequences composed only by gaps.
     */
    void removeAllGapsSeqsAndCols(bool seqs = true, bool cols = true);

     /**
      \brief Setter method to Terminal Only Flag.
      \param terminalOnly_
       New vlue of the Terminal Only Flag.
      */
    void setTrimTerminalGapsFlag(bool terminalOnly_);

     /**
     \warning Not In Use
     \brief Boundaries setter
     \param[in] boundaries
      New boundaries values.
     */
    void setBoundaries(int *boundaries);

     /**
      \brief Method to calculate identities between the sequences from the alignment.\n
      See Alignment::identities
      */
    void calculateSeqIdentity();

     /**
      \warning Not In Use
      \brief Method that makes a raw approximation of sequence identity computation.\n
      \note Designed for reducing comparisons for huge alignments.
      */
    void calculateRelaxedSeqIdentity();

     /**
        \brief Method to assign sequences to clusters.\n
        Clusters are calculated following this flow:
          - Select the longest sequence and
                use it as representative of the first cluster.\n
          - Using the previously computed identity values,
                check if second longest sequence should be part of
                cluster of first sequence.
            This is computed using the identity value
                in comparison with a threshold.\n
          - If the second sequence should be part of the existing cluster, we add it.
            Otherwise we create a new cluster and
            use this sequence as representative for this new cluster.\n
          - Continue with the rest of sequences,
                comparing them with the representatives of existing clusters.
            If a sequence can pertain to more than one cluster,
                we choose the one that maximizes the identity value
                with the representative sequence.
        \param maximumIdent
         Identity threshold used to decide if a sequence should be part of
         a cluster or create a new one.
        \return Vector that contains the clustering info.\n
        The first item in the vector contains it's size.
        */
    int *calculateRepresentativeSeq(float maximumIdent);

     /**
      \brief Method for computing the complementary alignment.\n
       Complementary alignment is an alignment containing all
       sequences and columns that the original alignment would reject.\n
      It inverses the saveResidues / saveSequences tags.
      \param residues Whether to reverse resiudes tags.
      \param sequences Whether to reverse sequences tags.
      */
    void computeComplementaryAlig(bool residues, bool sequences);

    /***
     * \brief Method to remove duplicate sequences. It will keep the latest ones.
     */
    void removeDuplicates();

private:

    friend class Alignment;

     /**
      \brief Pointer to the alignment that contains this object
      */
    Alignment *alig;

     /**
      \brief Class Constructor - Called by \link Alignment \endlink
      */
     explicit Cleaner(Alignment *parent);

     /**
     \brief Copy Constructor - Called by \link Alignment \endlink
     */
    Cleaner(Alignment *parent, Cleaner *mold);

};


#endif //TRIMAL_CLEANER_H
