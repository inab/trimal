//
// Created by bioinfo on 2/06/17.
//

#ifndef TRIMAL_CLEANER_H
#define TRIMAL_CLEANER_H

class newAlignment;
struct newValues;

//Class to make cleaning operations to an alignment
class Cleaner {

 public:
     
     /**
      \brief Flag for making only trimming on the terminal positions.
      */
     bool terminalGapOnly;
     /**
      \todo Give a good descriptiom
      */
     bool keepSequences;
    /**
      \brief Block size to use on the cleaning methods.
      */
     int blockSize;
     
     /**
      \brief Function that selects the best method based on statistics of the alignment.
      */
     int selectMethod(void);
     /**
      \brief Method to clean an alignment. \n
      It removes sequences that<b> overpass or equals </b>a certain threshold. \n
      The function detects if it would clean too many sequences, and relaxes the threshold until we have enough sequences to achieve the given percentage desired to keep. \n
      To achieve this, the method starts in the middle of the sequence and, alternating sides, adds column one by one.\n
      \param cut Cut value to use. If a column has a value lower or equal to the cut value, it is removed.
      \param baseLine Percent of sequences to keep
      \param gInCol Vector that contains the values that will be tested for each column.
      \param complementary Wheter or not to return the complementary version of the trimmed alignment.
      \return Pointer to the cleaned alignment.
      */
     newAlignment *cleanByCutValue(double cut, float baseLine, const int *gInCol, bool complementary);
     /**
      \brief Method to clean an alignment. \n
      It removes sequences that<b> fall behind </b>a certain threshold. \n
      The function detects if it would clean too many sequences, and relaxes the threshold until we have enough sequences to achieve the given percentage desired to keep. \n
      To achieve this, the method starts in the middle of the sequence and, alternating sides, adds column one by one.
      \param cut Gap cut value to use.
      \param baseLine Percent of sequences to keep
      \param ValueVect Vector that contains the gaps present on each column.
      \param complementary Wheter or not to return the complementary version of the trimmed alignment.
      \return Pointer to the cleaned alignment.
      */
     newAlignment *cleanByCutValue(float cut, float baseLine, const float *ValueVect, bool complementary);
           /**
      \brief Method to clean an alignment. \n
      It removes sequences that<b> overpass or equals </b> the gap threshold but <b> fall behind </b> the similarity threshold. \n
      The function detects if it would clean too many sequences, and relaxes both thresholds until we have enough sequences to achieve the given percentage desired to keep. \n
      To achieve this, the method starts in the middle of the sequence and, alternating sides, adds column one by one.
      \param cutGaps Gap cut value to use.
      \param baseLine Percent of sequences to keep
      \param gInCol Vector that contains the gaps present on each column.
      \param MDK_Win Vector that contains similarity value for each column.
      \param cutCons Minium convervation value to keep a column in the alignment.
      \param complementary Wheter or not to return the complementary version of the trimmed alignment.
      \return Pointer to the cleaned alignment.
      */
     newAlignment *cleanByCutValue(double cutGaps, const int * gInCol, float baseLine, float cutCons, const float *MDK_Win, bool complementary);
     /**
      \brief Method to clean an alignment. It carries out strict and strictplus.\n
      It removes sequences that<b> overpass </b>the gap threshold but<b> fall behind </b>the similarity threshold. \n
      The function recovers those columns that, by themselves would be rejected, but it's neighbours (3 of 4) don't.\n
      Column blocks that don't have a minimum size set by the method itself, will be removed too.
      \param gapCut Gap cut value to use.
      \param gInCol Vector that contains the gaps present on each column.
      \param simCut Minimim similarity value to keep a column.
      \param MDK_W Vector that contains the similarity of each column.
      \param complementary Wheter or not to return the complementary version of the trimmed alignment.
      \param variable Wheter to use a variable block length. 
      If false, block will be size 5. 
      Else, it will use 1% of the alignment length, with a minimum of 3 and maximum of 12. 
      This value will be overwritten if blockSize (of this object) is bigger than 0.
      \return Pointer to the cleaned alignment.
      */
     newAlignment *cleanStrict(int gapCut, const int * gInCol, float simCut, const float * MDK_W, bool complementary, bool variable);
/**
      \brief Method to clean sequences based on a minimum overlap threshold.\n
      The method selects a combination of parameters to maximize the final number of columns in the new alignment.\n
      \param minimumOverlap Min overlap to keep a sequence.
      \param overlapSeq Vector containing the overlap for each column.
      \param complementary Wheter or not to return the complementary version of the trimmed alignment.
      \return Pointer to the cleaned alignment.
*/
     newAlignment *cleanOverlapSeq(float minimumOverlap, float * overlapSeq, bool complementary);
/**
      \brief Method to clean an alignment based on the gap distribution values.
      Column blocks that don't have a minimum size set by the method itself, will be removed too.
      \param baseLine Minimim percentage of columns to conserve in the new alignment.
      \param gapsPct Maximum percentage of gaps per column.
      \param complementary Wheter to use a variable block length. 
      \return Pointer to the cleaned alignment.
*/
     newAlignment *cleanGaps(float baseLine, float gapsPct, bool complementary);
/**
      \brief Method to clean an alignment based on the similarity distribution values.
      \param baseLine Minimim percentage of columns to conserve in the new alignment.
      \param conservationPct Minimum value of similarity per column to keep the column.
      \param complementary Wheter to use a variable block length. 
      \return Pointer to the cleaned alignment.
*/
     newAlignment *cleanConservation(float baseLine, float conservationPct, bool complementary);
/**
      \brief Method to clean an alignment based on the similarity and gaps distribution values.\n
      \param baseLine Minimim percentage of columns to conserve in the new alignment.
      \param GapsPct Maximum percentage of gaps per column.
      \param conservationPct Minimum value of similarity per column to keep the column.
      \param complementary Wheter to use a variable block length. 
      \par "Take in mind" If baseLine is too strict, the other two will be relaxed to obtain the minimum percentage desired.
      \return Pointer to the cleaned alignment.
*/
     newAlignment *clean(float baseLine, float GapsPct, float conservationPct, bool complementary);
/**
      \brief Method to clean an alignment based on consistency values obtained from a dataset of alignments.\n
      The function computes the optimal paremeter combination values to trim an alignment based on the consistency value from the comparison among a dataset of alignments with the same sequences.
      \param cutpoint Lower limit (0-1) of comparefile value admits in the new alignment
      \todo cutpoint argument description is not clear.
      \param baseLine Minimim percentage of columns to conserve in the new alignment.
      \param vectValues Vector with alignment consistency values
      \param complementary Wheter to use a variable block length. 
      \return Pointer to the cleaned alignment.
*/
     newAlignment *cleanCompareFile(float cutpoint, float baseLine, float * vectValues, bool complementary);
/**
      \brief Method to compute the overlap for each sequence given a overlap value.\n
      This overlap sets the minimum fraction that has to be a position for the selected sequence to count as a hit.\n
      This proportion measures how much a position is similar (in terms of residues, indetermination and gaps) in comparison with the element on it's columns.
      \param overlap Overlap value to keep a residue.
      \param spuriousVector Pointer to the spuriousVector to fill.
      \return <b> True </b> if the calculation went ok.<b> False </b> otherwise.\n
      This should happen only if you pass a null pointer instead of a spuriousVector.
*/
     bool calculateSpuriousVector(float overlap, float *spuriousVector);
/**
      \todo Give a good description
      */
     newAlignment *cleanSpuriousSeq(float, float, bool);
/**
      \todo Give a good description
      */
     newAlignment *clean2ndSlope(bool);
/**
      \todo Give a good description
      */
     newAlignment *cleanCombMethods(bool, bool);
/**
      \todo Give a good description
      */
     newAlignment *cleanNoAllGaps(bool);
/**
      \todo Give a good description
      */
     newAlignment *removeColumns(int *, int, int, bool);
/**
      \todo Give a good description
      */
     newAlignment *removeSequences(int *, int, int, bool);
/**
      \todo Give a good description
      */
     newAlignment *getClustering(float identityThreshold);
/**
      \todo Give a good description
      */
     float getCutPointClusters(int);
/**
      \todo Give a good description
      */
     void removeSmallerBlocks(int);
/**
      \todo Give a good description
      */
     bool removeOnlyTerminal(void);
/**
      \todo Give a good description
      */
     newValues removeCols_SeqsAllGaps(void);
/**
      \todo Give a good description
      */
     void setTrimTerminalGapsFlag(bool);
    /**
      \todo Give a good description
      */
     void calculateSeqIdentity(void);
/**
      \todo Give a good description
      */
     void calculateRelaxedSeqIdentity(void);
  /**
      \todo Give a good description
      */
     int *calculateRepresentativeSeq(float maximumIdent);
    /**
      \todo Give a good description
      */
     void computeComplementaryAlig(bool, bool);
/**
      \todo Give a good description
      */
private:

     friend class newAlignment;
    /**
      \brief Pointer to the alignment that contains this object
      */
     newAlignment* _alignment;
/**
      \brief Class Constructor - Called by newAlignment
      */
     Cleaner(newAlignment* parent);
     /**
      \brief Copy Constructor - Called by newAlignment
      */
     Cleaner(newAlignment* parent, Cleaner* mold);

 };


#endif //TRIMAL_CLEANER_H
