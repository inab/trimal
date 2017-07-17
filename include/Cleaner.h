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
     
     bool terminalGapOnly;
     bool keepSequences;
     int blockSize;
     
     
     int selectMethod(void);
     
     newAlignment *cleanByCutValue(double, float, const int *, bool);

     newAlignment *cleanByCutValue(float, float, const float *, bool);

     newAlignment *cleanByCutValue(double, const int *, float, float, const float *, bool);

     newAlignment *cleanStrict(int, const int *, float, const float *, bool, bool);

     newAlignment *cleanOverlapSeq(float, float *, bool);

     newAlignment *cleanGaps(float, float, bool);

     newAlignment *cleanConservation(float, float, bool);

     newAlignment *clean(float, float, float, bool);

     newAlignment *cleanCompareFile(float, float, float *, bool);

     bool calculateSpuriousVector(float, float *);

     newAlignment *cleanSpuriousSeq(float, float, bool);

     newAlignment *clean2ndSlope(bool);

     newAlignment *cleanCombMethods(bool, bool);

     newAlignment *cleanNoAllGaps(bool);

     newAlignment *removeColumns(int *, int, int, bool);

     newAlignment *removeSequences(int *, int, int, bool);

     newAlignment *getClustering(float identityThreshold);

     float getCutPointClusters(int);

     void removeSmallerBlocks(int);

     bool removeOnlyTerminal(void);

     newValues removeCols_SeqsAllGaps(void);

     void trimTerminalGaps(bool);
    
     void calculateSeqIdentity(void);

     void calculateRelaxedSeqIdentity(void);
  
     int *calculateRepresentativeSeq(float maximumIdent);
    
     void computeComplementaryAlig(bool, bool);

private:

     friend class newAlignment;

     newAlignment* _alignment;

     Cleaner(newAlignment* parent);
     
     Cleaner(newAlignment* parent, Cleaner* mold);

 };


#endif //TRIMAL_CLEANER_H
