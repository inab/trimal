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
      \todo Give a good description
      */
     bool terminalGapOnly;
     /**
      \todo Give a good description
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
      \todo Give a good description
      */
     newAlignment *cleanByCutValue(double, float, const int *, bool);
      /**
      \todo Give a good description
      */
     newAlignment *cleanByCutValue(float, float, const float *, bool);
      /**
      \todo Give a good description
      */
     newAlignment *cleanByCutValue(double, const int *, float, float, const float *, bool);
      /**
      \todo Give a good description
      */
     newAlignment *cleanStrict(int, const int *, float, const float *, bool, bool);
     /**
      \todo Give a good description
      */
     newAlignment *cleanOverlapSeq(float, float *, bool);
/**
      \todo Give a good description
      */
     newAlignment *cleanGaps(float, float, bool);
/**
      \todo Give a good description
      */
     newAlignment *cleanConservation(float, float, bool);
/**
      \todo Give a good description
      */
     newAlignment *clean(float, float, float, bool);
/**
      \todo Give a good description
      */
     newAlignment *cleanCompareFile(float, float, float *, bool);
/**
      \todo Give a good description
      */
     bool calculateSpuriousVector(float, float *);
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
