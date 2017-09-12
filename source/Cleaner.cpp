//
// Created by bioinfo on 2/06/17.
//
#include "../include/utils.h"
#include "../include/values.h"

#include "../include/Cleaner.h"
#include "../include/StatisticsManager.h"
#include "../include/newAlignment.h"
#include "../include/defines.h"

#include "../include/reportsystem.h"

int Cleaner::selectMethod(void) {

    float mx, avg, maxSeq = 0, avgSeq = 0;
    int i, j;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Ask for the sequence identities assesment */
    if(_alignment -> identities == NULL)
        calculateSeqIdentity();
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Once we have the identities among all possible
     * combinations between each pair of sequence. We
     * compute the average identity as well as the
     * average identity for each sequence with its most
     * similar one */
    for(i = 0; i < _alignment -> sequenNumber; i++) {
        for(j = 0, mx = 0, avg = 0; j < _alignment -> sequenNumber; j++) {
            if(i != j) {
                mx  = mx < _alignment -> identities[i][j] ? _alignment -> identities[i][j] : mx;
                avg += _alignment -> identities[i][j];
            }
        }
        avgSeq += avg/(_alignment -> sequenNumber - 1);
        maxSeq += mx;
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    avgSeq = avgSeq/_alignment -> sequenNumber;
    maxSeq = maxSeq/_alignment -> sequenNumber;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* With the different parameters, we decide wich one
     * is the best automated method, based on a previous
     * simulated data benchmarks, to trim the alig */
    if(avgSeq >= 0.55)      return GAPPYOUT;
    else if(avgSeq <= 0.38) return STRICT;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Sometimes we need to use more parameters to select
     * the best automated method, always based on our
     * benchmarks, to trim the input alignment */
    else {
        if(_alignment -> sequenNumber <= 20) return GAPPYOUT;
        else {
            if((maxSeq >= 0.5) && (maxSeq <= 0.65)) return GAPPYOUT;
            else return STRICT;
        }
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

newAlignment* Cleaner::cleanByCutValue(double cut, float baseLine,
                                       const int *gInCol, bool complementary) {

    int i, j, k, jn, oth, pos, block, *vectAux;
    string *matrixAux, *newSeqsName;
    newAlignment *newAlig;
    newValues counter;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Select the columns with a gaps value less or equal
     * than the cut point. */
    for(i = 0, counter.residues = 0; i < _alignment->residNumber; i++)
        if(gInCol[i] <= cut) counter.residues++;
        else _alignment->saveResidues[i] = -1;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Compute, if it's necessary, the number of columns
     * necessary to achieve the minimum number of columns
     * fixed by coverage parameter. */
    oth = utils::roundInt((((baseLine/100.0) - (float) counter.residues/_alignment->residNumber)) * _alignment->residNumber);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* if it's necessary to recover some columns, we
     * applied this instructions to recover it */
    if(oth > 0) {
        counter.residues += oth;

        /* Allocate memory */
        vectAux = new int[_alignment->residNumber];

        /* Sort a copy of the gInCol vector, and take the value of the column that marks the % baseline */
        utils::copyVect((int *) gInCol, vectAux, _alignment->residNumber);
        utils::quicksort(vectAux, 0, _alignment->residNumber - 1);
        cut = vectAux[(int) ((float)(_alignment->residNumber - 1) * (baseLine)/100.0)];

        /* Deallocate memory */
        delete [] vectAux;
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Fixed the initial size of blocks as 0.5% of
     * _alignmentment's length */
    for(k = utils::roundInt(0.005 * _alignment->residNumber); (k >= 0) && (oth > 0); k--) {

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* We start in the _alignmentment middle then we move on
         * right and left side at the same time. */
        for(i = (_alignment->residNumber/2), j = (i + 1); (((i > 0) || (j < (_alignment->residNumber - 1))) && (oth > 0)); i--, j++) {

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* Left side. Here, we compute the block's size. */
            for(jn = i; ((_alignment->saveResidues[jn] != -1) && (jn >= 0) && (oth > 0)); jn--) ;

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* if block's size is greater or equal than the fixed
             * size then we save all columns that have not been
             * saved previously. */
            if((i - jn) >= k) {
                for( ; ((_alignment->saveResidues[jn] == -1) && (jn >= 0) && (oth > 0)); jn--) {
                    if(gInCol[jn] <= cut) {
                        _alignment->saveResidues[jn] = jn;
                        oth--;
                    } else
                        break;
                }
            }
            i = jn;
            /* ***** ***** ***** ***** ***** ***** ***** ***** */

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* Right side. Here, we compute the block's size. */
            for(jn = j; ((_alignment->saveResidues[jn] != -1) && (jn < _alignment->residNumber) && (oth > 0)); jn++) ;

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* if block's size is greater or equal than the fixed
             * size then we save all columns that have not been
             * saved previously. */
            if((jn - j) >= k) {
                for( ; ((_alignment->saveResidues[jn] == -1) && (jn < _alignment->residNumber) && (oth > 0)); jn++) {
                    if(gInCol[jn] <= cut) {
                        _alignment->saveResidues[jn] = jn;
                        oth--;
                    } else
                        break;
                }
            }
            j = jn;
            /* ***** ***** ***** ***** ***** ***** ***** ***** */
        }
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }

    /* Keep only columns blocks bigger than an input columns block size */
    if(blockSize != 0) {

        /* Traverse all _alignmentment looking for columns blocks greater than LONGBLOCK,
         * everytime than a column is not selected by the trimming method, check
         * whether the current block size until that point is big enough to be kept
         * or it should be removed from the final _alignmentment */
        for(i = 0, pos = 0, block = 0; i < _alignment->residNumber; i++) {
            if(_alignment->saveResidues[i] != -1)
                block++;
            else {
                /* Remove columns from blocks smaller than input blocks size */
                if(block < blockSize)
                    for(j = pos; j <= i; j++)
                        _alignment->saveResidues[j] = -1;
                pos = i + 1;
                block = 0;
            }
        }
        /* Check final block separately since it could happen than last block is not
         * big enough but because the loop end could provoke to ignore it */
        if(block < blockSize)
            for(j = pos; j < i; j++)
                _alignment->saveResidues[j] = -1;
    }

    /* If the flag -terminalony is set, apply a method to look for internal
     * boundaries and get back columns inbetween them, if they exist */
    if(terminalGapOnly == true)
        if(!_alignment->Cleaning->removeOnlyTerminal())
            return NULL;

    /* Once the columns/sequences selection is done, turn it around
     * if complementary flag is active */
    if(complementary == true)
        _alignment-> Cleaning -> computeComplementaryAlig(true, false);

    /* Check for any additional column/sequence to be removed */
    /* Compute new sequences and columns numbers */
    counter = _alignment->Cleaning->removeCols_SeqsAllGaps();


    /* Allocate memory  for selected sequences/columns */
    matrixAux = new string[counter.sequences];
    newSeqsName = new string[counter.sequences];

    /* Fill local allocated memory with previously selected data */
    _alignment->fillNewDataStructure(matrixAux, newSeqsName);

    /* When we have all parameters, we create the new _alignmentment */
    newAlig = new newAlignment(*_alignment);

    newAlig -> sequenNumber = counter.sequences;
    newAlig -> residNumber = counter.residues;

    if (newAlig-> sequences != NULL)
        delete[] newAlig->sequences;
    
    if (newAlig-> seqsName != NULL)
        delete[] newAlig->seqsName;
    
    newAlig -> sequences = matrixAux;
    newAlig -> seqsName = newSeqsName;
    
    delete [] counter.matrix;
    delete [] counter.seqsName;
    
    /* Return the new _alignmentment reference */
    return newAlig;
}

newAlignment* Cleaner::cleanByCutValue(float cut, float baseLine,
                                       const float *ValueVect, bool complementary) {

    int i, j, k, jn, oth, pos, block;
    string *matrixAux, *newSeqsName;
    newAlignment *newAlig;
    newValues counter;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Select the columns with a conservation's value
     * greater than the cut point. */
    for(i = 0, counter.residues = 0; i < _alignment->residNumber; i++)
        if(ValueVect[i] > cut) counter.residues++;
        else _alignment->saveResidues[i] = -1;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Compute, if it's necessary, the number of columns
     * necessary to achieve the minimum number of columns
     * fixed by coverage value. */
    oth = utils::roundInt((((baseLine/100.0) - (float) counter.residues/_alignment->residNumber)) * _alignment->residNumber);
    if(oth > 0) counter.residues += oth;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Fixed the initial size of blocks as 0.5% of
     * _alignmentment's length */
    for(k = utils::roundInt(0.005 * _alignment->residNumber); (k >= 0) && (oth > 0); k--) {

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* We start in the _alignmentment middle then we move on
         * right and left side at the same time. */
        for(i = (_alignment->residNumber/2), j = (i + 1); (((i > 0) || (j < (_alignment->residNumber - 1))) && (oth > 0)); i--, j++) {

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* Left side. Here, we compute the block's size. */
            for(jn = i; ((_alignment->saveResidues[jn] != -1) && (jn >= 0) && (oth > 0)); jn--) ;

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* if block's size is greater or equal than the fixed
             * size then we save all columns that have not been
             * saved previously. */
            /* Here, we only accept column with a conservation's
             * value equal to conservation cut point. */
            if((i - jn) >= k) {
                for( ; ((_alignment->saveResidues[jn] == -1) && (jn >= 0) && (oth > 0)); jn--) {
                    if(ValueVect[jn] == cut) {
                        _alignment->saveResidues[jn] = jn;
                        oth--;
                    } else
                        break;
                }
            }
            i = jn;
            /* ***** ***** ***** ***** ***** ***** ***** ***** */

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* Right side. Here, we compute the block's size. */
            for(jn = j; ((_alignment->saveResidues[jn] != -1) && (jn < _alignment->residNumber) && (oth > 0)); jn++) ;

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* if block's size is greater or equal than the fixed
             * size then we select the column and save the block's
             * size for the next iteraction. it's obvius that we
             * decrease the column's number needed to finish. */
            /* Here, we only accept column with a conservation's
             * value equal to conservation cut point. */
            if((jn - j) >= k) {
                for( ; ((_alignment->saveResidues[jn] == -1) && (jn < _alignment->residNumber) && (oth > 0)); jn++) {
                    if(ValueVect[jn] == cut) {
                        _alignment->saveResidues[jn] = jn;
                        oth--;
                    } else
                        break;
                }
            }
            j = jn;
            /* ***** ***** ***** ***** ***** ***** ***** ***** */
        }
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }

    /* Keep only columns blocks bigger than an input columns block size */
    if(blockSize != 0) {

        /* Traverse all _alignmentment looking for columns blocks greater than LONGBLOCK,
         * everytime than a column is not selected by the trimming method, check
         * whether the current block size until that point is big enough to be kept
         * or it should be removed from the final _alignmentment */
        for(i = 0, pos = 0, block = 0; i < _alignment->residNumber; i++) {
            if(_alignment->saveResidues[i] != -1)
                block++;
            else {
                /* Remove columns from blocks smaller than input blocks size */
                if(block < blockSize)
                    for(j = pos; j <= i; j++)
                        _alignment->saveResidues[j] = -1;
                pos = i + 1;
                block = 0;
            }
        }
        /* Check final block separately since it could happen than last block is not
         * big enough but because the loop end could provoke to ignore it */
        if(block < blockSize)
            for(j = pos; j < i; j++)
                _alignment->saveResidues[j] = -1;
    }

    /* If the flag -terminalony is set, apply a method to look for internal
     * boundaries and get back columns inbetween them, if they exist */
    if(terminalGapOnly == true)
        if(!_alignment->Cleaning->removeOnlyTerminal())
            return NULL;

    /* Once the columns/sequences selection is done, turn it around
     * if complementary flag is active */
    if(complementary == true)
        _alignment-> Cleaning -> computeComplementaryAlig(true, false);

    /* Check for any additional column/sequence to be removed */
    /* Compute new sequences and columns numbers */
    counter = _alignment->Cleaning->removeCols_SeqsAllGaps();

    /* Allocate memory  for selected sequences/columns */
    matrixAux = new string[counter.sequences];
    newSeqsName = new string[counter.sequences];

    /* Fill local allocated memory with previously selected data */
    _alignment->fillNewDataStructure(matrixAux, newSeqsName);

    newAlig = new newAlignment(*_alignment);
    
    newAlig -> sequenNumber = counter.sequences;
    newAlig -> residNumber = counter.residues;
    
    if (newAlig->sequences != NULL)
        delete[] newAlig -> sequences;
    
    if (newAlig -> seqsName != NULL)
        delete[] newAlig -> seqsName;
    
    newAlig -> sequences = matrixAux;
    newAlig -> seqsName = newSeqsName;

//     newAlig->fillMatrices(true);
    delete [] counter.matrix;
    delete [] counter.seqsName;
    /* Return the new _alignmentment reference */
    return newAlig;
}

newAlignment* Cleaner::cleanByCutValue(double cutGaps, const int *gInCol,
                                       float baseLine, float cutCons, const float *MDK_Win, bool complementary) {

    int i, j, k, oth, pos, block, jn, blGaps, *vectAuxGaps;
    string *matrixAux, *newSeqsName;
    float blCons, *vectAuxCons;
    newAlignment *newAlig;
    newValues counter;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Select the columns with a conservation's value
     * greater than the conservation cut point AND
     * less or equal than the gap cut point. */
    for(i = 0, counter.residues = 0; i < _alignment->residNumber; i++)
        if((MDK_Win[i] > cutCons) && (gInCol[i] <= cutGaps)) counter.residues++;
        else _alignment->saveResidues[i] = -1;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Compute, if it's necessary, the number of columns
     * necessary to achieve the minimum number of it fixed
     * by the coverage parameter. */
    oth = utils::roundInt((((baseLine/100.0) - (float) counter.residues/_alignment->residNumber)) * _alignment->residNumber);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If it's needed to add new columns, we compute the
     * news thresholds */
    if(oth > 0) {
        counter.residues += oth;

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* Allocate memory */
        vectAuxCons = new float[_alignment->residNumber];
        vectAuxGaps = new int[_alignment->residNumber];
        /* ***** ***** ***** ***** ***** ***** ***** ***** */

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* Sort a copy of the MDK_Win vector and of the gInCol
         * vector, and take the value of the column that marks
         * the % baseline */
        utils::copyVect((float *) MDK_Win, vectAuxCons, _alignment->residNumber);
        utils::copyVect((int *) gInCol, vectAuxGaps, _alignment->residNumber);

        utils::quicksort(vectAuxCons, 0, _alignment->residNumber-1);
        utils::quicksort(vectAuxGaps, 0, _alignment->residNumber-1);

        blCons = vectAuxCons[(int) ((float)(_alignment->residNumber - 1) * (100.0 - baseLine)/100.0)];
        blGaps = vectAuxGaps[(int) ((float)(_alignment->residNumber - 1) * (baseLine)/100.0)];
        /* ***** ***** ***** ***** ***** ***** ***** ***** */

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* Deallocate memory */
        delete [] vectAuxCons;
        delete [] vectAuxGaps;
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Fixed the initial size of blocks as 0.5% of
     * _alignmentment's length */
    for(k = utils::roundInt(0.005 * _alignment->residNumber); (k >= 0) && (oth > 0); k--) {

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* We start in the _alignmentment middle then we move on
         * right and left side at the same time. */
        for(i = (_alignment->residNumber/2), j = (i + 1); (((i > 0) || (j < (_alignment->residNumber - 1))) && (oth > 0)); i--, j++) {

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* Left side. Here, we compute the block's size. */
            for(jn = i; ((_alignment->saveResidues[jn] != -1) && (jn >= 0) && (oth > 0)); jn--) ;

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* if block's size is greater or equal than the fixed
             * size then we select the column and save the block's
             * size for the next iteraction. it's obvius that we
             * decrease the column's number needed to finish. */
            /* Here, we accept column with a conservation's value
             * greater or equal than the conservation cut point OR
             * less or equal than the gap cut point. */
            if((i - jn) >= k) {
                for( ; ((_alignment->saveResidues[jn] == -1) && (jn >= 0) && (oth > 0)); jn--) {
                    if((MDK_Win[jn] >= blCons) || (gInCol[jn] <= blGaps)) {
                        _alignment->saveResidues[jn] = jn;
                        oth--;
                    } else
                        break;
                }
            }
            i = jn;
            /* ***** ***** ***** ***** ***** ***** ***** ***** */

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* Right side. Here, we compute the block's size. */
            for(jn = j; ((_alignment->saveResidues[jn] != -1) && (jn < _alignment->residNumber) && (oth > 0)); jn++) ;

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* if block's size is greater or equal than the fixed
             * size then we select the column and save the block's
             * size for the next iteraction. it's obvius that we
             * decrease the column's number needed to finish. */
            /* Here, we accept column with a conservation's value
             * greater or equal than the conservation cut point OR
             * less or equal than the gap cut point. */
            if((jn - j) >= k) {
                for( ; ((_alignment->saveResidues[jn] == -1) && (jn < _alignment->residNumber) && (oth > 0)); jn++) {
                    if((MDK_Win[jn] >= blCons) || (gInCol[jn] <= blGaps)) {
                        _alignment->saveResidues[jn] = jn;
                        oth--;
                    } else
                        break;
                }
            }
            j = jn;
            /* ***** ***** ***** ***** ***** ***** ***** ***** */
        }
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }

    /* Keep only columns blocks bigger than an input columns block size */
    if(blockSize != 0) {

        /* Traverse all _alignmentment looking for columns blocks greater than LONGBLOCK,
         * everytime than a column is not selected by the trimming method, check
         * whether the current block size until that point is big enough to be kept
         * or it should be removed from the final _alignmentment */
        for(i = 0, pos = 0, block = 0; i < _alignment->residNumber; i++) {
            if(_alignment->saveResidues[i] != -1)
                block++;
            else {
                /* Remove columns from blocks smaller than input blocks size */
                if(block < blockSize)
                    for(j = pos; j <= i; j++)
                        _alignment->saveResidues[j] = -1;
                pos = i + 1;
                block = 0;
            }
        }
        /* Check final block separately since it could happen than last block is not
         * big enough but because the loop end could provoke to ignore it */
        if(block < blockSize)
            for(j = pos; j < i; j++)
                _alignment->saveResidues[j] = -1;
    }

    /* If the flag -terminalony is set, apply a method to look for internal
     * boundaries and get back columns inbetween them, if they exist */
    if(terminalGapOnly == true)
        if(!_alignment->Cleaning->removeOnlyTerminal())
            return NULL;

    /* Once the columns/sequences selection is done, turn it around
     * if complementary flag is active */
    if(complementary == true)
        _alignment->Cleaning->computeComplementaryAlig(true, false);

    /* Check for any additional column/sequence to be removed */
    /* Compute new sequences and columns numbers */
    counter = _alignment->Cleaning->removeCols_SeqsAllGaps();

    /* Allocate memory  for selected sequences/columns */
    matrixAux = new string[counter.sequences];
    newSeqsName = new string[counter.sequences];

    /* Fill local allocated memory with previously selected data */
    _alignment->fillNewDataStructure(matrixAux, newSeqsName);

    newAlig = new newAlignment(*_alignment);
    
    newAlig -> sequenNumber = counter.sequences;
    newAlig -> residNumber = counter.residues;
    
    if (newAlig->sequences != NULL)
        delete[] newAlig -> sequences;
    
    if (newAlig -> seqsName != NULL)
        delete[] newAlig -> seqsName;
    
    newAlig -> sequences = matrixAux;
    newAlig -> seqsName = newSeqsName;
    
    /* Deallocate local memory */

    delete [] counter.matrix;
    delete [] counter.seqsName;

    /* Return the new _alignmentment reference */
    return newAlig;
}

newAlignment* Cleaner::cleanStrict(int gapCut, const int *gInCol, float simCut,
                                   const float *MDK_W, bool complementary, bool variable) {

    int i, num, lenBlock;
    newAlignment *newAlig;

    /* Reject columns with gaps number greater than the gap threshold. */
    for(i = 0; i < _alignment->residNumber; i++)
        if(gInCol[i] > gapCut)
            _alignment->saveResidues[i] = -1;

    /* Reject columns with similarity score under the threshold. */
    for(i = 0; i < _alignment->residNumber; i++)
        if(MDK_W[i] < simCut)
            _alignment->saveResidues[i] = -1;

    /* Search for those columns that have been removed but have,
     * at least, 3 adjacent columns selected. */

    /* Special cases:
     * Second column */
    if((_alignment->saveResidues[0] != -1) && (_alignment->saveResidues[2] != -1)
            && (_alignment->saveResidues[3] != -1))
        _alignment->saveResidues[1] = 1;
    else
        _alignment->saveResidues[1] = -1;

    /* Second last column  */
    if((_alignment->saveResidues[_alignment->residNumber-1] != -1) && (_alignment->saveResidues[_alignment->residNumber-3] != -1)
            && (_alignment->saveResidues[_alignment->residNumber-4] != -1))
        _alignment->saveResidues[(_alignment->residNumber - 2)] = (_alignment->residNumber - 2);
    else
        _alignment->saveResidues[(_alignment->residNumber - 2)] = -1;

    /* Normal cases */
    for(i = 2, num = 0; i < (_alignment->residNumber - 2); i++, num = 0)
        if(_alignment->saveResidues[i] == -1) {
            if(_alignment->saveResidues[(i - 2)] != -1)
                num++;
            if(_alignment->saveResidues[(i - 1)] != -1)
                num++;
            if(_alignment->saveResidues[(i + 1)] != -1)
                num++;
            if(_alignment->saveResidues[(i + 2)] != -1)
                num++;
            _alignment->saveResidues[i] = (num >= 3) ? i : -1;
        }

    /* Select blocks size based on user input. It can be set either to 5 or to a
     * variable number between 3 and 12 depending on _alignmentment length (1% alig) */
    if(!variable)
        lenBlock = 5;
    else {
        lenBlock = utils::roundInt(_alignment->residNumber * 0.01F);
        lenBlock = lenBlock >  3 ? (lenBlock < 12 ? lenBlock : 12) : 3;
    }

    /* Allow to change minimal block size */
    blockSize = blockSize > 0 ? blockSize : lenBlock;

    /* Keep only columns blocks bigger than either a computed dinamically or
     * set by the user block size */
    _alignment->Cleaning->removeSmallerBlocks(blockSize);

    /* If the flag -terminalony is set, apply a method to look for internal
     * boundaries and get back columns inbetween them, if they exist */
    if(terminalGapOnly == true)
        if(!removeOnlyTerminal())
            return NULL;

    /* Once the columns/sequences selection is done, turn it around
     * if complementary flag is active */
    if(complementary == true)
        computeComplementaryAlig(true, false);

    /* Check for any additional column/sequence to be removed */
    /* Compute new sequences and columns numbers */

    newValues counter;
    counter = _alignment->Cleaning->removeCols_SeqsAllGaps();

    /* Allocate memory  for selected sequences/columns */
    //~ matrixAux = new string[counter.sequences];
    //~ newSeqsName = new string[counter.sequences];

    /* Fill local allocated memory with previously selected data */
    //~ _alignment->fillNewDataStructure(matrixAux, newSeqsName);

    _alignment->fillNewDataStructure(&counter);

    newAlig = new newAlignment(*_alignment);
    
    newAlig -> sequenNumber = counter.sequences;
    newAlig -> residNumber = counter.residues;
    
    if (newAlig->sequences != NULL)
        delete[] newAlig -> sequences;
    
    if (newAlig -> seqsName != NULL)
        delete[] newAlig -> seqsName;
    
    newAlig -> sequences = counter.matrix;
    newAlig -> seqsName = counter.seqsName;
    

    /* Return the new _alignmentment reference */
        
//     newAlig->fillMatrices(true);
    
    return newAlig;
}

newAlignment* Cleaner::cleanOverlapSeq(float minimumOverlap, float *overlapSeq,
                                       bool complementary) {

    string *matrixAux, *newSeqsName;
    newAlignment *newAlig;
    newValues counter;
    int i;

    /* Keep only those sequences with an overlap value equal or greater than
     * the minimum overlap value set by the user.  */
    for(i = 0; i < _alignment->sequenNumber; i++)
        if(overlapSeq[i] < minimumOverlap)
            _alignment->saveSequences[i] = -1;

    /* Once the columns/sequences selection is done, turn it around
     * if complementary flag is active */
    if(complementary == true)
        computeComplementaryAlig(false, true);

    /* Check for any additional column/sequence to be removed */
    /* Compute new sequences and columns numbers */
    counter = _alignment->Cleaning->removeCols_SeqsAllGaps();

    /* Allocate memory  for selected sequences/columns */
    matrixAux = new string[counter.sequences];
    newSeqsName = new string[counter.sequences];

    /* Fill local allocated memory with previously selected data */
    _alignment->fillNewDataStructure(matrixAux, newSeqsName);

    /* When we have all parameters, we create the new _alignmentment */

    newAlig = new newAlignment(*_alignment);
    
    newAlig -> sequenNumber = counter.sequences;
    newAlig -> residNumber = counter.residues;
    
    if (newAlig->sequences != NULL)
        delete[] newAlig -> sequences;
    
    if (newAlig -> seqsName != NULL)
        delete[] newAlig -> seqsName;
    
    newAlig -> sequences = matrixAux;
    newAlig -> seqsName = newSeqsName;
    
    delete [] counter.matrix;
    delete [] counter.seqsName;
//     newAlig->fillMatrices(true);

    /* Return the new _alignmentment reference */
    return newAlig;
}

newAlignment* Cleaner::cleanGaps(float baseLine, float gapsPct, bool complementary) {

    newAlignment *ret;
    double cut;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If gaps statistics are not calculated, we
     * calculate them */
    if(_alignment -> Statistics -> calculateGapStats() == false)
        return NULL;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Obtain the cut point using the given parameters */
    cut = _alignment -> sgaps -> calcCutPoint(baseLine, gapsPct);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Once we have the cut value proposed, we call the
     * appropiate method to clean the newAlignment and, then,
     * generate the new newAlignment. */
    ret = cleanByCutValue(cut, baseLine, _alignment -> sgaps -> getGapsWindow(), complementary);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Return a reference of the new newAlignment */
    return ret;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

newAlignment* Cleaner::cleanConservation(float baseLine, float conservationPct, bool complementary) {

    newAlignment *ret;
    float cut;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If conservation's statistics are not calculated,
     * we calculate them */
    if(_alignment->Statistics->calculateConservationStats() != true)
        return NULL;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Calculate the cut point using the given parameters */
    cut = (float) _alignment->scons -> calcCutPoint(baseLine, conservationPct);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Once we have the cut value, we call the appropiate
     * method to clean the newAlignment and, then, generate
       the new newAlignment */
    ret = cleanByCutValue(cut, baseLine, _alignment->scons -> getMdkwVector(), complementary);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Return a reference of the new newAlignment */
    return ret;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

newAlignment* Cleaner::clean(float baseLine, float GapsPct, float conservationPct, bool complementary) {

    newAlignment *ret;
    float cutCons;
    double cutGaps;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If gaps statistics are not calculated, we calculate
     *  them */
    if(_alignment->Statistics->calculateGapStats() != true)
        return NULL;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If conservation's statistics are not calculated,
     * we calculate them */
    if(_alignment->Statistics->calculateConservationStats() != true)
        return NULL;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Calculate the two cut points using the parameters */
    cutGaps = _alignment->sgaps->calcCutPoint(baseLine, GapsPct);
    cutCons = _alignment->scons->calcCutPoint(baseLine, conservationPct);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Clean the alingment using the two cut values, the
     * gapsWindow and MDK_Windows vectors and the baseline
     * value */
    ret = cleanByCutValue(cutGaps, _alignment->sgaps -> getGapsWindow(), baseLine, cutCons, _alignment->scons -> getMdkwVector(), complementary);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Return a reference of the clean newAlignment object */
    return ret;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

newAlignment* Cleaner::cleanCompareFile(float cutpoint, float baseLine, float *vectValues, bool complementary) {

    newAlignment *ret;
    float cut, *vectAux;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Allocate memory */
    vectAux = new float[_alignment->residNumber];
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Sort a copy of the vectValues vector, and take the
     *  value at 100% - baseline position. */
    utils::copyVect((float *) vectValues, vectAux, _alignment->residNumber);
    utils::quicksort(vectAux, 0, _alignment->residNumber-1);
    cut = vectAux[(int) ((float)(_alignment->residNumber - 1) * (100.0 - baseLine)/100.0)];
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* We have to decide which is the smallest value
     * between the cutpoint value and the value from
     * the minimum percentage threshold */
    cut = cutpoint < cut ? cutpoint : cut;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Deallocate memory */
    delete [] vectAux;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Clean the selected newAlignment using the input parameters. */
    ret = cleanByCutValue(cut, baseLine, vectValues, complementary);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Return a refernce of the new newAlignment */
    return ret;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

bool Cleaner::calculateSpuriousVector(float overlap, float *spuriousVector) {

    int i, j, k, seqValue, ovrlap, hit;
    float floatOverlap;
    char indet;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Compute the overlap */
    floatOverlap = overlap * float(_alignment -> sequenNumber-1);
    ovrlap = int(overlap * (_alignment -> sequenNumber-1));

    if(floatOverlap > float(ovrlap))
        ovrlap++;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If the spurious vectos is NULL, returns false. */
    if(spuriousVector == NULL)
        return false;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Depending on the kind of newAlignment, we have
     * different indetermination symbol */
    if(_alignment -> getAlignmentType() == SequenceTypes::AA)
        indet = 'X';
    else
        indet = 'N';
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* For each newAlignment's sequence, computes its overlap */
    for(i = 0, seqValue = 0; i < _alignment -> sequenNumber; i++, seqValue = 0) {

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* For each newAlignment's column, computes the overlap
         * between the selected sequence and the other ones */
        for(j = 0; j < _alignment -> residNumber; j++) {

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* For sequences are before the sequence selected */
            for(k = 0, hit = 0; k < i; k++) {

                /* ***** ***** ***** ***** ***** ***** ***** ***** */
                /* If the element of sequence selected is the same
                 * that the element of sequence considered, computes
                 * a hit */
                if(_alignment -> sequences[i][j] == _alignment -> sequences[k][j])
                    hit++;

                /* If the element of sequence selected isn't a 'X' nor
                 * 'N' (indetermination) or a '-' (gap) and the element
                 * of sequence considered isn't a  a 'X' nor 'N'
                 * (indetermination) or a '-' (gap), computes a hit */
                else if((_alignment -> sequences[i][j] != indet) && (_alignment -> sequences[i][j] != '-')
                        && (_alignment -> sequences[k][j] != indet) && (_alignment -> sequences[k][j] != '-'))
                    hit++;
            }
            /* ***** ***** ***** ***** ***** ***** ***** ***** */

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* For sequences are after the sequence selected */
            for(k = (i + 1); k < _alignment -> sequenNumber; k++) {

                /* ***** ***** ***** ***** ***** ***** ***** ***** */
                /* If the element of sequence selected is the same
                 * that the element of sequence considered, computes
                 * a hit */
                if(_alignment -> sequences[i][j] == _alignment -> sequences[k][j])
                    hit++;

                /* If the element of sequence selected isn't a 'X' nor
                 * 'N' (indetermination) or a '-' (gap) and the element
                 * of sequence considered isn't a  a 'X' nor 'N'
                 * (indetermination) or a '-' (gap), computes a hit */
                else if((_alignment -> sequences[i][j] != indet) && (_alignment -> sequences[i][j] != '-')
                        && (_alignment -> sequences[k][j] != indet) && (_alignment -> sequences[k][j] != '-'))
                    hit++;
            }
            /* ***** ***** ***** ***** ***** ***** ***** ***** */

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* Finally, if the hit's number divided by number of
             * sequences minus one is greater or equal than
             * overlap's value, computes a column's hit. */
            if(hit >= ovrlap)
                seqValue++;
            /* ***** ***** ***** ***** ***** ***** ***** ***** */
        }

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* For each newAlignment's sequence, computes its spurious's
         * or overlap's value as the column's hits -for that
           sequence- divided by column's number. */
        spuriousVector[i] = ((float) seqValue / _alignment -> residNumber);
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If there is not problem in the method, return true */
    return true;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}


newAlignment* Cleaner::cleanSpuriousSeq(float overlapColumn, float minimumOverlap, bool complementary) {

    float *overlapVector;
    newAlignment *newAlig;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Allocate local memory */
    overlapVector = new float[_alignment->sequenNumber];
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Compute the overlap's vector using the overlap
     * column's value */
    if(!calculateSpuriousVector(overlapColumn, overlapVector))
        return NULL;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Select and remove the sequences with a overlap less
     * than threshold's overlap and create a new _alignmentemnt */
    newAlig = cleanOverlapSeq(minimumOverlap, overlapVector, complementary);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Deallocate local memory */
    delete [] overlapVector;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Return a reference of the clean newAlignment object */
    newAlig->fillMatrices(false);
    return newAlig;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

newAlignment* Cleaner::clean2ndSlope(bool complementarity) {

    newAlignment *ret;
    int cut;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If gaps statistics are not calculated, we calculate
     *  them */
    if(_alignment->Statistics->calculateGapStats() != true)
        return NULL;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* We get the cut point using a automatic method for
     * this purpose. */
    cut = _alignment->sgaps -> calcCutPoint2ndSlope();
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Using the cut point calculates in last steps, we
     * clean the newAlignment and generate a new newAlignment */
    ret = cleanByCutValue(cut, 0, _alignment->sgaps->getGapsWindow(), complementarity);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Returns the new newAlignment. */
    return ret;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

newAlignment* Cleaner::cleanCombMethods(bool complementarity, bool variable) {

    float simCut, first20Point, last80Point, *simil, *vectAux;
    int i, j, acm, gapCut, *positions, *gaps;
    double inic, fin, vlr;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If gaps statistics are not calculated, we calculate
     *  them */
    if(_alignment->Statistics->calculateGapStats() != true)
        return NULL;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Computes the gap cut point using a automatic method
     * and at the same time, we get the gaps values from
     * the newAlignment. */
    gapCut = _alignment->sgaps -> calcCutPoint2ndSlope();
    gaps = _alignment->sgaps -> getGapsWindow();
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If conservation's statistics are not calculated,
     * we calculate them */
    if(_alignment->Statistics->calculateConservationStats() != true)
        return NULL;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Computes the conservations value for each column
     * in the newAlignment. At the same time, the method get
     * the vector with those values. */
    _alignment->scons -> calculateVectors(_alignment->sequences, _alignment->sgaps -> getGapsWindow());
    simil = _alignment->scons -> getMdkwVector();
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Allocate local memory and initializate it to -1 */
    positions = new int[_alignment->residNumber];
    utils::initlVect(positions, _alignment->residNumber, -1);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* The method only selects columns with gaps number
     * less or equal than the gap's cut point. Counts the
     * number of columns that have been selected */
    for(i = 0, acm = 0; i < _alignment->residNumber; i++) {
        if(gaps[i] <= gapCut) {
            positions[i] = i;
            acm++;
        }
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Allocate local memory and save the similaritys
     * values for the columns that have been selected */
    vectAux = new float[acm];
    for(i = 0, j = 0; i < _alignment->residNumber; i++)
        if(positions[i] != -1)
            vectAux[j++] = simil[i];
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Sort the conservation's value vector. */
    utils::quicksort(vectAux, 0, acm-1);

    /* ...and search for the vector points at the 20 and
     * 80% of length. */
    first20Point = 0;
    last80Point  = 0;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    for(i = acm - 1, j = 1; i >= 0; i--, j++) {
        if((((float) j/acm) * 100.0) <= 20.0)
            first20Point = vectAux[i];
        if((((float) j/acm) * 100.0) <= 80.0)
            last80Point = vectAux[i];
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Computes the logaritmic's values for those points.
     * Finally the method computes the similarity cut
     * point using these values. */
    inic = log10(first20Point);
    fin  = log10(last80Point);
    vlr  = ((inic - fin) / 10) + fin;
    simCut = (float) pow(10, vlr);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Clean the newAlignment and generate a new newAlignment
     * object using the gaps cut and the similaritys cut
     *  values */
    newAlignment *ret = cleanStrict(gapCut, _alignment -> sgaps -> getGapsWindow(),
                                    simCut, _alignment -> scons -> getMdkwVector(), 
                                    complementarity, variable);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Deallocate local memory */
    delete [] vectAux;
    delete [] positions;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Return a reference of the new newAlignment */
    return ret;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

newAlignment* Cleaner::cleanNoAllGaps(bool complementarity) {

    newAlignment *ret;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If gaps statistics are not calculated, we calculate
     *  them */
    if(_alignment->Statistics->calculateGapStats() != true)
        return NULL;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* We want to conserve the columns with gaps' number
     * less or equal than sequences' number - 1  */
    ret = cleanByCutValue((_alignment->sequenNumber - 1), 0, _alignment->sgaps->getGapsWindow(), complementarity);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Returns the new newAlignment. */
    return ret;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

float Cleaner::getCutPointClusters(int clusterNumber) {

    float max, min, avg, gMax, gMin, startingPoint, prevValue = 0, iter = 0;
    int i, j, clusterNum, *cluster, **seqs;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If the user wants only one cluster means that all
     * of sequences have to be in the same clauster.
     * Otherwise, if the users wants the maximum number of
     * clusters means that each sequence have to be in their
     * own cluster */
    if(clusterNumber == _alignment->sequenNumber) return 1;
    else if(clusterNumber == 1) return 0;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Ask for the sequence identities assesment */
    if(_alignment->identities == NULL)
        calculateSeqIdentity();
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Compute the maximum, the minimum and the average
     * identity values from the sequences */
    for(i = 0,gMax = 0, gMin = 1, startingPoint = 0; i < _alignment->sequenNumber; i++) {
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        for(j = 0, max = 0, avg = 0, min = 1; j < i; j++) {
            if(max < _alignment->identities[i][j])
                max  = _alignment->identities[i][j];
            if(min > _alignment->identities[i][j])
                min  = _alignment->identities[i][j];
            avg += _alignment->identities[i][j];
        }
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        for(j = i + 1; j < _alignment->sequenNumber; j++) {
            if(max < _alignment->identities[i][j])
                max  = _alignment->identities[i][j];
            if(min > _alignment->identities[i][j])
                min  = _alignment->identities[i][j];
            avg += _alignment->identities[i][j];
        }
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        startingPoint += avg / (_alignment->sequenNumber - 1);
        if(max > gMax) gMax = max;
        if(min < gMin) gMin = min;
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Take the starting point as the average value */
    startingPoint /= _alignment->sequenNumber;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Compute and sort the sequence length */
    seqs = new int*[_alignment->sequenNumber];
    for(i = 0; i < _alignment->sequenNumber; i++) {
        seqs[i] = new int[2];
        seqs[i][0] = utils::removeCharacter('-', _alignment->sequences[i]).size();
        seqs[i][1] = i;
    }
    utils::quicksort(seqs, 0, _alignment->sequenNumber-1);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Create the data structure to store the different
     * clusters for a given thresholds */
    cluster    = new int[_alignment->sequenNumber];
    cluster[0] = seqs[_alignment->sequenNumber - 1][1];
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /**** ***** ***** ***** ***** ***** ***** ***** */
    /* Look for the optimal identity value to get the
     * number of cluster set by the user. To do that, the
     * method starts in the average identity value and moves
     * to the right or lefe side of this value depending on
     * if this value is so strict or not. We set an flag to
     * avoid enter in an infinite loop, if we get the same
     * value for more than 10 consecutive point without get
     * the given cluster number means that it's impossible
     * to get this cut-off and we need to stop the search */
    do {
        clusterNum = 1;
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* Start the search */
        for(i = _alignment->sequenNumber - 2; i >= 0; i--) {
            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            for(j = 0; j < clusterNum; j++)
                if(_alignment->identities[seqs[i][1]][cluster[j]] > startingPoint)
                    break;
            /* ***** ***** ***** ***** ***** ***** ***** ***** */

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            if(j == clusterNum) {
                cluster[j] = seqs[i][1];
                clusterNum++;
            }
            /* ***** ***** ***** ***** ***** ***** ***** ***** */
        }
        /* ***** ***** ***** ***** ***** ***** ***** ***** */

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* Given the cutoff point, if we get the same cluster
         * number or we have been getting for more than 10
         * consecutive times the same cutoff point, we stop */
        if((clusterNum == clusterNumber) || (iter > 10))
            break;
        else {
            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* We have to move to the left side from the average
             * to decrease the cutoff value and get a smaller number
             * of clusters */
            if(clusterNum > clusterNumber) {
                gMax = startingPoint;
                startingPoint = (gMax + gMin)/2;
            } else {
                /* In the opposite side, we have to move to the right
                 * side from the cutting point to be more strict and get
                 * a bigger number of clusters */
                gMin = startingPoint;
                startingPoint = (gMax + gMin)/2;
            }
            /* ***** ***** ***** ***** ***** ***** ***** ***** */

            /* ***** ***** ***** ***** ***** ***** ***** ***** */
            /* If the clusters number from the previous iteration
             * is different from the current number, we put the
             * iteration number to 0 and store this new value */
            if(prevValue != clusterNum) {
                iter = 0;
                prevValue = clusterNum;
                /* Otherwise, we increase the iteration number to
                 * avoid enter in an infinitve loop */
            } else iter++;
            /* ***** ***** ***** ***** ***** ***** ***** ***** */
        }
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
    } while(true);

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Deallocate dinamic memory */
    for(i = 0; i < _alignment->sequenNumber; i++) delete [] seqs[i];
    delete [] seqs;
    delete [] cluster;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    return startingPoint;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Sets the condition to remove only terminal gaps after applying any
 * trimming method or not.   */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void Cleaner::setTrimTerminalGapsFlag(bool terminalOnly_)
{
    terminalGapOnly = terminalOnly_;
}

void Cleaner::setBoundaries(int * boundaries_)
{
    if (boundaries_ != NULL) {
        left_boundary = boundaries_[0];
        right_boundary = boundaries_[1];
    }
}

newAlignment * Cleaner::getClustering(float identityThreshold) {

    string *matrixAux, *newSeqsName;
    int i, j, *clustering;
    newAlignment *newAlig;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Get the representative member for each cluster
     * given a maximum identity threshold */
    clustering = calculateRepresentativeSeq(identityThreshold);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Put all sequences to be deleted and get back those
     * sequences that are representative for each cluster
     * */
    for(i = 0; i < _alignment -> sequenNumber; i ++)
        _alignment -> saveSequences[i] = -1;
    for(i = 1; i <= clustering[0]; i ++)
        _alignment -> saveSequences[clustering[i]] = clustering[i];
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* We allocate memory to save the sequences selected */
    matrixAux = new string[clustering[0]];
    newSeqsName = new string[clustering[0]];

    /* Copy to new structures the information that have
     * been selected previously. */
    for(i = 0, j = 0; i < _alignment -> sequenNumber; i++)
        if(_alignment -> saveSequences[i] != -1) {
            newSeqsName[j] = _alignment -> seqsName[i];
            matrixAux[j] = _alignment -> sequences[i];
            j++;
        }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* When we have all parameters, we create the new
     * alignment */

    newAlig = new newAlignment(*_alignment);
    
    newAlig -> sequenNumber = clustering[0];
    
    if (newAlig->sequences != NULL)
        delete[] newAlig -> sequences;
    
    if (newAlig -> seqsName != NULL)
        delete[] newAlig -> seqsName;
    
    newAlig -> sequenNumber = clustering[0];
    newAlig -> sequences = matrixAux;
    newAlig -> seqsName = newSeqsName;
    
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

     newAlig->fillMatrices(true);
     delete [] clustering;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Return the new alignment reference */
    return newAlig;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}



/* Remove those columns, expressed as range, set by the user. It can return
 * the complementary alignmnet if appropiate flags is set. */
newAlignment* Cleaner::removeColumns(int *columns, int init, int size,
                                     bool complementary) {

    string *matrixAux, *newSeqsName;
    newAlignment *newAlig;
    newValues counter;
    int i, j;

    /* Delete those range columns defines in the columns vector */
    for(i = init; i < size + init; i += 2)
        for(j = columns[i]; j <= columns[i+1]; j++)
            _alignment -> saveResidues[j] = -1;

    /* Once the columns/sequences selection is done, turn it around
     * if complementary flag is active */
    if(complementary == true)
        computeComplementaryAlig(true, false);

    /* Check for any additional column/sequence to be removed */
    /* Compute new sequences and columns numbers */
    counter = _alignment -> Cleaning -> removeCols_SeqsAllGaps();

    /* Allocate memory  for selected sequences/columns */
    matrixAux = new string[counter.sequences];
    newSeqsName = new string[counter.sequences];

    /* Fill local allocated memory with previously selected data */
    _alignment -> fillNewDataStructure(matrixAux, newSeqsName);

    newAlig = new newAlignment(*_alignment);
    
    newAlig -> sequenNumber = counter.sequences;
    newAlig -> residNumber = counter.residues;
    
    if (newAlig->sequences != NULL)
        delete[] newAlig -> sequences;
    
    if (newAlig -> seqsName != NULL)
        delete[] newAlig -> seqsName;
    
    newAlig -> sequences = matrixAux;
    newAlig -> seqsName = newSeqsName;

    /* Deallocate local memory */
//     delete[] matrixAux;
//     delete[] newSeqsName;
    delete [] counter.matrix;
    delete [] counter.seqsName;
//     newAlig->fillMatrices(true);

    /* Return the new alignment reference */
    return newAlig;
}

/* This method removes those sequences, expressed as range of sequences, set by
 * the user. The method also can return the complementary alignment, it means,
 * those sequences that originally should be removed from the input alignment */
newAlignment *Cleaner::removeSequences(int *seqs, int init, int size,
                                       bool complementary) {

    string *matrixAux, *newSeqsName;
    newAlignment *newAlig;
    newValues counter;
    int i, j;

    /* Delete those range of sequences defines by the seqs vector */
    for(i = init; i < size + init; i += 2)
        for(j = seqs[i]; j <= seqs[i+1]; j++)
            _alignment -> saveSequences[j] = -1;

    /* Once the columns/sequences selection is done, turn it around
     * if complementary flag is active */
    if(complementary == true)
        computeComplementaryAlig(false, true);

    /* Check for any additional column/sequence to be removed */
    /* Compute new sequences and columns numbers */
    counter = _alignment -> Cleaning -> removeCols_SeqsAllGaps();

    /* Allocate memory  for selected sequences/columns */
    matrixAux = new string[counter.sequences];
    newSeqsName = new string[counter.sequences];

    /* Fill local allocated memory with previously selected data */
    _alignment -> fillNewDataStructure(matrixAux, newSeqsName);

    /* When we have all parameters, we create the new alignment */

    newAlig = new newAlignment(*_alignment);
    
    newAlig -> sequenNumber = counter.sequences;
    newAlig -> residNumber = counter.residues;
    
    if (newAlig->sequences != NULL)
        delete[] newAlig -> sequences;
    
    if (newAlig -> seqsName != NULL)
        delete[] newAlig -> seqsName;
    
    newAlig -> sequences = matrixAux;
    newAlig -> seqsName = newSeqsName;
    
    delete [] counter.matrix;
    delete [] counter.seqsName;

    /* Return the new alignment reference */
    return newAlig;
}

/* Function designed for identifying right and left borders between central
 * and terminal regions in the newAlignment. The borders are those columns, first
 * and last, composed by only residues. Everything inbetween left and right
 * borders are keept independently of the applied methods */
bool Cleaner::removeOnlyTerminal(void) {

    int i;
  const int *gInCol;

  if((left_boundary == -1) and (right_boundary == -1)) {
    /* Get alignments gaps stats and copy it */
    if(_alignment->Statistics-> calculateGapStats() != true) {
      cerr << endl << "WARNING: Impossible to apply 'terminal-only' method"
        << endl << endl;
      return false;
    }
    gInCol = _alignment->sgaps -> getGapsWindow();

    /* Identify left and right borders. First and last columns with no gaps */
    for(i = 0; i < _alignment->residNumber && gInCol[i] != 0; i++) ;
    left_boundary = i;

    for(i = _alignment->residNumber - 1; i > -1 && gInCol[i] != 0; i--) ;
    right_boundary = i;
  }

  else if(left_boundary >= right_boundary) {
      ReportSystem::Report(ReportSystem::ErrorCode::LeftBoundaryBiggerThanRightBoundary, new std::string[2]{ std::to_string(left_boundary), std::to_string(right_boundary)} );
//     cerr << endl << "ERROR: Check your manually set left '"<< left_boundary
//       << "' and right '" << right_boundary << "' boundaries'" << endl << endl;
    return false;
  }

  /* We increase the right boundary in one position because we use lower strict
   * condition to get back columns inbetween boundaries removed by any method */
   right_boundary += 1;

  /* Once the interal boundaries have been established, if these limits exist
   * then retrieved all columns inbetween both boundaries. Columns out of these
   * limits will remain selected or not depending on the algorithm applied */
  for(i = left_boundary; i < right_boundary; i++)
    _alignment-> saveResidues[i] = i;

  return true;
}

/* Function designed to identify and remove those columns blocks smaller than
 * a given size */
void Cleaner::removeSmallerBlocks(int blockSize) {

    int i, j, pos, block;

    if(blockSize == 0)
        return ;

    /* Traverse the newAlignment looking for blocks greater than BLOCKSIZE, everytime
     * than a column hasn't been selected, check whether the current block is big
     * enough to be kept or it should be removed from the final newAlignment */
    for(i = 0, pos = 0, block = 0; i < _alignment -> residNumber; i++) {
        if(_alignment -> saveResidues[i] != -1)
            block++;
        else {
            /* Remove columns from blocks smaller than input blocks size */
            if(block < blockSize)
                for(j = pos; j <= i; j++)
                    _alignment -> saveResidues[j] = -1;
            pos = i + 1;
            block = 0;
        }
    }

    /* Check final block separately since it could happen than last block is not
     * big enough but because the loop end could provoke to ignore it */
    if(block < blockSize)
        for(j = pos; j < i; j++)
            _alignment -> saveResidues[j] = -1;
    return ;
}

/* Function designed to identify columns and sequences composed only by gaps.
 * Once these columns/sequences have been identified, they are removed from
 * final newAlignment. */
newValues Cleaner::removeCols_SeqsAllGaps(void) {
    int i, j, valid, gaps;
    bool warnings = false;
    newValues counter;

    /* Check all valid columns looking for those composed by only gaps */
    for(i = 0, counter.residues = 0; i < _alignment -> residNumber; i++) {
        if(_alignment -> saveResidues[i] == -1)
            continue;

        for(j = 0, valid = 0, gaps = 0; j < _alignment -> sequenNumber; j++) {
            if (_alignment -> saveSequences[j] == -1)
                continue;
            if (_alignment -> sequences[j][i] == '-')
                gaps ++;
            valid ++;
        }
        /* Once a column has been identified, warm about it and remove it */
        if(gaps == valid) {
            if(!warnings)
                cerr << endl;
            warnings = true;
            ReportSystem::Report(ReportSystem::WarningCode::RemovingOnlyGapsColumn);
//             cerr << "WARNING: Removing column '" << i << "' composed only by gaps"
//                  << endl;
            _alignment -> saveResidues[i] = -1;
        } else {
            counter.residues ++;
        }
    }

    /* Check for those selected sequences to see whether there is anyone with
     * only gaps */
    for(i = 0, counter.sequences = 0; i < _alignment -> sequenNumber; i++) {
        if(_alignment -> saveSequences[i] == -1)
            continue;

        for(j = 0, valid = 0, gaps = 0; j < _alignment -> residNumber; j++) {
            if(_alignment -> saveResidues[j] == -1)
                continue;
            if (_alignment -> sequences[i][j] == '-')
                gaps ++;
            valid ++;
        }
        /* Warm about it and remove each sequence composed only by gaps */
        if(gaps == valid) {
            if(!warnings)
                cerr << endl;
            warnings = true;

            if(keepSequences) {
                ReportSystem::Report(ReportSystem::WarningCode::KeepingOnlyGapsColumn, new std::string[1] { _alignment->seqsName[i] });
//                 cerr << "WARNING: Keeping sequence '" << _alignment -> seqsName[i]
//                      << "' composed only by gaps" << endl;
                counter.sequences ++;
            } else {
                ReportSystem::Report(ReportSystem::WarningCode::RemovingOnlyGapsColumn, new std::string[1] { _alignment->seqsName[i] });
//                 cerr << "WARNING: Removing sequence '" << _alignment -> seqsName[i]
//                      << "' composed only by gaps" << endl;
                _alignment -> saveSequences[i] = -1;
            }
        } else {
            counter.sequences ++;
        }
    }
    if(warnings)
        cerr << endl;

    if (counter.matrix != NULL) delete [] counter.matrix;
    if (counter.seqsName != NULL) delete [] counter.seqsName; 
    
    counter.matrix = new string[counter.sequences];
    counter.seqsName = new string[counter.sequences];

    return counter;
}

void Cleaner::calculateSeqIdentity(void) {

    int i, j, k, hit, dst;
    char indet;

    /* Depending on alignment type, indetermination symbol will be one or other */
    indet = _alignment -> getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

    /* Create identities matrix to store identities scores */
    _alignment -> identities = new float*[_alignment -> sequenNumber];

    /* For each seq, compute its identity score against the others in the MSA */
    for(i = 0; i < _alignment -> sequenNumber; i++) {
        _alignment -> identities[i] = new float[_alignment -> sequenNumber];

        /* It's a symmetric matrix, copy values that have been already computed */
        for(j = 0; j < i; j++)
            _alignment -> identities[i][j] = _alignment -> identities[j][i];
        _alignment -> identities[i][i] = 0;

        /* Compute identity scores for the current sequence against the rest */
        for(j = i + 1; j < _alignment -> sequenNumber; j++) {
            for(k = 0, hit = 0, dst = 0; k < _alignment -> residNumber; k++) {
                /* If one of the two positions is a valid residue,
                 * count it for the common length */
                if(((_alignment -> sequences[i][k] != indet) && (_alignment -> sequences[i][k] != '-')) ||
                        ((_alignment -> sequences[j][k] != indet) && (_alignment -> sequences[j][k] != '-'))) {
                    dst++;
                    /* If both positions are the same, count a hit */
                    if(_alignment -> sequences[i][k] == _alignment -> sequences[j][k])
                        hit++;
                }
            }

            /* Identity score between two sequences is the ratio of identical residues
             * by the total length (common and no-common residues) among them */
            _alignment -> identities[i][j] = (float) hit/dst;
        }
    }
}

void Cleaner::calculateRelaxedSeqIdentity(void) {
    /* Raw approximation of sequence identity computation designed for reducing
     * comparisons for huge alignemnts */

    int i, j, k, hit;

    /* Create identities matrix to store identities scores */
    _alignment -> identities = new float*[_alignment -> sequenNumber];

    /* For each seq, compute its identity score against the others in the MSA */
    for(i = 0; i < _alignment -> sequenNumber; i++) {
        _alignment -> identities[i] = new float[_alignment -> sequenNumber];

        /* It's a symmetric matrix, copy values that have been already computed */
        for(j = 0; j < i; j++)
            _alignment -> identities[i][j] = _alignment -> identities[j][i];
        _alignment -> identities[i][i] = 0;

        /* Compute identity score between the selected sequence and the others */
        for(j = i + 1; j < _alignment -> sequenNumber; j++) {
            for(k = 0, hit = 0; k < _alignment -> residNumber; k++) {
                /* If both positions are the same, count a hit */
                if(_alignment -> sequences[i][k] == _alignment -> sequences[j][k])
                    hit++;
            }
            /* Raw identity score is computed as the ratio of identical residues between
             * alignment length */
            _alignment -> identities[i][j] = (float) hit/_alignment -> residNumber;
        }
    }
}

int* Cleaner::calculateRepresentativeSeq(float maximumIdent) {

    int i, j, pos, clusterNum, **seqs;
    int *cluster;
    static int *repres;
    float max;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Ask for the sequence identities assesment */
    if(_alignment -> identities == NULL)
        _alignment ->  Cleaning ->  calculateSeqIdentity();
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    seqs = new int*[_alignment -> sequenNumber];
    for(i = 0; i < _alignment -> sequenNumber; i++) {
        seqs[i] = new int[2];
        seqs[i][0] = utils::removeCharacter('-', _alignment -> sequences[i]).size();
        seqs[i][1] = i;
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    utils::quicksort(seqs, 0, _alignment -> sequenNumber-1);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    cluster = new int[_alignment -> sequenNumber];
    cluster[0] = seqs[_alignment -> sequenNumber - 1][1];
    clusterNum = 1;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    for(i = _alignment -> sequenNumber - 2; i >= 0; i--) {

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        for(j = 0, max = 0, pos = -1; j < clusterNum; j++) {
            if(_alignment -> identities[seqs[i][1]][cluster[j]] > maximumIdent) {
                if(_alignment -> identities[seqs[i][1]][cluster[j]] > max) {
                    max = _alignment -> identities[seqs[i][1]][cluster[j]];
                    pos = j;
                }
            }
        }
        /* ***** ***** ***** ***** ***** ***** ***** ***** */

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        if(pos == -1) {
            cluster[j] = seqs[i][1];
            clusterNum++;
        }
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    repres = new int[clusterNum + 1];
    repres[0] = clusterNum;
    for(i = 0; i < clusterNum; i++)
        repres[i+1] = cluster[i];
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Deallocate dinamic memory */
    for(i = 0; i < _alignment -> sequenNumber; i++)
        delete [] seqs[i];

    delete [] cluster;
    delete [] seqs;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    return repres;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

void Cleaner::computeComplementaryAlig(bool residues, bool sequences) {
    int i;

    for(i = 0; i < _alignment -> residNumber && residues; i++)
        _alignment -> saveResidues[i] = (_alignment -> saveResidues[i] == -1) ? i : -1;

    for(i = 0; i < _alignment -> sequenNumber && sequences; i++)
        _alignment -> saveSequences[i] = (_alignment -> saveSequences[i] == -1) ? i : -1;
}


Cleaner::Cleaner(newAlignment *parent) {
    _alignment = parent;
    
    terminalGapOnly     = false;
    keepSequences       = false;
    
    blockSize =     0;
    
    left_boundary = -1;
    right_boundary = -1;
}

Cleaner::Cleaner(newAlignment *parent, Cleaner* mold) {
    _alignment = parent;
    
    terminalGapOnly     = mold -> terminalGapOnly;
    keepSequences       = mold -> keepSequences;
    blockSize           = mold -> blockSize;
    left_boundary       = mold -> left_boundary;
    right_boundary      = mold -> right_boundary;
}
