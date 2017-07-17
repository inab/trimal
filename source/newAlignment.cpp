/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated newAlignment trimming in large-scale
                 phylogenetics analyses.

    readAl v1.4: a tool for automated newAlignment conversion among different
                 formats.

    statAl v1.4: a tool for getting stats about multiple sequence newAlignments.


    2009-2012 Capella-Gutierrez S. and Gabaldon, T.
              [scapella, tgabaldon]@crg.es

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

***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
using namespace std;

#include <float.h>
#include "../include/newAlignment.h"
//#include "rwnewAlignment.cpp"
//#include "autnewAlignment.cpp"
#include "../include/Cleaner.h"
#include "../include/StatisticsManager.h"
#include "../include/ReadWriteManager.h"
#include <errno.h>
#include <ctype.h>
#include <string>
//#include <utils.h>
//#include <values.h>
#include "../include/defines.h"
//extern int errno;

using namespace std;


/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Class constructor */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

newAlignment::newAlignment(void) {

    Cleaning = new Cleaner(this);
    Statistics = new StatisticsManager(this);
    ReadWrite = new ReadWriteManager(this);

    /* newAlignment parameter */
    sequenNumber = 0;
    residNumber =  0;

    /* Are the input sequences aligned? */
    isAligned = false;

    /* Should the output file be reversed? */
    reverse   = false;

    /* Should be trimmed only terminal gaps? */
//     terminalGapOnly = false;

    /* Input and output formats */
//     ReadWrite -> iformat = 0;
//     ReadWrite -> oformat = 0;
//     ReadWrite -> shortNames = false;

//     forceCaps = false;
//     upperCase = false;
//     lowerCase = false;

    /* Indicate whether sequences composed only by gaps should be kept or not */
//     keepSequences = false;

    /* Indicate whether original header, they may include non-alphanumerical
     * characters, should be dumped into output stream without any preprocessing
     * step */
//     keepHeader = false;

//     gapSymbol = "-";

    /* Sequence datatype: DNA, RNA or Protein */
    dataType = 0;

    /* Window sizes to trim the input newAlignment */
//     ghWindow = 0;
//     shWindow = 0;

    /* Minimum block size in the new newAlignment */
//     blockSize = 0;

    /* Is this alignmnet new? */
//     oldnewAlignment  = false;

    /* Sequence residues number */
    residuesNumber = NULL;

    /* Columns and sequences that have been previously selected */
    saveResidues  = NULL;
    saveSequences = NULL;

    /* Input sequences as well other information such as sequences name, etc */
    sequences = NULL;
    seqsName  = NULL;
    seqsInfo  = NULL;

    /* Information about input newAlignment */
    filename = "";
//     aligInfo = "";

    /* Information computed from newAlignment */
    sgaps =     NULL;
    scons =     NULL;
    SequencesMatrix = NULL;
    identities = NULL;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Class constructor */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

// newAlignment::newAlignment(string o_filename, string o_aligInfo, string *o_sequences, string *o_seqsName,
//                            string *o_seqsInfo, int o_sequenNumber, int o_residNumber, int o_ReadWrite -> iformat, int o_ReadWrite -> oformat,
//                            bool o_ReadWrite -> shortNames, int o_dataType, int o_isAligned, bool o_reverse, bool o_terminalGapOnly,
//                            bool o_keepSeqs, bool o_keepHeader, int OldSequences, int OldResidues, int *o_residuesNumber,
//                            int *o_saveResidues, int *o_saveSequences, int o_ghWindow, int o_shWindow, int o_blockSize,
//                            float **o_identities) {
// 
//     Cleaning = new Cleaner(this);
//     Statistics = new StatisticsManager(this);
//     ReadWrite = new ReadWriteManager(this);
// 
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
//     int i, j, k, ll;
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
// 
//     oldnewAlignment = true;
// 
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
//     /* Assign the parameter values to the variables */
//     sequenNumber = o_sequenNumber;
//     residNumber  = o_residNumber;
// 
//     ReadWrite -> iformat = o_ReadWrite -> iformat;
//     ReadWrite -> oformat = o_ReadWrite -> oformat;
//     ReadWrite -> shortNames = o_ReadWrite -> shortNames;
// 
//     dataType = o_dataType;
// 
//     ghWindow = o_ghWindow;
//     shWindow = o_shWindow;
// 
//     blockSize = o_blockSize;
// 
//     isAligned = o_isAligned;
//     reverse   = o_reverse;
// 
//     terminalGapOnly = o_terminalGapOnly;
// 
//     filename = o_filename;
//     aligInfo = o_aligInfo;
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
// 
//     keepSequences = o_keepSeqs;
//     keepHeader = o_keepHeader;
// 
//     forceCaps = false;
//     upperCase = false;
//     lowerCase = false;
// 
//     gapSymbol = "-";
// 
// 
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
//     /* Basic information for the new newAlignment */
//     sequences = new string[sequenNumber];
//     for(i = 0; i < sequenNumber; i++)
//         sequences[i] = o_sequences[i];
// 
//     seqsName = new string[sequenNumber];
//     for(i = 0; i < sequenNumber; i++)
//         seqsName[i] = o_seqsName[i];
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
// 
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
//     residuesNumber = new int[sequenNumber];
//     if((isAligned) || (o_residuesNumber != NULL)) {
//         for(i = 0; i < sequenNumber; i++)
//             residuesNumber[i] = residNumber;
//     } else {
//         for(i = 0; i < sequenNumber; i++)
//             residuesNumber[i] = o_residuesNumber[i];
//     }
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
// 
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
//     if(o_seqsInfo != NULL) {
//         seqsInfo = new string[sequenNumber];
//         for(i = 0; i < sequenNumber; i++)
//             seqsInfo[i] = o_seqsInfo[i];
//     } else seqsInfo = NULL;
// 
//     saveResidues  = NULL;
//     if(o_saveResidues != NULL) {
//         saveResidues = new int[residNumber];
//         for(i = 0, j = 0; i < OldResidues; i++)
//             if(o_saveResidues[i] != -1) {
//                 saveResidues[j] = o_saveResidues[i];
//                 j++;
//             }
//     }
// 
//     saveSequences = NULL;
//     if(o_saveSequences != NULL) {
//         saveSequences = new int[sequenNumber];
//         for(i = 0, j = 0; i < OldSequences; i++)
//             if(o_saveSequences[i] != -1) {
//                 saveSequences[j] = o_saveSequences[i];
//                 j++;
//             }
//     }
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
// 
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
//     identities = NULL;
//     if(o_identities != NULL) {
//         identities = new float*[sequenNumber];
//         for(i = 0, j = 0; i < OldSequences; i++) {
//             if(o_saveSequences[i] != -1) {
//                 identities[j] = new float[sequenNumber];
//                 for(k = 0, ll = 0; k < OldSequences; k++) {
//                     if(o_saveSequences[k] != -1) {
//                         identities[j][ll] = o_identities[i][k];
//                         ll++;
//                     }
//                 }
//                 j++;
//             }
//         }
//     }
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
// 
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
//     /* Any structure associated to the new newAlignment is
//      * initialize to NULL. In this way, these structure,
//      * if it will be necessary, has to be computed */
//     sgaps  =     NULL;
//     scons  =     NULL;
//     SequencesMatrix =  NULL;
//     identities = NULL;
//     /* ***** ***** ***** ***** ***** ***** ***** ***** */
// }

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Overlapping operator = to use it as a kind of class constructor */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

newAlignment::newAlignment(newAlignment& old) {

    if(this != &old) {

        int i, j;

        sequenNumber = old.sequenNumber;
        residNumber =  old.residNumber;

        isAligned =  old.isAligned;
        reverse   =  old.reverse;

        dataType = old.dataType;

        sequences = new string[sequenNumber];
        for(i = 0; i < sequenNumber; i++)
            sequences[i] = old.sequences[i];

        seqsName = new string[sequenNumber];
        for(i = 0; i < sequenNumber; i++)
            seqsName[i] = old.seqsName[i];
        /* ***** ***** ***** ***** ***** ***** ***** ***** */

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        //delete [] residuesNumber;
        if(old.residuesNumber) {
            residuesNumber = new int[sequenNumber];
            for(i = 0; i < sequenNumber; i++)
                residuesNumber[i] = old.residuesNumber[i];
        }
        else residuesNumber = NULL;
        /* ***** ***** ***** ***** ***** ***** ***** ***** */

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        //delete [] seqsInfo;
        if(old.seqsInfo) {
            seqsInfo = new string[sequenNumber];
            for(i = 0; i < sequenNumber; i++)
                seqsInfo[i] = old.seqsInfo[i];
        } else seqsInfo = NULL;

        //delete [] saveResidues;
        if(old.saveResidues) {
            saveResidues = new int[residNumber];
            for(i = 0; i < residNumber; i++)
                saveResidues[i] = old.saveResidues[i];
        } else saveResidues = NULL;

        //delete [] saveSequences;
        if(old.saveSequences) {
            saveSequences = new int[sequenNumber];
            for(i = 0; i < sequenNumber; i++)
                saveSequences[i] = old.saveSequences[i];
        } else saveSequences = NULL;
        /* ***** ***** ***** ***** ***** ***** ***** ***** */

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        //delete [] identities;
        if(old.identities) {
            identities = new float*[sequenNumber];
            for(i = 0; i < sequenNumber; i++) {
                identities[i] = new float[sequenNumber];
                for(j = 0; j < sequenNumber; j++)
                    identities[i][j] = old.identities[i][j];
            }
        } else identities = NULL;
        /* ***** ***** ***** ***** ***** ***** ***** ***** */

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        //delete sgaps;
        sgaps = NULL;

        //delete scons;
        scons = NULL;

        //delete SequencesMatrix;
        SequencesMatrix = NULL;
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    this -> Cleaning = new Cleaner(this, old.Cleaning);
    
    this -> Statistics = new StatisticsManager(this, old.Statistics);
    
    this -> ReadWrite = new ReadWriteManager(this, old.ReadWrite);

}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Class destructor */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

newAlignment::~newAlignment(void) {
    int i;

    if(sequences != NULL)
        delete [] sequences;
    sequences = NULL;

    if(seqsName != NULL)
        delete [] seqsName;
    seqsName = NULL;

    if(residuesNumber != NULL)
        delete [] residuesNumber;
    residuesNumber = NULL;

    if(seqsInfo != NULL)
        delete [] seqsInfo;
    seqsInfo = NULL;

    if(saveResidues != NULL)
        delete[] saveResidues;
    saveResidues = NULL;

    if(saveSequences != NULL)
        delete[] saveSequences;
    saveSequences = NULL;

    if(identities != NULL) {
        for(i = 0; i < sequenNumber; i++)
            delete [] identities[i];
        delete [] identities;
    }
    identities = NULL;

    if(sgaps != NULL)
        delete sgaps;
    sgaps = NULL;

    if(scons != NULL)
        delete scons;
    scons = NULL;

    if(SequencesMatrix != NULL)
        delete SequencesMatrix;
    SequencesMatrix = NULL;

    sequenNumber = 0;
    residNumber  = 0;

    isAligned = false;
    reverse   = false;

    dataType = 0;
    
    delete Cleaning;
    delete Statistics;
    delete ReadWrite;
}



newAlignment *newAlignment::getTranslationCDS(int newResidues, int newSequences, int *ColumnsToKeep, string *oldSeqsName, sequencesMatrix *seqMatrix, 
                                              newAlignment *ProtAlig) {

    string *matrixAux;
    newAlignment *newAlig;
    int i, j, k, l, oldResidues;
    int *mappedSeqs, *tmpSequence, *selectedRes;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Map the selected protein sequences to the input
     * coding sequences */
    mappedSeqs = new int[newSequences];
    for(i = 0; i < sequenNumber; i++)
        for(j = 0; j < newSequences; j++)
            if(!seqsName[i].compare(oldSeqsName[j])) {
                mappedSeqs[j] = i;
                break;
            }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Get the original newAlignment size as well the original
     * selected/non-selected columns */
    oldResidues = seqMatrix -> getResidNumber();
    selectedRes = new int[oldResidues];

    for(i = 0; i < oldResidues; i++)
        selectedRes[i] = i;

    for(j = 0; j < ColumnsToKeep[0]; j++)
        selectedRes[j] = -1;

    for(i = 0; i < newResidues - 1; i++)
        for(j = ColumnsToKeep[i] + 1; j < ColumnsToKeep[i+1]; j++)
            selectedRes[j] = -1;

    for(j = ColumnsToKeep[newResidues - 1] + 1; j < oldResidues; j++)
        selectedRes[j] = -1;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* We allocate memory to save the columns selected */
    matrixAux = new string[newSequences];
    tmpSequence = new int[oldResidues];

    /* Using the information about which residues for each
     * sequence was selected by others function, we process
     * these residues to recover the corresponding codons */
    for(i = 0; i < newSequences; i++)
        if(seqMatrix -> getSequence(oldSeqsName[i], tmpSequence)) {
            for(j = 0; j < oldResidues; j++) {
                if((selectedRes[j] != -1) && (tmpSequence[j] != 0)) {

                    for(k = 3 * (tmpSequence[j] - 1), l = 0; l < 3; k++, l++) {
                        /* Check whether the nucleotide sequences end has been reached or not.
                         * If it has been reached, complete backtranslation using indetermination
                         * symbols 'N' */
                        if((int) sequences[mappedSeqs[i]].length() > k)
                            matrixAux[i].resize(matrixAux[i].size() + 1, sequences[mappedSeqs[i]][k]);
                        else
                            matrixAux[i].resize(matrixAux[i].size() + 1, 'N');
                    }
                } else if(selectedRes[j] != -1) {
                    matrixAux[i].resize(matrixAux[i].size() + 3, '-');
                }
            }
            /* If there is any problems with a sequence then
             * the function returns an error */
        } else {
            return NULL;
        }
    
    newAlig = new newAlignment(*this);
    newAlig -> ReadWrite -> aligInfo = "";
    newAlig -> seqsInfo = NULL;
    
    newAlig -> sequenNumber =   newSequences;
    newAlig -> residNumber =    newResidues * 3;
    
    newAlig -> sequences =  new string[newSequences];
    newAlig -> seqsName =   new string[newSequences];
    
    newAlig -> ReadWrite -> iformat =    ProtAlig -> ReadWrite->getInputFormat();
    newAlig -> ReadWrite -> oformat =    ProtAlig -> ReadWrite->getOutputFormat();
    newAlig -> ReadWrite -> shortNames = ProtAlig -> getShortNames();
    
    newAlig -> dataType =   DNAType;
    newAlig -> isAligned =  true;
    
    newAlig -> reverse =  ProtAlig -> getReverseFlag();
//     newAlig -> OldResidues = oldResidues * 3; TODO?
    newAlig -> residuesNumber = NULL;
    newAlig -> saveSequences = NULL;
    newAlig -> saveResidues = NULL;
    
//     newAlig -> ghWindow = 0;
//     newAlig -> shWindow = 0;
    
    newAlig -> Cleaning -> blockSize = ProtAlig -> getBlockSize();
    
    newAlig -> identities = NULL;
    
    newAlig -> saveResidues = new int[residNumber];

    
    
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Deallocated auxiliar memory */
    delete [] matrixAux;
    delete mappedSeqs;
    delete tmpSequence;
    delete selectedRes;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Return the new newAlignment reference */
    return newAlig;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}



/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Return the number of the sequences in the newAlignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int newAlignment::getNumSpecies(void) {
    return sequenNumber;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Return the number of residues in the newAlignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int newAlignment::getNumAminos(void) {
    return residNumber;
}


/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method stores the diferent windows values in the newAlignment object */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void newAlignment::setWindowsSize(int ghWindow_, int shWindow_) {
    Statistics -> ghWindow = ghWindow_;
    Statistics -> shWindow = shWindow_;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method lets to change the output format */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void newAlignment::setOutputFormat(int format_, bool shortNames_) {
    ReadWrite -> oformat    = format_;
    ReadWrite -> shortNames = shortNames_;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method sets a new block size value */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void newAlignment::setBlockSize(int blockSize_) {
    Cleaning -> blockSize = blockSize_;
}


/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Return true if the ReadWrite -> shortNames flag has been setted */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int newAlignment::getShortNames(void) {
    return ReadWrite -> shortNames;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Return true if the reverse flag has been setted. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int newAlignment::getReverseFlag(void) {
    return reverse;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method lets to change the output newAlignment orientation */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void newAlignment::setReverseFlag(void) {
    reverse = true;
}


/* Set appropiate flag to decide whether sequences composed only by gaps should
 * be kept or not */
void newAlignment::setKeepSequencesFlag(bool flag) {
    Cleaning -> keepSequences = flag;
}

/* Set appropiate flag to decide whether original sequences header should be
 * dumped into the output stream with or without any preprocessing step */
void newAlignment::setKeepSeqsHeaderFlag(bool flag) {
    ReadWrite -> keepHeader = flag;
}
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Return the block size value */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int newAlignment::getBlockSize(void) {
    return Cleaning -> blockSize;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method returns the sequences name */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void newAlignment::getSequences(string *Names) {
    for(int i = 0; i < sequenNumber; i++)
        Names[i] = seqsName[i];
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void newAlignment::getSequences(string *Names, int *lengths) {
    for(int i = 0; i < sequenNumber; i++) {
        lengths[i] = (int) utils::removeCharacter('-', sequences[i]).length();
        Names[i] = seqsName[i];
    }
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void newAlignment::getSequences(string *Names, string *Sequences, int *Lengths) {
    for(int i = 0; i < sequenNumber; i++) {
        Names[i] = seqsName[i];
        Sequences[i] = utils::removeCharacter('-', sequences[i]);
        Lengths[i] = (int) Sequences[i].length();
    }
}


/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* For a given set of sequences name, check for its correspondence with
 * the current sequences name. If there is not a complete correspondence
 * between both sets, the method return false, otherwise, return true */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool newAlignment::getSequenceNameOrder(string *names, int *orderVector) {
    int i, j, numNames;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* For each name in the current newAlignment, we look
     * for its correspondence in the input set */
    for(i = 0, numNames = 0; i < sequenNumber; i++) {
        for(j = 0; j < sequenNumber; j++) {
            if(seqsName[i] == names[j]) {
                orderVector[i] = j;
                numNames++;
                break;
            }
        }
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Depending on if we get the complete correspondence
     * between both sets of names, we return true or not */
    if(numNames == sequenNumber) return true;
    else return false;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
///* This method lets to build a sequence matrix. A sequence matrix contains
// * the residue position for each sequence without taking into account the
// * gaps in the sequence. This means that at position 10 we have the residue
// * 1 for that sequence */
///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
//void newAlignment::sequenMatrix(void) {
//  if(SequencesMatrix == NULL)
//    SequencesMatrix = new sequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
//}
//
///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
///* Use this method to destroy a given sequence matrix */
///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
//void newAlignment::destroySequenMatrix(void) {
//  if(SequencesMatrix != NULL)
//    delete SequencesMatrix;
//  SequencesMatrix = NULL;
//}
//
///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
///* Print the newAlignment's sequence matrix */
///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
//void newAlignment::printSequenMatrix(void) {
//  if(SequencesMatrix == NULL)
//    SequencesMatrix = new sequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
//  SequencesMatrix -> printMatrix();
//}
//
///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
///* Given an index, this method returns the sequence matrix column for that
// * index */
///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
//void newAlignment::getColumnSeqMatrix(int column, int *columnSeqMatrix) {
//  if(SequencesMatrix == NULL)
//    SequencesMatrix = new sequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
//  SequencesMatrix -> getColumn(column, columnSeqMatrix);
//}
//
///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
///* This function looks if a given value belongs a given row. In the affirmative
// * case, returns the columns value for that row and value combination */
///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
//void newAlignment::getColumnSeqMatrix(int value, int row, int *columnSeqMatrix) {
//  if(SequencesMatrix == NULL)
//    SequencesMatrix = new sequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
//  SequencesMatrix -> getColumn(value, row, columnSeqMatrix);
//}
//
///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
///* This function change the rows order in the sequence matrix */
///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
//void newAlignment::setSeqMatrixOrder(int *order) {
//  if(SequencesMatrix == NULL)
//    SequencesMatrix = new sequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
//  SequencesMatrix -> setOrder(order);
//}
//
///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
///* This function change the rows order in the sequence matrix */
///* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
//sequencesMatrix *newAlignment::getSeqMatrix(void) {
//  if(SequencesMatrix == NULL)
//    SequencesMatrix = new sequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
//  return SequencesMatrix;
//}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Computes, if it's necessary, and return the newAlignment's type */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int newAlignment::getAlignmentType(void) {
    if(dataType == 0)
        dataType = utils::checkAlignmentType(sequenNumber, residNumber, sequences);
    return dataType;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Returns the correspondence between the columns in the original and in the
 * trimmed newAlignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int *newAlignment::getCorrespResidues(void) {
    return saveResidues;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Returns the correspondence between the sequences in the original and in
 * the trimmed newAlignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int *newAlignment::getCorrespSequences(void) {
    return saveSequences;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Returns if the newAlignment is aligned or not */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool newAlignment::isFileAligned(void) {
    return isAligned;
}

/* *****************************************************************************
 * *****************************************************************************
 * *****************************************************************************
 * ************************************************************************** */

/* Function for computing the complementary newAlignment. It just turn around the
 * current columns/sequences selection */



/* Function for copying to previously allocated memory those data selected
 * for being in the final newAlignment */
void newAlignment::fillNewDataStructure(string *newMatrix, string *newNames) {
    int i, j, k;

    /* Copy only those sequences/columns selected */
    for(i = 0, j = 0; i < sequenNumber; i++) {
        if(saveSequences[i] == -1)
            continue;

        newNames[j] = seqsName[i];
        for(k = 0; k < residNumber; k++) {
            if(saveResidues[k] == -1)
                continue;
            newMatrix[j].resize(newMatrix[j].size() + 1, sequences[i][k]);
        }
        j++;
    }
}

/* Function for copying to previously allocated memory those data selected
 * for being in the final newAlignment */
void newAlignment::fillNewDataStructure(newValues *data) {
    int i, j, k;

    /* Copy only those sequences/columns selected */
    for(i = 0, j = 0; i < sequenNumber; i++) {
        if(saveSequences[i] == -1)
            continue;

        data -> seqsName[j] = seqsName[i];
        for(k = 0; k < residNumber; k++) {
            if(saveResidues[k] == -1)
                continue;
            data -> matrix[j].resize(data -> matrix[j].size() + 1, sequences[i][k]);
        }
        j++;
    }
    // cerr << data -> seqsName[j-1] << endl;
}

/* Check if CDS file is correct based on: Residues are DNA/RNA (at most). There
 * is not gaps in the whole dataset. Each sequence is multiple of 3. At the same
 * time, the function will remove stop codons if appropiate flags are used */
bool newAlignment::prepareCodingSequence(bool splitByStopCodon, bool ignStopCodon,\
        newAlignment *proteinAlig) {

    bool warning = false;
    size_t found;
    int i;

    /* New code: We now care about the presence of wildcards/indeterminations
     * characters such as 'X' or 'B' into protein sequences as well as about the
     * presence of Selenocysteines ('U') or Pyrrolysines ('O'). It only works, by
     * now, with the universal genetic code */
    int *protSeqsLengths, numbProtSeqs, current_prot;
    string *protSeqsNames, *protSequences;
    char aminoAcid;

    numbProtSeqs =  proteinAlig -> getNumSpecies();
    protSeqsNames = new string[numbProtSeqs];
    protSequences = new string[numbProtSeqs];
    protSeqsLengths = new int[numbProtSeqs];

    proteinAlig -> getSequences(protSeqsNames, protSequences, protSeqsLengths);

    /* Check read sequences are real DNA/RNA */
    if (getAlignmentType() == AAType) {
        cerr << endl << "ERROR: Check input CDS file. It seems to content protein "
             << "residues." << endl << endl;
        return false;
    }

    for(i = 0; i < sequenNumber; i++) {

        /* Get protein sequence to compare against any potential stop codon in the
         * coding sequence. If there is not protein sequence for current coding
         * sequence, skip its analysis */
        for(current_prot = 0; current_prot < numbProtSeqs; current_prot++)
            if(protSeqsNames[current_prot] == seqsName[i])
                break;
        if(current_prot == numbProtSeqs)
            continue;

        if(sequences[i].find("-") != string::npos) {
            if (!warning)
                cerr << endl;
            cerr << "ERROR: Sequence \"" << seqsName[i] << "\" has, at least, one gap"
                 << endl << endl;
            return false;
        }

        if((sequences[i].length() % 3) != 0) {
            if (!warning)
                cerr << endl;
            warning = true;
            cerr << "WARNING: Sequence length \"" << seqsName[i] << "\" is not "
                 << "multiple of 3 (length: " << sequences[i].length() << ")" << endl;
        }

        /* Ignore stop codons from the CDS if set by the user */
        if (ignStopCodon)
            continue;

        /* Detect universal stop codons in the CDS. Then, compare those stop codons
         * against the protein sequence to see whether they are real stop codons or
         * are representing rare amino-acids such as Selenocysteines or Pyrrolysines
         * It also allows stop-codons when there are wildcards/indet characters in
         * the protein sequence. CDS sequences could be splitted using stop codons
         * from the sequence itself */

        /* Initialize first appearence of a given stop codon to -1.
         * That means that it has not been found yet */
        found = -1;
        do {
            found = sequences[i].find("TGA", found + 1);

            /* If a stop codon has been found and its position is multiple of 3.
             * Analize it */
            if((found != string::npos) && (((int) found % 3) == 0)) {

                aminoAcid = (char) toupper(protSequences[current_prot][(int) found/3]);
                /* It may be a Selenocysteine ('TGA') which should be represented as 'U'
                 * or wildcard/indet characters such as 'X' or 'B' */
                //~ if ((aminoAcid == 'U') || (aminoAcid == 'X') || (aminoAcid == 'B'))
                /* If a rare amino-acids such as 'U'/'O' or a wildcard/indet character
                 * such as 'B'/'X' is present, skip current stop codon */
                if((aminoAcid == 'U') || (aminoAcid == 'O') || (aminoAcid == 'X') || \
                        (aminoAcid == 'B'))
                    continue;

                /* If split_by_stop_codon flag is activated then cut input CDS sequence
                 * up to first appearance of a stop codon */
                else if(splitByStopCodon) {
                    if (!warning)
                        cerr << endl;
                    warning = true;
                    cerr << "WARNING: Cutting sequence \"" << seqsName[i] << "\" at first"
                         << " appearance of stop codon \"TGA\" (residue \"" << aminoAcid
                         << "\") at position " << (int) found + 1 << " (length: "
                         << sequences[i].length() << ")" << endl;
                    sequences[i].resize((int) found);
                }
                /* Otherwise, warn about it and return an error */
                else {
                    if (!warning)
                        cerr << endl;
                    cerr << "ERROR: Sequence \"" << seqsName[i] << "\" has stop codon \""
                         << "TGA\" (residue \"" << aminoAcid << "\") at position "
                         << (int) found + 1 << " (length: " << sequences[i].length() << ")"
                         << endl << endl;
                    return false;
                }
            }
            /* Iterate over the CDS until not stop codon is found */
        } while(found != string::npos);

        /* Initialize first appearence of a given stop codon to -1.
         * That means that it has not been found yet */
        found = -1;
        do {
            found = sequences[i].find("TAA", found + 1);

            /* If a stop codon has been found and its position is multiple of 3.
             * Analize it */
            if((found != string::npos) && (((int) found % 3) == 0)) {

                aminoAcid = (char) toupper(protSequences[current_prot][(int) found/3]);
                /* Check if there is any wildcard/indet characters such as 'X' or 'B' */
                //~ if ((aminoAcid == 'X') || (aminoAcid == 'B'))
                /* If a rare amino-acids such as 'U'/'O' or a wildcard/indet character
                 * such as 'B'/'X' is present, skip current stop codon */
                if((aminoAcid == 'U') || (aminoAcid == 'O') || (aminoAcid == 'X') || \
                        (aminoAcid == 'B'))
                    continue;

                /* If split_by_stop_codon flag is activated then cut input CDS sequence
                 * up to first appearance of a stop codon */
                else if(splitByStopCodon) {
                    if (!warning)
                        cerr << endl;
                    warning = true;
                    cerr << "WARNING: Cutting sequence \"" << seqsName[i] << "\" at first"
                         << " appearance of stop codon \"TAA\" (residue \"" << aminoAcid
                         << "\") at position " << (int) found + 1 << " (length: "
                         << sequences[i].length() << ")" << endl;
                    sequences[i].resize((int) found);
                }
                /* Otherwise, warn about it and return an error */
                else {
                    if (!warning)
                        cerr << endl;
                    cerr << "ERROR: Sequence \"" << seqsName[i] << "\" has stop codon \""
                         << "TAA\" (residue \"" << aminoAcid << "\") at position "
                         << (int) found + 1 << " (length: " << sequences[i].length() << ")"
                         << endl << endl;
                    return false;
                }
            }
            /* Iterate over the CDS until not stop codon is found */
        } while(found != string::npos);

        /* Initialize first appearence of a given stop codon to -1.
         * That means that it has not been found yet */
        found = -1;
        do {
            found = sequences[i].find("TAG", found + 1);
            /* If a stop codon has been found and its position is multiple of 3.
             * Analize it */
            if((found != string::npos) && (((int) found % 3) == 0)) {

                aminoAcid = (char) toupper(protSequences[current_prot][(int) found/3]);
                /* It may be a Pyrrolysine ('TAG') which should be represented as 'O'
                 * or wildcard/indet characters such as 'X' or 'B' */
                //~ if ((aminoAcid == 'O') || (aminoAcid == 'X') || (aminoAcid == 'B'))
                /* If a rare amino-acids such as 'U'/'O' or a wildcard/indet character
                 * such as 'B'/'X' is present, skip current stop codon */
                if((aminoAcid == 'U') || (aminoAcid == 'O') || (aminoAcid == 'X') || \
                        (aminoAcid == 'B'))
                    continue;

                /* If split_by_stop_codon flag is activated then cut input CDS sequence
                 * up to first appearance of a stop codon */
                else if(splitByStopCodon) {
                    if (!warning)
                        cerr << endl;
                    warning = true;
                    cerr << "WARNING: Cutting sequence \"" << seqsName[i] << "\" at first"
                         << " appearance of stop codon \"TAG\" (residue \"" << aminoAcid
                         << "\") at position " << (int) found + 1 << " (length: "
                         << sequences[i].length() << ")" << endl;
                    sequences[i].resize((int) found);
                }
                /* Otherwise, warn about it and return an error */
                else {
                    if (!warning)
                        cerr << endl;
                    cerr << "ERROR: Sequence \"" << seqsName[i] << "\" has stop codon \""
                         << "TAG\" (residue \"" << aminoAcid << "\") at position "
                         << (int) found + 1 << " (length: " << sequences[i].length() << ")"
                         << endl << endl;
                    return false;
                }
            }
            /* Iterate over the CDS until not stop codon is found */
        } while(found != string::npos);
    }

    /* If everything was return an OK to informat about it. */
    return true;
}

/* Function designed to check whether input CDS file is correct or not based on
 * some features: Sequences are in both files (it could be more on CDS file),
 * they have (more or less) same ength. Otherwise, some nucleotides could be
 * excluded or some 'N's added to fit protein length. */
bool newAlignment::checkCorrespondence(string *names, int *lengths, int \
                                       totalInputSeqs, int multiple = 1) {

    int i, j, seqLength, indet;
    bool warnings = false;
    string tmp;

    /* For each sequence in the current protein newAlignment, look for its coding
     * DNA sequence checking that they have the same size. */
    for(i = 0; i < sequenNumber; i++) {

        /* Get protein sequence length removing any possible gap. Get as well last
         * residue from current sequence */

        tmp = utils::removeCharacter('-', sequences[i]);
        seqLength = tmp.length() * multiple;
        indet = ((int) tmp.length() - utils::min((int) tmp.find_last_not_of("X"), \
                 (int) tmp.find_last_not_of("x"))) - 1;

        /* Go through all available CDS looking for the one with the same ID */
        for(j = 0; j < totalInputSeqs; j++) {

            /* Once both ID matchs, compare its lengths */
            if(seqsName[i] == names[j]) {

                /* If both sequences have the same length, stop the search */
                if(seqLength == lengths[j])
                    break;

                /* If nucleotide sequence is larger than protein sequence, warn about
                 * it and continue the verification process. It will used the 'Nth'
                 * first nucleotides for the conversion */
                else if(seqLength < lengths[j]) {
                    if (!warnings)
                        cerr << endl;
                    warnings = true;
                    cerr << "WARNING: Sequence \"" << seqsName[i] << "\" will be cutted "
                         << "at position " << seqLength << " (length: "<< lengths[j] << ")"
                         << endl;
                    break;
                }

                /* It has been detected some indeterminations at the end of the protein
                 * sequence. That issue could be cause by some incomplete codons in the
                 * nucleotide sequences. This issue is solved adding as much 'N' symbols
                 * as it is needed to preserve the backtranslated newAlignment */
                else if((indet > 0) && (indet > (seqLength - lengths[j])/3)) {
                    if (!warnings)
                        cerr << endl;
                    warnings = true;
                    cerr << "WARNING: Sequence \"" << seqsName[i] << "\" has some inde"
                         << "termination symbols 'X' at the end of sequence. They will be"
                         << " included in the final newAlignment." << endl;
                    break;
                }

                /* If nucleotide sequence is shorter than protein sequence, return an
                 * error since it is not feasible to cut the input protein aligment to
                 * fit it into CDNA sequences size */
                else {
                    if (!warnings)
                        cerr << endl;
                    warnings = true;
                    cerr << "WARNING: Sequence \"" << seqsName[i] << "\" has less nucleo"
                         << "tides (" << lengths[j] << ") than expected (" << seqLength
                         << "). It will be added N's to complete the sequence"  << endl;
                    break;
                }
            }
        }

        /* Warn about a mismatch a sequences name level */
        if(j == totalInputSeqs) {
            cerr << endl << "ERROR: Sequence \"" << seqsName[i] << "\" is not in "
                 << "CDS file." << endl << endl;
            return false;
        }
    }

    /* If everything is OK, return an appropiate flag */
    return true;
}

bool newAlignment::fillMatrices(bool aligned) {
    /* Function to determine if a set of sequences, that can be aligned or not,
     * have been correctly load and are free of errors. */
    int i, j;

    /* Initialize some variables */
    residuesNumber = new int[sequenNumber];
    for(i = 0; i < sequenNumber; i++) {
        residuesNumber[i] = sequences[i].size();
    }

    /* Check whether there are any unknow/no allowed character in the sequences */
    for(i = 0; i < sequenNumber; i++)
        for(j = 0; j < residuesNumber[i]; j++)
            if((!isalpha(sequences[i][j])) && (!ispunct(sequences[i][j]))) {
                cerr << endl << "ERROR: The sequence \"" << seqsName[i] << "\" has an "
                     << "unknown (" << sequences[i][j] << ") character." << endl;
                return false;
            }

    /* Check whether all sequences have same size or not */
    for(i = 1; i < sequenNumber; i++)
        if(residuesNumber[i] != residuesNumber[i-1])
            break;
    /* Set an appropriate flag for indicating if sequences are aligned or not */
    isAligned = (i != sequenNumber) ? false : true;

    /* Warm about those cases where sequences should be aligned
     * and there are not */
    if (aligned and !isAligned) {
        cerr << endl << "ERROR: Sequences should be aligned (all with same length) "
             << "and there are not. Check your input alignment" << endl;
        return false;
    }

    /* Full-fill some information about input alignment */
    if(residNumber == 0)
        residNumber = residuesNumber[0];

    /* Check whether aligned sequences have the length fixed for the input alig */
    for(i = 0; (i < sequenNumber) and (aligned); i++) {
        if(residuesNumber[i] != residNumber) {
            cerr << endl << "ERROR: The sequence \"" << seqsName[i] << "\" ("
                 << residuesNumber[i] << ") does not have the same number of residues "
                 << "fixed by the alignment (" << residNumber << ")." << endl;
            return false;
        }
    }

    /* If the sequences are aligned, initialize some additional variables.
     * These variables will be useful for posterior analysis */
    if((aligned) || (isAligned)) {

        /* Asign its position to each column. That will be used to determine which
         * columns should be kept in output alignment after applying any method
         * and which columns should not */
        saveResidues = new int[residNumber];
        for(i = 0; i < residNumber; i++)
            saveResidues[i] = i;

        /* Asign its position to each sequence. Similar to the columns numbering
         * process, assign to each sequence its position is useful to know which
         * sequences will be in the output alignment */
        saveSequences = new int[sequenNumber];
        for(i = 0; i < sequenNumber; i++)
            saveSequences[i] = i;
    }

    /* Return an flag indicating that everything is fine */
    return true;
}



void newAlignment::printAlignmentInfo(ostream &file) {
    /* Print information about sequences number, average sequence length, maximum
     * and minimum sequences length, etc */

    int i, j, valid_res, max, min, max_pos, min_pos, total_res;

    /* Storage which sequences are the longest and shortest ones */
    max = 0;
    max_pos = 0;
    min_pos = 0;
    min = residuesNumber[0];

    for(i = 0, total_res = 0; i < sequenNumber; i++) {

        /* Discard gaps from current sequence and then compute real length */
        for(j = 0, valid_res = 0; j < residuesNumber[i]; j++)
            valid_res += (sequences[i][j] != '-' ? 1 : 0);

        /* Compute the total residues in the alignment to calculate avg. sequence
         * length */
        total_res += valid_res;

        /* Get values for the longest sequence */
        max_pos = (max > valid_res) ? max_pos : i;
        max = (max > valid_res) ? max : valid_res;
        /* Similarily, get values for the shortest sequence */
        min_pos = (min < valid_res) ? min_pos : i;
        min = (min < valid_res) ? min : valid_res;
    }

    file << "## Total sequences\t" << sequenNumber << endl;
    if (isFileAligned())
        file << "## Alignment length\t" << residNumber << endl;
    file  << "## Avg. sequence length\t" << (float) total_res / sequenNumber << endl
          << "## Longest seq. name\t'"   << seqsName[max_pos] << "'" << endl
          << "## Longest seq. length\t"  << max << endl
          << "## Shortest seq. name\t'"  << seqsName[min_pos] << "'" << endl
          << "## Shortest seq. length\t" << min << endl;
}






/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function computes some parameters from the input alignment such as
 * identity average, identity average from each sequence and its most similar
 * one, etc, to select which one is the best automated method to trim this
 * alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */


/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method prints different identity values computed from the alignment.
 * In this method, we asses the identity values matrix as well as diferent
 * average values. Moreover, the method computes which one is the most
 * similar sequence, in term of identity values, for each one in this alig */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void newAlignment::printSeqIdentity(void) {

    int i, j, k, pos, maxLongName;
    float mx, avg, maxAvgSeq = 0, maxSeq = 0, avgSeq = 0, **maxs;

    /* Ask for the sequence identities assesment */
    if(identities == NULL)
        Cleaning -> calculateSeqIdentity();

    /* For each sequence, we look for its most similar one */
    maxs = new float*[sequenNumber];

    for(i = 0; i < sequenNumber; i++) {
        maxs[i] = new float[2];

        /* Get the most similar sequence to the current one in term of identity */
        for(k = 0, mx = 0, avg = 0, pos = i; k < sequenNumber; k++) {
            if(i != k) {
                avg += identities[i][k];
                if(mx < identities[i][k]) {
                    mx = identities[i][k];
                    pos = k;
                }
            }
        }
        /* Update global average variables*/
        avgSeq += avg/(sequenNumber - 1);
        maxAvgSeq += mx;

        /* Save the maximum average identity value for each sequence */
        maxs[i][0] = mx;
        maxs[i][1] = pos;
    }

    /* Compute general averages */
    avgSeq = avgSeq/sequenNumber;
    maxAvgSeq = maxAvgSeq/sequenNumber;

    /* Compute longest sequences name */
    for(i = 0, maxLongName = 0; i < sequenNumber; i++)
        maxLongName = utils::max(maxLongName, seqsName[i].size());

    /* Once the method has computed all of different values, it prints it */
    cout.precision(4);
    cout << fixed;

    for(i = 0, maxSeq = 0; i < sequenNumber; i++)
        if(maxs[i][0] > maxSeq)
            maxSeq = maxs[i][0];

    cout << endl << "## MaxIdentity\t" << maxSeq;
    cout << endl << "#> MaxIdentity\tGet the maximum identity value for any pair "
         << "of sequences in the alignment" << endl;

    cout << endl << "## AverageIdentity\t" << avgSeq;
    cout << endl << "#> AverageIdentity\tAverage identity between all sequences";

    cout << endl << endl << "## Identity sequences matrix";
    for(i = 0; i < sequenNumber; i++) {
        cout << endl << setw(maxLongName + 2) << left << seqsName[i] << "\t";
        for(j = 0; j < i; j++)
            cout << setiosflags(ios::left) << setw(10) << identities[i][j] << "\t";
        cout << setiosflags(ios::left) << setw(10) << 1.00 << "\t";
        for(j = i + 1; j < sequenNumber; j++)
            cout << setiosflags(ios::left) << setw(10) << identities[i][j] << "\t";
    }
    cout << endl;

    cout << endl << "## AverageMostSimilarIdentity\t" << maxAvgSeq;
    cout << endl << "#> AverageMostSimilarIdentity\t Average identity between "
         << "most similar pair-wise sequences";

    cout << endl << endl << "## Identity for most similar pair-wise sequences "
         << "matrix" << endl;
    for(i = 0; i < sequenNumber; i++)
        cout << setw(maxLongName + 2) << left << seqsName[i]
             << "\t" << setiosflags(ios::left) << setw(5)
             << maxs[i][0] << "\t" << seqsName[(int) maxs[i][1]] << endl;
    cout << endl;
}

/* *** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *** */
/*                                                                           */
/*                             NEW CODE: feb/2012                            */
/*                                                                           */
/* *** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *** */
void newAlignment::calculateColIdentity(float *ColumnIdentities) {

    int i, j, counter, pos, max, columnLen;
    char letter, indet, gapSymbol;
    string column;

    /* Initialize some data for make computation more precise */
    indet = getAlignmentType() == AAType ? 'X' : 'N';
    gapSymbol = '-';

    /* Compute identity score for the most frequent residue, it can be as well
     * gaps and indeterminations, for each column */
    for(i = 0, max = 0; i < residNumber; i++, max = 0, column.clear()) {

        /* Get residues from each column in capital letters */
        for(j = 0; j < sequenNumber; j++)
            /* Discard gaps and indeterminations from calculations */
            if((toupper(sequences[j][i]) != indet) && (sequences[j][i] != gapSymbol))
                column += toupper(sequences[j][i]);
        columnLen = column.size();

        /* Count letter frequency. It only matter the frequency. Use some shorcuts
         * to speed-up the process */
        while (!column.empty()) {
            letter = column[0];
            counter = 0;
            pos = 0;

            do {
                counter += 1;
                column.erase(pos, 1);
                pos = column.find(letter, pos);
            } while(pos != (int) string::npos);

            /* Keep only the most frequent residue */
            if(counter > max)
                max = counter;
            /* If column size is smaller than the current max, stop the count */
            if((int) column.size() < max)
                break;
        }

        /* Store column identity values */
        if(columnLen != 0)
            ColumnIdentities[i] = float(max)/columnLen;
    }
}

void newAlignment::printColumnsIdentity_DescriptiveStats(void) {

    float *colIdentities, avg, std, max, min;
    int i, positions;

    /* Allocate local memory for the computation */
    colIdentities = new float[residNumber];

    utils::initlVect(colIdentities, residNumber, -1);
    calculateColIdentity(colIdentities);

    for(i = 0, max = 0, min = 1, avg = 0, positions = 0; i < residNumber; i++) {
        if(colIdentities[i] != -1) {
            /* Compute on-the-fly max and min scores. Store accumulative score */
            avg += colIdentities[i];
            max = (colIdentities[i] > max) ? colIdentities[i] : max;
            min = (colIdentities[i] < min) ? colIdentities[i] : min;
            /* Count how many columns have a value score */
            positions += 1;
        }
    }
    /* Compute average identity column score */
    avg /= positions;

    /* Compute standard desviation */
    for(i = 0, std = 0; i < residNumber; i++)
        if(colIdentities[i] != -1)
            std += pow((colIdentities[i] - avg), 2);
    std = sqrt(std/positions);

    /* Print general descriptive stats */
    cout << "#maxColIdentity\t" << max << endl;
    cout << "#minColIdentity\t" << min << endl;
    cout << "#avgColIdentity\t" << avg << endl;
    cout << "#stdColIdentity\t" << std << endl;
}




