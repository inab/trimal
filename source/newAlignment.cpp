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
// #include "../include/ReadWriteManager.h"
#include <errno.h>
#include <ctype.h>
#include <string>
//#include <utils.h>
#include "../include/values.h"
#include "../include/defines.h"
//extern int errno;

#include <sstream>
#include <map>
#include <vector>

using namespace std;


/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Clafile constructor */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

newAlignment::newAlignment(void) {

    Cleaning = new Cleaner(this);
    Statistics = new StatisticsManager(this);

    /* newAlignment parameter */
    sequenNumber = 0;
    residNumber =  0;

    /* Are the input sequences aligned? */
    isAligned = false;

    /* Sequence datatype: DNA, RNA or Protein */
    dataType = 0;

    /* Sequence residues number */
// // //     residuesNumber = NULL;

    /* Columns and sequences that have been previously selected */
    saveResidues  = NULL;
    saveSequences = NULL;

    /* Input sequences as well other information such as sequences name, etc */
    sequences = NULL;
    seqsName  = NULL;
    seqsInfo  = NULL;

    /* Information computed from newAlignment */
    sgaps =             NULL;
    scons =             NULL;
    SequencesMatrix =   NULL;
    identities =        NULL;

}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Overlapping operator = to use it as a kind of class constructor */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

newAlignment::newAlignment(newAlignment& originalAlignment) {

     if(this != &originalAlignment) {

        int i, j;
        
        aligInfo = originalAlignment.aligInfo;
        
        for ( i = 0, j = 0; i < originalAlignment.sequenNumber; i++)
            if (originalAlignment.saveSequences[i] != -1) j++;

        sequenNumber = j;
        
        for ( i = 0, j = 0; i < originalAlignment.residNumber; i++)
            if (originalAlignment.saveResidues[i] != -1) j++;
        
        residNumber =  j;
        
        isAligned =  originalAlignment.isAligned;
        
        dataType = originalAlignment.dataType;

        
        sequences = new string[sequenNumber];
        seqsName = new string[sequenNumber];
        saveSequences = new int[sequenNumber];
        for(i = 0, j = 0; i < originalAlignment.sequenNumber; i++)
            if (originalAlignment.saveSequences[i] != -1)
            {
                sequences[j] = originalAlignment.sequences[i];
                saveSequences[j] = originalAlignment.saveSequences[i];
                seqsName[j++] = originalAlignment.seqsName[i];
            }

        
        if(originalAlignment.seqsInfo) {
            seqsInfo = new string[sequenNumber];
            for(i = 0, j = 0; i < originalAlignment.sequenNumber; i++)
                if (originalAlignment.saveSequences[i] != -1)
                    seqsInfo[j++] = originalAlignment.seqsInfo[i];
        } else seqsInfo = NULL;

        
        if(originalAlignment.saveResidues) {
            saveResidues = new int[residNumber];
            for(i = 0, j = 0; i < originalAlignment.residNumber; i++)
            {
                if (originalAlignment.saveResidues[i] != -1)
                saveResidues[j++] = originalAlignment.saveResidues[i];
            }
        } else saveResidues = NULL;

        
//         if(originalAlignment.saveSequences) {
//             for(i = 0; i < sequenNumber; i++)
//             {
//                 saveSequences[i] = originalAlignment.saveSequences[i];
//             }
//         } else saveSequences = NULL;

        
//         if(originalAlignment.identities) {
//             identities = new float*[sequenNumber];
//             for(i = 0; i < sequenNumber; i++) {
//                 identities[i] = new float[sequenNumber];
//                 for(j = 0; j < sequenNumber; j++)
//                     identities[i][j] = originalAlignment.identities[i][j];
//             }
//         } else 
            identities = NULL;
        
        
//         this-> fillMatrices(false);

        //delete sgaps;
        sgaps = NULL;

        //delete scons;
        scons = NULL;

        //delete SequencesMatrix;
        SequencesMatrix = NULL;
        
        this -> Cleaning = new Cleaner(this, originalAlignment.Cleaning);
        
        this -> Statistics = new StatisticsManager(this, originalAlignment.Statistics);

        
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Clafile destructor */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

newAlignment::~newAlignment(void) {
    int i;
    
    if(sequences != NULL)
        delete [] sequences;
    sequences = NULL;

    if(seqsName != NULL)
        delete [] seqsName;
    seqsName = NULL;

//     if(residuesNumber != NULL)
//         delete [] residuesNumber;
//     residuesNumber = NULL;

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

    dataType = 0;
    
    delete Cleaning;
    delete Statistics;

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
     * sequence was selected by others function, we procefile
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
            delete[] matrixAux;
            delete[] tmpSequence;
            delete[] mappedSeqs;
            delete[] selectedRes;
            return NULL;
        }
    
    newAlig = new newAlignment(*this);
    newAlig -> aligInfo = "";
    newAlig -> seqsInfo = NULL;
    
    newAlig -> sequenNumber =   newSequences;
    newAlig -> residNumber =    newResidues * 3;
    
    newAlig -> sequences =  new string[newSequences];
    newAlig -> seqsName =   new string[newSequences];

    newAlig -> dataType =   SequenceTypes::DNA;
    newAlig -> isAligned =  true;
    
//     newAlig -> reverse =  ProtAlig -> getReverseFlag();
//     newAlig -> OldResidues = oldResidues * 3; TODO?
//     newAlig -> residuesNumber = NULL;
    newAlig -> saveSequences = NULL;
    newAlig -> saveResidues = NULL;
    
//     newAlig -> ghWindow = 0;
//     newAlig -> shWindow = 0;
    
    newAlig -> Cleaning -> blockSize = ProtAlig -> getBlockSize();
    
    newAlig -> identities = NULL;
    
    newAlig -> saveResidues = new int[residNumber];

    /* Deallocated auxiliar memory */
    delete [] matrixAux;
    delete[] mappedSeqs;
    delete[] tmpSequence;
    delete[] selectedRes;

    return newAlig;

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

void newAlignment::setBlockSize(int blockSize_) {
    Cleaning -> blockSize = blockSize_;
}

void newAlignment::setKeepSequencesFlag(bool flag) {
    Cleaning -> keepSequences = flag;
}

int newAlignment::getBlockSize(void) {
    return Cleaning -> blockSize;
}

void newAlignment::calculateSeqIdentity(void) {

  int i, j, k, hit, dst;
  char indet;

  /* Depending on alignment type, indetermination symbol will be one or other */
  indet = getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

  /* Create identities matrix to store identities scores */
  identities = new float*[sequenNumber];

  /* For each seq, compute its identity score against the others in the MSA */
  for(i = 0; i < sequenNumber; i++) {
    identities[i] = new float[sequenNumber];

    /* It's a symmetric matrix, copy values that have been already computed */
    for(j = 0; j < i; j++)
      identities[i][j] = identities[j][i];
    identities[i][i] = 0;

    /* Compute identity scores for the current sequence against the rest */
    for(j = i + 1; j < sequenNumber; j++) {
      for(k = 0, hit = 0, dst = 0; k < residNumber; k++) {
      /* If one of the two positions is a valid residue,
       * count it for the common length */
        if(((sequences[i][k] != indet) && (sequences[i][k] != '-')) ||
           ((sequences[j][k] != indet) && (sequences[j][k] != '-'))) {
          dst++;
          /* If both positions are the same, count a hit */
          if(sequences[i][k] == sequences[j][k])
            hit++;
        }
      }

      /* Identity score between two sequences is the ratio of identical residues
       * by the total length (common and no-common residues) among them */
      identities[i][j] = (float) hit/dst;
    }
  }
}

void newAlignment::calculateRelaxedSeqIdentity(void) {
  /* Raw approximation of sequence identity computation designed for reducing
   * comparisons for huge alignemnts */

  int i, j, k, hit;

  /* Create identities matrix to store identities scores */
  identities = new float*[sequenNumber];

  /* For each seq, compute its identity score against the others in the MSA */
  for(i = 0; i < sequenNumber; i++) {
    identities[i] = new float[sequenNumber];

    /* It's a symmetric matrix, copy values that have been already computed */
    for(j = 0; j < i; j++)
      identities[i][j] = identities[j][i];
    identities[i][i] = 0;

    /* Compute identity score between the selected sequence and the others */
    for(j = i + 1; j < sequenNumber; j++) {
      for(k = 0, hit = 0; k < residNumber; k++) {
        /* If both positions are the same, count a hit */
        if(sequences[i][k] == sequences[j][k])
          hit++;
      }
    /* Raw identity score is computed as the ratio of identical residues between
     * alignment length */
      identities[i][j] = (float) hit/residNumber;
    }
  }
}


void newAlignment::calculateSeqOverlap(void) {
  /* Compute the overlap between sequences taken each of them as the reference
   * to compute such scores. It will lead to a non-symmetric matrix. */

  int i, j, k, shared, referenceLength;
  char indet;

  /* Depending on alignment type, indetermination symbol will be one or other */
  indet = getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

  /* Create overlap matrix to store overlap scores */
  overlaps = new float*[sequenNumber];

  /* For each seq, compute its overlap score against the others in the MSA */
  for(i = 0; i < sequenNumber; i++) {
    overlaps[i] = new float[sequenNumber];

    for(j = 0; j < sequenNumber; j++) {
      for(k = 0, shared = 0, referenceLength = 0; k < residNumber; k++) {
        /* If there a valid residue for the reference sequence, then see if
         * there is a valid residue for the other sequence. */
        if((sequences[i][k] != indet) && (sequences[i][k] != '-')) {
          referenceLength++;
          if ((sequences[j][k] != indet) && (sequences[j][k] != '-'))
            shared++;
        }
      }
      /* Overlap score between two sequences is the ratio of shared valid
       * residues divided by the sequence length taken as reference. The
       * overlaps matrix, therefore, will be not symmetric. */
      overlaps[i][j] = (float) shared/referenceLength;
    }
  }
}

void newAlignment::getSequences(string *Names) {
    for(int i = 0; i < sequenNumber; i++)
        Names[i] = seqsName[i];
}

void newAlignment::getSequences(string *Names, int *lengths) {
    for(int i = 0; i < sequenNumber; i++) {
        lengths[i] = (int) utils::removeCharacter('-', sequences[i]).length();
        Names[i] = seqsName[i];
    }
}

void newAlignment::getSequences(string *Names, string *Sequences, int *Lengths) {
    for(int i = 0; i < sequenNumber; i++) {
        Names[i] = seqsName[i];
        Sequences[i] = utils::removeCharacter('-', sequences[i]);
        Lengths[i] = (int) Sequences[i].length();
    }
}

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

int newAlignment::getAlignmentType(void) {
    if(dataType == 0)
        dataType = utils::checkAlignmentType(sequenNumber, residNumber, sequences);
    return dataType;
}

int *newAlignment::getCorrespResidues(void) {
    return saveResidues;
}

int *newAlignment::getCorrespSequences(void) {
    return saveSequences;
}

bool newAlignment::isFileAligned(void) {
    return isAligned;
}

void newAlignment::fillNewDataStructure(string *newMatrix, string *newNames) {
    int i, j, k;

    /* Copy only those sequences/columns selected */
    for(i = 0, j = 0; i < sequenNumber; i++) {
        if(saveSequences[i] == -1)
            continue;

        newNames[j] = seqsName[i];
//         cout << "Sequen " << i << endl;
        for(k = 0; k < residNumber; k++) {
//             if (i == 2) cout << "Residue " << k << endl;
            if(saveResidues[k] == -1)
                continue;
            newMatrix[j].resize(newMatrix[j].size() + 1, sequences[i][k]);
        }
        j++;
    }
}

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

bool newAlignment::prepareCodingSequence(bool splitByStopCodon, bool ignStopCodon,\
        newAlignment *proteinAlig) {

    if (getAlignmentType() == SequenceTypes::AA) {
        cerr << endl << "ERROR: Check input CDS file. It seems to content protein "
             << "residues." << endl << endl;
        return false;
    }
    
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

bool newAlignment::checkCorrespondence(string *names, int *lengths, int \
                                       totalInputSeqs, int multiple = 1) {

    int i, j, seqLength, indet;
    bool warnings = false;
    string tmp;

    /* For each sequence in the current protein newAlignment, look for its coding
     * DNA sequence checking that they have the same size. */
    for(i = 0; i < sequenNumber; i++) {

        /* Get protein sequence length removing any pofileible gap. Get as well last
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
                 * it and continue the verification procefile. It will used the 'Nth'
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
                 * sequence. That ifileue could be cause by some incomplete codons in the
                 * nucleotide sequences. This ifileue is solved adding as much 'N' symbols
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
                    cerr << "WARNING: Sequence \"" << seqsName[i] << "\" has lefile nucleo"
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
//     residuesNumber = new int[sequenNumber];
//     for(i = 0; i < sequenNumber; i++) {
//         residuesNumber[i] = sequences[i].size();
//     }

    /* Check whether there are any unknow/no allowed character in the sequences */
    for(i = 0; i < sequenNumber; i++)
        for(j = 0; j < sequences[i].length(); j++)
            if((!isalpha(sequences[i][j])) && (!ispunct(sequences[i][j]))) {
                cerr << endl << "ERROR: The sequence \"" << seqsName[i] << "\" has an "
                     << "unknown (" << sequences[i][j] << ") character." << endl;
                return false;
            }

    /* Check whether all sequences have same size or not */
    for(i = 1; i < sequenNumber; i++)
    {
//         cout << sequences[i].length() << " " << sequences[i] << endl;
        if(sequences[i].length() != sequences[i-1].length())
            break;
    }
//     cout << endl;
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
        residNumber = sequences[0].length();

    /* Check whether aligned sequences have the length fixed for the input alig */
    for(i = 0; (i < sequenNumber) and (aligned); i++) {
        if(sequences[i].length() != residNumber) {
            cerr << endl << "ERROR: The sequence \"" << seqsName[i] << "\" ("
                 << sequences[i].length() << ") does not have the same number of residues "
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
         * procefile, afileign to each sequence its position is useful to know which
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
    min = sequences[0].length();

    for(i = 0, total_res = 0; i < sequenNumber; i++) {

        /* Discard gaps from current sequence and then compute real length */
        for(j = 0, valid_res = 0; j < sequences[i].length(); j++)
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

void newAlignment::printSeqIdentity(void) {

    int i, j, k, pos, maxLongName;
    float mx, avg, maxAvgSeq = 0, maxSeq = 0, avgSeq = 0, **maxs;

    /* Ask for the sequence identities afileesment */
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
    
    for(i = 0; i < sequenNumber; i++) {
        delete maxs[i]; 
    }
    delete[] maxs;
}

void newAlignment::printSeqOverlap()
{
  int i, j, k, pos, maxLongName;
  float mx, avg, maxAvgSeq = 0, maxSeq = 0, avgSeq = 0, **maxs;

  /* Ask for the sequence identities afileesment */
  if(overlaps == NULL)
    calculateSeqOverlap();

  /* For each sequence, we look for its most similar one */
  maxs = new float*[sequenNumber];

  for(i = 0; i < sequenNumber; i++) {
    maxs[i] = new float[2];

    /* Get the most similar sequence to the current one in term of overlap */
    for(k = 0, mx = 0, avg = 0, pos = i; k < sequenNumber; k++) {
      if(i != k) {
        avg += overlaps[i][k];
        if(mx < overlaps[i][k]) {
          mx = overlaps[i][k];
          pos = k;
        }
      }
    }
    /* Update global average variables*/
    avgSeq += avg/(sequenNumber - 1);
    maxAvgSeq += mx;

    /* Save the maximum average overlap value for each sequence */
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

  cout << "## MaxOverlap\t" << maxSeq;
  cout << endl << "#> MaxOverlap\tGet the maximum overlap value for any pair "
    << "of sequences in the alignment" << endl;

  cout << endl << "## AverageOverlap\t" << avgSeq;
  cout << endl << "#> AverageOverlap\tAverage overlap between all sequences";

  cout << endl << endl << "## Overlap sequences matrix";
  for(i = 0; i < sequenNumber; i++) {
    cout << endl << setw(maxLongName + 2) << left << seqsName[i] << "\t";
    for(j = 0; j < sequenNumber; j++)
      cout << setiosflags(ios::left) << setw(10) << overlaps[i][j] << "\t";
  }
  cout << endl;
  
  for(i = 0; i < sequenNumber; i++) {
        delete maxs[i]; 
  }
  delete[] maxs;
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
    indet = getAlignmentType() == SequenceTypes::AA ? 'X' : 'N';
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
         * to speed-up the procefile */
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

// bool newAlignment::alignmentSummaryHTML(char *destFile, int residues, int seqs,
//   int *selectedRes, int *selectedSeq, float *consValues, float blocks) {
// 
//     /* Generate an HTML file with a visual summary about which sequences/columns
//      * have been selected and which have not */
// 
//     int i, j, k, kj, upper, sequencesNamesLength, maxLongName, *gapsValues;
//     string tmpColumn;
//     float *simValues;
//     bool *res, *seq;
//     ofstream file;
//     char type;
// 
//     /* Allocate some local memory */
//     tmpColumn.reserve(sequenNumber);
// 
//     /* Check whether sequences in the alignment are aligned or not.
//      * Warn about it if there are not aligned. */
//     if (!isAligned) {
//         cerr << endl << "ERROR: Sequences are not aligned." << endl << endl;
//         return false;
//     }
// 
//     /* Open output file and check that file pointer is valid */
//     file.open(destFile);
//     if(!file)
//         return false;
// 
//     /* Compute maximum sequences name length. */
//     maxLongName = 0;
//     for(i = 0; i < sequenNumber; i++)
//         maxLongName = utils::max(maxLongName, seqsName[i].size());
// 
//     /* Compute HTML blank spaces */
//     sequencesNamesLength = utils::max(25, maxLongName + 10);
// 
//     /* Initialize local variables to control which columns/sequences
//      * will be kept in the output alignment */
//     
//     res = new bool[residNumber];
//     for(i = 0; i < residNumber; i++)
//         res[i] = false;
// 
//     seq = new bool[sequenNumber];
//     for(i = 0; i < sequenNumber; i++)
//         seq[i] = false;
// 
//     /* Record which columns/sequences from original alignment
//      * have been kept in the final one */
//     for(i = 0; i < residues; i++)
//         res[selectedRes[i]] = true;
//     for(i = 0; i < seqs; i++)
//         seq[selectedSeq[i]] = true;
// 
//     /* Recover some stats about different scores from current alignment */
//     gapsValues = NULL;
//     if (sgaps != NULL)
//         gapsValues = sgaps -> getGapsWindow();
//     simValues = NULL;
//     if (scons != NULL)
//         simValues = scons -> getMdkwVector();
// 
//     /* Print HTML header into output file */
//     file << "<!DOCTYPE html>" << endl << "<html><head>" << endl << "    <meta "
//          << "http-equiv=\"Content-Type\" content=\"text/html;charset=ISO-8859-1\" />"
//          << endl << "    <title>trimAl v1.4 Summary</title>" << endl
//          << "    <style type=\"text/cfile\" media=\"all\">" << endl
// 
//          << "    #b  { background-color: #3366ff; }\n"
//          << "    #r  { background-color: #cc0000; }\n"
//          << "    #g  { background-color: #33cc00; }\n"
//          << "    #p  { background-color: #ff6666; }\n"
//          << "    #m  { background-color: #cc33cc; }\n"
//          << "    #o  { background-color: #ff9900; }\n"
//          << "    #c  { background-color: #46C7C7; }\n"
//          << "    #y  { background-color: #FFFF00; }\n"
// 
//          << "    .sel  { background-color: #B9B9B9; }\n"
//          << "    .nsel { background-color: #E9E9E9; }\n"
// 
//          /* Sets of colors for high-lighting scores intervals */
//          << "    .c1   { background-color: #FFFBF2; }\n"
//          << "    .c2   { background-color: #FFF8CC; }\n"
//          << "    .c3   { background-color: #FAF0BE; }\n"
//          << "    .c4   { background-color: #F0EAD6; }\n"
//          << "    .c5   { background-color: #F3E5AB; }\n"
//          << "    .c6   { background-color: #F4C430; }\n"
//          << "    .c7   { background-color: #C2B280; color: white; }\n"
//          << "    .c8   { background-color: #DAA520; color: white; }\n"
//          << "    .c9   { background-color: #B8860B; color: white; }\n"
//          << "    .c10  { background-color: #918151; color: white; }\n"
//          << "    .c11  { background-color: #967117; color: white; }\n"
//          << "    .c12  { background-color: #6E5411; color: white; }\n"
// 
//          /* Other HTML elements */
//          << "    </style>\n  </head>\n\n" << "  <body>\n" << "  <pre>" << endl;
// 
//     /* Show information about how many sequences/residues have been selected */
//     file << "    <span class=sel>Selected Sequences: " << setw(5) << right << seqs
//          <<" /Selected Residues: " << setw(7) << right << residues << "</span>"
//          << endl << "    <span class=nsel>Deleted Sequences:  " << setw(5) << right
//          << sequenNumber - seqs << " /Deleted Residues:  " << setw(7) << right
//          << residNumber - residues << "</span>" << endl;
// 
//     /* Print headers for different scores derived from input alignment/s */
//     if (gapsValues != NULL)
//         file << endl << setw(sequencesNamesLength) << left << "    Gaps Scores:        "
//              << "<span  class=c1>  =0=  </span><span  class=c2> <.001 </span>"
//              << "<span  class=c3> <.050 </span><span  class=c4> <.100 </span>"
//              << "<span  class=c5> <.150 </span><span  class=c6> <.200 </span>"
//              << "<span  class=c7> <.250 </span><span  class=c8> <.350 </span>"
//              << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
//              << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";
// 
//     if (simValues != NULL)
//         file << endl << setw(sequencesNamesLength) << left << "    Similarity Scores:  "
//              << "<span  class=c1>  =0=  </span><span  class=c2> <1e-6 </span>"
//              << "<span  class=c3> <1e-5 </span><span  class=c4> <1e-4 </span>"
//              << "<span  class=c5> <.001 </span><span  class=c6> <.010 </span>"
//              << "<span  class=c7> <.100 </span><span  class=c8> <.250 </span>"
//              << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
//              << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";
// 
//     if (consValues != NULL)
//         file << endl << setw(sequencesNamesLength) << left << "    Consistency Scores: "
//              << "<span  class=c1>  =0=  </span><span  class=c2> <.001 </span>"
//              << "<span  class=c3> <.050 </span><span  class=c4> <.100 </span>"
//              << "<span  class=c5> <.150 </span><span  class=c6> <.200 </span>"
//              << "<span  class=c7> <.250 </span><span  class=c8> <.350 </span>"
//              << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
//              << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";
// 
//     if ((gapsValues != NULL) or (simValues == NULL) or (consValues == NULL))
//         file << endl;
// 
//     /* Print Sequences in block of BLOCK_SIZE */
//     for(j = 0, upper = HTMLBLOCKS; j < residNumber; j += HTMLBLOCKS, upper += \
//     HTMLBLOCKS) {
// 
//         /* Print main columns number */
//         file << endl << setw(sequencesNamesLength + 10) << right << (j + 10);
//         for(i = j + 20; ((i <= residNumber) && (i <= upper)); i += 10)
//             file << setw(10) << right << (i);
// 
//         /* Print special characters to delimit sequences blocks */
//         file << endl << setw(sequencesNamesLength + 1) << right;
//         for(i = j + 1; ((i <= residNumber) && (i <= upper)); i++)
//             file << (!(i % 10) ? "+" : "=");
//         file << endl;
// 
//         /* Print sequences name */
//         for(i = 0; i < sequenNumber; i++) {
//             file << "    <span class=" << ((seq[i]) ? "sel>" : "nsel>") << seqsName[i]
//                  << "</span>" << setw(sequencesNamesLength - 4 - seqsName[i].size()) << right << "";
// 
//             /* Print residues corresponding to current sequences block */
//             for(k = j; ((k < residNumber) && (k < upper)); k++) {
//                 for(kj = 0, tmpColumn.clear(); kj < sequenNumber; kj++)
//                     tmpColumn += sequences[kj][k];
//                 /* Determine residue color based on residues acrofile the alig column */
//                 type = utils::determineColor(sequences[i][k], tmpColumn);
//                 if (type == 'w')
//                     file << sequences[i][k];
//                 else
//                     file << "<span id=" << type << ">" << sequences[i][k] << "</span>";
//             }
//             file << endl;
//         }
// 
//         file << endl << setw(sequencesNamesLength) << left << "    Selected Cols:      ";
//         for(k = j; ((k < residNumber) && (k < (j + HTMLBLOCKS))); k++)
//             file << "<span class=" << (res[k] ? "sel" : "nsel") << "> </span>";
//         file << endl;
// 
//         /* If there is not any score to print, skip this part of the function */
//         if ((gapsValues == NULL) and (simValues == NULL) and (consValues == NULL))
//             continue;
// 
//         /* Print score colors according to certain predefined thresholds */
//         if (gapsValues != NULL) {
//             file << endl << setw(sequencesNamesLength) << left << "    Gaps Scores:        ";
//             for(k = j; ((k < residNumber) && (k < (j + HTMLBLOCKS))); k++)
//                 if(gapsValues[k] == 0)
//                     file << "<span class=c12> </span>";
//                 else if(gapsValues[k] == sequenNumber)
//                     file << "<span class=c1> </span>";
//                 else if(1 - (float(gapsValues[k])/sequenNumber) >= .750)
//                     file << "<span class=c11> </span>";
//                 else if(1 - (float(gapsValues[k])/sequenNumber) >= .500)
//                     file << "<span class=c10> </span>";
//                 else if(1 - (float(gapsValues[k])/sequenNumber) >= .350)
//                     file << "<span  class=c9> </span>";
//                 else if(1 - (float(gapsValues[k])/sequenNumber) >= .250)
//                     file << "<span  class=c8> </span>";
//                 else if(1 - (float(gapsValues[k])/sequenNumber) >= .200)
//                     file << "<span  class=c7> </span>";
//                 else if(1 - (float(gapsValues[k])/sequenNumber) >= .150)
//                     file << "<span  class=c6> </span>";
//                 else if(1 - (float(gapsValues[k])/sequenNumber) >= .100)
//                     file << "<span  class=c5> </span>";
//                 else if(1 - (float(gapsValues[k])/sequenNumber) >= .050)
//                     file << "<span  class=c4> </span>";
//                 else if(1 - (float(gapsValues[k])/sequenNumber) >= .001)
//                     file << "<span  class=c3> </span>";
//                 else
//                     file << "<span  class=c2> </span>";
//         }
//         if (simValues != NULL) {
//             file << endl << setw(sequencesNamesLength) << left << "    Similarity Scores:  ";
//             for(k = j; ((k < residNumber) && (k < (j + HTMLBLOCKS))); k++)
//                 if(simValues[k] == 1)
//                     file << "<span class=c12> </span>";
//                 else if(simValues[k] == 0)
//                     file << "<span class=c1> </span>";
//                 else if(simValues[k] >= .750)
//                     file << "<span class=c11> </span>";
//                 else if(simValues[k] >= .500)
//                     file << "<span class=c10> </span>";
//                 else if(simValues[k] >= .250)
//                     file << "<span  class=c9> </span>";
//                 else if(simValues[k] >= .100)
//                     file << "<span  class=c8> </span>";
//                 else if(simValues[k] >= .010)
//                     file << "<span  class=c7> </span>";
//                 else if(simValues[k] >= .001)
//                     file << "<span  class=c6> </span>";
//                 else if(simValues[k] >= 1e-4)
//                     file << "<span  class=c5> </span>";
//                 else if(simValues[k] >= 1e-5)
//                     file << "<span  class=c4> </span>";
//                 else if(simValues[k] >= 1e-6)
//                     file << "<span  class=c3> </span>";
//                 else
//                     file << "<span  class=c2> </span>";
//         }
//         if (consValues != NULL) {
//             file << endl << setw(sequencesNamesLength) << left << "    Consistency Scores: ";
//             for(k = j; ((k < residNumber) && (k < (j + HTMLBLOCKS))); k++)
//                 if(consValues[k] == 1)
//                     file << "<span class=c12> </span>";
//                 else if(consValues[k] == 0)
//                     file << "<span class=c1> </span>";
//                 else if(consValues[k] >= .750)
//                     file << "<span class=c11> </span>";
//                 else if(consValues[k] >= .500)
//                     file << "<span class=c10> </span>";
//                 else if(consValues[k] >= .350)
//                     file << "<span  class=c9> </span>";
//                 else if(consValues[k] >= .250)
//                     file << "<span  class=c8> </span>";
//                 else if(consValues[k] >= .200)
//                     file << "<span  class=c7> </span>";
//                 else if(consValues[k] >= .150)
//                     file << "<span  class=c6> </span>";
//                 else if(consValues[k] >= .100)
//                     file << "<span  class=c5> </span>";
//                 else if(consValues[k] >= .050)
//                     file << "<span  class=c4> </span>";
//                 else if(consValues[k] >= .001)
//                     file << "<span  class=c3> </span>";
//                 else
//                     file << "<span  class=c2> </span>";
//         }
//         file << endl;
//     }
// 
//     /* Print HTML footer into output file */
//     file << "    </pre>" << endl << "  </body>" << endl << "</html>" << endl;
// 
//     /* Close output file and deallocate local memory */
//     file.close();
// //     delete [] seq;
// //     delete [] res;
// 
//     return true;
// }

bool newAlignment::alignmentSummaryHTML(char *destFile, int residues, int seqs, int *selectedRes, int *selectedSeq, float *consValues, float blocks) {
    
    int i, j, k, kj, kn, upper;
    double H = 0.0;
    char type;
    
    int * gapsValues = NULL;
    if (sgaps != NULL)
        gapsValues = sgaps -> getGapsWindow();
    cout << "GAPS " << ((sgaps == NULL) ? "NULL" : "NOT NULL") << endl;
    float * simValues = NULL;
    if (scons != NULL)
        simValues = scons -> getMdkwVector();
    cout << "SIMS " << ((scons == NULL) ? "NULL" : "NOT NULL") << endl;
    
    cout << "CONS " << ((consValues == NULL) ? "NULL" : "NOT NULL") << endl;
    
    // Check if alignment is aligned;
    if (!isAligned) {
        cerr << endl << "ERROR: Sequences are not aligned." << endl << endl;
        return false;
    }

    // Get what residues and sequences have been kept
    bool * SEQS = new bool[sequenNumber + 1];
    std::fill(SEQS, SEQS + sequenNumber, false);
    SEQS[sequenNumber] = false;
    bool * RES = new bool[residNumber + 1];
    std::fill(RES, RES + residNumber, false);
    RES[residNumber] = false;
    
    for (int i = 0; i < residues; i++)
        if (selectedRes[i] != -1)
        {
            RES[selectedRes[i]] = true;
        }

    for (int i = 0; i < seqs ; i++)
        if (selectedSeq[i] != -1)
        {
            SEQS[selectedSeq[i]] = true;
        }
    
    // Calculate the blockSize;
    int blockSize = ((int)std::ceil((float)residNumber/ blocks * 0.1F)) * 10;
//     blockSize = 120;
    
    int fontSize = 15;
    
    // Allocate some local memory 
    string tmpColumn; 
    tmpColumn.reserve(sequenNumber);
    
    // Open the file;
    ofstream file;
    file.open(destFile);
    if(!file)
        return false;
    
    /* Compute HTML blank spaces */
    j = 0;
    for(i = 0; i < sequenNumber; i++)
        j = utils::max(j, seqsName[i].size());

    int sequencesNamesLength = utils::max(25, j + 20);
    
    bool textured = true;
    
    // Init Colors
    auto withTexture = []() { 
        
        return std::map<char, string>
        {
            {'o', "url(#colors-orange)"},
            {'y', "url(#colors-yellow)"},
            {'b', "url(#colors-blue)"},
            {'w', "lightgrey"},
            {'p', "darkviolet"},
            {'r', "url(#colors-red)"},
            {'g', "url(#colors-lime)"},
            {'m', "url(#colors-magenta)"},
            {'c', "url(#colors-light-blue)"}
        };
        
    };
    
    auto withoutTexture = []() { 
        
        return std::map<char, string>
        {
            {'o', "orange"},
            {'y', "yellow"},
            {'b', "royalblue"},
            {'w', "lightgrey"},
            {'p', "darkviolet"},
            {'r', "red"},
            {'g', "lime"},
            {'m', "magenta"},
            {'c', "aqua"}
        };
        
    };
    
    std::map<char, string> mappedColors = textured ? withTexture() : withoutTexture() ;
    
    // Start the html output
//     file    << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl
//             << "<title>trimAl v1.4 Summary</title>" << endl
//             << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
//     
    // Start the svg output
    file    << "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" version=\"1.1\" height=\"" <<
            125 +                           /*Header Height*/
            35  +
            (fontSize * 8 +                     /*Column Numbering*/ 
             fontSize * sequenNumber +        /*Sequences Height*/ 
             fontSize * ((gapsValues == NULL) ? 0 : 1) + 
             fontSize * ((simValues == NULL) ? 0 : 1) + 
             fontSize * ((consValues == NULL) ? 0 : 1) 
            )                        /*Selected Residues*/
            * std::ceil(sequences[0].length() / (float)blockSize) /*Blocks Number*/ <<"\"\
            width=\""<< std::max(1430.F, (( std::ceil(sequencesNamesLength / 2.F) + blockSize + 2) * fontSize)) << "px\">" << endl;
            
    // BEGIN defines
    file  << "<defs>" << endl;
    
    // Selected Block 
    file << "<pattern id=\"selected-no-focus\" patternUnits=\"userSpaceOnUse\" width=\"10\" height=\"10\"> " << endl;
    file << "<rect width=\"10\" height=\"10\" fill=\"green\"/>";
    file << "<path d=\"M-1,1 l2,-2 M0,10 l10,-10 M9,11 l2,-2\" stroke=\"white\" stroke-width=\"1\"/>";
    file << "</pattern> " << endl;
    
    file << "<pattern id=\"selected-focus\" patternUnits=\"userSpaceOnUse\" width=\"10\" height=\"10\"> " << endl;
    file << "<rect width=\"10\" height=\"10\" fill=\"mediumseagreen\"/>";
    file << "</pattern> " << endl;
    
    // Deleted Block
    file << "<pattern id=\"deleted-no-focus\" patternUnits=\"userSpaceOnUse\" width=\"10\" height=\"10\"> " << endl;
    file << "<rect width=\"10\" height=\"10\" fill=\"red\"/>";
    file << "<path d=\"M-1,1 l2,-2 M0,10 l10,-10 M9,11 l2,-2\" stroke=\"black\" stroke-width=\"2\"/>";
    file << "</pattern> " << endl;
    
    file << "<pattern id=\"deleted-focus\" patternUnits=\"userSpaceOnUse\" width=\"10\" height=\"10\"> " << endl;
    file << "<rect width=\"10\" height=\"10\" fill=\"black\"/>";
    file << "<path d=\"M-1,1 l2,-2 M0,10 l10,-10 M9,11 l2,-2\" stroke=\"red\" stroke-width=\"2\"/>";
    file << "</pattern> " << endl;
    
    // BEGIN COLORS
    
    // COLORS
    file << "<pattern id=\"colors-orange\" patternUnits=\"userSpaceOnUse\" width=\"10\" height=\"10\"> " << endl;
    file << "<rect width=\"10\" height=\"10\" fill=\"#F7BE81\"/>" ;
    file << "<path d=\"M 0 0 L 10 10 M 0 10 L 10 0\" stroke=\"#DF7401\" stroke-width=\"0.5\" fill=\"transparent\"/> ";
    file << "</pattern> " << endl;
    
    file << "<pattern id=\"colors-red\" patternUnits=\"userSpaceOnUse\" width=\"20\" height=\"12\"> " << endl;
    file << "<rect width=\"20\" height=\"12\" fill=\"red\"/>" ;
    file << "<path d=\"M6 12c0-.622-.095-1.221-.27-1.785A5.982 5.982 0 0 0 10 12c1.67 0 3.182-.683 4.27-1.785A5.998 5.998 0 0 0 14 12h2a4 4 0 0 1 4-4V6c-1.67 0-3.182.683-4.27 1.785C15.905 7.22 16 6.622 16 6c0-.622-.095-1.221-.27-1.785A5.982 5.982 0 0 0 20 6V4a4 4 0 0 1-4-4h-2c0 .622.095 1.221.27 1.785A5.982 5.982 0 0 0 10 0C8.33 0 6.818.683 5.73 1.785 5.905 1.22 6 .622 6 0H4a4 4 0 0 1-4 4v2c1.67 0 3.182.683 4.27 1.785A5.998 5.998 0 0 1 4 6c0-.622.095-1.221.27-1.785A5.982 5.982 0 0 1 0 6v2a4 4 0 0 1 4 4h2zm-4 0a2 2 0 0 0-2-2v2h2zm16 0a2 2 0 0 1 2-2v2h-2zM0 2a2 2 0 0 0 2-2H0v2zm20 0a2 2 0 0 1-2-2h2v2zm-10 8a4 4 0 1 0 0-8 4 4 0 0 0 0 8zm0-2a2 2 0 1 0 0-4 2 2 0 0 0 0 4z\" stroke=\"firebrick\" stroke-width=\"1\" fill=\"transparent\"/> ";
    file << "</pattern> " << endl;
    
    file << "<pattern id=\"colors-yellow\" patternUnits=\"userSpaceOnUse\" width=\"4\" height=\"4\"> " << endl;
    file << "<rect width=\"4\" height=\"4\" fill=\"yellow\"/>" ;
    file << "<path d=\"M 0 0 Q 0 4 4 0\" stroke=\"orange\" stroke-width=\"1\" fill=\"transparent\"/> ";
    file << "</pattern> " << endl;
    
    file << "<pattern id=\"colors-magenta\" patternUnits=\"userSpaceOnUse\" width=\"100\" height=\"100\"> " << endl;
    file << "<rect width=\"100\" height=\"100\" fill=\"magenta\"/>" ;
    file << "<path d=\"M11 18c3.866 0 7-3.134 7-7s-3.134-7-7-7-7 3.134-7 7 3.134 7 7 7zm48 25c3.866 0 7-3.134 7-7s-3.134-7-7-7-7 3.134-7 7 3.134 7 7 7zm-43-7c1.657 0 3-1.343 3-3s-1.343-3-3-3-3 1.343-3 3 1.343 3 3 3zm63 31c1.657 0 3-1.343 3-3s-1.343-3-3-3-3 1.343-3 3 1.343 3 3 3zM34 90c1.657 0 3-1.343 3-3s-1.343-3-3-3-3 1.343-3 3 1.343 3 3 3zm56-76c1.657 0 3-1.343 3-3s-1.343-3-3-3-3 1.343-3 3 1.343 3 3 3zM12 86c2.21 0 4-1.79 4-4s-1.79-4-4-4-4 1.79-4 4 1.79 4 4 4zm28-65c2.21 0 4-1.79 4-4s-1.79-4-4-4-4 1.79-4 4 1.79 4 4 4zm23-11c2.76 0 5-2.24 5-5s-2.24-5-5-5-5 2.24-5 5 2.24 5 5 5zm-6 60c2.21 0 4-1.79 4-4s-1.79-4-4-4-4 1.79-4 4 1.79 4 4 4zm29 22c2.76 0 5-2.24 5-5s-2.24-5-5-5-5 2.24-5 5 2.24 5 5 5zM32 63c2.76 0 5-2.24 5-5s-2.24-5-5-5-5 2.24-5 5 2.24 5 5 5zm57-13c2.76 0 5-2.24 5-5s-2.24-5-5-5-5 2.24-5 5 2.24 5 5 5zm-9-21c1.105 0 2-.895 2-2s-.895-2-2-2-2 .895-2 2 .895 2 2 2zM60 91c1.105 0 2-.895 2-2s-.895-2-2-2-2 .895-2 2 .895 2 2 2zM35 41c1.105 0 2-.895 2-2s-.895-2-2-2-2 .895-2 2 .895 2 2 2zM12 60c1.105 0 2-.895 2-2s-.895-2-2-2-2 .895-2 2 .895 2 2 2z\" stroke=\"mediumorchid\" stroke-width=\"1\" fill=\"mediumorchid\"/> ";
    file << "</pattern> " << endl;
    
    file << "<pattern id=\"colors-lime\" patternUnits=\"userSpaceOnUse\" width=\"28\" height=\"49\"> " << endl;
    file << "<rect width=\"28\" height=\"49\" fill=\"lime\"/>" ;
    file << "<path d=\"M13.99 9.25l13 7.5v15l-13 7.5L1 31.75v-15l12.99-7.5zM3 17.9v12.7l10.99 6.34 11-6.35V17.9l-11-6.34L3 17.9zM0 15l12.98-7.5V0h-2v6.35L0 12.69v2.3zm0 18.5L12.98 41v8h-2v-6.85L0 35.81v-2.3zM15 0v7.5L27.99 15H28v-2.31h-.01L17 6.35V0h-2zm0 49v-8l12.99-7.5H28v2.31h-.01L17 42.15V49h-2z\" stroke=\"green\" stroke-width=\"0.4\" fill=\"transparent\"/> ";
    file << "</pattern> " << endl;
    
    file << "<pattern id=\"colors-light-blue\" patternUnits=\"userSpaceOnUse\" width=\"100\" height=\"20\"> " << endl;
    file << "<rect width=\"100\" height=\"20\" fill=\"aqua\"/>" ;
    file << "<path d=\"M21.184 20c.357-.13.72-.264 1.088-.402l1.768-.661C33.64 15.347 39.647 14 50 14c10.271 0 15.362 1.222 24.629 4.928.955.383 1.869.74 2.75 1.072h6.225c-2.51-.73-5.139-1.691-8.233-2.928C65.888 13.278 60.562 12 50 12c-10.626 0-16.855 1.397-26.66 5.063l-1.767.662c-2.475.923-4.66 1.674-6.724 2.275h6.335zm0-20C13.258 2.892 8.077 4 0 4V2c5.744 0 9.951-.574 14.85-2h6.334zM77.38 0C85.239 2.966 90.502 4 100 4V2c-6.842 0-11.386-.542-16.396-2h-6.225zM0 14c8.44 0 13.718-1.21 22.272-4.402l1.768-.661C33.64 5.347 39.647 4 50 4c10.271 0 15.362 1.222 24.629 4.928C84.112 12.722 89.438 14 100 14v-2c-10.271 0-15.362-1.222-24.629-4.928C65.888 3.278 60.562 2 50 2 39.374 2 33.145 3.397 23.34 7.063l-1.767.662C13.223 10.84 8.163 12 0 12v2z \" stroke=\"blue\" stroke-width=\"0.3\" fill=\"transparent\"/> ";
    file << "</pattern> " << endl;
    
    file << "<pattern id=\"colors-blue\" patternUnits=\"userSpaceOnUse\" width=\"100\" height=\"18\"> " << endl;
    file << "<rect width=\"100\" height=\"18\" fill=\"deepskyblue\"/>" ;
    file << "<path d=\"M61.82 18c3.47-1.45 6.86-3.78 11.3-7.34C78 6.76 80.34 5.1 83.87 3.42 88.56 1.16 93.75 0 100 0v6.16C98.76 6.05 97.43 6 96 6c-9.59 0-14.23 2.23-23.13 9.34-1.28 1.03-2.39 1.9-3.4 2.66h-7.65zm-23.64 0H22.52c-1-.76-2.1-1.63-3.4-2.66C11.57 9.3 7.08 6.78 0 6.16V0c6.25 0 11.44 1.16 16.14 3.42 3.53 1.7 5.87 3.35 10.73 7.24 4.45 3.56 7.84 5.9 11.31 7.34zM61.82 0h7.66a39.57 39.57 0 0 1-7.34 4.58C57.44 6.84 52.25 8 46 8S34.56 6.84 29.86 4.58A39.57 39.57 0 0 1 22.52 0h15.66C41.65 1.44 45.21 2 50 2c4.8 0 8.35-.56 11.82-2z\" stroke=\"blue\" stroke-width=\"0.3\" fill=\"transparent\"/> ";
    file << "</pattern> " << endl;
    // END COLORS
    
    // BEGIN SCORES COLORS
    for (i = 0; i < 12; i++)
    {
        file << "<pattern id=\"score-" << i <<"\" patternUnits=\"userSpaceOnUse\" width=\"10\" height=\"10\"> " << endl;
        file << "<rect width=\"10\" height=\"10\" fill=\""; 
        switch(i) { 
            case 0:
                file << "#FFFBF2"; break;
            case 1:
                file << "#FFF8CC"; break;
            case 2:
                file << "#FAF0BE"; break;
            case 3:
                file << "#F0EAD6"; break;
            case 4:
                file << "#F3E5AB"; break;
            case 5:
                file << "#F4C430"; break;
            case 6:
                file << "#C2B280"; break;
            case 7:
                file << "#DAA520"; break;
            case 8:
                file << "#B8860B"; break;
            case 9:
                file << "#918151"; break;
            case 10:
                file << "#967117"; break;
            case 11:
                file << "#6E5411"; break;
        } 
        file << "\"/>";
        file << "<path d=\"M0,10 l10,-10 \" stroke=\""; 
        switch(i) { 
            case 0:
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
                file << "black"; break;
            case 6:
            case 7:
            case 8:
            case 9:
            case 10:
            case 11:
                file << "white"; break;
        } 
        file << "\" stroke-width=\"1\"/>";
        file << "</pattern> " << endl;
        
    }

    // END SCORES COLORS
    
    file << "</defs>" << endl;
    // END defines
    
    // BEGIN INFO
    // Calculate Info Header Lengths
    strtok(&filename[0], " ");
    string fname = "Filename: ";
    fname.append(strtok(nullptr, ";"));
    
    string sSequences = "Selected sequences: " + std::to_string(seqs) + " / " + std::to_string(sequenNumber);
    
    string sResidues =  "Selected residues:  " + std::to_string(residues) + " / "+  std::to_string(residNumber);
    
    string rSequences = "Deleted sequences:  " + std::to_string(sequenNumber - seqs) + " / " + std::to_string(sequenNumber);
    
    string rResidues =  "Deleted residues:   " + std::to_string(residNumber - residues) + " / "+  std::to_string(residNumber);
    
    
    int size = fname.length();
    size = std::max(size, (int)sSequences.length());
    size = std::max(size, (int)rSequences.length());
    size = std::max(size, (int)sResidues.length());
    size = std::max(size, (int)rResidues.length());
    
    size *= fontSize * 0.5F;
    
    //Filename
    file  << "<g class=\"bar\">" << endl;
    file  << "<rect style=\"fill:indianred\" height=\"20\" width=\"" << size <<"px\" x =\"0\" y=\"" << H << "\" dy=\".35em\" />" << endl;
    file  << "<text style=\"fill:black\" x =\"10\" y=\"" << (H+10) << "\" dy=\".35em\" text-anchor=\"start\" font-family = \"monospace\" xml:space=\"preserve\" kerning=\"0\" \
        lengthAdjust=\"spacingAndGlyphs\"\
        textLength='" << std::min(size - 10, (int)(fname.length() * fontSize * 0.5F)) << "' >"
        << fname << endl;
    file  << "</text>"<< endl;
    file  << "</g>"<< endl;;
    
    H += 25;
    
    // Selected Sequences
    file  << "<g class=\"bar\">"
        << "<rect style=\"fill:lightgrey\" height=\"20\" width=\"" << size << "px\" x =\"0\" y=\"" << (H) << "\" dy=\".35em\" />"
        << "<text style=\"fill:black\" width=\"200px\" x =\"10\" y=\"" << (H+10) << "\" dy=\".35em\" text-anchor=\"start\" font-family = \"monospace\" xml:space=\"preserve\" kerning=\"0\"\
        lengthAdjust=\"spacingAndGlyphs\"\
         textLength='" << std::min(size - 10, (int)(sSequences.length() * fontSize * 0.5F)) << "' >" 
        << sSequences
        << "</text>"
        << "</g>";
        
    H += 25;
    
    // Deleted Sequences
    file  << "<g class=\"bar\">"
        << "<rect style=\"fill:grey\" height=\"20\" width=\"" << size << "px\" x =\"0\" y=\"" << (H) << "\" dy=\".35em\" />"
        << "<text style=\"fill:black\" width=\"200px\" x =\"10\" y=\"" << (H+10) << "\" dy=\".35em\" text-anchor=\"start\" font-family = \"monospace\" xml:space=\"preserve\" kerning=\"0\"\
        lengthAdjust=\"spacingAndGlyphs\"\
         textLength='" << std::min(size - 10, (int)(rSequences.length() * fontSize * 0.5F))  << "' >" 
        << rSequences
        << "</text>"
        << "</g>";

    H += 25;
    
    // Selected Residues
    file  << "<g class=\"bar\">"
        << "<rect style=\"fill:lightgrey\" height=\"20\" width=\"" << size << "px\" x =\"0\" y=\"" << (H) << "\" dy=\".35em\" />"
        << "<text style=\"fill:black\" width=\"200px\" x =\"10\" y=\"" << (H+10) << "\" dy=\".35em\" text-anchor=\"start\" font-family = \"monospace\" xml:space=\"preserve\" kerning=\"0\"\
        lengthAdjust=\"spacingAndGlyphs\"\
         textLength='" << std::min(size - 10, (int)(sResidues.length() * fontSize * 0.5F)) << "' >" 
        << sResidues
        << "</text>"
        << "</g>";

    H += 25;

    // Deleted Residues
    file  << "<g class=\"bar\">"
        << "<rect style=\"fill:grey\" height=\"20\" width=\"" << size << "px\" x =\"0\" y=\"" << (H) << "\" dy=\".35em\" />"
        << "<text style=\"fill:black\" width=\"200px\" x =\"10\" y=\"" << (H+10) << "\" dy=\".35em\" text-anchor=\"start\" font-family = \"monospace\" xml:space=\"preserve\" kerning=\"0\" \
        lengthAdjust=\"spacingAndGlyphs\"\
         textLength='" << std::min(size - 10, (int)(rResidues.length() * fontSize * 0.5F)) << "' >" 
        << rResidues
        << "</text>"
        << "</g>";
        
    // END INFO
        
    H = 0;
        
    //BEGIN Legend
        file << "<rect \
                    style=\"fill:"<< mappedColors['b'] <<"\" \
                    height=\""<< 20 <<"\" \
                    width=\"" << 120 << "px\" \
                    x =\"" << size + 10 << "px\"\
                    y =\""<< (H) <<"\"/>" << endl;
        file << "<text \
                    width=\"100\" style=\"font-weight:bold\" \
                    text-anchor=\"middle\" \
                    x =\"" << size + 70 << "px\" \
                    y =\""<< (H + 15) <<"\" \
                    font-family = \"monospace\" \
                    kerning=\"0\" >" << endl
             << "Hidrophobic" << endl
             << "</text>" << endl;
             
        file << "<rect \
                    style=\"fill:"<< mappedColors['r'] <<"\" \
                    height=\""<< 20 <<"\" \
                    width=\"" << 120 << "px\" \
                    x =\"" << size + 10 + 120<< "px\"\
                    y =\""<< (H) <<"\"/>" << endl;
        file << "<text \
                    width=\"100\" style=\"font-weight:bold\"\
                    text-anchor=\"middle\" \
                    x =\"" << size + 70 + 120 << "px\" \
                    y =\""<< (H + 15) <<"\" \
                    font-family = \"monospace\" \
                    kerning=\"0\" >" << endl
             << "Positive Charge" << endl
             << "</text>" << endl;
             
        file << "<rect \
                    style=\"fill:"<< mappedColors['m'] <<"\" \
                    height=\""<< 20 <<"\" \
                    width=\"" << 120 << "px\" \
                    x =\"" << size + 10 + 120 * 2 << "px\"\
                    y =\""<< (H) <<"\"/>" << endl;
        file << "<text \
                    width=\"100\" style=\"font-weight:bold\"\
                    text-anchor=\"middle\" \
                    x =\"" << size + 70 + 120 * 2 << "px\" \
                    y =\""<< (H + 15) <<"\" \
                    font-family = \"monospace\" \
                    kerning=\"0\" >" << endl
             << "Negative Charge" << endl
             << "</text>" << endl;

        file << "<rect \
                    style=\"fill:"<< mappedColors['g'] <<"\" \
                    height=\""<< 20 <<"\" \
                    width=\"" << 120 << "px\" \
                    x =\"" << size + 10 + 120 * 3 << "px\"\
                    y =\""<< (H) <<"\"/>" << endl;
        file << "<text \
                    width=\"100\" style=\"font-weight:bold\"\
                    text-anchor=\"middle\" \
                    x =\"" << size + 70 + 120 * 3 << "px\" \
                    y =\""<< (H + 15) <<"\" \
                    font-family = \"monospace\" \
                    kerning=\"0\" >" << endl
             << "Polar" << endl
             << "</text>" << endl;
             
        file << "<rect \
                    style=\"fill:"<< mappedColors['p'] <<"\" \
                    height=\""<< 20 <<"\" \
                    width=\"" << 120 << "px\" \
                    x =\"" << size + 10 + 120 * 4 << "px\"\
                    y =\""<< (H) <<"\"/>" << endl;
        file << "<text \
                    width=\"100\" style=\"font-weight:bold\"\
                    text-anchor=\"middle\" \
                    x =\"" << size + 70 + 120 * 4 << "px\" \
                    y =\""<< (H + 15) <<"\" \
                    font-family = \"monospace\" \
                    kerning=\"0\" >" << endl
             << "Cysteine" << endl
             << "</text>" << endl;
             
        file << "<rect \
                    style=\"fill:"<< mappedColors['o'] <<"\" \
                    height=\""<< 20 <<"\" \
                    width=\"" << 120 << "px\" \
                    x =\"" << size + 10 + 120 * 5 << "px\"\
                    y =\""<< (H) <<"\"/>" << endl;
                    
            
        file << "<text \
                    width=\"100\" style=\"font-weight:bold\"\
                    text-anchor=\"middle\" \
                    x =\"" << size + 70 + 120 * 5 << "px\" \
                    y =\""<< (H + 15) <<"\" \
                    font-family = \"monospace\" \
                    kerning=\"0\" >" << endl
             << "Glycine" << endl
             << "</text>" << endl;

        file << "<rect \
                    style=\"fill:"<< mappedColors['y'] <<"\" \
                    height=\""<< 20 <<"\" \
                    width=\"" << 120 << "px\" \
                    x =\"" << size + 10 + 120 * 6 << "px\"\
                    y =\""<< (H) <<"\"/>" << endl;
        file << "<text \
                    width=\"100\" style=\"font-weight:bold\"\
                    text-anchor=\"middle\" \
                    x =\"" << size + 70 + 120 * 6 << "px\" \
                    y =\""<< (H + 15) <<"\" \
                    font-family = \"monospace\" \
                    kerning=\"0\" >" << endl
             << "Proline" << endl
             << "</text>" << endl;
             
        file << "<rect \
                    style=\"fill:"<< mappedColors['c'] <<"\" \
                    height=\""<< 20 <<"\" \
                    width=\"" << 120 << "px\" \
                    x =\"" << size + 10 + 120 * 7 << "px\"\
                    y =\""<< (H) <<"\"/>" << endl;
        file << "<text \
                    width=\"100\" style=\"font-weight:bold\"\
                    text-anchor=\"middle\" \
                    x =\"" << size + 70 + 120 * 7 << "px\" \
                    y =\""<< (H + 15) <<"\" \
                    font-family = \"monospace\" \
                    kerning=\"0\" >" << endl
             << "Aromatic" << endl
             << "</text>" << endl;
             
        file << "<rect \
                    style=\"fill:"<< mappedColors['w'] <<"\" \
                    height=\""<< 20 <<"\" \
                    width=\"" << 120 << "px\" \
                    x =\"" << size + 10 + 120 * 8 << "px\"\
                    y =\""<< (H) <<"\"/>" << endl;
        file << "<text \
                    width=\"100\" style=\"font-weight:bold\"\
                    text-anchor=\"middle\" \
                    x =\"" << size + 70 + 120 * 8 << "px\" \
                    y =\""<< (H + 15) <<"\" \
                    font-family = \"monospace\" \
                    kerning=\"0\">" << endl
             << "Unconserved" << endl
             << "</text>" << endl;
    H+= 25;

    if (gapsValues || consValues || simValues)
    {
        int width = 78;
        
        for (i = 0; i < 12; i++)
        {
            file << "<rect \
                        style=\"fill:"<< "url(#score-"<<i<<");stroke:black;stroke-width:1" <<"\" \
                        height=\""<< 20 <<"\" \
                        width=\"" << width << "px\" \
                        x =\"" << size + 10  + width * (2 + i) << "px\"\
                        y =\""<< (H) <<"\" />" << endl;
        }
        if (gapsValues) {
            std::array<std::string, 12> u = {{"0","&lt;.001","&lt;.050", "&lt;.100","&lt;.150", "&lt;.200", "&lt;.250", "&lt;.350", "&lt;.500", "&lt;.750", "&lt; 1.00", "1"}};
            file << "<text \
                        width=\"100\" style=\"font-weight:bold\" \
                        text-anchor=\"middle\" \
                        x =\"" << size + 10 + width << "px\" \
                        y =\""<< (H + 40) <<"\" \
                        font-family = \"monospace\" \
                        kerning=\"0\" \
                        textLength=\""<< width * 1.5F << "\">" << endl
                << "Gaps Scores" << endl
                << "</text>" << endl;
            for (i = 0; i < 12; i++)
                file << "<text \
                            width=\"100\" style=\"font-weight:bold\" \
                            text-anchor=\"middle\" \
                            x =\"" << size + 10 + width * (2.5F + i) << "px\" \
                            y =\""<< (H + 40) <<"\" \
                            font-family = \"monospace\" \
                            kerning=\"0\" \
                            textLength=\""<< width / 2 << "\">" << endl
                    << u[i] << endl
                    << "</text>" << endl;
            
            H+= 25;
        }
        if (simValues) {
            std::array<std::string, 12> u = {{"0", "&lt; 1e-6", "&lt; 1e-5", "&lt; 1e-4", "&lt;.001","&lt;.010", "&lt;.100", "&lt;.250", "&lt;.500", "&lt;.750", "&lt; 1.00", "1" }};
            file << "<text \
                        width=\"100\" style=\"font-weight:bold\" \
                        text-anchor=\"middle\" \
                        x =\"" << size + 10 + width << "px\" \
                        y =\""<< (H + 40) <<"\" \
                        font-family = \"monospace\" \
                        kerning=\"0\" \
                        textLength=\""<< width * 1.5F << "\">" << endl
                << "Similarity Scores" << endl
                << "</text>" << endl;
            for (i = 0; i < 12; i++)
                file << "<text \
                            width=\"100\" style=\"font-weight:bold\" \
                            text-anchor=\"middle\" \
                            x =\"" << size + 10 + width * (2.5F + i) << "px\" \
                            y =\""<< (H + 40) <<"\" \
                            font-family = \"monospace\" \
                            kerning=\"0\" \
                            textLength=\""<< width / 2 << "\">" << endl
                    << u[i] << endl
                    << "</text>" << endl;
            H+= 25;
        }
        if (consValues) {
            ;
            std::array<std::string, 12> u = {{"0"," <.001"," <.050"," <.100","  <.150"," <.200"," <.250"," <.350"," <.500"," <.750"," <1.00"," =1="}};
            file << "<text \
                        width=\"100\" style=\"font-weight:bold\" \
                        text-anchor=\"middle\" \
                        x =\"" << size + 10 + width << "px\" \
                        y =\""<< (H + 40) <<"\" \
                        font-family = \"monospace\" \
                        kerning=\"0\" \
                        textLength=\""<< width * 1.5F << "\">" << endl
                << "Similarity Scores" << endl
                << "</text>" << endl;
            for (i = 0; i < 12; i++)
                file << "<text \
                            width=\"100\" style=\"font-weight:bold\" \
                            text-anchor=\"middle\" \
                            x =\"" << size + 10 + width * (2.5F + i) << "px\" \
                            y =\""<< (H + 40) <<"\" \
                            font-family = \"monospace\" \
                            kerning=\"0\" \
                            textLength=\""<< width / 2 << "\">" << endl
                    << u[i] << endl
                    << "</text>" << endl;
        }
        
    }

    //END Legend
        
    H = 160;
    
    for(j = 0, upper = blockSize; 
        j < residNumber; 
        j += blockSize, upper += blockSize) {

        /* Print main columns number */
        H += fontSize;

        if ((j + blockSize) < residNumber)
        {
            file    << "<text font-family = \"monospace\" \
                    font-size=\"" << fontSize << "px\" dy=\".35em\"\
                    x =\"" << (sequencesNamesLength) * fontSize * 0.5F << "\" text-anchor=\"start\" y=\"" << (H) << "\" xml:space=\"preserve\" kerning=\"0\" \
                    textLength='"<< ((std::min((float)blockSize, residNumber - j - 0.25F) ) * fontSize )<<"' lengthAdjust=\"spacing\">";
        }
        else
        {
            file    << "<text font-family = \"monospace\" \
                    font-size=\"" << fontSize << "px\" dy=\".35em\" \
                    x =\"" << (sequencesNamesLength) * fontSize * 0.5F << "\" text-anchor=\"start\" y=\"" << (H) << "\" xml:space=\"preserve\" kerning=\"0\" \
                    textLength='"<< ((std::min((float)blockSize, residNumber - j - 0.25F) + std::to_string(j).length()) * fontSize )<<"' lengthAdjust=\"spacing\">";
        }
            
        for(i = j; ((i < residNumber) && (i < upper)); i += 10)            
            file << setw(10) << setfill(' ') << left << i;
        file << "</text>" << endl;
        
        H += fontSize;
        file    << "<text font-family = \"monospace\" \
                font-size=\"" << fontSize << "px\" dy=\".35em\" \
                x =\"" << (sequencesNamesLength) * fontSize * 0.5F << "\" text-anchor=\"start\" y=\"" << (H) << "\" xml:space=\"preserve\" kerning=\"0\" \
                textLength='"<< (std::min((float)blockSize, residNumber - j - 0.25F) * fontSize )<<"' lengthAdjust=\"spacing\">";
            
        for(i = j; ((i < residNumber) && (i < upper)); i += 10)
            
            file << setw(std::min(10, residNumber - i)) << setfill('~') << left << "+";
        file << "</text>" << endl;
        
//         H+=fontSize;

        char colors[sequenNumber][blockSize];
        for(i = 0; i < sequenNumber; i++) {
            
            /* Print residues corresponding to current sequences block */
            for(k = j; ((k < residNumber) && (k < upper)); k++) 
            {
                for(kj = 0, tmpColumn.clear(); kj < sequenNumber; kj++)
                    tmpColumn += sequences[kj][k];
                /* Determine residue color based on residues acrofile the alig column */
                type = utils::determineColor(sequences[i][k], tmpColumn);
                colors[i][k - j] = type;
            }

        }
        

        
        H += fontSize;
        
        float width =  (blockSize + 0.5F) * (fontSize) / blockSize;
        for (k = j; k < residNumber && k < upper; k++)
        {
            for(kj = 0; kj <= sequenNumber; kj++)
            {
                for(int kn = kj + 1; kn <= sequenNumber; kn++)
                {
                    if (colors[kj][k - j] != colors[kn][k - j] || kn + 1 == sequenNumber)
                    {
//                         if (colors[kj][k - j] != 'w')
                         file << "<rect " <<
                                    " style=\"fill:"<< mappedColors[colors[kj][k - j]] <<"\" " <<
                                    " height=\""<< fontSize * (kn - kj) <<"\" " <<
                                    " width=\"" << width << "px\" " <<
                                    " x =\"" << 
                            
                            sequencesNamesLength * fontSize * 0.5F  // Sequences Names
                            + (k - j - 0.25F) * width

                            << "\" " <<
                            " y =\"" << (H+(fontSize / 2)+(fontSize*kj)) << "\" />";
                        kj = kn;
                    }
                }
               
            }
        }
        
        /* Print sequences name */
        for(i = 0; i < sequenNumber; i++) {
            H += fontSize;
            file << "<text font-family = \"monospace\" font-size=\"" << fontSize << "px\" dy=\".35em\" x =\"0\" \
                text-anchor=\"start\" y=\"" << (H) << "\" kerning=\"0\" " << (SEQS[i] ? "style=\"font-weight:bold\"" : "text-decoration=\"line-through\"") << ">"  << 
                setfill(' ') << setw(sequencesNamesLength) << right << seqsName[i] << "</text>" << endl;
            
            file << "<text font-family = \"monospace\" font-size=\"" << fontSize << "px\" dy=\".35em\" x =\"" << sequencesNamesLength * fontSize * 0.5F << "px\" \
                text-anchor=\"start\" y=\"" << H << "\" xml:space=\"preserve\" kerning=\"0\" textLength='"<< (std::min(blockSize, residNumber - j) * fontSize) <<"' lengthAdjust=\"spacing\" " << (SEQS[i] ? "style=\"font-weight:bold\"" : "style=\"font-weight:100\"") << ">" ;
            
                /* Print residues corresponding to current sequences block */
            for(k = j; ((k < residNumber) && (k < upper)); k++) 
            {
                for(kj = 0, tmpColumn.clear(); kj < sequenNumber; kj++)
                    tmpColumn += sequences[kj][k];
                /* Determine residue color based on residues acrofile the alig column */
                type = utils::determineColor(sequences[i][k], tmpColumn);
                file << sequences[i][k];
            }
            file << "</text>" << endl;
            
        }
        
        H += fontSize * 2;

        // SELECTED OR REJECTED SEQUENCES AND RESIDUES
        bool accepted = RES[j];
        int oriPosi = j;
        for (k = j; k - j <= blockSize && k <= residNumber; k++)
        {
            if ( (k - j) == blockSize || ((k) == residNumber) || RES[k] != accepted)
            {
                    file << "<text font-family = \"monospace\" font-size=\"" << fontSize << "px\" dy=\".35em\" x =\"0\" \
                        text-anchor=\"start\" y=\"" << (H + fontSize / 2) << "\" kerning=\"0\"" << ">"  << 
                        setfill(' ') << setw(sequencesNamesLength) << right << "Selected Sequences" << "</text>" << endl;
                
                    file << 
                        "<rect style=\"fill:url(" << (accepted ? "#selected-no-focus" : "#deleted-no-focus") << ");stroke-width:1;stroke:black\" height=\"10\"" <<
                        " width=\"" << width * (k - oriPosi) << "px\"" <<
                        " x =\"" << ((sequencesNamesLength) * fontSize/2) + ((oriPosi - j - 0.25F) * width) << "px\" y=\"" << (H) << "\" dy=\".35em\" " <<
                        " onmouseover=\"evt.target.setAttribute('style', 'fill:url(" << (accepted ? "#selected-focus" : "#deleted-focus") << ");stroke-width:1;stroke:black');\" " <<
                        " onmouseout=\"evt.target.setAttribute('style', 'fill:url(" <<  (accepted ? "#selected-no-focus" : "#deleted-no-focus") << ");stroke-width:1;stroke:black');\"/>" << endl;
            
            }
            
            if (RES[k] != accepted)
            {
                accepted = RES[k];
                oriPosi = k;
            }

        }
        H += fontSize;
        
        // GAPS VALUES
        if (gapsValues)
        {
            float inverse = 1.F / sequenNumber;
            int step = utils::GetGapStep(&gapsValues[j], inverse), innerStep;
            oriPosi = j;
            
            for (k = j; k - j <= blockSize && k <= residNumber; k++)
            {
                if (((k) != residNumber))
                    innerStep = utils::GetGapStep(&gapsValues[k], inverse);
                if ( (k - j) == blockSize || ((k) == residNumber) || innerStep != step)
                {
                    file << "<text font-family = \"monospace\" font-size=\"" << fontSize << "px\" dy=\".35em\" x =\"0\" \
                            text-anchor=\"start\" y=\"" << (H + fontSize / 2) << "\" kerning=\"0\"" << ">"  << 
                            setfill(' ') << setw(sequencesNamesLength) << right << "Gaps Values" << "</text>" << endl;
                            
                    file << 
                            "<rect style=\"fill:url(#score-" << step << ");stroke-width:1.5;stroke:black\" height=\"10\"" <<
                            " width=\"" << width * (k - oriPosi) << "px\"" <<
                            " x =\"" << ((sequencesNamesLength) * fontSize/2) + ((oriPosi - j - 0.25F) * width) << "px\" y=\"" << (H) << "\" dy=\".35em\" " <<
                            " />" << endl;
                
                }
                
                if (innerStep != step)
                {
                    step = innerStep;
                    oriPosi = k;
                }
                /*
                */
            }
            
            H += fontSize;
        }
        
        // SIMILARITY VALUES
        if (simValues)
        {
            int step = utils::GetSimStep(&simValues[j]), innerStep;
            oriPosi = j;
            
            for (k = j; k - j <= blockSize && k <= residNumber; k++)
            {
                if (((k) != residNumber))
                    innerStep = utils::GetSimStep(&simValues[k]);
                if ( (k - j) == blockSize || ((k) == residNumber) || innerStep != step)
                {
                    file << "<text font-family = \"monospace\" font-size=\"" << fontSize << "px\" dy=\".35em\" x =\"0\" \
                            text-anchor=\"start\" y=\"" << (H + fontSize / 2) << "\" kerning=\"0\"" << ">"  << 
                            setfill(' ') << setw(sequencesNamesLength) << right << "Similarily Values" << "</text>" << endl;
                            
                    file << 
                            "<rect style=\"fill:url(#score-" << step << ");stroke-width:1.5;stroke:black\" height=\"10\"" <<
                            " width=\"" << width * (k - oriPosi) << "px\"" <<
                            " x =\"" << ((sequencesNamesLength) * fontSize/2) + ((oriPosi - j - 0.25F) * width) << "px\" y=\"" << (H) << "\" dy=\".35em\" " <<
                            " />" << endl;
                
                }
                
                if (innerStep != step)
                {
                    step = innerStep;
                    oriPosi = k;
                }
            }
            
            H += fontSize;
        }
        
        // CONSISTENCY VALUES
        if (consValues)
        {
            int step = utils::GetConsStep(&consValues[j]), innerStep;
            oriPosi = j;
            
            for (k = j; k - j <= blockSize && k <= residNumber; k++)
            {
                if (((k) != residNumber))
                    innerStep = utils::GetConsStep(&consValues[k]);
                if ( (k - j) == blockSize || ((k) == residNumber) || innerStep != step)
                {
                    file << "<text font-family = \"monospace\" font-size=\"" << fontSize << "px\" dy=\".35em\" x =\"0\" \
                            text-anchor=\"start\" y=\"" << (H + fontSize / 2) << "\" kerning=\"0\" style=\"font-weight:100\"" << ">"  << 
                            setfill(' ') << setw(sequencesNamesLength) << right << "Similarily Values" << "</text>" << endl;
                            
                    file << 
                            "<rect style=\"fill:url(#score-" << step << ");stroke-width:1.5;stroke:black\" height=\"10\"" <<
                            " width=\"" << width * (k - oriPosi) << "px\"" <<
                            " x =\"" << ((sequencesNamesLength) * fontSize/2) + ((oriPosi - j - 0.25F) * width) << "px\" y=\"" << (H) << "\" dy=\".35em\" " <<
                            " />" << endl;
                }
                
                if (innerStep != step)
                {
                    step = innerStep;
                    oriPosi = k;
                }
            }
            
            H += fontSize;
        }
        
        H += fontSize * 2;
    }
    file << "</svg>";
    delete [] SEQS;
    delete [] RES;
    file.close();
    return true;
}

bool newAlignment::alignmentColourHTML(ostream &file) {

    int i, j, kj, upper, k = 0, maxLongName = 0;
    string tmpColumn;
    char type;

    /* Allocate some local memory */
    tmpColumn.reserve(sequenNumber);

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!isAligned) {
        cerr << endl << "ERROR: Sequences are not aligned." << endl << endl;
        return false;
    }

    /* Compute maximum sequences name length */
    maxLongName = 0;
    for(i = 0; i < sequenNumber; i++)
        maxLongName = utils::max(maxLongName, seqsName[i].size());


    /* Print HTML header into output file */
    file << "<!DOCTYPE html>" << endl << "<html><head>" << endl << "    <meta "
         << "http-equiv=\"Content-Type\" content=\"text/html;charset=ISO-8859-1\" />"
         << endl << "    <title>readAl v1.4</title>" << endl
         << "    <style type=\"text/cfile\">" << endl
         << "    #b  { background-color: #3366ff; }\n"
         << "    #r  { background-color: #cc0000; }\n"
         << "    #g  { background-color: #33cc00; }\n"
         << "    #p  { background-color: #ff6666; }\n"
         << "    #m  { background-color: #cc33cc; }\n"
         << "    #o  { background-color: #ff9900; }\n"
         << "    #c  { background-color: #46C7C7; }\n"
         << "    #y  { background-color: #FFFF00; }\n"
         << "    </style>\n  </head>\n\n" << "  <body>\n  <pre>" << endl;

    /* Print sequences colored according to CLUSTAL scheme based on
     * physical-chemical properties */
    for(j = 0, upper = HTMLBLOCKS; j < residNumber; j += HTMLBLOCKS, upper += \
    HTMLBLOCKS) {

        file << endl;
        /* Print main columns number */
        file << setw(maxLongName + 19) << right << (j + 10);
        for(i = j + 20; ((i <= residNumber) && (i <= upper)); i += 10)
            file << setw(10) << right << i;

        /* Print special characters to delimit sequences blocks */
        file << endl << setw(maxLongName + 10);
        for(i = j + 1; ((i <= residNumber) && (i <= upper)); i++)
            file << (!(i % 10) ? "+" : "=");

        /* Print sequences themselves */
        for(i = 0; i < sequenNumber; i++) {

            /* Print sequences name */
            file << endl << setw(maxLongName + 9) << left << seqsName[i];

            /* Print residues corresponding to current sequences block */
            for(k = j; ((k < residNumber) && (k < upper)); k++) {
                for(kj = 0, tmpColumn.clear(); kj < sequenNumber; kj++)
                    tmpColumn += sequences[kj][k];
                /* Determine residue color based on residues acrofile the alig column */
                type = utils::determineColor(sequences[i][k], tmpColumn);
                if (type == 'w')
                    file << sequences[i][k];
                else
                    file << "<span id=" << type << ">" << sequences[i][k] << "</span>";
            }
        }
        file << endl;
    }

    /* Print HTML footer into output file */
    file << "    </pre>" << endl << "  </body>" << endl << "</html>" << endl;

    return true;
}


