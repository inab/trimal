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
#include <cmath>
#include "../include/reportsystem.h"

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
    dataType = SequenceTypes::NotDefined;

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

    SeqRef = new int(1);

}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Overlapping operator = to use it as a kind of class constructor */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

newAlignment::newAlignment(newAlignment& originalAlignment) {

    if(this != &originalAlignment) {

        int i, j;

        aligInfo = originalAlignment.aligInfo;



        for ( i = 0, j = 0; i < originalAlignment.originalSequenNumber; i++)
        {
            if (originalAlignment.saveSequences[i] != -1) j++;
        }

        sequenNumber = j;

        for ( i = 0, j = 0; i < originalAlignment.residNumber; i++)
            if (originalAlignment.saveResidues[i] != -1) j++;

        residNumber =  j;

        isAligned =  originalAlignment.isAligned;

        dataType = originalAlignment.dataType;

        originalSequenNumber = originalAlignment.originalSequenNumber;
        originalResidNumber = originalAlignment.originalResidNumber;

        sequences = originalAlignment.sequences;
        seqsName = originalAlignment.seqsName;
        seqsInfo = originalAlignment.seqsInfo;

        saveSequences = new int[originalSequenNumber];
        std::copy(originalAlignment.saveSequences, originalAlignment.saveSequences + originalAlignment.originalSequenNumber, saveSequences);

        saveResidues = new int[originalResidNumber];
        std::copy(originalAlignment.saveResidues, originalAlignment.saveResidues + originalAlignment.originalResidNumber, saveResidues);

        identities = NULL;

        //delete sgaps;
        sgaps = NULL;

        //delete scons;
        scons = NULL;

        //delete SequencesMatrix;
        SequencesMatrix = NULL;

        this -> Cleaning = new Cleaner(this, originalAlignment.Cleaning);

        this -> Statistics = new StatisticsManager(this, originalAlignment.Statistics);

        this -> SeqRef = originalAlignment.SeqRef;

        (*SeqRef)++;
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Clafile destructor */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

newAlignment::~newAlignment(void) {
    int i;

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

//     dataType = 0;

    delete Cleaning;

    delete Statistics;

    if (--(*SeqRef) == 0)
    {
        delete SeqRef;

        if(sequences != NULL)
            delete [] sequences;
        sequences = NULL;

        if(seqsName != NULL)
            delete [] seqsName;
        seqsName = NULL;

        if(seqsInfo != NULL)
            delete [] seqsInfo;
        seqsInfo = NULL;

    }

}

newAlignment *newAlignment::getTranslationCDS(/*int newResidues, int newSequences, int *ColumnsToKeep, string *oldSeqsName, sequencesMatrix *seqMatrix,*/
    newAlignment *ProtAlig) {

    int x, y, counter, * mappedSeqs;
    newAlignment *newAlig = new newAlignment();

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Map the selected protein sequences to the input
     * coding sequences */
    mappedSeqs = new int[ProtAlig->sequenNumber];

    for (x = 0; x < ProtAlig->originalSequenNumber; x++)
    {
        for (y = 0; y < originalSequenNumber; y++)
        {
            if (ProtAlig->seqsName[x] == seqsName[y])
            {
                mappedSeqs[x] = y;
                break;
            }
        }
//         if (y == originalSequenNumber)
//         {
//             Debug << "Name " << ProtAlig->seqsName[x] << " not found" << endl;
//         }
//         else
//         {
//             Debug << "Seq " << y << " ~ Seq " << x << "\t\t\t" << ProtAlig->seqsName[x] << "\t\t" << seqsName[y] << endl;
//         }
    }

    newAlig->sequences  = new std::string[ProtAlig->sequenNumber];
    newAlig->seqsInfo   = new std::string[ProtAlig->sequenNumber];
    newAlig->seqsName   = new std::string[ProtAlig->sequenNumber];

    for (x = 0; x < ProtAlig->originalSequenNumber; x++)
    {

        newAlig->sequences[x] = *new std::string[ProtAlig->sequences[x].size()];
        std::string & current = newAlig->sequences[x];

        if (ProtAlig->seqsInfo != NULL)
            newAlig->seqsInfo[x] = *new std::string(ProtAlig->seqsInfo[x]);
        newAlig->seqsName[x] = *new std::string(ProtAlig->seqsName[x]);


        for (y = 0, counter = 0; y < ProtAlig->sequences[x].size(); y++)
        {
            if (ProtAlig->saveResidues[y] == -1)
            {
                if (ProtAlig->sequences[x][y] != '-') counter++;
                continue;
            }
            if (ProtAlig->sequences[x][y] == '-') {
                current.push_back('-');
                current.push_back('-');
                current.push_back('-');
                continue;
            }
            current.push_back(sequences[mappedSeqs[x]][counter * 3]);
            current.push_back(sequences[mappedSeqs[x]][counter * 3 + 1]);
            current.push_back(sequences[mappedSeqs[x]][counter * 3 + 2]);
            counter ++;
        }
        current.shrink_to_fit();
    }

    newAlig->sequenNumber = ProtAlig->sequenNumber;
    newAlig->originalSequenNumber = ProtAlig->sequenNumber;

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
}

bool newAlignment::prepareCodingSequence(bool splitByStopCodon, bool ignStopCodon,\
        newAlignment *proteinAlig) {

    if (getAlignmentType() == SequenceTypes::AA) {
        Debug.Report(ErrorCode::CDScontainsProteinSequences);
//         cerr << endl << "ERROR: Check input CDS file. It seems to content protein "
//              << "residues." << endl << endl;
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
            Debug.Report(ErrorCode::SequenceContainsGap, new std::string[1] { seqsName[i] });
//             cerr << "ERROR: Sequence \"" << seqsName[i] << "\" has, at least, one gap"
//                  << endl << endl;
            return false;
        }

        if((sequences[i].length() % 3) != 0) {
            if (!warning)
                cerr << endl;
            warning = true;
            Debug.Report(ErrorCode::SequenceNotMultipleOfThree, new std::string [2] {seqsName[i], std::to_string(sequences[i].length()) });
//             cerr << "WARNING: Sequence length \"" << seqsName[i] << "\" is not "
//                  << "multiple of 3 (length: " << sequences[i].length() << ")" << endl;
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
                    warning = true;
                    Debug.Report(InfoCode::CuttingSequence, new std::string[5] { seqsName[i], "TGA", std::to_string(aminoAcid), std::to_string(found + 1), std::to_string(sequences[i].length()) });
                    sequences[i].resize((int) found);
                }
                /* Otherwise, warn about it and return an error */
                else {
                    Debug.Report(ErrorCode::SequenceHasStopCodon, new std::string[5] { seqsName[i], "TGA", std::to_string(aminoAcid), std::to_string(found + 1), std::to_string(sequences[i].length())});
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
                    warning = true;
                    Debug.Report(InfoCode::CuttingSequence, new std::string[5] { seqsName[i], "TAA", std::to_string(aminoAcid), std::to_string(found + 1), std::to_string(sequences[i].length()) });
                    sequences[i].resize((int) found);
                }
                /* Otherwise, warn about it and return an error */
                else {
                    Debug.Report(ErrorCode::SequenceHasStopCodon, new std::string[5] { seqsName[i], "TAA", std::to_string(aminoAcid), std::to_string(found + 1), std::to_string(sequences[i].length())});
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
                    warning = true;
                    Debug.Report(InfoCode::CuttingSequence, new std::string[5] { seqsName[i], "TAG", std::to_string(aminoAcid), std::to_string(found + 1), std::to_string(sequences[i].length()) });
                    sequences[i].resize((int) found);
                }
                /* Otherwise, warn about it and return an error */
                else {
                    Debug.Report(ErrorCode::SequenceHasStopCodon, new std::string[5] { seqsName[i], "TAG", std::to_string(aminoAcid), std::to_string(found + 1), std::to_string(sequences[i].length())});
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
                    Debug.Report(WarningCode::SequenceWillBeCutted, new std::string[3] { seqsName[i], std::to_string(seqLength), std::to_string(lengths[j])});
//                     cerr << "WARNING: Sequence \"" << seqsName[i] << "\" will be cutted "
//                          << "at position " << seqLength << " (length: "<< lengths[j] << ")"
//                          << endl;
                    break;
                }

                /* It has been detected some indeterminations at the end of the protein
                 * sequence. That ifileue could be cause by some incomplete codons in the
                 * nucleotide sequences. This issue is solved adding as much 'N' symbols
                 * as it is needed to preserve the backtranslated newAlignment */
                else if((indet > 0) && (indet > (seqLength - lengths[j])/3)) {
                    if (!warnings)
                        cerr << endl;
                    warnings = true;
                    Debug.Report(WarningCode::IncludingIndeterminationSymbols, new std::string[1] {seqsName[i]});
//                     cerr << "WARNING: Sequence \"" << seqsName[i] << "\" has some inde"
//                          << "termination symbols 'X' at the end of sequence. They will be"
//                          << " included in the final newAlignment." << endl;
                    break;
                }

                /* If nucleotide sequence is shorter than protein sequence, return an
                 * error since it is not feasible to cut the input protein aligment to
                 * fit it into CDNA sequences size */
                else {
                    if (!warnings)
                        cerr << endl;
                    warnings = true;
                    Debug.Report(WarningCode::LessNucleotidesThanExpected, new std::string[3] { seqsName[i], std::to_string(lengths[j]), std::to_string(seqLength)});
//                     cerr << "WARNING: Sequence \"" << seqsName[i] << "\" has less nucleo"
//                          << "tides (" << lengths[j] << ") than expected (" << seqLength
//                          << "). It will be added N's to complete the sequence"  << endl;
                    break;
                }
            }
        }

        /* Warn about a mismatch a sequences name level */
        if(j == totalInputSeqs) {
            Debug.Report(ErrorCode::SequenceNotPresentInCDS, new std::string[1] { seqsName[i] });
//             cerr << endl << "ERROR: Sequence \"" << seqsName[i] << "\" is not in "
//                  << "CDS file." << endl << endl;
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
                Debug.Report(ErrorCode::UnknownCharacter, new std::string[2] { seqsName[i], std::to_string(sequences[i][j]) });
//                 cerr << endl << "ERROR: The sequence \"" << seqsName[i] << "\" has an "
//                      << "unknown (" << sequences[i][j] << ") character." << endl;
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
        Debug.Report(ErrorCode::NotAligned, new std::string[1] { filename });
//         cerr << endl << "ERROR: Sequences should be aligned (all with same length) "
//              << "and there are not. Check your input alignment" << endl;
        return false;
    }

    /* Full-fill some information about input alignment */
    if(residNumber == 0)
        residNumber = sequences[0].length();

    /* Check whether aligned sequences have the length fixed for the input alig */
    for(i = 0; (i < sequenNumber) and (aligned); i++) {
        if(sequences[i].length() != residNumber) {
            Debug.Report(ErrorCode::SequencesNotSameSize, new std::string[3] { seqsName[i], std::to_string(sequences[i].length()), std::to_string(residNumber)});
//             cerr << endl << "ERROR: The sequence \"" << seqsName[i] << "\" ("
//                  << sequences[i].length() << ") does not have the same number of residues "
//                  << "fixed by the alignment (" << residNumber << ")." << endl;
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

// bool newAlignment::alignmentSummaryHTML(char *destFile, int residues, int seqs, \
//                                         int *selectedRes, int *selectedSeq, float *consValues) {
//
//     /* Generate an HTML file with a visual summary about which sequences/columns
//      * have been selected and which have not */
//
//     int i, j, k, kj, upper, minHTML, maxLongName, *gapsValues;
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
//         Debug.Report(ErrorCode::NotAligned, new std::string[1] { filename });
// //     cerr << endl << "ERROR: Sequences are not aligned." << endl << endl;
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
//     minHTML = utils::max(25, maxLongName + 10);
//
//     /* Initialize local variables to control which columns/sequences
//      * will be kept in the output alignment */
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
//     {
//         res[selectedRes[i]] = true;
//     }
//     for(i = 0; i < seqs; i++)
//     {
//
//         seq[selectedSeq[i]] = true;
//     }
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
//          << "    <style type=\"text/css\" media=\"all\">" << endl
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
//         file << endl << setw(minHTML) << left << "    Gaps Scores:        "
//              << "<span  class=c1>  =0=  </span><span  class=c2> <.001 </span>"
//              << "<span  class=c3> <.050 </span><span  class=c4> <.100 </span>"
//              << "<span  class=c5> <.150 </span><span  class=c6> <.200 </span>"
//              << "<span  class=c7> <.250 </span><span  class=c8> <.350 </span>"
//              << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
//              << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";
//
//     if (simValues != NULL)
//         file << endl << setw(minHTML) << left << "    Similarity Scores:  "
//              << "<span  class=c1>  =0=  </span><span  class=c2> <1e-6 </span>"
//              << "<span  class=c3> <1e-5 </span><span  class=c4> <1e-4 </span>"
//              << "<span  class=c5> <.001 </span><span  class=c6> <.010 </span>"
//              << "<span  class=c7> <.100 </span><span  class=c8> <.250 </span>"
//              << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
//              << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";
//
//     if (consValues != NULL)
//         file << endl << setw(minHTML) << left << "    Consistency Scores: "
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
//             HTMLBLOCKS) {
//
//         /* Print main columns number */
//         file << endl << setw(minHTML + 10) << right << (j + 10);
//         for(i = j + 20; ((i <= residNumber) && (i <= upper)); i += 10)
//             file << setw(10) << right << (i);
//
//         /* Print special characters to delimit sequences blocks */
//         file << endl << setw(minHTML + 1) << right;
//         for(i = j + 1; ((i <= residNumber) && (i <= upper)); i++)
//             file << (!(i % 10) ? "+" : "=");
//         file << endl;
//
//         /* Print sequences name */
//         for(i = 0; i < sequenNumber; i++) {
//             file << "    <span class=" << ((seq[i]) ? "sel>" : "nsel>") << seqsName[i]
//                  << "</span>" << setw(minHTML - 4 - seqsName[i].size()) << right << "";
//
//             /* Print residues corresponding to current sequences block */
//             for(k = j; ((k < residNumber) && (k < upper)); k++) {
//                 for(kj = 0, tmpColumn.clear(); kj < sequenNumber; kj++)
//                     tmpColumn += sequences[kj][k];
//                 /* Determine residue color based on residues across the alig column */
//                 type = utils::determineColor(sequences[i][k], tmpColumn);
//                 if (type == 'w')
//                     file << sequences[i][k];
//                 else
//                     file << "<span id=" << type << ">" << sequences[i][k] << "</span>";
//             }
//             file << endl;
//         }
//
//         file << endl << setw(minHTML) << left << "    Selected Cols:      ";
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
//             file << endl << setw(minHTML) << left << "    Gaps Scores:        ";
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
//             file << endl << setw(minHTML) << left << "    Similarity Scores:  ";
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
//             file << endl << setw(minHTML) << left << "    Consistency Scores: ";
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
//     delete [] seq;
//     delete [] res;
//
//     return true;
// }

bool newAlignment::alignmentSummaryHTML(newAlignment & _trimmedAlignment, char *destFile, float *consValues) {

    /* Generate an HTML file with a visual summary about which sequences/columns
     * have been selected and which have not */

    int i, j, k, kj, upper, minHTML, maxLongName, *gapsValues;
    string tmpColumn;
    float *simValues;
//     bool *res, *seq;
    ofstream file;
    char type;

    /* Allocate some local memory */
    tmpColumn.reserve(originalSequenNumber);

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!isAligned) {
        Debug.Report(ErrorCode::NotAligned, new std::string[1] { filename });
//     cerr << endl << "ERROR: Sequences are not aligned." << endl << endl;
        return false;
    }

    /* Open output file and check that file pointer is valid */
    file.open(destFile);
    if(!file)
        return false;

    /* Compute maximum sequences name length. */
    maxLongName = 0;
    for(i = 0; i < originalSequenNumber; i++)
        maxLongName = utils::max(maxLongName, seqsName[i].size());

    /* Compute HTML blank spaces */
    minHTML = utils::max(25, maxLongName + 10);

    /* Initialize local variables to control which columns/sequences
     * will be kept in the output alignment */
//     res = new bool[originalResidNumber];
//     for(i = 0; i < originalResidNumber; i++)
//         res[i] = false;
//
//     seq = new bool[originalSequenNumber];
//     for(i = 0; i < originalSequenNumber; i++)
//         seq[i] = false;

    /* Record which columns/sequences from original alignment
     * have been kept in the final one */
//     for(i = 0; i < originalResidNumber; i++)
//     {
//         res[selectedRes[i]] = true;
//     }
//     for(i = 0; i < seqs; i++)
//     {
//
//         seq[selectedSeq[i]] = true;
//     }

    /* Recover some stats about different scores from current alignment */
    gapsValues = NULL;
    if (sgaps != NULL)
        gapsValues = sgaps -> getGapsWindow();
    simValues = NULL;
    if (scons != NULL)
        simValues = scons -> getMdkwVector();

    /* Print HTML header into output file */
    file << "<!DOCTYPE html>" << endl << "<html><head>" << endl << "    <meta "
         << "http-equiv=\"Content-Type\" content=\"text/html;charset=ISO-8859-1\" />"
         << endl << "    <title>trimAl v1.4 Summary</title>" << endl
         << "    <style type=\"text/css\" media=\"all\">" << endl

         << "    #b  { background-color: #3366ff; }\n"
         << "    #r  { background-color: #cc0000; }\n"
         << "    #g  { background-color: #33cc00; }\n"
         << "    #p  { background-color: #ff6666; }\n"
         << "    #m  { background-color: #cc33cc; }\n"
         << "    #o  { background-color: #ff9900; }\n"
         << "    #c  { background-color: #46C7C7; }\n"
         << "    #y  { background-color: #FFFF00; }\n"

         << "    .sel  { background-color: #B9B9B9; }\n"
         << "    .nsel { background-color: #E9E9E9; }\n"

         /* Sets of colors for high-lighting scores intervals */
         << "    .c1   { background-color: #FFFBF2; }\n"
         << "    .c2   { background-color: #FFF8CC; }\n"
         << "    .c3   { background-color: #FAF0BE; }\n"
         << "    .c4   { background-color: #F0EAD6; }\n"
         << "    .c5   { background-color: #F3E5AB; }\n"
         << "    .c6   { background-color: #F4C430; }\n"
         << "    .c7   { background-color: #C2B280; color: white; }\n"
         << "    .c8   { background-color: #DAA520; color: white; }\n"
         << "    .c9   { background-color: #B8860B; color: white; }\n"
         << "    .c10  { background-color: #918151; color: white; }\n"
         << "    .c11  { background-color: #967117; color: white; }\n"
         << "    .c12  { background-color: #6E5411; color: white; }\n"

         /* Other HTML elements */
         << "    </style>\n  </head>\n\n" << "  <body>\n" << "  <pre>" << endl;

    /* Show information about how many sequences/residues have been selected */
    file << "    <span class=sel>Selected Sequences: " << setw(5) << right << _trimmedAlignment.sequenNumber
         <<" /Selected Residues: " << setw(7) << right << _trimmedAlignment.residNumber << "</span>"
         << endl << "    <span class=nsel>Deleted Sequences:  " << setw(5) << right
         << sequenNumber - _trimmedAlignment.sequenNumber << " /Deleted Residues:  " << setw(7) << right
         << residNumber - _trimmedAlignment.residNumber << "</span>" << endl;

    /* Print headers for different scores derived from input alignment/s */
    if (gapsValues != NULL)
        file << endl << setw(minHTML) << left << "    Gaps Scores:        "
             << "<span  class=c1>  =0=  </span><span  class=c2> <.001 </span>"
             << "<span  class=c3> <.050 </span><span  class=c4> <.100 </span>"
             << "<span  class=c5> <.150 </span><span  class=c6> <.200 </span>"
             << "<span  class=c7> <.250 </span><span  class=c8> <.350 </span>"
             << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
             << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";

    if (simValues != NULL)
        file << endl << setw(minHTML) << left << "    Similarity Scores:  "
             << "<span  class=c1>  =0=  </span><span  class=c2> <1e-6 </span>"
             << "<span  class=c3> <1e-5 </span><span  class=c4> <1e-4 </span>"
             << "<span  class=c5> <.001 </span><span  class=c6> <.010 </span>"
             << "<span  class=c7> <.100 </span><span  class=c8> <.250 </span>"
             << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
             << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";

    if (consValues != NULL)
        file << endl << setw(minHTML) << left << "    Consistency Scores: "
             << "<span  class=c1>  =0=  </span><span  class=c2> <.001 </span>"
             << "<span  class=c3> <.050 </span><span  class=c4> <.100 </span>"
             << "<span  class=c5> <.150 </span><span  class=c6> <.200 </span>"
             << "<span  class=c7> <.250 </span><span  class=c8> <.350 </span>"
             << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
             << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";

    if ((gapsValues != NULL) or (simValues == NULL) or (consValues == NULL))
        file << endl;

    /* Print Sequences in block of BLOCK_SIZE */
    for(j = 0, upper = HTMLBLOCKS; j < originalResidNumber; j += HTMLBLOCKS, upper += \
            HTMLBLOCKS) {

        /* Print main columns number */
        file << endl << setw(minHTML + 10) << right << (j + 10);
        for(i = j + 20; ((i <= originalResidNumber) && (i <= upper)); i += 10)
            file << setw(10) << right << (i);

        /* Print special characters to delimit sequences blocks */
        file << endl << setw(minHTML + 1) << right;
        for(i = j + 1; ((i <= originalResidNumber) && (i <= upper)); i++)
            file << (!(i % 10) ? "+" : "=");
        file << endl;

        /* Print sequences name */
        for(i = 0; i < originalSequenNumber; i++) {
            file << "    <span class=" << ((_trimmedAlignment.saveSequences[i] != -1 ) ? "sel>" : "nsel>") << seqsName[i]
                 << "</span>" << setw(minHTML - 4 - seqsName[i].size()) << right << "";

            /* Print residues corresponding to current sequences block */
            for(k = j; ((k < originalResidNumber) && (k < upper)); k++) {
                for(kj = 0, tmpColumn.clear(); kj < originalSequenNumber; kj++)
                    tmpColumn += sequences[kj][k];
                /* Determine residue color based on residues across the alig column */
                type = utils::determineColor(sequences[i][k], tmpColumn);

                if (type == 'w')
                    file << sequences[i][k];
                else
                    file << "<span id=" << type << ">" << sequences[i][k] << "</span>";
            }

            file << endl;
        }

        file << endl << setw(minHTML) << left << "    Selected Cols:      ";
        for(k = j; ((k < originalResidNumber) && (k < (j + HTMLBLOCKS))); k++)
            file << "<span class=" << (_trimmedAlignment.saveResidues[k] != -1 ? "sel" : "nsel") << "> </span>";
        file << endl;

        /* If there is not any score to print, skip this part of the function */
        if ((gapsValues == NULL) and (simValues == NULL) and (consValues == NULL))
            continue;

        /* Print score colors according to certain predefined thresholds */
        if (gapsValues != NULL) {
            file << endl << setw(minHTML) << left << "    Gaps Scores:        ";
            for(k = j; ((k < originalResidNumber) && (k < (j + HTMLBLOCKS))); k++)
                if(gapsValues[k] == 0)
                    file << "<span class=c12> </span>";
                else if(gapsValues[k] == originalSequenNumber)
                    file << "<span class=c1> </span>";
                else if(1 - (float(gapsValues[k])/originalSequenNumber) >= .750)
                    file << "<span class=c11> </span>";
                else if(1 - (float(gapsValues[k])/originalSequenNumber) >= .500)
                    file << "<span class=c10> </span>";
                else if(1 - (float(gapsValues[k])/originalSequenNumber) >= .350)
                    file << "<span  class=c9> </span>";
                else if(1 - (float(gapsValues[k])/originalSequenNumber) >= .250)
                    file << "<span  class=c8> </span>";
                else if(1 - (float(gapsValues[k])/originalSequenNumber) >= .200)
                    file << "<span  class=c7> </span>";
                else if(1 - (float(gapsValues[k])/originalSequenNumber) >= .150)
                    file << "<span  class=c6> </span>";
                else if(1 - (float(gapsValues[k])/originalSequenNumber) >= .100)
                    file << "<span  class=c5> </span>";
                else if(1 - (float(gapsValues[k])/originalSequenNumber) >= .050)
                    file << "<span  class=c4> </span>";
                else if(1 - (float(gapsValues[k])/originalSequenNumber) >= .001)
                    file << "<span  class=c3> </span>";
                else
                    file << "<span  class=c2> </span>";
        }
        if (simValues != NULL) {
            file << endl << setw(minHTML) << left << "    Similarity Scores:  ";
            for(k = j; ((k < originalResidNumber) && (k < (j + HTMLBLOCKS))); k++)
                if(simValues[k] == 1)
                    file << "<span class=c12> </span>";
                else if(simValues[k] == 0)
                    file << "<span class=c1> </span>";
                else if(simValues[k] >= .750)
                    file << "<span class=c11> </span>";
                else if(simValues[k] >= .500)
                    file << "<span class=c10> </span>";
                else if(simValues[k] >= .250)
                    file << "<span  class=c9> </span>";
                else if(simValues[k] >= .100)
                    file << "<span  class=c8> </span>";
                else if(simValues[k] >= .010)
                    file << "<span  class=c7> </span>";
                else if(simValues[k] >= .001)
                    file << "<span  class=c6> </span>";
                else if(simValues[k] >= 1e-4)
                    file << "<span  class=c5> </span>";
                else if(simValues[k] >= 1e-5)
                    file << "<span  class=c4> </span>";
                else if(simValues[k] >= 1e-6)
                    file << "<span  class=c3> </span>";
                else
                    file << "<span  class=c2> </span>";
        }
        if (consValues != NULL) {
            file << endl << setw(minHTML) << left << "    Consistency Scores: ";
            for(k = j; ((k < originalResidNumber) && (k < (j + HTMLBLOCKS))); k++)
                if(consValues[k] == 1)
                    file << "<span class=c12> </span>";
                else if(consValues[k] == 0)
                    file << "<span class=c1> </span>";
                else if(consValues[k] >= .750)
                    file << "<span class=c11> </span>";
                else if(consValues[k] >= .500)
                    file << "<span class=c10> </span>";
                else if(consValues[k] >= .350)
                    file << "<span  class=c9> </span>";
                else if(consValues[k] >= .250)
                    file << "<span  class=c8> </span>";
                else if(consValues[k] >= .200)
                    file << "<span  class=c7> </span>";
                else if(consValues[k] >= .150)
                    file << "<span  class=c6> </span>";
                else if(consValues[k] >= .100)
                    file << "<span  class=c5> </span>";
                else if(consValues[k] >= .050)
                    file << "<span  class=c4> </span>";
                else if(consValues[k] >= .001)
                    file << "<span  class=c3> </span>";
                else
                    file << "<span  class=c2> </span>";
        }
        file << endl;
    }

    /* Print HTML footer into output file */
    file << "    </pre>" << endl << "  </body>" << endl << "</html>" << endl;

    /* Close output file and deallocate local memory */
    file.close();

    return true;
}

bool newAlignment::alignmentSummarySVG(newAlignment & _trimmedAlignment, char *destFile, float *consValues, int blocks) {

    int i, j, k, counter = 0;

    // BEGIN Init variables
    int * gapsValues = NULL;
    if (sgaps != NULL)
        gapsValues = sgaps -> getGapsWindow();

    float * simValues = NULL;
    if (scons != NULL)
        simValues = scons -> getMdkwVector();

    // Check if alignment is aligned;
    if (!isAligned) {
        Debug.Report(ErrorCode::NotAligned, new std::string[1] { filename });
        return false;
    }

    // Calculate the blockSize;
    blocks = 120;

    int fontSize = 15;

    // Allocate some local memory
    char * tmpColumn = new char[originalSequenNumber + 1];

    // Open the file;
    ofstream file;
    file.open(destFile);
    if(!file)
        return false;

    /* Compute HTML blank spaces */
    j = 0;
    for(i = 0; i < sequenNumber; i++)
        j = utils::max(j, seqsName[i].size());

    int sequencesNamesLength = utils::max(25, j);

    // Init Colors
    std::map<char, string> mappedColors = std::map<char, string>
    {
        {'o', "orange"},    {'y', "yellow"},        {'b', "royalblue"},
        {'w', "lightgrey"}, {'p', "darkviolet"},    {'r', "red"},
        {'g', "lime"},      {'m', "magenta"},       {'c', "aqua"}
    };

    std::map<char, string> mappedColorsMeaning = std::map<char, string>
    {
        {'o', "Glycines"},      {'y', "Prolines"},          {'b', "Hidrophobic"},
        {'w', "Unconserved"},   {'p', "Prolines"},          {'r', "Positive Charge"},
        {'g', "Polar"},         {'m', "Negative Charge"},   {'c', "Aromatic"}
    };

    // END
    
    // BEGIN Init SVG
    float leftMargin          = 5,
          nameBlocksMargin    = 15,
          rightMargin         = 5,
          topMargin           = 5,
          bottomMargin        = 10,
          interBlocksMargin   = 20;

    float height = ((originalSequenNumber + 3 + (gapsValues != NULL || simValues != NULL ? 7 : 0) /* Colnumbering occupies 3 rows. Stats 5 */)
                    * fontSize + interBlocksMargin)                 /* Height of a block + interblock */
                   * std::ceil((float)originalResidNumber/blocks)  /* Number of blocks */
                   - interBlocksMargin /* Number of interblocks is number of blocks  - 1 */
                   + fontSize * 16
                   + topMargin
                   + bottomMargin;

    float width = sequencesNamesLength * fontSize
                  + blocks * fontSize * 0.75F
                  + leftMargin
                  + nameBlocksMargin
                  + rightMargin;

    float currentHeight = topMargin + fontSize;

    // Start the svg output
    file << "<svg \
        xmlns=\"http://www.w3.org/2000/svg\" \
        xmlns:xlink=\"http://www.w3.org/1999/xlink\" \
        version=\"1.1\" \
        height=\""  << height << "\" \
        width=\""   << width << "px\">" << endl;
    // END
        
    // BEGIN Functions
        file << "<script type=\"text/ecmascript\"> \n\
                <![CDATA[ \n\
                function onMouseOverLabel(evt) { \n\
                    var Line = document.childNodes[0].getElementById(evt.target.id.replace(\"Label\", \"Line\")); \n\
                    Line.setAttribute(\"style\", \"fill:none;stroke:black;stroke-width:3\"); \n\
                    evt.target.setAttribute(\"font-size\", 17); \n\
                } \n\
                function onMouseOutLabel(evt) { \n\
                    var Line = document.childNodes[0].getElementById(evt.target.id.replace(\"Label\", \"Line\")); \n\
                    Line.setAttribute(\"style\", \"fill:none;stroke:black;stroke-width:1\"); \n\
                    evt.target.setAttribute(\"font-size\", 15); \n\
                } \n\
                function onMouseOverLine(evt) { \n\
                    var Label = document.childNodes[0].getElementById(evt.target.id.replace(\"Line\", \"Label\")); \n\
                    evt.target.setAttribute(\"style\", \"fill:none;stroke:black;stroke-width:3\"); \n\
                    Label.setAttribute(\"font-size\", 17); \n\
                } \n\
                function onMouseOutLine(evt) { \n\
                    var Label = document.childNodes[0].getElementById(evt.target.id.replace(\"Line\", \"Label\")); \n\
                    evt.target.setAttribute(\"style\", \"fill:none;stroke:black;stroke-width:1\"); \n\
                    Label.setAttribute(\"font-size\", 15); \n\
                } \n\
                ]]> \n </script>" << endl;
    // END
        
    // BEGIN Patterns
    
    // Selected Block
    file << "<pattern id=\"selected\" patternUnits=\"userSpaceOnUse\" width=\"10\" height=\"10\"> " << endl;
    file << "<rect width=\"10\" height=\"10\" fill=\"green\" opacity=\"0.5\"/>";
    file << "<path d=\"M-1,1 l2,-2 M0,10 l10,-10 M9,11 l2,-2\" stroke=\"white\" stroke-width=\"1\"/>";
    file << "</pattern> " << endl;

//     file << "<pattern id=\"selected-focus\" patternUnits=\"userSpaceOnUse\" width=\"10\" height=\"10\"> " << endl;
//     file << "<rect width=\"10\" height=\"10\" fill=\"mediumseagreen\"/>";
//     file << "</pattern> " << endl;

    // Deleted Block
    file << "<pattern id=\"rejected\" patternUnits=\"userSpaceOnUse\" width=\"10\" height=\"10\"> " << endl;
    file << "<rect width=\"10\" height=\"10\" fill=\"red\" opacity=\"0.5\"/>";
    file << "<path d=\"M2,-2 l-1,1 M10,-10 l0,10 M2,-2 l9,11\" stroke=\"white\" stroke-width=\"1\"/>";
    file << "</pattern> " << endl;

//     file << "<pattern id=\"deleted-focus\" patternUnits=\"userSpaceOnUse\" width=\"10\" height=\"10\"> " << endl;
//     file << "<rect width=\"10\" height=\"10\" fill=\"black\"/>";
//     file << "<path d=\"M-1,1 l2,-2 M0,10 l10,-10 M9,11 l2,-2\" stroke=\"red\" stroke-width=\"2\"/>";
//     file << "</pattern> " << endl;
        
    // END
    
    // BEGIN Header
        
    strtok(&filename[0], " ");
    std::string fname = strtok(NULL, " ");
    
    file << "<rect \
                style=\"fill:none;stroke:lightgrey\" \
                height=\""<< fontSize * 6.5F <<"\" \
                width=\"" << fontSize * std::max(30.5F, fname.length() * 0.75F + 1.F) << "px\" \
                x =\"" << leftMargin + nameBlocksMargin - fontSize * 0.5F << "px\"\
                y =\""<< currentHeight - fontSize * 0.5F <<"\"/>" << endl;

    file << "<text \
                font-family = \"monospace\" \
                font-size=\"" << fontSize << "px\" \
                dy=\".35em\" \
                x =\"" << leftMargin + nameBlocksMargin << "\" \
                y=\"" << (currentHeight + fontSize / 2) << "\" \
                kerning=\"0\" \
                text-anchor=\"start\" \
                lengthAdjust=\"spacingAndGlyphs\"\
                textLength=\"" << fname.size() * 0.75F * fontSize << "\" \
                style=\"font-weight:100\" \
                xml:space=\"preserve\" >"

         << fname

         << "</text>" << endl;
    currentHeight += fontSize * 2;
    fname = "Selected Sequences \t" + std::to_string(_trimmedAlignment.sequenNumber) + " / " + std::to_string(originalSequenNumber);
    file << "<text \
                font-family = \"monospace\" \
                font-size=\"" << fontSize << "px\" \
                dy=\".35em\" \
                x =\"" << leftMargin + nameBlocksMargin << "\" \
                y=\"" << (currentHeight + fontSize / 2) << "\" \
                kerning=\"0\" \
                text-anchor=\"start\" \
                lengthAdjust=\"spacingAndGlyphs\"\
                textLength=\"" << std::min((float)fname.size() * 0.65F, (float)sequencesNamesLength) * fontSize << "\" \
                style=\"font-weight:100\" \
                xml:space=\"preserve\" >"

         << fname

         << "</text>" << endl;
    currentHeight += fontSize * 2;
    fname = "Selected Residues \t\t" + std::to_string(_trimmedAlignment.residNumber) + " / " + std::to_string(originalResidNumber);
    file << "<text \
                font-family = \"monospace\" \
                font-size=\"" << fontSize << "px\" \
                dy=\".35em\" \
                x =\"" << leftMargin + nameBlocksMargin << "\" \
                y=\"" << (currentHeight + fontSize / 2) << "\" \
                kerning=\"0\" \
                text-anchor=\"start\" \
                lengthAdjust=\"spacingAndGlyphs\"\
                textLength=\"" << std::min((float)fname.size() * 0.65F, (float)sequencesNamesLength) * fontSize << "\" \
                style=\"font-weight:100\" \
                xml:space=\"preserve\" >"

         << fname

         << "</text>" << endl;
    currentHeight += fontSize * 3;

    // END
    
    // BEGIN Legend
    file << "<rect \
                style=\"fill:none;stroke:lightgrey\" \
                height=\""<< fontSize * 6.5F <<"\" \
                width=\"" << fontSize * 30.5 << "px\" \
                x =\"" << leftMargin + nameBlocksMargin - fontSize * 0.5F << "px\"\
                y =\""<< currentHeight - fontSize * 0.5F <<"\"/>" << endl;
                    
    std::map<char,string>::iterator it = mappedColors.begin();
    for (i = 0 ; i < 3; i++)
    {
        for (j = 0; j < 3; j++, it++)
        {
            file << "<rect \
                    style=\"fill:"<< it->second <<";stroke:black;stroke-width:1px\" \
                    height=\""<< fontSize <<"\" \
                    width=\"" << fontSize << "px\" \
                    x =\"" << leftMargin + nameBlocksMargin * (j + 1) + (9 * fontSize * j) << "px\"\
                    y =\""<< currentHeight + fontSize * i * 1.2F <<"\"/>" << endl;

            file << "<text \
                font-family = \"monospace\" \
                font-size=\"" << fontSize * 0.75F << "px\" \
                dy=\".35em\" \
                x =\"" << leftMargin + nameBlocksMargin * (j + 3) + (9 * fontSize * j) << "\" \
                y =\"" << currentHeight + fontSize * i * 1.2F + fontSize * 0.5F << "\" \
                kerning=\"0\" \
                text-anchor=\"start\" \
                lengthAdjust=\"spacingAndGlyphs\"\
                textLength=\"" << mappedColorsMeaning[it->first].size() * 0.5F * fontSize << "\" \
                style=\"font-weight:100\" \
                xml:space=\"preserve\" >"

                 << mappedColorsMeaning[it->first]

                 << "</text>" << endl;
        }
        currentHeight += fontSize;
    }
    currentHeight += fontSize * 4;
    // END

    // Foreach block
    for (j = 0; j < originalResidNumber; j+= blocks)
    {
        // BEGIN Columms Numbering Numbers
        file << "<text \
                font-family = \"monospace\" \
                font-size=\"" << fontSize << "px\" \
                dy=\".35em\" \
                x =\"" << leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize << "\" \
                y=\"" << (currentHeight + fontSize / 2) << "\" \
                kerning=\"0\" \
                text-anchor=\"start\" \
                lengthAdjust=\"spacingAndGlyphs\"\
                textLength=\"" << std::min((int)blocks, originalResidNumber - j) * fontSize * 0.75F << "\" \
                style=\"font-weight:100\" \
                xml:space=\"preserve\" >";
        for (k = j; k < originalResidNumber && (k - j) < blocks; k+= 10)
        {
            file << std::left << std::setw(10) << std::setfill(' ') << k;
        }

        file << "</text>" << endl;

        currentHeight += fontSize;
        // END
        
        // BEGIN Columns Numbering Symbols

        file << "<text \
                font-family = \"monospace\" \
                font-size=\"" << fontSize << "px\" \
                dy=\".35em\" \
                x =\"" << leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize << "\" \
                y=\"" << (currentHeight + fontSize / 2) << "\" \
                kerning=\"0\" \
                text-anchor=\"start\" \
                lengthAdjust=\"spacingAndGlyphs\"\
                textLength=\"" << std::min((int)blocks, originalResidNumber - j) * fontSize * 0.75F << "\" \
                style=\"font-weight:100\" \
                xml:space=\"preserve\" >";
        for (k = j; k < originalResidNumber && (k - j) < blocks; k+= 10)
        {
            file << std::left << std::setw(10) << std::setfill(' ') << "|";
        }

        file << "</text>" << endl;

        file << "<line x1=\"" << leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize
             << "\" x2=\"" << leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize
             + std::min((int)blocks, originalResidNumber - j) * fontSize * 0.75F
             << "\" y1=\"" << currentHeight + fontSize / 2
             << "\" y2=\"" << currentHeight + fontSize / 2
             << "\" style=\"stroke:rgb(0,0,0);stroke-width:2\"/>" << endl;
        currentHeight += fontSize ;


        currentHeight += fontSize ;

        // END

        // BEGIN Color Matrix
        char lastColor, newColor;
        int lastPos = 0;
        for (i = 0; i < blocks && (i + j) < originalResidNumber; i++)
        {
            for (k = 0; k < originalSequenNumber; k ++)
            {
                tmpColumn[k]=(sequences[k][i + j]);
            }
            lastColor = utils::determineColor(tmpColumn[0], tmpColumn);
            lastPos = 0;
            for (k = 1; k < originalSequenNumber; k ++)
            {
                newColor = utils::determineColor(tmpColumn[k], tmpColumn);
                if (newColor != lastColor)
                {
                    if (lastColor != 'w')
                        file << "<rect "        <<
                             " style=\"fill:"   << mappedColors[lastColor]      <<"\" " <<
                             " height=\""       << fontSize * (k - lastPos)     <<"\" " <<
                             " width=\""        << fontSize * 0.75F << "px\" "  <<
                             " x =\""           <<

                             leftMargin +
                             nameBlocksMargin +
                             sequencesNamesLength * fontSize +
                             i * 0.75F * fontSize

                             << "\" "           <<
                             " y =\""           << currentHeight + fontSize * lastPos   << "\" />" << endl;
                    lastColor = newColor;
                    lastPos = k;
                }
            }
            if (lastColor != 'w')
                file << "<rect " <<
                     " style=\"fill:"<< mappedColors[lastColor] <<"\" " <<
                     " height=\""<< fontSize * (k - lastPos) <<"\" " <<
                     " width=\"" << fontSize * 0.75F << "px\" " <<
                     " x =\"" <<

                     leftMargin + nameBlocksMargin + (sequencesNamesLength) * fontSize + i * 0.75F * fontSize

                     << "\" " <<
                     " y =\"" << currentHeight + fontSize * lastPos << "\" />" << endl;
        }

        // END

        // BEGIN Sequences Names and Sequences

        for (i = 0; i < originalSequenNumber; i ++)
        {
            // Print names
            file << "<text \
                font-family = \"monospace\" \
                font-size=\"" << fontSize << "px\" \
                dy=\".35em\" \
                x =\"" << leftMargin << "\" \
                y=\"" << (currentHeight + fontSize / 2) << "\" \
                kerning=\"0\" \
                text-anchor=\"start\" \
                textLength=\"" << std::min((float)seqsName[i].size() * 0.75F, (float)sequencesNamesLength) * fontSize << "\" \
                style=\"font-weight:100\"" << ">"

                 << seqsName[i].substr(0, sequencesNamesLength)

                 << "</text>" << endl;
            // Print sequences
            file << "<text \
                font-family = \"monospace\" \
                font-size=\"" << fontSize << "px\" \
                dy=\".35em\" \
                x =\"" << leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize << "\" \
                y=\"" << (currentHeight + fontSize / 2) << "\" \
                kerning=\"0\" \
                text-anchor=\"start\" \
                lengthAdjust=\"spacingAndGlyphs\" \
                textLength=\"" << std::min((int)blocks, originalResidNumber - j) * fontSize * 0.75F << "\" \
                style=\"font-weight:100\"" << ">"

                 << sequences[i].substr(j, blocks)

                 << "</text>" << endl;
            if (_trimmedAlignment.saveSequences[i] == -1)
            {
                file << "<line x1=\""<< leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize <<"\"" << 
                             " x2=\""<< leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize + std::min(blocks, originalResidNumber - j) * fontSize * 0.75F <<"\"" << 
                             " y1=\""<< (currentHeight + fontSize / 2) <<"\"" <<
                             " y2=\""<< (currentHeight + fontSize / 2) <<"\"" << 
                             " style=\"stroke:rgb(0,0,0);stroke-width:2\""
                             "/>" << endl;
                file << "<line x1=\""<< leftMargin <<"\"" << 
                             " x2=\""<< leftMargin + std::min((float)seqsName[i].size() * 0.75F, (float)sequencesNamesLength) * fontSize <<"\"" << 
                             " y1=\""<< (currentHeight + fontSize / 2) <<"\"" <<
                             " y2=\""<< (currentHeight + fontSize / 2) <<"\"" << 
                             " style=\"stroke:rgb(0,0,0);stroke-width:2\""
                             "/>" << endl;
            }
            currentHeight += fontSize;
        }

        // END

        // BEGIN Selected / Rejected Boxes
            currentHeight += fontSize;

            // Big white box
            file << "<rect " <<
                 " style=\"fill:none;stroke:black\"" <<
                 " height=\""<< fontSize * 5 <<"\" " <<
                 " width=\"" << fontSize * 0.75F * std::min(blocks, originalResidNumber - j) << "px\" " <<
                 " x =\"" <<

                 leftMargin + nameBlocksMargin + (sequencesNamesLength) * fontSize

                 << "\" " <<
                 " y =\"" << currentHeight << "\" />" << endl;
            
            bool rejected = _trimmedAlignment.saveResidues[j] == -1;
            lastPos = j;
            for (k = j + 1; k < originalResidNumber && (k - j) < blocks; k++)
            {
                if ((_trimmedAlignment.saveResidues[k] == -1) != rejected)
                {
                    file << "<rect " <<
                            " style=\"fill:url(#" << (rejected ? "rejected" : "selected") << ");stroke:black\"" <<
                            " height=\""<< fontSize * 5 <<"\" " <<
                            " width=\"" << fontSize * 0.75F * (k - lastPos) << "px\" " <<
                            " opacity=\"0.5\"" <<

                            " x =\"" <<

                            leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize + fontSize * 0.75F * (lastPos - j) 

                        << "\" " <<
                        " y =\"" << currentHeight << "\" />" << endl;
                    rejected = !rejected;
                    lastPos = k;
                }
            }
            file << "<rect " <<
                    " style=\"fill:url(#" << (rejected ? "rejected" : "selected") << ");stroke:black\"" <<
                    " height=\""<< fontSize * 5 <<"\" " <<
                    " width=\"" << fontSize * 0.75F * (k - lastPos) << "px\" " <<
                    " opacity=\"0.5\""
                    " x =\"" <<

                    leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize + fontSize * 0.75F * (lastPos - j) 

                    << "\" " <<
                    " y =\"" << currentHeight << "\" />" << endl;
        // END
                    
        // BEGIN Stats report
            if (gapsValues || simValues || consValues)
            {
                file << "<text \
                        font-family = \"monospace\" \
                        font-size=\"" << fontSize << "px\" \
                        dy=\".35em\" \
                        x =\"" << leftMargin + (sequencesNamesLength) * fontSize << "\" \
                        y=\"" << currentHeight << "\" \
                        kerning=\"0\" \
                        text-anchor=\"end\" \
                        lengthAdjust=\"spacingAndGlyphs\"\
                        style=\"font-weight:100\" \
                        xml:space=\"preserve\" >"

                    << "1"

                    << "</text>" << endl;
                    
                file << "<text \
                        font-family = \"monospace\" \
                        font-size=\"" << fontSize << "px\" \
                        dy=\".35em\" \
                        x =\"" << leftMargin + (sequencesNamesLength) * fontSize << "\" \
                        y=\"" << currentHeight + fontSize * 5 << "\" \
                        kerning=\"0\" \
                        text-anchor=\"end\" \
                        lengthAdjust=\"spacingAndGlyphs\"\
                        style=\"font-weight:100\" \
                        xml:space=\"preserve\" >"

                    << "0"

                    << "</text>" << endl;
            }
            if (gapsValues)
            {
                file << "<polyline \
                        style=\"fill:none;stroke:black;stroke-width:1\" \
                        stroke-dasharray=\"10,4\" \
                        id =\"GapsLine" << counter << "\" \
                        onmouseover=\"onMouseOverLine(evt)\" \
                        onmouseout=\"onMouseOutLine(evt)\" \
                        points=\"";

                for (i = 0; i < blocks && (i+j) < originalResidNumber; i++)
                {
                    file << leftMargin + nameBlocksMargin + (sequencesNamesLength) * fontSize + (i + 0.5F) * 0.75F * fontSize << ","
                         << currentHeight + 10 + ( 1.F - (gapsValues[i + j] / (float)originalSequenNumber)) * fontSize * 4 << " ";
                }

                file << "\"/>" << endl;
                
                file << "<text \
                        font-family = \"monospace\" \
                        font-size=\"" << fontSize << "px\" \
                        dy=\".35em\" \
                        x =\"" << leftMargin << "\" \
                        y=\"" << (currentHeight + fontSize / 2) << "\" \
                        kerning=\"0\" \
                        text-anchor=\"start\" \
                        id =\"GapsLabel" << counter << "\" \
                        onmouseover=\"onMouseOverLabel(evt)\" \
                        onmouseout=\"onMouseOutLabel(evt)\" \
                        textLength=\"" << std::min(11 * 0.75F, (float)sequencesNamesLength) * fontSize << "\" \
                        style=\"font-weight:bold\"" << ">"

                        << "Gaps Values"

                        << "</text>" << endl;
                        
                   file << "<line x1=\"" << leftMargin 
                        << "\" x2=\"" << leftMargin + std::min(11 * 0.75F, (float)sequencesNamesLength) * fontSize
                        << "\" y1=\"" << currentHeight + fontSize * 1.25F
                        << "\" y2=\"" << currentHeight + fontSize * 1.25F
                        << "\" style=\"stroke:rgb(0,0,0);stroke-width:1\""
                        << " stroke-dasharray=\"10,4\""
                        << " id =\"GapsSubLine" << counter << "\""
                        << "/>" << endl;
            }

            if (simValues)
            {
                file << "<polyline \
                style=\"fill:none;stroke:black;stroke-width:1\" \
                stroke-dasharray=\"5,1\" \
                id =\"SimLine" << counter << "\" \
                onmouseover=\"onMouseOverLine(evt)\" \
                onmouseout=\"onMouseOutLine(evt)\" \
                points=\"";

                for (i = 0; i < blocks && (i+j) < originalResidNumber; i++)
                {
                    file << leftMargin + nameBlocksMargin + (sequencesNamesLength) * fontSize + (i + 0.5F) * 0.75F * fontSize << ","
                         << currentHeight + 10 + ( (1.F - simValues[i + j]) * fontSize * 4) << " ";
                }

                file << "\"/>";
                
                file << "<text \
                        font-family = \"monospace\" \
                        font-size=\"" << fontSize << "px\" \
                        dy=\".35em\" \
                        x =\"" << leftMargin << "\" \
                        y=\"" << (currentHeight + 2.5F * fontSize) << "\" \
                        kerning=\"0\" \
                        text-anchor=\"start\" \
                        id =\"SimLabel" << counter << "\" \
                        onmouseover=\"onMouseOverLabel(evt)\" \
                        onmouseout=\"onMouseOutLabel(evt)\" \
                        textLength=\"" << std::min(10 * 0.75F, (float)sequencesNamesLength) * fontSize << "\" \
                        style=\"font-weight:bold\"" << ">"

                        << "Sim Values"

                        << "</text>" << endl;
                        
                file << "<line x1=\"" << leftMargin 
                        << "\" x2=\"" << leftMargin + std::min(10 * 0.75F, (float)sequencesNamesLength) * fontSize
                        << "\" y1=\"" << currentHeight + fontSize * 3.25F
                        << "\" y2=\"" << currentHeight + fontSize * 3.25F
                        << "\" style=\"stroke:rgb(0,0,0);stroke-width:1\""
                        << " stroke-dasharray=\"5,1\""
                        << " id =\"SimSubLine" << counter << "\""
                        << "/>" << endl;
            }
            
            if (consValues)
            {
                file << "<polyline \
                style=\"fill:none;stroke:black;stroke-width:1\" \
                stroke-dasharray=\"2,2\" \
                id =\"ConsLine" << counter << "\" \
                onmouseover=\"onMouseOverLine(evt)\" \
                onmouseout=\"onMouseOutLine(evt)\" \
                points=\"";

                for (i = 0; i < blocks && (i+j) < originalResidNumber; i++)
                {
                    file << leftMargin + nameBlocksMargin + (sequencesNamesLength) * fontSize + (i + 0.5F) * 0.75F * fontSize << ","
                         << currentHeight + 10 + ( (1.F - consValues[i + j]) * fontSize * 4) << " ";
                }

                file << "\"/>";
                
                file << "<text \
                        font-family = \"monospace\" \
                        font-size=\"" << fontSize << "px\" \
                        dy=\".35em\" \
                        x =\"" << leftMargin << "\" \
                        y=\"" << (currentHeight + 4.5F * fontSize) << "\" \
                        kerning=\"0\" \
                        text-anchor=\"start\" \
                        id =\"ConsLabel" << counter << "\" \
                        onmouseover=\"onMouseOverLabel(evt)\" \
                        onmouseout=\"onMouseOutLabel(evt)\" \
                        textLength=\"" << std::min(10 * 0.75F, (float)sequencesNamesLength) * fontSize << "\" \
                        style=\"font-weight:bold\"" << ">"

                        << "Cons Values"

                        << "</text>" << endl;
                        
                file << "<line x1=\"" << leftMargin 
                        << "\" x2=\"" << leftMargin + std::min(10 * 0.75F, (float)sequencesNamesLength) * fontSize
                        << "\" y1=\"" << currentHeight + fontSize * 5.25F
                        << "\" y2=\"" << currentHeight + fontSize * 5.25F
                        << "\" style=\"stroke:rgb(0,0,0);stroke-width:1\""
                        << " stroke-dasharray=\"2,2\""
                        << " id =\"ConsSubLine" << counter << "\""
                        << "/>" << endl;
            }
            counter++;
            currentHeight += 6 * fontSize;
        
        // END

        currentHeight += interBlocksMargin;
    }
    file << "</svg>";

    delete [] tmpColumn;

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
        Debug.Report(ErrorCode::NotAligned, new std::string[1] { filename });
//         cerr << endl << "ERROR: Sequences are not aligned." << endl << endl;
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

void newAlignment::updateSequencesAndResiduesNums(bool countSequences, bool countResidues)
{
    int i;
    if (countSequences)
    {
        for (sequenNumber = 0, i = 0; i < originalSequenNumber; i++)
        {
            if (saveSequences[i] != -1) sequenNumber++;
        }
    }

    if (countResidues)
    {
        for (residNumber = 0, i = 0; i < originalResidNumber; i++)
        {
            if (saveResidues[i] != -1) residNumber++;
        }
    }
}
