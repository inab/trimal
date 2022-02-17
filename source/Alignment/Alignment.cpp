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
#include "Statistics/Similarity.h"
#include "Statistics/Consistency.h"
#include "InternalBenchmarker.h"
#include "Statistics/Manager.h"
#include "Alignment/sequencesMatrix.h"
#include "reportsystem.h"
#include "Cleaner.h"
#include "utils.h"
#include "Statistics/Similarity.h"
#include "Statistics/Consistency.h"
#include "InternalBenchmarker.h"
#include "Statistics/Manager.h"
#include "Alignment/sequencesMatrix.h"
#include "reportsystem.h"
#include "Cleaner.h"
#include "utils.h"
#include "omp.h"

Alignment::Alignment() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("Alignment::Alignment(void) ");

    // Submodules
    Cleaning = new Cleaner(this);
    Statistics = new statistics::Manager(this);

    // Sizes
    originalNumberOfResidues = 0;
    originalNumberOfSequences = 0;
    numberOfSequences = 0;
    numberOfResidues = 0;

    // Are the input sequences aligned?
    isAligned = false;

    // Sequence datatype: DNA, RNA or Protein
    dataType = SequenceTypes::NotDefined;

    // Columns and sequences that have been previously selected
    saveResidues = nullptr;
    saveSequences = nullptr;

    // Input sequences as well other information such as sequences name, etc 
    sequences = nullptr;
    seqsName = nullptr;
    seqsInfo = nullptr;

    // Information computed from Alignment
    SequencesMatrix = nullptr;
    identities = nullptr;

    // Pointer that helps to keep a count on how many items use the same information
    //  Upon destruction of the object, the counter is decreased, and, 
    //  when it arrives to zero, it removes the shared information from the memory.
    SeqRef = new int(1);

}

Alignment::Alignment(Alignment &originalAlignment) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("Alignment::Alignment(Alignment &originalAlignment) ");

    if (this != &originalAlignment) {
        filename = originalAlignment.filename;
        alignmentInfo = originalAlignment.alignmentInfo;
        dataType = originalAlignment.dataType;
        isAligned = originalAlignment.isAligned;

        seqsName = originalAlignment.seqsName;
        seqsInfo = originalAlignment.seqsInfo;
        sequences = originalAlignment.sequences;

        numberOfResidues = originalAlignment.numberOfResidues;
        numberOfSequences = originalAlignment.numberOfSequences;
        originalNumberOfResidues = originalAlignment.originalNumberOfResidues;
        originalNumberOfSequences = originalAlignment.originalNumberOfSequences;

        identities = nullptr;
        SequencesMatrix = nullptr;

        // Copy save(Sequences|Residues) vector to keep information of previous
        //  trimming steps
        saveSequences = new int[originalNumberOfSequences];
        if (originalAlignment.saveSequences != nullptr)
            std::copy(originalAlignment.saveSequences,
                      originalAlignment.saveSequences + originalAlignment.originalNumberOfSequences,
                      saveSequences);
        saveResidues = new int[originalNumberOfResidues];
        if (originalAlignment.saveResidues != nullptr)
            std::copy(originalAlignment.saveResidues,
                      originalAlignment.saveResidues + originalAlignment.originalNumberOfResidues,
                      saveResidues);

        // Submodules
        this->Cleaning = new Cleaner(this, originalAlignment.Cleaning);
        this->Statistics = new statistics::Manager(this, originalAlignment.Statistics);

        // Increase the number of elements using the shared information.
        this->SeqRef = originalAlignment.SeqRef;
        (*SeqRef)++;

    }
}

Alignment::~Alignment() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("Alignment::~Alignment) ");

    delete[] saveResidues;

    delete[] saveSequences;

    if (identities != nullptr) {
        for (int i = 0; i < numberOfSequences; i++)
            delete[] identities[i];
        delete[] identities;
    }

    delete SequencesMatrix;

    delete Cleaning;

    delete Statistics;

    if (--(*SeqRef) == 0) {
        delete SeqRef;

        delete[] sequences;

        delete[] seqsName;

        delete[] seqsInfo;

    } else if (*SeqRef < 0) {
        delete SeqRef;
    }
}

Alignment *Alignment::getTranslationCDS(Alignment *ProtAlig) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("Alignment *Alignment::getTranslationCDS(Alignment *ProtAlig) ");

    int x, y, counter = 0, *mappedSeqs;
    Alignment *newAlig = new Alignment();

    // Map the selected protein sequences to the input
    // coding sequences
    mappedSeqs = new int[ProtAlig->originalNumberOfSequences];

    for (x = 0; x < ProtAlig->originalNumberOfSequences; x++) {
        for (y = 0; y < originalNumberOfSequences; y++) {
            if (ProtAlig->seqsName[x] == seqsName[y]) {
                mappedSeqs[x] = y;
                break;
            }
        }
    }

    newAlig->sequences = new std::string[ProtAlig->originalNumberOfSequences];
    newAlig->seqsInfo = new std::string[ProtAlig->originalNumberOfSequences];
    newAlig->seqsName = new std::string[ProtAlig->originalNumberOfSequences];

    for (x = 0; x < ProtAlig->originalNumberOfSequences; x++) {

        newAlig->sequences[x] = std::string();
        std::string &current = newAlig->sequences[x];

        if (ProtAlig->seqsInfo != nullptr)
            newAlig->seqsInfo[x] = std::string(ProtAlig->seqsInfo[x]);
        newAlig->seqsName[x] = std::string(ProtAlig->seqsName[x]);

        for (y = 0, counter = 0; y < ProtAlig->sequences[x].size(); y++) {
            if (ProtAlig->saveResidues[y] == -1) {
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
            counter++;
        }
        current.shrink_to_fit();
    }

    newAlig->saveSequences = new int[ProtAlig->originalNumberOfSequences];
    std::copy(
            ProtAlig->saveSequences,
            ProtAlig->saveSequences + ProtAlig->originalNumberOfSequences,
            newAlig->saveSequences);

    newAlig->saveResidues = new int[ProtAlig->originalNumberOfResidues * 3];
    for (int i = 0; i < counter; i++) {
        newAlig->saveResidues[i] = i;
    }

    newAlig->numberOfSequences = ProtAlig->numberOfSequences;
    newAlig->originalNumberOfSequences = ProtAlig->originalNumberOfSequences;

    newAlig->numberOfResidues = counter;
    newAlig->originalNumberOfResidues = ProtAlig->originalNumberOfResidues * 3;

    delete[] mappedSeqs;

    return newAlig;

}

int Alignment::getNumSpecies(void) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("int Alignment::getNumSpecies(void) ");
    return numberOfSequences;
}

int Alignment::getNumAminos(void) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("int Alignment::getNumAminos(void) ");
    return numberOfResidues;
}

void Alignment::setWindowsSize(int ghWindow_, int shWindow_) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::setWindowsSize(int ghWindow_, int shWindow_) ");
    Statistics->ghWindow = ghWindow_;
    Statistics->shWindow = shWindow_;
}

void Alignment::setBlockSize(int blockSize_) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::setBlockSize(int blockSize_) ");
    Cleaning->blockSize = blockSize_;
}

void Alignment::setKeepSequencesFlag(bool flag) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::setKeepSequencesFlag(bool flag) ");
    Cleaning->keepSequences = flag;
}
/*
void Alignment::setKeepSeqsHeaderFlag(bool newFlagValue) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::setKeepSeqsHeaderFlag(bool flag) ");
    Cleaning->keepSequences = newFlagValue;
}
*/
/*
int Alignment::getBlockSize() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("int Alignment::getBlockSize() ");
    return Cleaning->blockSize;
}
 */
/*
void Alignment::calculateSeqIdentity() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::calculateSeqIdentity(void) ");

    int i, j, k, hit, dst;
    char indet;

    // Depending on alignment type, indetermination symbol will be one or other
    indet = getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

    // Create identities matrix to store identities scores
    identities = new float *[numberOfSequences];

    // For each seq, compute its identity score against the others in the MSA
    for (i = 0; i < numberOfSequences; i++) {
        identities[i] = new float[numberOfSequences];

        // It's a symmetric matrix, copy values that have been already computed
        for (j = 0; j < i; j++)
            identities[i][j] = identities[j][i];
        identities[i][i] = 0;

        // Compute identity scores for the current sequence against the rest
        for (j = i + 1; j < numberOfSequences; j++) {
            for (k = 0, hit = 0, dst = 0; k < numberOfResidues; k++) {
                // If one of the two positions is a valid residue,
                // count it for the common length
                if (((sequences[i][k] != indet) && (sequences[i][k] != '-')) ||
                    ((sequences[j][k] != indet) && (sequences[j][k] != '-'))) {
                    dst++;
                    // If both positions are the same, count a hit
                    if (sequences[i][k] == sequences[j][k])
                        hit++;
                }
            }

            // Identity score between two sequences is the ratio of identical residues
            // by the total length (common and no-common residues) among them
            identities[i][j] = (float) hit / dst;
        }
    }
}
 */
/*
void Alignment::calculateRelaxedSeqIdentity() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::calculateRelaxedSeqIdentity(void) ");
    // Raw approximation of sequence identity computation designed for reducing
    // comparisons for huge alignemnts

    int i, j, k, hit;

    // Create identities matrix to store identities scores
    identities = new float *[numberOfSequences];

    // For each seq, compute its identity score against the others in the MSA
    for (i = 0; i < numberOfSequences; i++) {
        identities[i] = new float[numberOfSequences];

        // It's a symmetric matrix, copy values that have been already computed
        for (j = 0; j < i; j++)
            identities[i][j] = identities[j][i];
        identities[i][i] = 0;

        // Compute identity score between the selected sequence and the others
        for (j = i + 1; j < numberOfSequences; j++) {
            for (k = 0, hit = 0; k < numberOfResidues; k++) {
                // If both positions are the same, count a hit
                if (sequences[i][k] == sequences[j][k])
                    hit++;
            }
            // Raw identity score is computed as the ratio of identical residues between
            // alignment length
            identities[i][j] = (float) hit / numberOfResidues;
        }
    }
}
*/

void Alignment::calculateSeqOverlap() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::calculateSeqOverlap(void) ");
    // Compute the overlap between sequences taken each of them as the reference
    // to compute such scores. It will lead to a non-symmetric matrix.

    int i, j, k, shared, referenceLength;
    char indet;

    // Depending on alignment type, indetermination symbol will be one or other
    indet = getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

    // Create overlap matrix to store overlap scores
    overlaps = new float *[numberOfSequences];

    #pragma omp parallel for num_threads(NUMTHREADS) if(numberOfSequences>MINPARALLELSIZE)
    for(i = 0; i < numberOfSequences; i++) overlaps[i] = new float[numberOfSequences];

    #pragma omp parallel for num_threads(NUMTHREADS) collapse(2) private(k, shared, referenceLength) if(numberOfSequences>MINPARALLELSIZE)
    // For each seq, compute its overlap score against the others in the MSA
    for (i = 0; i < numberOfSequences; i++) {
        for (j = 0; j < numberOfSequences; j++) {
            for (k = 0, shared = 0, referenceLength = 0; k < numberOfResidues; k++) {
                // If there a valid residue for the reference sequence, then see if
                // there is a valid residue for the other sequence.
                if ((sequences[i][k] != indet) && (sequences[i][k] != '-')) {
                    referenceLength++;
                    if ((sequences[j][k] != indet) && (sequences[j][k] != '-'))
                        shared++;
                }
            }
            // Overlap score between two sequences is the ratio of shared valid
            // residues divided by the sequence length taken as reference. The
            // overlaps matrix, therefore, will be not symmetric.
            overlaps[i][j] = (float) shared / referenceLength;
        }
    }
}

void Alignment::getSequences(string *Names) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::getSequences(string *Names) ");
    for (int i = 0; i < numberOfSequences; i++)
        Names[i] = seqsName[i];
}

void Alignment::getSequences(string *Names, int *lenghts) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::getSequences(string *Names, int *lengths) ");
    for (int i = 0; i < numberOfSequences; i++) {
        lenghts[i] = (int) utils::removeCharacter('-', sequences[i]).length();
        Names[i] = seqsName[i];
    }
}

void Alignment::getSequences(string *Names, string *Sequences, int *lengths) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::getSequences(string *Names, string *Sequences, int *lengths) ");
    for (int i = 0; i < numberOfSequences; i++) {
        Names[i] = seqsName[i];
        Sequences[i] = utils::removeCharacter('-', sequences[i]);
        lengths[i] = (int) Sequences[i].length();
    }
}

bool Alignment::getSequenceNameOrder(string *names, int *orderVector) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool Alignment::getSequenceNameOrder(string *names, int *orderVector) ");
    int i, j, numNames;

    // For each name in the current Alignment, we look
    // for its correspondence in the input set
    for (i = 0, numNames = 0; i < numberOfSequences; i++) {
        for (j = 0; j < numberOfSequences; j++) {
            if (seqsName[i] == names[j]) {
                orderVector[i] = j;
                numNames++;
                break;
            }
        }
    }

    // Depending on if we get the complete correspondence
    // between both sets of names, we return true or not
    return numberOfSequences == numNames ? true : false;

}

int Alignment::getAlignmentType() const {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("int Alignment::getAlignmentType(void) ");
    if (dataType == SequenceTypes::NotDefined)
        // Dropping constness,
        //      as the datatype is not a "modification" of the alignment
        const_cast<Alignment*>(this)->dataType = utils::checkAlignmentType(numberOfSequences, numberOfResidues, sequences);
    return dataType;
}
/*
int *Alignment::getCorrespResidues() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("int *Alignment::getCorrespResidues(void) ");
    return saveResidues;
}
*/
/*
int *Alignment::getCorrespSequences() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("int *Alignment::getCorrespSequences(void) ");
    return saveSequences;
}
*/
bool Alignment::isFileAligned() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool Alignment::isFileAligned(void) ");
    return isAligned;
}

/*
void Alignment::fillNewDataStructure(string *newMatrix, string *newNames) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::fillNewDataStructure(string *newMatrix, string *newNames) ");
    // TODO Remove method as it doesn't make sense now
    int i, j, k;

    // Copy only those sequences/columns selected
    for (i = 0, j = 0; i < numberOfSequences; i++) {
        if (saveSequences[i] == -1)
            continue;

        newNames[j] = seqsName[i];
        for (k = 0; k < numberOfResidues; k++) {
            if (saveResidues[k] == -1)
                continue;
            newMatrix[j].resize(newMatrix[j].size() + 1, sequences[i][k]);
        }
        j++;
    }
}
*/



bool Alignment::prepareCodingSequence(bool splitByStopCodon, bool ignStopCodon, Alignment *proteinAlig) {

#define deletePrepareCodingVariables() { \
        delete [] protSeqsNames; protSeqsNames = nullptr; \
        delete [] protSequences; protSequences = nullptr; \
        delete [] protSeqsLengths; protSeqsLengths = nullptr; \
    }

    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool Alignment::prepareCodingSequence(bool splitByStopCodon, bool ignStopCodon, Alignment *proteinAlig) ");

    if (getAlignmentType() & SequenceTypes::AA) {
        debug.report(ErrorCode::CDScontainsProteinSequences);
        return false;
    }

    bool warning = false;
    long found;
    int i;

    // New code: We now care about the presence of wildcards/indeterminations
    // characters such as 'X' or 'B' into protein sequences as well as about the
    // presence of Selenocysteines ('U') or Pyrrolysines ('O'). It only works, by
    // now, with the universal genetic code 
    int *protSeqsLengths, numbProtSeqs, current_prot;
    string *protSeqsNames, *protSequences;
    char aminoAcid;

    numbProtSeqs = proteinAlig->getNumSpecies();
    protSeqsNames = new string[numbProtSeqs];
    protSequences = new string[numbProtSeqs];
    protSeqsLengths = new int[numbProtSeqs];

    proteinAlig->getSequences(protSeqsNames, protSequences, protSeqsLengths);

    /* Check read sequences are real DNA/RNA */

    for (i = 0; i < numberOfSequences; i++) {

        // Get protein sequence to compare against any potential stop codon in the
        // coding sequence. If there is not protein sequence for current coding
        // sequence, skip its analysis
        for (current_prot = 0; current_prot < numbProtSeqs; current_prot++)
            if (protSeqsNames[current_prot] == seqsName[i])
                break;
        if (current_prot == numbProtSeqs)
            continue;

        if (sequences[i].find('-') != string::npos) {
            if (!warning)
                cerr << endl;
            debug.report(ErrorCode::SequenceContainsGap, new std::string[1]{seqsName[i]});
            deletePrepareCodingVariables()
            return false;
        }

        if ((sequences[i].length() % 3) != 0) {
            if (!warning)
                cerr << endl;
            warning = true;
            debug.report(ErrorCode::SequenceNotMultipleOfThree, new std::string[2]{seqsName[i], std::to_string(sequences[i].length())});
        }

        // Ignore stop codons from the CDS if set by the user
        if (ignStopCodon)
            continue;

        // Detect universal stop codons in the CDS. Then, compare those stop codons
        // against the protein sequence to see whether they are real stop codons or
        // are representing rare amino-acids such as Selenocysteines or Pyrrolysines
        // It also allows stop-codons when there are wildcards/indet characters in
        // the protein sequence. CDS sequences could be splitted using stop codons
        // from the sequence itself

        // Initialize first appearence of a given stop codon to -1.
        // That means that it has not been found yet
        found = -1;
        do {
            found = sequences[i].find("TGA", found + 1);

            // If a stop codon has been found and its position is multiple of 3.
            // Analize it
            if ((found != string::npos) && (((int) found % 3) == 0)) {

                aminoAcid = (char) toupper(protSequences[current_prot][(int) found / 3]);
                // It may be a Selenocysteine ('TGA') which should be represented as 'U'
                // or wildcard/indet characters such as 'X' or 'B'

                // If a rare amino-acids such as 'U'/'O' or a wildcard/indet character
                // such as 'B'/'X' is present, skip current stop codon
                if ((aminoAcid == 'U') || (aminoAcid == 'O') ||
                    (aminoAcid == 'X') || (aminoAcid == 'B'))
                    continue;

                // If split_by_stop_codon flag is activated then cut input CDS sequence
                // up to first appearance of a stop codon
                else if (splitByStopCodon) {
                    warning = true;
                    debug.report(InfoCode::CuttingSequence, new std::string[5]{seqsName[i], "TGA", std::to_string(aminoAcid), std::to_string(found + 1), std::to_string(sequences[i].length())});
                    sequences[i].resize(found);
                }
                // Otherwise, warn about it and return an error
                else {
                    debug.report(ErrorCode::SequenceHasStopCodon, new std::string[5]{seqsName[i], "TGA", std::to_string(aminoAcid), std::to_string(found + 1), std::to_string(sequences[i].length())});
                    deletePrepareCodingVariables()
                    return false;
                }
            }
            // Iterate over the CDS until not stop codon is found
        } while (found != string::npos);

        // Initialize first appearence of a given stop codon to -1.
        // That means that it has not been found yet
        found = -1;
        do {
            found = sequences[i].find("TAA", found + 1);

            // If a stop codon has been found and its position is multiple of 3.
            // Analize it
            if ((found != string::npos) && (((int) found % 3) == 0)) {

                aminoAcid = (char) toupper(protSequences[current_prot][(int) found / 3]);
                // Check if there is any wildcard/indet characters such as 'X' or 'B'
                // If a rare amino-acids such as 'U'/'O' or a wildcard/indet character
                // such as 'B'/'X' is present, skip current stop codon
                if ((aminoAcid == 'U') || (aminoAcid == 'O') ||
                    (aminoAcid == 'X') || (aminoAcid == 'B'))
                    continue;

                    // If split_by_stop_codon flag is activated then cut input CDS sequence
                    // up to first appearance of a stop codon
                else if (splitByStopCodon) {
                    warning = true;
                    debug.report(InfoCode::CuttingSequence, new std::string[5]{seqsName[i], "TAA", std::to_string(aminoAcid), std::to_string(found + 1), std::to_string(sequences[i].length())});
                    sequences[i].resize(found);
                }
                    // Otherwise, warn about it and return an error
                else {
                    debug.report(ErrorCode::SequenceHasStopCodon, new std::string[5]{seqsName[i], "TAA", std::to_string(aminoAcid), std::to_string(found + 1), std::to_string(sequences[i].length())});
                    deletePrepareCodingVariables()
                    return false;
                }
            }
            // Iterate over the CDS until not stop codon is found
        } while (found != string::npos);

        // Initialize first appearence of a given stop codon to -1.
        // That means that it has not been found yet
        found = -1;
        do {
            found = sequences[i].find("TAG", found + 1);
            // If a stop codon has been found and its position is multiple of 3.
            // Analize it
            if ((found != string::npos) && (((int) found % 3) == 0)) {

                aminoAcid = (char) toupper(protSequences[current_prot][(int) found / 3]);
                // It may be a Pyrrolysine ('TAG') which should be represented as 'O'
                // or wildcard/indet characters such as 'X' or 'B'
                // If a rare amino-acids such as 'U'/'O' or a wildcard/indet character
                // such as 'B'/'X' is present, skip current stop codon
                if ((aminoAcid == 'U') || (aminoAcid == 'O') || (aminoAcid == 'X') || \
                        (aminoAcid == 'B'))
                    continue;

                    // If split_by_stop_codon flag is activated then cut input CDS sequence
                    // up to first appearance of a stop codon
                else if (splitByStopCodon) {
                    warning = true;
                    debug.report(InfoCode::CuttingSequence, new std::string[5]{seqsName[i], "TAG", std::to_string(aminoAcid), std::to_string(found + 1), std::to_string(sequences[i].length())});
                    sequences[i].resize(found);
                }
                    // Otherwise, warn about it and return an error
                else {
                    debug.report(ErrorCode::SequenceHasStopCodon, new std::string[5]{seqsName[i], "TAG", std::to_string(aminoAcid), std::to_string(found + 1), std::to_string(sequences[i].length())});
                    deletePrepareCodingVariables()
                    return false;
                }
            }
            // Iterate over the CDS until not stop codon is found
        } while (found != string::npos);
    }

    deletePrepareCodingVariables()
    // If everything was return an OK to informat about it.
    return true;
}

bool Alignment::checkCorrespondence(string *names, int *lengths, int totalInputSeqs, int multiple = 1) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool Alignment::checkCorrespondence(string *names, int *lengths, int totalInputSeqs, int multiple = 1)");

    int i, j, seqLength, indet;
    bool warnings = false;
    string tmp;

    // For each sequence in the current protein Alignment, look for its coding
    // DNA sequence checking that they have the same size.
    for (i = 0; i < numberOfSequences; i++) {

        // Get protein sequence length removing any posible gap. Get as well last
        // residue from current sequence 

        tmp = utils::removeCharacter('-', sequences[i]);
        seqLength = tmp.length() * multiple;
        indet = ((int) tmp.length() - utils::min((int) tmp.find_last_not_of('X'), \
                 (int) tmp.find_last_not_of('x'))) - 1;

        // Go through all available CDS looking for the one with the same ID
        for (j = 0; j < totalInputSeqs; j++) {

            // Once both ID matchs, compare its lengths
            if (seqsName[i] == names[j]) {

                // If both sequences have the same length, stop the search 
                if (seqLength == lengths[j])
                    break;

                    // If nucleotide sequence is larger than protein sequence, warn about
                    // it and continue the verification procefile. It will used the 'Nth'
                    // first nucleotides for the conversion
                else if (seqLength < lengths[j]) {
                    if (!warnings)
                        cerr << endl;
                    warnings = true;
                    debug.report(WarningCode::SequenceWillBeCut, new std::string[3]{seqsName[i], std::to_string(seqLength), std::to_string(lengths[j])});
                    break;
                }

                    // It has been detected some indeterminations at the end of the protein
                    // sequence. That ifileue could be cause by some incomplete codons in the
                    // nucleotide sequences. This issue is solved adding as much 'N' symbols
                    // as it is needed to preserve the backtranslated Alignment
                else if ((indet > 0) && (indet > (seqLength - lengths[j]) / 3)) {
                    if (!warnings)
                        cerr << endl;
                    warnings = true;
                    debug.report(WarningCode::IncludingIndeterminationSymbols, new std::string[1]{seqsName[i]});
                    break;
                }

                    // If nucleotide sequence is shorter than protein sequence, return an
                    // error since it is not feasible to cut the input protein aligment to
                    // fit it into CDNA sequences size
                else {
                    if (!warnings)
                        cerr << endl;
                    warnings = true;
                    debug.report(WarningCode::LessNucleotidesThanExpected, new std::string[3]{seqsName[i], std::to_string(lengths[j]), std::to_string(seqLength)});
                    break;
                }
            }
        }

        // Warn about a mismatch a sequences name level
        if (j == totalInputSeqs) {
            debug.report(ErrorCode::SequenceNotPresentInCDS, new std::string[1]{seqsName[i]});
            return false;
        }
    }

    // If everything is OK, return an appropiate flag
    return true;
}

bool Alignment::fillMatrices(bool aligned, bool checkInvalidChars) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool Alignment::fillMatrices(bool aligned) ");
    // Function to determine if a set of sequences, that can be aligned or not,
    // have been correctly load and are free of errors.
    int i, j;

    // Initialize some variables

    // Check whether there are any unknow/no allowed character in the sequences
    if (checkInvalidChars)
        for (i = 0; i < numberOfSequences; i++)
            for (j = 0; j < sequences[i].length(); j++)
                if ((!isalpha(sequences[i][j])) && (!ispunct(sequences[i][j]))) {
                    debug.report(ErrorCode::UnknownCharacter, new std::string[2]{seqsName[i], std::to_string(sequences[i][j])});
                    return false;
                }

    // Check whether all sequences have same size or not
    for (i = 1; i < numberOfSequences; i++) {
        if (sequences[i].length() != sequences[i - 1].length())
            break;
    }
    // Set an appropriate flag for indicating if sequences are aligned or not
    isAligned = i == numberOfSequences;

    // Warm about those cases where sequences should be aligned
    // and there are not
    if (aligned and !isAligned) {
        debug.report(ErrorCode::NotAligned, new std::string[1]{filename});
        return false;
    }

    // Full-fill some information about input alignment
    if (numberOfResidues == 0)
        numberOfResidues = sequences[0].length();

    // Check whether aligned sequences have the length fixed for the input alig
    for (i = 0; (i < numberOfSequences) and (aligned); i++) {
        if (sequences[i].length() != numberOfResidues) {
            debug.report(ErrorCode::SequencesNotSameSize,
                    new std::string[3]{
                seqsName[i],
                std::to_string(sequences[i].length()),
                std::to_string(numberOfResidues)
            });
            return false;
        }
    }
    // If the sequences are aligned, initialize some additional variables.
    // These variables will be useful for posterior analysis
    if ((aligned) || (isAligned)) {

        // Assign its position to each column. That will be used to determine which
        // columns should be kept in output alignment after applying any method
        // and which columns should not
        saveResidues = new int[numberOfResidues];
        for (i = 0; i < numberOfResidues; i++)
            saveResidues[i] = i;

        // Assign its position to each sequence. Similar to the columns numbering
        // possible, assign to each sequence its position is useful to know which
        // sequences will be in the output alignment
        saveSequences = new int[numberOfSequences];
        for (i = 0; i < numberOfSequences; i++)
            saveSequences[i] = i;
    }

    // Return an flag indicating that everything is fine

    return true;
}

void Alignment::printAlignmentInfo(ostream &file) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::printAlignmentInfo(ostream &file) ");
    // Print information about sequences number, average sequence length, maximum
    // and minimum sequences length, etc

    int i, j, valid_res, max, min, max_pos, min_pos, total_res;

    // Storage which sequences are the longest and shortest ones
    max = 0;
    max_pos = 0;
    min_pos = 0;
    min = sequences[0].length();

    for (i = 0, total_res = 0; i < numberOfSequences; i++) {

        // Discard gaps from current sequence and then compute real length
        for (j = 0, valid_res = 0; j < sequences[i].length(); j++)
            valid_res += (sequences[i][j] != '-' ? 1 : 0);

        // Compute the total residues in the alignment to calculate avg. sequence
        // length
        total_res += valid_res;

        // Get values for the longest sequence
        max_pos = (max > valid_res) ? max_pos : i;
        max = (max > valid_res) ? max : valid_res;

        // Similarly, get values for the shortest sequence
        min_pos = (min < valid_res) ? min_pos : i;
        min = (min < valid_res) ? min : valid_res;
    }

    file << "## Total sequences\t" << numberOfSequences << endl;
    if (isFileAligned())
        file << "## Alignment length\t" << numberOfResidues << endl;
    file << "## Avg. sequence length\t" << (float) total_res / numberOfSequences << endl
         << "## Longest seq. name\t'" << seqsName[max_pos] << "'" << endl
         << "## Longest seq. length\t" << max << endl
         << "## Shortest seq. name\t'" << seqsName[min_pos] << "'" << endl
         << "## Shortest seq. length\t" << min << endl;
}

void Alignment::printSeqIdentity() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::printSeqIdentity(void) ");

    int i, j, k, pos, maxLongName;
    float mx, avg, maxAvgSeq = 0, maxSeq = 0, avgSeq = 0, **maxs;

    // Ask for the sequence identities assesment
    if (identities == nullptr)
        Cleaning->calculateSeqIdentity();

    // For each sequence, we look for its most similar one
    maxs = new float *[originalNumberOfSequences];

    #pragma omp parallel for private(k, mx, avg, pos) reduction(+: avgSeq, maxAvgSeq) num_threads(NUMTHREADS) if(originalNumberOfSequences>MINPARALLELSIZE)
    for (i = 0; i < originalNumberOfSequences; i++) {
        maxs[i] = new float[2];

        // Get the most similar sequence to the current one in term of identity
        for (k = 0, mx = 0, avg = 0, pos = i; k < originalNumberOfSequences; k++) {
            if (i != k) {
                avg += identities[i][k];
                if (mx < identities[i][k]) {
                    mx = identities[i][k];
                    pos = k;
                }
            }
        }
        // Update global average variables
        avgSeq += avg / (originalNumberOfSequences - 1);
        maxAvgSeq += mx;

        // Save the maximum average identity value for each sequence
        maxs[i][0] = mx;
        maxs[i][1] = pos;
    }

    // Compute general averages 
    avgSeq = avgSeq / originalNumberOfSequences;
    maxAvgSeq = maxAvgSeq / originalNumberOfSequences;

    // Compute longest sequences name
    for (i = 0, maxLongName = 0; i < originalNumberOfSequences; i++)
        maxLongName = utils::max(maxLongName, seqsName[i].size());

    // Once the method has computed all of different values, it prints it
    cout.precision(4);
    cout << fixed;

    for (i = 0, maxSeq = 0; i < originalNumberOfSequences; i++)
        if (maxs[i][0] > maxSeq)
            maxSeq = maxs[i][0];

    cout << endl << "## MaxIdentity\t" << maxSeq;
    cout << endl << "#> MaxIdentity\tGet the maximum identity value for any pair "
         << "of sequences in the alignment" << endl;

    cout << endl << "## AverageIdentity\t" << avgSeq;
    cout << endl << "#> AverageIdentity\tAverage identity between all sequences";

    cout << endl << endl << "## Identity sequences matrix";
    for (i = 0; i < numberOfSequences; i++) {
        cout << endl << setw(maxLongName + 2) << left << seqsName[i] << "\t";
        for (j = 0; j < i; j++)
            cout << setiosflags(ios::left) << setw(10) << identities[i][j] << "\t";
        cout << setiosflags(ios::left) << setw(10) << 1.00 << "\t";
        for (j = i + 1; j < numberOfSequences; j++)
            cout << setiosflags(ios::left) << setw(10) << identities[i][j] << "\t";
    }
    cout << endl;

    cout << endl << "## AverageMostSimilarIdentity\t" << maxAvgSeq;
    cout << endl << "#> AverageMostSimilarIdentity\t Average identity between "
         << "most similar pair-wise sequences";

    cout << endl << endl << "## Identity for most similar pair-wise sequences "
         << "matrix" << endl;
    for (i = 0; i < numberOfSequences; i++)
        cout << setw(maxLongName + 2) << left << seqsName[i]
             << "\t" << setiosflags(ios::left) << setw(5)
             << maxs[i][0] << "\t" << seqsName[(int) maxs[i][1]] << endl;
    cout << endl;

    for (i = 0; i < numberOfSequences; i++) {
        delete maxs[i];
    }
    delete[] maxs;
}

void Alignment::printSeqOverlap() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::printSeqOverlap() ");
    int i, j, k, pos, maxLongName;
    float mx, avg, maxAvgSeq = 0, maxSeq = 0, avgSeq = 0, **maxs;

    // Ask for the sequence identities assessment
    if (overlaps == nullptr)
        calculateSeqOverlap();

    // For each sequence, we look for its most similar one 
    maxs = new float *[numberOfSequences];

    for (i = 0; i < numberOfSequences; i++) {
        maxs[i] = new float[2];

        // Get the most similar sequence to the current one in term of overlap
        for (k = 0, mx = 0, avg = 0, pos = i; k < numberOfSequences; k++) {
            if (i != k) {
                avg += overlaps[i][k];
                if (mx < overlaps[i][k]) {
                    mx = overlaps[i][k];
                    pos = k;
                }
            }
        }
        // Update global average variables
        avgSeq += avg / (numberOfSequences - 1);
        maxAvgSeq += mx;

        // Save the maximum average overlap value for each sequence
        maxs[i][0] = mx;
        maxs[i][1] = pos;
    }

    // Compute general averages
    avgSeq = avgSeq / numberOfSequences;
//    maxAvgSeq = maxAvgSeq / numberOfSequences;

    // Compute longest sequences name 
    for (i = 0, maxLongName = 0; i < numberOfSequences; i++)
        maxLongName = utils::max(maxLongName, seqsName[i].size());

    // Once the method has computed all of different values, it prints it
    cout.precision(4);
    cout << fixed;

    for (i = 0, maxSeq = 0; i < numberOfSequences; i++)
        if (maxs[i][0] > maxSeq)
            maxSeq = maxs[i][0];

    cout << "## MaxOverlap\t" << maxSeq;
    cout << endl << "#> MaxOverlap\tGet the maximum overlap value for any pair "
         << "of sequences in the alignment" << endl;

    cout << endl << "## AverageOverlap\t" << avgSeq;
    cout << endl << "#> AverageOverlap\tAverage overlap between all sequences";

    cout << endl << endl << "## Overlap sequences matrix";
    for (i = 0; i < numberOfSequences; i++) {
        cout << endl << setw(maxLongName + 2) << left << seqsName[i] << "\t";
        for (j = 0; j < numberOfSequences; j++)
            cout << setiosflags(ios::left) << setw(10) << overlaps[i][j] << "\t";
    }
    cout << endl;

    for (i = 0; i < numberOfSequences; i++) {
        delete maxs[i];
    }
    delete[] maxs;
}

void Alignment::calculateColIdentity(float *ColumnIdentities) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::calculateColIdentity(float *ColumnIdentities) ");

    int i, j, counter, pos, max, columnLen;
    char letter, indet, gapSymbol;
    string column;

    // Initialize some data for make computation more precise
    indet = getAlignmentType() == SequenceTypes::AA ? 'X' : 'N';
    gapSymbol = '-';

    // Compute identity score for the most frequent residue, it can be as well
    // gaps and indeterminations, for each column 
    for (i = 0, max = 0; i < numberOfResidues; i++, max = 0, column.clear()) {

        // Get residues from each column in capital letters
        for (j = 0; j < numberOfSequences; j++)
            // Discard gaps and indeterminations from calculations
            if ((toupper(sequences[j][i]) != indet) && (sequences[j][i] != gapSymbol))
                column += toupper(sequences[j][i]);
        columnLen = column.size();

        // Count letter frequency. It only matter the frequency. Use some shorcuts
        // to speed-up the process
        while (!column.empty()) {
            letter = column[0];
            counter = 0;
            pos = 0;

            do {
                counter += 1;
                column.erase(pos, 1);
                pos = column.find(letter, pos);
            } while (pos != (int) string::npos);

            // Keep only the most frequent residue
            if (counter > max)
                max = counter;
            // If column size is smaller than the current max, stop the count
            if ((int) column.size() < max)
                break;
        }

        // Store column identity values
        if (columnLen != 0)
            ColumnIdentities[i] = float(max) / columnLen;
    }
}

/*
void Alignment::printColumnsIdentity_DescriptiveStats(void) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::printColumnsIdentity_DescriptiveStats(void) ");

    float *colIdentities, avg, std, max, min;
    int i, positions;

    // Allocate local memory for the computation
    colIdentities = new float[numberOfResidues];

    utils::initlVect(colIdentities, numberOfResidues, -1);
    calculateColIdentity(colIdentities);

    for (i = 0, max = 0, min = 1, avg = 0, positions = 0; i < numberOfResidues; i++) {
        if (colIdentities[i] != -1) {
            // Compute on-the-fly max and min scores. Store accumulative score
            avg += colIdentities[i];
            max = (colIdentities[i] > max) ? colIdentities[i] : max;
            min = (colIdentities[i] < min) ? colIdentities[i] : min;
            // Count how many columns have a value score
            positions += 1;
        }
    }
    // Compute average identity column score
    avg /= positions;

    // Compute standard deviation
    for (i = 0, std = 0; i < numberOfResidues; i++)
        if (colIdentities[i] != -1)
            std += pow((colIdentities[i] - avg), 2);
    std = sqrt(std / positions);

    // Print general descriptive stats
    cout << "#maxColIdentity\t" << max << endl;
    cout << "#minColIdentity\t" << min << endl;
    cout << "#avgColIdentity\t" << avg << endl;
    cout << "#stdColIdentity\t" << std << endl;
}
 */

bool Alignment::alignmentSummaryHTML(const Alignment &trimmedAlig, const char *const destFile) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool Alignment::alignmentSummaryHTML(Alignment &trimmedAlig, char *destFile, float *consValues) ");

    // Generate an HTML file with a visual summary about which sequences/columns
    // have been selected and which have not

    int i, j, k, kj, upper, minHTML, maxLongName, *gapsValues;
    string tmpColumn;
    ofstream file;
    char type;

    // Allocate some local memory
    tmpColumn.reserve(originalNumberOfSequences);

    // Check whether sequences in the alignment are aligned or not.
    // Warn about it if there are not aligned.
    if (!isAligned) {
        debug.report(ErrorCode::NotAligned, new std::string[1]{filename});
        return false;
    }

    // Open output file and check that file pointer is valid 
    file.open(destFile);
    if (!file)
        return false;

    // Compute maximum sequences name length.
    maxLongName = 0;
    for (i = 0; i < originalNumberOfSequences; i++)
        maxLongName = utils::max(maxLongName, seqsName[i].size());

    // Compute HTML blank spaces
    minHTML = utils::max(25, maxLongName + 10);

    // Recover some stats about different scores from current alignment
    gapsValues = nullptr;
    if (trimmedAlig.Statistics->gaps != nullptr)
        gapsValues = trimmedAlig.Statistics->gaps->getGapsWindow();

    float *simValues = nullptr;
    if (trimmedAlig.Statistics->similarity != nullptr)
        simValues = trimmedAlig.Statistics->similarity->getMdkWindowedVector();

    float *consValues = nullptr;
    if (trimmedAlig.Statistics->consistency != nullptr)
        consValues = trimmedAlig.Statistics->consistency->getValues();

    // Print HTML header into output file
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

         // Sets of colors for high-lighting scores intervals
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

         // Other HTML elements
         << "    </style>\n  </head>\n\n" << "  <body>\n" << "  <pre>" << endl;

    // Show information about how many sequences/residues have been selected
    file << "    <span class=sel>Selected Sequences: " << setw(5) << right << trimmedAlig.numberOfSequences
         << " /Selected Residues: " << setw(7) << right << trimmedAlig.numberOfResidues << "</span>"
         << endl << "    <span class=nsel>Deleted Sequences:  " << setw(5) << right
         << numberOfSequences - trimmedAlig.numberOfSequences << " /Deleted Residues:  " << setw(7) << right
         << numberOfResidues - trimmedAlig.numberOfResidues << "</span>" << endl;

    // Print headers for different scores derived from input alignment/s 
    if (gapsValues != nullptr)
        file << endl << setw(minHTML) << left << "    Gaps Scores:        "
             << "<span  class=c1>  =0=  </span><span  class=c2> <.001 </span>"
             << "<span  class=c3> <.050 </span><span  class=c4> <.100 </span>"
             << "<span  class=c5> <.150 </span><span  class=c6> <.200 </span>"
             << "<span  class=c7> <.250 </span><span  class=c8> <.350 </span>"
             << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
             << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";

    if (simValues != nullptr)
        file << endl << setw(minHTML) << left << "    Similarity Scores:  "
             << "<span  class=c1>  =0=  </span><span  class=c2> <1e-6 </span>"
             << "<span  class=c3> <1e-5 </span><span  class=c4> <1e-4 </span>"
             << "<span  class=c5> <.001 </span><span  class=c6> <.010 </span>"
             << "<span  class=c7> <.100 </span><span  class=c8> <.250 </span>"
             << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
             << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";

    if (consValues != nullptr)
        file << endl << setw(minHTML) << left << "    Consistency Scores: "
             << "<span  class=c1>  =0=  </span><span  class=c2> <.001 </span>"
             << "<span  class=c3> <.050 </span><span  class=c4> <.100 </span>"
             << "<span  class=c5> <.150 </span><span  class=c6> <.200 </span>"
             << "<span  class=c7> <.250 </span><span  class=c8> <.350 </span>"
             << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
             << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";

    if ((gapsValues != nullptr) or (simValues == nullptr) or (consValues == nullptr))
        file << endl;

    // Print Sequences in block of BLOCK_SIZE 
    for (j = 0, upper = HTMLBLOCKS; j < originalNumberOfResidues; j += HTMLBLOCKS, upper += \
            HTMLBLOCKS) {

        // Print main columns number
        file << endl << setw(minHTML + 10) << right << (j + 10);
        for (i = j + 20; ((i <= originalNumberOfResidues) && (i <= upper)); i += 10)
            file << setw(10) << right << (i);

        // Print special characters to delimit sequences blocks
        file << endl << setw(minHTML + 1) << right;
        for (i = j + 1; ((i <= originalNumberOfResidues) && (i <= upper)); i++)
            file << (!(i % 10) ? "+" : "=");
        file << endl;

        // Print sequences name
        for (i = 0; i < originalNumberOfSequences; i++) {
            file << "    <span class=" << ((trimmedAlig.saveSequences[i] != -1) ? "sel>" : "nsel>") << seqsName[i]
                 << "</span>" << setw(minHTML - 4 - seqsName[i].size()) << right << "";

            // Print residues corresponding to current sequences block
            for (k = j; ((k < originalNumberOfResidues) && (k < upper)); k++) {
                for (kj = 0, tmpColumn.clear(); kj < originalNumberOfSequences; kj++)
                    tmpColumn += sequences[kj][k];
                // Determine residue color based on residues across the alig column
                type = utils::determineColor(sequences[i][k], tmpColumn);

                if (type == 'w')
                    file << sequences[i][k];
                else
                    file << "<span id=" << type << ">" << sequences[i][k] << "</span>";
            }

            file << endl;
        }

        file << endl << setw(minHTML) << left << "    Selected Cols:      ";
        for (k = j; ((k < originalNumberOfResidues) && (k < (j + HTMLBLOCKS))); k++)
            file << "<span class=" << (trimmedAlig.saveResidues[k] != -1 ? "sel" : "nsel") << "> </span>";
        file << endl;

        // If there is not any score to print, skip this part of the function
        if ((gapsValues == nullptr) and (simValues == nullptr) and (consValues == nullptr))
            continue;

        // Print score colors according to certain predefined thresholds
        if (gapsValues != nullptr) {
            file << endl << setw(minHTML) << left << "    Gaps Scores:        ";
            for (k = j; ((k < originalNumberOfResidues) && (k < (j + HTMLBLOCKS))); k++)
                if (gapsValues[k] == 0)
                    file << "<span class=c12> </span>";
                else if (gapsValues[k] == originalNumberOfSequences)
                    file << "<span class=c1> </span>";
                else if (1 - (float(gapsValues[k]) / originalNumberOfSequences) >= .750)
                    file << "<span class=c11> </span>";
                else if (1 - (float(gapsValues[k]) / originalNumberOfSequences) >= .500)
                    file << "<span class=c10> </span>";
                else if (1 - (float(gapsValues[k]) / originalNumberOfSequences) >= .350)
                    file << "<span  class=c9> </span>";
                else if (1 - (float(gapsValues[k]) / originalNumberOfSequences) >= .250)
                    file << "<span  class=c8> </span>";
                else if (1 - (float(gapsValues[k]) / originalNumberOfSequences) >= .200)
                    file << "<span  class=c7> </span>";
                else if (1 - (float(gapsValues[k]) / originalNumberOfSequences) >= .150)
                    file << "<span  class=c6> </span>";
                else if (1 - (float(gapsValues[k]) / originalNumberOfSequences) >= .100)
                    file << "<span  class=c5> </span>";
                else if (1 - (float(gapsValues[k]) / originalNumberOfSequences) >= .050)
                    file << "<span  class=c4> </span>";
                else if (1 - (float(gapsValues[k]) / originalNumberOfSequences) >= .001)
                    file << "<span  class=c3> </span>";
                else
                    file << "<span  class=c2> </span>";
        }
        if (simValues != nullptr) {
            file << endl << setw(minHTML) << left << "    Similarity Scores:  ";
            for (k = j; ((k < originalNumberOfResidues) && (k < (j + HTMLBLOCKS))); k++)
                if (simValues[k] == 1)
                    file << "<span class=c12> </span>";
                else if (simValues[k] == 0)
                    file << "<span class=c1> </span>";
                else if (simValues[k] >= .750)
                    file << "<span class=c11> </span>";
                else if (simValues[k] >= .500)
                    file << "<span class=c10> </span>";
                else if (simValues[k] >= .250)
                    file << "<span  class=c9> </span>";
                else if (simValues[k] >= .100)
                    file << "<span  class=c8> </span>";
                else if (simValues[k] >= .010)
                    file << "<span  class=c7> </span>";
                else if (simValues[k] >= .001)
                    file << "<span  class=c6> </span>";
                else if (simValues[k] >= 1e-4)
                    file << "<span  class=c5> </span>";
                else if (simValues[k] >= 1e-5)
                    file << "<span  class=c4> </span>";
                else if (simValues[k] >= 1e-6)
                    file << "<span  class=c3> </span>";
                else
                    file << "<span  class=c2> </span>";
        }
        if (consValues != nullptr) {
            file << endl << setw(minHTML) << left << "    Consistency Scores: ";
            for (k = j; ((k < originalNumberOfResidues) && (k < (j + HTMLBLOCKS))); k++)
                if (consValues[k] == 1)
                    file << "<span class=c12> </span>";
                else if (consValues[k] == 0)
                    file << "<span class=c1> </span>";
                else if (consValues[k] >= .750)
                    file << "<span class=c11> </span>";
                else if (consValues[k] >= .500)
                    file << "<span class=c10> </span>";
                else if (consValues[k] >= .350)
                    file << "<span  class=c9> </span>";
                else if (consValues[k] >= .250)
                    file << "<span  class=c8> </span>";
                else if (consValues[k] >= .200)
                    file << "<span  class=c7> </span>";
                else if (consValues[k] >= .150)
                    file << "<span  class=c6> </span>";
                else if (consValues[k] >= .100)
                    file << "<span  class=c5> </span>";
                else if (consValues[k] >= .050)
                    file << "<span  class=c4> </span>";
                else if (consValues[k] >= .001)
                    file << "<span  class=c3> </span>";
                else
                    file << "<span  class=c2> </span>";
        }
        file << endl;
    }

    // Print HTML footer into output file
    file << "    </pre>" << endl << "  </body>" << endl << "</html>" << endl;

    // Close output file and deallocate local memory
    file.close();

    return true;
}

bool Alignment::statSVG(const char *const destFile) {

    // Init SVG size variables
    int
        whiteboxWidth = 1300,
        whiteboxHeight = 650,

        grayboxWidth = 1500,
        grayboxHeight = 900;

    // Init ratios and relative sizes
    float
        legendRatio = 0.175F,

        widthRatio = 0.5F,
        heightRatio = 0.75F,

        whiteboxDeltaHeight
        = whiteboxHeight * 0.05F,

        fontSize = whiteboxHeight * 0.02F;

    // Init chart sizes and origins
    float
        originX,
        originY,
        chartWidth,
        chartHeight;

    originX = (grayboxWidth - whiteboxWidth) * widthRatio + (whiteboxDeltaHeight * 0.5F);
    originY = (grayboxHeight - whiteboxHeight) * heightRatio + (whiteboxHeight - whiteboxDeltaHeight * 0.5F);
    chartWidth = (whiteboxWidth * (1.F - legendRatio) - whiteboxDeltaHeight);
    chartHeight = -(whiteboxHeight - whiteboxDeltaHeight);

    // Open the file;
    ofstream file;
    file.open(destFile);
    if (!file) {
        return false;
    }

    // svg header
    file << "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" version=\"1.1\" "
         << "height=\"" << grayboxHeight << "\" "
         << "width=\"" << grayboxWidth << "\">" << "\n";

    // White box
    file << "<rect "
         << "x=\"" << (grayboxWidth - whiteboxWidth) * widthRatio << "\" "
         << "width=\"" << whiteboxWidth * (1.F - legendRatio) << "\" "
         << "y=\"" << (grayboxHeight - whiteboxHeight) * heightRatio << "\" "
         << "height=\"" << whiteboxHeight << "\" "
         << "style=\"fill:white; stroke:black; stroke-width:2\" "
         << "/>" << "\n";

    // Header text
    file << "<text text-anchor=\"middle\" "
         << "x=\"" << grayboxWidth * 0.5F << "\" "
         << "y=\"" << (grayboxHeight - whiteboxHeight) * heightRatio * 0.75F << "\" "
         << "font-size=\"" << (grayboxHeight - whiteboxHeight) * heightRatio * 10.F / filename.size() << "\" "
         << ">"
         << filename
         << "</text>" << "\n";
    {
        float   halfDeltaHeight = (whiteboxDeltaHeight * 0.5F),

                x1 = (grayboxWidth - whiteboxWidth) * widthRatio,
                x2 = (grayboxWidth - whiteboxWidth) * widthRatio + whiteboxWidth * (1.F - legendRatio),
                y1 = (grayboxHeight - whiteboxHeight) * heightRatio + halfDeltaHeight,
                y2 = (grayboxHeight - whiteboxHeight) * heightRatio + halfDeltaHeight,

                x3 = (whiteboxWidth * (1.F - legendRatio) - whiteboxDeltaHeight) * 0.1F;
        for (int xx = 0; xx < 11; xx++) {
            // Vertical lines
            file << "<line "
                 << "x1=\"" << x1 << "\" "
                 << "y1=\"" << y1 + (whiteboxHeight - whiteboxDeltaHeight) * (xx * 0.1F) << "\" "
                 << "x2=\"" << x2 << "\" "
                 << "y2=\"" << y2 + (whiteboxHeight - whiteboxDeltaHeight) * (xx * 0.1F) << "\" "
                 << "style=\"stroke:black;stroke-width:1\" "
                 << "stroke-dasharray=\"1, 1\" "
                 << "opacity=\"0.5\"/>" << "\n";

            // Labels
            file << "<text "
                 << "x=\"" << x1 * 0.95F << "\" "
                 << "y=\"" << y2 + (whiteboxHeight - whiteboxDeltaHeight) * (xx * 0.1F) + (0.25F * fontSize)  << "\" "
                 << "text-anchor=\"end\" "
                 << "xml:space=\"preserve\" "
                 << "font-size=\"" << fontSize << "\">"
                 << (10 - xx) / 10.F
                 << "</text>" << "\n";

            // Horizontal lines
            file << "<line "
                 << "x1=\"" << x1 + halfDeltaHeight + x3 * xx  << "\" "
                 << "y1=\"" << (grayboxHeight - whiteboxHeight) * heightRatio << "\" "
                 << "x2=\"" << x1 + halfDeltaHeight + x3 * xx  << "\" "
                 << "y2=\"" << (grayboxHeight - whiteboxHeight) * heightRatio + whiteboxHeight << "\" "
                 << "style=\"stroke:black;stroke-width:1\" "
                 << "stroke-dasharray=\"1, 1\" "
                 << "opacity=\"0.5\"/>" << "\n";

            // Labels
            file << "<text "
                 << "x=\"" << x1 + (whiteboxWidth * (1.F - legendRatio) - whiteboxDeltaHeight) * (xx * 0.1F) + (whiteboxDeltaHeight * 0.5F) << "\" "
                 << "y=\"" << (grayboxHeight - whiteboxHeight) * heightRatio + whiteboxHeight + fontSize * 1.5F << "\" "
                 << "text-anchor=\"middle\" "
                 << "xml:space=\"preserve\" "
                 << "font-size=\"" << fontSize << "\">"
                 << xx * 10 << " %"
                 << "</text>" << "\n";
        }

        // Print Lines
        int statsAddedCounter = 0;
        {
            float deltaHeigth = whiteboxHeight / (float) (2 + 1);
            deltaHeigth = std::min(whiteboxHeight * 0.12F, deltaHeigth);
            float height = whiteboxWidth * legendRatio * 0.1F;

            // Lambda method that accepts an array of values to order and plot,
            //      the name of the stat, and the color of the lines and dots.
            //
            // The method does not perform a copy of the values to re-order them
            //      and thus, this has to be performed before the call.
            //
            // This allows to transform the values and to prevent multiple
            //      allocations and deallocations of big arrays.
            //
            // The values on the array should be normalized prior to the call,
            //      otherwise, the chart can appear deformed.
            //
            // It is a lambda due to the fact it uses most of the variables
            //      that are internal to this method: sizes and ratios.
            auto addStat = [&] (float * values,
                                const string & name,
                                const string & color) -> void {

                // Sort the provided array
                utils::quicksort(values, 0, originalNumberOfResidues - 1);

                // Legend
                {
                    // Color box
                    file << "<rect "
                         << "x=\"" << (grayboxWidth - whiteboxWidth) * widthRatio + whiteboxWidth * (1.F - legendRatio) + whiteboxWidth * legendRatio * 0.1F << "\" "
                         << "y=\"" << (grayboxHeight - whiteboxHeight) * heightRatio + deltaHeigth * (statsAddedCounter + 1) + deltaHeigth * 0.5F - height * 0.5F - fontSize * 0.25F << "\" "
                         << "width=\"" << height << "\" "
                         << "height=\"" << height << "\" "
                         << "style=\"fill:" << color << "; stroke:black; stroke-width:2\" "
                         << "fill-opacity=\"0.75\" "
                         << "/>" << "\n";

                    // Stat name
                    file << "<text "
                         << "x=\"" << (grayboxWidth - whiteboxWidth) * widthRatio + whiteboxWidth * (1.F - legendRatio) + whiteboxWidth * legendRatio * 0.5F << "\" "
                         << "y=\"" << (grayboxHeight - whiteboxHeight) * heightRatio + deltaHeigth * (statsAddedCounter + 1) + deltaHeigth * 0.5F << "\" "
                         << "text-anchor=\"middle\" "
                         << "xml:space=\"preserve\" "
                         << "font-size=\"" << fontSize << "\">"
                         << name
                         << "</text>" << "\n";
                }

                // Start a polyline
                file << "<polyline stroke-linecap=\"round\" "
                     << "style=\"fill:none;stroke:" << color << ";stroke-width:0.8\" opacity=\"0.8\" points=\"";

                // Add each value individually
                for (int X = 0; X < numberOfResidues; X++)
                {
                    file << originX + ((float)X / numberOfResidues) * chartWidth << ",\t"
                         << originY + values[X] * chartHeight << " \n";
                }
                // Finish the polyline
                file << "\"/>" << "\n";

                // Add a dot for each value
                for (int X = 0; X < numberOfResidues; X++)
                {
                    file << "<circle cx=\""
                    << originX + ((float)X / numberOfResidues) * chartWidth
                    << "\" cy=\""
                    << (originY + values[X] * chartHeight)
                    << "\" r=\"2\" stroke=\"black\" stroke-width=\"0.1\" fill=\""
                    << color
                    << "\" />\n";
                }

                // Increase the number of stats added.
                statsAddedCounter ++;
            };


            // Create a vector to contain all the values for each column, for each stat.
            //      This allows to prevent multiple allocations for the copy of each stat.
            float * vectAux = new float[numberOfResidues];

            // Add the gap stat if calculated.
            if (Statistics->gaps != nullptr)
            {
                // Copy each value, and normalize using the original number of sequences.
                //      This is needed due to gapWindow being the total number of gaps per column.
                for (int i = 0; i < originalNumberOfResidues; i++)
                    vectAux[i] = Statistics->gaps->getGapsWindow()[i] / (float)originalNumberOfSequences;

                // Add the stat.
                addStat(vectAux, "Gaps", "Red");
            }

            // Add the similarity stat if calculated
            if (Statistics->similarity != nullptr)
            {
                // Make a copy of the values, to allow reordering without modification.
                for (int i = 0; i < originalNumberOfResidues; i++)
                    vectAux[i] = Statistics->similarity->getMdkWindowedVector()[i];

                addStat(vectAux, "Similarity", "Blue");
            }

            // Add the consistency stat if calculated
            if (Statistics->consistency != nullptr)
            {
                // Make a copy of the values, to allow reordering without modification.
                for (int i = 0; i < originalNumberOfResidues; i++)
                    vectAux[i] = Statistics->consistency->getValues()[i];

                addStat(vectAux, "Consistency", "Green");
            }

            // Delete the temporal values array
            delete [] vectAux;

        }

        // Legend
        {
            float deltaHeigth = whiteboxHeight / (float) (statsAddedCounter + 1);
            deltaHeigth = std::min(whiteboxHeight * 0.12F, deltaHeigth);

            // Legend box
            file << "<rect "
                 << "x=\"" << (grayboxWidth - whiteboxWidth) * widthRatio + whiteboxWidth * (1.F - legendRatio) + 6 << "\" "
                 << "width=\"" << whiteboxWidth * legendRatio << "\" "
                 << "y=\"" << (grayboxHeight - whiteboxHeight) * heightRatio << "\" "
                 << "height=\"" << deltaHeigth * (statsAddedCounter + 1) << "\" "
                 << "style=\"fill:white; stroke:black; stroke-width:2\" "
                 << "fill-opacity=\"0.25\" "
                 << "/>" << "\n";

            file << "<text "
                 << "x=\"" << (grayboxWidth - whiteboxWidth) * widthRatio + whiteboxWidth * (1.F - legendRatio) + whiteboxWidth * legendRatio * 0.5F << "\" "
                 << "y=\"" << (grayboxHeight - whiteboxHeight) * heightRatio + deltaHeigth * 0 + deltaHeigth * 0.5F << "\" "
                 << "text-anchor=\"middle\" "
                 << "xml:space=\"preserve\" "
                 << "font-size=\"" << fontSize * 2 << "\">"
                 << "statistics"
                 << "</text>" << "\n";

            file << "<line "
                 << "x1=\"" << (grayboxWidth - whiteboxWidth) * widthRatio + whiteboxWidth * (1.F - legendRatio) + 12 << "\" "
                 << "x2=\"" << (grayboxWidth - whiteboxWidth) * widthRatio + whiteboxWidth * (1.F - legendRatio) + whiteboxWidth * legendRatio * 1.F << "\" "
                 << "y1=\"" << (grayboxHeight - whiteboxHeight) * heightRatio + deltaHeigth * 0.3F + deltaHeigth * 0.5F << "\" "
                 << "y2=\"" << (grayboxHeight - whiteboxHeight) * heightRatio + deltaHeigth * 0.3F + deltaHeigth * 0.5F << "\" "
                 << "style=\"stroke:black;stroke-width:2\" />" << "\n";
        }

    }


    file << "</svg>";
    file.close();

    return true;
}

bool Alignment::alignmentSummarySVG(Alignment &trimmedAlig, const char *const destFile, int blocks) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool Alignment::alignmentSummarySVG(Alignment &trimmedAlig, char *destFile, float *consValues, int blocks) ");

    int i, j, k, counter = 0;

    // BEGIN Init variables
    int *gapsValues = nullptr;
    if (Statistics->gaps != nullptr)
        gapsValues = Statistics->gaps->getGapsWindow();

    float *simValues = nullptr;
    if (Statistics->similarity != nullptr)
        simValues = Statistics->similarity->getMdkWindowedVector();

    float *consValues = nullptr;
    if (Statistics->consistency != nullptr)
        consValues = Statistics->consistency->getValues();

    // Check if alignment is aligned;
    if (!isAligned) {
        debug.report(ErrorCode::NotAligned, new std::string[1]{filename});
        return false;
    }

    // Calculate the blockSize;
    blocks = 120;

    int fontSize = 15;

    // Open the file;
    ofstream file;
    file.open(destFile);
    if (!file) {
        return false;
    }

    // Allocate some local memory
    char *tmpColumn = new char[originalNumberOfSequences + 1];
    tmpColumn[originalNumberOfSequences] = '\0';

    // Compute HTML blank spaces
    j = 0;
    for (i = 0; i < numberOfSequences; i++)
        j = utils::max(j, seqsName[i].size());

    int sequencesNamesLength = utils::max(25, j);

    // Init Colors
    std::map<char, string> mappedColors = std::map<char, string>{
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

    std::map<char, string> mappedColorsMeaning = std::map<char, string>{
            {'o', "Glycines"},
            {'y', "Prolines"},
            {'b', "Hidrophobic"},
            {'w', "Unconserved"},
            {'p', "Prolines"},
            {'r', "Positive Charge"},
            {'g', "Polar"},
            {'m', "Negative Charge"},
            {'c', "Aromatic"}
    };

    // END

    // BEGIN Init SVG
    float leftMargin = 5, nameBlocksMargin = 15,
            rightMargin = 5, topMargin = 5,
            bottomMargin = 10, interBlocksMargin = 20;

    float height = ((originalNumberOfSequences + 3 + 7) // ((gapsValues != nullptr || simValues != nullptr) ? 7 : 0) /* Col numbering occupies 3 rows. Stats 5 */)
                    * fontSize + interBlocksMargin)                 // Height of a block + interblock
                   * std::ceil((float) originalNumberOfResidues / blocks)  // Number of blocks
                   - interBlocksMargin // Number of interblocks is number of blocks  - 1
                   + fontSize * 16
                   + topMargin
                   + bottomMargin;

    float width = sequencesNamesLength * fontSize
                  + std::min(blocks, trimmedAlig.originalNumberOfResidues) * fontSize * 0.75F
                  + leftMargin
                  + nameBlocksMargin
                  + rightMargin;

    float currentHeight = topMargin + fontSize;

    // Start the svg output
    file << "<svg xmlns=\"http://www.w3.org/2000/svg\" \
            xmlns:xlink=\"http://www.w3.org/1999/xlink\" \
            version=\"1.1\" \
            viewbox=\"0 0 " << width << " " << height << "\">" << endl;
    // END

    // BEGIN Functions
    file << "<script type=\"text/ecmascript\"> \n\
                <![CDATA[ \n\
                function onMouseOverLabel(evt) { \n\
                    var Line = document.getElementById(evt.target.id.replace(\"Label\", \"Line\")); \n\
                    Line.setAttribute(\"style\", \"fill:none;stroke:black;stroke-width:3\"); \n\
                    evt.target.setAttribute(\"font-size\", 17); \n\
                } \n\
                function onMouseOutLabel(evt) { \n\
                    var Line = document.getElementById(evt.target.id.replace(\"Label\", \"Line\")); \n\
                    Line.setAttribute(\"style\", \"fill:none;stroke:black;stroke-width:1\"); \n\
                    evt.target.setAttribute(\"font-size\", 15); \n\
                } \n\
                function onMouseOverLine(evt) { \n\
                    var Label = document.getElementById(evt.target.id.replace(\"Line\", \"Label\")); \n\
                    evt.target.setAttribute(\"style\", \"fill:none;stroke:black;stroke-width:3\"); \n\
                    Label.setAttribute(\"font-size\", 17); \n\
                } \n\
                function onMouseOutLine(evt) { \n\
                    var Label = document.getElementById(evt.target.id.replace(\"Line\", \"Label\")); \n\
                    evt.target.setAttribute(\"style\", \"fill:none;stroke:black;stroke-width:1\"); \n\
                    Label.setAttribute(\"font-size\", 15); \n\
                } \n\
                ]]> \n </script>" << endl;
    // END

    // BEGIN Patterns

    // Selected Block
    file << "<pattern id=\"selected_pattern\" patternUnits=\"userSpaceOnUse\" width=\"10\" height=\"10\"> " << endl;
    file << "<rect width=\"10\" height=\"10\" fill=\"green\" opacity=\"0.5\"/>";
    file << "<path d=\"M-1,1 l2,-2 M0,10 l10,-10 M9,11 l2,-2\" stroke=\"white\" stroke-width=\"1\"/>";
    file << "</pattern> " << endl;

    // Deleted Block
    file << "<pattern id=\"rejected\" patternUnits=\"userSpaceOnUse\" width=\"10\" height=\"10\"> " << endl;
    file << "<rect width=\"10\" height=\"10\" fill=\"red\" opacity=\"0.5\"/>";
    file << "<path d=\"M2,-2 l-1,1 M10,-10 l0,10 M2,-2 l9,11\" stroke=\"white\" stroke-width=\"1\"/>";
    file << "</pattern> " << endl;

    // END

    // BEGIN Header



    file << "<rect \
                style=\"fill:none;stroke:lightgrey\" \
                height=\"" << fontSize * 6.5F << "\" \
                width=\"" << fontSize * std::max(30.5F, filename.length() * 0.75F + 1.F) << "px\" \
                x =\"" << leftMargin + nameBlocksMargin - fontSize * 0.5F << "px\"\
                y =\"" << currentHeight - fontSize * 0.5F << "\"/>" << endl;

    file << "<text \
                font-family = \"monospace\" \
                font-size=\"" << fontSize << "px\" \
                dy=\".35em\" \
                x =\"" << leftMargin + nameBlocksMargin << "\" \
                y=\"" << (currentHeight + fontSize / 2) << "\" \
                kerning=\"0\" \
                text-anchor=\"start\" \
                lengthAdjust=\"spacingAndGlyphs\"\
                textLength=\"" << filename.size() * 0.75F * fontSize << "\" \
                style=\"font-weight:100\" \
                xml:space=\"preserve\" >"

         << filename

         << "</text>" << endl;
    currentHeight += fontSize * 2;
    filename = "Selected Sequences \t" + std::to_string(trimmedAlig.numberOfSequences) + " / " + std::to_string(originalNumberOfSequences);
    file << "<text \
                font-family = \"monospace\" \
                font-size=\"" << fontSize << "px\" \
                dy=\".35em\" \
                x =\"" << leftMargin + nameBlocksMargin << "\" \
                y=\"" << (currentHeight + fontSize / 2) << "\" \
                kerning=\"0\" \
                text-anchor=\"start\" \
                lengthAdjust=\"spacingAndGlyphs\"\
                textLength=\"" << std::min((float) filename.size() * 0.65F, (float) sequencesNamesLength) * fontSize << "\" \
                style=\"font-weight:100\" \
                xml:space=\"preserve\" >"

         << filename

         << "</text>" << endl;
    currentHeight += fontSize * 2;
    filename = "Selected Residues \t\t" + std::to_string(trimmedAlig.numberOfResidues) + " / " + std::to_string(originalNumberOfResidues);
    file << "<text \
                font-family = \"monospace\" \
                font-size=\"" << fontSize << "px\" \
                dy=\".35em\" \
                x =\"" << leftMargin + nameBlocksMargin << "\" \
                y=\"" << (currentHeight + fontSize / 2) << "\" \
                kerning=\"0\" \
                text-anchor=\"start\" \
                lengthAdjust=\"spacingAndGlyphs\"\
                textLength=\"" << std::min((float) filename.size() * 0.65F, (float) sequencesNamesLength) * fontSize << "\" \
                style=\"font-weight:100\" \
                xml:space=\"preserve\" >"

         << filename

         << "</text>" << endl;
    currentHeight += fontSize * 3;

    // END

    // BEGIN Legend
    file << "<rect \
                style=\"fill:none;stroke:lightgrey\" \
                height=\"" << fontSize * 6.5F << "\" \
                width=\"" << fontSize * 30.5 << "px\" \
                x =\"" << leftMargin + nameBlocksMargin - fontSize * 0.5F << "px\"\
                y =\"" << currentHeight - fontSize * 0.5F << "\"/>" << endl;

    std::map<char, string>::iterator it = mappedColors.begin();
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++, it++) {
            file << "<rect \
                    style=\"fill:" << it->second << ";stroke:black;stroke-width:1px\" \
                    height=\"" << fontSize << "\" \
                    width=\"" << fontSize << "px\" \
                    x =\"" << leftMargin + nameBlocksMargin * (j + 1) + (9 * fontSize * j) << "px\"\
                    y =\"" << currentHeight + fontSize * i * 1.2F << "\"/>" << endl;

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
    for (j = 0; j < originalNumberOfResidues; j += blocks) {
        // BEGIN Columns Numbering Numbers
        file << "<text \
                font-family = \"monospace\" \
                font-size=\"" << fontSize << "px\" \
                dy=\".35em\" \
                x =\"" << leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize << "\" \
                y=\"" << (currentHeight + fontSize / 2) << "\" \
                kerning=\"0\" \
                text-anchor=\"start\" \
                lengthAdjust=\"spacingAndGlyphs\"\
                textLength=\"" << blocks * fontSize * 0.75F << "\" \
                style=\"font-weight:100\" \
                xml:space=\"preserve\" >";
        for (k = j; k < originalNumberOfResidues && (k - j) < blocks; k += 10) {
            file << std::left << std::setw(10) << std::setfill(' ') << k;
        }
        for (; (k - j) < blocks; k += 10)
            file << std::left << std::setw(10) << std::setfill(' ') << " ";

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
                textLength=\"" << blocks * fontSize * 0.75F << "\" \
                style=\"font-weight:100\" \
                xml:space=\"preserve\" >";
        for (k = j; k < originalNumberOfResidues && (k - j) < blocks; k += 10) {
            file << std::left << std::setw(10) << std::setfill(' ') << "|";
        }
        for (; (k - j) < blocks; k += 10)
            file << std::left << std::setw(10) << std::setfill(' ') << " ";

        file << "</text>" << endl;

        file << "<line x1=\"" << leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize
             << "\" x2=\"" << leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize
                              + std::min((int) blocks, originalNumberOfResidues - j) * fontSize * 0.75F
             << "\" y1=\"" << currentHeight + fontSize / 2
             << "\" y2=\"" << currentHeight + fontSize / 2
             << "\" style=\"stroke:rgb(0,0,0);stroke-width:2\"/>" << endl;
        currentHeight += fontSize;


        currentHeight += fontSize;

        // END

        // BEGIN Color Matrix
        char lastColor, newColor;
        int lastPos = 0;
        for (i = 0; i < blocks && (i + j) < originalNumberOfResidues; i++) {
            for (k = 0; k < originalNumberOfSequences; k++) {
                tmpColumn[k] = (sequences[k][i + j]);
            }
            lastColor = utils::determineColor(tmpColumn[0], tmpColumn);
            lastPos = 0;
            for (k = 1; k < originalNumberOfSequences; k++) {
                newColor = utils::determineColor(tmpColumn[k], tmpColumn);
                if (newColor != lastColor) {
                    if (lastColor != 'w')
                        file << "<rect " <<
                             " style=\"fill:" << mappedColors[lastColor] << "\" " <<
                             " height=\"" << fontSize * (k - lastPos) << "\" " <<
                             " width=\"" << fontSize * 0.75F << "px\" " <<
                             " x =\"" <<

                             leftMargin +
                             nameBlocksMargin +
                             sequencesNamesLength * fontSize +
                             i * 0.75F * fontSize

                             << "\" " <<
                             " y =\"" << currentHeight + fontSize * lastPos << "\" />" << endl;
                    lastColor = newColor;
                    lastPos = k;
                }
            }
            if (lastColor != 'w')
                file << "<rect " <<
                     " style=\"fill:" << mappedColors[lastColor] << "\" " <<
                     " height=\"" << fontSize * (k - lastPos) << "\" " <<
                     " width=\"" << fontSize * 0.75F << "px\" " <<
                     " x =\"" <<

                     leftMargin + nameBlocksMargin + (sequencesNamesLength) * fontSize + i * 0.75F * fontSize

                     << "\" " <<
                     " y =\"" << currentHeight + fontSize * lastPos << "\" />" << endl;
        }

        // END

        // BEGIN Sequences Names and Sequences

        for (i = 0; i < originalNumberOfSequences; i++) {
            // Print names
            file << "<text \
                font-family = \"monospace\" \
                font-size=\"" << fontSize << "px\" \
                dy=\".35em\" \
                x =\"" << leftMargin << "\" \
                y=\"" << (currentHeight + fontSize / 2) << "\" \
                kerning=\"0\" \
                text-anchor=\"start\" \
                textLength=\"" << std::min((float) seqsName[i].size() * 0.75F, (float) sequencesNamesLength) * fontSize << "\" \
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
                textLength=\"" << std::min((int) blocks, originalNumberOfResidues - j) * fontSize * 0.75F << "\" \
                style=\"font-weight:100\"" << ">"

                 << sequences[i].substr(j, blocks)

                 << "</text>" << endl;
            if (trimmedAlig.saveSequences[i] == -1) {
                file << "<line x1=\"" << leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize << "\"" <<
                     " x2=\"" << leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize + std::min(blocks, originalNumberOfResidues - j) * fontSize * 0.75F << "\"" <<
                     " y1=\"" << (currentHeight + fontSize / 2) << "\"" <<
                     " y2=\"" << (currentHeight + fontSize / 2) << "\"" <<
                     " style=\"stroke:rgb(0,0,0);stroke-width:2\""
                     "/>" << endl;
                file << "<line x1=\"" << leftMargin << "\"" <<
                     " x2=\"" << leftMargin + std::min((float) seqsName[i].size() * 0.75F, (float) sequencesNamesLength) * fontSize << "\"" <<
                     " y1=\"" << (currentHeight + fontSize / 2) << "\"" <<
                     " y2=\"" << (currentHeight + fontSize / 2) << "\"" <<
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
             " height=\"" << fontSize * 5 << "\" " <<
             " width=\"" << fontSize * 0.75F * std::min(blocks, originalNumberOfResidues - j) << "px\" " <<
             " x =\"" <<

             leftMargin + nameBlocksMargin + (sequencesNamesLength) * fontSize

             << "\" " <<
             " y =\"" << currentHeight << "\" />" << endl;

        bool rejected = trimmedAlig.saveResidues[j] == -1;
        lastPos = j;
        for (k = j + 1; k < originalNumberOfResidues && (k - j) < blocks; k++) {
            if ((trimmedAlig.saveResidues[k] == -1) != rejected) {
                file << "<rect " <<
                     " style=\"fill:url(#" << (rejected ? "rejected" : "selected_pattern") << ");stroke:black\"" <<
                     " height=\"" << fontSize * 5 << "\" " <<
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
             " style=\"fill:url(#" << (rejected ? "rejected" : "selected_pattern") << ");stroke:black\"" <<
             " height=\"" << fontSize * 5 << "\" " <<
             " width=\"" << fontSize * 0.75F * (k - lastPos) << "px\" " <<
             " opacity=\"0.5\""
             " x =\"" <<

             leftMargin + nameBlocksMargin + sequencesNamesLength * fontSize + fontSize * 0.75F * (lastPos - j)

             << "\" " <<
             " y =\"" << currentHeight << "\" />" << endl;
        // END

        // BEGIN Stats report
        if (gapsValues || simValues || consValues) {
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
        if (gapsValues) {
            file << "<polyline \
                        style=\"fill:none;stroke:black;stroke-width:1\" \
                        stroke-dasharray=\"10,4\" \
                        id =\"GapsLine" << counter << "\" \
                        onmouseover=\"onMouseOverLine(evt)\" \
                        onmouseout=\"onMouseOutLine(evt)\" \
                        points=\"";

            for (i = 0; i < blocks && (i + j) < originalNumberOfResidues; i++) {
                file << leftMargin + nameBlocksMargin + (sequencesNamesLength) * fontSize + (i + 0.5F) * 0.75F * fontSize << ","
                     << currentHeight + 10 + (1.F - (gapsValues[i + j] / (float) originalNumberOfSequences)) * fontSize * 4 << " ";
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
                        textLength=\"" << std::min(11 * 0.75F, (float) sequencesNamesLength) * fontSize << "\" \
                        style=\"font-weight:bold\"" << ">"

                 << "Gaps Values"

                 << "</text>" << endl;

            file << "<line x1=\"" << leftMargin
                 << "\" x2=\"" << leftMargin + std::min(11 * 0.75F, (float) sequencesNamesLength) * fontSize
                 << "\" y1=\"" << currentHeight + fontSize * 1.25F
                 << "\" y2=\"" << currentHeight + fontSize * 1.25F
                 << "\" style=\"stroke:rgb(0,0,0);stroke-width:1\""
                 << " stroke-dasharray=\"10,4\""
                 << " id =\"GapsSubLine" << counter << "\""
                 << "/>" << endl;
        }

        if (simValues) {
            file << "<polyline \
                style=\"fill:none;stroke:black;stroke-width:1\" \
                stroke-dasharray=\"5,1\" \
                id =\"SimLine" << counter << "\" \
                onmouseover=\"onMouseOverLine(evt)\" \
                onmouseout=\"onMouseOutLine(evt)\" \
                points=\"";

            for (i = 0; i < blocks && (i + j) < originalNumberOfResidues; i++) {
                file << leftMargin + nameBlocksMargin + (sequencesNamesLength) * fontSize + (i + 0.5F) * 0.75F * fontSize << ","
                     << currentHeight + 10 + ((1.F - simValues[i + j]) * fontSize * 4) << " ";
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
                        textLength=\"" << std::min(10 * 0.75F, (float) sequencesNamesLength) * fontSize << "\" \
                        style=\"font-weight:bold\"" << ">"

                 << "Sim Values"

                 << "</text>" << endl;

            file << "<line x1=\"" << leftMargin
                 << "\" x2=\"" << leftMargin + std::min(10 * 0.75F, (float) sequencesNamesLength) * fontSize
                 << "\" y1=\"" << currentHeight + fontSize * 3.25F
                 << "\" y2=\"" << currentHeight + fontSize * 3.25F
                 << "\" style=\"stroke:rgb(0,0,0);stroke-width:1\""
                 << " stroke-dasharray=\"5,1\""
                 << " id =\"SimSubLine" << counter << "\""
                 << "/>" << endl;
        }

        if (consValues) {
            file << "<polyline \
                style=\"fill:none;stroke:black;stroke-width:1\" \
                stroke-dasharray=\"2,2\" \
                id =\"ConsLine" << counter << "\" \
                onmouseover=\"onMouseOverLine(evt)\" \
                onmouseout=\"onMouseOutLine(evt)\" \
                points=\"";

            for (i = 0; i < blocks && (i + j) < originalNumberOfResidues; i++) {
                file << leftMargin + nameBlocksMargin + (sequencesNamesLength) * fontSize + (i + 0.5F) * 0.75F * fontSize << ","
                     << currentHeight + 10 + ((1.F - consValues[i + j]) * fontSize * 4) << " ";
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
                        textLength=\"" << std::min(10 * 0.75F, (float) sequencesNamesLength) * fontSize << "\" \
                        style=\"font-weight:bold\"" << ">"

                 << "Cons Values"

                 << "</text>" << endl;

            file << "<line x1=\"" << leftMargin
                 << "\" x2=\"" << leftMargin + std::min(10 * 0.75F, (float) sequencesNamesLength) * fontSize
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

    delete[] tmpColumn;

    return true;
}

void Alignment::updateSequencesAndResiduesNums(bool countSequences, bool countResidues) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void Alignment::updateSequencesAndResiduesNums(bool countSequences, bool countResidues) ");
    int i;
    if (countSequences) {
        for (numberOfSequences = 0, i = 0; i < originalNumberOfSequences; i++) {
            if (saveSequences[i] != -1) numberOfSequences++;
        }
    }

    if (countResidues) {
        for (numberOfResidues = 0, i = 0; i < originalNumberOfResidues; i++) {
            if (saveResidues[i] != -1) numberOfResidues++;
        }
    }
}

