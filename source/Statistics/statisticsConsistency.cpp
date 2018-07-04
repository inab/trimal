#include <TimerFactory.h>
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    2009-2015 Capella-Gutierrez S. and Gabaldon, T.
              [scapella, tgabaldon]@crg.es

    This file is part of trimAl.

    trimAl is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl. If not, see <http://www.gnu.org/licenses/>.

***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#include "Statistics/statisticsConsistency.h"
#include "reportsystem.h"
#include <sstream>

#define LONG 80



void statisticsConsistency::perform(char *comparesetFilePath,
                                    ReadWriteMS &ReadWriteMachine,
                                    trimAlManager &manager,
                                    char *forceFile) {

    line = new char[1024];

    // Open the file that contains the paths of the files.
    compare.open(comparesetFilePath, ifstream::in);
    while (compare.getline(line, 1024)) numFiles++;
    compare.close();
    compare.open(comparesetFilePath);

    compareAlignmentsArray = new newAlignment *[numFiles];
    filesToCompare = new char *[numFiles];

    // Load all the alignments to compare
    // Check if they: are aligned and the type is the same
    // Store the maximum number of amino acids present on all alignments
    for (i = 0; i < numFiles; i++) {


        // Search for end of line.
        for (nline.clear(), compare.read(&c, 1);
             (c != '\n') && ((!compare.eof()));
             compare.read(&c, 1))
        {
            nline += c;
        }

        // Save the alignment path
        filesToCompare[i] = new char[nline.size() + 1];
        strcpy(filesToCompare[i], nline.c_str());

        // Load the alignment
        compareAlignmentsArray[i]
                = ReadWriteMachine.loadAlignment(filesToCompare[i]);

        // Check if alignment could be loaded
        if (compareAlignmentsArray[i] == nullptr) {
            appearErrors = true;
        } else {
            // Check if alignment is not aligned
            if (!compareAlignmentsArray[i]->isFileAligned()) {
                debug.report(ErrorCode::NotAligned,
                             new std::string[1]{filesToCompare[i]});
                appearErrors = true;

            } else {
                compareAlignmentsArray[i]->SequencesMatrix
                        = new sequencesMatrix(compareAlignmentsArray[i]);

                // Store maximum number of aminoacids
                if (compareAlignmentsArray[i]->getNumAminos() > maxAminos)
                    maxAminos = compareAlignmentsArray[i]->getNumAminos();

                // Check if alignment type is the same as last one.
                if (prevType == -1)
                    prevType = compareAlignmentsArray[i]->getAlignmentType();
                else if (compareAlignmentsArray[i]->getAlignmentType() != prevType) {
                    debug.report(ErrorCode::AlignmentTypesNotMatching);
                    appearErrors = true;
                }
            }
        }
    }

    // If the analysis couldn't be performed, stop the program.
    if (appearErrors)
    {
        debug.report(ErrorCode::ComparesetFailedAlignmentMissing);
        delete_variables();
        delete [] values;
        exit(ErrorCode::ComparesetFailedAlignmentMissing);
    }

    else {
        // If no alignment is forced to be selected, select one of them
        if (forceFile == nullptr) {
            values = new float[maxAminos];
            // Perform stat calculation and
            //  choice the best scoring alignment
            referFile = statisticsConsistency::compareAndChoose(
                    compareAlignmentsArray,
                    filesToCompare,
                    values,
                    numFiles,
                    // Verbosity
                    manager.stats >= 0 && manager.outfile != nullptr);

            // If no alignment could be selected, stop the program
            if (referFile == -1)
            {
                delete_variables();
                delete [] values;
                exit(-1);
            }

            // Specify the selected alignment as origAlig
            //  (as if it was fed thought -in argument
            manager.origAlig
                    = new newAlignment(*compareAlignmentsArray[referFile]);
        }
            // If there is an alignment that has been forcibly selected
        else {
            values = new float[manager.origAlig->getNumAminos()];
            appearErrors = !statisticsConsistency::forceComparison(
                    compareAlignmentsArray,
                    numFiles,
                    manager.origAlig,
                    values);

            // If forcing comparison failed, exit the program
            if (appearErrors)
            {
                delete_variables();
                delete [] values;
                exit(-1);
            }
        }

        // Store the cross reference between alignment
        //  and this, it's consistency stat
        manager.origAlig->Statistics->consistency = this;
        _alignment = manager.origAlig;
        residues = _alignment->originalResidNumber;


        // Apply window sizes
        if (manager.windowSize != -1)
            appearErrors += !applyWindow(manager.windowSize);
        else if (manager.consistencyWindow != -1)
            appearErrors += !applyWindow(manager.consistencyWindow);

        if (appearErrors)
        {
            exit(-1);
        }

        // If no output format is provided
        //  we'll use the selected alignment format
        if (manager.oformats.empty()) {
            manager.oformats.emplace_back(
                    ReadWriteMachine.getFileFormatName(
                            forceFile == nullptr ? filesToCompare[referFile] : forceFile)
            );
        }

        delete_variables();
    }
}


// This method compares a set of alignment in order to select the most
// consistent one respect of the other ones. To compute the consistency
// values we use the proportion of residue pairs per column in the aligs
// to compare 

int statisticsConsistency::compareAndChoose(newAlignment **vectAlignments,
                                            char **fileNames,
                                            float *columnsValue,
                                            int numAlignments,
                                            bool verbosity) {
	 // Create a timerLevel that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("int statisticsConsistency::compareAndChoose("
                        "newAlignment **vectAlignments, "
                        "char **fileNames, "
                        "float *columnsValue, "
                        "int numAlignments, "
                        "bool verbosity) ");

    int *numResiduesAlig, *correspNames, *columnSeqMatrix, *columnSeqMatrixAux;
    int i, j, k, l, m, numSeqs, pairRes, hits, alig = 0;
    float max = 0, value = 0, **vectHits;
    bool appearErrors = false;
    string *names;


    // Get some parameters from the alignment that has
    // been selected */
    numSeqs = vectAlignments[0]->getNumSpecies();

    // Allocate dinamic local memory
    names               = new string[numSeqs];
    correspNames        = new int[numSeqs];
    numResiduesAlig     = new int[numAlignments];
    columnSeqMatrix     = new int[numSeqs];
    vectHits            = new float *[numAlignments];
    columnSeqMatrixAux  = new int[numSeqs];

    // Check that all of alignment has the same number of
    // sequence as well as there exists a correspondence
    // between the names for each pars of aligs.
    for (i = 1; i < numAlignments; i++) {

        if (numSeqs != vectAlignments[i]->getNumSpecies()) {
            debug.report(ErrorCode::DifferentNumberOfSequencesInCompareset);
            appearErrors = true;
            break;
        }

        vectAlignments[i]->getSequences(names);
        if (!vectAlignments[0]->getSequenceNameOrder(names, correspNames)) {
            debug.report(ErrorCode::DifferentSeqsNamesInCompareset);
            appearErrors = true;
            break;
        }
    }

    if (!appearErrors)
    {

        // Changes the order in sequences number matrix
        // according to the order in the selected alignment
        for (i = 1; ((i < numAlignments) && (!appearErrors)); i++) {
            vectAlignments[i]->getSequences(names);
            vectAlignments[0]->getSequenceNameOrder(names, correspNames);
            vectAlignments[i]->SequencesMatrix->setOrder(correspNames);
        }

        // Get back the residues number for each alignment
        for (i = 0; ((i < numAlignments) && (!appearErrors)); i++)
            numResiduesAlig[i] = vectAlignments[i]->getNumAminos();

        // Start the comparison among the alignments
        for (i = 0; ((i < numAlignments) && (!appearErrors)); i++) {
            value = 0;
            // If it's necessary, we print some information
            if (verbosity)
                cout << endl;

            // Initialize the hits vector for each alignment
            vectHits[i] = new float[numResiduesAlig[i]];
            utils::initlVect(vectHits[i], numResiduesAlig[i], 0);

            for (j = 0, pairRes = 0, hits = 0; j < numResiduesAlig[i]; j++, pairRes = 0, hits = 0) {


                // Get back each column from the current selected
                // alignment
                vectAlignments[i]->SequencesMatrix->getColumn(j, columnSeqMatrix);

                // For each position from the previous recovered
                // columns, we carry out the analysis to detect the
                // same residues pair in the rest of the alignmetns
                for (k = 0; k < numSeqs; k++) {

                    // If there is a valid residue, we go ahead with
                    // the analysis
                    if (columnSeqMatrix[k] != 0) {
                        for (l = 0; l < i; l++) {
                            // Recover the residue pairs from the others aligs
                            vectAlignments[l]->SequencesMatrix->getColumn(columnSeqMatrix[k], k, columnSeqMatrixAux);
                            // and count the similar residue pairs
                            for (m = k + 1; m < numSeqs; m++)
                                if (columnSeqMatrix[m] != 0) {
                                    if (columnSeqMatrix[m] == columnSeqMatrixAux[m])
                                        hits++;
                                    pairRes++;
                                }
                        }

                        for (l = i + 1; l < numAlignments; l++) {
                            // Recover the residue pairs from the others aligs
                            vectAlignments[l]->SequencesMatrix->getColumn(columnSeqMatrix[k], k, columnSeqMatrixAux);
                            // and count the similar residue pairs
                            for (m = k + 1; m < numSeqs; m++)
                                if (columnSeqMatrix[m] != 0) {
                                    if (columnSeqMatrix[m] == columnSeqMatrixAux[m])
                                        hits++;
                                    pairRes++;
                                }
                        }
                    }
                }

                // For each column, compute the hits proportion for
                // every residue pair against the rest of alignments
                if (pairRes != 0) {
                    vectHits[i][j] += ((1.0 * hits) / pairRes);
                    value += vectHits[i][j];
                }
            }

            // The method can offer some information about the
            // comparison progression
            if (verbosity)
                cout << "File:\t\t" << fileNames[i] << endl << "Values:\t\tSequences: " << numSeqs
                     << "\tResidues: " << numResiduesAlig[i] << "\tPond. Hits: " << setw(8)
                     << value << "\t%Consistency: " << value / numResiduesAlig[i] << endl;
            // Keep the alignment with higher consistency value
            if ((value / numResiduesAlig[i]) > max) {
                alig = i;
                max = value / numResiduesAlig[i];
            }
        }

        // Prints the alignment that have been selected
        if ((verbosity) && (!appearErrors)) {
            cout << "\t\t\t\t\t--------------" << endl;
            cout << endl << "File Selected:\t" << fileNames[alig] << endl << "Value:\t\t" << max << endl << endl;
        }

        // The method returns a vector with the consistency
        // value for each column in the selected alignment
        if ((columnsValue != nullptr) && (!appearErrors)) {
            utils::initlVect(columnsValue, numResiduesAlig[alig], -1);
            for (i = 0; i < numResiduesAlig[alig]; i++)
                columnsValue[i] = vectHits[alig][i];
        }
    }

    // Deallocate memory
    for (i = 0; ((i < numAlignments) && (!appearErrors)); i++)
        delete[] vectHits[i];

    delete[] vectHits;
    delete[] names;
    delete[] correspNames;
    delete[] numResiduesAlig;
    delete[] columnSeqMatrix;
    delete[] columnSeqMatrixAux;

    // Return the selected alignment index or an error
    // flag otherwise
    if (appearErrors) return -1;
    else return alig;
}

// This method returns the consistency value vector for a given alignment
// against a set of alignments with the same sequences
bool statisticsConsistency::forceComparison(newAlignment **vectAlignments,
                                            int numAlignments,
                                            newAlignment *selected,
                                            float *columnsValue) {
	 // Create a timerLevel that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("bool statisticsConsistency::forceComparison("
                        "newAlignment **vectAlignments, "
                        "int numAlignments, "
                        "newAlignment *selected, "
                        "float *columnsValue) ");

    int *correspNames, *columnSeqMatrix, *columnSeqMatrixAux;
    int i, j, k, ll, numResidues, numSeqs, pairRes, hit;
    bool appearErrors = false;
    string *names;

    // Get some parameters from the alignment that has
    // been selected
    numResidues = selected->getNumAminos();
    numSeqs = selected->getNumSpecies();

    // Initialize the vector where we are going to store
    // the proportion of hits for each column in the
    // selected alignment
    utils::initlVect(columnsValue, numResidues, 0);

    // Allocate dinamic local memory
    names = new string[numSeqs];
    correspNames = new int[numSeqs];
    columnSeqMatrix = new int[numSeqs];
    columnSeqMatrixAux = new int[numSeqs];
    // Check that all of alignment has the same number of
    // sequence as well as there exists a correspondence
    // between the names for each pars of aligs.
    for (i = 0; i < numAlignments; i++) {

        if (numSeqs != vectAlignments[i]->getNumSpecies()) {
            debug.report(ErrorCode::DifferentNumberOfSequencesInCompareset);
            appearErrors = true;
            break;
        }

        vectAlignments[i]->getSequences(names);
        if (!selected->getSequenceNameOrder(names, correspNames)) {
            debug.report(ErrorCode::DifferentSeqsNamesInCompareset);
            appearErrors = true;
            break;
        }
    }
    // Changes the order in sequences number matrix
    // according to the order in the selected alignment
    for (i = 0; i < numAlignments; i++) {
        vectAlignments[i]->getSequences(names);
        selected->getSequenceNameOrder(names, correspNames);
        vectAlignments[i]->SequencesMatrix->setOrder(correspNames);
    }
    // Do the same analysis for each column
    for (i = 0, pairRes = 0, hit = 0;
         i < numResidues && !appearErrors;
         i++, pairRes = 0, hit = 0) {

        // We get back the sequence position for each residue
        // from every column in the selected alignment
        utils::initlVect(columnSeqMatrix, numSeqs, 0);
        selected->SequencesMatrix->getColumn(i, columnSeqMatrix);
        // For each residue pairs, we look for it in the rest
        // of alignments
        for (j = 0; j < numSeqs; j++) {

            // If there is a residue, not a gap, we carry out the
            // assesment
            if (columnSeqMatrix[j] != 0) {
                for (k = 0; k < numAlignments; k++) {

                    // We look for the same residue in the same row in
                    // the rest of alignments
                    utils::initlVect(columnSeqMatrixAux, numSeqs, 0);
                    vectAlignments[k]->SequencesMatrix->getColumn(columnSeqMatrix[j], j, columnSeqMatrixAux);
                    // We count when we get the same residue pairs in the
                    // rest of alignments
                    for (ll = j + 1; ll < numSeqs; ll++)
                        if (columnSeqMatrix[ll] != 0) {
                            if (columnSeqMatrix[ll] == columnSeqMatrixAux[ll])
                                hit++;
                            pairRes++;
                        }
                }
            }
        }
        // Store the hits proportion for each column
        if (pairRes != 0) columnsValue[i] += ((1.0 * hit) / pairRes);
    }

    // Deallocate dynamic memory
    delete[] names;
    delete[] correspNames;
    delete[] columnSeqMatrix;
    delete[] columnSeqMatrixAux;

    // Return if the process was perform without errors
    return !appearErrors;

}

// This method applies a specific windows size to a selected alignment 
bool statisticsConsistency::applyWindow(int _halfWindow) {
	 // Create a timerLevel that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("bool statisticsConsistency::applyWindow(int columns, int halfWindowApplied, float *columnsValue) ");

    if (_halfWindow > residues / 4)
    {
        debug.report(ErrorCode::ConsistencyWindowTooBig);
        exit(-1);
    }

    // If the current half window is the same as the last one, don't do anything
    if (halfWindow == _halfWindow) return true;

    // Save the requested half window. This is useful when making a copy of the
    // alignment, as the window values are not valid anymore but don't want to
    // calculate them if not needed anymore
    halfWindow = _halfWindow;

    // If the half window requested is 0 or a negative number
    // we simply delete the window values.
    if (_halfWindow < 1) {
        delete[] values_windowed;
        values_windowed = nullptr;
        return true;
    }

    // Initialize the values used in the calculation
    int i, j, window;

    // Initialize the consistency window array if it's null
    if (values_windowed == nullptr)
        values_windowed = new float[residues];


    window = 2 * halfWindow + 1;

    // Do the average window calculations
    for (i = 0; i < residues; i++) {
        values_windowed[i] = 0.F;
        for (j = i - halfWindow; j <= i + halfWindow; j++) {
            if (j < 0)
                values_windowed[i] += values[-j];
            else if (j >= residues)
                values_windowed[i] += values[((2 * residues - j) - 2)];
            else
                values_windowed[i] += values[j];
        }

        // Calculate the similarity value for the i column
        values_windowed[i] = values_windowed[i] / window;
    }

    return true;
}

bool statisticsConsistency::isDefinedWindow() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool statisticsConservation::isDefinedWindow(void) ");

    return (halfWindow != -1);
}

float *statisticsConsistency::getValues() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("float *statisticsConservation::getMdkWindowedVector(void) ");

    // If a window is defined
    if (isDefinedWindow()) {
        // Check if the window has been applied
        if (values_windowed == nullptr)
            applyWindow(halfWindow);
        // Return the windowed value
        return values_windowed;
    }
        // Return the original values
    else return values;
}

// Print the consistency value for each column from the selected alignment 
void statisticsConsistency::printStatisticsFileColumns(newAlignment &_alignment,
                                                       float *compareVect) {
	 // Create a timerLevel that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("void statisticsConsistency::printStatisticsFileColumns("
                        "newAlignment &_alignment, "
                        "float *values) ");

    int size = 20;

    std::string fname = _alignment.filename.substr(6, _alignment.filename.size() - 7);

    cout
            << std::setw(fname.length() + 7)
            << std::setfill(' ')
            << std::left << "" << endl;

    cout << "#\33[0;31m File :\33[0;1m" << fname << "\33[0m";

    fname = std::to_string(size);

    cout
            << std::setw(fname.length() + 7)
            << std::setfill(' ')
            << std::left << "" << endl;

    cout << "#\33[0;36m BlockSize : \33[0;1m" << fname << "\33[0m" << endl;

    fname = " Conservation per Column";

    cout << "#\33[0;32m Statistic :\33[0;1m" << fname << "\33[0m" << endl;

    cout << std::setw(_alignment.filename.substr(6, _alignment.filename.size() - 7).length() + 7)
         << std::setfill('-')
         << std::left << ""
         << std::setfill(' ')
         << std::fixed
         << endl;

    cout.precision(10);

    cout << "\33[0;33;1m"
         << std::setw(size) << std::left << " Residue" << std::setw(size) << std::left << " Consistency " << endl
         << std::setw(size) << std::left << " Number" << std::setw(size) << std::left << " Value " << endl
         << std::setfill('-')
         << "\33[0m"
         << std::setw(size) << std::right << "  " << std::setw(size) << std::right << "  " << endl
         << std::setfill(' ');

    // Print the consistency values for each column from
    // the selected alignment
    for (int i = 0; i < _alignment.residNumber; i++)
        cout << setw(size) << std::left << i + 1
             << setw(size) << std::left << compareVect[i]
             << endl;

}

// Print the consistency values accumulative distribution for the selected
// alignment
void statisticsConsistency::printStatisticsFileAcl(newAlignment &_alignment,
                                                   float *compareVect) {
	 // Create a timerLevel that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("void statisticsConsistency::printStatisticsFileAcl("
                        "newAlignment &_alignment, "
                        "float *values) ");

    int size = 20;
    float refer, *vectAux;
    int i, num;

    std::string fname = _alignment.filename.substr(6, _alignment.filename.size() - 7);

    cout
            << std::setw(fname.length() + 7)
            << std::setfill(' ')
            << std::left << "" << endl;

    cout << "#\33[0;31m File :\33[0;1m" << fname << "\33[0m";

    fname = std::to_string(size);

    cout
            << std::setw(fname.length() + 7)
            << std::setfill(' ')
            << std::left << "" << endl;

    cout << "#\33[0;36m BlockSize : \33[0;1m" << fname << "\33[0m" << endl;

    fname = " Conservation Total";

    cout << "#\33[0;32m Statistic :\33[0;1m" << fname << "\33[0m" << endl;

    cout << std::setw(_alignment.filename.substr(6, _alignment.filename.size() - 7).length() + 7)
         << std::setfill('-')
         << std::left << ""
         << std::setfill(' ')
         << std::fixed
         << endl;

    cout.precision(10);


    // Allocate dinamic memory to copy the input vector
    // and sort it
    vectAux = new float[_alignment.residNumber];
    utils::copyVect(compareVect, vectAux, _alignment.residNumber);
    utils::quicksort(vectAux, 0, _alignment.residNumber - 1);

    // Set the output precision and print the header
    std::stringstream firstLine;
    std::stringstream secondLine;
    std::stringstream thirdLine;

    firstLine << std::setw(size) << std::left << "";
    secondLine << std::setw(size) << std::left << " Number of";
    thirdLine << std::setw(size) << std::left << " residues";

    firstLine << std::setw(size) << std::left << "";
    secondLine << std::setw(size) << std::left << " Percentage";
    thirdLine << std::setw(size) << std::left << " of alignment";

    firstLine << std::setw(size) << std::left << " Accumulative";
    secondLine << std::setw(size) << std::left << " Number of";
    thirdLine << std::setw(size) << std::left << " Residues";

    firstLine << std::setw(size) << std::left << " Accumulative";
    secondLine << std::setw(size) << std::left << " Percentage of";
    thirdLine << std::setw(size) << std::left << " Residues";

    firstLine << std::setw(size) << std::left << "";
    secondLine << std::setw(size) << std::left << " Consistency";
    thirdLine << std::setw(size) << std::left << " Value";

    cout << "\33[0;33;1m"
         << firstLine.rdbuf() << endl
         << secondLine.rdbuf() << endl
         << thirdLine.rdbuf() << endl
         << "\33[0;m"
         << std::setfill('-');

    for (i = 0; i < 5; i++) {
        cout << std::setw(size) << std::right << "   ";
    }

    cout << std::endl
         << std::fixed
         << std::setfill(' ');

    cout.precision(10);
    // Fix the initial values to count how many columns
    // has the same consistency value
    refer = vectAux[0];
    num = 1;
    // Print the accumulative distribution
    for (i = 1; i < _alignment.residNumber; i++) {

        // When the method detects a new consistency value
        // print the previous value as well as its frequency
        // and starts to count how many columns are for this
        // new value

        if (refer != vectAux[i]) {
            // 
            cout << std::setw(size) << std::left << num;

            // 
            cout << std::setw(size) << std::left
                 << std::setw(size - 6) << std::right << ((float) num / _alignment.residNumber * 100.0F)
                 << std::setw(6) << std::left << " ";

            // 
            cout << std::setw(size) << std::left << i;

            // 
            cout << std::setw(size) << std::left
                 << std::setw(size - 6) << std::right << ((float) i / _alignment.residNumber * 100.0F)
                 << std::setw(6) << std::left << " ";

            // 
            cout << std::setw(size) << std::left << refer;

            // End line
            cout << endl;

            refer = vectAux[i];
            num = 1;
        } else num++;
    }

    // Print the last consistency value as well as its
    // frequency

    cout << std::setw(size) << std::left << num;

    cout << std::setw(size) << std::left
         << std::setw(size - 6) << std::right << ((float) num / _alignment.residNumber * 100.0F)
         << std::setw(6) << std::left << " ";

    cout << std::setw(size) << std::left << i;

    cout << std::setw(size) << std::left
         << std::setw(size - 6) << std::right << ((float) i / _alignment.residNumber * 100.0F)
         << std::setw(6) << std::left << " ";

    cout << std::setw(size) << std::left << refer;

    cout << endl;

    // Deallocate dynamic memory
    delete[] vectAux;
}

void statisticsConsistency::delete_variables() {

    delete [] compareAlignmentsArray;
    for (i = 0; i < numFiles; i++) {
        delete [] filesToCompare[i];
    }
    delete [] filesToCompare;

    delete [] line;
}

statisticsConsistency::~statisticsConsistency() {
    if (--(*refCounter) == 0)
    {
        delete [] values;
        delete [] values_windowed;
    }
    _alignment = nullptr;
}

statisticsConsistency::statisticsConsistency(newAlignment *pAlignment,
                                             statisticsConsistency *pConsistency) {
    _alignment      = pAlignment;
    values          = pConsistency->values;
    values_windowed = pConsistency->values_windowed;

    refCounter      = pConsistency->refCounter;
    (*refCounter)++;
}

statisticsConsistency::statisticsConsistency() {
    refCounter = new int(1);
}

