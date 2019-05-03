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

#include "Alignment/sequencesMatrix.h"
#include "Statistics/Consistency.h"
#include "Alignment/Alignment.h"
#include "InternalBenchmarker.h"
#include "Statistics/Manager.h"
#include "trimalManager.h"
#include "reportsystem.h"
#include "defines.h"
#include "utils.h"

namespace statistics {
    
#define LONG 80


    bool Consistency::perform(char *comparesetFilePath,
                              FormatHandling::FormatManager &formatManager,
                              trimAlManager &manager,
                              char *forceFile) {

        std::unique_ptr<char[]> line(new char[1024]);

        // Open the file that contains the paths of the files.
        std::ifstream compare;
        compare.open(comparesetFilePath, std::ifstream::in);
        while (compare.getline(line.get(), 1024)) numFiles++;
        compare.close();
        compare.open(comparesetFilePath);

        compareAlignmentsArray = new Alignment *[numFiles];
        char **filesToCompare = new char *[numFiles];
        for (i = 0; i < numFiles; i++) filesToCompare[i] = nullptr;


        int prevType = SequenceTypes::NotDefined;

        char c;

        string nline;

        // Load all the alignments to compare
        // Check if they: are aligned and the type is the same
        // Store the maximum number of amino acids present on all alignments
        for (i = 0; i < numFiles; i++) {
            // Search for end of line.
            for (nline.clear(), compare.read(&c, 1);
                 (c != '\n') && ((!compare.eof()));
                 compare.read(&c, 1)) {
                nline += c;
            }

            // Save the alignment path
            filesToCompare[i] = new char[nline.size() + 1];
            strcpy(filesToCompare[i], nline.c_str());

            // Load the alignment
            compareAlignmentsArray[i]
                    = formatManager.loadAlignment(filesToCompare[i]);

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
                            = new Alignment::sequencesMatrix(compareAlignmentsArray[i]);

                    // Store maximum number of aminoacids
                    if (compareAlignmentsArray[i]->getNumAminos() > maxResidues)
                        maxResidues = compareAlignmentsArray[i]->getNumAminos();

                    // Check if alignment type is the same as last one.
                    if (prevType == SequenceTypes::NotDefined)
                        prevType = compareAlignmentsArray[i]->getAlignmentType();
                    else if (compareAlignmentsArray[i]->getAlignmentType() != prevType) {
                        debug.report(ErrorCode::AlignmentTypesNotMatching);
                        appearErrors = true;
                    }
                }
            }
        }

        // If the analysis couldn't be performed, stop the program.
        if (appearErrors) {
            debug.report(ErrorCode::ComparesetFailedAlignmentMissing);
            for (int x = 0; x < numFiles; x++)
            {
                delete [] filesToCompare[x];
                delete compareAlignmentsArray[x];
            }
            delete[] compareAlignmentsArray;
            delete[] filesToCompare;
            delete[] values;
            return appearErrors;
        } else {
            int referFile = -1;
            // If no alignment is forced to be selected, select one of them
            if (forceFile == nullptr) {
                values = new float[maxResidues];
                // Perform stat calculation and
                //  choice the best scoring alignment
                referFile = Consistency::compareAndChoose(
                        compareAlignmentsArray,
                        filesToCompare,
                        values,
                        numFiles,
                        // Verbosity
                        manager.stats >= 0 && manager.outfile != nullptr);

                // If no alignment could be selected, stop the program
                if (referFile == -1) {
                    for (int x = 0; x < numFiles; x++)
                    {
                        delete [] filesToCompare[x];
                        delete compareAlignmentsArray[x];
                    }
                    delete[] filesToCompare;
                    delete[] compareAlignmentsArray;
                    delete[] values;
                    manager.appearErrors = true;
                    return true;
                }

                // Specify the selected alignment as origAlig
                //  (as if it was fed thought -in argument)
                manager.origAlig
                        = new Alignment(*compareAlignmentsArray[referFile]);
            }

                // If there is an alignment that has been forcibly selected
            else {
                values = new float[manager.origAlig->getNumAminos()];
                appearErrors = !Consistency::forceComparison(
                        compareAlignmentsArray,
                        numFiles,
                        manager.origAlig,
                        values);

                // If forcing comparison failed, exit the program
                if (appearErrors) {
                    for (int x = 0; x < numFiles; x++)
                    {
                        delete [] filesToCompare[x];
                        delete compareAlignmentsArray[x];
                    }
                    delete[] filesToCompare;
                    delete[] compareAlignmentsArray;
                    delete[] values;
                    manager.appearErrors = true;
                    return true;
                }
            }

            // Store the cross reference between alignment
            //  and this, it's consistency stat
            manager.origAlig->Statistics->consistency = this;
            manager.CS = nullptr;
            alig = manager.origAlig;
            residues = alig->originalNumberOfResidues;


            // Apply window sizes
            if (manager.windowSize != -1)
                appearErrors += !applyWindow(manager.windowSize);
            else if (manager.consistencyWindow != -1)
                appearErrors += !applyWindow(manager.consistencyWindow);

            if (appearErrors) {
                for (int x = 0; x < numFiles; x++)
                {
                    delete [] filesToCompare[x];
                    delete compareAlignmentsArray[x];
                }
                delete[] filesToCompare;
                delete[] compareAlignmentsArray;
                delete[] values;
                manager.appearErrors = true;
                return appearErrors;
            }

            // If no output format is provided
            //  we'll use the selected alignment format
            if (manager.oformats.empty()) {
                manager.oformats.emplace_back(
                        formatManager.getFileFormatName(
                                forceFile == nullptr ? filesToCompare[referFile] : forceFile)
                );
            }
            for (int x = 0; x < numFiles; x++)
            {
                delete [] filesToCompare[x];
                delete compareAlignmentsArray[x];
            }
            delete[] filesToCompare;
            delete[] compareAlignmentsArray;
        }
        return appearErrors;
    }


// This method compares a set of alignment in order to select the most
// consistent one respect of the other ones. To compute the consistency
// values we use the proportion of residue pairs per column in the aligs
// to compare 

    int Consistency::compareAndChoose(Alignment **vectAlignments,
                                                char **fileNames,
                                                float *columnsValue,
                                                int numAlignments,
                                                bool verbosity) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("int Consistency::compareAndChoose("
                    "Alignment **vectAlignments, "
                    "char **fileNames, "
                    "float *columnsValue, "
                    "int numAlignments, "
                    "bool verbosity) ");

        int *numResiduesAlig, *correspNames, *columnSeqMatrix, *columnSeqMatrixAux;
        int i, k, l, j, m, numSeqs, pairRes, hits, alignmentIndex = 0;
        float max = 0, value = 0, **vectHits;
        bool appearErrors = false;
        string *names;


        // Get some parameters from the alignment that has
        // been selected */
        numSeqs = vectAlignments[0]->getNumSpecies();

        // Allocate dinamic local memory
        names = new string[numSeqs];
        correspNames = new int[numSeqs];
        numResiduesAlig = new int[numAlignments];
        columnSeqMatrix = new int[numSeqs];
        vectHits = new float *[numAlignments];
        columnSeqMatrixAux = new int[numSeqs];

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

        if (!appearErrors) {

            max = 0;

            // Changes the order in sequences number matrix
            // according to the order in the selected alignment
            for (i = 1; i < numAlignments; i++) {
                vectAlignments[i]->getSequences(names);
                vectAlignments[0]->getSequenceNameOrder(names, correspNames);
                vectAlignments[i]->SequencesMatrix->setOrder(correspNames);
            }

            // Get back the residues number for each alignment
            for (i = 0; i < numAlignments; i++)
                numResiduesAlig[i] = vectAlignments[i]->getNumAminos();

            // Start the comparison among the alignments
            for (i = 0; i < numAlignments; i++) {
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
                    alignmentIndex = i;
                    max = value / numResiduesAlig[i];
                }
            }

            // Prints the alignment that have been selected
            if (verbosity) {
                cout << "\t\t\t\t\t--------------" << endl;
                cout << endl << "File Selected:\t" << fileNames[alignmentIndex] << endl << "Value:\t\t" << max << endl << endl;
            }

            // The method returns a vector with the consistency
            // value for each column in the selected alignment
            if (columnsValue != nullptr) {
                utils::initlVect(columnsValue, numResiduesAlig[alignmentIndex], -1);
                for (i = 0; i < numResiduesAlig[alignmentIndex]; i++)
                    columnsValue[i] = vectHits[alignmentIndex][i];
            }
        }

        // Deallocate memory
        for (i = 0; i < numAlignments; i++)
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
        else return alignmentIndex;
    }

// This method returns the consistency value vector for a given alignment
// against a set of alignments with the same sequences
    bool Consistency::forceComparison(Alignment **vectAlignments,
                                                int numAlignments,
                                                Alignment *selected,
                                                float *columnsValue) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Consistency::forceComparison("
                    "Alignment **vectAlignments, "
                    "int numAlignments, "
                    "Alignment *selected, "
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

        // Allocate dynamic local memory
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
        
        if (selected->SequencesMatrix == nullptr) 
            selected->SequencesMatrix = new Alignment::sequencesMatrix(selected);
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
    bool Consistency::applyWindow(int halfW) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Consistency::applyWindow(int columns, int halfWindowApplied, float *columnsValue) ");

        if (halfW > residues / 4) {
            debug.report(ErrorCode::ConsistencyWindowTooBig);
            return false;
        }

        // If the current half window is the same as the last one, don't do anything
        if (halfWindow == halfW) return true;

        // Save the requested half window. This is useful when making a copy of the
        // alignment, as the window values are not valid anymore but don't want to
        // calculate them if not needed anymore
        halfWindow = halfW;

        // If the half window requested is 0 or a negative number
        // we simply delete the window values.
        if (halfW < 1) {
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

    bool Consistency::isWindowDefined() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Similarity::isWindowDefined(void) ");

        return (halfWindow != -1);
    }

    float *Consistency::getValues() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("float *Similarity::getMdkWindowedVector(void) ");

        // If a window is defined
        if (isWindowDefined()) {
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
    void Consistency::printStatisticsFileColumns(Alignment &alig,
                                                           float *compareVect) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Consistency::printStatisticsFileColumns("
                    "Alignment &alig, "
                    "float *values) ");

        int size = 20;

        std::string fname = alig.filename;

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

        fname = " Similarity per Column";

        cout << "#\33[0;32m Statistic :\33[0;1m" << fname << "\33[0m" << endl;

        cout << std::setw(alig.filename.size())
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
        for (int i = 0; i < alig.numberOfResidues; i++)
            cout << setw(size) << std::left << i + 1
                 << setw(size) << std::left << compareVect[i]
                 << endl;

    }

    // Print the consistency values accumulative distribution for the selected
    // alignment
    void Consistency::printStatisticsFileAcl(
            Alignment &alig, float *compareVect) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Consistency::printStatisticsFileAcl("
                    "Alignment &alig, "
                    "float *compareVect) ");

        int size = 20;
        float refer, *vectAux;
        int i, num;

        std::string fname = alig.filename;

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

        fname = " Similarity Total";

        cout << "#\33[0;32m Statistic :\33[0;1m" << fname << "\33[0m" << endl;

        cout << std::setw(alig.filename.size())
             << std::setfill('-')
             << std::left << ""
             << std::setfill(' ')
             << std::fixed
             << endl;

        cout.precision(10);


        // Allocate dinamic memory to copy the input vector
        // and sort it
        vectAux = new float[alig.numberOfResidues];
        utils::copyVect(compareVect, vectAux, alig.numberOfResidues);
        utils::quicksort(vectAux, 0, alig.numberOfResidues - 1);

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
        for (i = 1; i < alig.numberOfResidues; i++) {

            // When the method detects a new consistency value
            // print the previous value as well as its frequency
            // and starts to count how many columns are for this
            // new value

            if (refer != vectAux[i]) {
                //
                cout << std::setw(size) << std::left << num;

                //
                cout << std::setw(size) << std::left
                     << std::setw(size - 6) << std::right << ((float) num / alig.numberOfResidues * 100.0F)
                     << std::setw(6) << std::left << " ";

                //
                cout << std::setw(size) << std::left << i;

                //
                cout << std::setw(size) << std::left
                     << std::setw(size - 6) << std::right << ((float) i / alig.numberOfResidues * 100.0F)
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
             << std::setw(size - 6) << std::right << ((float) num / alig.numberOfResidues * 100.0F)
             << std::setw(6) << std::left << " ";

        cout << std::setw(size) << std::left << i;

        cout << std::setw(size) << std::left
             << std::setw(size - 6) << std::right << ((float) i / alig.numberOfResidues * 100.0F)
             << std::setw(6) << std::left << " ";

        cout << std::setw(size) << std::left << refer;

        cout << endl;

        // Deallocate dynamic memory
        delete[] vectAux;
    }

    Consistency::~Consistency() {
        if (--(*refCounter) == 0) {
            delete[] values;
            delete[] values_windowed;
        }
        alig = nullptr;
        delete refCounter;
    }

    Consistency::Consistency(Alignment *pAlignment,
                                                 Consistency *pConsistency) {
        alig = pAlignment;
        values = pConsistency->values;
        values_windowed = pConsistency->values_windowed;

        refCounter = pConsistency->refCounter;
        (*refCounter)++;
    }

    Consistency::Consistency() {
        refCounter = new int(1);
    }

}