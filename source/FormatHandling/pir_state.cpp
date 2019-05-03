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

#include "FormatHandling/pir_state.h"

#include "FormatHandling/FormatManager.h"
#include "defines.h"
#include "utils.h"

namespace FormatHandling {
int pir_state::CheckAlignment(std::istream *origin) {
    char *line;
    origin->seekg(0);
    line = utils::readLine(*origin);
    if (line == nullptr) return 0;
    if (strlen(line) > 4) {
        if (line[0] == '>')
            if (line[3] == ';')
            {
                delete [] line;
                return 2;
            }
    }
    delete[] line;
    return 0;
}

Alignment *pir_state::LoadAlignment(const std::string &filename) {
    /* NBRF/PIR file format parser */

    Alignment *alig = new Alignment();

    bool seqIdLine, seqLines;
    char *str, *line = nullptr;
    std::ifstream file;
    int i;

    /* Check the file and its content */
    file.open(filename, std::ifstream::in);
    if (!utils::checkFile(file))
        return nullptr;

    /* Store input file name for posterior uses in other formats */
    // alig->filename.append("!Title ");
    alig->filename.append(filename);
    alig->filename.append(";");

    /* Compute how many sequences are in the input alignment */
    alig->numberOfSequences = 0;
    while (!file.eof()) {

        /* Deallocate previously used dinamic memory */
        if (line != nullptr)
            delete[] line;

        /* Read lines in a safe way */
        line = utils::readLine(file);
        if (line == nullptr)
            continue;

        /* It the line starts by ">" means that a new sequence has been found */
        str = strtok(line, DELIMITERS);
        if (str == nullptr)
            continue;

        /* If a sequence name flag is detected, increase sequences counter */
        if (str[0] == '>')
            alig->numberOfSequences++;
    }

    /* Finish to preprocess the input file. */
    file.clear();
    file.seekg(0);

    /* Allocate memory for the input alignmet */
    alig->sequences = new std::string[alig->numberOfSequences];
    alig->seqsName = new std::string[alig->numberOfSequences];
    alig->seqsInfo = new std::string[alig->numberOfSequences];

    /* Initialize some local variables */
    seqIdLine = true;
    seqLines = false;
    i = -1;

    /* Read the entire input file */
    while (!file.eof()) {

        /* Deallocate local memory */
        delete[] line;

        /* Read lines in a safe way */
        line = utils::readLine(file);
        if (line == nullptr)
            continue;

        /* Sequence ID line.
         * Identification of these kind of lines is based on presence of ">" and ";"
         * symbols at positions 0 and 3 respectively */
        if ((line[0] == '>') && (line[3] == ';') && (seqIdLine)) {
            seqIdLine = false;
            i += 1;

            /* Skip information about sequence datatype - Better to check */
            strtok(line, ">;");

            /* and the sequence identifier itself */
            str = strtok(nullptr, ">;");
            alig->seqsName[i].append(str, strlen(str));
        }

            /* Line just after sequence Id line contains a textual description of
             * the sequence. */
        else if ((!seqIdLine) && (!seqLines)) {
            seqLines = true;
            alig->seqsInfo[i].append(line, strlen(line));
        }

            /* Sequence lines itself */
        else if (seqLines) {

            /* Check whether a sequence end symbol '*' exists in current line.
             * In that case, set appropriate flags to read a new sequence */
            if (line[strlen(line) - 1] == '*') {
                seqLines = false;
                seqIdLine = true;
            }

            /* Process line */
            str = strtok(line, OTHDELIMITERS);
            while (str != nullptr) {
                if (str[strlen(str) - 1] == '*')
                    alig->sequences[i].append(str, strlen(str) - 1);
                else
                    alig->sequences[i].append(str, strlen(str));
                str = strtok(nullptr, OTHDELIMITERS);
            }
        }
    }
    /* Close the input file */
    file.close();

    /* Deallocate dinamic memory */
    delete[] line;

    /* Check the matrix's content */
    alig->fillMatrices(true);
    alig->originalNumberOfSequences = alig->numberOfSequences;
    alig->originalNumberOfResidues = alig->numberOfResidues;


    return alig;
}

bool pir_state::SaveAlignment(const Alignment &alignment, std::ostream *output) {

    /* Generate output alignment in NBRF/PIR format. Sequences can be unaligned */

    int i, j, k, l;
    std::string alg_datatype, *tmpMatrix;

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    if (Machine->reverse)
    {
        /* Allocate local memory for generating output alignment */
        tmpMatrix = new std::string[alignment.originalNumberOfSequences];
        for(i = 0; i < alignment.originalNumberOfSequences; i++)
            tmpMatrix[i] = utils::getReverse(alignment.sequences[i]);
    }
    else tmpMatrix = alignment.sequences;

    /* Compute output file datatype */
    alignment.getAlignmentType();
    if (alignment.getAlignmentType() & SequenceTypes::DNA)
        alg_datatype = "DL";
    else if (alignment.getAlignmentType() & SequenceTypes::RNA)
        alg_datatype = "RL";
    else if (alignment.getAlignmentType() & SequenceTypes::AA)
        alg_datatype = "P1";




    /* Print alignment */
    for (i = 0; i < alignment.originalNumberOfSequences; i++) {
        if (alignment.saveSequences && alignment.saveSequences[i] == -1) continue;

        /* Print sequence datatype and its name */
        if ((alignment.seqsInfo != nullptr) /*&& (iformat == oformat)*/)
            (*output) << ">" << alg_datatype << ";" << alignment.seqsName[i]
                      << "\n" << alignment.seqsInfo[i] << "\n";
        else
            (*output) << ">" << alg_datatype << ";" << alignment.seqsName[i] << "\n"
                      << alignment.seqsName[i] << " " << alignment.sequences[i].length() << " bases\n";

        for (j = 0, k = 0, l = 0; j < alignment.sequences[i].length(); j++) {
            if (alignment.saveResidues != nullptr && alignment.saveResidues[j] == -1) {
//                 if (j == alignment->sequences[i].length() -1 ) 
//                     (*output) << "\n";
            } else {
                if (k % 10 == 0) (*output) << " ";
                (*output) << tmpMatrix[i][j];
                k++;
                if (j == alignment.sequences[i].length() - 1);
                else if (k % 50 == 0) (*output) << "\n";
            }
        }
        if (k % 50 == 0) (*output) << "\n";
        if (k % 10 == 0) (*output) << " ";
        (*output) << "*\n\n";
    }

    /* Deallocate local memory */
    if (Machine->reverse)
        delete [] tmpMatrix;

    return true;
}

bool pir_state::RecognizeOutputFormat(const std::string &FormatName) {
    if (BaseFormatHandler::RecognizeOutputFormat(FormatName))
        return true;
    return FormatName == "pir" || FormatName == "nbrf" ||
           FormatName == "PIR" || FormatName == "NBRF";
}

}