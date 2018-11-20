#include "../../include/ReadWriteMS/pir_state.h"

#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/defines.h"
#include "../../include/utils.h"

int pir_state::CheckAlignment(std::istream *origin) {
    char *line;
    origin->seekg(0);
    line = utils::readLine(*origin);
    if (line == nullptr) return 0;
    if (strlen(line) > 4) {
        if (line[0] == '>')
            if (line[3] == ';')
                return 2;
    }
    delete[] line;
    return 0;
}

newAlignment *pir_state::LoadAlignment(std::string filename) {
    /* NBRF/PIR file format parser */

    newAlignment *_alignment = new newAlignment();

    bool seqIdLine, seqLines;
    char *str, *line = nullptr;
    std::ifstream file;
    int i;

    /* Check the file and its content */
    file.open(filename, std::ifstream::in);
    if (!utils::checkFile(file))
        return nullptr;

    /* Store input file name for posterior uses in other formats */
    _alignment->filename.append("!Title ");
    _alignment->filename.append(filename);
    _alignment->filename.append(";");

    /* Compute how many sequences are in the input alignment */
    _alignment->sequenNumber = 0;
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
            _alignment->sequenNumber++;
    }

    /* Finish to preprocess the input file. */
    file.clear();
    file.seekg(0);

    /* Allocate memory for the input alignmet */
    _alignment->sequences = new std::string[_alignment->sequenNumber];
    _alignment->seqsName = new std::string[_alignment->sequenNumber];
    _alignment->seqsInfo = new std::string[_alignment->sequenNumber];

    /* Initialize some local variables */
    seqIdLine = true;
    seqLines = false;
    i = -1;

    /* Read the entire input file */
    while (!file.eof()) {

        /* Deallocate local memory */
        if (line != nullptr)
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

            /* Store information about sequence datatype */
            str = strtok(line, ">;");
//             _alignment->seqsInfo[i].append(str, strlen(str));

            /* and the sequence identifier itself */
            str = strtok(nullptr, ">;");
            _alignment->seqsName[i].append(str, strlen(str));
        }

            /* Line just after sequence Id line contains a textual description of
             * the sequence. */
        else if ((!seqIdLine) && (!seqLines)) {
            seqLines = true;
            _alignment->seqsInfo[i].append(line, strlen(line));
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
                    _alignment->sequences[i].append(str, strlen(str) - 1);
                else
                    _alignment->sequences[i].append(str, strlen(str));
                str = strtok(nullptr, OTHDELIMITERS);
            }
        }
    }
    /* Close the input file */
    file.close();

    /* Deallocate dinamic memory */
    delete[] line;

    /* Check the matrix's content */
    _alignment->fillMatrices(true);
    _alignment->originalSequenNumber = _alignment->sequenNumber;
    _alignment->originalResidNumber = _alignment->residNumber;


    return _alignment;
}

bool pir_state::SaveAlignment(newAlignment *alignment, std::ostream *output, std::string *FileName) {

    /* Generate output alignment in NBRF/PIR format. Sequences can be unaligned */

    int i, j, k, l;
    std::string alg_datatype, *tmpMatrix;

    /* Allocate local memory for generating output alignment */
    tmpMatrix = new std::string[alignment->originalSequenNumber];

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    for (i = 0; i < alignment->originalSequenNumber; i++)
        tmpMatrix[i] = (!Machine->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    /* Compute output file datatype */
    alignment->getAlignmentType();
    if (alignment->getAlignmentType() & SequenceTypes::DNA)
        alg_datatype = "DL";
    else if (alignment->getAlignmentType() & SequenceTypes::RNA)
        alg_datatype = "RL";
    else if (alignment->getAlignmentType() & SequenceTypes::AA)
        alg_datatype = "P1";




    /* Print alignment */
    for (i = 0; i < alignment->originalSequenNumber; i++) {
        if (alignment->saveSequences && alignment->saveSequences[i] == -1) continue;

        /* Print sequence datatype and its name */
        if ((alignment->seqsInfo != nullptr) /*&& (iformat == oformat)*/)
            (*output) << ">" << alg_datatype << ";" << alignment->seqsName[i]
                      << "\n" << alignment->seqsInfo[i] << "\n";
        else
            (*output) << ">" << alg_datatype << ";" << alignment->seqsName[i] << "\n"
                      << alignment->seqsName[i] << " " << alignment->sequences[i].length() << " bases\n";

        for (j = 0, k = 0, l = 0; j < alignment->sequences[i].length(); j++) {
            if (alignment->saveResidues != nullptr && alignment->saveResidues[j] == -1) {
//                 if (j == alignment->sequences[i].length() -1 ) 
//                     (*output) << "\n";
            } else {
                if (k % 10 == 0) (*output) << " ";
                (*output) << tmpMatrix[i][j];
                k++;
                if (j == alignment->sequences[i].length() - 1);
                else if (k % 50 == 0) (*output) << "\n";
            }
        }
        if (k % 50 == 0) (*output) << "\n";
        if (k % 10 == 0) (*output) << " ";
        (*output) << "*\n\n";
    }

    /* Deallocate local memory */
    delete[] tmpMatrix;

    return true;
}

bool pir_state::RecognizeOutputFormat(std::string FormatName) {
    if (ReadWriteBaseState::RecognizeOutputFormat(FormatName))
        return true;
    return FormatName == "pir" || FormatName == "nbrf" ||
           FormatName == "PIR" || FormatName == "NBRF";
}

