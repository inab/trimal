#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/ReadWriteMS/pir_state.h"
#include "../../include/defines.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

int PirState::CheckAlignment(istream* origin)
{
    char * line;
    origin->seekg(0);
    line = utils::readLine(*origin);
    if (strlen(line) > 4)
    {
        if (line[0] == '>')
            if (line[3] == ';')
                return 2;
    }
    return 0;
}

newAlignment* PirState::LoadAlignment(std::__cxx11::string filename)
{
    /* NBRF/PIR file format parser */
    
    newAlignment* _alignment = new newAlignment();

    bool seqIdLine, seqLines;
    char *str, *line = NULL;
    ifstream file;
    int i;

    /* Check the file and its content */
    file.open(filename, ifstream::in);
    if(!utils::checkFile(file))
        return nullptr;

    /* Store input file name for posterior uses in other formats */
    filename.append("!Title ");
    filename.append(filename);
    filename.append(";");

    /* Compute how many sequences are in the input alignment */
    _alignment->sequenNumber = 0;
    while(!file.eof()) {

        /* Deallocate previously used dinamic memory */
        if (line != NULL)
            delete [] line;

        /* Read lines in a safe way */
        line = utils::readLine(file);
        if (line == NULL)
            continue;

        /* It the line starts by ">" means that a new sequence has been found */
        str = strtok(line, DELIMITERS);
        if (str == NULL)
            continue;

        /* If a sequence name flag is detected, increase sequences counter */
        if(str[0] == '>')
            _alignment->sequenNumber++;
    }

    /* Finish to preprocess the input file. */
    file.clear();
    file.seekg(0);

    /* Allocate memory for the input alignmet */
    _alignment->sequences = new string[_alignment->sequenNumber];
    _alignment->seqsName  = new string[_alignment->sequenNumber];
    _alignment->seqsInfo  = new string[_alignment->sequenNumber];

    /* Initialize some local variables */
    seqIdLine = true;
    seqLines = false;
    i = -1;

    /* Read the entire input file */
    while(!file.eof()) {

        /* Deallocate local memory */
        if (line != NULL)
            delete [] line;

        /* Read lines in a safe way */
        line = utils::readLine(file);
        if (line == NULL)
            continue;

        /* Sequence ID line.
         * Identification of these kind of lines is based on presence of ">" and ";"
         * symbols at positions 0 and 3 respectively */
        if((line[0] == '>') && (line[3] == ';') && (seqIdLine)) {
            seqIdLine = false;
            i += 1;

            /* Store information about sequence datatype */
            str = strtok(line, ">;");
            _alignment->seqsInfo[i].append(str, strlen(str));

            /* and the sequence identifier itself */
            str = strtok(NULL, ">;");
            _alignment->seqsName[i].append(str, strlen(str));
        }

            /* Line just after sequence Id line contains a textual description of
             * the sequence. */
        else if((!seqIdLine) && (!seqLines)) {
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
            while (str != NULL) {
                _alignment->sequences[i].append(str, strlen(str));
                str = strtok(NULL, OTHDELIMITERS);
            }
        }
    }
    /* Close the input file */
    file.close();

    /* Deallocate dinamic memory */
    if (line != NULL)
        delete [] line;

    /* Check the matrix's content */
    _alignment->fillMatrices(true);
    return _alignment; 
}

bool PirState::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
   /* Generate output alignment in NBRF/PIR format. Sequences can be unaligned */

    int i, j, k;
    string alg_datatype, *tmpMatrix;

    /* Allocate local memory for generating output alignment */
    tmpMatrix = new string[alignment->sequenNumber];

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    for(i = 0; i < alignment->sequenNumber; i++)
        tmpMatrix[i] = (!alignment->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    /* Compute output file datatype */
    alignment->getAlignmentType();
    if (alignment->dataType == DNAType)
        alg_datatype = "DL";
    else if (alignment->dataType == RNAType)
        alg_datatype = "RL";
    else if (alignment->dataType == AAType)
        alg_datatype = "P1";

    /* Print alignment */
    for(i = 0; i < alignment->sequenNumber; i++) {

        /* Print sequence datatype and its name */
        if((alignment->seqsInfo != NULL) /*&& (iformat == oformat)*/)
            (*output) << ">" << alignment->seqsInfo[i].substr(0, 2) << ";" << alignment->seqsName[i]
                 << endl << alignment->seqsInfo[i].substr(2) << endl;
        else
            (*output) << ">" << alg_datatype << ";" << alignment->seqsName[i] << endl
                 << alignment->seqsName[i] << " " << alignment->residuesNumber[i] << " bases" << endl;

        /* Write the sequence */
        for(j = 0; j < alignment->residNumber; j += 50) {
            for(k = j; (k < alignment->residNumber) && (k < (j + 50)); k += 10)
                (*output) << " " << tmpMatrix[i].substr(k, 10);
            if((j + 50) >= alignment->residNumber)
                (*output) << "*";
            (*output) << endl;
        }
        (*output) << endl;
    }

    /* Deallocate local memory */
    delete [] tmpMatrix;
    
    return true;
}

bool PirState::RecognizeOutputFormat(std::string FormatName)
{
    if (ReadWriteBaseState::RecognizeOutputFormat(FormatName)) return true;
    if (FormatName == "pir" || FormatName == "PIR" ||
        FormatName == "nbrf" || FormatName == "NBRF")
        return true;
    return false;
}

