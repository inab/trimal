#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/ReadWriteMS/phylip40_state.h"
#include "../../include/defines.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include "../../include/newAlignment.h"

using namespace std;

int Phylip40State::CheckAlignment(istream* origin)
{
        origin->seekg(0);
    origin->clear();
    char *firstWord = NULL, *line = NULL;
    int blocks = 0;
    string nline;

    /* Read first valid line in a safer way */
    do {
        line = utils::readLine(*origin);
    } while ((line == NULL) && (!origin->eof()));

    /* If the file end is reached without a valid line, warn about it */
    if (origin->eof())
        return false;

    /* Otherwise, split line */
    firstWord = strtok(line, OTHDELIMITERS);

    /* Phylip Format */
    {
        /* Determine specific phylip format: sequential or interleaved. */

        /* Get number of sequences and residues */
        int sequenNumber = atoi(firstWord);
        int residNumber = 0;
        firstWord = strtok(NULL, DELIMITERS);
        if(firstWord != NULL)
            residNumber = atoi(firstWord);
        else return 0;

        /* If there is only one sequence, use by default sequential format since
         * it is impossible to determine exactly which phylip format is */
        if((sequenNumber == 1) && (residNumber != 0))
            return 1;

            /* If there are more than one sequence, analyze sequences distribution to
             * determine its format. */
        else if((sequenNumber != 0) && (residNumber != 0)) {
            blocks = 0;

            /* Read line in a safer way */
            do {
                if (line != NULL)
                    delete [] line;
                line = utils::readLine(*origin);
            } while ((line == NULL) && (!origin->eof()));

            /* If the file end is reached without a valid line, warn about it */
            if (origin->eof())
                return 0;

            firstWord = strtok(line, DELIMITERS);
            while(firstWord != NULL) {
                blocks++;
                firstWord = strtok(NULL, DELIMITERS);
            }

            /* Read line in a safer way */
            do {
                if (line != NULL)
                    delete [] line;
                line = utils::readLine(*origin);
            } while ((line == NULL) && (!origin->eof()));

            firstWord = strtok(line, DELIMITERS);
            while(firstWord != NULL) {
                blocks--;
                firstWord = strtok(NULL, DELIMITERS);
            }

            /* If the file end is reached without a valid line, warn about it */
            if (origin->eof())
                return false;

            /* Phylip Interleaved (12) or Sequential (11) */
            return (!blocks) ? 1 : 0;
        }
    }
    return 0;
}

newAlignment* Phylip40State::LoadAlignment(std::string filename)
{
    /* PHYLIP/PHYLIP 4 (Sequential) file format parser */
    newAlignment * _alignment = new newAlignment();
    char *str, *line = NULL;
    ifstream file;
    int i;

    /* Check the file and its content */
    file.open(filename, ifstream::in);
    if(!utils::checkFile(file))
        return nullptr;

    /* Store some data about filename for possible uses in other formats */
    filename.append("!Title ");
    filename.append(filename);
    filename.append(";");

    /* Read first valid line in a safer way */
    do {
        line = utils::readLine(file);
    } while ((line == NULL) && (!file.eof()));

    /* If the file end is reached without a valid line, warn about it */
    if (file.eof())
        return nullptr;

    /* Read the input sequences and residues for each sequence numbers */
    str = strtok(line, DELIMITERS);
    _alignment->sequenNumber = 0;
    if(str != NULL)
        _alignment->sequenNumber = atoi(str);

    str = strtok(NULL, DELIMITERS);
    _alignment->residNumber = 0;
    if(str != NULL)
        _alignment->residNumber = atoi(str);

    /* If something is wrong about the sequences or/and residues number,
     * return an error to warn about that */
    if((_alignment->sequenNumber == 0) || (_alignment->residNumber == 0))
        return nullptr;

    /* Allocate memory  for the input data */
    _alignment->sequences  = new string[_alignment->sequenNumber];
    _alignment->seqsName   = new string[_alignment->sequenNumber];

    /* Read the lines block containing the sequences name + first fragment */
    i = 0;
    while((i < _alignment->sequenNumber) && (!file.eof())){

        /* Read lines in a safer way. Destroy previous stored information */
        if (line != NULL)
            delete [] line;
        line = utils::readLine(file);

        /* It the input line/s are blank lines, skip the loop iteration  */
        if(line == NULL)
            continue;

        /* First token: Sequence name */
        str = strtok(line, DELIMITERS);
        _alignment->seqsName[i].append(str, strlen(str));

        /* Trim the rest of the line from blank spaces, tabs, etc and store it */
        str = strtok(NULL, DELIMITERS);
        while(str != NULL) {
            _alignment->sequences[i].append(str, strlen(str));
            str = strtok(NULL, DELIMITERS);
        }
        i++;
    }

    /* Read the rest of the input file */
    while(!file.eof()) {

        /* Try to get for each sequences its corresponding residues */
        i = 0;
        while((i < _alignment->sequenNumber) && (!file.eof())) {
            /* Read lines in a safer way. Destroy previous stored information */
            if (line != NULL)
                delete [] line;

            line = utils::readLine(file);
            /* It the input line/s are blank lines, skip the loop iteration  */
            if(line == NULL)
                continue;

            /* Remove from the current line non-printable characters and add fragments
             * to previous stored sequence */
            str = strtok(line, DELIMITERS);
            while(str != NULL) {
                _alignment->sequences[i].append(str, strlen(str));
                str = strtok(NULL, DELIMITERS);
            }
            i++;
        }
    }

    /* Close the input file and delete dinamic memory */
    file.close();
    if (line != NULL)
        delete [] line;

    /* Check the matrix's content */
    _alignment->fillMatrices(true);
    return _alignment;
}

void Phylip40State::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
   /* Generate output alignment in PHYLIP/PHYLIP 4 format (sequential) */

    int i, j, maxLongName;
    string *tmpMatrix;

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!alignment->isAligned) {
        cerr << endl << "ERROR: Sequences are not aligned. Format (phylip) "
             << "not compatible with unaligned sequences." << endl << endl;
        return ;
    }

    /* Allocate local memory for generating output alignment */
    tmpMatrix = new string[alignment->sequenNumber];

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    for(i = 0; i < alignment->sequenNumber; i++)
        tmpMatrix[i] = (!alignment->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    /* Depending on if short name flag is activated (limits sequence name up to
     * 10 characters) or not, get maximum sequence name length */
    maxLongName = PHYLIPDISTANCE;
    for(i = 0; (i < alignment->sequenNumber) && (!Machine->shortNames); i++)
        maxLongName = utils::max(maxLongName, alignment->seqsName[i].size());

    /* Generating output alignment */
    /* First Line: Sequences Number & Residued Number */
    (*output) << " " << alignment->sequenNumber << " " << alignment->residNumber << endl;

    /* First Block: Sequences Names & First 60 residues */
    for(i = 0; i < alignment->sequenNumber; i++)
        (*output) << setw(maxLongName + 3) << left << alignment->seqsName[i].substr(0, maxLongName)
             << tmpMatrix[i].substr(0, 60) << endl;
    (*output) << endl;

    /* Rest of blocks: Print 60 residues per each blocks of sequences */
    for(i = 60; i < alignment->residNumber; i += 60) {
        for(j = 0; j < alignment->sequenNumber; j++)
            (*output) << tmpMatrix[j].substr(i, 60) << endl;
        (*output) << endl;
    }
    (*output) << endl;

    /* Deallocate local memory */
    delete [] tmpMatrix;
}

bool Phylip40State::RecognizeOutputFormat(std::string FormatName)
{
    if (FormatName == "phylip" || FormatName == "phylip40") return true;
    return false;
}

