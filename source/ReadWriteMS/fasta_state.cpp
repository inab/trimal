#include "../../include/ReadWriteMS/fasta_state.h"
#include <iostream>
#include "../../include/defines.h"
#include <stdio.h>
#include <string>
#include "../../include/ReadWriteMS/ReadWriteMachineState.h"

using namespace std;

int FastaState::CheckAlignment(istream* origin)
{
    char c;
    origin->seekg(0);
    origin->get(c);
    if (!strcmp(&c, ">"))
        return 1;
    return 0;
}

newAlignment* FastaState::LoadAlignment(std::__cxx11::string filename)
{
    /* FASTA file format parser */
    newAlignment* _alignment = new newAlignment();
    char *str, *line = NULL;
    ifstream file;
    int i;

    /* Check the file and its content */
    file.open(filename, ifstream::in);
    if(!utils::checkFile(file))
        return nullptr;

    /* Store input file name for posterior uses in other formats */
    _alignment->filename.append("!Title ");
    _alignment->filename.append(filename);
    _alignment->filename.append(";");

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
    _alignment->seqsName  = new string[_alignment->sequenNumber];
    _alignment->sequences = new string[_alignment->sequenNumber];
    _alignment->seqsInfo  = NULL;//new string[_alignment->sequenNumber];

    for(i = -1; (i < _alignment->sequenNumber) && (!file.eof()); ) {

        /* Deallocate previously used dinamic memory */
        if (line != NULL)
            delete [] line;

        /* Read lines in a safe way */
        line = utils::readLine(file);
        if (line == NULL)
            continue;

        /* Store original header fom input sequences including non-standard
         * characters */
//         if (line[0] == '>')
//             _alignment->seqsInfo[i+1].append(&line[1], strlen(line) - 1);

        /* Cut the current line and check whether there are valid characters */
        str = strtok(line, OTHDELIMITERS);
        if (str == NULL)
            continue;

        /* Check whether current line belongs to the current sequence
         * or it is a new one. In that case, store the sequence name */
        if(str[0] == '>') {
            /* Move sequence name pointer until a valid string name is obtained */
            do {
                str = str + 1;
            } while(strlen(str) == 0);
            _alignment->seqsName[++i].append(str, strlen(str));
            continue;
        }

        /* Sequence */
        while(str != NULL) {
            _alignment->sequences[i].append(str, strlen(str));
            str = strtok(NULL, DELIMITERS);
        }
    }

    /* Close the input file */
    file.close();

    /* Deallocate previously used dinamic memory */
    if (line != NULL)
        delete [] line;
        
    /* Check the matrix's content */
    _alignment->fillMatrices(false);
    _alignment->originalSequenNumber = _alignment-> sequenNumber;
    _alignment->originalResidNumber = _alignment->residNumber;
    return _alignment; 
}

bool FastaState::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
    /* Generate output alignment in FASTA format. Sequences can be unaligned. */

    int i, j, k, maxLongName;
    string *tmpMatrix;
    bool lastcharIsnewline = false;

    /* Allocate local memory for generating output alignment */
    tmpMatrix = new string[alignment->originalSequenNumber];

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    for(i = 0; i < alignment->originalSequenNumber; i++)
        tmpMatrix[i] = (!Machine->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    /* Depending on if short name flag is activated (limits sequence name up to
     * 10 characters) or not, get maximum sequence name length. Consider those
     * cases when the user has asked to keep original sequence header */
    maxLongName = 0;
    for(i = 0; i < alignment->originalSequenNumber; i++)
    {
        if (alignment->saveSequences && alignment->saveSequences[i] == -1) continue;
        if (!Machine->keepHeader)
            maxLongName = utils::max(maxLongName, alignment->seqsName[i].size());
        else if (alignment->seqsInfo != NULL)
            maxLongName = utils::max(maxLongName, alignment->seqsInfo[i].size());
    }

    if (Machine->shortNames && maxLongName > PHYLIPDISTANCE) {
        maxLongName = PHYLIPDISTANCE;
        if (!Machine->keepHeader)
            Debug.Report(WarningCode::HeaderWillBeCut);
    }
    /* Print alignment. First, sequences name id and then the sequences itself */
    for(i = 0; i < alignment->originalSequenNumber; i++) {
        
        if (alignment->saveSequences != NULL && alignment->saveSequences[i] == -1) continue;
        
        if (!Machine->keepHeader)
            (*output) << ">" << alignment->seqsName[i].substr(0, maxLongName) << endl;
        
        else if (alignment->seqsInfo != NULL)
            (*output) << ">" << alignment->seqsInfo[i].substr(0, maxLongName) << endl;
        
        
        for (j = 0, k = 0; j < alignment->sequences[i].length(); j++)
        {
            if (alignment->saveResidues != NULL && alignment->saveResidues[j] == -1) 
            {
                if (!lastcharIsnewline && j == alignment->sequences[i].length() -1 ) 
                {
                    (*output) << endl;
                    lastcharIsnewline = true;
                }
            }
            else
            {
                (*output) << tmpMatrix[i][j];
                k++;
                lastcharIsnewline = false;
                if (!lastcharIsnewline && (k % 60 == 0 || j == alignment->sequences[i].length() -1 ))
                {
                    (*output) << endl;
                    lastcharIsnewline = true;
                }
            }
        }
    }

    /* Deallocate local memory */
    delete [] tmpMatrix;
    
    return true;
}

bool FastaState::RecognizeOutputFormat(std::string FormatName)
{
    if (ReadWriteBaseState::RecognizeOutputFormat(FormatName)) return true;
    if (FormatName == "fasta") return true;
    return false;
}
