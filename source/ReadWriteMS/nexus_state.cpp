#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/ReadWriteMS/nexus_state.h"
#include "../../include/defines.h"
#include "../../include/values.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include "../../include/newAlignment.h"

using namespace std;

int NexusState::CheckAlignment(istream* origin)
{
    origin->seekg(0);
    origin->clear();
    char *firstWord = NULL, *line = NULL;
    
    /* Read first valid line in a safer way */
    do {
        line = utils::readLine(*origin);
    } while ((line == NULL) && (!origin->eof()));

    /* If the file end is reached without a valid line, warn about it */
    if (origin->eof())
        return false;

    /* Otherwise, split line */
    firstWord = strtok(line, OTHDELIMITERS);

    /* Clustal Format */
    if((!strcmp(firstWord, "#NEXUS")) || (!strcmp(firstWord, "#nexus")))
        return 1;

    delete[] line;
 
    return 0;
}

newAlignment* NexusState::LoadAlignment(std::string filename)
{
    newAlignment* _alignment = new newAlignment();
    /* NEXUS file format parser */
    char *frag = NULL, *str = NULL, *line = NULL;
    int i, pos, state, firstBlock;
    ifstream file;

    /* Check the file and its content */
    file.open(filename, ifstream::in);
    if(!utils::checkFile(file))
        return nullptr;

    /* Store input file name for posterior uses in other formats */
    /* We store the file name */
    _alignment->filename.append("!Title ");
    _alignment->filename.append(filename);
    _alignment->filename.append(";");

    state = false;
    do {

        /* Destroy previous assigned memory */
        if (line != NULL)
            delete [] line;

        /* Read line in a safer way */
        line = utils::readLine(file);
        if (line == NULL)
            continue;

        /* Discard line where there is not information */
        str = strtok(line, DELIMITERS);
        if(str == NULL)
            continue;

        /* If the line has any kind of information, try to catch it */
        /* Firstly, convert to capital letters the input line */
        for(i = 0; i < (int) strlen(str); i++)
            str[i] = toupper(str[i]);

        /* ... and then compare it again specific tags */
        if(!strcmp(str, "BEGIN"))
            state = true;

        else if(!strcmp(str, "MATRIX"))
            break;

            /* Store information about input format file */
        else if(!strcmp(str, "FORMAT")) {
            str = strtok(NULL, DELIMITERS);
            while(str != NULL) {
                _alignment-> aligInfo.append(str, strlen(str));
                _alignment-> aligInfo.append(" ", strlen(" "));
                str = strtok(NULL, DELIMITERS);
            }
        }

            /* In this case, try to get matrix dimensions */
        else if((!strcmp(str, "DIMENSIONS")) && state) {
            str = strtok(NULL, DELIMITERS);
            frag = strtok(NULL, DELIMITERS);
            str = strtok(str, "=;");
            _alignment->sequenNumber = atoi(strtok(NULL, "=;"));
            frag = strtok(frag, "=;");
            _alignment->residNumber = atoi(strtok(NULL, "=;"));
        }
    } while(!file.eof());

    /* Check all parameters */
    if(strcmp(str, "MATRIX") || (_alignment->sequenNumber == 0) || (_alignment->residNumber == 0))
        return nullptr;

    /* Allocate memory for the input alignmet */
    _alignment->seqsName  = new string[_alignment->sequenNumber];
    _alignment->sequences = new string[_alignment->sequenNumber];

    pos = 0;
    state = false;
    firstBlock = true;

    while(!file.eof()) {
        /* Destroy previous assigned memory */
        if (line != NULL)
            delete [] line;

        /* Read line in a safer way */
        line = utils::readLine(file);
        if (line == NULL)
            continue;

        /* Discard any comments from input file */
        for(i = 0; i < (int) strlen(line); i++) {
            if (line[i] == '[')
                state = true;
            else if (line[i] == ']' && state) {
                state = false;
                break;
            }
        }

        /* If there is a multi-line comments, skip it as well */
        if ((state) || (not state && i != (int) strlen(line)))
            continue;

        /* Check for a specific tag indicating matrix end */
        if((!strncmp(line, "end;", 4)) || (!strncmp(line, "END;", 4)))
            break;

        /* Split input line and check it if it is valid */
        str = strtok(line, OTH2DELIMITERS);
        if (str == NULL)
            continue;

        /* Store the sequence name, only from the first block */
        if(firstBlock)
            _alignment->seqsName[pos].append(str, strlen(str));

        /* Store rest of line as part of sequence */
        str = strtok(NULL, OTH2DELIMITERS);
        while(str != NULL) {
            _alignment->sequences[pos].append(str, strlen(str));
            str = strtok(NULL, OTH2DELIMITERS);
        }

        /* Move sequences pointer to next one. It if it is last one, move it to
         * the beginning and set the first block to false for avoiding to rewrite
         * sequences name */
        pos = (pos + 1) % _alignment->sequenNumber;
        if (not pos)
            firstBlock = false;
    }

    /* Deallocate memory */
    if (line != NULL)
        delete [] line;

    /* Close the input file */
    file.close();

    /* Check the matrix's content */
    _alignment->fillMatrices(true);
    return _alignment;
}

bool NexusState::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
    /* Generate output alignment in NEXUS format setting only alignment block */

    int i, j, k, maxLongName = 0;
    string *tmpMatrix;

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!alignment->isAligned) {
        cerr << endl << "ERROR: Sequences are not aligned. Format (phylip) "
             << "not compatible with unaligned sequences." << endl << endl;
        return false;
    }

    /* Allocate local memory for generating output alignment */
    tmpMatrix = new string[alignment->sequenNumber];

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    for(i = 0; i < alignment->sequenNumber; i++)
        tmpMatrix[i] = (!Machine->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    /* Compute maximum sequences name length */
    for(i = 0; (i < alignment->sequenNumber) && (!Machine->shortNames); i++)
        maxLongName = utils::max(maxLongName, alignment->seqsName[i].size());

    /* Compute output file datatype */
    alignment->getAlignmentType();

    /* Remove characters like ";" from input alignment information line */
    while((int) alignment -> aligInfo.find(";") != (int) string::npos)
        alignment ->aligInfo.erase(alignment -> aligInfo.find(";"), 1);

    /* Print Alignment header */
    *output << "#NEXUS" << endl << "BEGIN DATA;" << endl << " DIMENSIONS NTAX="
         << alignment->sequenNumber << " NCHAR=" << alignment->residNumber <<";" << endl;

    /* Print alignment datatype */
    if (alignment->dataType & SequenceTypes::DNA)
        *output << "FORMAT DATATYPE=DNA INTERLEAVE=yes GAP=-";
    else if (alignment->dataType & SequenceTypes::RNA)
        *output << "FORMAT DATATYPE=RNA INTERLEAVE=yes GAP=-";
    else if (alignment->dataType & SequenceTypes::AA)
        *output << "FORMAT DATATYPE=PROTEIN INTERLEAVE=yes GAP=-";

    i = 0;
    /* Using information from input alignment. Use only some tags. */
    while((j = alignment ->aligInfo.find(" ", i)) != (int) string::npos) {

        if((alignment ->aligInfo.substr(i, j - i)).compare(0, 7, "MISSING") == 0 ||
           (alignment ->aligInfo.substr(i, j)).compare(0, 7, "missing") == 0)
            *output << " " << (alignment ->aligInfo.substr(i, j - i));

        else if((alignment ->aligInfo.substr(i, j)).compare(0, 9, "MATCHCHAR") == 0 ||
                (alignment ->aligInfo.substr(i, j)).compare(0, 9, "matchchar") == 0)
            *output << " " << (alignment ->aligInfo.substr(i, j - i));

        i = j + 1;
    }
    *output << ";" << endl;

    /* Print sequence name and sequence length */
    for(i = 0; i < alignment->sequenNumber; i++)
        *output << "[Name: " << setw(maxLongName + 4) << left << alignment->seqsName[i] << "Len: "
             << alignment->residNumber << "]" << endl;
    *output << endl << "MATRIX" << endl;

    /* Print alignment itself. Sequence name and 50 residues blocks */
    for(j = 0; j < alignment->residNumber; j += 50) {
        for(i = 0; i < alignment->sequenNumber; i++) {
            *output << setw(maxLongName + 4) << left << alignment->seqsName[i];
            for(k = j; k < (j + 50) && k < alignment->residNumber; k += 10)
                *output << " " << alignment->sequences[i].substr(k, 10);
            *output << endl;
        }
        *output << endl;
    }
    *output << ";" << endl << "END;" << endl;

    /* Deallocate local memory */
    delete [] tmpMatrix;
    
    return true;
}

bool NexusState::RecognizeOutputFormat(std::string FormatName)
{
    if (ReadWriteBaseState::RecognizeOutputFormat(FormatName)) return true;
    if (FormatName == "nexus") return true;
    return false;
}

