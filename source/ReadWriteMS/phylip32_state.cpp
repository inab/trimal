#include "../../include/ReadWriteMS/phylip32_state.h"

#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/defines.h"
#include "../../include/utils.h"

int phylip32_state::CheckAlignment(std::istream* origin)
{
    origin->seekg(0);
    origin->clear();
    char *firstWord = nullptr, *line = nullptr;
    int blocks = 0;
    std::string nline;
    
    /* Read first valid line in a safer way */
    do {
        line = utils::readLine(*origin);
    } while ((line == nullptr) && (!origin->eof()));

    /* If the file end is reached without a valid line, warn about it */
    if (origin->eof())
    {
        delete [] line;
        return false;
    }

    /* Otherwise, split line */
    firstWord = strtok(line, OTHDELIMITERS);

    /* Phylip Format */
    {
        /* Determine specific phylip format: sequential or interleaved. */

        /* Get number of sequences and residues */
        int sequenNumber = atoi(firstWord);
        int residNumber = 0;
        firstWord = strtok(nullptr, DELIMITERS);
        if(firstWord != nullptr)
            residNumber = atoi(firstWord);
        else 
        {
            delete [] line;
            return false;
        }

        /* If there is only one sequence, use by default sequential format since
         * it is impossible to determine exactly which phylip format is */
        if((sequenNumber == 1) && (residNumber != 0))
        {
            delete [] line;
            return false;
        }

            /* If there are more than one sequence, analyze sequences distribution to
             * determine its format. */
        else if((sequenNumber != 0) && (residNumber != 0)) {
            blocks = 0;

            /* Read line in a safer way */
            delete [] line;
            do {
                line = utils::readLine(*origin);
            } while ((line == nullptr) && (!origin->eof()));

            /* If the file end is reached without a valid line, warn about it */
            if (origin->eof())
                return 0;

            firstWord = strtok(line, DELIMITERS);
            while(firstWord != nullptr) {
                blocks++;
                firstWord = strtok(nullptr, DELIMITERS);
            }

            delete [] line;
            /* Read line in a safer way */
            do {
                line = utils::readLine(*origin);
            } while ((line == nullptr) && (!origin->eof()));

            firstWord = strtok(line, DELIMITERS);
            while(firstWord != nullptr) {
                blocks--;
                firstWord = strtok(nullptr, DELIMITERS);
            }
            delete [] line;
            /* If the file end is reached without a valid line, warn about it */
            if (origin->eof())
            {
                return false;
            }

            /* Phylip Interleaved (12) or Sequential (11) */
            return (!blocks) ? 0 : 1;
        }
    }
    delete[] line;
    return 0;
}

newAlignment* phylip32_state::LoadAlignment(std::string filename)
{
    /* PHYLIP 3.2 (Interleaved) file format parser */
    newAlignment* _alignment = new newAlignment();
    
    int i, blocksFirstLine, firstLine = true;
    char *str, *line = nullptr;
    std::ifstream file;

    /* Check the file and its content */
    file.open(filename, std::ifstream::in);
    if(!utils::checkFile(file))
        return nullptr;

    /* Store the file name for futher format conversion*/
    _alignment->filename.append("!Title ");
    _alignment->filename.append(filename);
    _alignment->filename.append(";");

    /* Read first valid line in a safer way */
    do {
        line = utils::readLine(file);
    } while ((line == nullptr) && (!file.eof()));

    /* If the file end is reached without a valid line, warn about it */
    if (file.eof())
        return nullptr;

    /* Get the sequences and residues numbers. If there is any mistake,
     * return a FALSE value to warn about the possible error */
    str = strtok(line, DELIMITERS);
    _alignment->sequenNumber = 0;
    if(str != nullptr)
        _alignment->sequenNumber = atoi(str);

    str = strtok(nullptr, DELIMITERS);
    _alignment->residNumber = 0;
    if(str != nullptr)
        _alignment->residNumber = atoi(str);

    if((_alignment->sequenNumber == 0) || (_alignment->residNumber == 0))
        return nullptr;

    /* Reserve memory according to the input parameters */
    _alignment->sequences  = new std::string[_alignment->sequenNumber];
    _alignment->seqsName   = new std::string[_alignment->sequenNumber];

    /* Point to the first sequence in the alignment. Since the alignment could not
     * have blank lines to separate the different sequences. Store the blocks size
     * for the first line including a sequence identifier */
    i = 0;
    blocksFirstLine = 0;

    do {
        /* Read lines in a safer way. Destroy previous stored information */
        delete [] line;

        line = utils::readLine(file);
        /* If there is nothing in the input line, skip the loop instructions */
        if(line == nullptr)
            continue;

        str = strtok(line, OTHDELIMITERS);
        /* First block: Sequence Name + Sequence fragment. Count how many blocks
         * the first sequence line is divided. It could help to identify the
         * different sequences from the input file */
        if(firstLine) {
            _alignment->seqsName[i].append(str, strlen(str));
            str = strtok(nullptr, OTHDELIMITERS);
            firstLine = 1;
        }

        /* Sequence fragment */
        while(str != nullptr) {
            _alignment->sequences[i].append(str, strlen(str));
            str = strtok(nullptr, OTHDELIMITERS);
            /* Count the blocks number for the sequences first line */
            if (firstLine)
                firstLine += 1;
        }

        /* Store the blocks number for the first sequence including the name */
        if ((blocksFirstLine == 0) && firstLine)
            blocksFirstLine = firstLine;

        /* If a false positive new sequence was detected, add stored information for
         * the current sequence to the previous one and clear the data structure
         * for the current sequence. Finally, move the sequence pointer to the
         * previous one. */
        if ((firstLine ) && (firstLine != blocksFirstLine)) {
            _alignment->sequences[i-1].append(_alignment->seqsName[i]);
            _alignment->seqsName[i].clear();
            _alignment->sequences[i-1].append(_alignment->sequences[i]);
            _alignment->sequences[i].clear();
            i --;
        }

        firstLine = false;
        /* There are many ways to detect a new sequence. */
        /* One of them -experimental- is just to detect if the residues number for
         * the current entry is equal to the residues number for the whole align */
        if ((int) _alignment->sequences[i].size() == _alignment->residNumber) {
            firstLine = true;
            i++;
        }
    } while(!file.eof());

    /* Close the input file and delete dinamic memory */
    file.close();
    delete [] line;

    /* Check the matrix's content */
    _alignment->fillMatrices(true);
    _alignment->originalSequenNumber = _alignment-> sequenNumber;
    _alignment->originalResidNumber = _alignment->residNumber;
    return _alignment; 
}

bool phylip32_state::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
    /* Generate output alignment in PHYLIP 3.2 format (interleaved) */

    int i, j, k, maxLongName;
    std::string *tmpMatrix;

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!alignment->isAligned) {
        debug.report(ErrorCode::UnalignedAlignmentToAlignedFormat, new std::string[1] { this->name });
        return false;
    }

    /* Allocate local memory for generating output alignment */
    tmpMatrix = new std::string[alignment->originalSequenNumber];

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    for(i = 0; i < alignment->originalSequenNumber; i++)
        tmpMatrix[i] = (!Machine->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    maxLongName = PHYLIPDISTANCE;
    for(i = 0; (i < alignment->originalSequenNumber); i++)
        if  (alignment->saveSequences[i] != -1)
            maxLongName = utils::max(maxLongName, alignment->seqsName[i].size());

    /* Generating output alignment */
    /* First Line: Sequences Number & Residued Number */
    (*output) << " " << alignment->sequenNumber << " " << alignment->residNumber;

    /* Alignment */
    /* For each sequence, print its identifier and then the sequence itself in
     * blocks of 50 residues */
    for(i = 0; i < alignment->originalSequenNumber; i++) {
        /* Sequence Name */
        if (alignment->saveSequences[i] == -1) continue;
        (*output) << "\n" << std::setw(maxLongName + 3) << std::left << alignment->seqsName[i].substr(0, maxLongName);
        /* Sequence. Each line contains a block of 5 times 10 residues. */
        
        for (j = 0, k = 0; j < alignment->originalResidNumber; j++)
        {
            if (alignment->saveResidues[j] == -1) continue;
            if (k == 50)
            {
                *output << "\n" << std::setw(maxLongName + 3) << std::left << " " ;
                k = 0;
            }
            *output << alignment->sequences[i][j];
            k++;
            if (k % 10 == 0) 
                *output << " ";
        }
        if (k % 10 != 0)
            *output << " ";
        
        /* Print a blank line to mark sequences separation */
//         if (k % 50 != 0)
            (*output) << "\n";
    }
    *output << "\n";

    /* Deallocate local memory */
    delete [] tmpMatrix;
    
    return true;
}

bool phylip32_state::RecognizeOutputFormat(std::string FormatName)
{
    if (ReadWriteBaseState::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "phylip32";
}

