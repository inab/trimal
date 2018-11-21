#include "../../include/ReadWriteMS/nexus_state.h"

#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/defines.h"
#include "../../include/utils.h"

int nexus_state::CheckAlignment(std::istream* origin)
{
    origin->seekg(0);
    origin->clear();
    char *firstWord = nullptr, *line = nullptr;
    
    /* Read first valid line in a safer way */
    do {
        line = utils::readLine(*origin);
    } while ((line == nullptr) && (!origin->eof()));

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

newAlignment* nexus_state::LoadAlignment(std::string& filename)
{
    newAlignment* _alignment = new newAlignment();
    /* NEXUS file format parser */
    char *frag = nullptr, *str = nullptr, *line = nullptr;
    int i, pos, state, firstBlock;
    std::ifstream file;

    /* Check the file and its content */
    file.open(filename, std::ifstream::in);
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
        delete [] line;

        /* Read line in a safer way */
        line = utils::readLine(file);
        if (line == nullptr)
            continue;

        /* Discard line where there is not information */
        str = strtok(line, DELIMITERS);
        if(str == nullptr)
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
            str = strtok(nullptr, DELIMITERS);
            while(str != nullptr) {
                _alignment-> aligInfo.append(str, strlen(str));
                _alignment-> aligInfo.append(" ", strlen(" "));
                str = strtok(nullptr, DELIMITERS);
            }
        }

            /* In this case, try to get matrix dimensions */
        else if((!strcmp(str, "DIMENSIONS")) && state) {
            str = strtok(nullptr, DELIMITERS);
            frag = strtok(nullptr, DELIMITERS);
            str = strtok(str, "=;");
            _alignment->sequenNumber = atoi(strtok(nullptr, "=;"));
            frag = strtok(frag, "=;");
            _alignment->residNumber = atoi(strtok(nullptr, "=;"));
        }
    } while(!file.eof());

    /* Check all parameters */
    if(strcmp(str, "MATRIX") || (_alignment->sequenNumber == 0) || (_alignment->residNumber == 0))
        return nullptr;

    /* Allocate memory for the input alignment */
    _alignment->seqsName  = new std::string[_alignment->sequenNumber];
    _alignment->sequences = new std::string[_alignment->sequenNumber];

    pos = 0;
    state = false;
    firstBlock = true;

    while(!file.eof()) {
        /* Destroy previous assigned memory */
        delete [] line;

        /* Read line in a safer way */
        line = utils::readLine(file);
        if (line == nullptr)
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
        if ((state) || (i != (int) strlen(line)))
            continue;

        /* Check for a specific tag indicating matrix end */
        if((!strncmp(line, "end;", 4)) || (!strncmp(line, "END;", 4)))
            break;

        /* Split input line and check it if it is valid */
        str = strtok(line, OTH2DELIMITERS);
        if (str == nullptr)
            continue;

        /* Store the sequence name, only from the first block */
        if(firstBlock)
            _alignment->seqsName[pos].append(str, strlen(str));

        /* Store rest of line as part of sequence */
        str = strtok(nullptr, OTH2DELIMITERS);
        while(str != nullptr) {
            _alignment->sequences[pos].append(str, strlen(str));
            str = strtok(nullptr, OTH2DELIMITERS);
        }

        /* Move sequences pointer to next one. It if it is last one, move it to
         * the beginning and set the first block to false for avoiding to rewrite
         * sequences name */
        pos = (pos + 1) % _alignment->sequenNumber;
        if (not pos)
            firstBlock = false;
    }

    /* Deallocate memory */
    delete [] line;

    /* Close the input file */
    file.close();

    /* Check the matrix's content */
    _alignment->fillMatrices(true);
    _alignment->originalSequenNumber = _alignment-> sequenNumber;
    _alignment->originalResidNumber =_alignment->residNumber;
    return _alignment;
}

bool nexus_state::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
    /* Generate output alignment in NEXUS format setting only alignment block */

    int i, j, k, l, maxLongName = 0;
    std::string *tmpMatrix;

    // Check whether sequences in the alignment are aligned or not.
    // Warn about it if there are not aligned.
    if (!alignment->isAligned) {
        debug.report(ErrorCode::UnalignedAlignmentToAlignedFormat, new std::string[1] { this->name });
        return false;
    }

    // Allocate local memory for generating output alignment
    tmpMatrix = new std::string[alignment->originalSequenNumber];

    // Depending on alignment orientation: forward or reverse. Copy directly
    // sequence information or get firstly the reversed sequences and then
    // copy it into local memory
    for(i = 0; i < alignment->originalSequenNumber; i++)
        tmpMatrix[i] = (!Machine->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    // Compute maximum sequences name length
    for (i = 0; (i < alignment->originalSequenNumber); i++)
        if (alignment->saveSequences[i] != -1)
            maxLongName = utils::max(maxLongName, alignment->seqsName[i].size());

    // Compute output file datatype
    alignment->getAlignmentType();

    // Remove characters like ";" from input alignment information line
    while((int) alignment -> aligInfo.find(';') != (int) std::string::npos)
        alignment ->aligInfo.erase(alignment -> aligInfo.find(';'), 1);

    // Print Alignment header
    *output << "#NEXUS\nBEGIN DATA;\n DIMENSIONS NTAX="
         << alignment->sequenNumber << " NCHAR=" << alignment->residNumber <<";\n";

    // Print alignment datatype
    if (alignment->getAlignmentType() & SequenceTypes::DNA)
        *output << "FORMAT DATATYPE=DNA INTERLEAVE=yes GAP=-";
    else if (alignment->getAlignmentType() & SequenceTypes::RNA)
        *output << "FORMAT DATATYPE=RNA INTERLEAVE=yes GAP=-";
    else if (alignment->getAlignmentType() & SequenceTypes::AA)
        *output << "FORMAT DATATYPE=PROTEIN INTERLEAVE=yes GAP=-";

    i = 0;
    // Using information from input alignment. Use only some tags.
    while((j = alignment ->aligInfo.find(' ', i)) != (int) std::string::npos) {

        if((alignment ->aligInfo.substr(i, j - i)).compare(0, 7, "MISSING") == 0 ||
           (alignment ->aligInfo.substr(i, j)).compare(0, 7, "missing") == 0)
            *output << " " << (alignment ->aligInfo.substr(i, j - i));

        else if((alignment ->aligInfo.substr(i, j)).compare(0, 9, "MATCHCHAR") == 0 ||
                (alignment ->aligInfo.substr(i, j)).compare(0, 9, "matchchar") == 0)
            *output << " " << (alignment ->aligInfo.substr(i, j - i));

        i = j + 1;
    }
    *output << ";\n";

    // Add a header indicating the number of residues of each sequence.
    for(i = 0; i < alignment->originalSequenNumber; i++)
    {
        if (alignment->saveSequences[i] == -1 )
            continue;

        *output << "[Name: "
                << std::setw(maxLongName + 4) << std::left << alignment->seqsName[i]
                << "Len: " << alignment->residNumber << "]\n";
    }
    *output << "\nMATRIX";

    // Start filling the file with sequence names and sequences.
    for (j = 0, k = 0; j < alignment->originalResidNumber;)
    {
        // Move until next not rejected residue
        if (alignment->saveResidues[j] == -1)
        {
            j++;
            continue;
        }

        // Iterate over the sequences
        for (i = 0; i < alignment->originalSequenNumber; i++)
        {
            // Skip rejected sequences
            if (alignment->saveSequences[i] == -1) continue;
            // Output sequence name
            *output << "\n" << std::setw(maxLongName + 5) << std::left << alignment->seqsName[i];
            // Add residues per block.
            // k = residue position;
            // l = residues added on current line
            for (k = j, l = 0; k < alignment->originalResidNumber && l < 50;)
            {
                // Don't save residues that have been rejected
                if (alignment->saveResidues[k] == -1)
                {
                    k++;
                    continue;
                }
                // Add the residue at position K
                *output << alignment->sequences[i][k];
                k++; l++;
                // If a block of 10 residues have been added, append a space to sequence
                if (l % 10 == 0 && l != 50) *output << " ";
            }
        }
        // Add the line to split blocks
        *output << "\n";
        j = k;
    }

    // End of MATRIX
    *output << "\n;\nEND;\n";

    /* Deallocate local memory */
    delete [] tmpMatrix;
    
    return true;
}

bool nexus_state::RecognizeOutputFormat(std::string& FormatName)
{
    if (ReadWriteBaseState::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "nexus";
}

