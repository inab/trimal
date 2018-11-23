#include "../../include/FormatHandling/mega_sequential_state.h"

#include "../../include/FormatHandling/FormatManager.h"
#include "../../include/defines.h"
#include "../../include/utils.h"

int mega_sequential_state::CheckAlignment(std::istream* origin)
{
    origin->seekg(0);
    origin->clear();
 
    char c, *firstWord = nullptr, *line = nullptr;
    int blocks = 0;
    std::string nline;

    /* Read first valid line in a safer way */
    do {
        line = utils::readLine(*origin);
    } while ((line == nullptr) && (!origin->eof()));

    /* If the file end is reached without a valid line, warn about it */
    if (origin->eof())
        return 0;

    /* Otherwise, split line */
    firstWord = strtok(line, OTHDELIMITERS);

        /* Mega Format */
    if((!strcmp(firstWord, "#MEGA")) || (!strcmp(firstWord, "#mega"))) {

        /* Determine specific mega format: sequential or interleaved.
         * Counting the number of blocks (set of lines starting by "#") in
         * the input file. */
        blocks = 0;
        do {
            origin->read(&c, 1);
        } while((c != '#') && (!origin->eof()));

        do {
            while((c != '\n') && (!origin->eof()))
                origin->read(&c, 1);
            origin->read(&c, 1);
            if(c == '#')
                blocks++;
        } while((c != '\n') && (!origin->eof()));

        /* MEGA Sequential (22) or Interleaved (21) */
        return (!blocks) ? 1 : 0;
    }
    delete[] line;
    return 0;
}

Alignment* mega_sequential_state::LoadAlignment(std::string& filename)
{
    Alignment * _alignment = new Alignment();
   /* MEGA sequential file format parser */

    char *frag = nullptr, *str = nullptr, *line = nullptr;
    std::ifstream file;
    int i;

    /* Check the file and its content */
    file.open(filename, std::ifstream::in);
    if(!utils::checkFile(file))
        return nullptr;

    /* Filename is stored as a title for MEGA input alignment.
     * If it is detected later a label "TITLE" in input file, this information
     * will be replaced for that one */
    _alignment->filename.append("!Title ");
    _alignment->filename.append(filename);
    _alignment->filename.append(";");

    /* Skip first valid line */
    do {
        line = utils::readLine(file);
    } while ((line == nullptr) && (!file.eof()));

    /* If the file end is reached without a valid line, warn about it */
    if (file.eof())
        return nullptr;

    /* Try to get input alignment information */
    while(!file.eof()) {

        /* Destroy previously allocated memory */
        if (line != nullptr)
            delete [] line;

        /* Read a new line in a safe way */
        line = utils::readLine(file);
        if (line == nullptr)
            continue;

        /* If a sequence name flag is found, go out from getting information loop */
        if(!strncmp(line, "#", 1))
            break;

        /* Destroy previously allocated memory */
        if (frag != nullptr)
            delete [] frag;

        /* Create a local copy from input line */
        frag = new char[strlen(line) + 1];
        strcpy(frag, line);

        /* Split input line copy into pieces and analize it
         * looking for specific labels */
        str = strtok(frag, "!: ");
        for(i = 0; i < (int) strlen(str); i++)
            str[i] = toupper(str[i]);

        /* If TITLE label is found, replace previously stored information with
         * this info */
        if(!strcmp(str, "TITLE")) {
             _alignment->filename.clear();
            if(strncmp(line, "!", 1))
                 _alignment->filename += "!";
             _alignment->filename += line;
        }

            /* If FORMAT label is found, try to get some details from input file */
        else if(!strcmp(str, "FORMAT"))
            _alignment -> alignmentInfo.append(line, strlen(line));
    }

    /* Deallocate local memory */
    delete [] frag;

    /* Count how many sequences are in input alignment file */
    do {

        /* Check whether input line is valid or not */
        if (line == nullptr) {
            line = utils::readLine(file);
            continue;
        }

        /* If current line starts by a # means that it is a sequence name */
        if (!strncmp(line, "#", 1))
            _alignment->numberOfSequences++;

        /* Destroy previously allocated memory */
        delete [] line;

        /* Read a new line in a safe way */
        line = utils::readLine(file);

    } while(!file.eof());

    /* Move file pointer to the beginner */
    file.clear();
    file.seekg(0);

    /* Allocate memory */
    _alignment->seqsName  = new std::string[_alignment->numberOfSequences];
    _alignment->sequences = new std::string[_alignment->numberOfSequences];

    /* Skip first line */
    line = utils::readLine(file);

    /* Skip lines until first sequence name is found */
    while(!file.eof()) {

        /* Destroy previously allocated memory */
        delete [] line;

        /* Read a new line in a safe way */
        line = utils::readLine(file);
        if (line == nullptr)
            continue;

        /* If sequence name label is found, go out from loop */
        if (!strncmp(line, "#", 1))
            break;
    }

    /* First sequence is already detected so its name should be stored */
    i = -1;

    /* This loop is a bit tricky because first sequence name has been already
     * detected, so it is necessary to process it before moving to next line.
     * That implies that lines are read at loop ends */
    while(!file.eof()) {

        /* Skip blank lines */
        if (line == nullptr) {
            line = utils::readLine(file);
            continue;
        }

        /* Skip lines with comments */
        if (!strncmp(line, "!", 1)) {
            /* Deallocate memory and read a new line */
            delete [] line;
            line = utils::readLine(file);
            continue;
        }

        /* Remove comments inside of sequences and split input line */
        frag = utils::trimLine(line);

        /* Skip lines with only comments */
        if (frag == nullptr) {

            /* Deallocate memory and read a new line */
            delete [] line;
            line = utils::readLine(file);
            continue;
        }

        /* Otherwise, split it into fragments */
        str = strtok(frag, " #\n");

        /* Sequence Name */
        if (!strncmp(line, "#", 1)) {
            i += 1;
            _alignment->seqsName[i].append(str, strlen(str));
            str = strtok(nullptr, " #\n");
        }

        /* Sequence itself */
        while(str != nullptr) {
            _alignment->sequences[i].append(str, strlen(str));
            str = strtok(nullptr, " \n");
        }

        /* Deallocate dynamic memory */
        delete [] frag;

        delete [] line;

        /* Read a new line in a safe way */
        line = utils::readLine(file);
    }

    /* Close input file */
    file.close();

    /* Deallocate dynamic memory */
    delete [] line;

    /* Check the matrix's content */
    _alignment->fillMatrices(true);
    _alignment->originalNumberOfSequences = _alignment-> numberOfSequences;
    _alignment->originalNumberOfResidues = _alignment->numberOfResidues;
    return _alignment;
}

bool mega_sequential_state::SaveAlignment(Alignment* alignment, std::ostream* output, std::string* FileName)
{
    /* Generate output alignment in MEGA format */

    int i, j, k, l, m;
    std::string *tmpMatrix;

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!alignment->isAligned || !alignment->saveResidues || !alignment->saveSequences) {
        debug.report(ErrorCode::UnalignedAlignmentToAlignedFormat, new std::string[1] { this->name });
        return false;
    }

    /* Allocate local memory for generating output alignment */
    tmpMatrix = new std::string[alignment->originalNumberOfSequences];

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    for(i = 0; i < alignment->originalNumberOfSequences; i++)
        tmpMatrix[i] = (!Machine->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    /* Compute output file datatype */
    alignment->getAlignmentType();

    /* Print output alignment header */
    *output << "#MEGA\n" << alignment->filename << "\n";

    /* Print alignment datatype */
    if (alignment->getAlignmentType() & SequenceTypes::DNA)
        *output << "!Format DataType=DNA ";
    else if (alignment->getAlignmentType() & SequenceTypes::RNA)
        *output << "!Format DataType=RNA ";
    else if (alignment->getAlignmentType() & SequenceTypes::AA)
        *output << "!Format DataType=protein ";

    /* Print number of sequences and alignment length */
    *output << "NSeqs=" << alignment->numberOfSequences << " Nsites=" << alignment->numberOfResidues
         << " indel=- CodeTable=Standard;\n";

    /* Print sequences name and sequences divided into blocks of 50 residues */
    for(i = 0; i < alignment->originalNumberOfSequences; i++) {
        if (alignment->saveSequences[i] != -1)
        {
            *output << "\n#" << alignment->seqsName[i] << "\n";
            
            for(j = 0, l = 0; j < alignment->sequences[i].length(); j ++) 
            {
                if (alignment->saveResidues[j] != -1)
                {
                    *output << alignment->sequences[i][j];
                    if (++l % 10 == 0) *output << " ";
                    if (l == 50)
                    {
                        l = 0;
                        *output << "\n";
                    }
                }
            }
            if (l % 10 != 0) *output << " ";
            if (l != 0) *output << "\n";
        }
    }
    *output << "\n";

    /* Deallocate local memory */
    delete [] tmpMatrix;
    
    return true;
}

bool mega_sequential_state::RecognizeOutputFormat(std::string& FormatName)
{
    if (BaseFormatHandler::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "mega";
}
