#include "../../include/ReadWriteMS/mega_interleaved_state.h"

#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include "../../include/defines.h"

using namespace std;

int MegaInterleavedState::CheckAlignment(istream* origin)
{
    origin->seekg(0);
    origin->clear();
    
    char c, *firstWord = NULL, *line = NULL;
    int format = 0, blocks = 0;
    string nline;
    
    /* Read first valid line in a safer way */
    do {
        line = utils::readLine(*origin);
    } while ((line == NULL) && (!origin->eof()));

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
        return (!blocks) ? 0 : 1;
    }
    delete[] line;
    return 0;
}

newAlignment* MegaInterleavedState::LoadAlignment(std::__cxx11::string filename)
{
    newAlignment * _alignment = new newAlignment();
   /* MEGA interleaved file format parser */

    char *frag = NULL, *str = NULL, *line = NULL;
    int i, firstBlock = true;
    ifstream file;
    
    /* Check the file and its content */
    file.open(filename, ifstream::in);
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
    } while ((line == NULL) && (!file.eof()));

    /* If the file end is reached without a valid line, warn about it */
    if (file.eof())
        return nullptr;

    /* Try to get input alignment information */
    while(!file.eof()) {

        /* Destroy previously allocated memory */
        if (line != NULL)
            delete [] line;

        /* Read a new line in a safe way */
        line = utils::readLine(file);
        if (line == NULL)
            continue;

        /* If a sequence name flag is found, go out from getting information loop */
        if(!strncmp(line, "#", 1))
            break;

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
            filename.clear();
            if(strncmp(line, "!", 1))
                filename += "!";
            filename += line;
        }

            /* If FORMAT label is found, try to get some details from input file */
        else if(!strcmp(str, "FORMAT"))
             _alignment->aligInfo.append(line, strlen(line));

        /* Destroy previously allocated memory */
        if (frag != NULL)
            delete [] frag;
    }

    /* Count how many sequences are in input file */
    while(!file.eof()) {

        /* If a sequence name flag has been detected, increase counter */
        if(!strncmp(line, "#", 1))
            _alignment->sequenNumber++;

        /* Deallocate dynamic memory */
        if(line != NULL)
            delete [] line;

        /* Read lines in a safe way */
        line = utils::readLine(file);

        /* If a blank line is detected means first block of sequences is over */
        /* Then, break counting sequences loop */
        if (line == NULL)
            break;
    }

    /* Finish to preprocess the input file. */
    file.clear();
    file.seekg(0);

    /* Allocate memory */
    _alignment->seqsName  = new string[_alignment->sequenNumber];
    _alignment->sequences = new string[_alignment->sequenNumber];

    /* Skip first line */
    line = utils::readLine(file);

    /* Skip lines until first # flag is reached */
    while(!file.eof()) {

        /* Deallocate previously used dynamic memory */
        if (line != NULL)
            delete [] line;

        /* Read line in a safer way */
        line = utils::readLine(file);

        /* Determine whether a # flag has been found in current string */
        if (line != NULL)
            if(!strncmp(line, "#", 1))
                break;
    }

    /* Read sequences and get from first block, the sequences names */
    i = 0;
    firstBlock = true;

    while(!file.eof()) {

        if (line == NULL) {
            /* Read line in a safer way */
            line = utils::readLine(file);
            continue;
        }

        if (!strncmp(line, "!", 1)) {
            /* Deallocate memory and read a new line */
            delete [] line;
            line = utils::readLine(file);
            continue;
        }

        /* Trim lines from any kind of comments and split it */
        frag = utils::trimLine(line);
        str = strtok(frag, " #\n");

        /* Check whether a line fragment is valid or not */
        if (str == NULL)
            continue;

        /* Store sequences names if firstBlock flag is TRUE */
        if(firstBlock)
            _alignment->seqsName[i].append(str, strlen(str));

        /* Store sequence */
        str = strtok(NULL, " \n");
        while(str != NULL) {
            _alignment->sequences[i].append(str, strlen(str));
            str = strtok(NULL, " \n");
        }

        /* Deallocate previously used dynamic memory */
        if (frag != NULL)
            delete [] frag;

        if (line != NULL)
            delete [] line;

        /* Read line in a safer way */
        line = utils::readLine(file);

        i = (i + 1) % _alignment->sequenNumber;
        if (i == 0)
            firstBlock = false;
    }

    /* Close input file */
    file.close();

    /* Deallocate local memory */
    if (line != NULL)
        delete [] line;

    /* Check the matrix's content */
    _alignment->fillMatrices(true);
    return _alignment;
}

bool MegaInterleavedState::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
    return false;
}

bool MegaInterleavedState::RecognizeOutputFormat(std::string FormatName)
{
    return false;
}
