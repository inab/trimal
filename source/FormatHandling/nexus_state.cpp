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

#include "FormatHandling/nexus_state.h"

#include "FormatHandling/FormatManager.h"
#include "defines.h"
#include "utils.h"

namespace FormatHandling {
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
    {
        delete [] line;
        return false;
    }

    /* Otherwise, split line */
    firstWord = strtok(line, OTHDELIMITERS);

    /* Clustal Format */
    if((!strcmp(firstWord, "#NEXUS")) || (!strcmp(firstWord, "#nexus")))
    {
        delete[] line;
        return 1;
    }

    delete[] line;
    return 0;
}

Alignment* nexus_state::LoadAlignment(const std::string &filename)
{
    Alignment* alig = new Alignment();
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
    // alig->filename.append("!Title ");
    alig->filename.append(filename);
    alig->filename.append(";");

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
                alig-> alignmentInfo.append(str, strlen(str));
                alig-> alignmentInfo.append(" ", strlen(" "));
                str = strtok(nullptr, DELIMITERS);
            }
        }

            /* In this case, try to get matrix dimensions */
        else if((!strcmp(str, "DIMENSIONS")) && state) {
            str = strtok(nullptr, DELIMITERS);
            frag = strtok(nullptr, DELIMITERS);
            str = strtok(str, "=;");
            alig->numberOfSequences = atoi(strtok(nullptr, "=;"));
            frag = strtok(frag, "=;");
            alig->numberOfResidues = atoi(strtok(nullptr, "=;"));
        }
    } while(!file.eof());

    /* Check all parameters */
    if(strcmp(str, "MATRIX") || (alig->numberOfSequences == 0) || (alig->numberOfResidues == 0))
        return nullptr;

    /* Allocate memory for the input alignment */
    alig->seqsName  = new std::string[alig->numberOfSequences];
    alig->sequences = new std::string[alig->numberOfSequences];

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
            alig->seqsName[pos].append(str, strlen(str));

        /* Store rest of line as part of sequence */
        str = strtok(nullptr, OTH2DELIMITERS);
        while(str != nullptr) {
            alig->sequences[pos].append(str, strlen(str));
            str = strtok(nullptr, OTH2DELIMITERS);
        }

        /* Move sequences pointer to next one. It if it is last one, move it to
         * the beginning and set the first block to false for avoiding to rewrite
         * sequences name */
        pos = (pos + 1) % alig->numberOfSequences;
        if (not pos)
            firstBlock = false;
    }

    /* Deallocate memory */
    delete [] line;

    /* Close the input file */
    file.close();

    /* Check the matrix's content */
    alig->fillMatrices(true);
    alig->originalNumberOfSequences = alig-> numberOfSequences;
    alig->originalNumberOfResidues =alig->numberOfResidues;
    return alig;
}

bool nexus_state::SaveAlignment(const Alignment &alignment, std::ostream *output)
{
    /* Generate output alignment in NEXUS format setting only alignment block */

    int i, j, k, l, maxLongName = 0;
    std::string *tmpMatrix;

    // Check whether sequences in the alignment are aligned or not.
    // Warn about it if there are not aligned.
    if (!alignment.isAligned) {
        debug.report(ErrorCode::UnalignedAlignmentToAlignedFormat, new std::string[1] { this->name });
        return false;
    }

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

    // Compute maximum sequences name length
    for (i = 0; (i < alignment.originalNumberOfSequences); i++)
        if (alignment.saveSequences[i] != -1)
            maxLongName = utils::max(maxLongName, alignment.seqsName[i].size());

    // Compute output file datatype
    alignment.getAlignmentType();

    std::string alignmentInfo {alignment.alignmentInfo};
    // Remove characters like ";" from input alignment information line
    while((int) alignment. alignmentInfo.find(';') != (int) std::string::npos)
        alignmentInfo.erase(alignment. alignmentInfo.find(';'), 1);

    // Print Alignment header
    *output << "#NEXUS\nBEGIN DATA;\n DIMENSIONS NTAX="
         << alignment.numberOfSequences << " NCHAR=" << alignment.numberOfResidues <<";\n";

    // Print alignment datatype
    if (alignment.getAlignmentType() & SequenceTypes::DNA)
        *output << "FORMAT DATATYPE=DNA INTERLEAVE=yes GAP=-";
    else if (alignment.getAlignmentType() & SequenceTypes::RNA)
        *output << "FORMAT DATATYPE=RNA INTERLEAVE=yes GAP=-";
    else if (alignment.getAlignmentType() & SequenceTypes::AA)
        *output << "FORMAT DATATYPE=PROTEIN INTERLEAVE=yes GAP=-";

    i = 0;
    // Using information from input alignment. Use only some tags.
    while((j = alignmentInfo.find(' ', i)) != (int) std::string::npos) {

        if((alignmentInfo.substr(i, j - i)).compare(0, 7, "MISSING") == 0 ||
           (alignmentInfo.substr(i, j)).compare(0, 7, "missing") == 0)
            *output << " " << (alignmentInfo.substr(i, j - i));

        else if((alignmentInfo.substr(i, j)).compare(0, 9, "MATCHCHAR") == 0 ||
                (alignmentInfo.substr(i, j)).compare(0, 9, "matchchar") == 0)
            *output << " " << (alignmentInfo.substr(i, j - i));

        i = j + 1;
    }
    *output << ";\n";

    // Add a header indicating the number of residues of each sequence.
    for(i = 0; i < alignment.originalNumberOfSequences; i++)
    {
        if (alignment.saveSequences[i] == -1 )
            continue;

        *output << "[Name: "
                << std::setw(maxLongName + 4) << std::left << alignment.seqsName[i]
                << "Len: " << alignment.numberOfResidues << "]\n";
    }
    *output << "\nMATRIX";

    // Start filling the file with sequence names and sequences.
    for (j = 0, k = 0; j < alignment.originalNumberOfResidues;)
    {
        // Move until next not rejected residue
        if (alignment.saveResidues[j] == -1)
        {
            j++;
            continue;
        }

        // Iterate over the sequences
        for (i = 0; i < alignment.originalNumberOfSequences; i++)
        {
            // Skip rejected sequences
            if (alignment.saveSequences[i] == -1) continue;
            // Output sequence name
            *output << "\n" << std::setw(maxLongName + 5) << std::left << alignment.seqsName[i];
            // Add residues per block.
            // k = residue position;
            // l = residues added on current line
            for (k = j, l = 0; k < alignment.originalNumberOfResidues && l < 50;)
            {
                // Don't save residues that have been rejected
                if (alignment.saveResidues[k] == -1)
                {
                    k++;
                    continue;
                }
                // Add the residue at position K
                *output << alignment.sequences[i][k];
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
    if (Machine->reverse)
        delete [] tmpMatrix;
    
    return true;
}

bool nexus_state::RecognizeOutputFormat(const std::string &FormatName)
{
    if (BaseFormatHandler::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "nexus";
}

}