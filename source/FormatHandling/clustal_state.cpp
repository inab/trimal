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

#include "FormatHandling/clustal_state.h"

#include "FormatHandling/FormatManager.h"
#include "defines.h"
#include "utils.h"

namespace FormatHandling {
int clustal_state::CheckAlignment(std::istream* origin)
{
    origin->seekg(0);
    origin->clear();
    char *firstWord = nullptr, *line = nullptr;

    
    /* Read first valid line in a safer way */
    do {
        delete[] line;
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
    if((!strcmp(firstWord, "CLUSTAL")) || (!strcmp(firstWord, "clustal")))
    {
        delete [] line;
        return 1;
    }

    delete[] line;
    
    return 0;
}

Alignment* clustal_state::LoadAlignment(const std::string &filename)
{
    Alignment* alignment = new Alignment();
    int i, seqLength, pos, firstBlock;
    char *str, *line = nullptr;
    std::ifstream file;
    file.open(filename, std::ifstream::in);
    
    /* Store some details about input file to be used in posterior format
     * conversions */
    // alignment.filename.append("!Title ");
    alignment->filename.append(filename);
    alignment->filename.append(";");

    /* The first valid line corresponding to CLUSTAL label is ignored */
    do {
        delete [] line;
        line = utils::readLine(file);
    } while ((line == nullptr) && (!file.eof()));

    /* If the file end is reached without a valid line, warn about it */
    if (file.eof())
        return nullptr;

    /* Ignore blank lines before first sequence block starts */
    while(!file.eof()) {

        /* Deallocate previously used dynamic memory */
        delete [] line;

        /* Read lines in safe way */
        line = utils::readLine(file);

        if (line != nullptr)
            break;
    }

    /* The program in only interested in the first blocks of sequences since
     * it wants to know how many sequences are in the input file */
    alignment->numberOfSequences = 0;
    while(!file.eof()) {

        /* If a new line without any valid character is detected
         * means the first block is over */
        if (line == nullptr)
            break;

        /* Count how many times standard characters as well as
         * gap symbol "-" is detected in current line. */
        seqLength = (int) strlen(line);
        for(pos = 0; pos < seqLength; pos++)
            if((isalpha(line[pos])) || (line[pos] == '-'))
                break;

        /* If not standard characters are detected in current line means that
         * the program has found typical line in clustal alignment files to mark
         * some scores for some columns. In that case, the first block is over */
        if(pos == seqLength)
            break;
        alignment->numberOfSequences++;

        /* Deallocate previously used dynamic memory */
        delete [] line;

        /* Read lines in safe way */
        line = utils::readLine(file);
    }

    delete [] line;

    /* Finish to preprocess the input file. */
    file.clear();
    file.seekg(0);

    /* Allocate memory for the input alignmet */
    alignment->seqsName  = new std::string[alignment->numberOfSequences];
    alignment->sequences = new std::string[alignment->numberOfSequences];

    /* Read the title line and store it */
    line = utils::readLine(file);
    if (line == nullptr)
        return nullptr;
    alignment->alignmentInfo.append(line, strlen(line));

    /* Ignore blank lines before first sequence block starts */
    while(!file.eof()) {

        /* Deallocate previously used dynamic memory */
        delete [] line;

        /* Read lines in safe way */
        line = utils::readLine(file);

        if (line != nullptr)
            break;
    }

    /* Set-up sequences pointer to the first one and the flag to indicate
     * the first blocks. That flag implies that sequences names have to be
     * stored */
    i = 0;
    firstBlock = true;

    while(!file.eof()) {

        if (line == nullptr) {
            /* Sometimes, clustalw files does not have any marker after first block
             * to indicate similarity between its columns residues. In that cases,
             * mark the end of first block */
            if (i == 0)
                firstBlock = false;
            /* Read current line and analyze it*/
            line = utils::readLine(file);
            continue;
        }

        /* Check whteher current line is a standard line or it is a line to mark
         * quality scores for that alignment column */
        seqLength = (int) strlen(line);
        for(pos = 0; pos < seqLength; pos++)
            if((isalpha(line[pos])) || (line[pos] == '-'))
                break;

        /* Start a new block in the input alignment */
        if (pos == seqLength) {
            firstBlock = false;

            /* Deallocate dinamic memory if it has been used before */
            delete [] line;

            /* Read current line and analyze it*/
            line = utils::readLine(file);

            continue;
        }

        /* If it is a standard line, split it into two parts. The first one contains
         * sequence name and the second one the residues. If the "firstBlock" flag
         * is active then store the sequence name */
        str = strtok(line, OTHDELIMITERS);
        if(str != nullptr) {
            if(firstBlock)
                alignment->seqsName[i].append(str, strlen(str));
            str = strtok(nullptr, OTHDELIMITERS);
            if(str != nullptr)
                alignment->sequences[i].append(str, strlen(str));

            /* Move sequences pointer in a circular way */
            i = (i + 1) % alignment->numberOfSequences;
        }

        /* Deallocate dinamic memory if it has been used before */
        delete [] line;

        /* Read current line and analyze it*/
        line = utils::readLine(file);
    }

    /* Close the input file */
    file.close();

    /* Deallocate dinamic memory */
    delete [] line;

    /* Check the matrix's content */
    alignment->fillMatrices(true);
    alignment->originalNumberOfSequences = alignment->numberOfSequences;
    alignment->originalNumberOfResidues = alignment->numberOfResidues;
    return alignment; 
}

bool clustal_state::SaveAlignment(const Alignment &alignment, std::ostream *output)
{
    /* Generate output alignment in CLUSTAL format */

    int i, j, k, l, maxLongName = 0;
    std::string *tmpMatrix;

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
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

    /* Print alignment header */
    if(!alignment.alignmentInfo.empty() && alignment.alignmentInfo.substr(0,7) == "CLUSTAL")
        (*output) << alignment.alignmentInfo << "\n\n";
    else
        (*output) << "CLUSTAL multiple sequence alignment\n\n";

    /* Print alignment itself */
    /* Print as many blocks as it is needed of lines composed
     * by sequences name and 60 residues */
    
    for(j = 0, k = 0; j < alignment.originalNumberOfResidues; j = k)
    {
        for(i = 0; i < alignment.originalNumberOfSequences; i++)
        {
            if (alignment.saveSequences != nullptr && alignment.saveSequences[i] != -1)
            {
                (*output) << std::setw(maxLongName + 5) << std::left << alignment.seqsName[i];
                for (k = j, l = 0; l < 60 && k < alignment.originalNumberOfResidues; k++)
                {
                    if (alignment.saveResidues != nullptr && alignment.saveResidues[k] != -1)
                    {
                        (*output) << tmpMatrix[i][k];
                        l++;
                    }
                }
                (*output) << "\n";
            }
        }
        (*output) << "\n\n";
    }

    /* Deallocate local memory */
    if (Machine->reverse)
        delete [] tmpMatrix;
    
    return true;
}

bool clustal_state::RecognizeOutputFormat(const std::string &FormatName)
{
    if (BaseFormatHandler::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "clustal";
}

}