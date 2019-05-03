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

#include "FormatHandling/fasta_state.h"

#include "FormatHandling/FormatManager.h"
#include "defines.h"
#include "utils.h"

namespace FormatHandling {
int fasta_state::CheckAlignment(std::istream* origin)
{
    char c;
    origin->seekg(0);
    origin->get(c);
    if (c == '>')
        return 1;
    return 0;
}

Alignment* fasta_state::LoadAlignment(const std::string &filename)
{
    /* FASTA file format parser */
    Alignment* alig = new Alignment();
    char *str, *line = nullptr;
    std::ifstream file;
    int i;

    /* Check the file and its content */
    file.open(filename, std::ifstream::in);
    if(!utils::checkFile(file))
        return nullptr;

    /* Store input file name for posterior uses in other formats */
    // alig->filename.append("!Title ");
    alig->filename.append(filename);
    alig->filename.append(";");

    /* Compute how many sequences are in the input alignment */
    alig->numberOfSequences = 0;
    while(!file.eof()) {

        /* Deallocate previously used dinamic memory */
        delete [] line;

        /* Read lines in a safe way */
        line = utils::readLine(file);
        if (line == nullptr)
            continue;

        /* It the line starts by ">" means that a new sequence has been found */
        str = strtok(line, DELIMITERS);
        if (str == nullptr)
            continue;

        /* If a sequence name flag is detected, increase sequences counter */
        if(str[0] == '>')
            alig->numberOfSequences++;
    }

    /* Finish to preprocess the input file. */
    file.clear();
    file.seekg(0);

    /* Allocate memory for the input alignmet */
    alig->seqsName  = new std::string[alig->numberOfSequences];
    alig->sequences = new std::string[alig->numberOfSequences];
    alig->seqsInfo  = nullptr;//new std::string[alig->numberOfSequences];

    for(i = -1; (i < alig->numberOfSequences) && (!file.eof()); ) {

        /* Deallocate previously used dinamic memory */
        delete [] line;

        /* Read lines in a safe way */
        line = utils::readLine(file);
        if (line == nullptr)
            continue;

        /* Store original header fom input sequences including non-standard
         * characters */
//         if (line[0] == '>')
//             alig->seqsInfo[i+1].append(&line[1], strlen(line) - 1);

        /* Cut the current line and check whether there are valid characters */
        str = strtok(line, OTHDELIMITERS);
        if (str == nullptr)
            continue;

        /* Check whether current line belongs to the current sequence
         * or it is a new one. In that case, store the sequence name */
        if(str[0] == '>') {
            /* Move sequence name pointer until a valid std::string name is obtained */
            do {
                str = str + 1;
            } while(strlen(str) == 0);
            alig->seqsName[++i].append(str, strlen(str));
            continue;
        }

        /* Sequence */
        while(str != nullptr) {
            alig->sequences[i].append(str, strlen(str));
            str = strtok(nullptr, DELIMITERS);
        }
    }

    /* Close the input file */
    file.close();

    /* Deallocate previously used dinamic memory */
    if (line != nullptr)
        delete [] line;
        
    /* Check the matrix's content */
    alig->fillMatrices(false);
    alig->originalNumberOfSequences = alig-> numberOfSequences;
    alig->originalNumberOfResidues = alig->numberOfResidues;
    return alig; 
}

bool fasta_state::SaveAlignment(const Alignment &alignment,
                                std::ostream *output)
{
    /* Generate output alignment in FASTA format. Sequences can be unaligned. */

    int i, j, k, maxLongName;
    std::string *tmpMatrix;
    bool lastcharIsnewline = true;

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

    /* Depending on if short name flag is activated (limits sequence name up to
     * 10 characters) or not, get maximum sequence name length. Consider those
     * cases when the user has asked to keep original sequence header */
    maxLongName = 0;
    for(i = 0; i < alignment.originalNumberOfSequences; i++)
    {
        if (alignment.saveSequences && alignment.saveSequences[i] == -1)
            continue;
        if (!Machine->keepHeader)
            maxLongName = utils::max(maxLongName, alignment.seqsName[i].size());
        else if (alignment.seqsInfo != nullptr)
            maxLongName = utils::max(maxLongName, alignment.seqsInfo[i].size());
    }

    /* Print alignment. First, sequences name id and then the sequences itself */
    for(i = 0; i < alignment.originalNumberOfSequences; i++) {
        
        if (alignment.saveSequences != nullptr &&
            alignment.saveSequences[i] == -1)
            continue;
        
        if (!Machine->keepHeader)
            (*output) << ">" << alignment.seqsName[i].substr(0, maxLongName) << "\n";
        
        else if (alignment.seqsInfo != nullptr)
            (*output) << ">" << alignment.seqsInfo[i].substr(0, maxLongName) << "\n";
        
        
        for (j = 0, k = 0; j < alignment.sequences[i].length(); j++)
        {
            if (alignment.saveResidues != nullptr && alignment.saveResidues[j] == -1) 
            {
                if (!lastcharIsnewline && j == alignment.sequences[i].length() -1 ) 
                {
                    (*output) << "\n";
                    lastcharIsnewline = true;
                }
            }
            else
            {
                (*output) << tmpMatrix[i][j];
                k++;
                lastcharIsnewline = false;
                if (k % 60 == 0 || j == alignment.sequences[i].length() -1 )
                {
                    (*output) << "\n";
                    lastcharIsnewline = true;
                }
            }
        }
    }

    /* Deallocate local memory */
    if (Machine->reverse)
        delete [] tmpMatrix;
    
    return true;
}

bool fasta_state::RecognizeOutputFormat(const std::string &FormatName)
{
    if (BaseFormatHandler::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "fasta";
}
}