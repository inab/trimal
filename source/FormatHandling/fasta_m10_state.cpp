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

#include "FormatHandling/fasta_m10_state.h"

#include "FormatHandling/FormatManager.h"
#include "utils.h"

namespace FormatHandling {
int fasta_m10_state::CheckAlignment(std::istream* origin)
{
    return 0;
}

Alignment* fasta_m10_state::LoadAlignment(const std::string &filename)
{
    return nullptr;
}

bool fasta_m10_state::SaveAlignment(const Alignment &alignment, std::ostream *output)
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
        if (alignment.saveSequences && alignment.saveSequences[i] == -1) continue;
        if (!Machine->keepHeader)
            maxLongName = utils::max(maxLongName, alignment.seqsName[i].size());
        else if (alignment.seqsInfo != nullptr)
            maxLongName = utils::max(maxLongName, alignment.seqsInfo[i].size());
    }

    if (maxLongName > PHYLIPDISTANCE) {
        maxLongName = PHYLIPDISTANCE;
        debug.report(WarningCode::HeaderWillBeCut, new std::string[1]{std::string(name)});
    }
    /* Print alignment. First, sequences name id and then the sequences itself */
    for(i = 0; i < alignment.originalNumberOfSequences; i++) {
        
        if (alignment.saveSequences != nullptr && alignment.saveSequences[i] == -1) continue;
        
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
                if ((k % 60 == 0 || j == alignment.sequences[i].length() - 1))
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

bool fasta_m10_state::RecognizeOutputFormat(const std::string &FormatName)
{
    if (BaseFormatHandler::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "fasta_m10";
}
}