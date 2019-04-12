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

#include "FormatHandling/phylip40_m10_state.h"

#include "FormatHandling/FormatManager.h"
#include "utils.h"

namespace FormatHandling {
int phylip40_m10_state::CheckAlignment(std::istream* origin)
{
    return 0;
}

Alignment* phylip40_m10_state::LoadAlignment(const std::string &filename)
{
    return nullptr;
}

bool phylip40_m10_state::SaveAlignment(const Alignment &alignment, std::ostream *output)
{
  
    
   /* Generate output alignment in PHYLIP/PHYLIP 4 format (sequential) */

    int i, j, k, l, maxLongName;
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

    /* Depending on if short name flag is activated (limits sequence name up to
     * 10 characters) or not, get maximum sequence name length */
    maxLongName = PHYLIPDISTANCE;
    for(i = 0; (i < alignment.originalNumberOfSequences); i++)
        maxLongName = utils::max(maxLongName, alignment.seqsName[i].size());

    if (maxLongName > PHYLIPDISTANCE) {
        maxLongName = PHYLIPDISTANCE;
        debug.report(WarningCode::HeaderWillBeCut, new std::string[1]{std::string(name)});
    }

    /* Generating output alignment */
    /* First Line: Sequences Number & Residued Number */
    (*output) << " " << alignment.numberOfSequences << " " << alignment.numberOfResidues;

    /* First Block: Sequences Names & First 60 residues */
    for(i = 0; i < alignment.originalNumberOfSequences; i++)
    {
        if (alignment.saveSequences[i] == -1) continue;
        (*output) << "\n" << std::setw(maxLongName + 3) << std::left << alignment.seqsName[i].substr(0, maxLongName);
            
        for (k = 0, l = 0; k < alignment.originalNumberOfResidues && l < 60; k++)
        {
            if (alignment.saveResidues[k] == -1) continue;
            *output << alignment.sequences[i][k];
            l++;
        }
    }


    for (i = k; i < alignment.originalNumberOfResidues; i=k)
    {
        if (alignment.saveResidues[i] == -1) {
            k++;
            continue;
        }
        *output << "\n";
        for (j = 0; j < alignment.originalNumberOfSequences; j++)
        {
            if (alignment.saveSequences[j] == -1) continue;
            *output << "\n";
            for (k = i, l = 0; k < alignment.originalNumberOfResidues && l < 60; k++)
            {
                if (alignment.saveResidues[k] == -1) continue;
                *output << alignment.sequences[j][k];
                l++;
            }
            
        }
    }
    
    *output << "\n\n\n";

    /* Deallocate local memory */
    if (Machine->reverse)
        delete [] tmpMatrix;
    
    return true;
}

bool phylip40_m10_state::RecognizeOutputFormat(const std::string &FormatName)
{
    if (BaseFormatHandler::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "phylip_m10" || FormatName == "phylip40_m10";
}

}