#include "FormatHandling/fasta_m10_state.h"

#include "FormatHandling/FormatManager.h"
#include "defines.h"
#include "utils.h"

namespace FormatHandling {
int fasta_m10_state::CheckAlignment(std::istream* origin)
{
    return 0;
}

Alignment* fasta_m10_state::LoadAlignment(std::string& filename)
{
    return nullptr;
}

bool fasta_m10_state::SaveAlignment(Alignment* alignment, std::ostream* output)
{
    /* Generate output alignment in FASTA format. Sequences can be unaligned. */

    int i, j, k, maxLongName;
    std::string *tmpMatrix;
    bool lastcharIsnewline = false;

    /* Allocate local memory for generating output alignment */
    tmpMatrix = new std::string[alignment->originalNumberOfSequences];

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    for(i = 0; i < alignment->originalNumberOfSequences; i++)
        tmpMatrix[i] = (!Machine->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    /* Depending on if short name flag is activated (limits sequence name up to
     * 10 characters) or not, get maximum sequence name length. Consider those
     * cases when the user has asked to keep original sequence header */
    maxLongName = 0;
    for(i = 0; i < alignment->originalNumberOfSequences; i++)
    {
        if (alignment->saveSequences && alignment->saveSequences[i] == -1) continue;
        if (!Machine->keepHeader)
            maxLongName = utils::max(maxLongName, alignment->seqsName[i].size());
        else if (alignment->seqsInfo != nullptr)
            maxLongName = utils::max(maxLongName, alignment->seqsInfo[i].size());
    }

    if (maxLongName > PHYLIPDISTANCE) {
        maxLongName = PHYLIPDISTANCE;
        if (!Machine->keepHeader)
            debug.report(WarningCode::HeaderWillBeCut);
    }
    /* Print alignment. First, sequences name id and then the sequences itself */
    for(i = 0; i < alignment->originalNumberOfSequences; i++) {
        
        if (alignment->saveSequences != nullptr && alignment->saveSequences[i] == -1) continue;
        
        if (!Machine->keepHeader)
            (*output) << ">" << alignment->seqsName[i].substr(0, maxLongName) << "\n";
        
        else if (alignment->seqsInfo != nullptr)
            (*output) << ">" << alignment->seqsInfo[i].substr(0, maxLongName) << "\n";
        
        
        for (j = 0, k = 0; j < alignment->sequences[i].length(); j++)
        {
            if (alignment->saveResidues != nullptr && alignment->saveResidues[j] == -1) 
            {
                if (!lastcharIsnewline && j == alignment->sequences[i].length() -1 ) 
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
                if (!lastcharIsnewline && (k % 60 == 0 || j == alignment->sequences[i].length() -1 ))
                {
                    (*output) << "\n";
                    lastcharIsnewline = true;
                }
            }
        }
    }

    /* Deallocate local memory */
    delete [] tmpMatrix;
    
    return true;
}

bool fasta_m10_state::RecognizeOutputFormat(std::string& FormatName)
{
    if (BaseFormatHandler::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "fasta_m10";
}
}