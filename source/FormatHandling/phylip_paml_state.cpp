#include "../../include/FormatHandling/phylip_paml_state.h"

#include "../../include/FormatHandling/FormatManager.h"
#include "../../include/defines.h"
#include "../../include/utils.h"

int phylip_paml_state::CheckAlignment(std::istream* origin)
{
    return 0;
}

Alignment* phylip_paml_state::LoadAlignment(std::string& filename)
{
    return nullptr;
}

bool phylip_paml_state::SaveAlignment(Alignment* alignment, std::ostream* output, std::string* FileName)
{
    /* Generate output alignment in PHYLIP format compatible with PAML program */

    int i, maxLongName;
    std::string *tmpMatrix;

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!alignment->isAligned) {
        debug.report(ErrorCode::UnalignedAlignmentToAlignedFormat, new std::string[1] { this->name });
        return false;
    }

    /* Allocate local memory for generating output alignment */
    tmpMatrix = new std::string[alignment->numberOfSequences];

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    for(i = 0; i < alignment->numberOfSequences; i++)
        tmpMatrix[i] = (!Machine->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    maxLongName = PHYLIPDISTANCE;
    for(i = 0; (i < alignment->numberOfSequences); i++)
        maxLongName = utils::max(maxLongName, alignment->seqsName[i].size());

    /* Generating output alignment */
    /* First Line: Sequences Number & Residued Number */
    *output << " " << alignment->numberOfSequences << " " << alignment->numberOfResidues << "\n";

    /* Print alignment */
    /* Print sequences name follow by the sequence itself in the same line */
    for(i = 0; i < alignment->numberOfSequences; i++)
        *output << std::setw(maxLongName + 3) << std::left << alignment->seqsName[i].substr(0, maxLongName)
             << alignment->sequences[i] << "\n";
    *output << "\n";

    /* Deallocate local memory */
    delete [] tmpMatrix;
    
    return true;
}

bool phylip_paml_state::RecognizeOutputFormat(std::string& FormatName)
{
    if (BaseFormatHandler::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "phylippaml";
}

