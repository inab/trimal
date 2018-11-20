#include "../../include/ReadWriteMS/phylip32_m10_state.h"

#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/defines.h"
#include "../../include/utils.h"

int phylip32_m10_state::CheckAlignment(std::istream* origin)
{
    return 0;
}

newAlignment* phylip32_m10_state::LoadAlignment(std::string filename)
{
    return nullptr;
}

bool phylip32_m10_state::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
    /* Generate output alignment in PHYLIP 3.2 format (interleaved) */

    int i, j, k, maxLongName;
    std::string *tmpMatrix;

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!alignment->isAligned) {
        debug.report(ErrorCode::UnalignedAlignmentToAlignedFormat, new std::string[1] { this->name });
        return false;
    }

    /* Allocate local memory for generating output alignment */
    tmpMatrix = new std::string[alignment->originalSequenNumber];

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    for(i = 0; i < alignment->originalSequenNumber; i++)
        tmpMatrix[i] = (!Machine->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    /* Depending on if short name flag is activated (limits sequence name up to
     * 10 characters) or not, get maximum sequence name length */
    maxLongName = PHYLIPDISTANCE;
    for(i = 0; (i < alignment->originalSequenNumber); i++)
        if  (alignment->saveSequences[i] != -1)
            maxLongName = utils::max(maxLongName, alignment->seqsName[i].size());

    maxLongName = std::min(maxLongName, PHYLIPDISTANCE);

    /* Generating output alignment */
    /* First Line: Sequences Number & Residued Number */
    (*output) << " " << alignment->sequenNumber << " " << alignment->residNumber;

    /* Alignment */
    /* For each sequence, print its identifier and then the sequence itself in
     * blocks of 50 residues */
    for(i = 0; i < alignment->originalSequenNumber; i++) {
        /* Sequence Name */
        if (alignment->saveSequences[i] == -1) continue;
        (*output) << "\n" << std::setw(maxLongName + 3) << std::left << alignment->seqsName[i].substr(0, maxLongName);
        /* Sequence. Each line contains a block of 5 times 10 residues. */
        
        for (j = 0, k = 0; j < alignment->originalResidNumber; j++)
        {
            if (alignment->saveResidues[j] == -1) continue;
            if (k == 50)
            {
                *output << "\n" << std::setw(maxLongName + 3) << std::left << " " ;
                k = 0;
            }
            *output << alignment->sequences[i][j];
            k++;
            if (k % 10 == 0) 
                *output << " ";
        }
        if (k % 10 != 0)
            *output << " ";
        
        /* Print a blank line to mark sequences separation */
//         if (k % 50 != 0)
        (*output) << "\n";
    }
    *output << "\n";

    /* Deallocate local memory */
    delete [] tmpMatrix;
    
    return true;
}

bool phylip32_m10_state::RecognizeOutputFormat(std::string FormatName)
{
    if (ReadWriteBaseState::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "phylip32_m10";
}

