#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/ReadWriteMS/phylip40_m10_state.h"
#include "../../include/defines.h"
#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include "../../include/newAlignment.h"

using namespace std;

int phylip40_m10_state::CheckAlignment(istream* origin)
{
    return 0;
}

newAlignment* phylip40_m10_state::LoadAlignment(std::string filename)
{
    return nullptr;
}

bool phylip40_m10_state::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
  
    
   /* Generate output alignment in PHYLIP/PHYLIP 4 format (sequential) */

    int i, j, k, l, maxLongName;
    string *tmpMatrix;

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!alignment->isAligned) {
        debug.report(ErrorCode::UnalignedAlignmentToAlignedFormat, new std::string[1] { this->name });
        return false;
    }

    /* Allocate local memory for generating output alignment */
    tmpMatrix = new string[alignment->originalSequenNumber];

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
        maxLongName = utils::max(maxLongName, alignment->seqsName[i].size());

    maxLongName = std::min(maxLongName, PHYLIPDISTANCE);

    /* Generating output alignment */
    /* First Line: Sequences Number & Residued Number */
    (*output) << " " << alignment->sequenNumber << " " << alignment->residNumber;

    /* First Block: Sequences Names & First 60 residues */
    for(i = 0; i < alignment->originalSequenNumber; i++)
    {
        if (alignment->saveSequences[i] == -1) continue;
        (*output) << endl << setw(maxLongName + 3) << left << alignment->seqsName[i].substr(0, maxLongName);
            
        for (k = 0, l = 0; k < alignment->originalResidNumber && l < 60; k++)
        {
            if (alignment->saveResidues[k] == -1) continue;
            *output << alignment->sequences[i][k];
            l++;
        }
    }


    for (i = k; i < alignment->originalResidNumber; i=k)
    {
        if (alignment->saveResidues[i] == -1) {
            k++;
            continue;
        }
        *output << endl;
        for (j = 0; j < alignment->originalSequenNumber; j++)
        {
            if (alignment->saveSequences[j] == -1) continue;
            *output << endl;
            for (k = i, l = 0; k < alignment->originalResidNumber && l < 60; k++)
            {
                if (alignment->saveResidues[k] == -1) continue;
                *output << alignment->sequences[j][k];
                l++;
            }
            
        }
    }
    
    *output << endl << endl << endl;

    /* Deallocate local memory */
    delete [] tmpMatrix;
    
    return true;
}

bool phylip40_m10_state::RecognizeOutputFormat(std::string FormatName)
{
    if (ReadWriteBaseState::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "phylip_m10" || FormatName == "phylip40_m10";
}

