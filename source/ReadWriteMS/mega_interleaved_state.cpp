#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/ReadWriteMS/mega_interleaved_state.h"
#include "../../include/defines.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

int MegaInterleavedState::CheckAlignment(istream* origin)
{
    string str;
     origin->seekg(0);
     getline( (*origin), str );
     str = str.substr(0, 5);
     if (strcmp(str.c_str(), "#MEGA"))
         return 0;
    
     
     
    return 0;
}

newAlignment* MegaInterleavedState::LoadAlignment(istream* origin)
{
    // TODO Not Implemented
    origin->clear();
    origin->seekg(0);

    newAlignment* _alignment = new newAlignment();
    
    return _alignment;
}

void MegaInterleavedState::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
    (*output) << "#MEGA" << endl << "!Title " << *FileName << ";" << endl << "!Format DataType=";
    
    alignment->getAlignmentType();
        
    if (alignment->dataType == DNAType)
        (*output) << "dna";
    else if (alignment->dataType == RNAType)
        (*output) << "rna";
    else if (alignment->dataType == AAType)
        (*output) << "protein";
    else if (alignment->dataType == DNADeg)
        (*output) << "dna";
    else if (alignment->dataType == RNADeg)
        (*output) << "rna";
    
    (*output) << " NSeqs=" << alignment->sequenNumber << " Nsites=" << alignment->residNumber << " indel=- CodeTable=Standard" << endl << endl;
    
    int residCount = 0, y;
    bool continueTag = true;
    
    int max = 0;
    
    for (int x = 0; x < alignment->sequenNumber; x++)
    {
        max = std::max(max, (int)alignment->seqsName[x].length());
    }
    
    max += 1;
    
    while (continueTag)
    {
        continueTag = false;
        
        for (int x = 0; x < alignment->sequenNumber; x++)
        {
            (*output) << "#" << setw(max) << left << alignment->seqsName[x] ;
            
            for (y = residCount; y < std::min((int)alignment->sequences[x].length(), residCount + 60); y+=10)
            {
                (*output) << alignment->sequences[x].substr(y, 10) << " ";
            }
            
            (*output) << endl;
            
            if (residCount + 60 < alignment->sequences[x].length())
                continueTag = true;
        }
        
        (*output) << endl;
        residCount += 60;
    }
}

bool MegaInterleavedState::RecognizeOutputFormat(std::string FormatName)
{
    if (FormatName == "megai") return true;
    return false;
}

