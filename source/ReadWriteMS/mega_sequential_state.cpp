#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/ReadWriteMS/mega_sequential_state.h"
#include "../../include/defines.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

int MegaSequentialState::CheckAlignment(istream* origin)
{
    string str;
     origin->seekg(0);
     getline( (*origin), str );
     str = str.substr(0, 5);
     if (!strcmp(str.c_str(), "#MEGA"))
         return 1;
 
    return 0;
}

newAlignment* MegaSequentialState::LoadAlignment(istream* origin)
{
    origin->clear();
    origin->seekg(0);

    newAlignment* _alignment = new newAlignment();
    _alignment->sequenNumber = 0;
    
    string line;
    char * temp;
    getline(*origin , line);
    
    // Read the sequences number from Format Comentary.
    for (;getline(*origin, line); )
    {
        if (line[0] == '#') break;
        if (line.substr(0, 7) == "!Format")
        {
            temp = strtok(&line[0], " ");
            
            while (temp != NULL)
            {
                temp = strtok(nullptr, "=");
                if (temp != nullptr)
                {
                    if (!strcmp(temp, "NSeqs"))
                    {
                        _alignment->sequenNumber = stoi(strtok(nullptr, " "));
                        continue;
                    }
                    temp = strtok(nullptr, " ");
                }
            }
        }
    }
    // Calculate the number of sequences from file.
    if (_alignment->sequenNumber == 0)
    {
        if (line[0] == '#') _alignment->sequenNumber++;
        
        for (;getline(*origin, line); )
            if (line[0] == '#') _alignment->sequenNumber++;
            
        origin->clear();
        origin->seekg(0);
        
        for (;getline(*origin, line); )
        {
            if (line[0] == '#') break;
        }
    }
    
    _alignment->seqsName = new string[_alignment->sequenNumber];
    _alignment->sequences = new string[_alignment->sequenNumber];
    
    int x = -1;
    while (!origin->eof())
    {
        if (line[0] == '#')
        {
            _alignment->seqsName[++x] = line.substr(1);
            _alignment->sequences[x] = "";
        }
        else
        {
            for (char c : line)
            {
                if (c == ' ') continue;
                _alignment->sequences[x].append(1, c);
            }
        }
        getline(*origin, line);
    }
    _alignment->fillMatrices(false);
    
    return _alignment;
}

void MegaSequentialState::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
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
    
    for (int i = 0; i < alignment->sequenNumber; i++)
    {
        if (Machine->shortNames)
            (*output) << "#" << alignment->seqsName[i].substr(0, std::min(10, (int)alignment->seqsName[i].length()));
        else
            (*output) << "#" << alignment->seqsName[i];
        
        
        for (int x = 0; x < alignment->sequences[i].length(); x+= 10)
        {
            if (x % 50 == 0) 
                (*output) << endl;
            
            (*output) << alignment->sequences[i].substr(x, 10) << " ";
        }
        
        (*output) << endl ;
    }
}

bool MegaSequentialState::RecognizeOutputFormat(std::string FormatName)
{
    if (FormatName == "meganon") return true;
    return false;
}

