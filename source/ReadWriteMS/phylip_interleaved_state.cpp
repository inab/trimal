#include "../../include/ReadWriteMS/phylip_interleaved_state.h"
#include <iostream>
#include "../../include/defines.h"
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

int PhylipInterleavedState::CheckAlignment(istream* origin)
{
    string str;
    origin->seekg(0);
    getline( (*origin), str );
    
    //Check if first line is composed by two numbers
    string str2 = strtok(&str[0], " \t");
    // First number
    for (int i = 0; i < str2.length(); i++)
    {
        if (isdigit(str2[i]))
            continue;
        else return 0;
    }
    
    str2 = strtok(nullptr, " \t");
    // Second number
    for (int i = 0; i < str2.length(); i++)
    {
        if (isdigit(str2[i]))
            continue;
        else return 0;
    }
    
    // Check if the first two lines have the same length.
    // If not, it may be the sequential or relaxed format.
    getline( (*origin), str );
    str = str.substr(str.find_first_not_of(" "), str.find_last_not_of(" \t\n"));
    int size = str.length();
    getline( (*origin), str );
    str = str.substr(str.find_first_not_of(" "), str.find_last_not_of(" \t\n"));
    
    // Check interleaved or sequential
    if (size != str.length())
        return 0;
    
//      cout << "Input Format: Phylip interleaved" << endl;
    return 1;
}

newAlignment* PhylipInterleavedState::LoadAlignment(istream* origin)
{
    origin->clear();
    origin->seekg(0);

    newAlignment* _alignment = new newAlignment();
    
    std::string line;
    // Get rid of the first line.
    getline( (*origin), line );

    
    /* Compute how many sequences are in the input alignment */
    _alignment->sequenNumber = stoi(strtok(&line[0u], " "));
    /* Allocate memory for the input alignmet */
    _alignment->seqsName  = new string[_alignment->sequenNumber];
    _alignment->sequences = new string[_alignment->sequenNumber];
    _alignment->seqsInfo  = NULL;
    
    for (int x = 0; getline( (*origin), line ), x < _alignment->sequenNumber; x++)
    {
        _alignment->seqsName[x] = std::string().append(&line[0u], 10);
        
        char * substring = strtok(&line[10], " ");
        
        do {
            _alignment->sequences[x].append(substring);
            substring = strtok(nullptr, " ");
            
        } while (substring != nullptr);
    }
    
    while (!origin->eof())
    {
        for (int x = 0; getline( (*origin), line ), x < _alignment->sequenNumber;)
        {
            if (line.empty()) 
            {
                if (origin->eof()) break;
                else continue;
            }

            char * substring = strtok(&line[0u], " ");
            
            do {
                _alignment->sequences[x].append(substring);
                substring = strtok(nullptr, " ");
                
            } while (substring != nullptr);
            
            x++;
        }
    }
    
    /* Check the matrix's content */
    _alignment->fillMatrices(false);
    
    return _alignment;
}

void PhylipInterleavedState::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
    setfill(' ');
    
    (*output) << " " << alignment->sequenNumber << " " << alignment->residNumber << endl << " I";
    
    int y = 0;
    for (int x = 0; x < alignment->sequenNumber; x++)
    {
        (*output) 
                << setw(13) << left
                << alignment->seqsName[x].substr(0, std::min(10, (int)alignment->seqsName[x].length()));
        
        for (y = 0; y < std::min((int)alignment->sequences[x].length(), 60); y+=10)
        {
            (*output) << alignment->sequences[x].substr(y, 10) << " ";
        }
        
        (*output) << endl;
        
    }
    
    int residCount = 60;
    bool continueTag = true;
    
    while (continueTag)
    {
        continueTag = false;
        
        (*output) << endl;
        
        for (int x = 0; x < alignment->sequenNumber; x++)
        {
            for (y = residCount; y < std::min((int)alignment->sequences[x].length(), residCount + 60); y+=10)
            {
                (*output) << alignment->sequences[x].substr(y, 10) << " ";
            }
            
            (*output) << endl;
            
            if (residCount + 60 < alignment->sequences[x].length())
                continueTag = true;
        }
        residCount += 60;
    }
}

bool PhylipInterleavedState::RecognizeOutputFormat(std::string FormatName)
{
    if (FormatName == "phylipinterleaved" || FormatName == "phylip_interleaved" ||
        FormatName == "PhylipInterleaved" || FormatName == "Phylip_Interleaved")
            return true;
    return false;
}

