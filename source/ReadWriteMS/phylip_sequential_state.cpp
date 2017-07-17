#include "../../include/ReadWriteMS/phylip_sequential_state.h"
#include <iostream>
#include "../../include/defines.h"
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

int PhylipSequentialState::CheckAlignment(istream* origin)
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
    if (size == str.length())
        return 0;
//          cout << "Input Format: Phylip sequential" << endl;
    return 1;
}

newAlignment* PhylipSequentialState::LoadAlignment(istream* origin)
{
    origin->clear();
    origin->seekg(0);

    newAlignment* _alignment = new newAlignment();
    
    std::string line;
    // Get rid of the first line.
    getline( (*origin), line );
    
    /* Compute how many sequences are in the input alignment */
    _alignment->sequenNumber = stoi(strtok(&line[0u], " "));
    _alignment->residNumber = stoi(strtok(nullptr, " "));
    /* Allocate memory for the input alignmet */
    _alignment->seqsName  = new string[_alignment->sequenNumber];
    _alignment->sequences = new string[_alignment->sequenNumber];
    _alignment->seqsInfo  = NULL;
    
    int sequenceNumber = 0, residNumber = 0;
    while (!origin -> eof())
    {
        getline( (*origin), line);
        if (line.empty()) continue;
        
        _alignment->seqsName[sequenceNumber] = line.substr(0, 10);
        _alignment->sequences[sequenceNumber] = std::string();
        residNumber = 0;
        
        for (int x = 10; x < line.length(); x++)
        {
            if (line[x] == ' ') continue;
            _alignment->sequences[sequenceNumber].append(1, line[x]);
            residNumber++;
        }
        
        while (residNumber < _alignment->residNumber)
        {
            getline( (*origin), line);
            for (char c : line)
            {
                if (c == ' ') continue;
                _alignment->sequences[sequenceNumber].append(1, c);
                residNumber++;
            }
        }
        
        sequenceNumber += 1;
        
    }
    
    /* Check the matrix's content */
    _alignment->fillMatrices(false);
    
    return _alignment;
}

void PhylipSequentialState::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
    setfill(' ');
    
    (*output) << " " << alignment->sequenNumber << " " << alignment->residuesNumber << endl;
    
    for (int i = 0; i < alignment->sequenNumber; i++)
    {
        (*output) 
                << setw(13) << left
                << alignment->seqsName[i].substr(0, std::min(10, (int)alignment->seqsName[i].length()));
                
        for (int x = 0; x < alignment->sequences[i].length(); x+= 10)
        {
            if (x % 60 == 0) 
                (*output) << endl << setw(13) << left << "";
            
            (*output) << alignment->sequences[i].substr(x, 10) << " ";
        }
        
        (*output) << endl ;
    }
}

bool PhylipSequentialState::RecognizeOutputFormat(std::string FormatName)
{
    if (FormatName == "phylipsequential" || FormatName == "phylip_sequential" ||
        FormatName == "PhylipSequential" || FormatName == "Phylip_Sequential")
            return true;
    return false;
}

