#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/ReadWriteMS/clustal_state.h"
#include "../../include/defines.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

int ClustalState::CheckAlignment(istream* origin)
{
    string str;
     origin->seekg(0);
     getline( (*origin), str );
     str = str.substr(0, 7);
     if (!strcmp(str.c_str(), "CLUSTAL"))
         return 1;
 
    return 0;
}

newAlignment* ClustalState::LoadAlignment(istream* origin)
{
    origin->clear();
    origin->seekg(0);

    newAlignment* _alignment = new newAlignment();

    /* Compute how many sequences are in the input alignment */
    _alignment->sequenNumber = 0;
    
    std::string line;
    // Get rid of the first line.
    getline( (*origin), line );
    // Get rid of blank lines.
    for( ; getline( (*origin), line ); )
    {
        if (!line.empty())
        {
            break;
        }
    }
    
    // Start a names vector.
    std::vector<std::string> names;
    do
    {
        // Read until a blank 
        if (line.empty() || line[0] == ' ')
            break;
        else
        {
            // We also get the names from first iteration.
            names.push_back(std::string().append(strtok(&line[0u], " ")));
            
        }
        
        getline( (*origin), line );
        
    } while (!origin->eof());
    
    _alignment -> sequenNumber = names.size();

    /* Allocate memory for the input alignmet */
    _alignment->seqsName  = new string[_alignment->sequenNumber];
    _alignment->sequences = new string[_alignment->sequenNumber];
    _alignment->seqsInfo  = NULL;
    
    // As we got the names previously, we can copy them from the vector
    std::copy(names.begin(), names.end(), _alignment->seqsName);
    names.clear();
    
    int i = 0;
    for ( ; i < _alignment->sequenNumber; i++)
    {
        _alignment->sequences[i] = "";
    }
    
    
    /* Finish to preprocess the input file. */
    origin->clear();
    origin->seekg(0);
    getline( (*origin), line );
    for( i = -1; getline( (*origin), line ); )
    {
        if (line.empty() || line[0] == ' ')
            i = -1;
        else 
        {
            strtok(&line[0u], " ");
            string strr = strtok(NULL, " ");
            _alignment->sequences[++i].append(strr);
        }
    }

    
    /* Check the matrix's content */
    _alignment->fillMatrices(false);
    
    return _alignment;
}

void ClustalState::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
    int max = 0;
    
    if (Machine -> shortNames)
    {
        max = 15;
    }
    else 
    {
        for (int x = 0; x < alignment->sequenNumber; x++)
        {
            max = std::max((int)alignment->seqsName[x].length(), max);
        }
        
        max += 5;
    }
    
    setfill(' ');
    
    (*output) << "CLUSTAL W" << endl << endl;
    
    for (int x = 0; x < alignment->residNumber; x+= 60)
    {
        for (int i = 0; i < alignment->sequenNumber; i++)
        {
            if (Machine -> shortNames)
                (*output) << setw(max) << left << alignment->seqsName[i].substr(0, 10);
            else
                (*output) << setw(max) << left << alignment->seqsName[i];
            (*output) << alignment->sequences[i].substr(x, 60) << endl;
        }
        
        if (x + 60 < alignment->residNumber)
            (*output) << endl << endl;
    }
}

bool ClustalState::RecognizeOutputFormat(std::string FormatName)
{
    if (FormatName == "clustal") return true;
    return false;
}

