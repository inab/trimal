#include "../../include/ReadWriteMS/fasta_state.h"
#include <iostream>
#include <../../home/vfernandez/git/trimal/include/defines.h>
#include <stdio.h>
#include <string>
#include "../../include/ReadWriteMS/ReadWriteMachineState.h"

using namespace std;

int FastaState::CheckAlignment(istream* origin)
{
    char c;
    origin->seekg(0);
    origin->get(c);
    if (!strcmp(&c, ">"))
        return 1;
    return 0;
}

newAlignment* FastaState::LoadAlignment(istream* origin)
{
    origin->clear();
    origin->seekg(0);

    newAlignment* _alignment = new newAlignment();

    /* Compute how many sequences are in the input alignment */
    _alignment->sequenNumber = 0;
    
    for( std::string line; getline( (*origin), line ); )
    {
        if (line[0] == '>')
            _alignment -> sequenNumber++;
    }
    
    /* Finish to preprocess the input file. */
    origin->clear();
    origin->seekg(0);

    /* Allocate memory for the input alignmet */
    _alignment->seqsName  = new string[_alignment->sequenNumber];
    _alignment->sequences = new string[_alignment->sequenNumber];
    _alignment->seqsInfo  = NULL;
    
    int i = -1;
    for( std::string line; getline( (*origin), line ); )
    {
        if (line[0] == '>')
        {
            int last_not = line.find_last_not_of("\t\n ");
            
            _alignment->seqsName[++i] = std::string().append(&line[1], last_not);
            _alignment->sequences[i] = "";
//             _alignment->seqsInfo[i] = std::string().append(&line[1], last_not);
        }
        else 
        {
            _alignment->sequences[i].append(line);
        }
    }

    /* Check the matrix's content */
    _alignment->fillMatrices(false);
    return _alignment;
}

void FastaState::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
    for (int i = 0; i < alignment->sequenNumber; i++)
    {
        if (Machine -> shortNames)
            (*output) << ">" << alignment->seqsName[i].substr(0, 10) << endl;
        else
            (*output) << ">" << alignment->seqsName[i] << endl;
        for (int x = 0; x < alignment->sequences[i].length(); x+= 60)
        {
            (*output) << alignment->sequences[i].substr(x, 60) << endl;
        }
    }
}

bool FastaState::RecognizeOutputFormat(std::string FormatName)
{
    if (FormatName == "fasta") return true;
    return false;
}
