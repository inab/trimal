#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/ReadWriteMS/pir_state.h"
#include "../../include/defines.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

int PirState::CheckAlignment(istream* origin)
{
    string str;
     origin->seekg(0);
     getline( (*origin), str );
     if (str[0] == '>' && str[3] == ';')
         return 2;
    return 0;
}

newAlignment* PirState::LoadAlignment(istream* origin)
{
    origin->clear();
    origin->seekg(0);

    newAlignment* alignment = new newAlignment();

    /* Compute how many sequences are in the input alignment */
    alignment->sequenNumber = 0;
    
    for( std::string line; getline( (*origin), line ); )
    {
        if (line[0] == '>')
            alignment -> sequenNumber++;
    }
    
    /* Finish to preprocess the input file. */
    origin->clear();
    origin->seekg(0);

    /* Allocate memory for the input alignmet */
    alignment->seqsName  = new string[alignment->sequenNumber];
    alignment->sequences = new string[alignment->sequenNumber];
    alignment->seqsInfo  = new string[alignment->sequenNumber];
    
    int i = -1;
    char * str_c = NULL;
    for( std::string line; getline( (*origin), line ); )
    {
        if (line[0] == '>')
        {
            strtok(&line[0u], ";");
            
            alignment->seqsName[++i] = std::string().append(strtok(NULL, ";"));
            getline( (*origin), line );
            alignment->seqsInfo[i] = std::string().append(&line[0], line.find_first_not_of(" "), line.length());
            alignment->sequences[i] = "";
            
            
        }
        else 
        {
            for (int x = 0; x < line.length(); x++)
            {
                if (line[x] == ' ' || line[x] == '*') continue;
                alignment->sequences[i].append(&line[x], 1);
            }
            
        }
    }
    
    delete str_c;
    
//     for (int y = 0; y < alignment->sequenNumber; y++)
//     {
//         cout << alignment->seqsName[y] << "_" << alignment->sequences[y] << endl;
//     }

    /* Check the matrix's content */
    alignment->fillMatrices(false);
    return alignment;
}

void PirState::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
   int i, j, k;
    string alg_datatype, *tmpMatrix;

    /* Allocate local memory for generating output alignment */
    tmpMatrix = new string[alignment->sequenNumber];

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    for(i = 0; i < alignment->sequenNumber; i++)
        tmpMatrix[i] = (!alignment->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    /* Compute output file datatype */
    alignment->getAlignmentType();
    if (alignment->dataType == DNAType)
        alg_datatype = "DL";
    else if (alignment->dataType == RNAType)
        alg_datatype = "RL";
    else if (alignment->dataType == AAType)
        alg_datatype = "P1";

    /* Print alignment */
    for(i = 0; i < alignment->sequenNumber; i++) {

        /* Print sequence datatype and its name */
        
        (*output) << ">" << alg_datatype << ";" ;
        
        if (Machine -> shortNames)
            (*output) << alignment->seqsName[i].substr(0, 10); 
        else
            (*output) << alignment->seqsName[i]; 
        
        if((alignment->seqsInfo != NULL))
            (*output) << endl << alignment->seqsInfo[i] << endl;
        else
            (*output) << endl << " " << alignment->residuesNumber[i] << " bases" << endl;

        /* Write the sequence */
        for(j = 0; j < alignment->residNumber; j += 50) {
            (*output) << " ";
            for(k = j; (k < alignment->residNumber) && (k < (j + 50)); k += 10)
                (*output) << " " << tmpMatrix[i].substr(k, 10);
            if((j + 50) >= alignment->residNumber)
                (*output) << "*";
            (*output) << endl;
        }
        (*output) << endl;
    }

    /* Deallocate local memory */
    delete [] tmpMatrix;
}

bool PirState::RecognizeOutputFormat(std::string FormatName)
{
    if (FormatName == "pir" || FormatName == "PIR" ||
        FormatName == "nbrf" || FormatName == "NBRF")
        return true;
    return false;
}

