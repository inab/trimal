#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/ReadWriteMS/readwrites.h"
#include "../../include/newAlignment.h"

#include "../../include/ReadWriteMS/fasta_state.h"
#include "../../include/ReadWriteMS/clustal_state.h"
#include "../../include/ReadWriteMS/pir_state.h"
#include "../../include/ReadWriteMS/phylip_interleaved_state.h"
#include "../../include/ReadWriteMS/phylip_sequential_state.h"
#include "../../include/ReadWriteMS/mega_interleaved_state.h"
#include "../../include/ReadWriteMS/mega_sequential_state.h"


#include <iostream>
#include <fstream>

ReadWriteMS::ReadWriteMS()
{
    addState(new FastaState(this));
    addState(new ClustalState(this));
    addState(new PirState(this));
    addState(new PhylipInterleavedState(this));
    addState(new PhylipSequentialState(this));
    addState(new MegaInterleavedState(this));
    addState(new MegaSequentialState(this));
};

ReadWriteMS::~ReadWriteMS()
{
    for(readwrites* child : available_states) {
        delete child;
    }
}

void ReadWriteMS::addState(readwrites* newState)
{
    this -> available_states.push_back(newState);
}

void ReadWriteMS::processFile(std::string inFile, std::string outFile, std::string outFormat)
{
    inState = NULL;
    outState = NULL;
    
    ifstream inFileHandler;
    inFileHandler.open(inFile);
    if (!inFileHandler.is_open())
    {
        cerr << "Couldnt open input file" << endl;
        return;
    }
    
    // Recognize format. 
    //      As some formats have similar structure (Pir vs Fasta), 
    //      we use the one with a higher coincidence.
    {
        int format_value = 0;
        int temp_value = 0;
        
        for(readwrites* state : available_states)
        {
            temp_value = state -> CheckAlignment(&inFileHandler);
            
            if (temp_value > format_value)
            {
                temp_value = format_value;
                inState = state;
            }
            
            if (state->RecognizeOutputFormat(outFormat))
            {
                outState = state;
            }
        }
    }
    
    
    if (inState == NULL || outState == NULL)
    {
        if (inState == NULL)
        {
            cerr << "Cannot recognize input format" << endl;
        }
        if (outState == NULL)
        {
            cerr << "Cannot recognize output format" << endl;
        }
        
        inFileHandler.close();
        return;
    }
    
    cout << "Loading Alignment" << endl;
    newAlignment* alignment = inState->LoadAlignment(&inFileHandler);
    inFileHandler.close();
    
    cout << "Saving Alignment" << endl;
    ofstream outFileHandler;
    outFileHandler.open(outFile);
    outState -> SaveAlignment(alignment, &outFileHandler, &inFile);
    outFileHandler.close();
    
    cout << "Alignment Saved" << endl;
    
    delete alignment;
}

// void ReadWriteMS::processFile(std::string inFile, std::string outPattern, std::string outFormats[]){}
//     
// void ReadWriteMS::processFile(std::string inFiles[], std::string outPattern, std::string outFormat) {}
// 
// void ReadWriteMS::processFile(std::string inFiles[], std::string outPattern, std::string outFormats[]){}
