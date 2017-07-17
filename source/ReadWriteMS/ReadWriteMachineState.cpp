#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/ReadWriteMS/ReadWriteBaseState.h"
#include "../../include/newAlignment.h"

#include "../../include/ReadWriteMS/clustal_state.h"
#include "../../include/ReadWriteMS/fasta_state.h"
#include "../../include/ReadWriteMS/pir_state.h"
#include "../../include/ReadWriteMS/phylip32_state.h"
#include "../../include/ReadWriteMS/phylip40_state.h"
#include "../../include/ReadWriteMS/phylip_paml_state.h"
#include "../../include/ReadWriteMS/nexus_state.h"
#include "../../include/ReadWriteMS/mega_interleaved_state.h"
#include "../../include/ReadWriteMS/mega_sequential_state.h"
#include "../../include/ReadWriteMS/htmlreport_state.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <array>

ReadWriteMS::ReadWriteMS()
{
    addState(new ClustalState(this));
    addState(new FastaState(this));
    addState(new PirState(this));
    addState(new Phylip32State(this));
    addState(new Phylip40State(this));
    addState(new PhylipPamlState(this));
    addState(new NexusState(this));
    addState(new MegaInterleavedState(this));
    addState(new MegaSequentialState(this));
    addState(new HTMLState(this));
};

ReadWriteMS::~ReadWriteMS()
{
    for(ReadWriteBaseState* child : available_states) {
        delete child;
    }
}

void ReadWriteMS::addState(ReadWriteBaseState* newState)
{
    this -> available_states.push_back(newState);
}

newAlignment* ReadWriteMS::loadAlignment(std::string inFile)
{
    ReadWriteBaseState* inState = NULL;
    
    ifstream inFileHandler;
    inFileHandler.open(inFile);
    if (!inFileHandler.is_open())
    {
        cerr << "Couldn't open input file" << endl;
        return nullptr;
    }

    int format_value = 0;
    int temp_value = 0;
    
    for(ReadWriteBaseState* state : available_states)
    {
        temp_value = state -> CheckAlignment(&inFileHandler);
        
        if (temp_value > format_value)
        {
            temp_value = format_value;
            inState = state;
        }
    }
    
    if (inState == NULL)
    {
        cerr << "Cannot recognize input format" << endl;
        inFileHandler.close();
        return nullptr;
    }
    
    newAlignment* alignment = inState->LoadAlignment(inFile);
    inFileHandler.close();
    return alignment;
}

void ReadWriteMS::saveAlignment(std::string outFile, std::string outFormat, newAlignment* alignment)
{
    ofstream output;
    output.open(outFile);
    for(ReadWriteBaseState* state : available_states)
    {
        if (state->RecognizeOutputFormat(outFormat))
        {
            state->SaveAlignment(alignment, &output, &alignment->filename);
            return;
        }
    }
    
}

void ReadWriteMS::processFile(std::string inFile, std::string outFile, std::string outFormat)
{
    ReadWriteBaseState* inState = NULL;
    ReadWriteBaseState* outState = NULL;
    
    ifstream inFileHandler;
    inFileHandler.open(inFile);
    if (!inFileHandler.is_open())
    {
        cerr << "Couldnt open input file" << endl;
        return;
    }
    
    {
        int format_value = 0;
        int temp_value = 0;
        
        for(ReadWriteBaseState* state : available_states)
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
    newAlignment* alignment = inState->LoadAlignment(inFile);
    inFileHandler.close();
    
    cout << "Saving Alignment in " << outState->name << " format" << endl;
    ofstream outFileHandler;
    outFileHandler.open(outFile);
    outState -> SaveAlignment(alignment, &outFileHandler, &inFile);
    outFileHandler.close();
    
    cout << "Alignment Saved" << endl;
    
    delete alignment;
}

void ReadWriteMS::processFile(std::string inFile, std::string outPattern, std::vector< std::string > outFormats[])
{
    ReadWriteBaseState* inState = NULL;
    std::vector<ReadWriteBaseState*> outStates = std::vector<ReadWriteBaseState*>();
    
    ifstream inFileHandler;
    inFileHandler.open(inFile);
    if (!inFileHandler.is_open())
    {
        cerr << "Couldn't open input file" << endl;
        return;
    }
    
    
    int format_value = 0;
    int temp_value = 0;
    
    for(ReadWriteBaseState* state : available_states)
    {
        temp_value = state -> CheckAlignment(&inFileHandler);
        
        if (temp_value > format_value)
        {
            temp_value = format_value;
            inState = state;
        }
    }
    
    if (inState == NULL)
    {
        cerr << "Cannot recognize input format" << endl;
        inFileHandler.close();
        return;
    }
    
    bool recognized;
    for (int formatID = 0; formatID < outFormats->size(); formatID++)
    {
        recognized = false;
        for(ReadWriteBaseState* state : available_states)
        {
            if (state->RecognizeOutputFormat(outFormats->at(formatID)))
            {
                outStates.push_back(state);
                recognized = true;
                break;
            }
        }
        
        if (!recognized)
        {
            cerr << "Cannot recognize output format" << endl;
            inFileHandler.close();
            return;
        }
    }
    
    cout << "Loading Alignment" << endl;
    newAlignment* alignment = inState->LoadAlignment(inFile);
    inFileHandler.close();
    {
        string filename;
        ofstream outFileHandler;
        for (ReadWriteBaseState * state : outStates)
        {
            cout << "Saving Alignment in " << state->name << " format" << endl;
            filename = outPattern;
            outFileHandler.open(filename.replace(filename.find('$'), 1, state->name));
            state -> SaveAlignment(alignment, &outFileHandler, &inFile);
            outFileHandler.close();
            cout << "Alignment Saved" << endl;
            
        }
    }
    
    delete alignment;
    
}
