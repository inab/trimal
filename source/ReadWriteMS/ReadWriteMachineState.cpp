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
#include <sstream>

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

bool ReadWriteMS::saveAlignment(std::string outFile, std::string outFormat, newAlignment* alignment)
{
    ofstream output;
    output.open(outFile);
    for(ReadWriteBaseState* state : available_states)
    {
        if (state->RecognizeOutputFormat(outFormat))
        {
            return state->SaveAlignment(alignment, &output, &alignment->filename);
        }
    }

}

void ReadWriteMS::processFile(
    std::vector< std::string >* inFiles,
    std::string* outPattern,
    std::vector< std::string >* outFormats)
{

    // Get all the output format states needed.
    std::vector<ReadWriteBaseState*> outStates = std::vector<ReadWriteBaseState*>();
    {
        bool recognized;
        for (std::string outFormat : *outFormats)
        {
            recognized = false;
            for(ReadWriteBaseState* state : available_states)
            {
                if (state->RecognizeOutputFormat(outFormat))
                {
                    outStates.push_back(state);
                    recognized = true;
                    break;
                }
            }

            if (!recognized)
            {
                cerr << "Cannot recognize output format " << outFormat << endl;
                return;
            }
        }
    }


    // Process input files one by one.
    ReadWriteBaseState* inState = NULL;
    ifstream inFileHandler;
    int format_value = 0;
    int temp_value = 0;

    for (std::string inFile : *inFiles)
    {
        // Open file.
        inFileHandler.open(inFile);
        if (!inFileHandler.is_open())
        {
            cerr << "Couldn't open input file " << inFile << endl;
            return;
        }

        inState = NULL;
        format_value = 0;
        temp_value = 0;

        // Get format State that can handle this alignment.
        for(ReadWriteBaseState* state : available_states)
        {
            temp_value = state -> CheckAlignment(&inFileHandler);

            if (temp_value > format_value)
            {
                temp_value = format_value;
                inState = state;
            }
        }

        // Check if there is a format State to handle the alignment.
        if (inState == NULL)
        {
            cerr << "Cannot recognize input format" << endl;
            inFileHandler.close();
            return;
        }

        // Load alignment one by one and store it on each of the formats specified.
        newAlignment* alignment = inState->LoadAlignment(inFile);
        if (reverse) alignment->setReverseFlag(true);
        int start, end;
        inFileHandler.close();
        {
            string filename;
            ofstream outFileHandler;
            for (ReadWriteBaseState * state : outStates)
            {
                start = inFile.find_last_of("/");
                end = inFile.find_last_of(".");
                filename = utils::ReplaceString(*outPattern, "[in]", inFile.substr(start, end-start));
                utils::ReplaceStringInPlace(filename, "[extension]", state->extension);
                utils::ReplaceStringInPlace(filename, "[format]", state->name);

                outFileHandler.open(filename);
                if (this->hasOutputFile)
                    state -> SaveAlignment(alignment, &outFileHandler, &inFile);
                else
                    state -> SaveAlignment(alignment, &cout, &inFile);
                outFileHandler.close();
            }
        }
        if (alignment != nullptr)
            delete alignment;
    }
}

std::string ReadWriteMS::getInputStateName(std::string inFile)
{
    ReadWriteBaseState* inState = NULL;
    ifstream inFileHandler;
    int format_value = 0;
    int temp_value = 0;

    // Open file.
    inFileHandler.open(inFile);
    if (!inFileHandler.is_open())
    {
        cerr << "Couldn't open input file " << inFile << endl;
        return "Unknown";
    }

    inState = NULL;
    format_value = 0;
    temp_value = 0;

    // Get format State that can handle this alignment.
    for(ReadWriteBaseState* state : available_states)
    {
        temp_value = state -> CheckAlignment(&inFileHandler);

        if (temp_value > format_value)
        {
            temp_value = format_value;
            inState = state;
        }
    }

    // Check if there is a format State to handle the alignment.
    if (inState == NULL)
    {
        inFileHandler.close();
        return "Unknown";
    }

    return inState->name;
}

std::string ReadWriteMS::getFormatsAvailable()
{
    std::stringstream ss("");

    for (ReadWriteBaseState* state : available_states)
    {
        ss << state->name << ", " ;
    }
    ss.seekp(-2, std::ios_base::end);
    ss << ". ";

    return ss.str();

}
