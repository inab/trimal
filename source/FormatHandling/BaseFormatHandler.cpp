#include "../../include/FormatHandling/FormatManager.h"
#include "../../include/FormatHandling/BaseFormatHandler.h"
#include "../../include/FormatHandling/formats_header.h"
#include "../../include/TimerFactory.h"
#include "../../include/newAlignment.h"
#include "../../include/utils.h"

FormatManager::~FormatManager()
{
    for(BaseFormatHandler* child : available_states) {
        delete child;
    }
}

void FormatManager::addState(BaseFormatHandler* newState)
{
    this -> available_states.push_back(newState);
}

newAlignment* FormatManager::loadAlignment(std::string inFile)
{
    // Check input file.
    std::ifstream inFileHandler;
    inFileHandler.open(inFile);
    if (!inFileHandler.is_open())
    {
        debug.report(ErrorCode::CantOpenFile, &inFile[0]);
        return nullptr;
    }
    else if (inFileHandler.peek() == std::ifstream::traits_type::eof())
    {
        debug.report(ErrorCode::FileIsEmpty, &inFile[0]);
        return nullptr;
    }

    BaseFormatHandler* inState = nullptr;
    int format_value = 0;
    int temp_value = 0;

    for(BaseFormatHandler* state : available_states)
    {
        temp_value = state -> CheckAlignment(&inFileHandler);

        if (temp_value > format_value)
        {
            format_value = temp_value;
            inState = state;
        }
    }

    if (inState == nullptr)
    {
        debug.report(ErrorCode::AlignmentFormatNotRecognized, &inFile[0]);
        inFileHandler.close();
        return nullptr;
    }

    newAlignment* alignment = inState->LoadAlignment(inFile);
    inFileHandler.close();
    return alignment;
}

bool FormatManager::saveAlignment(std::string outPattern, std::vector< std::string >* outFormats, newAlignment* alignment)
{
    StartTiming("bool FormatManager::saveAlignment()");
    if (alignment->residNumber == 0 || alignment->sequenNumber == 0)
    {
        debug.report(ErrorCode::AlignmentIsEmpty);
        return false;
    }
    std::string filename;
    unsigned long start;
    unsigned long end;
    if (alignment->filename.empty())
    {
        filename = utils::ReplaceString(outPattern, "[in]", "NoInputFileName");
    }
    else
    {
        start = std::max((int)alignment->filename.find_last_of('/'), 0);
        end = alignment->filename.find_last_of('.');
        filename = utils::ReplaceString(outPattern, "[in]", alignment->filename.substr(start, end-start));
    }

    
    if (outPattern.empty())
    {
        if (outFormats->empty())
        {
            outFormats->push_back("fasta");
        }
        if (outFormats->size() == 1)
        {
            for(BaseFormatHandler* state : available_states)
            {
                if (state->RecognizeOutputFormat( outFormats->at(0) ))
                {
                    return state->SaveAlignment(alignment, &std::cout, &filename);
                }
            }
            debug.report(ErrorCode::OutputFormatNotRecognized, &outFormats->at(0)[0]);
            return false;
        }
        else
        {
            debug.report(ErrorCode::OnlyOneFormatOnConsoleOutput);
            return false;
        }
    }
    else
    {
        if (outFormats->empty())
        {
            outFormats->push_back("fasta");
        }
        if (outFormats->size() == 1)
        {
            for(BaseFormatHandler* state : available_states)
            {
                if (state->RecognizeOutputFormat( outFormats->at(0) ))
                {
                    utils::ReplaceStringInPlace(filename, "[extension]", state->extension);
                    utils::ReplaceStringInPlace(filename, "[format]", state->name);
                    
                    std::ofstream outFileHandler;
                    outFileHandler.open(filename);
                    
                    return state->SaveAlignment(alignment, &outFileHandler, &alignment->filename);
                }
            }
            debug.report(ErrorCode::OutputFormatNotRecognized, &outFormats->at(0)[0]);
            return false;
        }
        else 
        {
            // Tranform the list of std::string states to a BaseFormatHandler list.
            bool recognized;
            std::vector<BaseFormatHandler*> outStates = std::vector<BaseFormatHandler*>();
            {
                for (std::string outFormat : *outFormats)
                {
                    recognized = false;
                    for(BaseFormatHandler* state : available_states)
                    {
                        if (state->RecognizeOutputFormat(outFormat))
                        {
                            outStates.push_back(state);
                            recognized = true;
                        }
                    }
                    if (!recognized)
                    {
                        debug.report(ErrorCode::OutputFormatNotRecognized, &outFormat[0]);
                        return false; 
                    }
                }
            }
            std::ofstream outFileHandler;
            bool isCorrect = true;
            std::string filename_2;
            for (BaseFormatHandler * state : outStates)
            {
                filename_2 = utils::ReplaceString(filename, "[extension]", state->extension);
                utils::ReplaceStringInPlace(filename_2, "[format]", state->name);
                outFileHandler.open(filename_2);
                if(!state -> SaveAlignment(alignment, &outFileHandler, &alignment->filename))
                {
                    debug.report(ErrorCode::AlignmentNotSaved, &state->name[0]);
                    isCorrect = false;
                }
                
                outFileHandler.close();
            }
            if (!isCorrect)
                debug.report(ErrorCode::ImpossibleToGenerate, new std::string[1]{"the output file"});
            return isCorrect;
        }
    }
}

void FormatManager::loadAndSaveMultipleAlignments(
    std::vector< std::string >* inFiles,
    std::string* outPattern,
    std::vector< std::string >* outFormats)
{

    // Get all the output format states needed.
    std::vector<BaseFormatHandler*> outStates = std::vector<BaseFormatHandler*>();
    {
        bool recognized;
        for (std::string outFormat : *outFormats)
        {
            recognized = false;
            for(BaseFormatHandler* state : available_states)
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
                debug.report(ErrorCode::OutputFormatNotRecognized, &outFormat[0]);
                return;
            }
        }
    }

    // Process input files one by one.
    BaseFormatHandler* inState = nullptr;
    std::ifstream inFileHandler;
    int format_value = 0;
    int temp_value = 0;

    for (std::string inFile : *inFiles)
    {
        // Open file.
        inFileHandler.open(inFile);
        if (!inFileHandler.is_open())
        {
            debug.report(ErrorCode::CantOpenFile, &inFile[0]);
            return;
        }
        else if (inFileHandler.peek() == std::ifstream::traits_type::eof())
        {
            debug.report(ErrorCode::FileIsEmpty, &inFile[0]);
            return;
        }

        inState = nullptr;
        format_value = 0;

        // Get format State that can handle this alignment.
        for(BaseFormatHandler* state : available_states)
        {
            temp_value = state -> CheckAlignment(&inFileHandler);

            if (temp_value > format_value)
            {
                format_value = temp_value;
                inState = state;
            }
        }

        // Check if there is a format State to handle the alignment.
        if (inState == nullptr)
        {
            debug.report(ErrorCode::AlignmentFormatNotRecognized, &inFile[0]);
            inFileHandler.close();
            return;
        }

        // Load alignment one by one and store it on each of the formats specified.
        newAlignment* alignment = inState->LoadAlignment(inFile);
        unsigned long start;
        unsigned long end;
        inFileHandler.close();
        {
            std::string filename;
            std::ofstream outFileHandler;
            for (BaseFormatHandler * state : outStates)
            {
                start = std::max((int)inFile.find_last_of('/'), 0);
                end = inFile.find_last_of('.');
                filename = utils::ReplaceString(*outPattern, "[in]", inFile.substr(start, end-start));
                utils::ReplaceStringInPlace(filename, "[extension]", state->extension);
                utils::ReplaceStringInPlace(filename, "[format]", state->name);

                outFileHandler.open(filename);
                if (this->hasOutputFile)
                    state -> SaveAlignment(alignment, &outFileHandler, &inFile);
                else
                    state -> SaveAlignment(alignment, &std::cout, &inFile);
                outFileHandler.close();
            }
        }

        delete alignment;
    }
}

std::string FormatManager::getFileFormatName(std::string inFile)
{
    std::ifstream inFileHandler;

    // Open file.
    inFileHandler.open(inFile);
    if (!inFileHandler.is_open())
    {
        debug.report(ErrorCode::CantOpenFile, &inFile[0]);
        return "Unknown";
    }
    else if (inFileHandler.peek() == std::ifstream::traits_type::eof())
    {
        debug.report(ErrorCode::FileIsEmpty, &inFile[0]);
        return "None";
    }

    BaseFormatHandler* inState = nullptr;
    int format_value = 0;
    int temp_value = 0;

    // Get format State that can handle this alignment.
    for(BaseFormatHandler* state : available_states)
    {
        temp_value = state -> CheckAlignment(&inFileHandler);

        if (temp_value > format_value)
        {
            format_value = temp_value;
            inState = state;
        }
    }

    // Check if there is a format State to handle the alignment.
    if (inState == nullptr)
    {
        inFileHandler.close();
        return "Unknown";
    }

    return inState->name;
}

std::string FormatManager::getInputFormatsAvailable()
{
    std::stringstream ss("");

    for (BaseFormatHandler* state : available_states)
    {
        if (state->canLoad)
            ss << state->name << ", " ;
    }
    ss.seekp(-2, std::ios_base::end);
    ss << "  ";

    return ss.str();

}

std::string FormatManager::getOutputFormatsAvailable()
{
    std::stringstream ss("");

    for (BaseFormatHandler* state : available_states)
    {
        if (state->canSave)
            ss << state->name << ", " ;
    }
    ss.seekp(-2, std::ios_base::end);
    ss << "  ";

    return ss.str();

}

std::vector<newAlignment*> FormatManager::splitAlignmentKeeping(newAlignment& alignment)
{
    std::vector<newAlignment *> splitted = std::vector<newAlignment *>(alignment.originalSequenNumber);

    for (int i = 0; i < alignment.originalSequenNumber; i++)
    {
        newAlignment * tempAlignment = new newAlignment();
        tempAlignment->sequences = new std::string[1];
        tempAlignment->sequences[0] = std::string(alignment.sequences[i]);
        tempAlignment->seqsName = new std::string[1] { alignment.seqsName[i] };
        tempAlignment->sequenNumber = 1;
        tempAlignment->originalSequenNumber = 1;
        tempAlignment->residNumber = tempAlignment->sequences[0].size();
        tempAlignment->originalResidNumber = tempAlignment->residNumber;
        tempAlignment->filename = tempAlignment->seqsName[0];
        splitted[i] = tempAlignment;
    }

    return splitted;
}
