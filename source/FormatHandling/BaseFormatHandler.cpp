/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2009-2019
        Fernandez-Rodriguez V.  (victor.fernandez@bsc.es)
        Capella-Gutierrez S.    (salvador.capella@bsc.es)
        Gabaldon, T.            (tgabaldon@crg.es)

    This file is part of trimAl/readAl.

    trimAl/readAl are free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl/readAl are distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl/readAl. If not, see <http://www.gnu.org/licenses/>.

***************************************************************************** */

#include <limits>

#include <FormatHandling/FormatManager.h>
#include "FormatHandling/FormatManager.h"
#include "FormatHandling/BaseFormatHandler.h"
#include "FormatHandling/formats_header.h"
#include "InternalBenchmarker.h"
#include "Alignment/Alignment.h"
#include "utils.h"


namespace FormatHandling {

    FormatManager::~FormatManager()
    {
        for(BaseFormatHandler* child : available_states) {
            delete child;
        }
    }

    void FormatManager::addState(FormatHandling::BaseFormatHandler* newState)
    {
        this -> available_states.emplace_back(newState);
    }

    Alignment* FormatManager::loadAlignment(const std::string &inFile)
    {
        auto inState = getFormatFromFile(inFile);

        if (inState == nullptr)
        {
            debug.report(ErrorCode::AlignmentFormatNotRecognized, &inFile[0]);
            return nullptr;
        }

        return inState->LoadAlignment(inFile);
    }

    std::string FormatManager::replaceINtag(
            const Alignment & alignment,
            const std::string & outPattern)
    {
        // Replace the [in] tag with correspondent information.
        unsigned long start, end;
        // Although an alignment file should have a filename
        //      we make a guard to prevent issues.
        if (alignment.filename.empty())
        {
            return utils::ReplaceString(outPattern, "[in]", "NoInputFileName");
        }
            // If the alignment contains a filename
            //      we replace the tag with it's name
        else
        {
            start = std::max((int)alignment.filename.find_last_of('/') + 1, 0);
            end = alignment.filename.find_last_of('.');
            return utils::ReplaceString(outPattern, "[in]", alignment.filename.substr(start, end-start));
        }
    }

    bool FormatManager::saveAlignment(
            const std::string &outPattern,
            const std::vector<std::string> &outFormats,
            Alignment &alignment)
    {
        const std::vector<Alignment*>tmpVector{ &alignment };
        return saveAlignments(outPattern, outFormats, tmpVector);
    }

    bool FormatManager::saveAlignments(
            const std::string &outPattern,
            const std::vector<std::string> &outFormats,
            const std::vector<Alignment *> &alignmentVector)
    {
        StartTiming("bool FormatManager::saveAlignment()");

        bool returnValue = true;

        std::vector<FormatHandling::BaseFormatHandler*> formats{};

        for(const std::string & outFormat : outFormats)
        {
            // Get the format desired
            BaseFormatHandler* state = getFormatFromToken(outFormat);
            if (state == nullptr)
            {
                debug.report(ErrorCode::OutputFormatNotRecognized, &outFormats.at(0)[0]);
                returnValue = false;
                continue;
            }
            formats.emplace_back(state);
        }


        // Check if any format has been requested.
        //      If not, the default is fasta.
        if (formats.empty())
            formats.emplace_back(getFormatFromToken("fasta"));


        // If we are outputting to terminal
        if (outPattern.empty())
        {
            if (formats.size() != 1)
            {
                // Report an error if more than one format is requested on terminal
                debug.report(ErrorCode::OnlyOneFormatOnConsoleOutput);
                return false;
            }

            return formats[0]->SaveAlignment(*alignmentVector[0], &std::cout);
        }

        // Iterate over all alignments
        for (const Alignment * alignment : alignmentVector)
        {
            // Check if is possible to save the alignment
            if (alignment->numberOfResidues  == 0 ||
                alignment->numberOfSequences == 0)
            {
                debug.report(ErrorCode::AlignmentIsEmpty);
                returnValue = false;
                continue;
            }

            // Retrieve the IN parameter
            std::string filename = replaceINtag(*alignment, outPattern);

            for(FormatHandling::BaseFormatHandler * state : formats)
            {
                // Replace the [extension] and [format] tokens
                std::string finalFilename =
                        utils::ReplaceString(filename, "[extension]", state->extension);
                utils::ReplaceStringInPlace(finalFilename, "[format]", state->name);
// Macro passed as compile time. See FormatHandlerOverwritePolicy.cmake
#if FormatHandlerOverwrites && false
                if (utils::fileExists(finalFilename) && (openmode & std::ofstream::app) == 0)
                {
                    debug.report(OverwrittingFile, new std::string[2]{state->name, finalFilename} );
                }
#else
                if((openmode & std::ofstream::app) == 0)
                {

                    uint i;
                    for (i = 0; i < std::numeric_limits<uint>::max(); i++)
                    {
                        if (!utils::fileExists(finalFilename + "." + std::to_string(i)))
                        {
                            debug.report(RenamingOutputPreventOverride,
                                         new std::string[3]{state->name, finalFilename, finalFilename + "." + std::to_string(i)} );
                            finalFilename += "." + std::to_string(i);
                            break;
                        }
                    }
                    if (i == std::numeric_limits<uint>::max())
                    {
// Macro passed as compile time. See FormatHandlerOverwritePolicy.cmake
#if FormatHandlerOverwritesOriginal
                        debug.report(TriedRenamingOutputPreventOverride,
                            new std::string[3]{state->name, finalFilename, finalFilename} );
#else
                        debug.report(TriedRenamingOutputPreventOverride,
                                new std::string[3]{state->name, finalFilename, finalFilename + "." + std::to_string(i)} );
                        finalFilename += "." + std::to_string(i);
#endif
                    }
#endif
                }

                // Open the file handler
                std::ofstream outFileHandler(finalFilename, openmode);

                // Save the alignment
                if (!state->SaveAlignment(*alignment, &outFileHandler))
                {
                    debug.report(ErrorCode::AlignmentNotSaved, &state->name[0]);
                    returnValue = false;
                }
            }
        }

        return returnValue;
    }

    void FormatManager::loadAndSaveMultipleAlignments(
            const std::vector<std::string> &inFiles,
            const std::string &outPattern,
            const std::vector<std::string> &outFormats)
    {
        // Store all the alignments
        std::vector<Alignment *> alignments {};
        // Iterate over all of them
        for(const std::string & inFile : inFiles)
        {
            Alignment * alignment = loadAlignment(inFile);
            if (alignment == nullptr) continue;
            alignments.emplace_back(alignment);
        }
        // Save them
        saveAlignments(outPattern, outFormats, alignments);

        for(Alignment * alignment : alignments)
        {
            delete alignment;
            alignment = nullptr;
        }
    }

    std::string FormatManager::getFileFormatName(const std::string &inFile)
    {
        FormatHandling::BaseFormatHandler * format = getFormatFromFile(inFile);
        return format == nullptr ? "Unknown" : format->name;
    }

    std::string FormatManager::getInputFormatsAvailable()
    {
        std::stringstream ss("");
        for (BaseFormatHandler* state : available_states)
            if (state->canLoad)
                ss << state->name << ", " ;
        ss.seekp(-2, std::ios_base::end);
        ss << "  ";
        return ss.str();
    }

    std::string FormatManager::getOutputFormatsAvailable()
    {
        std::stringstream ss("");
        for (BaseFormatHandler* state : available_states)
            if (state->canSave)
                ss << state->name << ", " ;
        ss.seekp(-2, std::ios_base::end);
        ss << "  ";
        return ss.str();
    }

    std::vector<Alignment*> FormatManager::splitAlignmentKeeping(
            const Alignment &alignment)
    {
        std::vector<Alignment *> splitted =
                std::vector<Alignment *>(alignment.originalNumberOfSequences);

        for (int i = 0; i < alignment.originalNumberOfSequences; i++)
        {
            Alignment * tempAlignment = new Alignment();
            tempAlignment->sequences = new std::string[1]{alignment.sequences[i]};
            tempAlignment->seqsName = new std::string[1]{ alignment.seqsName[i] };
            tempAlignment->numberOfSequences
                    = tempAlignment->originalNumberOfSequences
                    = 1;
            tempAlignment->numberOfResidues
                    = tempAlignment->originalNumberOfResidues
                    = (int) tempAlignment->sequences[0].size();
            tempAlignment->filename = tempAlignment->seqsName[0];
            splitted[i] = tempAlignment;
        }

        return splitted;
    }

    FormatHandling::BaseFormatHandler * FormatManager::getFormatFromToken(
            const std::string &token)
    {
        // Search for a format able to recognize the token
        for(BaseFormatHandler* state : available_states)
            if (state->RecognizeOutputFormat( token ))
                return state;

        return nullptr;
    }

    FormatHandling::BaseFormatHandler * FormatManager::getFormatFromFile(
            const std::string &filename)
    {
        // Open the file
        std::ifstream * input = getNonEmptyFile(filename);

        if (input == nullptr)
        {
            delete input;
            return nullptr;
        }

        // Temporal value
        BaseFormatHandler* inState = nullptr;

        // Compare variables
        int format_value = 0, temp_value = 0;

        // Iterate over all available states
        for(BaseFormatHandler* state : available_states)
        {
            // Store the score that state produces
            temp_value = state -> CheckAlignment(input);

            // If the score is better than last maximum, this format
            //      is better at recognizing the input file
            if (temp_value > format_value)
            {
                format_value = temp_value;
                inState = state;
            }
        }

        // We don't need the input open anymore. Destructor deletes it
        delete input;

        // Check if any input recognized the format
        if (inState == nullptr)
        {
            debug.report(ErrorCode::AlignmentFormatNotRecognized, &filename[0]);
            return nullptr;
        }

        // Return the recognized format
        return inState;
    }

    std::ifstream * FormatManager::getNonEmptyFile(
            const std::string &filename)
    {
        auto inFileHandler = new std::ifstream(filename);
        // Input file must exist
        if (!inFileHandler->is_open())
        {
            debug.report(ErrorCode::CantOpenFile, &filename[0]);
            delete inFileHandler;
            return nullptr;
        }
        // Input file cannot be empty
        if (inFileHandler->peek() == std::ifstream::traits_type::eof())
        {
            debug.report(ErrorCode::FileIsEmpty, &filename[0]);
            delete inFileHandler;
            return nullptr;
        }


        return inFileHandler;
    }
}