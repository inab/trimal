#include <iomanip>
#include "../include/reportsystem.h"

ReportSystem Debug = ReportSystem();

// Out of line initializations.

ReportSystem::VerboseLevel ReportSystem::Level = VerboseLevel::ERROR;

const std::map<ReportSystem::InfoCode, const char *> ReportSystem::InfoMessages = 
{
    { ReportSystem::InfoCode::CuttingSequence, 
        "Cutting sequence \"[ŧag]\" at first appearance of stop codon \"[tag]\" (residue \"[tag]\") at position [tag] (length: [tag] \")" }, 
        
    { ReportSystem::InfoCode::WindowSizeCompareset, 
        "Try with specific comparison file window value. Parameter -cw" }, 
};
    
const std::map<ReportSystem::WarningCode, const char *> ReportSystem::WarningMessages = 
{
    { ReportSystem::WarningCode::RemovingOnlyGapsColumn, 
        "Removing column '[tag]' composed only by gaps" }, 
        
    { ReportSystem::WarningCode::KeepingOnlyGapsColumn, 
        "Keeping column '[tag]' composed only by gaps" }, 
        
    { ReportSystem::WarningCode::SequenceWillBeCutted, 
        "Sequence \"[tag]\" will be cutted at position [tag] (length:[tag])" }, 
        
    { ReportSystem::WarningCode::IncludingIndeterminationSymbols,
        "Sequence \"[tag]\" has some indetermination symbols 'X' at the end of sequence. They will be included in the final newAlignment." },
        
    { ReportSystem::WarningCode::LessNucleotidesThanExpected,
        "Sequence \"[tag]\" has less nucleotides ([tag]) than expected ([tag]). It will be added N's to complete the sequence" },
        
    { ReportSystem::WarningCode::HeaderWillBeCut,
        "Original sequence header will be cutted by 10 characters" }
        
};
    
const std::map<ReportSystem::ErrorCode, const char *> ReportSystem::ErrorMessages = 
{
    { ReportSystem::ErrorCode::AlignmentNotLoaded, 
        "Alignment not loaded: \" [tag] \" Check the file's content" }, 
        
    { ReportSystem::ErrorCode::NoFormatsSpecified, 
        "You must specify at least one format after the '-formats' argument" }, 
        
    { ReportSystem::ErrorCode::AlternativeMatrixNotRecognized, 
        "Alternative matrix \" [tag] \" not recognized" }, 
        
    { ReportSystem::ErrorCode::ReferenceFileNotLoaded, 
        "Reference file \" [tag] \" not loaded" }, 
        
    { ReportSystem::ErrorCode::GapThresholdOutOfRange, 
        "The gap threshold value should be between 0 and 1" }, 
        
    { ReportSystem::ErrorCode::GapThresholdNotRecognized, 
        "The gap threshold value should be a positive real number between 0 and 1" }, 
        
    { ReportSystem::ErrorCode::SimilarityThresholdOutOfRange, 
        "The similarity threshold value should be between 0 and 1" }, 
        
    { ReportSystem::ErrorCode::SimilarityThresholdNotRecognized, 
        "The similarity threshold value should be a positive real number between 0 and 1" }, 
        
    { ReportSystem::ErrorCode::ConsistencyThresholdOutOfRange, 
        "The consistency threshold value should be between 0 and 1" }, 
        
    { ReportSystem::ErrorCode::ConsistencyThresholdNotRecognized, 
        "The consistency threshold value should be a positive real number between 0 and 1" }, 
        
    { ReportSystem::ErrorCode::ConservationThresholdOutOfRange, 
        "The conservation threshold value should be between 0 and 100" }, 
        
    { ReportSystem::ErrorCode::ConservationThresholdNotRecognized, 
        "The conservation threshold value should be a positive real number between 0 and 100" }, 
        
    { ReportSystem::ErrorCode::ResidueOverlapOutOfRange, 
        "The residue overlap value should be between 0 and 100" }, 
        
    { ReportSystem::ErrorCode::ResidueOverlapNotRecognized, 
        "The residue overlap value should be a positive real number between 0 and 100" }, 
        
    { ReportSystem::ErrorCode::SequencesOverlapOutOfRange, 
        "The sequences overlap value should be between 0 and 100" }, 
        
    { ReportSystem::ErrorCode::SequencesOverlapNotRecognized, 
        "The sequences overlap value should be a positive real number between 0 and 100" }, 
        
    { ReportSystem::ErrorCode::MaxIdentityOutOfRange, 
        "The max identity value should be between 0 and 1" }, 
        
    { ReportSystem::ErrorCode::MaxIdentityNotRecognized, 
        "The max identity value should be a positive real number between 0 and 1" }, 
        
    { ReportSystem::ErrorCode::ClustersValueOutOfRange, 
        "The clusters value should be greater than 0" }, 
        
    { ReportSystem::ErrorCode::ClustersValueNotRecognized, 
        "The clusters value should be a positive integer number greater than 0" }, 
        
    { ReportSystem::ErrorCode::WindowValueOutOfRange, 
        "The window value should be greater than 0" }, 
        
    { ReportSystem::ErrorCode::WindowValueNotRecognized, 
        "The window value should be a positive integer number greater than 0" }, 
        
    { ReportSystem::ErrorCode::SelectSeqsNotRecognized, 
        "Could not parse the -selectseqs ranges" }, 
        
    { ReportSystem::ErrorCode::SelectColsNotRecognized, 
        "Could not parse the -selectres ranges" }, 
        
    { ReportSystem::ErrorCode::GapWindowValueOutOfRange, 
        "The gap window value should be greater than 0" }, 
        
    { ReportSystem::ErrorCode::GapWindowValueNotRecognized, 
        "The gap window value should be a positive integer number grater than 0" }, 
        
    { ReportSystem::ErrorCode::SimilarityWindowValueOutOfRange, 
        "The similarity window value should be greater than 0" }, 
        
    { ReportSystem::ErrorCode::SimilarityWindowValueNotRecognized, 
        "The similarity window value should be a positive integer number grater than 0" }, 
        
    { ReportSystem::ErrorCode::ConsistencyThresholdOutOfRange, 
        "The consistency window value should be greater than 0" }, 
        
    { ReportSystem::ErrorCode::ConsistencyThresholdNotRecognized, 
        "The consistency window value should be a positive integer number grater than 0" }, 
        
    { ReportSystem::ErrorCode::BlockSizeOutOfRange, 
        "The consistency window value should be greater than 0" }, 
        
    { ReportSystem::ErrorCode::BlockSizeNotRecognized, 
        "The consistency window value should be a positive integer number grater than 0" }, 
        
    { ReportSystem::ErrorCode::InFileComparisonStatistics, 
        "Option -in not valid in combination with file comparision options: -sft || -sfc || -ct" }, 
        
    { ReportSystem::ErrorCode::IncompatibleArguments, 
        "Argument [tag] is not compatible with argument [tag]" }, 
        
    { ReportSystem::ErrorCode::SelectSeqsResAndThresholdIncompatibilities, 
        "Options -selectCols and -selectSeqs are not allowed in combination with threshold options: (-gt || -st || -ct || -cons)" }, 
    
    { ReportSystem::ErrorCode::SelectSeqsResAndAutomathedMethodsIncompatibilities, 
        "Options -selectCols and -selectSeqs are not allowed in combination with automated trimming options: (-nogaps || -noallgaps || -gappyout || -strict || -strictplus || -automated1)." }, 
        
    { ReportSystem::ErrorCode::SelectSeqsResAndWindowIncompatibilities, 
        "Options -selectCols and -selectSeqs are not allowed in combination with window options: (-w || -sw || -gw || -wc)" }, 
        
    { ReportSystem::ErrorCode::SelectSeqsResAndOverlapIncompatibilites, 
        "Options -selectCols and -selectSeqs are not allowed in combination of overlap options (-resoverlap || -seqoverlap)" }, 
        
    { ReportSystem::ErrorCode::OnlyOneSequencesSelectionMethodAllowed, 
        "Only one method to chose sequences can be applied: (-selectseqs || -clusters || -maxIdentity" }, 
        
    { ReportSystem::ErrorCode::CombinationAmongTrimmingMethods, 
        "Only one trimming method can be used at the same time, either automatic or manual" }, 
        
    { ReportSystem::ErrorCode::AutomathicMethodAndBlock, 
        "Combination between automatic methods and -block options is not allowed" }, 
        
    { ReportSystem::ErrorCode::WindowAndArgumentIncompatibilities, 
        "Combination between general and specific windows is not allowed" }, 
        
    { ReportSystem::ErrorCode::CombinationAmongThresholdsMethods, 
        "Combination among thresholds are not allowed" }, 
        
    { ReportSystem::ErrorCode::GeneralAndSpecificWindows, 
        "General window (-w) is not compatible with specific window options: (-cw, -gw, -sw)" }, 
    
    { ReportSystem::ErrorCode::StatisticsArgumentIncompatibilities, 
        "Parameter [tag] is not valid when statistics' parameters are defined" }, 
        
    { ReportSystem::ErrorCode::TrimmingMethodNeeded, 
        "Parameter [tag] can only be used with either an automatic or a manual method" }, 
        
    { ReportSystem::ErrorCode::ForceFileWithoutCompareDataset, 
        "You can not force the alignment selection without setting an alignment dataset to compare against" }, 
        
    { ReportSystem::ErrorCode::BacktranslationWithoutMainAlignment, 
        "You need to specify an input alignment (-in or -forcefile) to use a Coding Sequence File (-backtranslation) to apply the back translation method" }, 
        
    { ReportSystem::ErrorCode::NotAligned, 
        "The sequences in the input alignment [tag] should be aligned in order to use any trimming method or statistics" }, 
        
    { ReportSystem::ErrorCode::MatrixGivenWithNoMethodToUseIt, 
        "The Similarity Matrix can only be used with methods that use this matrix" }, 
        
    { ReportSystem::ErrorCode::SameNameOutput, 
        "The [tag] and [tag] files can't be the same" }, 
        
    { ReportSystem::ErrorCode::SequenceAndResiduesOverlapMutuallyNeeded, 
        "Sequence and residues overlap values are mutually needed. You only specified [tag]" }, 
        
    { ReportSystem::ErrorCode::OutFileNeededWhenPrintingStatistics, 
        "Out file must be specified in order to print any statistics" }, 
        
    { ReportSystem::ErrorCode::AlignmentTypesNotMatching, 
        "The alignments' datatypes are different. Check your dataset" }, 
        
    { ReportSystem::ErrorCode::BlocksizeTooBig, 
        "The block size value is too big. Please, choose another one smaller than residues number / 4. In this case, the limit is: [tag]" }, 
        
    { ReportSystem::ErrorCode::ParemeterOnlyOnBacktranslation, 
        "The [tag] parameter can be only set up with backtranslation functionality" }, 
        
    { ReportSystem::ErrorCode::ProteinAlignmentMustBeAligned, 
        "The input protein file has to be aligned to carry out the backtranslation process" }, 
                
    { ReportSystem::ErrorCode::BacktransAlignIsDNA, 
        "Check your Coding sequences file. It has been detected other kind of biological sequences" }, 
        
    { ReportSystem::ErrorCode::ImpossibleToGenerate, 
        "Impossible to generate [tag]" }, 
        
    { ReportSystem::ErrorCode::ImpossibleToProcessMatrix, 
        "It's impossible to process the Similarity Matrix" }, 
        
    { ReportSystem::ErrorCode::SelectOnlyAccepts, 
        "Option [tag] only accepts ranges from 0 to number of [tag] - 1" }, 
        
    { ReportSystem::ErrorCode::MoreClustersThanSequences, 
        "The number of clusters from the Alignment can not be larger than the number of sequences from that alignment" }, 
        
    { ReportSystem::ErrorCode::LeftBoundaryBiggerThanRightBoundary, 
        "Check your manually set left '[tag]' and right '[tag]' boundaries'" }, 
        
    { ReportSystem::ErrorCode::DifferentNumberOfSequencesInCompareset, 
        "The files to compare do not have the same number of sequences" }, 
    
    { ReportSystem::ErrorCode::DifferentSeqsNamesInCompareset, 
        "Sequences names differ in compareset files" }, 

    { ReportSystem::ErrorCode::CDScontainsProteinSequences, 
        "Check input CDS file. It seems to content protein residues." }, 
        
    { ReportSystem::ErrorCode::SequenceContainsGap, 
        "Sequence \"[tag]\" has, at least, one gap" }, 
        
    { ReportSystem::ErrorCode::SequenceNotMultipleOfThree, 
        "Sequence length \"[tag]\" is not multiple of 3 (length: [tag])" }, 
        
    { ReportSystem::ErrorCode::SequenceHasStopCodon, 
        "Sequence \"[tag]\" has stop codon \"[tag]\" (residue \"[tag]\") at position [tag] (length: [tag])" }, 
        
    { ReportSystem::ErrorCode::SequenceNotPresentInCDS, 
        "Sequence \"[tag]\" is not in CDS file." }, 
        
    { ReportSystem::ErrorCode::UnknownCharacter, 
        "The sequence \"[tag]\" has an unknown ([tag]) character." }, 
        
    { ReportSystem::ErrorCode::SequencesNotSameSize, 
        "The sequence \"[tag]\" ([tag]) does not have the same number of residues fixed by the alignment ([tag])" }, 
        
    { ReportSystem::ErrorCode::IncorrectSymbol, 
        "the symbol '[tag]' is incorrect" }, 
        
    { ReportSystem::ErrorCode::UndefinedSymbol, 
        "the symbol '[tag]' accesing the matrix is not defined in this object" }, 
        
    { ReportSystem::ErrorCode::ParameterNotFoundOrRepeated, 
        "Parameter \"[tag]\" not valid or repeated." }, 
        
    { ReportSystem::ErrorCode::SimilarityMatrixNotCompatibleWindow, 
        "The Similarity Matrix can only be used with general/similarity windows size." }, 
        
    { ReportSystem::ErrorCode::PossibleMissmatch, 
        "Possible (\") mismatch for comments" }, 
        
    { ReportSystem::ErrorCode::BracketsMissmatchFound, 
        "Brackets (]) mismatch found" }, 
        
    { ReportSystem::ErrorCode::UnalignedAlignmentToAlignedFormat, 
        "Sequences are not aligned. Format ([tag]) not compatible with unaligned sequences." }, 
};

void ReportSystem::PrintCodesAndMessages()
{
    switch(ReportSystem::Level)
    {
        case VerboseLevel::NONE:
            std::cout << "[VerboseLevel] None" << std::endl;
            break;
        case VerboseLevel::INFO:
            std::cout << "[VerboseLevel] Info" << std::endl;
            break;
        case VerboseLevel::WARNING:
            std::cout << "[VerboseLevel] Warning" << std::endl;
            break;
        case VerboseLevel::ERROR:
            std::cout << "[VerboseLevel] Error" << std::endl;
            break;
    }

    for (int i = 1; i < ReportSystem::InfoCode::__MAXINFO ; i++)
    {
//         std::cout << "Info Code:    " << std::setw(3) << std::setfill('0') << std::right << i << " -> ";
        ReportSystem::Report((ReportSystem::InfoCode)i);
    }
    
    for (int i = 1; i < ReportSystem::WarningCode::__MAXWARNING ; i++)
    {
//         std::cout << "Warning Code: " << std::setw(3) << std::setfill('0') << std::right << i << " -> ";
        ReportSystem::Report((ReportSystem::WarningCode)i);
    }  
    
    for (int i = 1; i < ReportSystem::ErrorCode::__MAXERROR ; i++)
    {
//         std::cerr << "Error Code:   " << std::setw(3) << std::setfill('0') << std::right << i << " -> ";
        ReportSystem::Report((ReportSystem::ErrorCode)i);
    }
}

void ReportSystem::Report(ReportSystem::ErrorCode message, std::string * vars)
{
    if (Level < VerboseLevel::ERROR)
    {
        if (vars != NULL)
            delete [] vars;
    }
    else 
    {
        if (vars == NULL)
        {
            std::cerr << "[ERROR "<< std::setw(3) << std::setfill('0') << message << "] " << ErrorMessages.at(message) << std::endl;
            return;
        }
        
        std::string s(ErrorMessages.at(message));

        std::string FindWord = "[tag]";
        
        int counter = 0;

        std::size_t index;
            while ((index = s.find(FindWord)) != std::string::npos)
                s.replace(index, FindWord.length(), vars[counter++]);

        std::cerr << "[ERROR "<< std::setw(3) << std::setfill('0') << message << "] " << s << std::endl;
        
        delete [] vars;
    }
}

void ReportSystem::Report(ReportSystem::ErrorCode message, char * vars)
{
    if (Level < VerboseLevel::ERROR) return;
    
    if (vars == NULL)
    {
        std::cerr << "[ERROR "<< std::setw(3) << std::setfill('0') << message << "] " << ErrorMessages.at(message) << std::endl;
        return;
    }
    
    std::string s(ErrorMessages.at(message));

    std::string FindWord = "[tag]";
    
    std::string Vars = vars;

    std::size_t index;
        while ((index = s.find(FindWord)) != std::string::npos)
            s.replace(index, FindWord.length(), Vars);

    std::cerr << "[ERROR "<< std::setw(3) << std::setfill('0') << message << "] " << s << std::endl;

}

void ReportSystem::Report(ReportSystem::WarningCode message, std::string * vars)
{
    if (Level < VerboseLevel::WARNING)
    {
        if (vars != NULL)
            delete [] vars;
    }
    else 
    {
        if (vars == NULL)
        {
            std::cout << "[WARNING "<< std::setw(3) << std::setfill('0') << message << "] " << WarningMessages.at(message) << std::endl;
            return;
        }
        
        std::string s(WarningMessages.at(message));

        std::string FindWord = "[tag]";
        
        int counter = 0;

        std::size_t index;
            while ((index = s.find(FindWord)) != std::string::npos)
                s.replace(index, FindWord.length(), vars[counter++]);

        std::cout << "[WARNING "<< std::setw(3) << std::setfill('0') << message << "] " << s << std::endl;
        
        delete [] vars;
    }
}

void ReportSystem::Report(ReportSystem::WarningCode message, char * vars)
{
    if (Level < VerboseLevel::WARNING) return;
    
    if (vars == NULL)
    {
        std::cout << "[WARNING "<< std::setw(3) << std::setfill('0') << message << "] " << WarningMessages.at(message) << std::endl;
        return;
    }
    
    std::string s(WarningMessages.at(message));

    std::string FindWord = "[tag]";
    
    std::string Vars = vars;

    std::size_t index;
        while ((index = s.find(FindWord)) != std::string::npos)
            s.replace(index, FindWord.length(), Vars);

    std::cout << "[WARNING "<< std::setw(3) << std::setfill('0') << message << "] " << s << std::endl;

}

void ReportSystem::Report(ReportSystem::InfoCode message, std::string * vars)
{
    if (Level < VerboseLevel::INFO)
    {
        if (vars != NULL)
            delete [] vars;
    }
    else 
    {
        if (vars == NULL)
        {
            std::cout << "[INFO "<< std::setw(3) << std::setfill('0') << message << "] " << InfoMessages.at(message) << std::endl;
            return;
        }
        
        std::string s(InfoMessages.at(message));

        std::string FindWord = "[tag]";
        
        int counter = 0;

        std::size_t index;
            while ((index = s.find(FindWord)) != std::string::npos)
                s.replace(index, FindWord.length(), vars[counter++]);

        std::cout << "[INFO "<< std::setw(3) << std::setfill('0') << message << "] " << s << std::endl;
        
        delete [] vars;
    }
}

void ReportSystem::Report(ReportSystem::InfoCode message, char * vars)
{
    if (Level < VerboseLevel::INFO) return;
    
    if (vars == NULL)
    {
        std::cout << "[INFO "<< std::setw(3) << std::setfill('0') << message << "] " << InfoMessages.at(message) << std::endl;
        return;
    }
    
    std::string s(InfoMessages.at(message));

    std::string FindWord = "[tag]";
    
    std::string Vars = vars;

    std::size_t index;
        while ((index = s.find(FindWord)) != std::string::npos)
            s.replace(index, FindWord.length(), Vars);

    std::cout << "[INFO "<< std::setw(3) << std::setfill('0') << message << "] " << s << std::endl;

}

void ReportSystem::Debug(std::string debugMessage)
{
    if (Level == VerboseLevel::DEBUG)
    {
        std::cout << debugMessage;
    }
}
