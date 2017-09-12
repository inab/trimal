#include <iomanip>
#include "../include/verbosemanager.h"
#include <iostream>

// Out of line initializations.

VerboseManager::VerboseLevel VerboseManager::Level = VerboseLevel::Errors;

const std::map<VerboseManager::InfoCode, const char *> VerboseManager::InfoMessages = 
{
    { VerboseManager::InfoCode::CuttingSequence, 
        "WARNING: Cutting sequence \"[ŧag]\" at first appearance of stop codon \"[tag]\" (residue \"[tag]\") at position [tag] (length: [tag] \")" }, 
};
    
const std::map<VerboseManager::WarningCode, const char *> VerboseManager::WarningMessages = 
{
    { VerboseManager::WarningCode::RemovingOnlyGapsColumn, 
        "WARNING: Removing column '[tag]' composed only by gaps" }, 
        
    { VerboseManager::WarningCode::KeepingOnlyGapsColumn, 
        "WARNING: Keeping column '[tag]' composed only by gaps" }, 
        
    { VerboseManager::WarningCode::SequenceWillBeCutted, 
        "WARNING: Sequence \"[tag]\" will be cutted at position [tag] (length:[tag])" }, 
        
    { VerboseManager::WarningCode::IncludingIndeterminationSymbols,
        "WARNING: Sequence \"[tag]\" has some indetermination symbols 'X' at the end of sequence. They will be included in the final newAlignment." },
        
    { VerboseManager::WarningCode::LessNucleotidesThanExpected,
        "WARNING: Sequence \"[tag]\" has less nucleotides ([tag]) than expected ([tag]). It will be added N's to complete the sequence" }
};
    
const std::map<VerboseManager::ErrorCode, const char *> VerboseManager::ErrorMessages = 
{
    { VerboseManager::ErrorCode::AlignmentNotLoaded, 
        "Alignment not loaded: \" [tag] \" Check the file's content" }, 
        
    { VerboseManager::ErrorCode::NoFormatsSpecified, 
        "You must specify at least one format after the '-formats' argument" }, 
        
    { VerboseManager::ErrorCode::AlternativeMatrixNotRecognized, 
        "Alternative matrix \" [tag] \" not recognized" }, 
        
    { VerboseManager::ErrorCode::ReferenceFileNotLoaded, 
        "Reference file \" [tag] \" not loaded" }, 
        
    { VerboseManager::ErrorCode::GapThresholdOutOfRange, 
        "The gap threshold value should be between 0 and 1" }, 
        
    { VerboseManager::ErrorCode::GapThresholdNotRecognized, 
        "The gap threshold value should be a positive real number between 0 and 1" }, 
        
    { VerboseManager::ErrorCode::SimilarityThresholdOutOfRange, 
        "The similarity threshold value should be between 0 and 1" }, 
        
    { VerboseManager::ErrorCode::SimilarityThresholdNotRecognized, 
        "The similarity threshold value should be a positive real number between 0 and 1" }, 
        
    { VerboseManager::ErrorCode::ConsistencyThresholdOutOfRange, 
        "The consistency threshold value should be between 0 and 1" }, 
        
    { VerboseManager::ErrorCode::ConsistencyThresholdNotRecognized, 
        "The consistency threshold value should be a positive real number between 0 and 1" }, 
        
    { VerboseManager::ErrorCode::ConservationThresholdOutOfRange, 
        "The conservation threshold value should be between 0 and 100" }, 
        
    { VerboseManager::ErrorCode::ConservationThresholdNotRecognized, 
        "The conservation threshold value should be a positive real number between 0 and 100" }, 
        
    { VerboseManager::ErrorCode::ResidueOverlapOutOfRange, 
        "The residue overlap value should be between 0 and 100" }, 
        
    { VerboseManager::ErrorCode::ResidueOverlapNotRecognized, 
        "The residue overlap value should be a positive real number between 0 and 100" }, 
        
    { VerboseManager::ErrorCode::SequencesOverlapOutOfRange, 
        "The sequences overlap value should be between 0 and 100" }, 
        
    { VerboseManager::ErrorCode::SequencesOverlapNotRecognized, 
        "The sequences overlap value should be a positive real number between 0 and 100" }, 
        
    { VerboseManager::ErrorCode::MaxIdentityOutOfRange, 
        "The max identity value should be between 0 and 1" }, 
        
    { VerboseManager::ErrorCode::MaxIdentityNotRecognized, 
        "The max identity value should be a positive real number between 0 and 1" }, 
        
    { VerboseManager::ErrorCode::ClustersValueOutOfRange, 
        "The clusters value should be greater than 0" }, 
        
    { VerboseManager::ErrorCode::ClustersValueNotRecognized, 
        "The clusters value should be a positive integer number greater than 0" }, 
        
    { VerboseManager::ErrorCode::WindowValueOutOfRange, 
        "The window value should be greater than 0" }, 
        
    { VerboseManager::ErrorCode::WindowValueNotRecognized, 
        "The window value should be a positive integer number greater than 0" }, 
        
    { VerboseManager::ErrorCode::SelectSeqsNotRecognized, 
        "Could not parse the -selectseqs ranges" }, 
        
    { VerboseManager::ErrorCode::SelectColsNotRecognized, 
        "Could not parse the -selectres ranges" }, 
        
    { VerboseManager::ErrorCode::GapWindowValueOutOfRange, 
        "The gap window value should be greater than 0" }, 
        
    { VerboseManager::ErrorCode::GapWindowValueNotRecognized, 
        "The gap window value should be a positive integer number grater than 0" }, 
        
    { VerboseManager::ErrorCode::SimilarityWindowValueOutOfRange, 
        "The similarity window value should be greater than 0" }, 
        
    { VerboseManager::ErrorCode::SimilarityWindowValueNotRecognized, 
        "The similarity window value should be a positive integer number grater than 0" }, 
        
    { VerboseManager::ErrorCode::ConsistencyThresholdOutOfRange, 
        "The consistency window value should be greater than 0" }, 
        
    { VerboseManager::ErrorCode::ConsistencyThresholdNotRecognized, 
        "The consistency window value should be a positive integer number grater than 0" }, 
        
    { VerboseManager::ErrorCode::BlockSizeOutOfRange, 
        "The consistency window value should be greater than 0" }, 
        
    { VerboseManager::ErrorCode::BlockSizeNotRecognized, 
        "The consistency window value should be a positive integer number grater than 0" }, 
        
    { VerboseManager::ErrorCode::InFileComparisonStatistics, 
        "Option -in not valid in combination with file comparision options: -sft || -sfc || -ct" }, 
        
    { VerboseManager::ErrorCode::IncompatibleArguments, 
        "Argument [tag] is not compatible with argument [tag]" }, 
        
    { VerboseManager::ErrorCode::SelectSeqsResAndThresholdIncompatibilities, 
        "Options -selectCols and -selectSeqs are not allowed in combination with threshold options: (-gt || -st || -ct || -cons)" }, 
    
    { VerboseManager::ErrorCode::SelectSeqsResAndAutomathedMethodsIncompatibilities, 
        "Options -selectCols and -selectSeqs are not allowed in combination with automated trimming options: (-nogaps || -noallgaps || -gappyout || -strict || -strictplus || -automated1)." }, 
        
    { VerboseManager::ErrorCode::SelectSeqsResAndWindowIncompatibilities, 
        "Options -selectCols and -selectSeqs are not allowed in combination with window options: (-w || -sw || -gw || -wc)" }, 
        
    { VerboseManager::ErrorCode::SelectSeqsResAndOverlapIncompatibilites, 
        "Options -selectCols and -selectSeqs are not allowed in combination of overlap options (-resoverlap || -seqoverlap)" }, 
        
    { VerboseManager::ErrorCode::OnlyOneSequencesSelectionMethodAllowed, 
        "Only one method to chose sequences can be applied: (-selectseqs || -clusters || -maxIdentity" }, 
        
    { VerboseManager::ErrorCode::CombinationAmongTrimmingMethods, 
        "Only one trimming method can be used at the same time, either automatic or manual" }, 
        
    { VerboseManager::ErrorCode::AutomathicMethodAndBlock, 
        "Combination between automatic methods and -block options is not allowed" }, 
        
    { VerboseManager::ErrorCode::WindowAndArgumentIncompatibilities, 
        "Combination between general and specific windows is not allowed" }, 
        
    { VerboseManager::ErrorCode::CombinationAmongThresholdsMethods, 
        "Combination among thresholds are not allowed" }, 
        
    { VerboseManager::ErrorCode::GeneralAndSpecificWindows, 
        "General window (-w) is not compatible with specific window options: (-cw, -gw, -sw)" }, 
    
    { VerboseManager::ErrorCode::StatisticsArgumentIncompatibilities, 
        "Parameter [tag] is not valid when statistics' parameters are defined" }, 
        
    { VerboseManager::ErrorCode::TrimmingMethodNeeded, 
        "Parameter [tag] can only be used with either an automatic or a manual method" }, 
        
    { VerboseManager::ErrorCode::ForceFileWithoutCompareDataset, 
        "You can not force the alignment selection without setting an alignment dataset to compare against" }, 
        
    { VerboseManager::ErrorCode::BacktranslationWithoutMainAlignment, 
        "You need to specify an input alignment (-in or -forcefile) to use a Coding Sequence File (-backtranslation) to apply the back translation method" }, 
        
    { VerboseManager::ErrorCode::NotAligned, 
        "The sequences in the input alignment [tag] should be aligned in order to use any trimming method or statistics" }, 
        
    { VerboseManager::ErrorCode::MatrixGivenWithNoMethodToUseIt, 
        "The Similarity Matrix can only be used with methods that use this matrix" }, 
        
    { VerboseManager::ErrorCode::SameNameOutput, 
        "The [tag] and [tag] files can't be the same" }, 
        
    { VerboseManager::ErrorCode::SequenceAndResiduesOverlapMutuallyNeeded, 
        "Sequence and residues overlap values are mutually needed. You only specified [tag]" }, 
        
    { VerboseManager::ErrorCode::OutFileNeededWhenPrintingStatistics, 
        "Out file must be specified in order to print any statistics" }, 
        
    { VerboseManager::ErrorCode::AlignmentTypesNotMatching, 
        "The alignments' datatypes are different. Check your dataset" }, 
        
    { VerboseManager::ErrorCode::BlocksizeTooBig, 
        "The block size value is too big. Please, choose another one smaller than residues number / 4. In this case, the limit is: [tag]" }, 
        
    { VerboseManager::ErrorCode::ParemeterOnlyOnBacktranslation, 
        "The [tag] parameter can be only set up with backtranslation functionality" }, 
        
    { VerboseManager::ErrorCode::ProteinAlignmentMustBeAligned, 
        "The input protein file has to be aligned to carry out the backtranslation process" }, 
                
    { VerboseManager::ErrorCode::BacktransAlignIsDNA, 
        "Check your Coding sequences file. It has been detected other kind of biological sequences" }, 
                        
    { VerboseManager::ErrorCode::WindowSizeCompareset, 
        "Try with specific comparison file window value. Parameter -cw" }, 
        
    { VerboseManager::ErrorCode::ImpossibleToGenerate, 
        "Impossible to generate [tag]" }, 
        
    { VerboseManager::ErrorCode::ImpossibleToProcessMatrix, 
        "It's impossible to process the Similarity Matrix" }, 
        
    { VerboseManager::ErrorCode::SelectOnlyAccepts, 
        "Option [tag] only accepts ranges from 0 to number of [tag] - 1" }, 
        
    { VerboseManager::ErrorCode::MoreClustersThanSequences, 
        "The number of clusters from the Alignment can not be larger than the number of sequences from that alignment" }, 
        
    { VerboseManager::ErrorCode::LeftBoundaryBiggerThanRightBoundary, 
        "Check your manually set left '[tag]' and right '[tag]' boundaries'" }, 
        
    { VerboseManager::ErrorCode::DifferentNumberOfSequencesInCompareset, 
        "The files to compare do not have the same number of sequences" }, 
    
    { VerboseManager::ErrorCode::DifferentSeqsNamesInCompareset, 
        "Sequences names differ in compareset files" }, 

    { VerboseManager::ErrorCode::CDScontainsProteinSequences, 
        "Check input CDS file. It seems to content protein residues." }, 
        
    { VerboseManager::ErrorCode::SequenceContainsGap, 
        "Sequence \"[tag]\" has, at least, one gap" }, 
        
    { VerboseManager::ErrorCode::SequenceNotMultipleOfThree, 
        "Sequence length \"[tag]\" is not multiple of 3 (length: [tag])" }, 
        
    { VerboseManager::ErrorCode::SequenceHasStopCodon, 
        "Sequence \"[tag]\" has stop codon \"[tag]\" (residue \"[tag]\") at position [tag] (length: [tag])" }, 
        
    { VerboseManager::ErrorCode::SequenceNotPresentInCDS, 
        "Sequence \"[tag]\" is not in CDS file." }, 
        
    { VerboseManager::ErrorCode::UnknownCharacter, 
        "The sequence \"[tag]\" has an unknown ([tag]) character." }, 
        
    { VerboseManager::ErrorCode::SequencesNotSameSize, 
        "The sequence \"[tag]\" ([tag]) does not have the same number of residues fixed by the alignment ([tag])" }, 
        
    { VerboseManager::ErrorCode::IncorrectSymbol, 
        "the symbol '[tag]' is incorrect" }, 
        
    { VerboseManager::ErrorCode::UndefinedSymbol, 
        "the symbol '[tag]' accesing the matrix is not defined in this object" }, 
        
    { VerboseManager::ErrorCode::ParameterNotFoundOrRepeated, 
        "Parameter \"[tag]\" not valid or repeated." }, 
        
    { VerboseManager::ErrorCode::SimilarityMatrixNotCompatibleWindow, 
        "The Similarity Matrix can only be used with general/similarity windows size." }, 
        
    { VerboseManager::ErrorCode::PossibleMissmatch, 
        "Possible (\") mismatch for comments" }, 
        
    { VerboseManager::ErrorCode::BracketsMissmatchFound, 
        "Brackets (]) mismatch found" }, 
};

void VerboseManager::PrintCodesAndMessages()
{
    switch(VerboseManager::Level)
    {
        case VerboseLevel::None:
            std::cout << "[VerboseLevel] None" << std::endl;
        case VerboseLevel::Info:
            std::cout << "[VerboseLevel] Info" << std::endl;
        case VerboseLevel::Warnings:
            std::cout << "[VerboseLevel] Warning" << std::endl;
        case VerboseLevel::Errors:
            std::cout << "[VerboseLevel] Error" << std::endl;
    }

    for (int i = 1; i < VerboseManager::InfoCode::__MAXINFO ; i++)
    {
        std::cout << "Info Code:    " << std::setw(3) << std::setfill('0') << std::right << i << " -> ";
        VerboseManager::Report((VerboseManager::InfoCode)i);
    }
    
    for (int i = 1; i < VerboseManager::WarningCode::__MAXWARNING ; i++)
    {
        std::cout << "Warning Code: " << std::setw(3) << std::setfill('0') << std::right << i << " -> ";
        VerboseManager::Report((VerboseManager::WarningCode)i);
    }  
    
    for (int i = 1; i < VerboseManager::ErrorCode::__MAXERROR ; i++)
    {
        std::cerr << "Error Code:   " << std::setw(3) << std::setfill('0') << std::right << i << " -> ";
        VerboseManager::Report((VerboseManager::ErrorCode)i);
    }
}

void VerboseManager::Report(VerboseManager::ErrorCode message, std::string * vars)
{
    if (Level < VerboseLevel::Errors)
    {
        if (vars != NULL)
            delete [] vars;
    }
    else 
    {
        if (vars == NULL)
        {
            std::cerr << "[ERROR] " << ErrorMessages.at(message) << std::endl;
            return;
        }
        
        std::string s(ErrorMessages.at(message));

        std::string FindWord = "[tag]";
        
        int counter = 0;

        std::size_t index;
            while ((index = s.find(FindWord)) != std::string::npos)
                s.replace(index, FindWord.length(), vars[counter++]);

        std::cerr << "[ERROR] " << s << std::endl;
        
        delete [] vars;
    }
}

void VerboseManager::Report(VerboseManager::ErrorCode message, char * vars)
{
    if (Level < VerboseLevel::Errors) return;
    
    if (vars == NULL)
    {
        std::cerr << "[ERROR] " << ErrorMessages.at(message) << std::endl;
        return;
    }
    
    std::string s(ErrorMessages.at(message));

    std::string FindWord = "[tag]";
    
    std::string Vars = vars;

    std::size_t index;
        while ((index = s.find(FindWord)) != std::string::npos)
            s.replace(index, FindWord.length(), Vars);

    std::cerr << "[ERROR] " << s << std::endl;

}

void VerboseManager::Report(VerboseManager::WarningCode message, std::string * vars)
{
    if (Level < VerboseLevel::Warnings)
    {
        if (vars != NULL)
            delete [] vars;
    }
    else 
    {
        if (vars == NULL)
        {
            std::cout << "[WARNING] " << WarningMessages.at(message) << std::endl;
            return;
        }
        
        std::string s(WarningMessages.at(message));

        std::string FindWord = "[tag]";
        
        int counter = 0;

        std::size_t index;
            while ((index = s.find(FindWord)) != std::string::npos)
                s.replace(index, FindWord.length(), vars[counter++]);

        std::cout << "[WARNING] " << s << std::endl;
        
        delete [] vars;
    }
}

void VerboseManager::Report(VerboseManager::WarningCode message, char * vars)
{
    if (Level < VerboseLevel::Warnings) return;
    
    if (vars == NULL)
    {
        std::cout << "[WARNING] " << WarningMessages.at(message) << std::endl;
        return;
    }
    
    std::string s(WarningMessages.at(message));

    std::string FindWord = "[tag]";
    
    std::string Vars = vars;

    std::size_t index;
        while ((index = s.find(FindWord)) != std::string::npos)
            s.replace(index, FindWord.length(), Vars);

    std::cout << "[WARNING] " << s << std::endl;

}

void VerboseManager::Report(VerboseManager::InfoCode message, std::string * vars)
{
    if (Level < VerboseLevel::Info)
    {
        if (vars != NULL)
            delete [] vars;
    }
    else 
    {
        if (vars == NULL)
        {
            std::cout << "[INFO] " << InfoMessages.at(message) << std::endl;
            return;
        }
        
        std::string s(InfoMessages.at(message));

        std::string FindWord = "[tag]";
        
        int counter = 0;

        std::size_t index;
            while ((index = s.find(FindWord)) != std::string::npos)
                s.replace(index, FindWord.length(), vars[counter++]);

        std::cout << "[INFO] " << s << std::endl;
        
        delete [] vars;
    }
}

void VerboseManager::Report(VerboseManager::InfoCode message, char * vars)
{
    if (Level < VerboseLevel::Info) return;
    
    if (vars == NULL)
    {
        std::cout << "[INFO] " << InfoMessages.at(message) << std::endl;
        return;
    }
    
    std::string s(InfoMessages.at(message));

    std::string FindWord = "[tag]";
    
    std::string Vars = vars;

    std::size_t index;
        while ((index = s.find(FindWord)) != std::string::npos)
            s.replace(index, FindWord.length(), Vars);

    std::cout << "[INFO] " << s << std::endl;

}
