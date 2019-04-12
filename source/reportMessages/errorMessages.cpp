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

#include "reportsystem.h"

const std::map<ErrorCode, const char *> reporting::reportManager::ErrorMessages =
{
    {ErrorCode::AlignmentNotLoaded,
            "Alignment not loaded: \" [tag] \" Check the file's content"},

    {ErrorCode::NoFormatsSpecified,
            "You must specify at least one format after the '-formats' argument"},

    {ErrorCode::AlternativeMatrixNotRecognized,
            "Alternative matrix \" [tag] \" not recognized"},

    {ErrorCode::ReferenceFileNotLoaded,
            "Reference file \" [tag] \" not loaded"},

    {ErrorCode::GapThresholdOutOfRange,
            "The gap threshold value should be between 0 and 1"},

    {ErrorCode::GapThresholdNotRecognized,
            "The gap threshold value should be a positive real number between 0 and 1"},

    {ErrorCode::SimilarityThresholdOutOfRange,
            "The similarity threshold value should be between 0 and 1"},

    {ErrorCode::SimilarityThresholdNotRecognized,
            "The similarity threshold value should be a positive real number between 0 and 1"},

    {ErrorCode::ConsistencyThresholdOutOfRange,
            "The consistency threshold value should be between 0 and 1"},

    {ErrorCode::ConsistencyThresholdNotRecognized,
            "The consistency threshold value should be a positive real number between 0 and 1"},

    {ErrorCode::ConservationThresholdOutOfRange,
            "The similarity threshold value should be between 0 and 100"},

    {ErrorCode::ConservationThresholdNotRecognized,
            "The similarity threshold value should be a positive real number between 0 and 100"},

    {ErrorCode::ResidueOverlapOutOfRange,
            "The residue overlap value should be between 0 and 1"},

    {ErrorCode::ResidueOverlapNotRecognized,
            "The residue overlap value should be a positive real number between 0 and 100"},

    {ErrorCode::SequencesOverlapOutOfRange,
            "The sequences overlap value should be between 0 and 100"},

    {ErrorCode::SequencesOverlapNotRecognized,
            "The sequences overlap value should be a positive real number between 0 and 100"},

    {ErrorCode::MaxIdentityOutOfRange,
            "The max identity value should be between 0 and 1"},

    {ErrorCode::MaxIdentityNotRecognized,
            "The max identity value should be a positive real number between 0 and 1"},

    {ErrorCode::ClustersValueOutOfRange,
            "The clusters value should be greater than 0"},

    {ErrorCode::ClustersValueNotRecognized,
            "The clusters value should be a positive integer number greater than 0"},

    {ErrorCode::WindowValueOutOfRange,
            "The window value should be greater than 0"},

    {ErrorCode::WindowValueNotRecognized,
            "The window value should be a positive integer number greater than 0"},

    {ErrorCode::SelectSeqsNotRecognized,
            "Could not parse the -selectseqs ranges"},

    {ErrorCode::SelectColsNotRecognized,
            "Could not parse the -selectres ranges"},

    {ErrorCode::GapWindowValueOutOfRange,
            "The gap window value should be greater than 0"},

    {ErrorCode::GapWindowValueNotRecognized,
            "The gap window value should be a positive integer number grater than 0"},

    {ErrorCode::SimilarityWindowValueOutOfRange,
            "The similarity window value should be greater than 0"},

    {ErrorCode::SimilarityWindowValueNotRecognized,
            "The similarity window value should be a positive integer number grater than 0"},

    {ErrorCode::ConsistencyThresholdOutOfRange,
            "The consistency window value should be greater than 0"},

    {ErrorCode::ConsistencyThresholdNotRecognized,
            "The consistency window value should be a positive integer number grater than 0"},

    {ErrorCode::BlockSizeOutOfRange,
            "The consistency window value should be greater than 0"},

    {ErrorCode::BlockSizeNotRecognized,
            "The consistency window value should be a positive integer number grater than 0"},

    {ErrorCode::InFileComparisonStatistics,
            "Option -in not valid in combination with file comparision options: -sft || -sfc || -ct"},

    {ErrorCode::IncompatibleArguments,
            "Argument [tag] is not compatible with argument [tag]"},

    {ErrorCode::SelectSeqsResAndThresholdIncompatibilities,
            "Options -selectCols and -selectSeqs are not allowed in combination with threshold options: (-gt || -st || -ct || -cons)"},

    {ErrorCode::SelectSeqsResAndAutomathedMethodsIncompatibilities,
            "Options -selectCols and -selectSeqs are not allowed in combination with automated trimming options: (-nogaps || -noallgaps || -gappyout || -strict || -strictplus || -automated1)."},

    {ErrorCode::SelectSeqsResAndWindowIncompatibilities,
            "Options -selectCols and -selectSeqs are not allowed in combination with window options: (-w || -sw || -gw || -wc)"},

    {ErrorCode::SelectSeqsResAndOverlapIncompatibilites,
            "Options -selectCols and -selectSeqs are not allowed in combination of overlap options (-resoverlap || -seqoverlap)"},

    {ErrorCode::OnlyOneSequencesSelectionMethodAllowed,
            "Only one method to chose sequences can be applied: (-selectseqs || -clusters || -maxIdentity"},

    {ErrorCode::CombinationAmongTrimmingMethods,
            "Only one trimming method can be used at the same time, either automatic or manual. Manual argument: \"[tag]\" + \"[tag]\""},

    {ErrorCode::AutomathicMethodAndBlock,
            "Combination between automatic methods and -block options is not allowed"},

    {ErrorCode::WindowAndArgumentIncompatibilities,
            "Combination between general and specific windows is not allowed"},

    {ErrorCode::CombinationAmongThresholdsMethods,
            "Combination of -ct + -cons + either -gt or -st is not allowed"},

    {ErrorCode::GeneralAndSpecificWindows,
            "General window (-w) is not compatible with specific window options: (-cw, -gw, -sw)"},

    {ErrorCode::StatisticsArgumentIncompatibilities,
            "Parameter [tag] is not valid when statistics' parameters are defined"},

    {ErrorCode::TrimmingMethodNeeded,
            "Parameter [tag] can only be used with either an automatic or a manual method"},

    {ErrorCode::ForceFileWithoutCompareDataset,
            "You can not force the alignment selection without setting an alignment dataset to compare against"},

    {ErrorCode::BacktranslationWithoutMainAlignment,
            "You need to specify an input alignment (-in or -forcefile) to use a Coding Sequence File (-backtranslation) to apply the back translation method"},

    {ErrorCode::NotAligned,
            "The sequences in the input alignment [tag] should be aligned in order to use any trimming method or statistics"},

    {ErrorCode::MatrixGivenWithNoMethodToUseIt,
            "The Similarity Matrix can only be used with methods that use this matrix"},

    {ErrorCode::SameNameOutput,
            "The [tag] and [tag] files can't be the same"},

    {ErrorCode::SequenceAndResiduesOverlapMutuallyNeeded,
            "Sequence and residues overlap values are mutually needed. You only specified [tag]"},

    {ErrorCode::OutFileNeededWhenPrintingStatistics,
            "Out file must be specified in order to print any statistics"},

    {ErrorCode::AlignmentTypesNotMatching,
            "The alignments' datatypes are different. Check your dataset"},

    {ErrorCode::BlocksizeTooBig,
            "The block size value is too big. Please, choose another one smaller than residues number / 4. In this case, the limit is: [tag]"},

    {ErrorCode::ParemeterOnlyOnBacktranslation,
            "The [tag] parameter can be only set up with backtranslation functionality"},

    {ErrorCode::ProteinAlignmentMustBeAligned,
            "The input protein file has to be aligned to carry out the backtranslation process"},

    {ErrorCode::BacktransAlignIsDNA,
            "Check your Coding sequences file. It has been detected other kind of biological sequences"},

    {ErrorCode::ImpossibleToGenerate,
            "Impossible to generate [tag]"},

    {ErrorCode::ImpossibleToProcessMatrix,
            "It's impossible to process the Similarity Matrix"},

    {ErrorCode::SelectOnlyAccepts,
            "Option [tag] only accepts ranges from 0 to number of [tag] - 1"},

    {ErrorCode::MoreClustersThanSequences,
            "The number of clusters from the Alignment can not be larger than the number of sequences from that alignment"},

    {ErrorCode::LeftBoundaryBiggerThanRightBoundary,
            "Check your manually set left '[tag]' and right '[tag]' boundaries'"},

    {ErrorCode::DifferentNumberOfSequencesInCompareset,
            "The files to compare do not have the same number of sequences"},

    {ErrorCode::DifferentSeqsNamesInCompareset,
            "Sequences names differ in compareset files"},

    {ErrorCode::CDScontainsProteinSequences,
            "Check input CDS file. It seems to content protein residues."},

    {ErrorCode::SequenceContainsGap,
            "Sequence \"[tag]\" has, at least, one gap"},

    {ErrorCode::SequenceNotMultipleOfThree,
            "Sequence length \"[tag]\" is not multiple of 3 (length: [tag])"},

    {ErrorCode::SequenceHasStopCodon,
            "Sequence \"[tag]\" has stop codon \"[tag]\" (residue \"[tag]\") at position [tag] (length: [tag])"},

    {ErrorCode::SequenceNotPresentInCDS,
            "Sequence \"[tag]\" is not in CDS file"},

    {ErrorCode::UnknownCharacter,
            "The sequence \"[tag]\" has an unknown ([tag]) character"},

    {ErrorCode::SequencesNotSameSize,
            "The sequence \"[tag]\" ([tag]) does not have the same number of residues fixed by the alignment ([tag])"},

    {ErrorCode::IncorrectSymbol,
            "the symbol '[tag]' is incorrect"},

    {ErrorCode::UndefinedSymbol,
            "the symbol '[tag]' accesing the matrix is not defined in this object"},

    {ErrorCode::ParameterNotFoundOrRepeated,
            "Parameter \"[tag]\" not valid or repeated"},

    {ErrorCode::SimilarityMatrixNotCompatibleWindow,
            "The Similarity Matrix can only be used with general/similarity windows size"},

    {ErrorCode::PossibleMissmatch,
            "Possible (\") mismatch for comments"},

    {ErrorCode::BracketsMissmatchFound,
            "Brackets (]) mismatch found"},

    {ErrorCode::UnalignedAlignmentToAlignedFormat,
            "Sequences are not aligned. Format ([tag]) not compatible with unaligned sequences."},

    {ErrorCode::CantOpenFile,
            "File [tag] not found or impossible to open"},

    {ErrorCode::FileIsEmpty,
            "File [tag] is empty"},

    {ErrorCode::AlignmentFormatNotRecognized,
            "Alignment not loaded: \"[tag]\". Format couldn't be recognized."},

    {ErrorCode::OutputFormatNotRecognized,
            "Output format [tag] not recognized"},

    {ErrorCode::OnlyOneFormatOnConsoleOutput,
            "You must specify only one output format if you don't provide an output file pattern"},

    {ErrorCode::AlignmentNotSaved,
            "Alignment couldn't be saved on [tag] format"},

    {ErrorCode::VerboseLevelNotRecognized,
            "Verbose Level specified ([tag]) wasn't recognized. Current level is: [tag]"},

    {ErrorCode::NeedToSpecifyVerboseLevel,
            "Verbose Level has to be specified after the [tag] argument. Acceptable values are: 'error', 'warning', 'info', 'none' and their numerical equivalents '3', '2', '1' and '0'. Current level is [tag]"},

    {ErrorCode::NoReferenceSequenceForContig,
            "No reference sequence found for contig \"[tag]\""},

    {ErrorCode::SNPoutOfBounds,
            "SNP at positon \"[tag]\" in file \"[tag]\" cannot be applied as sequence has a length of \"[tag]\""},

    {ErrorCode::NoInputFile,
            "An MSA input file has to be provided"},

    {ErrorCode::ComparesetFailedAlignmentMissing,
            "Compareset couldn't be performed as some alignments are missing"},

    {ErrorCode::GapWindowTooBig,
            "Gap window size (-gw) provided is too big, please specify a window lesser than 1/4 of residues"},

    {ErrorCode::SimilarityWindowTooBig,
            "Similarity window size (-sw) provided is too big, please specify a window lesser than 1/4 of residues"},

    {ErrorCode::ConsistencyWindowTooBig,
            "Consistency window size (-cw) provided is too big, please specify a window lesser than 1/4 of residues"},

    {ErrorCode::WindowTooBig,
            "Window size (-w) provided is too big, please specify a window lesser than 1/4 of residues"},

    {ErrorCode::AlignmentIsEmpty,
            "Couldn't output any alignment as it doesn't contain any sequence."},

    {ErrorCode::AlignmentTypeIsUnknown,
            "Couldn't perform the trimming step, alignment type unkwnown. "},

    {ErrorCode::MultipleOutputFormatsSameName,
            "Multiple formats have been requested, but output pattern does not contain tags \"[extension]\" or \"[format]\".\n"
            "            Output files will be overwritten. Please, provide an output pattern containing tags to prevent overwritting."},

    {ErrorCode::MultipleInputs,
            "Multiple input files have been provided (Ej: \"-in X Y\"). This is not supported."},

    {ErrorCode::SomethingWentWrong_reportToDeveloper,
            "Something went wrong. Please, report to the developer with this message: \"[tag]\""},

    {ErrorCode::ReferenceNucleotideNotCorresponding,
             "Failed to apply SNP to \"[tag]\":\"[tag]\" at position \"[tag]\".\n"
             "\tCharacter expected at reference: \"[tag]\"\n"
             "\tCharacter at reference:          \"[tag]\""},

    {ErrorCode::OverwrittingSNP,
            "Overwriting SNP to \"[tag]\":\"[tag]\" at position \"[tag]\".\n"
            "\tCharacter found at reference:      \"[tag]\"\n"
            "\tCharacter found at destination:    \"[tag]\"\n"
            "\tCharacter applied to destination:  \"[tag]\""},

    {ErrorCode::MoreDonorsOnLineThanPresented,
            "There are more donors on the line than present on header. Please check your VCF: \"[tag]\""},

    {ErrorCode::MinQualityLesserThan0,
            "Min Quality should be equal or greater than 0."},

    {ErrorCode::MinQualityNotRecognized,
            "Min Quality provided not recognized"},

    {ErrorCode::MinCoverageLesserThan0,
            "Min Coverage should be equal or greater than 0."},

    {ErrorCode::MinCoverageNotRecognized,
            "Min Coverage provided not recognized"},

    {ErrorCode::OnlyValidWithVCF,
            "[tag] only valid in combination with -vcf argument."},

    {ErrorCode::TriedRenamingOutputPreventOverride,
            "[[tag]] Tried to prevent overriding file [tag] but found no luck adding suffixes.\n"
            "Overwritting file [tag]"
            "You should check your output folder."},
    {ErrorCode::AbsoluteAndRelativeGapThreshold,
            "Combination among absolute (-gat) and relative (-gt) gap thresholds is not allowed."},
    {ErrorCode::AbsoluteGapThresholdLessThanZero,
            "The absolute gap window value should be greater or equal to 0"},
    {ErrorCode::AbsoluteGapThresholdBiggerThanNumberOfSequences,
            "The absolute gap window value should be lesser than the number of sequences.\n"
            "\tNumber of sequences: \"[tag]\".\n"
            "\tAbsolute Gap Threshold provided: \"[tag]\"."},
    {ErrorCode::AbsoluteGapThresholdNotRecognized,
            "Absolute gap (-gat) not recognized"},

    {ErrorCode::MoreThanOneAutomatedMethod,
            "More than one automated method has been requested."},

    {ErrorCode::ForceSelectAndInArgumentsProvided,
            "Arguments -in <x> and -forceselect <x> are incompatible"},

    {ErrorCode::ComparesetAndInArgumentsProvided,
            "Arguments -in <x> and -compareset <x> are incompatible"},

};