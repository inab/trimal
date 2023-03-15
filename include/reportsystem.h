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

#ifndef VERBOSEMANAGER_H
#define VERBOSEMANAGER_H

#include <iostream>
#include <cstring>
#include <iomanip>
#include <vector>
#include <array>
#include <map>


/** 
 * \brief VerboseLevel used to report messages. 
 */
enum VerboseLevel {

    /// 1 = Info, warning and error messages
    INFO        = 1,
    /// 2 = Error and warning messages
    WARNING    = 2,
    /// 3 = Only error messages
    ERROR      = 3,
    /// 4 = No output messages
    NONE        = 4,
};

enum ErrorCode {

    SomethingWentWrong_reportToDeveloper                = 0,

    AlignmentNotLoaded                                  = 1,

    NoFormatsSpecified                                  = 2,

    AlternativeMatrixNotRecognized                      = 3,

    ReferenceFileNotLoaded                              = 4,

    GapThresholdOutOfRange                              = 5,
    GapThresholdNotRecognized                           = 6,

    SimilarityThresholdOutOfRange                       = 7,
    SimilarityThresholdNotRecognized                    = 8,

    ConsistencyThresholdOutOfRange                      = 9,
    ConsistencyThresholdNotRecognized                   = 10,

    ConservationThresholdOutOfRange                     = 11,
    ConservationThresholdNotRecognized                  = 12,

    ResidueOverlapOutOfRange                            = 13,
    ResidueOverlapNotRecognized                         = 14,

    SequencesOverlapOutOfRange                          = 15,
    SequencesOverlapNotRecognized                       = 16,

    MaxIdentityOutOfRange                               = 17,
    MaxIdentityNotRecognized                            = 18,

    ClustersValueOutOfRange                             = 19,
    ClustersValueNotRecognized                          = 20,

    WindowValueOutOfRange                               = 21,
    WindowValueNotRecognized                            = 22,

    SelectSeqsNotRecognized                             = 23,
    SelectColsNotRecognized                             = 24,

    GapWindowValueOutOfRange                            = 25,
    GapWindowValueNotRecognized                         = 26,

    SimilarityWindowValueOutOfRange                     = 27,
    SimilarityWindowValueNotRecognized                  = 28,

    ConsistencyWindowValueOutOfRange                    = 27,
    ConsistencyWindowValueNotRecognized                 = 28,

    BlockSizeOutOfRange                                 = 29,
    BlockSizeNotRecognized                              = 30,

    InFileComparisonStatistics                          = 31,

    IncompatibleArguments                               = 32,

    SelectSeqsResAndThresholdIncompatibilities          = 33,
    SelectSeqsResAndAutomathedMethodsIncompatibilities  = 34,
    SelectSeqsResAndWindowIncompatibilities             = 35,
    SelectSeqsResAndOverlapIncompatibilites             = 36,

    OnlyOneSequencesSelectionMethodAllowed              = 37,

    CombinationAmongTrimmingMethods                     = 38,
    AutomathicMethodAndBlock                            = 39,

    WindowAndArgumentIncompatibilities                  = 40,
    CombinationAmongThresholdsMethods                   = 41,
    GeneralAndSpecificWindows                           = 42,
    StatisticsArgumentIncompatibilities                 = 43,

    TrimmingMethodNeeded                                = 44,
    ForceFileWithoutCompareDataset                      = 45,

    BacktranslationWithoutMainAlignment                 = 46,

    NotAligned                                          = 47,

    MatrixGivenWithNoMethodToUseIt                      = 48,

    SameNameOutput                                      = 49,

    SequenceAndResiduesOverlapMutuallyNeeded            = 50,

    OutFileNeededWhenPrintingStatistics                 = 51,

    AlignmentTypesNotMatching                           = 52,

    BlocksizeTooBig                                     = 53,

    ParemeterOnlyOnBacktranslation                      = 54,

    ProteinAlignmentMustBeAligned                       = 55,

    BacktransAlignIsDNA                                 = 56,



    ImpossibleToGenerate                                = 57,

    ImpossibleToProcessMatrix                           = 58,

    SelectOnlyAccepts                                   = 59,

    MoreClustersThanSequences                           = 60,

    LeftBoundaryBiggerThanRightBoundary                 = 61,

    DifferentNumberOfSequencesInCompareset              = 62,

    DifferentSeqsNamesInCompareset                      = 63,

    CDScontainsProteinSequences                         = 64,

    SequenceContainsGap                                 = 65,

    SequenceNotMultipleOfThree                          = 66,

    SequenceHasStopCodon                                = 67,

    SequenceNotPresentInCDS                             = 68,

    UnknownCharacter                                    = 69,

    SequencesNotSameSize                                = 70,

    IncorrectSymbol                                     = 71,

    UndefinedSymbol                                     = 72,

    ParameterNotFoundOrRepeated                         = 73,

    SimilarityMatrixNotCompatibleWindow                 = 74,

    PossibleMissmatch                                   = 75,

    BracketsMissmatchFound                              = 76,

    UnalignedAlignmentToAlignedFormat                   = 77,
    
    CantOpenFile                                        = 78,
    
    FileIsEmpty                                         = 79,
    
    AlignmentFormatNotRecognized                        = 80,
    
    OutputFormatNotRecognized                           = 81,
    
    OnlyOneFormatOnConsoleOutput                        = 82,
    
    AlignmentNotSaved                                   = 83,
    
    VerboseLevelNotRecognized                           = 84,
    
    NeedToSpecifyVerboseLevel                           = 85,
    
    NoReferenceSequenceForContig                        = 86,
    
    SNPoutOfBounds                                      = 87,

    NoInputFile                                         = 88,

    ComparesetFailedAlignmentMissing                    = 89,

    GapWindowTooBig                                     = 90,
    SimilarityWindowTooBig                              = 91,
    ConsistencyWindowTooBig                             = 92,
    WindowTooBig                                        = 93,

    AlignmentIsEmpty                                    = 94,

    AlignmentTypeIsUnknown                              = 95,

    MultipleOutputFormatsSameName                       = 96,
    MultipleInputs                                      = 97, // Not in use

    ReferenceNucleotideNotCorresponding                 = 98,

    OverwrittingSNP                                     = 99,


    MoreDonorsOnLineThanPresented                       = 100,

    MinQualityLesserThan0                               = 101,
    MinQualityNotRecognized                             = 102,

    MinCoverageLesserThan0                              = 103,
    MinCoverageNotRecognized                            = 104,

    OnlyValidWithVCF                                    = 105,

    TriedRenamingOutputPreventOverride                  = 106,

    AbsoluteAndRelativeGapThreshold                     = 107,
    AbsoluteGapThresholdLessThanZero                    = 108,
    AbsoluteGapThresholdBiggerThanNumberOfSequences     = 109,
    AbsoluteGapThresholdNotRecognized                   = 110,

    MoreThanOneAutomatedMethod                          = 111,


    ForceSelectAndInArgumentsProvided                   = 112,
    ComparesetAndInArgumentsProvided                    = 113,

    NoResidueSequences                                  = 114,


    __MAXERROR,
};

enum WarningCode {
    RemovingOnlyGapsSequence                    = 1,

    KeepingOnlyGapsSequence                     = 2,

    SequenceWillBeCut                           = 3,

    IncludingIndeterminationSymbols             = 4,

    LessNucleotidesThanExpected                 = 5,

    HeaderWillBeCut                             = 6,
    
    DonorAlreadyAdded                           = 7,

    SNPAlreadApplied                            = 8,

    OverwrittingFile                            = 9,

    RenamingOutputPreventOverride               = 10,


    __MAXWARNING
};

enum InfoCode {
    CuttingSequence                             = 1,

    WindowSizeCompareset                        = 2,

    AddingSNP                                   = 3,

    RemovingDuplicateSequences                  = 4,


    __MAXINFO
};

/// \brief Internal classes to handle reporting to user in several ways.\n
/// The reporting system is made so a developer that wants to use the system
///     doesn't need to enter to the namespace, but use the global variable #debug
namespace reporting {

/**
 \brief Internal class used by reporting::reportManager \n
 This class serves as a proxy for message reporting. \n
 It will only output to std::cout if
  reporting::reportWrapper::CanReport it's true
 */
class reportWrapper {
    /// \brief Variable that specifies if the object is able to output to cout or not.
    bool CanReport;
    
    public:
    
    /// \brief Constructor
    explicit reportWrapper(bool CanReport) : CanReport{CanReport} { }
    
    /// \brief Overloaded operator that allows to use this object as
    /// it was an iostream, as cout.
    template <typename T>
    reportWrapper &operator<<(const T &a) {
        if (CanReport)
            std::cout<<a;
        return *this;
    }
    
    /// \brief Overloaded operator that allows to use this object as
    /// it was an iostream, as cout.
    reportWrapper &operator<<(std::ostream& (*pf) (std::ostream&)) {
        if (CanReport)
            std::cout<<pf;
        return *this;
    }
};

/**
 \brief Class that allows us to centralize all the reporting messages that
  should be used to inform the user of Errors, Warnings and Info.\n
 The object can also be used as a substitute to cout for temporal messages using
  the << overloaded operator or the \link reportManager::log log \endlink method,
  and the messages will be behind the IsDebug variable,
  which should be set to false in Release.\n
 This allows us to protect the user from receiving Debug/Developing messages
  that shouldn't bother them, in case any of them is not removed by error.
 \warning <b> THIS CLASS SHOULDN'T BE USED DIRECTLY </b> \n
 Use instanced global variable \link debug debug \endlink instead.
 */
class reportManager
{
private:
    
    /** \brief Object that will be returned by
         \link reporting::reportManager::log log \endlink
         if it allows to output the message
     */
    reportWrapper canReport = reportWrapper(true);
    /** \brief Object that will be returned by
     *   \link reporting::reportManager::log log \endlink
     *   if it doesn't allow to output the message
     * */
    reportWrapper cantReport = reportWrapper(false);

    static const std::map<InfoCode, const char *>       InfoMessages ;
    static const std::map<WarningCode, const char *>    WarningMessages ;
    static const std::map<ErrorCode, const char *>      ErrorMessages ;
    
public:
    /** \brief Level of Verbosity.\n
     *   The report system won't output messages
     *    that are lower than the current level */
    VerboseLevel Level = VerboseLevel::INFO;
    
    /** \brief Protection variable.
      This variable should be set to true while in development,
      and false on Release. */
    bool IsDebug = false;

    /**
     \brief Method to print all Info, Warning and Error codes
      and their respective message.\n
     This method is useful to check whether all codes
      have a message linked to them.
     It will use the current VerboseLevel of the object,
      which defaults to ERROR level.
     */
    void PrintCodesAndMessages();

    /**
     \brief Method to report an Error. \n
     It will be displayed if
      \link reporting::reportManager::Level Level \endlink
      is equal or higher to VerboseLevel::ERROR
     \param message
      Code to report.
     \param vars
      Array of strings that will replace the '[tags]' in the message.\n
      The number of elements in the array must be the same as
      [tag] occurrences on the message.\n
     <b> The method will take care of destroying the array </b>
     */
    void report(ErrorCode message, std::string * vars = nullptr);

    /**
     \brief Method to report an Error.\n
     It will be displayed if
      \link reporting::reportManager::Level Level \endlink
      is equal or higher to VerboseLevel::ERROR
     \param message
      Code to report.
     \param vars
      Array of chars that will replace the first occurrence of
      '[tags]' in the message. \n
     <b> The method wont take care of destroying the array </b>\n
      This allows to reuse the parameter passed.
     */
    void report(ErrorCode message, const char *vars);
    /**
     \brief Method to report a Warning. \n
     It will be displayed if
      \link reporting::reportManager::Level Level \endlink
      is equal or higher to VerboseLevel::WARNING
     \param message
      Code to report.
     \param vars
      Array of strings that will replace the '[tags]' in the message.\n
      The number of elements in the array must be the same as
      [tag] occurrences on the message.\n
     <b> The method will take care of destroying the array </b>
     */
    void report(WarningCode message, std::string * vars = nullptr);
    /**
     \brief Method to report a Warning.\n
     It will be displayed if
      \link reporting::reportManager::Level Level \endlink
      is equal or higher to VerboseLevel::WARNING
     \param message
      Code to report.
     \param vars
      Array of chars that will replace the first occurrence of
      '[tags]' in the message. \n
     <b> The method wont take care of destroying the array </b>\n
      This allows to reuse the parameter passed.
     */
    void report(WarningCode message, const char *vars);
    /**
     \brief Method to report an Info message. \n
     It will be displayed if
      \link reporting::reportManager::Level Level \endlink
      is equal or higher to VerboseLevel::INFO
     \param message
      Code to report.
     \param vars
      Array of strings that will replace the '[tags]' in the message.\n
      The number of elements in the array must be the same as
      [tag] occurrences on the message.\n
     <b> The method will take care of destroying the array </b>
     */
    void report(InfoCode message, std::string * vars = nullptr);
    /**
     \brief Method to report an Info message.\n
     It will be displayed if
      \link reporting::reportManager::Level Level \endlink
      is equal or higher to VerboseLevel::WARNING
     \param message
      Code to report.
     \param vars
      Array of chars that will replace the first occurrence of
      '[tags]' in the message. \n
     <b> The method wont take care of destroying the array </b>\n
      This allows to reuse the parameter passed.
     */
    void report(InfoCode message, const char *vars);

    /**
     \brief Method to output a message behind two checks:
      \link reporting::reportManager::IsDebug IsDebug \endlink
      and VerboseLevel passed.\n
     If \link reporting::reportManager::IsDebug IsDebug \endlink is false,
      it will ignore the message returning the
      \link reporting::reportManager::cantReport cantReport \endlink,
      otherwise, VerboseLevel level will be checked.
     \param level
      Level to use for this message.\n
     If \link reporting::reportManager::Level Level \endlink is higher
      than this argument, message will be ignored,
      otherwise, message will be outputted to cout.
     */
    reportWrapper log(VerboseLevel level)
    {
        if (!IsDebug) return cantReport;
        
        if (level >= Level)
            return canReport;
        else
            return cantReport;
    }
    
    /**
     \brief Overloaded operator that allows us to debug some messages behind
      \link reporting::reportManager::IsDebug IsDebug \endlink
     \param a
      Message to be outputed if
      \link reporting::reportManager::IsDebug IsDebug \endlink
     */
    template <typename T>
    reportManager &operator<<(const T &a) {
        if (IsDebug)
            std::cout<<a;
        return *this;
    }

    /**
     \brief Overloaded operator that allows us to debug some messages behind
      \link reporting::reportManager::IsDebug IsDebug \endlink
     \param pf
      Object that overloads the '<<' operator
     */
    reportManager &operator<<(std::ostream& (*pf) (std::ostream&)) {
        if (IsDebug)
            std::cout<<pf;
        return *this;
    }
};

}
/** \brief <b> This instance is the one that should be used</b> \n
 * It's use is similar to a singleton, without the need to obtain
 *  the instance every time it's needed.\n
 * Instead, it's a global instance of the reporting::reportManager \n
 * This allows us to have the '<<' operator overloaded,
 *  as it can't be statically overloaded.
    \relates reporting::reportManager
 */
extern reporting::reportManager debug;

#endif // VERBOSEMANAGER_H

