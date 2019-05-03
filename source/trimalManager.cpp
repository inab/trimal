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


#include <trimalManager.h>

#include "FormatHandling/BaseFormatHandler.h"
#include "Statistics/similarityMatrix.h"
#include "Statistics/Consistency.h"
#include "Statistics/Similarity.h"
#include "InternalBenchmarker.h"
#include "Alignment/Alignment.h"
#include "Statistics/Manager.h"
#include "trimalManager.h"
#include "reportsystem.h"
#include "VCFHandler.h"
#include "Cleaner.h"
#include "defines.h"
#include "values.h"
#include "utils.h"

// Macros defined on this file will be defined
//  before the method that uses them.
// This is due to macros that are only used on
//  one method.

trimAlManager::trimAlManager() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("trimAlManager::trimAlManager() ");
}

trimAlManager::~trimAlManager() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("trimAlManager::~trimAlManager() ");

    delete_variables();
}

#define checkArgument(argument) { \
    auto x = argument(&argc, argv, &i); \
    if (x == Recognized) continue; \
    if (x == Errored) return x;\
    if (x == Final) return x;\
}

int trimAlManager::parseArguments(int argc, char **argv) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void trimAlManager::parseArguments(int argc, char *argv[]) ");

    // Call to verbosity argument has to be done
    //  before any other argument.
    verbosity_argument(&argc, argv);

    // Check if there is only an argument passed to trimAl,
    //  if so, show the menu and the examples
    if (argc == 1) {
        menu();
        examples();
        return Final;
    }

    // Check every provided argument
    //  with each argument method
    for (int i = 1; i < argc; i++) {
        if (appearErrors) break;

        checkArgument(help_arguments);

        // Check arguments. Makes use of the 
        //  macro checkArgument (defined before the method)
        //  to make code more understandable.
        //
        // Argument passed to checkArgument is,
        //  indeed, a method with signature:
        //  x_argument(&argc, argv, &i)

        // In-Out
        {
            // Alignment in-out
            checkArgument(in_argument)
            checkArgument(out_argument)
            checkArgument(out_format_arguments)

            // Trimming Reports
            checkArgument(html_out_argument)
            checkArgument(svg_out_argument)

            // Compare set
            checkArgument(compareset_argument)
            checkArgument(force_select_argument)

            // Back translation
            checkArgument(back_trans_argument)

            // Stats report
            checkArgument(stats_arguments)
            checkArgument(svg_stats_argument)

            // Time Tracker out
            checkArgument(timetracker_out_argument)

            // Custom Matrix argument
            checkArgument(matrix_argument)

            // Column Numbering
            checkArgument(col_numbering_argument)
        }

        // Stats
        {
            // Thresholds
            checkArgument(gap_threshold_argument)
            checkArgument(similarity_threshold_argument)
            checkArgument(consistency_threshold_argument)
            checkArgument(conservation_threshold_argument)

            // Windows
            checkArgument(window_argument)
            checkArgument(gap_window_argument)
            checkArgument(similarity_window_argument)
            checkArgument(consistency_window_argument)
        }

        // Trimming methods
        {
            // Gaps based
            checkArgument(no_gaps_argument)
            checkArgument(no_all_gaps_argument)
            // Automated
            checkArgument(gappy_out_argument)
            checkArgument(strict_argument)
            checkArgument(strict_plus_argument)
            checkArgument(automated1_argument)
            // Overlap
            checkArgument(residue_overlap_argument)
            checkArgument(sequence_overlap_argument)
            // Manual
            checkArgument(seqs_select_argument)
            checkArgument(select_cols_argument)
            // Duplicates
            checkArgument(remove_duplicates_argument)
            // Identity
            checkArgument(max_identity_argument)
            checkArgument(clusters_argument)
        }

        // Modifiers
        {
            checkArgument(keep_seqs_argument)
            checkArgument(keep_header_argument)

            // Back translation Modifiers
            checkArgument(ignore_stop_codon_argument)
            checkArgument(split_by_stop_codon_argument)

            checkArgument(block_argument)

            checkArgument(complementary_argument)

            checkArgument(terminal_only_argument)
        }

        // VCF
        {
            checkArgument(vcf_argument)
            checkArgument(ignore_filter_argument)
            checkArgument(min_quality_argument)
            checkArgument(min_coverage_argument)
        }

        // Skip the verbosity option, as it was checked before the loop
        if (!strcmp(argv[i], "--verbosity") || !strcmp(argv[i], "-v")) {
            i++;
            continue;
        }

        debug.report(ErrorCode::ParameterNotFoundOrRepeated, argv[i]);
        appearErrors = true;
        break;
    }

    if (performCompareset())
        appearErrors = true;

    // Check if we've provided multiple correct arguments
    // but no input file is provided.
    if (infile == nullptr && compareset == nullptr) {
        this->appearErrors = true;
        debug.report(ErrorCode::NoInputFile);
        return Errored;
    }
    return this->appearErrors ? Errored : Recognized;
}

#undef checkArgument

/**inline**/ bool trimAlManager::check_arguments_incompatibilities() {
    StartTiming("bool trimAlManager::check_arguments_incompatibilities()");
    // The incompatibilities are checked only once,
    // so there are arguments with no function to check it's incompatibilities although they have.
    // These are checked within other functions.
    // So if argument A is incompatible with B,
    // A may have this checked in it's incompatibilities function,
    // but B may have no function to check them.

    check_inFile_incompatibilities();
    check_select_cols_and_seqs_incompatibilities();
    check_thresholds_incompatibilities();
    check_automated_methods_incompatibilities();
    check_max_identity_incompatibilities();
    check_clusters_incompatibilities();
    check_windows_incompatibilities();
    check_stats_incompatibilities();
    check_codon_behaviour_incompatibility();
    check_combinations_among_thresholds_incompatibility();

    check_vcf_incompatibility();

    return appearErrors;
}

/**inline**/ void trimAlManager::verbosity_argument(const int *argc, char *argv[]) {
    for (int i = 1; i < *argc; i++) {
        if (!strcmp(argv[i], "--verbosity") || !strcmp(argv[i], "-v")) {
            if ((i + 1) != *argc) {
                if (!strcmp(argv[i + 1], "error") || !strcmp(argv[i + 1], "3")) {
                    debug.Level = VerboseLevel::ERROR;
                    return;
                }
                if (!strcmp(argv[i + 1], "warning") || !strcmp(argv[i + 1], "2")) {
                    debug.Level = VerboseLevel::WARNING;
                    return;
                }
                if (!strcmp(argv[i + 1], "info") || !strcmp(argv[i + 1], "1")) {
                    debug.Level = VerboseLevel::INFO;
                    return;
                }
                if (!strcmp(argv[i + 1], "none") || !strcmp(argv[i + 1], "0")) {
                    debug.Level = VerboseLevel::NONE;
                    return;
                }

                debug.report(ErrorCode::VerboseLevelNotRecognized, new std::string[2]{argv[i + 1], std::to_string(debug.Level)});
            } else
                debug.report(ErrorCode::NeedToSpecifyVerboseLevel, new std::string[2]{argv[i], std::to_string(debug.Level)});
        }
    }
}

/**inline**/ trimAlManager::argumentReport trimAlManager::help_arguments(const int *argc, char **argv, int *currentArg) {
    if (!strcmp(argv[*currentArg], "-h") || !strcmp(argv[*currentArg], "-help")) {
        menu();
        examples();
        return Final; // We don't want to continue if we show the help.
    }

    if (!strcmp(argv[*currentArg], "--version")) {
        std::cout << VERSION << std::endl;
        return Final; // We don't want to continue if we show the version.
    }

    if (!strcmp(argv[*currentArg], "-lf") || !strcmp(argv[*currentArg], "--listformats")) {
        std::cout << "Input Formats:  \t" << formatManager.getInputFormatsAvailable() << "\n\n";
        std::cout << "Output Formats: \t" << formatManager.getOutputFormatsAvailable() << "\n";
        return Final;
    }

    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::in_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-in") && ((*i) + 1 != *argc) && infile == nullptr)
    {
        // Need to check this here, as later origAlig is set and checked
        //      the same happens on force_select_argument, and thus,
        //      one alignment would be leaked.
        if (forceFile != nullptr)
        {
            debug.report(ForceSelectAndInArgumentsProvided);
            return Errored;
        }
        argumentLength = strlen(argv[++*i]);
        infile = new char[argumentLength + 1];
        strcpy(infile, argv[*i]);
        if ((origAlig = formatManager.loadAlignment(infile)) == nullptr) {
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::vcf_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-vcf") && ((*i) + 1 != *argc)) {
        vcfs = new std::vector<std::string>();
        while (((*i) + 1 != *argc)) {
            ++*i;
            if (argv[*i][0] == '-') {
                --*i;
                break;
            }
            vcfs->emplace_back(argv[*i]);
        }
        return Recognized;

    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::out_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-out")) && ((*i) + 1 != *argc) && (outfile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        outfile = new char[argumentLength + 1];
        strcpy(outfile, argv[*i]);
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::html_out_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-htmlout")) && ((*i) + 1 != *argc) && (htmlOutFile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        htmlOutFile = new char[argumentLength + 1];
        strcpy(htmlOutFile, argv[*i]);
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::timetracker_out_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-timetrackerout")) && ((*i) + 1 != *argc)) {
        (*i) += 1;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::svg_out_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-svgout")) && ((*i) + 1 != *argc) && (svgOutFile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        svgOutFile = new char[argumentLength + 1];
        strcpy(svgOutFile, argv[*i]);
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::svg_stats_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-svgstats")) && ((*i) + 1 != *argc) && (svgStatsOutFile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        svgStatsOutFile = new char[argumentLength + 1];
        strcpy(svgStatsOutFile, argv[*i]);
        return Recognized;
    }
    return NotRecognized;
}

// Macro to allow using legacy format arguments
//  with the new system.
#define LegacyFormatArgumentWrapper(arg, format) \
    if (!strcmp(argv[*i], arg)) { \
        oformats.emplace_back(format); \
        return Recognized; \
    }

/**inline**/ trimAlManager::argumentReport trimAlManager::out_format_arguments(const int *argc, char *argv[], int *i) {
    // Detect formats using the ReadWriteMachineState (RWMS) method
    //  Store the formats desired on 'oformats' list, which
    //  will be passed to the RWMS
    if (!strcmp(argv[*i], "-formats")) {
        if ((*i + 1) == *argc) {
            debug.report(ErrorCode::NoFormatsSpecified);
            appearErrors = true;
            return Errored;
        }
        while (++(*i) != *argc && argv[*i][0] != '-')
        {
            if (formatManager.getFormatFromToken(argv[*i]) == nullptr)
            {
                debug.report(ErrorCode::OutputFormatNotRecognized, argv[*i]);
                appearErrors = true;
                return Errored;
            } else oformats.emplace_back(argv[*i]);
        }
        (*i)--;
        return Recognized;
    }

    // Legacy formats handling code.
    //  Makes use of the macro LegacyFormatArgumentWrapper 
    //  for the sake of understandability.
    // This transforms the old arguments to the RWMS format,
    //  emplacing the IDS for each format on oformats list.
    LegacyFormatArgumentWrapper("-clustal",         "clustal")
    LegacyFormatArgumentWrapper("-fasta",           "fasta")
    LegacyFormatArgumentWrapper("-fasta_m10",       "fasta_m10")
    LegacyFormatArgumentWrapper("-nbrf",            "pir")
    LegacyFormatArgumentWrapper("-nexus",           "nexus")
    LegacyFormatArgumentWrapper("-mega",            "mega")
    LegacyFormatArgumentWrapper("-phylip3.2",       "phylip32")
    LegacyFormatArgumentWrapper("-phylip3.2_m10",   "phylip32_m10")
    LegacyFormatArgumentWrapper("-phylip",          "phylip40")
    LegacyFormatArgumentWrapper("-phylip_m10",      "phylip40_m10")
    LegacyFormatArgumentWrapper("-phylip_paml",     "phylip_paml")
    LegacyFormatArgumentWrapper("-phylip_paml_m10", "phylip_paml_m10")

    return NotRecognized;
}

#undef LegacyFormatArgumentWrapper

/**inline**/ trimAlManager::argumentReport trimAlManager::matrix_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-matrix") && ((*i) + 1 != *argc) && (matrixFile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        matrixFile = new char[argumentLength + 1];
        strcpy(matrixFile, argv[*i]);
        return Recognized;
    } else if (!strcmp(argv[*i], "--alternative_matrix") && ((*i) + 1 != *argc) && (alternative_matrix == -1)) {
        (*i)++;
        if (!strcmp(argv[*i], "degenerated_nt_identity"))
        {
            alternative_matrix = 1;
            return Recognized;
        }
        else {
            debug.report(ErrorCode::AlternativeMatrixNotRecognized, argv[*i]);
            appearErrors = true;
        }
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::compareset_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-compareset") && ((*i) + 1 != *argc) && (compareset == nullptr)) {
        // Check if file can be opened
        compare.open(argv[++*i], std::ifstream::in);
        if (!compare) {
            debug.report(ErrorCode::ReferenceFileNotLoaded, argv[*i]);
            appearErrors = true;
        }
        compare.close();
        compareset = argv[*i];
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::force_select_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-forceselect") && ((*i) + 1 != *argc) && forceFile == nullptr)
    {
        // Need to check this here, as later origAlig is set and checked
        //      the same happens on in_argument, and thus,
        //      one alignment would be leaked.
        if (infile != nullptr)
        {
            debug.report(ForceSelectAndInArgumentsProvided);
            return Errored;
        }
        argumentLength = strlen(argv[++*i]);
        forceFile = new char[argumentLength + 1];
        strcpy(forceFile, argv[*i]);

        if ((origAlig = formatManager.loadAlignment(forceFile)) == nullptr) {
            debug.report(ErrorCode::AlignmentNotLoaded, forceFile);
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::back_trans_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-backtrans") && ((*i) + 1 != *argc) && (backtransFile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        backtransFile = new char[argumentLength + 1];
        strcpy(backtransFile, argv[*i]);

        if ((backtranslationAlig = formatManager.loadAlignment(backtransFile)) == nullptr) {
            debug.report(ErrorCode::AlignmentNotLoaded, backtransFile);
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::gap_threshold_argument(const int *argc, char *argv[], int *i) {

    if ((!strcmp(argv[*i], "-gapthreshold") || !strcmp(argv[*i], "-gt"))
        && ((*i) + 1 != *argc) && (gapThreshold == -1)) {

        if (gapAbsoluteThreshold != -1) {
            debug.report(ErrorCode::AbsoluteAndRelativeGapThreshold);
            appearErrors = true;
        }

        if (utils::isNumber(argv[++*i])) {
            gapThreshold = 1.0F - (float)atof(argv[*i]);
            if ((gapThreshold < 0) || (gapThreshold > 1)) {
                debug.report(ErrorCode::GapThresholdOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::GapThresholdNotRecognized);
            appearErrors = true;
        }
        return Recognized;
    }

    if ((!strcmp(argv[*i], "-gapabsolutethreshold") || !strcmp(argv[*i], "-gat"))
        && ((*i) + 1 != *argc) && (gapAbsoluteThreshold == -1)) {

        if (gapThreshold != -1) {
            debug.report(ErrorCode::AbsoluteAndRelativeGapThreshold);
            appearErrors = true;
        }

        if (utils::isNumber(argv[++*i])) {
            gapAbsoluteThreshold = atoi(argv[*i]);
            if (gapAbsoluteThreshold < 0) {
                debug.report(ErrorCode::AbsoluteGapThresholdLessThanZero);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::AbsoluteGapThresholdNotRecognized);
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;

}

/**inline**/ trimAlManager::argumentReport trimAlManager::similarity_threshold_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-simthreshold") || !strcmp(argv[*i], "-st")) && ((*i) + 1 != *argc) && (similarityThreshold == -1)) {
        if (utils::isNumber(argv[++*i])) {
            similarityThreshold = atof(argv[*i]);
            if ((similarityThreshold < 0) || (similarityThreshold > 1)) {
                debug.report(ErrorCode::SimilarityThresholdOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::SimilarityThresholdNotRecognized);
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::consistency_threshold_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-conthreshold") || !strcmp(argv[*i], "-ct")) && ((*i) + 1 != *argc) && (consistencyThreshold == -1)) {
        if (utils::isNumber(argv[++*i])) {
            consistencyThreshold = atof(argv[*i]);
            if ((consistencyThreshold < 0) || (consistencyThreshold > 1)) {
                debug.report(ErrorCode::ConsistencyThresholdOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::ConsistencyThresholdNotRecognized);
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::conservation_threshold_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-cons")) && ((*i) + 1 != *argc) && (conservationThreshold == -1)) {
        if (utils::isNumber(argv[++*i])) {
            conservationThreshold = atof(argv[*i]);
            if ((conservationThreshold < 0) || (conservationThreshold > 100)) {
                debug.report(ErrorCode::ConservationThresholdOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::ConservationThresholdNotRecognized);
            appearErrors = true;
        }

        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::select_cols_argument(const int *argc, char *argv[], int *i) {
    StartTiming("trimAlManager::select_cols_argument");
    if ((!strcmp(argv[*i], "-selectcols")) &&
        !selectCols &&
        (*i + 3) < *argc &&
        !strcmp(argv[++(*i)], "{") &&
        !strcmp(argv[(*i) + 2], "}")) {
        if ((delColumns = utils::readNumbers(argv[++(*i)])) == nullptr) {
            debug.report(ErrorCode::SelectColsNotRecognized);
            appearErrors = true;
        } else selectCols = true;
        (*i)++;

        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::no_gaps_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-nogaps") && (!nogaps)) {
        nogaps = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::remove_duplicates_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-noduplicateseqs") && (!removeDuplicates)) {
        removeDuplicates = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::no_all_gaps_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-noallgaps") && (!noallgaps)) {
        noallgaps = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::keep_seqs_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-keepseqs") && (!keepSeqs)) {
        keepSeqs = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::keep_header_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-keepheader") && (!formatManager.keepHeader)) {
        formatManager.keepHeader = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::gappy_out_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-gappyout") && (!gappyout)) {
        gappyout = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::strict_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-strict") && (!strict)) {
        strict = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::strict_plus_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-strictplus")) && (!strictplus)) {
        strictplus = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::automated1_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-automated1")) && (!automated1)) {
        automated1 = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::residue_overlap_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-resoverlap")) && ((*i) + 1 != *argc) && (residuesOverlap == -1)) {
        if (utils::isNumber(argv[++*i])) {
            residuesOverlap = atof(argv[*i]);
            if ((residuesOverlap < 0) || (residuesOverlap > 1)) {
                debug.report(ErrorCode::ResidueOverlapOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::ResidueOverlapNotRecognized);
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::sequence_overlap_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-seqoverlap")) && ((*i) + 1 != *argc) && (sequenceOverlap == -1)) {
        if (utils::isNumber(argv[++*i])) {
            sequenceOverlap = atof(argv[*i]);
            if ((sequenceOverlap < 0) || (sequenceOverlap > 100)) {
                debug.report(ErrorCode::SequencesOverlapOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::SequencesOverlapNotRecognized);
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::seqs_select_argument(const int *argc, char *argv[], int *i) {
    StartTiming("trimAlManager::select_seqs_argument");
    if ((!strcmp(argv[*i], "-selectseqs")) &&
         !selectSeqs &&
         ((*i + 3) < *argc) &&
         (!strcmp(argv[++*i], "{")) &&
         (!strcmp(argv[*i + 2], "}"))) {
        if ((delSequences = utils::readNumbers(argv[++*i])) == nullptr) {
            debug.report(ErrorCode::SelectSeqsNotRecognized);
            appearErrors = true;
        } else selectSeqs = true;
        (*i)++;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::max_identity_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-maxidentity")) && ((*i) + 1 != *argc) && (maxIdentity == -1)) {

        if (utils::isNumber(argv[++*i])) {
            maxIdentity = atof(argv[*i]);
            if ((maxIdentity < 0) || (maxIdentity > 1)) {
                debug.report(ErrorCode::MaxIdentityOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::MaxIdentityNotRecognized);
            appearErrors = true;
        }

        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::clusters_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-clusters")) && ((*i) + 1 != *argc) && (clusters == -1)) {
        if (utils::isNumber(argv[++*i])) {
            clusters = atoi(argv[*i]);
            if (clusters < 1) {
                debug.report(ErrorCode::ClustersValueOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::ClustersValueNotRecognized);
            appearErrors = true;
        }

        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::terminal_only_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-terminalonly")) && (!terminalOnly)) {
        terminalOnly = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::window_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-w") && ((*i) + 1 != *argc) && (windowSize == -1)) {
        if (utils::isNumber(argv[*i + 1])) {
            windowSize = atoi(argv[++*i]);
            if (windowSize <= 0) {
                debug.report(ErrorCode::WindowValueOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::WindowValueNotRecognized);
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::gap_window_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-gw") && ((*i) + 1 != *argc) && (gapWindow == -1)) {
        if (utils::isNumber(argv[*i + 1])) {
            gapWindow = atoi(argv[++*i]);
            if (gapWindow <= 0) {
                debug.report(ErrorCode::GapWindowValueOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::GapWindowValueNotRecognized);
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::similarity_window_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-sw") && ((*i) + 1 != *argc) && (similarityWindow == -1)) {
        if (utils::isNumber(argv[*i + 1])) {
            similarityWindow = atoi(argv[++*i]);
            if (similarityWindow <= 0) {
                debug.report(ErrorCode::SimilarityWindowValueOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::SimilarityWindowValueNotRecognized);
            appearErrors = true;
        }

        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::consistency_window_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-cw") && ((*i) + 1 != *argc) && (consistencyWindow == -1)) {
        if (utils::isNumber(argv[*i + 1])) {
            consistencyWindow = atoi(argv[++*i]);
            if (consistencyWindow <= 0) {
                debug.report(ErrorCode::ConsistencyWindowValueOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::ConsistencyWindowValueNotRecognized);
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::block_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-block") && ((*i) + 1 != *argc) && (blockSize == -1)) {
        if (utils::isNumber(argv[*i + 1])) {
            blockSize = atoi(argv[++*i]);
            if (blockSize <= 0) {
                debug.report(ErrorCode::BlockSizeOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::BlockSizeNotRecognized);
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;
}

// Macros to make the code on stats_arguments easier to understand.
//  The macro needs the stat variable name to be the same as the argument used. 

// This macro also tightens the relation between a stats argument and 
//  the variable name where that information is stored.

// This is an intentional coupling that would help when the number of variables
//  arises, as their names may be not fully informative.

#define stat_check(stat) \
    if (!strcmp(argv[*i], "-" #stat )) { \
        if (!stat) { \
            stat = true; \
            stats--; \
            return Recognized; \
        } else return NotRecognized; \
    }

/**inline**/ trimAlManager::argumentReport trimAlManager::stats_arguments(const int *argc, char *argv[], int *i) {
    stat_check(sgc)
    stat_check(sgt)
    stat_check(ssc)
    stat_check(sst)

    stat_check(sident)
    stat_check(soverlap)

    stat_check(sfc)
    stat_check(sft)

    // As macro stat_check is able to make the method to return values 
    //  under certain circumstances, althought the method (stats_arguments) 
    //  seem to always return false, this is incorrect and the method can
    //  return True or False, depeding if a stat argument has been recognized.
    return NotRecognized;
}

#undef stat_check

/**inline**/ trimAlManager::argumentReport trimAlManager::complementary_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-complementary")) && !getComplementary) {
        getComplementary = true;
        return Recognized;
    }
    if ((!strcmp(argv[*i], "-complementaryseq")) && !getComplementarySeq) {
        getComplementarySeq = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::col_numbering_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-colnumbering")) && !columnNumbering) {
        columnNumbering = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::split_by_stop_codon_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-splitbystopcodon")) && !splitByStopCodon) {
        splitByStopCodon = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::ignore_stop_codon_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-ignorestopcodon")) && !ignoreStopCodon) {
        ignoreStopCodon = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::ignore_filter_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-ignorefilter")) && !ignoreFilter) {
        ignoreFilter = true;
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::min_quality_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-minquality") && ((*i) + 1 != *argc) && (minQuality == -1)) {
        if (utils::isNumber(argv[*i + 1])) {
            minQuality = atof(argv[++*i]);
            if (minQuality < 0) {
                debug.report(ErrorCode::MinQualityLesserThan0);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::MinQualityNotRecognized);
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;
}

/**inline**/ trimAlManager::argumentReport trimAlManager::min_coverage_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-mincoverage") && ((*i) + 1 != *argc) && (minCoverage == -1)) {
        if (utils::isNumber(argv[*i + 1])) {
            minCoverage = atof(argv[++*i]);
            if (minCoverage < 0) {
                debug.report(ErrorCode::MinCoverageLesserThan0);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::MinCoverageNotRecognized);
            appearErrors = true;
        }
        return Recognized;
    }
    return NotRecognized;
}

bool trimAlManager::processArguments(char *argv[]) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool trimAlManager::processArguments(char *argv[]) ");
    if (!appearErrors) {

        // This allows us to have a count of all the automatic methods 
        //  that have been requested.
        // We can use this information to prevent performing more than one method,
        //  and using as a bool, to check if any automatic method has been used.
        automatedMethodCount =
                nogaps      + noallgaps +
                gappyout    + strict    +
                strictplus  + automated1 + removeDuplicates;

        check_arguments_incompatibilities();
        check_arguments_needs(argv);
    }
    return appearErrors;
}

/**inline**/ bool trimAlManager::check_inFile_incompatibilities() {
    if (infile != nullptr) {
        if ((sfc) || (sft) || (consistencyThreshold != -1)) {
            debug.report(ErrorCode::InFileComparisonStatistics);
            appearErrors = true;
            i++;
        }
        if (compareset != nullptr) {
            debug.report(ErrorCode::IncompatibleArguments, new std::string[2]{"-in", "-compareset"});
            appearErrors = true;
        }
        if (forceFile != nullptr) {
            debug.report(ErrorCode::IncompatibleArguments, new std::string[2]{"-in", "-forceselect"});
            appearErrors = true;
        }
    }
    return appearErrors;
}

/**inline**/ bool trimAlManager::check_select_cols_and_seqs_incompatibilities() {
    if (selectCols || selectSeqs) {
        if ((clusters != -1) || (maxIdentity != -1)) {
            debug.report(ErrorCode::OnlyOneSequencesSelectionMethodAllowed);
            appearErrors = true;
        }

        if (selectCols)
        {
            if (blockSize != -1) {
                debug.report(ErrorCode::IncompatibleArguments,
                             new std::string[2]{"-selectcols", "-block"});
                appearErrors = true;
            }

            for (int i = 1; i <= delColumns[0]; i++)
            {
                if (delColumns[i] >= origAlig->numberOfResidues) {
                    debug.report(ErrorCode::SelectOnlyAccepts,
                                 new std::string[2]{"-selectcols", "columns"});
                    appearErrors = true;
                    break;
                }
            }
        }

        if (selectSeqs)

            for (int i = 1; i <= delSequences[0]; i++)
            {
                if (delSequences[i] >= origAlig->numberOfSequences) {
                    debug.report(ErrorCode::SelectOnlyAccepts,
                                 new std::string[2]{"-selectseqs", "sequences"});
                    appearErrors = true;
                    break;
                }
            }
    }
    return appearErrors;
}

/**inline**/ bool trimAlManager::check_thresholds_incompatibilities() {
    if ((gapThreshold != -1)            || 
        (similarityThreshold != -1)     || 
        (consistencyThreshold != -1)    || 
        (conservationThreshold != -1)
        ) 
    {


        if (automatedMethodCount) {

            const char *autom;
            if (strict)         autom = "-strict";
            if (strictplus)     autom = "-strictplus";
            if (gappyout)       autom = "-gappyout";
            if (nogaps)         autom = "-nogaps";
            if (noallgaps)      autom = "-noallgaps";
            if (automated1)     autom = "-automated1";


            if (gapThreshold != -1)
                debug.report(ErrorCode::CombinationAmongTrimmingMethods,  new std::string[2]{autom,"-gt"});
            if (similarityThreshold != -1)
                debug.report(ErrorCode::CombinationAmongTrimmingMethods,  new std::string[2]{autom,"-st"});
            if (consistencyThreshold != -1)
                debug.report(ErrorCode::CombinationAmongTrimmingMethods,  new std::string[2]{autom,"-ct"});
            if (conservationThreshold != -1)
                debug.report(ErrorCode::CombinationAmongTrimmingMethods,  new std::string[2]{autom,"-cons"});

            appearErrors = true;
        }

        if (conservationThreshold != -1) {
            if (blockSize != -1) {
                debug.report(ErrorCode::IncompatibleArguments,
                             new std::string[2]{"-conthreshold", "-block"});
                appearErrors = true;
            }
        }
    }
    return appearErrors;
}

/**inline**/ bool trimAlManager::check_automated_methods_incompatibilities() {
    if (automatedMethodCount) {

        const char *autom;
        if (strict)         autom = "-strict";
        if (strictplus)     autom = "-strictplus";
        if (gappyout)       autom = "-gappyout";
        if (nogaps)         autom = "-nogaps";
        if (noallgaps)      autom = "-noallgaps";
        if (automated1)     autom = "-automated1";


        if ((windowSize != -1)) {
            debug.report(ErrorCode::CombinationAmongTrimmingMethods, new std::string[2]{autom, "-w"});
            appearErrors = true;
        }
        if ((gapWindow != -1)) {
            debug.report(ErrorCode::CombinationAmongTrimmingMethods, new std::string[2]{autom,"-gw"});
            appearErrors = true;
        }
        if ((similarityWindow != -1)) {
            debug.report(ErrorCode::CombinationAmongTrimmingMethods, new std::string[2]{autom,"-sw"});
            appearErrors = true;
        }

        if (automatedMethodCount > 1) {
            debug.report(ErrorCode::MoreThanOneAutomatedMethod);
            appearErrors = true;
        }
    }
    return appearErrors;
}

/**inline**/ bool trimAlManager::check_max_identity_incompatibilities() {
    if (maxIdentity != -1) {
        if ((windowSize != -1) || (gapWindow != -1) || (similarityWindow != -1) || (consistencyWindow != -1)) {
            debug.report(ErrorCode::WindowAndArgumentIncompatibilities, new std::string[1]{"-maxIdentity"});
            appearErrors = true;
        }
        if (clusters != -1) {
            debug.report(ErrorCode::OnlyOneSequencesSelectionMethodAllowed);
            appearErrors = true;
        }
    }
    return appearErrors;
}

/**inline**/ bool trimAlManager::check_clusters_incompatibilities() {
    if (clusters != -1) {
        if (clusters > origAlig->getNumSpecies()) {
            debug.report(ErrorCode::MoreClustersThanSequences);
            appearErrors = true;
        }
    }
    return appearErrors;
}

/**inline**/ bool trimAlManager::check_windows_incompatibilities() {
    if (windowSize != -1) {
        if (consistencyWindow != -1 || gapWindow != -1 || similarityWindow != -1) {
            debug.report(ErrorCode::GeneralAndSpecificWindows);
            appearErrors = true;
        }
    }
    return appearErrors;
}

/**inline**/ bool trimAlManager::check_stats_incompatibilities() {
    if (stats < 0) {
        if (columnNumbering) {
            debug.report(ErrorCode::StatisticsArgumentIncompatibilities, new std::string[1]{"-colnumbering"});
            appearErrors = true;
        }
    }
    return appearErrors;
}

/**inline**/ bool trimAlManager::check_arguments_needs(char *argv[]) {
    StartTiming("bool trimAlManager::check_arguments_needs(char *argv[])");

    check_absolute_gap_theshold();

    check_force_selection();
    check_input_file_with_coding_sequences_argument();
    check_file_aligned();
    check_similarity_matrix();
    check_outputs_coincidence();
    check_col_numbering();
    check_residue_and_sequence_overlap();
    check_output_relevance();
    check_output_file_with_statistics();
    check_automated_manual_incompatibilities();
    check_block_size();
    check_backtranslations();
    check_coding_sequences_type();
    check_and_prepare_coding_sequence();
    check_backtranslation_infile_names_correspondence();
    check_compareset_window_argument();
    check_output_format();

    check_thresholds_dependencies();
    return appearErrors;
}

/**inline**/ bool trimAlManager::check_codon_behaviour_incompatibility() {
    if ((!appearErrors) && (ignoreStopCodon) && (splitByStopCodon)) {
        debug.report(ErrorCode::IncompatibleArguments, new std::string[2]{"-ignorestopcodon", "-splitbystopcodon"});
        appearErrors = true;
        return true;
    }
    return false;
}

/**inline**/ bool trimAlManager::check_combinations_among_thresholds_incompatibility()
{
    if ((consistencyThreshold != -1) 
        && (conservationThreshold != -1) 
        && (!appearErrors)) 
    {
        if ((gapThreshold != -1) || (similarityThreshold != -1)) {
            debug.report(ErrorCode::CombinationAmongThresholdsMethods);
            appearErrors = true;
            return true;
        }
    }
    return false;
}

/**inline**/ bool trimAlManager::check_vcf_incompatibility()
{
    bool returnValue = false;
    if (vcfs == nullptr)
    {
        if (minQuality != -1)
        {
            appearErrors = true;
            debug.report(ErrorCode::OnlyValidWithVCF, "MinQuality");
        }

        if (minCoverage != -1)
        {
            appearErrors = true;
            debug.report(ErrorCode::OnlyValidWithVCF, "MinCoverage");
        }

        if (ignoreFilter)
        {
            appearErrors = true;
            debug.report(ErrorCode::OnlyValidWithVCF, "IgnoreFilter");
        }
    }
    return returnValue;
}


/**inline**/ bool trimAlManager::check_automated_manual_incompatibilities() {
    if ((getComplementary) && (!appearErrors))
        if (!automatedMethodCount && // Are we not using an automated method?
            (gapThreshold == -1) && (consistencyThreshold == -1) && (similarityThreshold == -1) && // Neither a threshold method.
            (!selectCols) && (!selectSeqs) && (residuesOverlap == -1) && (sequenceOverlap == -1) && // Neither a sequence and residues semimanual selection methods
            (maxIdentity == -1) && (clusters == -1)) // Or complex selection of sequences.
        {
            debug.report(ErrorCode::TrimmingMethodNeeded, new std::string[1]{"-complementary"});
            appearErrors = true;
            return true;
        }

    if ((terminalOnly) && (!appearErrors))
        if (!automatedMethodCount && // Are we not using an automated method?
            (gapThreshold == -1) && (conservationThreshold == -1) && (similarityThreshold == -1) && // Neither a threshold method.
            (!selectCols) && (!selectSeqs) && (residuesOverlap == -1) && (sequenceOverlap == -1) && // Neither a sequence and residues semimanual selection methods
            (maxIdentity == -1) && (clusters == -1)) // Or complex selection of sequences.
        {
            debug.report(ErrorCode::TrimmingMethodNeeded, new std::string[1]{"-terminalonly"});
            appearErrors = true;
            return true;
        }
    return false;
}

inline bool trimAlManager::check_absolute_gap_theshold() {
    if (!appearErrors) {
        if (gapAbsoluteThreshold != -1)
        {
            gapThreshold = ((float)gapAbsoluteThreshold) / origAlig->originalNumberOfSequences;
            if (gapThreshold >= 1.0F) {
                appearErrors = true;
                debug.report(
                        ErrorCode::AbsoluteGapThresholdBiggerThanNumberOfSequences,
                        new std::string[2] {
                    std::to_string(origAlig->originalNumberOfSequences),
                    std::to_string(gapAbsoluteThreshold)
                    });
                gapThreshold = -1;
                return true;
            }
        }
    }
    return false;
}

/**inline**/ bool trimAlManager::check_force_selection() {
    if (!appearErrors) {
        if ((compareset == nullptr) && (forceFile != nullptr)) {
            debug.report(ErrorCode::ForceFileWithoutCompareDataset);
            appearErrors = true;
            return true;
        }
    }
    return false;
}

/**inline**/ bool trimAlManager::check_input_file_with_coding_sequences_argument() {
    if ((!appearErrors) && (infile == nullptr) && (compareset == nullptr) && (forceFile == nullptr) && (backtransFile != nullptr)) {
        debug.report(ErrorCode::BacktranslationWithoutMainAlignment);
        appearErrors = true;
        return true;
    }
    return false;
}

/**inline**/ bool trimAlManager::check_file_aligned() {
    if ((!appearErrors) && (infile != nullptr)) {

        if (// Are we requesting an automated method ? or...
                (automatedMethodCount || 
            // Are we requesting any manual threshold method ? or...
                (gapThreshold != -1) || (consistencyThreshold != -1) || (similarityThreshold != -1) || 
            // Are we selecting columns or sequences ? or...
                (selectCols) || (selectSeqs) || 
            // Are we using max overlap between sequences ? or...
                (residuesOverlap != -1) || (sequenceOverlap != -1) || 
            // Are we asking for any stats?
                (stats < 0)) 
            &&
            // Then we need the alignment to be aligned; 
            //  If not, we should report the error
            (!origAlig->isFileAligned())
            &&
            // We are not doing VCF)
            (vcfs == nullptr))
        {
            debug.report(ErrorCode::NotAligned, new std::string[1]{infile});
            appearErrors = true;
            return true;
        }
    }
    return false;
}

/**inline**/ bool trimAlManager::check_similarity_matrix() {
    if ((matrixFile != nullptr || alternative_matrix != -1) && (!appearErrors)) {
        if ((!strict) && (!strictplus) && (!automated1) && (similarityThreshold == -1) && (!ssc) && (!sst)) {
            debug.report(ErrorCode::MatrixGivenWithNoMethodToUseIt);
            appearErrors = true;
            return true;
        }

        // TODO this is an incompatibility.
        if ((gapWindow != -1) || ((compareset == nullptr) && (consistencyWindow != -1))) {
            debug.report(ErrorCode::SimilarityMatrixNotCompatibleWindow);
            appearErrors = true;
            return true;
        }
    }
    return false;
}

/**inline**/ bool trimAlManager::check_outputs_coincidence() {
    // Check dynamically that every output path is unique
    // We use two arrays: 
    //      outFiles, which contains the variables
    //      outFilesNames, which is used to identify each method

    std::array<char *, 4> outFiles
    {{
        htmlOutFile,
        outfile,
        svgOutFile,
        svgStatsOutFile
    }};

    std::array<std::string, 4> outFilesNames
    {{
        "html report (-htmlout)",
        "output alignment (-out)",
        "svg report (-svgout)",
        "svg stats (-svgstats)"
    }};


    for (int i = 0, x = 0; i < outFiles.size(); i++) {
        if (outFiles.at(i) != nullptr)
            for (x = i + 1; x < outFiles.size(); x++) {
                if (outFiles.at(x) != nullptr)
                    if (!strcmp(outFiles.at(i), outFiles.at(x))) {
                        debug.report(ErrorCode::SameNameOutput, new std::string[2]{outFilesNames.at(i), outFilesNames.at(x)});
                        appearErrors = true;
                    }
            }
    }

    return false;
}

/**inline**/ bool trimAlManager::check_col_numbering() {

    // As colnumbering doesn't make sense if we don't trim the alignment...
    if ((columnNumbering) && (!appearErrors)) {
        if (
            // Are we not using any automated method? and...
                (!automatedMethodCount) && 
            // Are we not using any threshold? and...
                (gapThreshold == -1) && 
                (conservationThreshold == -1) && 
                (similarityThreshold == -1) && 
                (consistencyThreshold == -1) && 
            // Neither selecting any column or sequence?
                (!selectCols) && (!selectSeqs))
        {
            // If all that happens, we are not trimming the alignment and should report it.
            debug.report(ErrorCode::TrimmingMethodNeeded, new std::string[1]{"-colnumbering"});
            appearErrors = true;
            return true;
        }
    }
    return false;
}

/**inline**/ bool trimAlManager::check_residue_and_sequence_overlap() {
    if (!appearErrors) {
        if ((residuesOverlap != -1) && (sequenceOverlap == -1)) {
            debug.report(ErrorCode::SequenceAndResiduesOverlapMutuallyNeeded, new std::string[1]{"residues overlap"});
            appearErrors = true;
            return true;
        } else if ((residuesOverlap == -1) && (sequenceOverlap != -1)) {
            debug.report(ErrorCode::SequenceAndResiduesOverlapMutuallyNeeded, new std::string[1]{"sequences overlap"});
            appearErrors = true;
            return true;
        }
    }
    return false;
}

/**inline**/ bool trimAlManager::check_output_relevance() {
    if (((htmlOutFile != nullptr) ||
         (svgOutFile != nullptr) ||
         (svgStatsOutFile != nullptr)) 
         && (!appearErrors)) {
        if (
            // Are we not using any automated method? and...
                !automatedMethodCount && 
            // Neither using thresholds? and...
                (gapThreshold == -1) && (conservationThreshold == -1) && (similarityThreshold == -1) && (consistencyThreshold == -1) && 
            // Neither selecting columns or sequences? and...
                (!selectCols) && (!selectSeqs) && (residuesOverlap == -1) && (sequenceOverlap == -1) && 
            // Neither using other selecting methods? and...
                (maxIdentity == -1) && (clusters == -1))
        {
            // Then, we are not trimming the alignment and should report it.
            if (htmlOutFile != nullptr)
                debug.report(ErrorCode::TrimmingMethodNeeded, new std::string[1]{"-htmlout"});
            if (svgOutFile != nullptr)
                debug.report(ErrorCode::TrimmingMethodNeeded, new std::string[1]{"-svgout"});
            if (svgStatsOutFile != nullptr)
                debug.report(ErrorCode::TrimmingMethodNeeded, new std::string[1]{"-svgstats"});
            appearErrors = true;
            return true;
        }
    }
    return false;
}

/**inline**/ bool trimAlManager::check_output_file_with_statistics() {
    // TODO: Create a method to check if a trimming method has been requested.
    // TODO:    This will help to maintain coherence when scaling up 
    // TODO:    With more trimming methods. 

    if ((stats < 0) && (!appearErrors)) {
        stats--; 

        if (
            // If we are using an automated method
                ((automatedMethodCount) ||
            // Or a manual threshold
                (gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1)) 
            // We need the outFile to be specified.
            && (outfile == nullptr)) 
        {
            // If so, report the error to the user.
            debug.report(ErrorCode::OutFileNeededWhenPrintingStatistics);
            appearErrors = true;
            return true;
        }
    }
    return false;
}

/**inline**/ bool trimAlManager::check_block_size() {
    if ((!appearErrors) && (blockSize > origAlig->getNumAminos() / 4)) {
        debug.report(ErrorCode::BlocksizeTooBig, new std::string[1]{std::to_string(origAlig->getNumAminos() / 4)});
        appearErrors = true;
        return true;
    }
    return false;
}

/**inline**/ bool trimAlManager::check_backtranslations() {
    if (!appearErrors) {
        if (backtransFile == nullptr) {

            // Crossed dependency between Backtranslation and SplitByCodon

            if (splitByStopCodon) {
                debug.report(ErrorCode::ParemeterOnlyOnBacktranslation, new std::string[1]{"-splitbystopcodon"});
                appearErrors = true;
                return true;
            }

            // Crossed dependency between Backtranslation and IgnoreStopCodon
            if (ignoreStopCodon) {
                debug.report(ErrorCode::ParemeterOnlyOnBacktranslation, new std::string[1]{"-ignorestopcodon"});
                appearErrors = true;
                return true;
            }
        } 
        // If doing backtranslation, we need the input file to be aligned.
        else if (!origAlig->isFileAligned()) {
            debug.report(ErrorCode::ProteinAlignmentMustBeAligned);
            appearErrors = true;
            return true;
        }
    }
    return false;
}

/**inline**/ bool trimAlManager::check_coding_sequences_type() {
    if (
        // Is there a backtranslation file?
            (!appearErrors) && (backtransFile != nullptr) && 
        // If so, is it from DNA? If not...
            ( (!backtranslationAlig->getAlignmentType()) & SequenceTypes::DNA)) 
    {
        // When doing backtranslation, 
        //  a DNA alignment (without gaps) is needed.
        // Report the error to the user.
        debug.report(ErrorCode::BacktransAlignIsDNA);
        appearErrors = true;
        return true;
    }
    return false;
}

/**inline**/ bool trimAlManager::check_and_prepare_coding_sequence() {
    if ((!appearErrors) && 
        // If we are going to do a back translation
            (backtransFile != nullptr) &&
        // Perform it and save if it has errored.
            (!backtranslationAlig->prepareCodingSequence(
                splitByStopCodon,
                ignoreStopCodon,
                origAlig)
            )
        ) 
    {
        // Error reporting is made by prepareCodingSequence function.
        appearErrors = true;
        return true;
    }
    return false;
}

/**inline**/ bool trimAlManager::check_backtranslation_infile_names_correspondence() {
    //NOTE Maybe we don't need to copy the names and lengths to two new arrays 
    // as we could pass the original names and lengths to the check checkCorrespondence function, 
    // which doesn't modify the pointers passed to them

    if ((!appearErrors) && (backtransFile != nullptr)) {
        sequencesNames = new std::string[backtranslationAlig->getNumSpecies()];
        sequencesLengths = new int[backtranslationAlig->getNumSpecies()];
        backtranslationAlig->getSequences(sequencesNames, sequencesLengths);

        if (!origAlig->checkCorrespondence(
                sequencesNames,
                sequencesLengths,
                backtranslationAlig->getNumSpecies(),
                3
        )) {
            appearErrors = true;
        }
    }
    delete [] sequencesNames;
    delete [] sequencesLengths;

    return appearErrors;
}

/**inline**/ void trimAlManager::check_compareset_window_argument() {
    if ((!appearErrors) && (windowSize != -1) && (compareset != nullptr))
        debug.report(InfoCode::WindowSizeCompareset);
}

/**inline**/ void trimAlManager::check_output_format() {
    if (oformats.empty() && infile) {
        oformats.emplace_back(formatManager.getFileFormatName(infile));
    }
}

/**inline**/ void trimAlManager::check_thresholds_dependencies() {
    bool generalMakesSense = false;
    for(auto& thresholdDependency : std::vector<std::tuple<const char *, float&, int&>>{
            {"gaps",        gapThreshold,         gapWindow},
            {"similarity",  similarityThreshold,  similarityWindow},
            {"consistency", consistencyThreshold, consistencyWindow},
    })
    {
        if (std::get<2>(thresholdDependency) != -1) // Window is defined
        {
            if (std::get<1>(thresholdDependency) == -1.0F) // Threshold is not defined
            {
                this->appearErrors = true;
                std::cout << std::get<0>(thresholdDependency) << " window provided, but stat not requested.";
            }
        } else if (windowSize != -1) { // Specific window is not defined, but general is
            if (std::get<1>(thresholdDependency) != -1.0F) // Threshold is defined
            {
                generalMakesSense = true;
            }
        }
    }
    if (windowSize != -1) // General window size specified
    {
        if (!generalMakesSense) // And there isn't a stat with no specific window
        {
            this->appearErrors = true;
            std::cout << "General window provided, but no stat requested without specific window";
        }
    }
}

int trimAlManager::perform() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("int trimAlManager::perform() ");

    // We don't perform any manipulation in previous steps an error has
    // been detected
    if (appearErrors) 
        return 1;

    // If similarity threshold is -1, it hasn't been specified, and we
    // use 0 as default value, which makes no effect
    if (conservationThreshold == -1)
        conservationThreshold = 0;

    // In case we are not doing VCF transformations, the perform
    //      should be applied normally and once.
    if (vcfs == nullptr)
        return innerPerform();

    // In case we are doing VCF transformations, the perform should
    // be applied in all the collection.
    else
    {
        int returnValue = 0;
        auto XX = formatManager.splitAlignmentKeeping(*origAlig);
        char replacement = '-';
        ngs::readVCF(
                /* Dataset          */ XX,
                /* VCF Collection   */ *vcfs,
                /* min Quality      */ minQuality,
                /* min Coverage     */ minCoverage,
                /* ignore Filters   */ ignoreFilter,
                /* replacement char */ &replacement
        );

        // If no output filename has been provided,
        //      or the filename does not contain the [contig] tag,
        //      the sequence names should contain the name of the contig
        //      to prevent repeated sequences names.
        if (outfile == nullptr || std::string(outfile).find("[contig]") == std::string::npos)
        {
            // On each alignment, add the prefix to each non-reference sequence -> [1:]
            for (Alignment* &i : XX)
            {
                for (int x = 1; x < i->originalNumberOfSequences; x++)
                    i->seqsName[x] = i->seqsName[0] + "." + i->seqsName[x];
            }

            // We check if outfile is not null because checking
            //      if it contains the contig tag
            //      it needs to be converted to string
            //      and it needs to check if outfile is not null previously
            // As this two conditions have previously checked,
            //      if outfile is not null, it doesn't contain the tag either.
            if (outfile != nullptr)
            {
                // Set the mode of opening files of the manager to append.
                formatManager.openmode = std::ofstream::out | std::ofstream::app;

                // Store the output filename
                std::string fname = formatManager.replaceINtag(*origAlig, outfile);

                // Iterate over all formats requested
                for (const std::string & token : oformats) {

                    // Get the handler for each format handler
                    auto * handler = formatManager.getFormatFromToken(token);
                    if (handler != nullptr)
                    {
                        // Open the destination file to empty it
                        std::string newName =
                                utils::ReplaceString(
                                        fname, "[format]", handler->name);
                        utils::ReplaceStringInPlace(
                                newName, "[extension]", handler->extension);
                        std::ofstream tmp(newName);
                    }
                }
            }
        }

        // To allow users to use [in] and [contig] tags, we should recover the original filename.
        std::string originalName = origAlig->filename;

        // We can get rid of the original alignment.
        delete origAlig;
        origAlig = nullptr;

        // Save the original filename pattern containing the tags
        std::string originalOutFile = outfile == nullptr ? "" : outfile;

        // Free the outfile array. It will be replenished before it's use.
        delete [] outfile;

        for (Alignment* &i : XX) {

            // If original filename wasn't empty,
            //      reset it and apply the [contig] tag
            std::string newOutFile =
                    utils::ReplaceString(
                            originalOutFile, "[contig]", i->filename);

            if (!originalOutFile.empty())
            {
                outfile = &newOutFile[0];
                i->filename = originalName;
            }

            // To use innerPerform, alignment must be stored on origAlig
            origAlig = i;
            returnValue = std::max(innerPerform(), returnValue);

            // Delete and nullify tempAlig
            if (tempAlig != origAlig &&
                tempAlig != singleAlig)
                delete tempAlig;
            tempAlig = nullptr;

            // Delete and nullify singleAlig
            if (singleAlig != origAlig)
                delete singleAlig;
            singleAlig = nullptr;

            // Delete and nullify origAlig
            delete origAlig;
            origAlig = nullptr;

        }
        outfile = nullptr;
        origAlig = nullptr;
        return returnValue;
    }
}

/**inline**/ int trimAlManager::innerPerform() {
    origAlig->Cleaning->setTrimTerminalGapsFlag(terminalOnly);
    origAlig->setKeepSequencesFlag(keepSeqs);

    set_window_size();

    if (blockSize != -1)
        origAlig->setBlockSize(blockSize);

    if (!create_or_use_similarity_matrix())
        return 2;

    print_statistics();

    clean_alignment();

    if (singleAlig == nullptr) {
        singleAlig = origAlig;
        origAlig = nullptr;
    }

    postprocess_alignment();

    output_reports();

    save_alignment();

    if (svgStatsOutFile != nullptr)
        svg_stats_out();

    if ((columnNumbering) && (!appearErrors))
        singleAlig->Statistics->printCorrespondence();

    return 0;

}

/**inline**/ void trimAlManager::save_alignment()
{
    if ((outfile != nullptr) && (!appearErrors)) {
        std::string outFileString = std::string(outfile);
        if (singleAlig != nullptr)
        {
            if (!formatManager.saveAlignment(outFileString, oformats, *singleAlig)) {
                appearErrors = true;
            }
        } else debug.report(
                ErrorCode::SomethingWentWrong_reportToDeveloper,
                "Trying to save a nullptr alignment on save_alignment");

    } else if ((stats >= 0) && (!appearErrors))
    {
        std::string emptyString;
        if (singleAlig != nullptr)
            formatManager.saveAlignment(emptyString, oformats, *singleAlig);
        else debug.report(
                ErrorCode::SomethingWentWrong_reportToDeveloper,
                "Trying to save a nullptr alignment on save_alignment");
    }
}

/**inline**/ bool trimAlManager::performCompareset() {

    if (compareset != nullptr && !appearErrors) {

        // Need to check this here, as inside Consistency::perform
        //      origAlig is set and checked
        //      (if no forceselect applied),
        //      and thus, one alignment would be leaked.
        if (infile != nullptr)
        {
            debug.report(ComparesetAndInArgumentsProvided);
            return appearErrors = true;
        }
        CS = new statistics::Consistency();
        if (CS->perform(compareset, formatManager, *this, forceFile))
            appearErrors = true;
    }
    return appearErrors;
}

/**inline**/ void trimAlManager::output_reports()
{
    if ((svgOutFile != nullptr) && (!appearErrors))
    {
        size_t start = origAlig->filename.find_last_of('/');
        if (start == std::string::npos) start = 0;
        size_t end = origAlig->filename.find_last_of('.');
        if (end == std::string::npos) end = origAlig->filename.size();
        std::string filename = utils::ReplaceString(svgOutFile, "[in]", origAlig->filename.substr(start, end-start));
        utils::ReplaceStringInPlace(filename, "[extension]", "svg");
        if (!origAlig->alignmentSummarySVG(*singleAlig, filename.c_str(), 0)) {
            debug.report(ErrorCode::ImpossibleToGenerate, new std::string[1]{"the SVG output file"});
            appearErrors = true;
        }
    }

    if ((htmlOutFile != nullptr) && (!appearErrors))
    {
        size_t start = origAlig->filename.find_last_of('/');
        if (start == std::string::npos) start = 0;
        size_t end = origAlig->filename.find_last_of('.');
        if (end == std::string::npos) end = 0;
        std::string filename = utils::ReplaceString(htmlOutFile, "[in]", origAlig->filename.substr(start, end-start));
        utils::ReplaceStringInPlace(filename, "[extension]", "svg");
        if (!origAlig->alignmentSummaryHTML(*singleAlig, filename.c_str())) {
            debug.report(ErrorCode::ImpossibleToGenerate, new std::string[1]{"the HTML output file"});
            appearErrors = true;
        }
    }
}

/**inline**/ void trimAlManager::print_statistics() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("/**inline**/ void trimAlManager::print_statistics() ");

    if (sgc) {
        origAlig->Statistics->printStatisticsGapsColumns();
        stats++;
    }

    if (sgt) {
        origAlig->Statistics->printStatisticsGapsTotal();
        stats++;
    }

    if (ssc) {
        origAlig->Statistics->printStatisticsConservationColumns();
        stats++;
    }

    if (sst) {
        origAlig->Statistics->printStatisticsConservationTotal();
        stats++;
    }

    if (sident) {
        origAlig->printSeqIdentity();
        stats++;
    }

    if (soverlap) {
        origAlig->printSeqOverlap();
        stats++;
    }

    if (compareset != nullptr) {
        if (sfc)
            statistics::Consistency::printStatisticsFileColumns(*origAlig, origAlig->Statistics->consistency->getValues());
        if (sft)
            statistics::Consistency::printStatisticsFileAcl(*origAlig, origAlig->Statistics->consistency->getValues());
    }
}

/**inline**/ void trimAlManager::svg_stats_out()
{
    origAlig->statSVG(svgStatsOutFile);
}

/**inline**/ bool trimAlManager::create_or_use_similarity_matrix() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("/**inline**/ bool trimAlManager::create_or_use_similarity_matrix() ");
    if ((strict) || (strictplus) || (automated1) || (similarityThreshold != -1.0) || (ssc == 1) || (sst == 1)) {
        similMatrix = new statistics::similarityMatrix();

    // Load Matrix
        if (matrixFile != nullptr)
            similMatrix->loadSimMatrix(matrixFile);

    // Alternative Default Matrix
        else if (alternative_matrix != -1) {
            similMatrix->alternativeSimilarityMatrices(alternative_matrix, origAlig->getAlignmentType());
        }

    // Default Matrices
        else {
            int alignDataType = origAlig->getAlignmentType();
            if (alignDataType & SequenceTypes::AA)
                similMatrix->defaultAASimMatrix();
            else if ((alignDataType == SequenceTypes::DNA) || (alignDataType == SequenceTypes::RNA))
                similMatrix->defaultNTSimMatrix();
            else if ((alignDataType == (SequenceTypes::DNA | SequenceTypes::DEG)) ||
                     (alignDataType == (SequenceTypes::RNA | SequenceTypes::DEG)))
                similMatrix->defaultNTDegeneratedSimMatrix();
        }

    // Check if Matrix has been loaded
        if (!origAlig->Statistics->setSimilarityMatrix(similMatrix)) {
            debug.report(ErrorCode::ImpossibleToProcessMatrix);
            return false;
        }
    }
    return true;
}

/**inline**/ void trimAlManager::clean_alignment() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("inline void trimAlManager::clean_alignment() ");

    // Check alignment is aligned
    if (!origAlig->isFileAligned()) {
        debug.report(ErrorCode::NotAligned, infile);
        exit(ErrorCode::NotAligned);
    }

    // We apply first the cleaning methods that remove sequences
    CleanSequences();

    // We apply the cleaning methods that remove columns / residues
    // Combination between Automatic and Non Automatic are not allowed
    if (automatedMethodCount)
        CleanResiduesAuto();
    else
        CleanResiduesNonAuto();

}

/**inline**/ void trimAlManager::postprocess_alignment() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("/**inline**/ void trimAlManager::postprocess_alignment()");
    
    // Only terminal
    if (terminalOnly)
        singleAlig->Cleaning->removeOnlyTerminal();

    // Complementary
    if (getComplementary)
        singleAlig->Cleaning->computeComplementaryAlig(true, false);

    // ComplementarySeq
    if (getComplementarySeq)
        singleAlig->Cleaning->computeComplementaryAlig(false, true);

    // Backtranslate
    if (backtransFile != nullptr) {
        tempAlig = backtranslationAlig->getTranslationCDS(singleAlig);
        if (singleAlig != origAlig)
            delete singleAlig;
        singleAlig = tempAlig;
        tempAlig = nullptr;
    }

}

/**inline**/ void trimAlManager::CleanSequences() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("/**inline**/ void trimAlManager::CleanSequences() ");

    // Clustering by number of clusters
    if (clusters != -1) {
        tempAlig = origAlig->Cleaning->getClustering(
                origAlig->Cleaning->getCutPointClusters(clusters)
        );

    // Clustering by Max Identity
    } else if (maxIdentity != -1) {
        tempAlig = origAlig->Cleaning->getClustering(maxIdentity);

    // Remove Sequences
    } else if (delSequences != nullptr) {
        tempAlig = origAlig->Cleaning->removeSequences(
                delSequences,       // Informative array
                1,                  // Position to start. 
                                    // Position on number 0 indicates the number of elements.
                delSequences[0],    // Num of sequences in array
                false               // Complementary? Legacy
        );

    // Remove by Overlap
    } else if ((residuesOverlap != -1) && (sequenceOverlap != -1)) {
        tempAlig = origAlig->Cleaning->cleanSpuriousSeq(
                residuesOverlap,
                sequenceOverlap / 100.0F,
                /* getComplementary*/ false
        );
    } else if (removeDuplicates) {
        origAlig->Cleaning->removeDuplicates();
    }


    // We'll use singleAlig as input for the next cleaning step.
    // If tempAlig is not null it means that we've performed a sequence trim
    // If not, we move the origAlig to singleAlig
    if (tempAlig) {
        singleAlig = tempAlig->Cleaning->cleanNoAllGaps(false);

        delete tempAlig;
        tempAlig = nullptr;

        delete singleAlig->Statistics->gaps;
        singleAlig->Statistics->gaps = nullptr;

        delete singleAlig->Statistics->similarity;
        singleAlig->Statistics->similarity = nullptr;

    } else {
        singleAlig = origAlig;
    }
}

/**inline**/ void trimAlManager::CleanResiduesAuto() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("/**inline**/ void trimAlManager::CleanResiduesAuto() ");
    // Here we use singleAlig as source alignment as the previous step,
    // CleanSequences initializes singleAlig.
    // singleAlig can be a derived alignment from origAlig or be origAlig itself

    if (automated1) {
        if (singleAlig->Cleaning->selectMethod() == GAPPYOUT)
            gappyout = true;
        else
            strict = true;
    }
    if (nogaps) {
        tempAlig = singleAlig->Cleaning->cleanGaps(0, 0, /* getComplementary*/ false);
    } else if (noallgaps) {
        tempAlig = singleAlig->Cleaning->cleanNoAllGaps(/* getComplementary*/ false);
    } else if (gappyout) {
        tempAlig = singleAlig->Cleaning->clean2ndSlope(/* getComplementary*/ false);
    } else if (strict) {
        tempAlig = singleAlig->Cleaning->cleanCombMethods(/* getComplementary*/ false, false);
    } else if (strictplus) {
        tempAlig = singleAlig->Cleaning->cleanCombMethods(/* getComplementary*/ false, true);
    }

    // Move the new formed alignment to the variable singleAlig
    if (tempAlig) {
        // If singleAlig was an created, we delete it from memory.
        if (singleAlig != origAlig)
            delete singleAlig;

        singleAlig = tempAlig;
        tempAlig = nullptr;
    }
}

/**inline**/ void trimAlManager::CleanResiduesNonAuto() {
    // Here we use singleAlig as source alignment as the previous step,
    // CleanSequences initializes singleAlig.
    // singleAlig can be a derived alignment from origAlig or origAlig itself

    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("/**inline**/ void trimAlManager::CleanResiduesNonAuto() ");

    if (delColumns != nullptr) {
        for (int i = 0 ; i <= delColumns[0]; i++)
        {
            if (delColumns[i] >= singleAlig->getNumAminos()) {
                debug.report(
                        ErrorCode::SelectOnlyAccepts,
                        new std::string[2]{"-selectcols", "residues"}
                );
                appearErrors = true;
            }
        }
         if (!appearErrors)
            tempAlig = singleAlig->Cleaning->removeColumns(
                    delColumns,
                    1,
                    delColumns[0],
                    /* getComplementary*/ false
            );

    } else {

        // Consistency will be applied prior to gap or similarity
        if (consistencyThreshold != -1.0F) {
            tempAlig = singleAlig->Cleaning->cleanCompareFile(
                    consistencyThreshold,
                    conservationThreshold,
                    origAlig->Statistics->consistency->getValues(),
                    /* getComplementary*/ false
            );
            if (singleAlig != origAlig) delete singleAlig;
            singleAlig = tempAlig;
            tempAlig = nullptr;
        }

        if (similarityThreshold != -1.0F) {

            // Clean by similarity AND gaps
            if (gapThreshold != -1.0F) {
                tempAlig = singleAlig->Cleaning->clean(
                        conservationThreshold,
                        gapThreshold,
                        similarityThreshold,
                        /* getComplementary*/ false
                );
                // Clean only by similarity
            } else {
                tempAlig = singleAlig->Cleaning->cleanConservation(
                        conservationThreshold,
                        similarityThreshold,
                        /* getComplementary*/ false
                );
            }
            // Clean only by gaps
        } else if (gapThreshold != -1.0F) {

            tempAlig = singleAlig->Cleaning->cleanGaps(
                    conservationThreshold,
                    gapThreshold,
                    /* getComplementary*/ false
            );

        }
    }

    if (tempAlig) {
        if (singleAlig != origAlig)
            delete singleAlig;
        singleAlig = tempAlig;
        tempAlig = nullptr;
    }
}

/**inline**/ void trimAlManager::set_window_size() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("inline void trimAlManager::set_window_size() ");
    if (windowSize != -1) {
        gapWindow = windowSize;
        similarityWindow = windowSize;
    } else {
        if (gapWindow == -1)
            gapWindow = 0;
        if (similarityWindow == -1)
            similarityWindow = 0;
    }
    origAlig->setWindowsSize(gapWindow, similarityWindow);
}

/**inline**/ void trimAlManager::delete_variables() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("inline void trimAlManager::delete_variables() ");

    if (tempAlig == origAlig || tempAlig == singleAlig)
        tempAlig = nullptr;
    if (singleAlig == origAlig)
        singleAlig = nullptr;

    delete tempAlig;
    tempAlig = nullptr;

    delete singleAlig;
    singleAlig = nullptr;

    delete origAlig;
    origAlig = nullptr;

    delete[] compareAlignmentsArray;
    compareAlignmentsArray = nullptr;

    delete similMatrix;
    similMatrix = nullptr;

    delete[] delColumns;
    delColumns = nullptr;

    delete[] delSequences;
    delSequences = nullptr;

    delete[] filesToCompare;
    filesToCompare = nullptr;

    delete[] outfile;
    outfile = nullptr;
    delete[] htmlOutFile;
    htmlOutFile = nullptr;
    delete[] svgOutFile;
    svgOutFile = nullptr;
    delete[] svgStatsOutFile;
    svgStatsOutFile = nullptr;

    delete[] infile;
    infile = nullptr;
    delete[] matrixFile;
    matrixFile = nullptr;

    delete [] forceFile;
    forceFile = nullptr;
    delete [] backtransFile;
    backtransFile = nullptr;
    delete backtranslationAlig;
    backtranslationAlig = nullptr;

    delete vcfs;
    vcfs = nullptr;

    delete CS;
    CS = nullptr;

}

/**inline**/ void trimAlManager::menu() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void trimAlManager::menu(void) ");

    // Simple idiom to allow including a text file into a char array at compile time
    // Explanation: https://stackoverflow.com/a/25021520
    std::string menu = {
        #include "RawText/menu.txt"
    };

    utils::ReplaceStringInPlace(menu, "[iformat]",  formatManager.getInputFormatsAvailable());
    utils::ReplaceStringInPlace(menu, "[oformat]",  formatManager.getOutputFormatsAvailable());
    utils::ReplaceStringInPlace(menu, "[version]",  VERSION);
    utils::ReplaceStringInPlace(menu, "[revision]", REVISION);
    utils::ReplaceStringInPlace(menu, "[build]",    BUILD);
    utils::ReplaceStringInPlace(menu, "[authors]",  AUTHORS);

    using namespace utils::TerminalColors;

    utils::ReplaceStringInPlace(menu, "[b][u]", colors[BOLD] + colors[UNDERLINE]);
    utils::ReplaceStringInPlace(menu, "[b]",    colors[BOLD]);
    utils::ReplaceStringInPlace(menu, "[u]",    colors[UNDERLINE]);
    utils::ReplaceStringInPlace(menu, "[r]",    colors[RESET]);

    std::cout << menu;
}

/**inline**/ void trimAlManager::examples() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void trimAlManager::examples(void) ");

    // Simple idiom to allow including a text file into a char array at compile time
    // Explanation: https://stackoverflow.com/a/25021520
    std::string examples = {
        #include "RawText/examples.txt"
    };

    using namespace utils::TerminalColors;

    utils::ReplaceStringInPlace(examples, "[b][u]", colors[BOLD] + colors[UNDERLINE]);
    utils::ReplaceStringInPlace(examples, "[b]",    colors[BOLD]);
    utils::ReplaceStringInPlace(examples, "[u]",    colors[UNDERLINE]);
    utils::ReplaceStringInPlace(examples, "[r]",    colors[RESET]);

    std::cout << examples;
}

FormatHandling::FormatManager &trimAlManager::getFormatManager() {
    return this->formatManager;
}



