#include <TimerFactory.h>
//
// Created by bioinfo on 8/06/17.
//

#include <string>

#include "../include/trimalManager.h"
#include "Statistics/statisticsConsistency.h"
#include "defines.h"
#include "values.h"
#include "reportsystem.h"
#include "Statistics/statisticsConservation.h"

#include "vcf_statish.h"

//void menu();
//
//void examples();

trimAlManager::trimAlManager() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("trimAlManager::trimAlManager() ");

}

trimAlManager::~trimAlManager() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("trimAlManager::~trimAlManager() ");
    delete_variables();
}

void trimAlManager::parseArguments(int argc, char *argv[]) {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void trimAlManager::parseArguments(int argc, char *argv[]) ");

    // Call to verbosity argument has to be done
    //  before any other argument.
    verbosity_argument(&argc, argv);

    // Check if there is only a call to trimAl,
    //  if so, show the menu and the examples
    if (argc == 1) {
        menu();
        examples();
        exit(0);
    }

    // Check every provided argument
    //  with each argument method
    for (int i = 1; i < argc; i++) {
        if (appearErrors) break;

        help_arguments(&argc, argv, &i);

        if (in_argument(&argc, argv, &i)) continue;
        if (vcf_argument(&argc, argv, &i)) continue;
        if (out_argument(&argc, argv, &i)) continue;
        if (html_out_argument(&argc, argv, &i)) continue;
        if (timetracker_out_argument(&argc, argv, &i)) continue;
        if (svg_out_argument(&argc, argv, &i)) continue;
        if (svg_stats_argument(&argc, argv, &i)) continue;
        if (out_format_arguments(&argc, argv, &i)) continue;
        if (matrix_argument(&argc, argv, &i)) continue;

        if (compareset_argument(&argc, argv, &i)) continue;
        if (force_select_argument(&argc, argv, &i)) continue;
        if (back_trans_argument(&argc, argv, &i)) continue;

        if (gap_threshold_argument(&argc, argv, &i)) continue;
        if (similarity_threshold_argument(&argc, argv, &i)) continue;
        if (consistency_threshold_argument(&argc, argv, &i)) continue;

        if (conservation_threshold_argument(&argc, argv, &i)) continue;
        if (select_cols_argument(&argc, argv, &i)) continue;

        if (no_gaps_argument(&argc, argv, &i)) continue;
        if (no_all_gaps_argument(&argc, argv, &i)) continue;

        if (keep_seqs_argument(&argc, argv, &i)) continue;
        if (keep_header_argument(&argc, argv, &i)) continue;

        if (gappy_out_argument(&argc, argv, &i)) continue;
        if (strict_argument(&argc, argv, &i)) continue;
        if (strict_plus_argument(&argc, argv, &i)) continue;
        if (automated1_argument(&argc, argv, &i)) continue;

        if (residue_overlap_argument(&argc, argv, &i)) continue;
        if (sequence_overlap_argument(&argc, argv, &i)) continue;

        if (seqs_select_argument(&argc, argv, &i)) continue;

        if (max_identity_argument(&argc, argv, &i)) continue;
        if (clusters_argument(&argc, argv, &i)) continue;

        if (terminal_only_argument(&argc, argv, &i)) continue;

        if (window_argument(&argc, argv, &i)) continue;
        if (gap_window_argument(&argc, argv, &i)) continue;
        if (similarity_window_argument(&argc, argv, &i)) continue;
        if (consistency_window_argument(&argc, argv, &i)) continue;

        if (block_argument(&argc, argv, &i)) continue;
        if (stats_arguments(&argc, argv, &i)) continue;

        if (complementary_argument(&argc, argv, &i)) continue;

        if (col_numbering_argument(&argc, argv, &i)) continue;
        if (split_by_stop_codon_argument(&argc, argv, &i)) continue;
        if (ignore_stop_codon_argument(&argc, argv, &i)) continue;

        // Skip the verbosity option, as it was checked before the loop
        if (!strcmp(argv[i], "--verbosity") || !strcmp(argv[i], "-v")) {
            i++;
            continue;
        }

        debug.report(ErrorCode::ParameterNotFoundOrRepeated, argv[i]);
        appearErrors = true;
        break;
    }

    // Check if we've provided multiple correct arguments
    // but no input file is provided.
    if (infile == nullptr && compareset == -1) {
        debug.report(ErrorCode::NoInputFile);
        exit(0);
    }
}

inline bool trimAlManager::check_arguments_incompatibilities() {
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

    return appearErrors;
}

inline void trimAlManager::verbosity_argument(const int *argc, char *argv[]) {
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

inline void trimAlManager::help_arguments(const int *argc, char **argv, int *i) {
    if (!strcmp(argv[*i], "-h") || !strcmp(argv[*i], "-help")) {
        menu();
        examples();
        exit(0); // We don't want to continue if we show the help.
    }

//    if (!strcmp(argv[*i], "-examples")) {
//        examples();
//        delete_variables();
//        exit(0); // We don't want to continue if we show the help.
//    }

    if (!strcmp(argv[*i], "--version")) {
        std::cout << VERSION << std::endl;
        exit(0); // We don't want to continue if we show the examples.
    }

    if (!strcmp(argv[*i], "-lf") || !strcmp(argv[*i], "--listformats")) {
        cout << "Input Formats:  " << ReadWriteMachine.getInputFormatsAvailable() << endl;
        cout << "Output Formats: " << ReadWriteMachine.getOutputFormatsAvailable() << endl;
        exit(0);
    }
}

inline bool trimAlManager::in_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-in") && ((*i) + 1 != *argc) && (infile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        infile = new char[argumentLength + 1];
        strcpy(infile, argv[*i]);
        if ((origAlig = ReadWriteMachine.loadAlignment(infile)) == nullptr) {
            appearErrors = true;
        }
        return true;
    }
    return false;
}

inline bool trimAlManager::vcf_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-vcf") && ((*i) + 1 != *argc) && (infile == nullptr)) {
        vcfs = new std::vector<std::string>();
        while (((*i) + 1 != *argc)) {
            ++*i;
            if (argv[*i][0] == '-') {
                --*i;
                break;
            }
            vcfs->emplace_back(argv[*i]);
        }
        return true;

    }
    return false;
}

inline bool trimAlManager::out_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-out")) && ((*i) + 1 != *argc) && (outfile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        outfile = new char[argumentLength + 1];
        strcpy(outfile, argv[*i]);
        return true;
    }
    return false;
}

inline bool trimAlManager::html_out_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-htmlout")) && ((*i) + 1 != *argc) && (htmlOutFile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        htmlOutFile = new char[argumentLength + 1];
        strcpy(htmlOutFile, argv[*i]);
        return true;
    }
    return false;
}

inline bool trimAlManager::timetracker_out_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-timetrackerout")) && ((*i) + 1 != *argc)) {
        (*i) += 1;
        return true;
    }
    return false;
}

inline bool trimAlManager::svg_out_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-svgout")) && ((*i) + 1 != *argc) && (svgOutFile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        svgOutFile = new char[argumentLength + 1];
        strcpy(svgOutFile, argv[*i]);
        return true;
    }
    return false;
}

inline bool trimAlManager::svg_stats_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-svgstats")) && ((*i) + 1 != *argc) && (svgStatsOutFile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        svgStatsOutFile = new char[argumentLength + 1];
        strcpy(svgStatsOutFile, argv[*i]);
        return true;
    }
    return false;
}

bool trimAlManager::out_format_arguments(const int *argc, char *argv[], int *i) {
    {
        if (!strcmp(argv[*i], "-formats")) {
            if ((*i + 1) == *argc) {
                debug.report(ErrorCode::NoFormatsSpecified);
            }
            while (++(*i) != *argc && argv[*i][0] != '-')
                oformats.emplace_back(argv[*i]);
            (*i)--;
            return true;
        }
            // Option -clustal
        else if (!strcmp(argv[*i], "-clustal")) {
            oformats.emplace_back("clustal");
            return true;
        }

            // Option -fasta
        else if (!strcmp(argv[*i], "-fasta")) {
            oformats.emplace_back("fasta");
            return true;
        }

            // Option -fasta-m10
        else if (!strcmp(argv[*i], "-fasta_m10")) {
            oformats.emplace_back("fasta_m10");
            return true;
        }

            // Option -nbrf
        else if (!strcmp(argv[*i], "-nbrf")) {
            oformats.emplace_back("pir");
            return true;
        }

            // Option -nexus
        else if (!strcmp(argv[*i], "-nexus")) {
            oformats.emplace_back("nexus");
            return true;
        }

            // Option -mega
        else if (!strcmp(argv[*i], "-mega")) {
            oformats.emplace_back("mega");
            return true;
        }

            // Option -phylip3.2
        else if (!strcmp(argv[*i], "-phylip3.2")) {
            oformats.emplace_back("phylip32");
            return true;
        }

            // Option -phylip3.2-m10
        else if (!strcmp(argv[*i], "-phylip3.2_m10")) {
            oformats.emplace_back("phylip32_m10");
            return true;
        }

            // Option -phylip
        else if (!strcmp(argv[*i], "-phylip")) {
            oformats.emplace_back("phylip40");
            return true;
        }

            // Option -phylip-m10
        else if (!strcmp(argv[*i], "-phylip_m10")) {
            oformats.emplace_back("phylip40_m10");
            return true;
        }

            // Option -phylip_paml
        else if (!strcmp(argv[*i], "-phylip_paml")) {
            oformats.emplace_back("phylip_paml");
            return true;
        }

            // Option -phylip_paml-m10
        else if (!strcmp(argv[*i], "-phylip_paml_m10")) {
            oformats.emplace_back("phylip_paml_m10");
            return true;
        }
    }
    return false;
}

inline bool trimAlManager::matrix_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-matrix") && ((*i) + 1 != *argc) && (matrixFile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        matrixFile = new char[argumentLength + 1];
        strcpy(matrixFile, argv[*i]);
        return true;
    } else if (!strcmp(argv[*i], "--alternative_matrix") && ((*i) + 1 != *argc) && (alternative_matrix == -1)) {
        i++;
        if (!strcmp(argv[*i], "degenerated_nt_identity"))
            alternative_matrix = 1;
        else {
            debug.report(ErrorCode::AlternativeMatrixNotRecognized, argv[*i]);
            appearErrors = true;
        }
    }
    return false;
}

inline bool trimAlManager::compareset_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-compareset") && ((*i) + 1 != *argc) && (compareset == -1)) {
        compare.open(argv[++*i], ifstream::in);
        if (!compare) {
            debug.report(ErrorCode::ReferenceFileNotLoaded, argv[*i]);
            appearErrors = true;
        }
        compare.close();

        compareset = *i;

        return true;
    }
    return false;
}

inline bool trimAlManager::force_select_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-forceselect") && ((*i) + 1 != *argc) && (forceFile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        forceFile = new char[argumentLength + 1];
        strcpy(forceFile, argv[*i]);
        if ((origAlig = ReadWriteMachine.loadAlignment(forceFile)) == nullptr) {
            debug.report(ErrorCode::AlignmentNotLoaded, forceFile);
            appearErrors = true;
        }
        return true;
    }
    return false;
}

inline bool trimAlManager::back_trans_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-backtrans") && ((*i) + 1 != *argc) && (backtransFile == nullptr)) {
        argumentLength = strlen(argv[++*i]);
        backtransFile = new char[argumentLength + 1];
        strcpy(backtransFile, argv[*i]);

        if ((backtranslationAlig = ReadWriteMachine.loadAlignment(backtransFile)) == nullptr) {
            debug.report(ErrorCode::AlignmentNotLoaded, backtransFile);
            appearErrors = true;
        }
        return true;
    }
    return false;
}

inline bool trimAlManager::gap_threshold_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-gapthreshold") || !strcmp(argv[*i], "-gt")) && ((*i) + 1 != *argc) && (gapThreshold == -1)) {
        if (utils::isNumber(argv[++*i])) {
            gapThreshold = 1. - atof(argv[*i]);
            if ((gapThreshold < 0) || (gapThreshold > 1)) {
                debug.report(ErrorCode::GapThresholdOutOfRange);
                appearErrors = true;
            }
        } else {
            debug.report(ErrorCode::GapThresholdNotRecognized);
            appearErrors = true;
        }
        return true;
    }
    return false;
}

inline bool trimAlManager::similarity_threshold_argument(const int *argc, char *argv[], int *i) {
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
        return true;
    }
    return false;
}

inline bool trimAlManager::consistency_threshold_argument(const int *argc, char *argv[], int *i) {
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
        return true;
    }
    return false;
}

inline bool trimAlManager::conservation_threshold_argument(const int *argc, char *argv[], int *i) {
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

        return true;
    }
    return false;
}

bool trimAlManager::select_cols_argument(const int *argc, char *argv[], int *i) {

    if (!strcmp(argv[*i], "-selectcols") &&
        !selectCols &&
        (*i + 3) < *argc &&
        !strcmp(argv[++(*i)], "{") &&
        !strcmp(argv[(*i) + 2], "}")) {
        if ((delColumns = utils::readNumbers(argv[++(*i)])) == nullptr) {
            debug.report(ErrorCode::SelectColsNotRecognized);
            appearErrors = true;
        } else selectCols = true;
        (*i)++;

        return true;
    }
    return false;
}

inline bool trimAlManager::no_gaps_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-nogaps") && (!nogaps)) {
        nogaps = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::no_all_gaps_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-noallgaps") && (!noallgaps)) {
        noallgaps = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::keep_seqs_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-keepseqs") && (!keepSeqs)) {
        keepSeqs = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::keep_header_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-keepheader") && (!ReadWriteMachine.keepHeader)) {
        ReadWriteMachine.keepHeader = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::gappy_out_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-gappyout") && (!strict)) {
        gappyout = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::strict_argument(const int *argc, char *argv[], int *i) {
    if (!strcmp(argv[*i], "-strict") && (!strict)) {
        strict = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::strict_plus_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-strictplus")) && (!strictplus)) {
        strictplus = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::automated1_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-automated1")) && (!automated1)) {
        automated1 = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::residue_overlap_argument(const int *argc, char *argv[], int *i) {
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
        return true;
    }
    return false;
}

inline bool trimAlManager::sequence_overlap_argument(const int *argc, char *argv[], int *i) {
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
        return true;
    }
    return false;
}

bool trimAlManager::seqs_select_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-selectseqs")) && !selectSeqs && ((*i + 3) < *argc) && (!strcmp(argv[++*i], "{")) && (!strcmp(argv[*i + 2], "}"))) {
        if ((delSequences = utils::readNumbers(argv[++*i])) == nullptr) {
            debug.report(ErrorCode::SelectSeqsNotRecognized);
            appearErrors = true;
        } else selectSeqs = true;
        (*i)++;
        return true;
    }
    return false;
}

inline bool trimAlManager::max_identity_argument(const int *argc, char *argv[], int *i) {
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

        return true;
    }
    return false;
}

inline bool trimAlManager::clusters_argument(const int *argc, char *argv[], int *i) {
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

        return true;
    }
    return false;
}

inline bool trimAlManager::terminal_only_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-terminalonly")) && (!terminalOnly)) {
        terminalOnly = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::window_argument(const int *argc, char *argv[], int *i) {
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
        return true;
    }
    return false;
}

inline bool trimAlManager::gap_window_argument(const int *argc, char *argv[], int *i) {
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
        return true;
    }
    return false;
}

inline bool trimAlManager::similarity_window_argument(const int *argc, char *argv[], int *i) {
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

        return true;
    }
    return false;
}

inline bool trimAlManager::consistency_window_argument(const int *argc, char *argv[], int *i) {
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
        return true;
    }
    return false;
}

inline bool trimAlManager::block_argument(const int *argc, char *argv[], int *i) {
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
        return true;
    }
    return false;
}

inline bool trimAlManager::stats_arguments(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-sgc"))) {
        if (!sgc) {
            sgc = true;
            stats--;
        } else return false;
    } else if ((!strcmp(argv[*i], "-sgt"))) {
        if (!sgt) {
            sgt = true;
            stats--;
        } else return false;
    } else if ((!strcmp(argv[*i], "-ssc"))) {
        if (!scc) {
            scc = true;
            stats--;
        } else return false;
    } else if ((!strcmp(argv[*i], "-sst"))) {
        if (!sct) {
            sct = true;
            stats--;
        } else return false;
    } else if ((!strcmp(argv[*i], "-sident")) && (!sident)) {
        sident = true;
        stats--;
    } else if ((!strcmp(argv[*i], "-soverlap")) && (!soverlap)) {
        soverlap = true;
        stats--;
    } else if ((!strcmp(argv[*i], "-sfc"))) {
        if (!sfc) {
            sfc = true;
            stats--;
        } else return false;
    } else if ((!strcmp(argv[*i], "-sft"))) {
        if (!sft) {
            sft = true;
            stats--;
        } else return false;
    } else
        return false;
    return true;
}

inline bool trimAlManager::complementary_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-complementary")) && !getComplementary) {
        getComplementary = true;
        return true;
    }
    return false;

}

inline bool trimAlManager::col_numbering_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-colnumbering")) && !columnNumbering) {
        columnNumbering = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::split_by_stop_codon_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-splitbystopcodon")) && !splitByStopCodon) {
        splitByStopCodon = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::ignore_stop_codon_argument(const int *argc, char *argv[], int *i) {
    if ((!strcmp(argv[*i], "-ignorestopcodon")) && !ignoreStopCodon) {
        ignoreStopCodon = true;
        return true;
    }
    return false;
}

bool trimAlManager::processArguments(char *argv[]) {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool trimAlManager::processArguments(char *argv[]) ");
    if (!appearErrors) {
        automatedMethodCount = nogaps + noallgaps + gappyout + strict + strictplus + automated1;

        check_arguments_incompatibilities();

        check_arguments_needs(argv);

    }

    return appearErrors;
}

inline bool trimAlManager::check_inFile_incompatibilities() {
    if (infile != nullptr) {
        if ((sfc) || (sft) || (consistencyThreshold != -1)) {
            debug.report(ErrorCode::InFileComparisonStatistics);
            appearErrors = true;
            i++;
        }
        if (compareset != -1) {
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

inline bool trimAlManager::check_select_cols_and_seqs_incompatibilities() {
    if (selectCols || selectSeqs) {
        if ((clusters != -1) || (maxIdentity != -1)) {
            debug.report(ErrorCode::OnlyOneSequencesSelectionMethodAllowed);
            appearErrors = true;
        }

        if (selectCols)
            if (blockSize != -1) {
                debug.report(ErrorCode::IncompatibleArguments,
                             new std::string[2]{"-selectcols", "-block"});
                appearErrors = true;
            }

        if (selectSeqs)

            for (int i = 0; i < delSequences[0]; i++)
            {
                if (delSequences[i] >= origAlig->getNumSpecies()) {
                    debug.report(ErrorCode::SelectOnlyAccepts,
                                 new string[2]{"-selectseqs", "sequences"});
                    appearErrors = true;
                    break;
                }
            }
    }
    return appearErrors;
}

inline bool trimAlManager::check_thresholds_incompatibilities() {
    //NOTE consistencyThreshold was not present on this Incompatibilities. Is this correct?
    if ((gapThreshold != -1) || (similarityThreshold != -1) || (consistencyThreshold != -1)) {

        if (automatedMethodCount) {
            debug.report(ErrorCode::CombinationAmongTrimmingMethods);
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

inline bool trimAlManager::check_automated_methods_incompatibilities() {
    if (automatedMethodCount) {
        if ((windowSize != -1) || (gapWindow != -1) || (similarityWindow != -1)) {
            debug.report(ErrorCode::CombinationAmongTrimmingMethods);
            appearErrors = true;
        }

        if (automatedMethodCount > 1) {
            debug.report(ErrorCode::CombinationAmongTrimmingMethods);
            appearErrors = true;
        }
    }
    return appearErrors;
}

inline bool trimAlManager::check_max_identity_incompatibilities() {
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

inline bool trimAlManager::check_clusters_incompatibilities() {
    if (clusters != -1) {
        if (clusters > origAlig->getNumSpecies()) {
            debug.report(ErrorCode::MoreClustersThanSequences);
            appearErrors = true;
        }
    }
    return appearErrors;
}

inline bool trimAlManager::check_windows_incompatibilities() {
    if (windowSize != -1) {
        if (consistencyWindow != -1 || gapWindow != -1 || similarityWindow != -1) {
            debug.report(ErrorCode::GeneralAndSpecificWindows);
            appearErrors = true;
        }
    }
    return appearErrors;
}

inline bool trimAlManager::check_stats_incompatibilities() {
    if (stats < 0) {
        if (columnNumbering) {
            debug.report(ErrorCode::StatisticsArgumentIncompatibilities, new string[1]{"-colnumbering"});
            appearErrors = true;
        }
    }
    return appearErrors;
}

inline bool trimAlManager::check_arguments_needs(char *argv[]) {
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
    check_multiple_files_comparison(argv);
    check_block_size();
    check_backtranslations();
    check_coding_sequences_type();
    check_and_prepare_coding_sequence();
    check_backtranslation_infile_names_corresponde();
    check_compareset_window_argument();
    check_output_format();
    return appearErrors;
}

inline bool trimAlManager::check_codon_behaviour_incompatibility() {
    if ((!appearErrors) && (ignoreStopCodon) && (splitByStopCodon)) {
        debug.report(ErrorCode::IncompatibleArguments, new std::string[2]{"-ignorestopcodon", "-splitbystopcodon"});
        appearErrors = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::check_combinations_among_thresholds_incompatibility() // TODO is this ok?
{
    if ((consistencyThreshold != -1) && (conservationThreshold != -1) && (!appearErrors)) {

        if ((gapThreshold != -1) || (similarityThreshold != -1)) {
            debug.report(ErrorCode::CombinationAmongThresholdsMethods);
            appearErrors = true;
            return true;
        }
    }
    return false;
}

inline bool trimAlManager::check_automated_manual_incompatibilities() {
    if ((getComplementary) && (!appearErrors))
        if (!automatedMethodCount && // Are we not using an automated method?
            (gapThreshold == -1) && (conservationThreshold == -1) && (similarityThreshold == -1) && // Neither a threshold method.
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

inline bool trimAlManager::check_force_selection() {
    if (!appearErrors) {
        if ((compareset == -1) && (forceFile != nullptr)) {
            debug.report(ErrorCode::ForceFileWithoutCompareDataset);
            appearErrors = true;
            return true;
        }
    }
    return false;
}

inline bool trimAlManager::check_input_file_with_coding_sequences_argument() {
    if ((!appearErrors) && (infile == nullptr) && (compareset == -1) && (forceFile == nullptr) && (backtransFile != nullptr)) {
        debug.report(ErrorCode::BacktranslationWithoutMainAlignment);
        appearErrors = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::check_file_aligned() {
    if ((!appearErrors) && (infile != nullptr)) {

        if ((automatedMethodCount || // Are we using an automated method?
             (gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1) || // Are we using any threshold?
             (selectCols) || (selectSeqs) || (residuesOverlap != -1) || (sequenceOverlap != -1) || // Are we selecting columns or sequences?
             (stats < 0))  // Are we asking for any stats?
            &&
            (!origAlig->isFileAligned())) // Then we need the alignment to be aligned;
        {
            debug.report(ErrorCode::NotAligned, new std::string[1]{infile});
            appearErrors = true;
            return true;
        }
    }
    return false;
}

inline bool trimAlManager::check_similarity_matrix() {
    if ((matrixFile != nullptr) && (!appearErrors)) {
        if ((!strict) && (!strictplus) && (!automated1) && (similarityThreshold == -1) && (!scc) && (!sct)) {
            debug.report(ErrorCode::MatrixGivenWithNoMethodToUseIt);
            appearErrors = true;
            return true;
        }

        // TODO this is an incompatibility.
        if ((gapWindow != -1) || ((compareset == -1) && (consistencyWindow != -1))) {
            debug.report(ErrorCode::SimilarityMatrixNotCompatibleWindow);
            appearErrors = true;
            return true;
        }
    }
    return false;
}

inline bool trimAlManager::check_outputs_coincidence() {

    std::array<char *, 4> outFiles
            {{
                     htmlOutFile,
                     outfile,
                     svgOutFile,
                     svgStatsOutFile
             }
            };
    std::array<string, 4> outFilesNames
            {{
                     "html (-htmlout)",
                     "output (-out)",
                     "svg (-svgout)",
                     "svg stats (-svgstats)"
             }
            };

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

inline bool trimAlManager::check_col_numbering() {
    if ((columnNumbering) && (!appearErrors)) {
        if ((!automatedMethodCount) && // Are we not using any automated method?
            (gapThreshold == -1) && (conservationThreshold == -1) && (similarityThreshold == -1) && (consistencyThreshold == -1) && // Are we not using any threshold?
            (!selectCols) && (!selectSeqs)) // Neither selecting any column or sequence?
        {
            debug.report(ErrorCode::TrimmingMethodNeeded, new std::string[1]{"-colnumbering"});
            appearErrors = true;
            return true;
        }
    }
    return false;
}

inline bool trimAlManager::check_residue_and_sequence_overlap() {
    if (!appearErrors) {
        if ((residuesOverlap != -1) && (sequenceOverlap == -1)) {
            debug.report(ErrorCode::SequenceAndResiduesOverlapMutuallyNeeded, new string[1]{"residues overlap"});
            appearErrors = true;
            return true;
        } else if ((residuesOverlap == -1) && (sequenceOverlap != -1)) {
            debug.report(ErrorCode::SequenceAndResiduesOverlapMutuallyNeeded, new string[1]{"sequences overlap"});
            appearErrors = true;
            return true;
        }
    }
    return false;
}

inline bool trimAlManager::check_output_relevance() {
    if (((htmlOutFile != nullptr) ||
         (svgOutFile != nullptr) ||
         (svgStatsOutFile != nullptr)) && (!appearErrors)) {
        if (!automatedMethodCount && // Are we not using any automated method?
            (gapThreshold == -1) && (conservationThreshold == -1) && (similarityThreshold == -1) && (consistencyThreshold == -1) && // Neither using thresholds
            (!selectCols) && (!selectSeqs) && (residuesOverlap == -1) && (sequenceOverlap == -1) && // Neither selecting columns or sequences?
            (maxIdentity == -1) && (clusters == -1)) // Neither using other selecting methods?
        {
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

inline bool trimAlManager::check_output_file_with_statistics() {
    // NOTE should this function only take in mind if statistics are requested and outFile is not defined?
    // Either we are using 'some' trimming methods or not, the incompatibility to output multiple info (alignment + stats) is present.
    if ((stats < 0) && (!appearErrors)) {
        stats--; //NOTE Why?

        if (((automatedMethodCount) || // If we are using an automated method
             (gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1)) // Or a manual threshold

            && (outfile == nullptr)) // We need the outFile to be specified.
        {
            debug.report(ErrorCode::OutFileNeededWhenPrintingStatistics);
            appearErrors = true;
            return true;
        }
    }
    return false;
}

inline bool trimAlManager::check_multiple_files_comparison(char *argv[]) {
    if ((compareset != -1) && (!appearErrors)) {

        statisticsConsistency *CS = new statisticsConsistency();
        CS->perform(argv[compareset], ReadWriteMachine, *this, forceFile);
    }
}

inline bool trimAlManager::check_block_size() {
    if ((!appearErrors) && (origAlig->getNumAminos() < (blockSize / 4))) {
        debug.report(ErrorCode::BlocksizeTooBig, new string[1]{std::to_string(origAlig->getNumAminos() / 4)});
        appearErrors = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::check_backtranslations() {
    if (!appearErrors) {
        if (backtransFile == nullptr) {
            if (splitByStopCodon) {
                debug.report(ErrorCode::ParemeterOnlyOnBacktranslation, new string[1]{"-splitbystopcodon"});
                appearErrors = true;
                return true;
            }
            if (ignoreStopCodon) {
                debug.report(ErrorCode::ParemeterOnlyOnBacktranslation, new string[1]{"-ignorestopcodon"});
                appearErrors = true;
                return true;
            }
        } else if (!origAlig->isFileAligned()) {
            debug.report(ErrorCode::ProteinAlignmentMustBeAligned);
            appearErrors = true;
            return true;
        }
    }
    return false;
}

inline bool trimAlManager::check_coding_sequences_type() {
    if ((!appearErrors) && (backtransFile != nullptr) && // Is there a backtranslation file?
        (!backtranslationAlig->getAlignmentType() & SequenceTypes::DNA)) // If so, is it from DNA?
    {
        debug.report(ErrorCode::BacktransAlignIsDNA);
        appearErrors = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::check_and_prepare_coding_sequence() {
    if ((!appearErrors) && (backtransFile != nullptr) &&
        (!backtranslationAlig->prepareCodingSequence(
                splitByStopCodon,
                ignoreStopCodon,
                origAlig
        ))) {
        // Error reporting is made by prepareCodingSequence function.
        appearErrors = true;
        return true;
    }
    return false;
}

inline bool trimAlManager::check_backtranslation_infile_names_corresponde() {
    //NOTE Maybe we don't need to copy the names and lengths to two new arrays 
    // as we could pass the original names and lengths to the check checkCorrespondence function, 
    // which doesn't modify the pointers passed to them
    if ((!appearErrors) && (backtransFile != nullptr)) {
        sequencesNames = new string[backtranslationAlig->getNumSpecies()];
        sequencesLengths = new int[backtranslationAlig->getNumSpecies()];
        backtranslationAlig->getSequences(sequencesNames, sequencesLengths);

        if (!origAlig->checkCorrespondence(
                sequencesNames,
                sequencesLengths,
                backtranslationAlig->getNumSpecies(),
                3
        )) {
            appearErrors = true;
            return true;
        }
    }
    return false;
}

inline void trimAlManager::check_compareset_window_argument() {
    if ((!appearErrors) && (windowSize != -1) && (compareset != -1))
        debug.report(InfoCode::WindowSizeCompareset);
}

inline void trimAlManager::check_output_format() {
    if (oformats.empty() && infile) {
        oformats.emplace_back(ReadWriteMachine.getFileFormatName(infile));
    }
}

int trimAlManager::perform() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("int trimAlManager::perform() ");

    if (appearErrors) return -1;

    if (conservationThreshold == -1)
        conservationThreshold = 0;

    if (vcfs != nullptr) {
        auto XX = ReadWriteMachine.splitAlignmentKeeping(*origAlig);

        char replacement = '-';

        ngs::readVCF(
                /* Dataset          */ XX,
                /* VCF Collection   */ *vcfs,
                /* min Quality      */ 0,
                /* min Coverage     */ 30,
                /* ignore Filters   */ false,
                /* replacement char */ &replacement
        );

        for (auto &i : XX) {
            delete origAlig;
            origAlig = i;

            origAlig->Cleaning->setTrimTerminalGapsFlag(terminalOnly);
            origAlig->setKeepSequencesFlag(keepSeqs);

            set_window_size();

            if (blockSize != -1)
                origAlig->setBlockSize(blockSize);

            if (!create_or_use_similarity_matrix())
                return -2;

            print_statistics();

            clean_alignment();

            postprocess_alignment();

            if ((outfile != nullptr) && (!appearErrors)) {
                std::string outFileString = std::string(outfile);
                if (!ReadWriteMachine.saveAlignment(outFileString, &oformats, origAlig)) {
                    debug.report(ErrorCode::ImpossibleToGenerate, new string[1]{"the output file"});
                    appearErrors = true;
                }
            }
        }
        return 0;
    }

    origAlig->Cleaning->setTrimTerminalGapsFlag(terminalOnly);
    origAlig->setKeepSequencesFlag(keepSeqs);

    set_window_size();

    if (blockSize != -1)
        origAlig->setBlockSize(blockSize);

    if (!create_or_use_similarity_matrix())
        return -2;

    print_statistics();

    clean_alignment();

    if (singleAlig == nullptr) {
        singleAlig = origAlig;
        origAlig = nullptr;
    }

    postprocess_alignment();

    if ((svgOutFile != nullptr) && (!appearErrors))
        if (!origAlig->
                alignmentSummarySVG(*singleAlig, svgOutFile, 0)) {
            debug.report(ErrorCode::ImpossibleToGenerate, new string[1]{"the SVG output file"});
            appearErrors = true;
        }

    if ((htmlOutFile != nullptr) && (!appearErrors))
        if (!origAlig->
                alignmentSummaryHTML(*singleAlig, htmlOutFile)) {
            debug.report(ErrorCode::ImpossibleToGenerate, new string[1]{"the HTML output file"});
            appearErrors = true;
        }

    if ((outfile != nullptr) && (!appearErrors)) {
        std::string outFileString = std::string(outfile);
        if (!ReadWriteMachine.saveAlignment(outFileString, &oformats, singleAlig)) {
            appearErrors = true;
        }

    } else if ((stats >= 0) && (!appearErrors))
        ReadWriteMachine.saveAlignment("", &oformats, singleAlig);

    if ((columnNumbering) && (!appearErrors))
        singleAlig->Statistics->printCorrespondence();

    return 0;
}

inline void trimAlManager::print_statistics() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("inline void trimAlManager::print_statistics() ");
    if (svgStatsOutFile != nullptr) {
        std::string title = infile;
        std::string filename = svgStatsOutFile;
        std::string linename = "";
        std::string color = "";

        utils::streamSVG(nullptr, nullptr, 0, nullptr, nullptr, &title, &filename);

        if (origAlig->Statistics->calculateGapStats()) {
            float acm = 0.0F;
            float x = 0;
            float y = 1.F;
            color = "red";
            linename = "gaps";
            utils::streamSVG(&x, &y, 0, &linename, &color, nullptr, nullptr);

            for (i = 0, acm = 0; i <= origAlig->Statistics->gaps->maxGaps; i++) {
                // If the columns' number with this gaps' number is not equal to zero, we will count the columns. 
                if (origAlig->Statistics->gaps->numColumnsWithGaps[i] != 0) {
                    // Compute and prints the accumulative values for the gaps in the alignment. 
                    acm += origAlig->Statistics->gaps->numColumnsWithGaps[i];
                    x = acm / origAlig->residNumber;
                    y = 1.F - ((i * 1.0F) / origAlig->sequenNumber);
                    utils::streamSVG(&x, &y, 0, &linename, &color, nullptr, nullptr);
                }
            }
        }

        if (origAlig->Statistics->calculateConservationStats()) {
            color = "blue";
            linename = "conservation";
            float x = 0;
            float y = 1.F;
            utils::streamSVG(&x, &y, 0, &linename, &color, nullptr, nullptr);
            float refer, *vectAux;
            int i, num, acm;

            // Allocate memory 
            vectAux = new float[origAlig->residNumber];

            // Select the conservation's value source and copy that vector in a auxiliar vector 
            if (origAlig->Statistics->conservation->MDK_Window != nullptr)
                utils::copyVect(origAlig->Statistics->conservation->MDK_Window, vectAux, origAlig->residNumber);
            else
                utils::copyVect(origAlig->Statistics->conservation->MDK, vectAux, origAlig->residNumber);

            // Sort the auxiliar vector. 
            utils::quicksort(vectAux, 0, origAlig->residNumber - 1);

            // Initializate some values 
            refer = vectAux[origAlig->residNumber - 1];
            acm = 0;
            num = 1;

            // Count the columns with the same conservation's value and compute this information to shows the accunulative
            // statistics in the alignment. 
            for (i = origAlig->residNumber - 2; i >= 0; i--) {
                acm++;

                if (refer != vectAux[i]) {
                    x = ((float) acm / origAlig->residNumber);
                    y = refer;
                    utils::streamSVG(&x, &y, 0, &linename, &color, nullptr, nullptr);
                    refer = vectAux[i];
                    num = 1;
                } else num++;
            }
            acm++;
            x = ((float) acm / origAlig->residNumber);
            y = refer;
            utils::streamSVG(&x, &y, 0, &linename, &color, nullptr, nullptr);

            // Deallocate the reserved memory. 
            delete[] vectAux;
        }

        utils::streamSVG(nullptr, nullptr, 0, nullptr, nullptr, nullptr, nullptr);
    }


    {
        if (sgc) {
            origAlig->Statistics->printStatisticsGapsColumns();
            stats++;
        }

        if (sgt) {
            origAlig->Statistics->printStatisticsGapsTotal();
            stats++;
        }

        if (scc) {
            origAlig->Statistics->printStatisticsConservationColumns();
            stats++;
        }

        if (sct) {
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

        if (compareset != -1) {
            if (sfc)
                statisticsConsistency::printStatisticsFileColumns(*origAlig, origAlig->Statistics->consistency->getValues());
            if (sft)
                statisticsConsistency::printStatisticsFileAcl(*origAlig, origAlig->Statistics->consistency->getValues());
        }
    }
}

inline bool trimAlManager::create_or_use_similarity_matrix() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("inline bool trimAlManager::create_or_use_similarity_matrix() ");
    if ((strict) || (strictplus) || (automated1) || (similarityThreshold != -1.0) || (scc == 1) || (sct == 1)) {
        similMatrix = new similarityMatrix();

        // Simple Matrix
        if (matrixFile != nullptr)
            similMatrix->loadSimMatrix(matrixFile);

            // Alternative Matrix
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

        if (!origAlig->Statistics->setSimilarityMatrix(similMatrix)) {
            debug.report(ErrorCode::ImpossibleToProcessMatrix);
            return false;
        }
    }
    return true;
}

inline void trimAlManager::clean_alignment() {
    // Create a timerLevel that will report times upon its destruction
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

inline void trimAlManager::postprocess_alignment() {

    // only terminal
    if (terminalOnly)
        singleAlig->Cleaning->removeOnlyTerminal();

    // complementary
    if (getComplementary)
        singleAlig->Cleaning->computeComplementaryAlig(true, false);

    // backtranslate
    if (backtransFile != nullptr) {
        tempAlig = backtranslationAlig->getTranslationCDS(singleAlig);
        delete singleAlig;
        singleAlig = tempAlig;
        tempAlig = nullptr;
    }

}

inline void trimAlManager::CleanSequences() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("inline void trimAlManager::CleanSequences() ");
    tempAlig = nullptr;

    if (clusters != -1) {

        tempAlig = origAlig->Cleaning->getClustering(
                origAlig->Cleaning->getCutPointClusters(clusters)
        );
    } else if (maxIdentity != -1) {
        tempAlig = origAlig->Cleaning->getClustering(maxIdentity);
    } else if (delSequences != nullptr) {

        {
            tempAlig = origAlig->Cleaning->removeSequences(
                    delSequences,
                    1,
                    delSequences[0], // Num of sequences in array
                    false
            );
        }
    } else if ((residuesOverlap != -1) && (sequenceOverlap != -1)) {
        tempAlig = origAlig->Cleaning->cleanSpuriousSeq(
                residuesOverlap,
                sequenceOverlap / 100.0F,
                /* getComplementary*/ false
        );
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

        delete singleAlig->Statistics->conservation;
        singleAlig->Statistics->conservation = nullptr;

    } else {
        singleAlig = origAlig;
    }
}

inline void trimAlManager::CleanResiduesAuto() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("inline void trimAlManager::CleanResiduesAuto() ");
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

inline void trimAlManager::CleanResiduesNonAuto() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("inline void trimAlManager::CleanResiduesNonAuto() ");
    // Here we use singleAlig as source alignment as the previous step,
    // CleanSequences initializes singleAlig.
    // singleAlig can be a derived alignment from origAlig or origAlig itself

    if (delColumns != nullptr) {
        for (int i = 0 ; i < delColumns[0]; i++)
        {
            if (delColumns[i] >= singleAlig->getNumAminos()) {
                debug.report(
                        ErrorCode::SelectOnlyAccepts,
                        new string[2]{"-selectcols", "residues"}
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

inline void trimAlManager::set_window_size() {
    // Create a timerLevel that will report times upon its destruction
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

inline void trimAlManager::delete_variables() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("inline void trimAlManager::delete_variables() ");

    if (singleAlig == origAlig)
        singleAlig = nullptr;
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

    delete[] infile;
    infile = nullptr;
    delete[] matrixFile;
    matrixFile = nullptr;

    delete forceFile;
    forceFile = nullptr;
    delete backtransFile;
    backtransFile = nullptr;
    delete backtranslationAlig;
    backtranslationAlig = nullptr;

    delete vcfs;
    vcfs = nullptr;
}

void trimAlManager::menu() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void trimAlManager::menu(void) ");

    // Simple idiom to allow including a text file into a char array at compile time
    // Explanation: https://stackoverflow.com/a/25021520
    std::string menu = {
#include "menu.txt"
    };

    utils::ReplaceStringInPlace(menu, "[iformat]", ReadWriteMachine.getInputFormatsAvailable());
    utils::ReplaceStringInPlace(menu, "[oformat]", ReadWriteMachine.getOutputFormatsAvailable());
    utils::ReplaceStringInPlace(menu, "[version]", VERSION);
    utils::ReplaceStringInPlace(menu, "[revision]", REVISION);
    utils::ReplaceStringInPlace(menu, "[build]", BUILD);
    utils::ReplaceStringInPlace(menu, "[authors]", AUTHORS);

    std::cout << menu;
}

void trimAlManager::examples() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void trimAlManager::examples(void) ");

    // Simple idiom to allow including a text file into a char array at compile time
    // Explanation: https://stackoverflow.com/a/25021520
    std::string examples = {
#include "examples.txt"
    };

    std::cout << examples;
}



