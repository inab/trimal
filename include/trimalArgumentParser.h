//
// Created by bioinfo on 8/06/17.
//

#ifndef TRIMAL_TRIMALARGUMENTPARSER_H
#define TRIMAL_TRIMALARGUMENTPARSER_H

#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include "newAlignment.h"

using namespace std;

class trimalArgumentParser
{

    bool appearErrors = false,
         getComplementary = false, columnNumbering = false,
         nogaps = false, noallgaps = false,
         gappyout = false, strict = false,
         strictplus = false, automated1 = false,
         sgc = false, sgt = false, scc = false, sct = false,
         sfc = false, sft = false, sident = false, selectSeqs = false,
         selectCols = false, shortNames = false, splitByStopCodon = false,
         terminalOnly = false, keepSeqs = false,
         keepHeader = false, ignoreStopCodon = false;

    float conservationThreshold = -1,
          gapThreshold = -1,
          similarityThreshold = -1,
          consistencyThreshold = -1,
          residuesOverlap = -1,
          sequenceOverlap = -1,
          maxIdentity = -1;
    int outformat = -1,
        prevType = -1,
        compareset = -1,
        stats = 0,
        windowSize = -1, gapWindow = -1, similarityWindow = -1,
        consistencyWindow = -1, blockSize = -1, clusters = -1;

    /* Others variables */
    ifstream compare;
    float *compareVect = NULL;
    newAlignment  **compareAlignmentsArray  = NULL;
    string nline, *sequencesNames = NULL;
    sequencesMatrix *seqMatrix = NULL;
    similarityMatrix *similMatrix = NULL;
    newAlignment  *origAlig = NULL, *singleAlig = NULL, *backtranslationAlig = NULL;

    int i = 1, argumentLength, num = 0, maxAminos = 0, numfiles = 0, referFile = 0, *delColumns = NULL, *delSequences = NULL, *sequencesLengths = NULL;
    char c, *forceFile = NULL, *infile = NULL, *backtransFile = NULL, *outfile = NULL, *htmlOutFile = NULL, *matrixFile = NULL,
             **filesToCompare = NULL, line[256];

public:

    void parseArguments(int argc, char *argv[]);

    bool info_arguments(int* argc, char* argv[], int* i);
    bool in_argument(int* argc, char* argv[], int* i);
    bool out_argument(int* argc, char* argv[], int* i);
    bool html_out_argument(int* argc, char* argv[], int* i);
    bool out_format_arguments(int* argc, char* argv[], int* i);
    bool matrix_argument(int* argc, char* argv[], int* i);
    bool compareset_argument(int* argc, char* argv[], int* i);
    bool force_select_argument(int* argc, char* argv[], int* i);
    bool back_trans_argument(int* argc, char* argv[], int* i);
    bool gap_threshold_argument(int* argc, char* argv[], int* i);
    bool similarity_threshold_argument(int* argc, char* argv[], int* i);
    bool consistency_threshold_argument(int* argc, char* argv[], int* i);
    bool conservation_argument(int* argc, char* argv[], int* i);
    bool select_cols_argument(int* argc, char* argv[], int* i);
    bool no_gaps_argument(int* argc, char* argv[], int* i);
    bool no_all_gaps_argument(int* argc, char* argv[], int* i);
    bool keep_seqs_argument(int* argc, char* argv[], int* i);
    bool keep_header_argument(int* argc, char* argv[], int* i);
    bool gappy_out_argument(int* argc, char* argv[], int* i);
    bool strict_argument(int* argc, char* argv[], int* i);
    bool strict_plus_argument(int* argc, char* argv[], int* i);
    bool automated1_argument(int* argc, char* argv[], int* i);
    bool residue_overlap_argument(int* argc, char* argv[], int* i);
    bool sequence_overlap_argument(int* argc, char* argv[], int* i);
    bool seqs_select_argument(int* argc, char* argv[], int* i);
    bool max_identity_argument(int* argc, char* argv[], int* i);
    bool clusters_argument(int* argc, char* argv[], int* i);
    bool terminal_only_argument(int* argc, char* argv[], int* i);
    bool window_argument(int* argc, char* argv[], int* i);
    bool gap_window_argument(int* argc, char* argv[], int* i);
    bool similarity_window_argument(int* argc, char* argv[], int* i);
    bool consistency_window_argument(int* argc, char* argv[], int* i);
    bool block_argument(int* argc, char* argv[], int* i);
    bool stats_arguments(int* argc, char* argv[], int* i);
    bool complementary_argument(int* argc, char* argv[], int* i);
    bool col_numbering_argument(int* argc, char* argv[], int* i);
    bool split_by_stop_codon_argument(int* argc, char* argv[], int* i);
    bool ignore_stop_codon_argument(int* argc, char* argv[], int* i);

    // TODO: Names of this functions should be more informative.
    bool post_process(char* argv[]);
    
    bool check_argument_incompatibilities();

    bool check_force_selection();
    bool check_input_file_with_coding_sequences_argument();
    bool check_file_aligned();
    bool check_similarity_matrix();
    bool check_outputs_coincidence();
    bool check_col_numbering();
    bool check_residue_and_sequence_overlap();
    bool check_html_output_interest();
    bool check_output_file_with_statistics();
    bool check_combinations_among_thresholds();
    bool check_automated_manual_incompatibilities();
    bool check_multiple_files_comparison(char* argv[]);
    bool check_block_size();
    bool check_backtranslations();
    bool check_coding_sequences_type();
    bool check_ignore_or_splitby_stop_codon();
    bool check_and_prepare_coding_sequence();
    bool check_correspondence();
    void check_cw_argument(); // TODO <- HAS TO CHANGE ITS NAME

    int perform();

    void print_statistics();
    bool create_or_use_similarity_matrix();
    void clean_alignment();
    void set_window_size();



    void menu();
    void examples();

};


#endif //TRIMAL_TRIMALARGUMENTPARSER_H
