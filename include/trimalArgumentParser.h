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
#include "ReadWriteMS/ReadWriteMachineState.h"

using namespace std;

class trimAlManager
{

    bool 
    appearErrors        = false,
    getComplementary    = false, 
    columnNumbering     = false,
    nogaps              = false, 
    noallgaps           = false,
    gappyout            = false, 
    strict              = false,
    strictplus          = false, 
    automated1          = false,
    sgc                 = false, 
    sgt                 = false, 
    scc                 = false, 
    sct                 = false,
    sfc                 = false, 
    sft                 = false, 
    sident              = false, 
    selectSeqs          = false,
    selectCols          = false, 
    shortNames          = false, 
    splitByStopCodon    = false,
    terminalOnly        = false, 
    keepSeqs            = false,
    ignoreStopCodon     = false;

    float 
    conservationThreshold   = -1,
    gapThreshold            = -1,
    similarityThreshold     = -1,
    consistencyThreshold    = -1,
    residuesOverlap         = -1,
    sequenceOverlap         = -1,
    maxIdentity             = -1;
    
    int 
    stats                   = 0,
    prevType                = -1,
    compareset              = -1,
    windowSize              = -1, 
    gapWindow               = -1, 
    similarityWindow        = -1,
    consistencyWindow       = -1, 
    blockSize               = -1, 
    clusters                = -1,
    argumentLength          = -1, 
    i                       = 1, 
    num                     = 0, 
    maxAminos               = 0, 
    numfiles                = 0, 
    referFile               = 0, 
    automatedMethodCount    = -1,
    
    *delColumns         = NULL, 
    *delSequences       = NULL, 
    *sequencesLengths   = NULL;
        
    string 
    nline, 
    *sequencesNames     = NULL;
    

    /* Others variables */
    ifstream compare;
    float *compareVect              = NULL;
    sequencesMatrix *seqMatrix      = NULL;
    similarityMatrix *similMatrix   = NULL;
    
    newAlignment  
    *origAlig                   = NULL, 
    *singleAlig                 = NULL, 
    **compareAlignmentsArray    = NULL,
    *backtranslationAlig        = NULL;

    char 
    c, 
    *forceFile          = NULL, 
    *infile             = NULL, 
    *backtransFile      = NULL, 
    *outfile            = NULL, 
    *htmlOutFile        = NULL, 
    *matrixFile         = NULL,
             
    **filesToCompare    = NULL, 
    line[256];
             
    std::vector<std::string> oformats;

public:
///\addtogroup TrimalParsingArguments
///@{ 
    ///\addtogroup ParsingArgumentsFunctions 
    // read arguments.
    ///@{ 
    void parseArguments(int argc, char *argv[]);
private:
        void info_arguments(int* argc, char* argv[], int* i);
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
        bool conservation_threshold_argument(int* argc, char* argv[], int* i);
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

    ///@}
public:
    bool process_arguments(char* argv[]);
    
private:
        ///\addtogroup CheckIncompatibilities 
    // Check incompatibilities between arguments
        ///@{ 
        bool check_arguments_incompatibilities();
            bool check_inFile_incompatibilities();
            bool check_select_cols_and_seqs_incompatibilities();
            bool check_thresholds_incompatibilities();
            bool check_automated_methods_incompatibilities();
            bool check_max_identity_incompatibilities();
            bool check_clusters_incompatibilities();
            bool check_windows_incompatibilities();
            bool check_stats_incompatibilities();
            bool check_codon_behaviour_incompatibility();
            bool check_combinations_among_thresholds();
        ///@}

        ///\addtogroup CheckNeeds 
    // Check cross arguments needs.
        ///@{ 
        bool check_arguments_needs(char* argv[]);

            bool check_force_selection();
            bool check_input_file_with_coding_sequences_argument();
            bool check_file_aligned();
            bool check_similarity_matrix();
            bool check_outputs_coincidence();
            bool check_col_numbering();
            bool check_residue_and_sequence_overlap();
            bool check_html_output_interest();
            bool check_output_file_with_statistics();
            bool check_automated_manual_incompatibilities();
            bool check_multiple_files_comparison(char* argv[]);
            bool check_block_size();
            bool check_backtranslations();
            bool check_coding_sequences_type();
            bool check_and_prepare_coding_sequence();
            bool check_backtranslation_infile_names_corresponde();
            void check_cw_argument(); // TODO <- HAS TO CHANGE ITS NAME
            void check_output_format();
        ///@}
            
    
    ///\addtogroup PerformAlgorithm 
// Arguments are valid, perform.
    ///@{ 
public:
    int perform();
private:
        void print_statistics();
        bool create_or_use_similarity_matrix();
        void clean_alignment();
        void set_window_size();
        
    ///@}
    ///\addtogroup OtherMethods
    ///@{ 
    void delete_variables();
    
    // NON COMPLEX OPTIONS
    void menu();
    void examples();
    ///@}

private:
    ReadWriteMS ReadWriteMachine;
};
///@}


#endif //TRIMAL_TRIMALARGUMENTPARSER_H
