//
// Created by bioinfo on 8/06/17.
//

#ifndef TRIMAL_TRIMALARGUMENTPARSER_H
#define TRIMAL_TRIMALARGUMENTPARSER_H

#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <iosfwd>
#include "newAlignment.h"
#include "ReadWriteMS/ReadWriteMachineState.h"

using namespace std;

class trimAlManager
{
public:
    std::vector<std::string> * vcfs = nullptr;
    
    
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
    ssc                 = false, 
    sst                 = false,
    sfc                 = false, 
    sft                 = false, 
    sident              = false, 
    soverlap            = false, 
    selectSeqs          = false,
    selectCols          = false, 
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
    i                       = 1,
    stats                   = 0,
//    num                     = 0,
    compareset              = -1,
    windowSize              = -1, 
    gapWindow               = -1, 
    similarityWindow        = -1,
    consistencyWindow       = -1, 
    blockSize               = -1, 
    clusters                = -1
    ,
    automatedMethodCount    = -1,
    alternative_matrix      = -1,
    
    *delColumns         = nullptr,
    *delSequences       = nullptr,
    *sequencesLengths   = nullptr;
    
    size_t
            argumentLength = size_t(-1);

    string 
    *sequencesNames     = nullptr;
    

    /* Others variables */
    ifstream compare;
//    float *values              = nullptr;
    similarityMatrix *similMatrix   = nullptr;
    
    newAlignment  
    *origAlig                   = nullptr,
    *singleAlig                 = nullptr,
    *tempAlig                   = nullptr,
    **compareAlignmentsArray    = nullptr,
    *backtranslationAlig        = nullptr;

    char 
    c, 
    *forceFile          = nullptr,
    
    *infile             = nullptr,
    
    *backtransFile      = nullptr,
    
    *outfile            = nullptr,
    *htmlOutFile        = nullptr,
    *svgOutFile         = nullptr,
    *svgStatsOutFile    = nullptr,
    
    *matrixFile         = nullptr,
             
    **filesToCompare    = nullptr,
    line[256];
             
    std::vector<std::string> oformats;

public:
    trimAlManager();
    ~trimAlManager();
    void parseArguments(int argc, char *argv[]);
private:
        void verbosity_argument(const int* argc, char* argv[]);
        void help_arguments(const int *argc, char **argv, int *i);
        bool in_argument(const int* argc, char* argv[], int* i);
        bool vcf_argument(const int* argc, char* argv[], int* i);
        bool out_argument(const int* argc, char* argv[], int* i);
        bool html_out_argument(const int* argc, char* argv[], int* i);
        bool timetracker_out_argument(const int* argc, char* argv[], int* i);
        bool svg_out_argument(const int* argc, char* argv[], int* i);
        bool svg_stats_argument(const int* argc, char* argv[], int* i);
        bool out_format_arguments(const int* argc, char* argv[], int* i);
        bool matrix_argument(const int* argc, char* argv[], int* i);
        bool compareset_argument(const int* argc, char* argv[], int* i);
        bool force_select_argument(const int* argc, char* argv[], int* i);
        bool back_trans_argument(const int* argc, char* argv[], int* i);
        bool gap_threshold_argument(const int* argc, char* argv[], int* i);
        bool similarity_threshold_argument(const int* argc, char* argv[], int* i);
        bool consistency_threshold_argument(const int* argc, char* argv[], int* i);
        bool conservation_threshold_argument(const int* argc, char* argv[], int* i);
        bool select_cols_argument(const int* argc, char* argv[], int* i);
        bool no_gaps_argument(const int* argc, char* argv[], int* i);
        bool no_all_gaps_argument(const int* argc, char* argv[], int* i);
        bool keep_seqs_argument(const int* argc, char* argv[], int* i);
        bool keep_header_argument(const int* argc, char* argv[], int* i);
        bool gappy_out_argument(const int* argc, char* argv[], int* i);
        bool strict_argument(const int* argc, char* argv[], int* i);
        bool strict_plus_argument(const int* argc, char* argv[], int* i);
        bool automated1_argument(const int* argc, char* argv[], int* i);
        bool residue_overlap_argument(const int* argc, char* argv[], int* i);
        bool sequence_overlap_argument(const int* argc, char* argv[], int* i);
        bool seqs_select_argument(const int* argc, char* argv[], int* i);
        bool max_identity_argument(const int* argc, char* argv[], int* i);
        bool clusters_argument(const int* argc, char* argv[], int* i);
        bool terminal_only_argument(const int* argc, char* argv[], int* i);
        bool window_argument(const int* argc, char* argv[], int* i);
        bool gap_window_argument(const int* argc, char* argv[], int* i);
        bool similarity_window_argument(const int* argc, char* argv[], int* i);
        bool consistency_window_argument(const int* argc, char* argv[], int* i);
        bool block_argument(const int* argc, char* argv[], int* i);
        bool stats_arguments(const int* argc, char* argv[], int* i);
        bool complementary_argument(const int* argc, char* argv[], int* i);
        bool col_numbering_argument(const int* argc, char* argv[], int* i);
        bool split_by_stop_codon_argument(const int* argc, char* argv[], int* i);
        bool ignore_stop_codon_argument(const int* argc, char* argv[], int* i);

public:
    bool processArguments(char* argv[]);
    
private:

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
            bool check_combinations_among_thresholds_incompatibility();
            bool check_automated_manual_incompatibilities();

        bool check_arguments_needs(char* argv[]);

            bool check_force_selection();
            bool check_input_file_with_coding_sequences_argument();
            bool check_file_aligned();
            bool check_similarity_matrix();
            bool check_outputs_coincidence();
            bool check_col_numbering();
            bool check_residue_and_sequence_overlap();
            bool check_output_relevance();
            bool check_output_file_with_statistics();
            bool check_multiple_files_comparison(char* argv[]);
            bool check_block_size();
            bool check_backtranslations();
            bool check_coding_sequences_type();
            bool check_and_prepare_coding_sequence();
            bool check_backtranslation_infile_names_corresponde();
            void check_compareset_window_argument(); // TODO <- HAS TO CHANGE ITS NAME
            void check_output_format();

public:
    int perform();
    void delete_variables();
private:
    void output_reports();
    void save_alignment();
    int perform_VCF();
    void svg_stats_out();
    void print_statistics();
    bool create_or_use_similarity_matrix();
    void clean_alignment();
    void postprocess_alignment();

        void set_window_size();

    // NON COMPLEX OPTIONS
    void menu();
    void examples();


    void CleanSequences();
    void CleanResiduesAuto();
    void CleanResiduesNonAuto();


private:
    ReadWriteMS ReadWriteMachine;
};



#endif //TRIMAL_TRIMALARGUMENTPARSER_H
