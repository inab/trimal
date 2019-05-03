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

#ifndef TRIMAL_TRIMALARGUMENTPARSER_H
#define TRIMAL_TRIMALARGUMENTPARSER_H


#include "FormatHandling/FormatManager.h"
#include "Alignment/Alignment.h"

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <iosfwd>
#include <string>
#include <functional>

//#include "../tests/catch.hpp";

namespace Catch
{
    struct SectionInfo;
}


// Forward declatarion
class similarityMatrix;
namespace statistics
{
    class Consistency;
    class similarityMatrix;
}

class trimAlManager
{
public:

    friend struct Catch::SectionInfo;

    std::vector<std::string> * vcfs = nullptr;

    bool 
        appearErrors        = false,

        getComplementary    = false, 
        getComplementarySeq = false,
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
        ignoreStopCodon     = false,
        ignoreFilter        = false,
        removeDuplicates    = false;

    float 
        conservationThreshold   = -1,
        gapThreshold            = -1,
        similarityThreshold     = -1,
        consistencyThreshold    = -1,
        residuesOverlap         = -1,
        sequenceOverlap         = -1,
        maxIdentity             = -1,
        minCoverage             = -1,
        minQuality              = -1;
    
    int
        i                       = 1,
        stats                   = 0,
        windowSize              = -1, 
        gapWindow               = -1, 
        similarityWindow        = -1,
        consistencyWindow       = -1, 
        blockSize               = -1, 
        clusters                = -1
        ,
        automatedMethodCount    = -1,
        alternative_matrix      = -1,
        gapAbsoluteThreshold    = -1,

        *delColumns         = nullptr,
        *delSequences       = nullptr,
        *sequencesLengths   = nullptr;
    
    size_t
            argumentLength = size_t(-1);

    std::string 
        *sequencesNames     = nullptr;
    

    /* Others variables */
    std::ifstream compare;
    statistics::similarityMatrix *similMatrix   = nullptr;
    
    Alignment
        *origAlig                   = nullptr,
        *singleAlig                 = nullptr,
        *tempAlig                   = nullptr,
        *backtranslationAlig        = nullptr,
        **compareAlignmentsArray    = nullptr;

    char 
        c, 
        *forceFile          = nullptr,
        
        *infile             = nullptr,
        
        *backtransFile      = nullptr,
        
        *outfile            = nullptr,
        *htmlOutFile        = nullptr,
        *svgOutFile         = nullptr,
        *svgStatsOutFile    = nullptr,

        *compareset         = nullptr,

        
        *matrixFile         = nullptr,
                
        **filesToCompare    = nullptr,
        line[256];
             
    std::vector<std::string> oformats;

    /// Consistency Manager. 
    ///     We have to save the reference, 
    ///     while it is still not part of an alignment
    statistics::Consistency *CS = nullptr;

public:
    trimAlManager();
    ~trimAlManager();
    int parseArguments(int argc, char **argv);

    enum argumentReport {
        NotRecognized   = 0,
        Recognized      = 1,
        Errored         = 2,
        Final           = 3
    };

private: // Parse Arguments Methods
        void verbosity_argument             (const int* argc, char* argv[]);
        
        argumentReport help_arguments(const int *argc, char **argv, int *currentArg);

        argumentReport in_argument                    (const int* argc, char* argv[], int* currentArg);
        argumentReport vcf_argument                   (const int* argc, char* argv[], int* currentArg);
        argumentReport out_argument                   (const int* argc, char* argv[], int* currentArg);
        argumentReport html_out_argument              (const int* argc, char* argv[], int* currentArg);
        argumentReport timetracker_out_argument       (const int* argc, char* argv[], int* currentArg);
        argumentReport svg_out_argument               (const int* argc, char* argv[], int* currentArg);
        argumentReport svg_stats_argument             (const int* argc, char* argv[], int* currentArg);
        argumentReport out_format_arguments           (const int* argc, char* argv[], int* currentArg);
        argumentReport matrix_argument                (const int* argc, char* argv[], int* currentArg);
        argumentReport compareset_argument            (const int* argc, char* argv[], int* currentArg);
        argumentReport force_select_argument          (const int* argc, char* argv[], int* currentArg);
        argumentReport back_trans_argument            (const int* argc, char* argv[], int* currentArg);
        argumentReport gap_threshold_argument         (const int* argc, char* argv[], int* currentArg);
        argumentReport similarity_threshold_argument  (const int* argc, char* argv[], int* currentArg);
        argumentReport consistency_threshold_argument (const int* argc, char* argv[], int* currentArg);
        argumentReport conservation_threshold_argument(const int* argc, char* argv[], int* currentArg);
        argumentReport select_cols_argument           (const int* argc, char* argv[], int* currentArg);
        argumentReport no_gaps_argument               (const int* argc, char* argv[], int* currentArg);
        argumentReport no_all_gaps_argument           (const int* argc, char* argv[], int* currentArg);
        argumentReport keep_seqs_argument             (const int* argc, char* argv[], int* currentArg);
        argumentReport keep_header_argument           (const int* argc, char* argv[], int* currentArg);
        argumentReport gappy_out_argument             (const int* argc, char* argv[], int* currentArg);
        argumentReport strict_argument                (const int* argc, char* argv[], int* currentArg);
        argumentReport strict_plus_argument           (const int* argc, char* argv[], int* currentArg);
        argumentReport automated1_argument            (const int* argc, char* argv[], int* currentArg);
        argumentReport residue_overlap_argument       (const int* argc, char* argv[], int* currentArg);
        argumentReport sequence_overlap_argument      (const int* argc, char* argv[], int* currentArg);
        argumentReport seqs_select_argument           (const int* argc, char* argv[], int* currentArg);
        argumentReport max_identity_argument          (const int* argc, char* argv[], int* currentArg);
        argumentReport clusters_argument              (const int* argc, char* argv[], int* currentArg);
        argumentReport terminal_only_argument         (const int* argc, char* argv[], int* currentArg);
        argumentReport window_argument                (const int* argc, char* argv[], int* currentArg);
        argumentReport gap_window_argument            (const int* argc, char* argv[], int* currentArg);
        argumentReport similarity_window_argument     (const int* argc, char* argv[], int* currentArg);
        argumentReport consistency_window_argument    (const int* argc, char* argv[], int* currentArg);
        argumentReport block_argument                 (const int* argc, char* argv[], int* currentArg);
        argumentReport stats_arguments                (const int* argc, char* argv[], int* currentArg);
        argumentReport complementary_argument         (const int* argc, char* argv[], int* currentArg);
        argumentReport col_numbering_argument         (const int* argc, char* argv[], int* currentArg);
        argumentReport split_by_stop_codon_argument   (const int* argc, char* argv[], int* currentArg);
        argumentReport ignore_stop_codon_argument     (const int* argc, char* argv[], int* currentArg);
        argumentReport ignore_filter_argument         (const int* argc, char* argv[], int* currentArg);
        argumentReport min_quality_argument           (const int* argc, char* argv[], int* currentArg);
        argumentReport min_coverage_argument          (const int* argc, char* argv[], int* currentArg);
        argumentReport remove_duplicates_argument     (const int* argc, char* argv[], int* currentArg);

public:
    bool processArguments(char* argv[]);
private: // Process Arguments Methods
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
        bool check_vcf_incompatibility();

    bool check_arguments_needs(char* argv[]);

        bool check_absolute_gap_theshold();

        bool check_force_selection();
        bool check_input_file_with_coding_sequences_argument();
        bool check_file_aligned();
        bool check_similarity_matrix();
        bool check_outputs_coincidence();
        bool check_col_numbering();
        bool check_residue_and_sequence_overlap();
        bool check_output_relevance();
        bool check_output_file_with_statistics();
        bool check_block_size();
        bool check_backtranslations();
        bool check_coding_sequences_type();
        bool check_and_prepare_coding_sequence();
        bool check_backtranslation_infile_names_correspondence();
        void check_compareset_window_argument();
        void check_output_format();
        void check_thresholds_dependencies();

public:
    int perform();
    void delete_variables();

private: // General, private, methods
    int innerPerform();
    bool performCompareset();

    void output_reports();
    void save_alignment();

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
    FormatHandling::FormatManager formatManager;
public:
    FormatHandling::FormatManager& getFormatManager();


public:

    // Used for testing to allow exposing private methods and variables
    // Each implementation of access_by contains different methods and variables.
    // This structs can access all private methods of trimAlManager
    //
    // Original implementation: https://softwareengineering.stackexchange.com/a/257721
    template <uint T> struct access_by;
    template <uint T> friend struct access_by;

// Defines to simplify the creation and usage of access_by structs

// Takes one argument, a const char * that will be used as UUID
// This const char * is converted to a hash using a hasher (hasher.h) in compile time
// Then, the struct must be defined. In the method definitions of the struct,
//      we can access the private methods and variables of trimAlManager
#define trimAlManagerPrivateExposerCreator(TOKEN) template <> struct trimAlManager::access_by< TOKEN""_hash>

// Allows the easy access to access_by structs.
// This structs can access all private methods of trimAlManager
// It allows to create an instance of the struct: trimAlManagerPrivateExposerCaller(X)().call()
// It also allows to access the struct statically: trimAlManagerPrivateExposerCaller(X)::call()
#define trimAlManagerPrivateExposerCaller(TOKEN) trimAlManager::access_by<TOKEN""_hash>

};




#endif //TRIMAL_TRIMALARGUMENTPARSER_H
