//
// Created by bioinfo on 8/06/17.
//

#include "trimalArgumentParser.h"
#include <string.h>
void menu();
void examples();
void trimalArgumentParser::parseArguments(int argc, char *argv[]) {

    origAlig = new newAlignment();

    for(int i = 0; i < argc; i++ )
    {
        if (appearErrors) break;

        if (info_arguments(&argc, argv, &i)) continue;
        if (in_argument(&argc, argv, &i)) continue;
        if (out_argument(&argc, argv, &i)) continue;
        if (html_out_argument(&argc, argv, &i)) continue;
        if (out_format_arguments(&argc, argv, &i)) continue;
        if (matrix_argument(&argc, argv, &i)) continue;

        if (compareset_argument(&argc, argv, &i)) continue;
        if (force_select_argument(&argc, argv, &i)) continue;
        if (back_trans_argument(&argc, argv, &i)) continue;

        if (gap_threshold_argument(&argc, argv, &i)) continue;
        if (sim_threshold_argument(&argc, argv, &i)) continue;
        if (cont_threshold_argument(&argc, argv, &i)) continue;

        if (cons_argument(&argc, argv, &i)) continue;
        if (select_cols_argument(&argc, argv, &i)) continue;

        if (no_gaps_argument(&argc, argv, &i)) continue;
        if (no_all_gaps_argument(&argc, argv, &i)) continue;

        if (keep_seqs_argument(&argc, argv, &i)) continue;
        if (keep_header_argument(&argc, argv, &i)) continue;

        if (gappy_out_argument(&argc, argv, &i)) continue;
        if (strict_argument(&argc, argv, &i)) continue;
        if (strict_plus_argument(&argc, argv, &i)) continue;
        if (automated1_argument(&argc, argv, &i)) continue;

        if (col_overlap_argument(&argc, argv, &i)) continue;
        if (seq_overlap_argument(&argc, argv, &i)) continue;

        if (seqs_select_argument(&argc, argv, &i)) continue;

        if (max_identity_argument(&argc, argv, &i)) continue;
        if (clusters_argument(&argc, argv, &i)) continue;

        if (terminal_only_argument(&argc, argv, &i)) continue;

        if (window_argument(&argc, argv, &i)) continue;
        if (gap_window_argument(&argc, argv, &i)) continue;
        if (sim_window_argument(&argc, argv, &i)) continue;
        if (con_window_argument(&argc, argv, &i)) continue;

        if (block_argument(&argc, argv, &i)) continue;
        if (stats_arguments(&argc, argv, &i)) continue;

        if (complementary_argument(&argc, argv, &i)) continue;

        if (col_numbering_argument(&argc, argv, &i)) continue;
        if (split_by_stop_codon_argument(&argc, argv, &i)) continue;
        if (ignore_stop_codon_argument(&argc, argv, &i)) continue;
    }
}

bool trimalArgumentParser::info_arguments(int *argc, char *argv[], int *i) {

    if (!strcmp(argv[*i], "-h"))
    {
        cout << "Print help" << endl;
        return true;
    }

    if(!strcmp(argv[*i], "--version")) {

        cout << "Print version" << endl;
        return true;
    }

    return false;
}

bool trimalArgumentParser::in_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-in") && (i+1 != argc) && (infile == NULL)) {

        if((sfc) || (sft) || (comThreshold != -1)) {
            cerr << endl << "ERROR: Not allowed in combination of file comparision." << endl << endl;
            appearErrors = true;
            i++;
        }

        else if((compareset == -1) || (forceFile != NULL)) {
            lng = strlen(argv[++*i]);
            infile = new char[lng + 1];
            strcpy(infile, argv[*i]);

            if(!origAlig -> ReadWrite -> loadAlignment(infile)) {
                cerr << endl << "ERROR: newAlignment  not loaded: \"" << infile << "\" Check the file's content." << endl << endl;
                appearErrors = true;
            }
        }

        else {
            if(compareset != -1)
                cerr << endl << "ERROR: Option \"" << argv[*i] << "\" not valid. A reference file exists with alignments to compare." << endl << endl;
            if(forceFile != NULL)
                cerr << endl << "ERROR: Option \"" << argv[*i] << "\" not valid. A newAlignment  file has been setting up to be compare with a set of alignmets." << endl << endl;
            appearErrors = true;
            i++;
        }
    }
}

bool trimalArgumentParser::out_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-out")) && (i+1 != argc) && (outfile == NULL)) {
        lng = strlen(argv[++*i]);
        outfile = new char[lng + 1];
        strcpy(outfile, argv[*i]);
    }
}

bool trimalArgumentParser::html_out_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-htmlout")) && (i+1 != argc) && (outhtml == NULL)) {
        lng = strlen(argv[++*i]);
        outhtml = new char[lng + 1];
        strcpy(outhtml, argv[*i]);
    }
}

bool trimalArgumentParser::out_format_arguments(int *argc, char *argv[], int *i) {

    if (outformat == -1)
    {
        /* Option -clustal -------------------------------------------------------------------------------------- */
        if(!strcmp(argv[*i], "-clustal") ) {
            outformat = 1;
        }

            /* Option -fasta -------------------------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-fasta") )
            outformat = 8;

            /* Option -fasta-m10 -------------------------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-fasta_m10") ) {
            outformat = 8; shortNames = true;
        }

            /* Option -nbrf ------------------------------------------------------------------------------------ */
        else if(!strcmp(argv[*i], "-nbrf") )
            outformat = 3;

            /* Option -nexus ------------------------------------------------------------------------------------ */
        else if(!strcmp(argv[*i], "-nexus") )
            outformat = 17;

            /* Option -mega ------------------------------------------------------------------------------------ */
        else if(!strcmp(argv[*i], "-mega") )
            outformat = 21;

            /* Option -phylip3.2 --------------------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-phylip3.2") )
            outformat = 11;

            /* Option -phylip3.2-m10 ----------------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-phylip3.2_m10") ) {
            outformat = 11; shortNames = true;
        }

            /* Option -phylip --------------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-phylip") )
            outformat = 12;

            /* Option -phylip-m10 ----------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-phylip_m10") ) {
            outformat = 12; shortNames = true;
        }

            /* Option -phylip_paml ---------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-phylip_paml") )
            outformat = 13;

            /* Option -phylip_paml-m10 ------------------------------------------------------------------ */
        else if(!strcmp(argv[*i], "-phylip_paml_m10") ) {
            outformat = 13; shortNames = true;
        }
    }

}

bool trimalArgumentParser::matrix_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-matrix") && (i+1 != argc) && (matrix == NULL)) {
        lng = strlen(argv[++*i]);
        matrix = new char[lng + 1];
        strcpy(matrix, argv[*i]);
    }
}

bool trimalArgumentParser::compareset_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-compareset") && (i+1 != argc) && (compareset == -1)) {

        if(infile == NULL) {
            compare.open(argv[++*i], ifstream::in);
            if(!compare) {
                cerr << endl << "ERROR: Check the reference file with the alignments to compare." << endl << endl;
                appearErrors = true;
            }

            while(compare.getline(line, 256)) numfiles++;
            compare.close();

            compareset = *i;
        }

        else {
            cerr << endl << "ERROR: Option \"" << argv[*i] << "\" not valid. A single newAlignment  file has been set by the user." << endl << endl;
            appearErrors = true;
            i++;
        }
    }
}

bool trimalArgumentParser::force_select_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-forceselect") && (i+1 != argc) && (forceFile == NULL)) {

        if(infile == NULL) {
            lng = strlen(argv[++*i]);
            forceFile = new char[lng + 1];
            strcpy(forceFile, argv[*i]);
            if(!origAlig -> ReadWrite -> loadAlignment(forceFile)) {
                cerr << endl << "ERROR: newAlignment  not loaded: \"" << forceFile << "\" Check the file's content." << endl << endl;
                appearErrors = true;
            }
        }

        else {
            cerr << endl << "ERROR: Option \"" << argv[*i] << "\" not valid. A single newAlignment  file has been setting it up" << endl << endl;
            appearErrors = true;
            i++;
        }
    }
}

bool trimalArgumentParser::back_trans_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-backtrans") && (i+1 != argc) && (backtransFile == NULL)) {

        lng = strlen(argv[++*i]);
        backtransFile = new char[lng + 1];
        strcpy(backtransFile, argv[*i]);

        backtranslation = new newAlignment;
        if(!backtranslation -> ReadWrite -> loadAlignment(backtransFile)) {
            cerr << endl << "ERROR: newAlignment  not loaded: \"" << backtransFile << "\" Check the file's content." << endl << endl;
            appearErrors = true;
        }
    }
}

bool trimalArgumentParser::gap_threshold_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-gapthreshold") || !strcmp(argv[*i], "-gt")) && (i+1 != argc) && (gapThreshold == -1)) {

        if((selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else {
            if(utils::isNumber(argv[++*i])) {
                gapThreshold = 1 - atof(argv[*i]);
                if((gapThreshold < 0) || (gapThreshold > 1)) {
                    cerr << endl << "ERROR: The gap threshold value should be between 0 and 1." << endl << endl;
                    appearErrors = true;
                }
            }
            else {
                cerr << endl << "ERROR: The gap threshold value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
    }
}

bool trimalArgumentParser::sim_threshold_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-simthreshold") || !strcmp(argv[*i], "-st")) && (i+1 != argc) && (simThreshold == -1)) {

        if((selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else {
            if(utils::isNumber(argv[++*i])) {
                simThreshold = atof(argv[*i]);
                if((simThreshold < 0) || (simThreshold > 1)) {
                    cerr << endl << "ERROR: The similarity threshold value should be between 0 and 1." << endl << endl;
                    appearErrors = true;
                }
            }
            else {
                cerr << endl << "ERROR: The similarity threshold value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
    }
}

bool trimalArgumentParser::cont_threshold_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-conthreshold") || !strcmp(argv[*i], "-ct")) && (i+1 != argc) && (comThreshold == -1)) {

        if((selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

            //~ else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
            //~ cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            //~ appearErrors = true;
            //~ }

        else if(infile != NULL) {
            cerr << endl << "ERROR: Not allowed in combination with -in option." << endl << endl;
            appearErrors = true;

        }

        else {
            if(utils::isNumber(argv[++*i])) {
                comThreshold = atof(argv[*i]);
                if((comThreshold < 0) || (comThreshold > 1)) {
                    cerr << endl << "ERROR: The consistency threshold value should be between 0 and 1." << endl << endl;
                    appearErrors = true;
                }
            }
            else {
                cerr << endl << "ERROR: The consistency threshold value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
    }
}

bool trimalArgumentParser::cons_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-cons")) && (i+1 != argc) && (conserve == -1)) {

        if((selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else if(blockSize != -1) {
            cerr << endl << "ERROR: Not allowed in combination of column block size value." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus)  || (automated1)) {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else {
            if(utils::isNumber(argv[++*i])) {
                conserve = atof(argv[*i]);
                if((conserve < 0) || (conserve > 100)) {
                    cerr << endl << "ERROR: The minimal positions value should be between 0 and 100." << endl << endl;
                    appearErrors = true;
                }
            }
            else {
                cerr << endl << "ERROR: The minimal positions value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
    }
}

bool trimalArgumentParser::select_cols_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-selectcols")) && (selectCols == false) &&
            ((i+3) < argc) && (!strcmp(argv[++*i], "{")) && (!strcmp(argv[*i+2], "}"))) {

        if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed." << endl << endl;
            appearErrors = true;
        }

        else if(blockSize != -1) {
            cerr << endl << "ERROR: Not allowed in combination of column block size value." << endl << endl;
            appearErrors = true;
        }

        else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || (comThreshold != -1)) {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods." << endl << endl;
            appearErrors = true;
        }

        else if((windowSize != -1) || (gapWindow != -1)|| (simWindow != -1)) {
            cerr << endl << "ERROR: It's imposible to use this windows size in combination of selection method." << endl << endl;
            appearErrors = true;
        }

        else if((delColumns = utils::readNumbers(argv[++*i])) == NULL) {
            cerr << endl << "ERROR: Impossible to parser the sequences number" << endl << endl;
            appearErrors = true;
        }

        else selectCols = true;
        i++;
    }
}

bool trimalArgumentParser::no_gaps_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-nogaps") && (!nogaps)) {

        if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
            cerr << endl << "ERROR: Not allowed in combination of window values." << endl << endl;
            appearErrors = true;
        }

        else if(blockSize != -1) {
            cerr << endl << "ERROR: Not allowed in combination of column block size value." << endl << endl;
            appearErrors = true;
        }

        else if((noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
            cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
            appearErrors = true;
        }

            //~ else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
            //~ (comThreshold != -1) || (selectCols) || (selectSeqs)) {
            //~ cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            //~ appearErrors = true;
            //~ }
        else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
                (selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
            nogaps = true;
    }
}

bool trimalArgumentParser::no_all_gaps_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-noallgaps") && (!noallgaps)) {

        if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
            cerr << endl << "ERROR: Not allowed in combination of window values." << endl << endl;
            appearErrors = true;
        }

        else if(blockSize != -1) {
            cerr << endl << "ERROR: Not allowed in combination of column block size value." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
            cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
            appearErrors = true;
        }

            //~ else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
            //~ (comThreshold != -1) || (selectCols) || (selectSeqs)) {
            //~ cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            //~ appearErrors = true;
            //~ }
        else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
                (selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
            noallgaps = true;
    }
}

bool trimalArgumentParser::keep_seqs_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-keepseqs") && (!keepSeqs)) {
        keepSeqs = true;
    }
}

bool trimalArgumentParser::keep_header_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-keepheader") && (!keepHeader)) {
        keepHeader = true;
    }
}

bool trimalArgumentParser::gappy_out_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-gappyout") && (!strict)) {

        if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
            cerr << endl << "ERROR: Not allowed in combination of window values." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (strict) || (strictplus) || (automated1)) {
            cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
            appearErrors = true;
        }

            //~ else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
            //~ (comThreshold != -1) || (selectCols) || (selectSeqs)) {
            //~ cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            //~ appearErrors = true;
            //~ }
        else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
                (selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
            gappyout = true;
    }
}

bool trimalArgumentParser::strict_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-strict") && (!strict)) {

        if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
            cerr << endl << "ERROR: Not allowed in combination of window values." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strictplus) || (automated1)) {
            cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
            appearErrors = true;
        }

        else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
                (selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
            strict = true;
    }
}

bool trimalArgumentParser::strict_plus_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-strictplus")) && (!strictplus)) {

        if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
            cerr << endl << "ERROR: Not allowed in combination with this window value." << endl << endl;
            appearErrors = true;
        }

            //~ else if(blockSize != -1) {
            //~ cerr << endl << "ERROR: Not allowed in combination of column block size value." << endl << endl;
            //~ appearErrors = true;
            //~ }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (automated1)) {
            cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
            appearErrors = true;
        }

            //~ else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
            //~ (comThreshold != -1) || (selectCols) || (selectSeqs)) {
            //~ cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            //~ appearErrors = true;
            //~ }
        else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
                (selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
            strictplus = true;
    }
}

bool trimalArgumentParser::automated1_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-automated1")) && (!automated1)) {

        if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
            cerr << endl << "ERROR: Not allowed in combination with this window value." << endl << endl;
            appearErrors = true;
        }

        else if(blockSize != -1) {
            cerr << endl << "ERROR: Not allowed in combination of column block size value." << endl << endl;
            appearErrors = true;
        }

        else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
                (comThreshold != -1) || (selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
            appearErrors = true;
        }

        else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
                (selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
            automated1 = true;
    }
}

bool trimalArgumentParser::col_overlap_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-resoverlap")) && (i+1 != argc) && (resOverlap == -1)) {

        if((selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Not allowed in combination of methods such as manual selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else {
            if(utils::isNumber(argv[++*i])) {
                resOverlap = atof(argv[*i]);
                if((resOverlap < 0) || (resOverlap > 1)) {
                    cerr << endl << "ERROR: The residue overlap value should be between 0 and 1." << endl << endl;
                    appearErrors = true;
                }
            }
            else {
                cerr << endl << "ERROR: The residue overlap value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
    }
 }

bool trimalArgumentParser::seq_overlap_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-seqoverlap")) && (i+1 != argc) && (seqOverlap == -1)) {

        if((selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Not allowed in combination of methods such as manual selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else {
            if(utils::isNumber(argv[++*i])) {
                seqOverlap = atof(argv[*i]);
                if((seqOverlap < 0) || (seqOverlap > 100)) {
                    cerr << endl << "ERROR: The sequences overlap value should be between 0 and 100." << endl << endl;
                    appearErrors = true;
                }
            }
            else {
                cerr << endl << "ERROR: The minimal positions value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
    }
}

bool trimalArgumentParser::seqs_select_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-selectseqs")) && (selectSeqs == false) && ((i+3) < argc) && (!strcmp(argv[++*i], "{")) && (!strcmp(argv[*i+2], "}"))) {

        if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed." << endl << endl;
            appearErrors = true;
        }

        else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || (comThreshold != -1)) {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods." << endl << endl;
            appearErrors = true;
        }

        else if((windowSize != -1) || (gapWindow != -1)|| (simWindow != -1)) {
            cerr << endl << "ERROR: It's imposible to use this windows size in combination of selection method." << endl << endl;
            appearErrors = true;
        }

        else if((clusters != -1) || (maxIdentity != -1)) {
            cerr << endl << "ERROR: Only one method to chose sequences can be applied." << endl << endl;
            appearErrors = true;
        }

        else if((delSequences = utils::readNumbers(argv[++*i])) == NULL) {
            cerr << endl << "ERROR: Impossible to parser the sequences number" << endl << endl;
            appearErrors = true;
        }

        else selectSeqs = true;
        i++;
    }
}

bool trimalArgumentParser::max_identity_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-maxidentity")) && (i+1 != argc) && (maxIdentity == -1)) {

        if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
           (comThreshold != -1) || (selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual "
                 << "selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus)  || (automated1)) {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
            cerr << endl << "ERROR: Not allowed in combination with window values." << endl << endl;
            appearErrors = true;
        }

        else if(clusters != -1) {
            cerr << endl << "ERROR: Only one method to chose representative sequences can be applied." << endl << endl;
            appearErrors = true;
        }

        else {
            if(utils::isNumber(argv[++*i])) {
                maxIdentity = atof(argv[*i]);
                if((maxIdentity < 0) || (maxIdentity > 1)) {
                    cerr << endl << "ERROR: The maximum identity threshold should be between 0 and 1." << endl << endl;
                    appearErrors = true;
                }
            }
            else {
                cerr << endl << "ERROR: The minimal positions value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
    }
}

bool trimalArgumentParser::clusters_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-clusters")) && (i+1 != argc) && (clusters == -1)) {

        if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
           (comThreshold != -1) || (selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual "
                 << "selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus)  || (automated1)) {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
            cerr << endl << "ERROR: Not allowed in combination with window values." << endl << endl;
            appearErrors = true;
        }

        else if(maxIdentity != -1) {
            cerr << endl << "ERROR: Only one method to chose representative sequences can be applied." << endl << endl;
            appearErrors = true;
        }

        else {
            if(utils::isNumber(argv[++*i])) {
                clusters = atoi(argv[*i]);
                if(clusters < 1) {
                    cerr << endl << "ERROR: There is a problem with the given clusters number." << endl << endl;
                    appearErrors = true;
                }
            }
            else {
                cerr << endl << "ERROR: The clusters number should be a positive integer number." << endl << endl;
                appearErrors = true;
            }
        }
    }
}

bool trimalArgumentParser::terminal_only_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-terminalonly")) && (!terminal)) {
        terminal = true;
    }
}

bool trimalArgumentParser::window_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-w") && (i+1 != argc) && (windowSize == -1)){

        if((gapWindow != -1) || (simWindow != -1) || (conWindow != -1)) {
            cerr << endl << "ERROR: Not allowed in combination with this specific window value." << endl << endl;
            appearErrors = true;
        }

        if((selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: It's imposible to use this windows size in combination of manual selection method." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
            cerr << endl << "ERROR: Not allowed in combination of automatic methods." << endl << endl;
            appearErrors = true;
        }

        else {
            if(utils::isNumber(argv[*i+1])) {
                windowSize = atoi(argv[++*i]);
                if(windowSize <= 0){
                    cerr << endl << "ERROR: The window value should be a positive integer number." << endl << endl;
                    appearErrors = true;
                }
            }
            else {
                cerr << endl << "ERROR: The window value should be a number." << endl << endl;
                appearErrors = true;
            }
        }
        //~ i++;
    }
}

bool trimalArgumentParser::gap_window_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-gw") && (i+1 != argc) && (gapWindow == -1)){

        if(windowSize != -1) {
            cerr << endl << "ERROR: Not allowed in combination of general window value." << endl << endl;
            appearErrors = true;
        }

        if((selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: It's imposible to use this windows size in combination of manual selection method." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
            cerr << endl << "ERROR: Not allowed in combination of automatic methods." << endl << endl;
            appearErrors = true;
        }

        else {
            if(utils::isNumber(argv[*i+1])) {
                gapWindow = atoi(argv[++*i]);
                if(gapWindow <= 0){
                    cerr << endl << "ERROR: The window value should be a positive integer number." << endl << endl;
                    appearErrors = true;
                }
            }
            else {
                cerr << endl << "ERROR: The window value should be a number." << endl << endl;
                appearErrors = true;
            }
        }
        //~ i++;
    }
}

bool trimalArgumentParser::sim_window_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-sw") && (i+1 != argc) && (simWindow == -1)){

        if(windowSize != -1) {
            cerr << endl << "ERROR: Not allowed in combination of general window value." << endl << endl;
            appearErrors = true;
        }

        if((selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: It's imposible to use this windows size in combination of manual selection method." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
            cerr << endl << "ERROR: Not allowed in combination of automatic methods." << endl << endl;
            appearErrors = true;
        }

        else {
            if(utils::isNumber(argv[*i+1])) {
                simWindow = atoi(argv[++*i]);
                if(simWindow <= 0){
                    cerr << endl << "ERROR: The window value should be a positive integer number." << endl << endl;
                    appearErrors = true;
                }
            }
            else {
                cerr << endl << "ERROR: The window value should be a number." << endl << endl;
                appearErrors = true;
            }
        }
        //~ i++;
    }
}

bool trimalArgumentParser::con_window_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-cw") && (i+1 != argc) && (conWindow == -1)){

        if(windowSize != -1) {
            cerr << endl << "ERROR: Not allowed in combination of general window value." << endl << endl;
            appearErrors = true;
        }

        if((selectCols) || (selectSeqs)) {
            cerr << endl << "ERROR: It's imposible to use this windows size in combination of manual selection method." << endl << endl;
            appearErrors = true;
        }

        else {
            if(utils::isNumber(argv[*i+1])) {
                conWindow = atoi(argv[++*i]);
                if(conWindow <= 0){
                    cerr << endl << "ERROR: The window value should be a positive integer number." << endl << endl;
                    appearErrors = true;
                }
            }
            else {
                cerr << endl << "ERROR: The window value should be a number." << endl << endl;
                appearErrors = true;
            }
        }
        //~ i++;
    }
}

bool trimalArgumentParser::block_argument(int *argc, char *argv[], int *i) {
    if(!strcmp(argv[*i], "-block") && (i+1 != argc) && (blockSize == -1)){

        if(selectCols) {
            cerr << endl << "ERROR: It's imposible to set a block size value in combination with a column manual selection" << endl << endl;
            appearErrors = true;
        }

        else if(conserve != -1) {
            cerr << endl << "ERROR: It's imposible to ask for a minimum percentage of the input newAlignment  in combination with column block size" << endl << endl;
            appearErrors = true;
        }

            //~ else if((nogaps) || (noallgaps) || (strict) || (strictplus) || (automated1)) {
        else if((nogaps) || (noallgaps)) {
            cerr << endl << "ERROR: Not allowed in combination of automatic methods." << endl << endl;
            appearErrors = true;
        }

        else {
            if(utils::isNumber(argv[*i+1])) {
                blockSize = atoi(argv[++*i]);
                if(blockSize <= 0){
                    cerr << endl << "ERROR: The block size value should be a positive integer number." << endl << endl;
                    appearErrors = true;
                }
            }
            else {
                cerr << endl << "ERROR: The block size value should be a number." << endl << endl;
                appearErrors = true;
            }
        }
    }
}

bool trimalArgumentParser::stats_arguments(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-sgc")) && (!sgc)) {
        sgc = true;
        stats--;
    }
    else if((!strcmp(argv[*i], "-sgt")) && (!sgt)) {
        sgt = true;
        stats--;
    }
    else if((!strcmp(argv[*i], "-ssc")) && (!scc)) {
        scc = true;
        stats--;
    }
    else if((!strcmp(argv[*i], "-sst")) && (!sct)) {
        sct = true;
        stats--;
    }
    else if((!strcmp(argv[*i], "-sident")) && (!sident)) {
        sident = true;
        stats--;
    }
    else if((!strcmp(argv[*i], "-sfc")) && (!sfc)) {

        if(infile != NULL) {
            cerr << endl << "ERROR: Not allowed in combination with -in option." << endl << endl;
            appearErrors = true;
            i++;
        }

        else {
            sfc = true;
            stats--;
        }
    }
    else if((!strcmp(argv[*i], "-sft")) && (!sft)) {

        if(infile != NULL) {
            cerr << endl << "ERROR: Not allowed in combination with -in option." << endl << endl;
            appearErrors = true;
            i++;
        }

        else {
            sft = true;
            stats--;
        }
    }

}

bool trimalArgumentParser::complementary_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-complementary")) && (complementary == false)) {
        complementary = true;
    }

}

bool trimalArgumentParser::col_numbering_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-colnumbering")) && (colnumbering == false)) {
        colnumbering = true;
    }
}

bool trimalArgumentParser::split_by_stop_codon_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-splitbystopcodon")) && (splitbystop == false)) {
        splitbystop = true;
    }
}

bool trimalArgumentParser::ignore_stop_codon_argument(int *argc, char *argv[], int *i) {
    if((!strcmp(argv[*i], "-ignorestopcodon")) && (ignorestop == false)) {
        ignorestop = true;
    }
}


