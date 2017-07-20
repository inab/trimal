//
// Created by bioinfo on 8/06/17.
//

#include "../include/trimalArgumentParser.h"
#include <string.h>
#include "../include/compareFiles.h"
#include "../include/defines.h"
void menu();
void examples();
void trimalArgumentParser::parseArguments(int argc, char *argv[])
{

    origAlig = new newAlignment();

    if (argc == 1)
    {
        menu();
        return;
    }

    for(int i = 1; i < argc; i++ )
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
        if (similarity_threshold_argument(&argc, argv, &i)) continue;
        if (consistency_threshold_argument(&argc, argv, &i)) continue;

        if (conservation_argument(&argc, argv, &i)) continue;
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

        cerr << endl << "ERROR: Parameter \"" << argv[i] << "\" not valid or repeated." << endl << endl;
        appearErrors = true;
        break;
    }
}

bool trimalArgumentParser::info_arguments(int *argc, char *argv[], int *i)
{

    if (!strcmp(argv[*i], "-h"))
    {
//         cout << "Print help" << endl;
        menu();
        return true;
    }

    if(!strcmp(argv[*i], "--version"))
    {

//         cout << "Print version" << endl;
        examples();
        return true;
    }

    return false;
}

bool trimalArgumentParser::in_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-in") && (i+1 != argc) && (infile == NULL))
    {

        if((sfc) || (sft) || (consistencyThreshold != -1))
        {
            cerr << endl << "ERROR: Not allowed in combination of file comparision." << endl << endl;
            appearErrors = true;
            i++;
        }

        else if((compareset == -1) || (forceFile != NULL))
        {
            argumentLength = strlen(argv[++*i]);
            infile = new char[argumentLength + 1];
            strcpy(infile, argv[*i]);

            if(!origAlig -> ReadWrite -> loadAlignment(infile))
            {
                cerr << endl << "ERROR: newAlignment  not loaded: \"" << infile << "\" Check the file's content." << endl << endl;
                appearErrors = true;
            }
        }

        else
        {
            if(compareset != -1)
                cerr << endl << "ERROR: Option \"" << argv[*i] << "\" not valid. A reference file exists with alignments to compare." << endl << endl;
            if(forceFile != NULL)
                cerr << endl << "ERROR: Option \"" << argv[*i] << "\" not valid. An alignment file has been setting up to be compare with a set of alignmets." << endl << endl;
            appearErrors = true;
            i++;
        }
        return true;
    }
    return false;
}

bool trimalArgumentParser::out_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-out")) && (i+1 != argc) && (outfile == NULL))
    {
        argumentLength = strlen(argv[++*i]);
        outfile = new char[argumentLength + 1];
        strcpy(outfile, argv[*i]);
        return true;
    }
    return false;
}

bool trimalArgumentParser::html_out_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-htmlout")) && (i+1 != argc) && (htmlOutFile == NULL))
    {
        argumentLength = strlen(argv[++*i]);
        htmlOutFile = new char[argumentLength + 1];
        strcpy(htmlOutFile, argv[*i]);
        return true;
    }
    return false;
}

bool trimalArgumentParser::out_format_arguments(int *argc, char *argv[], int *i)
{

    if (outformat == -1)
    {
        /* Option -clustal -------------------------------------------------------------------------------------- */
        if(!strcmp(argv[*i], "-clustal") )
        {
            outformat = 1;
            return true;
        }

        /* Option -fasta -------------------------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-fasta") )
        {
            outformat = 8;
            return true;
        }

        /* Option -fasta-m10 -------------------------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-fasta_m10") )
        {
            outformat = 8;
            shortNames = true;
            return true;
        }

        /* Option -nbrf ------------------------------------------------------------------------------------ */
        else if(!strcmp(argv[*i], "-nbrf") )
        {
            outformat = 3;
            return true;
        }

        /* Option -nexus ------------------------------------------------------------------------------------ */
        else if(!strcmp(argv[*i], "-nexus") )
        {
            outformat = 17;
            return true;
        }

        /* Option -mega ------------------------------------------------------------------------------------ */
        else if(!strcmp(argv[*i], "-mega") )
        {
            outformat = 21;
            return true;
        }

        /* Option -phylip3.2 --------------------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-phylip3.2") )
        {
            outformat = 11;
            return true;
        }

        /* Option -phylip3.2-m10 ----------------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-phylip3.2_m10") )
        {
            outformat = 11;
            shortNames = true;
            return true;
        }

        /* Option -phylip --------------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-phylip") )
        {
            outformat = 12;
            return true;
        }

        /* Option -phylip-m10 ----------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-phylip_m10") )
        {
            outformat = 12;
            shortNames = true;
            return true;
        }

        /* Option -phylip_paml ---------------------------------------------------------------------- */
        else if(!strcmp(argv[*i], "-phylip_paml") )
        {
            outformat = 13;
            return true;
        }

        /* Option -phylip_paml-m10 ------------------------------------------------------------------ */
        else if(!strcmp(argv[*i], "-phylip_paml_m10") )
        {
            outformat = 13;
            shortNames = true;
            return true;
        }
    }
    return false;

}

bool trimalArgumentParser::matrix_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-matrix") && (i+1 != argc) && (matrixFile == NULL))
    {
        argumentLength = strlen(argv[++*i]);
        matrixFile = new char[argumentLength + 1];
        strcpy(matrixFile, argv[*i]);
        return true;
    }
    return false;
}

bool trimalArgumentParser::compareset_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-compareset") && (i+1 != argc) && (compareset == -1))
    {

        if(infile == NULL)
        {
            compare.open(argv[++*i], ifstream::in);
            if(!compare)
            {
                cerr << endl << "ERROR: Check the reference file with the alignments to compare." << endl << endl;
                appearErrors = true;
            }

            while(compare.getline(line, 256)) numfiles++;
            compare.close();

            compareset = *i;
        }

        else
        {
            cerr << endl << "ERROR: Option \"" << argv[*i] << "\" not valid. A single alignment file has been set by the user." << endl << endl;
            appearErrors = true;
            i++;
        }
        return true;
    }
    return false;
}

bool trimalArgumentParser::force_select_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-forceselect") && (i+1 != argc) && (forceFile == NULL))
    {

        if(infile == NULL)
        {
            argumentLength = strlen(argv[++*i]);
            forceFile = new char[argumentLength + 1];
            strcpy(forceFile, argv[*i]);
            if(!origAlig -> ReadWrite -> loadAlignment(forceFile))
            {
                cerr << endl << "ERROR: alignment not loaded: \"" << forceFile << "\" Check the file's content." << endl << endl;
                appearErrors = true;
            }
        }

        else
        {
            cerr << endl << "ERROR: Option \"" << argv[*i] << "\" not valid. A single alignment file has been setting it up" << endl << endl;
            appearErrors = true;
            i++;
        }
        return true;
    }
    return false;
}

bool trimalArgumentParser::back_trans_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-backtrans") && (i+1 != argc) && (backtransFile == NULL))
    {

        argumentLength = strlen(argv[++*i]);
        backtransFile = new char[argumentLength + 1];
        strcpy(backtransFile, argv[*i]);

        backtranslationAlig = new newAlignment();
        if(!backtranslationAlig -> ReadWrite -> loadAlignment(backtransFile))
        {
            cerr << endl << "ERROR: alignment not loaded: \"" << backtransFile << "\" Check the file's content." << endl << endl;
            appearErrors = true;
        }
        return true;
    }
    return false;
}

bool trimalArgumentParser::gap_threshold_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-gapthreshold") || !strcmp(argv[*i], "-gt")) && (i+1 != argc) && (gapThreshold == -1))
    {

        if((selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1))
        {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
        {
            if(utils::isNumber(argv[++*i]))
            {
                gapThreshold = 1. - atof(argv[*i]);
                if((gapThreshold < 0) || (gapThreshold > 1))
                {
                    cerr << endl << "ERROR: The gap threshold value should be between 0 and 1." << endl << endl;
                    appearErrors = true;
                }
            }
            else
            {
                cerr << endl << "ERROR: The gap threshold value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
        return true;
    }
    return false;
}

bool trimalArgumentParser::similarity_threshold_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-simthreshold") || !strcmp(argv[*i], "-st")) && (i+1 != argc) && (similarityThreshold == -1))
    {

        if((selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1))
        {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
        {
            if(utils::isNumber(argv[++*i]))
            {
                similarityThreshold = atof(argv[*i]);
                if((similarityThreshold < 0) || (similarityThreshold > 1))
                {
                    cerr << endl << "ERROR: The similarity threshold value should be between 0 and 1." << endl << endl;
                    appearErrors = true;
                }
            }
            else
            {
                cerr << endl << "ERROR: The similarity threshold value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
        return true;
    }
    return false;
}

bool trimalArgumentParser::consistency_threshold_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-conthreshold") || !strcmp(argv[*i], "-ct")) && (i+1 != argc) && (consistencyThreshold == -1))
    {

        if((selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else if(infile != NULL)
        {
            cerr << endl << "ERROR: Not allowed in combination with -in option." << endl << endl;
            appearErrors = true;

        }

        else
        {
            if(utils::isNumber(argv[++*i]))
            {
                consistencyThreshold = atof(argv[*i]);
                if((consistencyThreshold < 0) || (consistencyThreshold > 1))
                {
                    cerr << endl << "ERROR: The consistency threshold value should be between 0 and 1." << endl << endl;
                    appearErrors = true;
                }
            }
            else
            {
                cerr << endl << "ERROR: The consistency threshold value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
        return true;
    }
    return false;
}

bool trimalArgumentParser::conservation_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-cons")) && (i+1 != argc) && (conservationThreshold == -1))
    {

        if((selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else if(blockSize != -1)
        {
            cerr << endl << "ERROR: Not allowed in combination of column block size value." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus)  || (automated1))
        {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
        {
            if(utils::isNumber(argv[++*i]))
            {
                conservationThreshold = atof(argv[*i]);
                if((conservationThreshold < 0) || (conservationThreshold > 100))
                {
                    cerr << endl << "ERROR: The minimal positions value should be between 0 and 100." << endl << endl;
                    appearErrors = true;
                }
            }
            else
            {
                cerr << endl << "ERROR: The minimal positions value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
        return true;
    }
    return false;
}
//TODO
bool trimalArgumentParser::select_cols_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-selectcols")) && (selectCols == false) &&
            ((i+3) < argc) && (!strcmp(argv[++*i], "{")) && (!strcmp(argv[*i+2], "}")))
    {

        if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1))
        {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed." << endl << endl;
            appearErrors = true;
        }

        else if(blockSize != -1)
        {
            cerr << endl << "ERROR: Not allowed in combination of column block size value." << endl << endl;
            appearErrors = true;
        }

        else if((gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1) || (consistencyThreshold != -1))
        {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods." << endl << endl;
            appearErrors = true;
        }

        else if((windowSize != -1) || (gapWindow != -1)|| (similarityWindow != -1))
        {
            cerr << endl << "ERROR: It's imposible to use windows size in combination of selection method." << endl << endl;
            appearErrors = true;
        }

        else if((delColumns = utils::readNumbers(argv[++*i])) == NULL)
        {
            cerr << endl << "ERROR: Impossible to parse the sequences number" << endl << endl;
            appearErrors = true;
        }

        else selectCols = true;
        i++;

        return true;
    }
    return false;
}

bool trimalArgumentParser::no_gaps_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-nogaps") && (!nogaps))
    {

        if((windowSize != -1) || (gapWindow != -1) || (similarityWindow != -1))
        {
            cerr << endl << "ERROR: Not allowed in combination of window values." << endl << endl;
            appearErrors = true;
        }

        else if(blockSize != -1)
        {
            cerr << endl << "ERROR: Not allowed in combination of column block size value." << endl << endl;
            appearErrors = true;
        }

        else if((noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1))
        {
            cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
            appearErrors = true;
        }

        else if((gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1) ||
                (selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
            nogaps = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::no_all_gaps_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-noallgaps") && (!noallgaps))
    {

        if((windowSize != -1) || (gapWindow != -1) || (similarityWindow != -1))
        {
            cerr << endl << "ERROR: Not allowed in combination of window values." << endl << endl;
            appearErrors = true;
        }

        else if(blockSize != -1)
        {
            cerr << endl << "ERROR: Not allowed in combination of column block size value." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (gappyout) || (strict) || (strictplus) || (automated1))
        {
            cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
            appearErrors = true;
        }

        else if((gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1) ||
                (selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
            noallgaps = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::keep_seqs_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-keepseqs") && (!keepSeqs))
    {
        keepSeqs = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::keep_header_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-keepheader") && (!keepHeader))
    {
        keepHeader = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::gappy_out_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-gappyout") && (!strict))
    {

        if((windowSize != -1) || (gapWindow != -1) || (similarityWindow != -1))
        {
            cerr << endl << "ERROR: Not allowed in combination of window values." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (strict) || (strictplus) || (automated1))
        {
            cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
            appearErrors = true;
        }

        else if((gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1) ||
                (selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
            gappyout = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::strict_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-strict") && (!strict))
    {

        if((windowSize != -1) || (gapWindow != -1) || (similarityWindow != -1))
        {
            cerr << endl << "ERROR: Not allowed in combination of window values." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strictplus) || (automated1))
        {
            cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
            appearErrors = true;
        }

        else if((gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1) ||
                (selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
            strict = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::strict_plus_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-strictplus")) && (!strictplus))
    {

        if((windowSize != -1) || (gapWindow != -1) || (similarityWindow != -1))
        {
            cerr << endl << "ERROR: Not allowed in combination with this window value." << endl << endl;
            appearErrors = true;
        }

        //~ else if(blockSize != -1) {
        //~ cerr << endl << "ERROR: Not allowed in combination of column block size value." << endl << endl;
        //~ appearErrors = true;
        //~ }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (automated1))
        {
            cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
            appearErrors = true;
        }

        //~ else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
        //~ (comThreshold != -1) || (selectCols) || (selectSeqs)) {
        //~ cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        //~ appearErrors = true;
        //~ }
        else if((gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1) ||
                (selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
            strictplus = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::automated1_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-automated1")) && (!automated1))
    {

        if((windowSize != -1) || (gapWindow != -1) || (similarityWindow != -1))
        {
            cerr << endl << "ERROR: Not allowed in combination with this window value." << endl << endl;
            appearErrors = true;
        }

        else if(blockSize != -1)
        {
            cerr << endl << "ERROR: Not allowed in combination of column block size value." << endl << endl;
            appearErrors = true;
        }

        else if((gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1) ||
                (consistencyThreshold != -1) || (selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
            appearErrors = true;
        }

        else if((gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1) ||
                (selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else
            automated1 = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::residue_overlap_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-resoverlap")) && (i+1 != argc) && (residuesOverlap == -1))
    {

        if((selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Not allowed in combination of methods such as manual selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else
        {
            if(utils::isNumber(argv[++*i]))
            {
                residuesOverlap = atof(argv[*i]);
                if((residuesOverlap < 0) || (residuesOverlap > 1))
                {
                    cerr << endl << "ERROR: The residue overlap value should be between 0 and 1." << endl << endl;
                    appearErrors = true;
                }
            }
            else
            {
                cerr << endl << "ERROR: The residue overlap value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
        return true;
    }
    return false;
}

bool trimalArgumentParser::sequence_overlap_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-seqoverlap")) && (i+1 != argc) && (sequenceOverlap == -1))
    {

        if((selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Not allowed in combination of methods such as manual selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else
        {
            if(utils::isNumber(argv[++*i]))
            {
                sequenceOverlap = atof(argv[*i]);
                if((sequenceOverlap < 0) || (sequenceOverlap > 100))
                {
                    cerr << endl << "ERROR: The sequences overlap value should be between 0 and 100." << endl << endl;
                    appearErrors = true;
                }
            }
            else
            {
                cerr << endl << "ERROR: The minimal positions value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
        return true;
    }
    return false;
}

bool trimalArgumentParser::seqs_select_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-selectseqs")) && (selectSeqs == false) && ((i+3) < argc) && (!strcmp(argv[++*i], "{")) && (!strcmp(argv[*i+2], "}")))
    {

        if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1))
        {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed." << endl << endl;
            appearErrors = true;
        }

        else if((gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1) || (consistencyThreshold != -1))
        {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods." << endl << endl;
            appearErrors = true;
        }

        else if((windowSize != -1) || (gapWindow != -1)|| (similarityWindow != -1))
        {
            cerr << endl << "ERROR: It's imposible to use this windows size in combination of selection method." << endl << endl;
            appearErrors = true;
        }

        else if((clusters != -1) || (maxIdentity != -1))
        {
            cerr << endl << "ERROR: Only one method to chose sequences can be applied." << endl << endl;
            appearErrors = true;
        }

        else if((delSequences = utils::readNumbers(argv[++*i])) == NULL)
        {
            cerr << endl << "ERROR: Impossible to parser the sequences number" << endl << endl;
            appearErrors = true;
        }

        else selectSeqs = true;
        i++;
        return true;
    }
    return false;
}

bool trimalArgumentParser::max_identity_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-maxidentity")) && (i+1 != argc) && (maxIdentity == -1))
    {

        if((gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1) ||
                (consistencyThreshold != -1) || (selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual "
                 << "selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus)  || (automated1))
        {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else if((windowSize != -1) || (gapWindow != -1) || (similarityWindow != -1))
        {
            cerr << endl << "ERROR: Not allowed in combination with window values." << endl << endl;
            appearErrors = true;
        }

        else if(clusters != -1)
        {
            cerr << endl << "ERROR: Only one method to chose representative sequences can be applied." << endl << endl;
            appearErrors = true;
        }

        else
        {
            if(utils::isNumber(argv[++*i]))
            {
                maxIdentity = atof(argv[*i]);
                if((maxIdentity < 0) || (maxIdentity > 1))
                {
                    cerr << endl << "ERROR: The maximum identity threshold should be between 0 and 1." << endl << endl;
                    appearErrors = true;
                }
            }
            else
            {
                cerr << endl << "ERROR: The minimal positions value should be a positive real number." << endl << endl;
                appearErrors = true;
            }
        }
        return true;
    }
    return false;
}

bool trimalArgumentParser::clusters_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-clusters")) && (i+1 != argc) && (clusters == -1))
    {

        if((gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1) ||
                (consistencyThreshold != -1) || (selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual "
                 << "selection of sequences/columns." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus)  || (automated1))
        {
            cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
            appearErrors = true;
        }

        else if((windowSize != -1) || (gapWindow != -1) || (similarityWindow != -1))
        {
            cerr << endl << "ERROR: Not allowed in combination with window values." << endl << endl;
            appearErrors = true;
        }

        else if(maxIdentity != -1)
        {
            cerr << endl << "ERROR: Only one method to chose representative sequences can be applied." << endl << endl;
            appearErrors = true;
        }

        else
        {
            if(utils::isNumber(argv[++*i]))
            {
                clusters = atoi(argv[*i]);
                if(clusters < 1)
                {
                    cerr << endl << "ERROR: There is a problem with the given clusters number." << endl << endl;
                    appearErrors = true;
                }
            }
            else
            {
                cerr << endl << "ERROR: The clusters number should be a positive integer number." << endl << endl;
                appearErrors = true;
            }
        }
        return true;
    }
    return false;
}

bool trimalArgumentParser::terminal_only_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-terminalonly")) && (!terminalOnly))
    {
        terminalOnly = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::window_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-w") && (i+1 != argc) && (windowSize == -1))
    {

        if((gapWindow != -1) || (similarityWindow != -1) || (consistencyWindow != -1))
        {
            cerr << endl << "ERROR: Not allowed in combination with this specific window value." << endl << endl;
            appearErrors = true;
        }

        if((selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: It's imposible to use this windows size in combination of manual selection method." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1))
        {
            cerr << endl << "ERROR: Not allowed in combination of automatic methods." << endl << endl;
            appearErrors = true;
        }

        else
        {
            if(utils::isNumber(argv[*i+1]))
            {
                windowSize = atoi(argv[++*i]);
                if(windowSize <= 0)
                {
                    cerr << endl << "ERROR: The window value should be a positive integer number." << endl << endl;
                    appearErrors = true;
                }
            }
            else
            {
                cerr << endl << "ERROR: The window value should be a number." << endl << endl;
                appearErrors = true;
            }
        }
        //~ i++;
        return true;
    }
    return false;
}

bool trimalArgumentParser::gap_window_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-gw") && (i+1 != argc) && (gapWindow == -1))
    {

        if(windowSize != -1)
        {
            cerr << endl << "ERROR: Not allowed in combination of general window value." << endl << endl;
            appearErrors = true;
        }

        if((selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: It's imposible to use this windows size in combination of manual selection method." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1))
        {
            cerr << endl << "ERROR: Not allowed in combination of automatic methods." << endl << endl;
            appearErrors = true;
        }

        else
        {
            if(utils::isNumber(argv[*i+1]))
            {
                gapWindow = atoi(argv[++*i]);
                if(gapWindow <= 0)
                {
                    cerr << endl << "ERROR: The window value should be a positive integer number." << endl << endl;
                    appearErrors = true;
                }
            }
            else
            {
                cerr << endl << "ERROR: The window value should be a number." << endl << endl;
                appearErrors = true;
            }
        }
        //~ i++;
        return true;
    }
    return false;
}

bool trimalArgumentParser::similarity_window_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-sw") && (i+1 != argc) && (similarityWindow == -1))
    {

        if(windowSize != -1)
        {
            cerr << endl << "ERROR: Not allowed in combination of general window value." << endl << endl;
            appearErrors = true;
        }

        if((selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: It's imposible to use this windows size in combination of manual selection method." << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1))
        {
            cerr << endl << "ERROR: Not allowed in combination of automatic methods." << endl << endl;
            appearErrors = true;
        }

        else
        {
            if(utils::isNumber(argv[*i+1]))
            {
                similarityWindow = atoi(argv[++*i]);
                if(similarityWindow <= 0)
                {
                    cerr << endl << "ERROR: The window value should be a positive integer number." << endl << endl;
                    appearErrors = true;
                }
            }
            else
            {
                cerr << endl << "ERROR: The window value should be a number." << endl << endl;
                appearErrors = true;
            }
        }
        //~ i++;
        return true;
    }
    return false;
}

bool trimalArgumentParser::consistency_window_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-cw") && (i+1 != argc) && (consistencyWindow == -1))
    {

        if(windowSize != -1)
        {
            cerr << endl << "ERROR: Not allowed in combination of general window value." << endl << endl;
            appearErrors = true;
        }

        if((selectCols) || (selectSeqs))
        {
            cerr << endl << "ERROR: It's imposible to use this windows size in combination of manual selection method." << endl << endl;
            appearErrors = true;
        }

        else
        {
            if(utils::isNumber(argv[*i+1]))
            {
                consistencyWindow = atoi(argv[++*i]);
                if(consistencyWindow <= 0)
                {
                    cerr << endl << "ERROR: The window value should be a positive integer number." << endl << endl;
                    appearErrors = true;
                }
            }
            else
            {
                cerr << endl << "ERROR: The window value should be a number." << endl << endl;
                appearErrors = true;
            }
        }
        //~ i++;
        return true;
    }
    return false;
}

bool trimalArgumentParser::block_argument(int *argc, char *argv[], int *i)
{
    if(!strcmp(argv[*i], "-block") && (i+1 != argc) && (blockSize == -1))
    {

        if(selectCols)
        {
            cerr << endl << "ERROR: It's imposible to set a block size value in combination with a column manual selection" << endl << endl;
            appearErrors = true;
        }

        else if(conservationThreshold != -1)
        {
            cerr << endl << "ERROR: It's imposible to ask for a minimum percentage of the input newAlignment  in combination with column block size" << endl << endl;
            appearErrors = true;
        }

        else if((nogaps) || (noallgaps))
        {
            cerr << endl << "ERROR: Not allowed in combination of automatic methods." << endl << endl;
            appearErrors = true;
        }

        else
        {
            if(utils::isNumber(argv[*i+1]))
            {
                blockSize = atoi(argv[++*i]);
                if(blockSize <= 0)
                {
                    cerr << endl << "ERROR: The block size value should be a positive integer number." << endl << endl;
                    appearErrors = true;
                }
            }
            else
            {
                cerr << endl << "ERROR: The block size value should be a number." << endl << endl;
                appearErrors = true;
            }
        }
        return true;
    }
    return false;
}

bool trimalArgumentParser::stats_arguments(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-sgc")) && (!sgc))
    {
        sgc = true;
        stats--;
    }
    else if((!strcmp(argv[*i], "-sgt")) && (!sgt))
    {
        sgt = true;
        stats--;
    }
    else if((!strcmp(argv[*i], "-ssc")) && (!scc))
    {
        scc = true;
        stats--;
    }
    else if((!strcmp(argv[*i], "-sst")) && (!sct))
    {
        sct = true;
        stats--;
    }
    else if((!strcmp(argv[*i], "-sident")) && (!sident))
    {
        sident = true;
        stats--;
    }
    else if((!strcmp(argv[*i], "-sfc")) && (!sfc))
    {

        if(infile != NULL)
        {
            cerr << endl << "ERROR: Not allowed in combination with -in option." << endl << endl;
            appearErrors = true;
            i++;
        }

        else
        {
            sfc = true;
            stats--;
        }
    }
    else if((!strcmp(argv[*i], "-sft")) && (!sft))
    {

        if(infile != NULL)
        {
            cerr << endl << "ERROR: Not allowed in combination with -in option." << endl << endl;
            appearErrors = true;
            i++;
        }

        else
        {
            sft = true;
            stats--;
        }
        return true;
    }
    return false;
}

bool trimalArgumentParser::complementary_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-complementary")) && (getComplementary == false))
    {
        getComplementary = true;
        return true;
    }
    return false;

}

bool trimalArgumentParser::col_numbering_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-colnumbering")) && (columnNumbering == false))
    {
        columnNumbering = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::split_by_stop_codon_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-splitbystopcodon")) && (splitByStopCodon == false))
    {
        splitByStopCodon = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::ignore_stop_codon_argument(int *argc, char *argv[], int *i)
{
    if((!strcmp(argv[*i], "-ignorestopcodon")) && (ignoreStopCodon == false))
    {
        ignoreStopCodon = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::post_process(char* argv[])
{
    if (appearErrors) return true;

    check_force_selection();
    check_input_file_with_coding_sequences_argument();
    check_file_aligned();
    check_similarity_matrix();
    check_outputs_coincidence();
    check_col_numbering();
    check_residue_and_sequence_overlap();
    check_html_output_interest();
    check_output_file_with_statistics();
    check_combinations_among_thresholds();
    check_automated_manual_incompatibilities();
    check_multiple_files_comparison(argv);
    check_block_size();
    check_backtranslations();
    check_coding_sequences_type();
    check_ignore_or_splitby_stop_codon();
    check_and_prepare_coding_sequence();
    check_correspondence();
    check_cw_argument();

    if(appearErrors)
    {

        delete singleAlig;
        delete origAlig;
        delete[] compareAlignmentsArray;

        delete similMatrix;
        delete[] delColumns;

        delete[] filesToCompare;
        delete[] compareVect;

        delete[] outfile;
        delete[] htmlOutFile;

        delete[] infile;
        delete[] matrixFile;

        if(forceFile != NULL) delete forceFile;
        if(backtransFile != NULL) delete backtransFile;
        if(backtranslationAlig != NULL) delete backtranslationAlig;

        return true;
    }
    return false;
}

bool trimalArgumentParser::check_force_selection()
{
    if (!appearErrors)
    {
        if((infile != NULL) && (forceFile != NULL))
        {
            cerr << endl << "ERROR: You can not use a single alignment at the same "
                 << "time that you force the alignment selection." << endl << endl;
            appearErrors = true;
            return true;
        }
        /* ------------------------------------------------------------------------------------------------------ */
        if((compareset == -1) && (forceFile != NULL))
        {
            cerr << endl << "ERROR: You can not force the alignment selection without set"
                 << " an alignment  dataset against to compare it." << endl << endl;
            appearErrors = true;
            return true;
        }
    }
    return false;
}
//TODO TODO
bool trimalArgumentParser::check_input_file_with_coding_sequences_argument()
{
    if((!appearErrors) && (infile == NULL) && (compareset == -1) && (forceFile == NULL) && (backtransFile != NULL))
    {
        cerr << endl << "ERROR: It is impossible to use a Coding Sequences file to apply the back translation method"
             << " without define an input alignment." << endl << endl;
        appearErrors = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::check_file_aligned()
{
    if((!appearErrors) && (infile != NULL))
    {

        if(((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1) ||
                (gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1) || (selectCols) || (selectSeqs) ||
                (residuesOverlap != -1) || (sequenceOverlap != -1) || (stats < 0)) &&
                (!origAlig -> isFileAligned()))
        {
            cerr << endl << "ERROR: The sequences in the input alignment should be aligned in order to use trimming method." << endl << endl;
            appearErrors = true;
            return true;
        }
    }
    return false;
}

bool trimalArgumentParser::check_similarity_matrix()
{
    if((matrixFile != NULL) && (!appearErrors))
    {
        if((!strict) && (!strictplus) && (!automated1) && (similarityThreshold == -1) && (!scc) && (!sct))
        {
            cerr << endl << "ERROR: The Similarity Matrix can only be used with methods that use this matrix." << endl << endl;
            appearErrors = true;
            return true;
        }

        if((gapWindow != -1) ||((compareset == -1) && (consistencyWindow != -1)))
        {
            cerr << endl << "ERROR: The Similarity Matrix can only be used with general/similarity windows size." << endl << endl;
            appearErrors = true;
            return true;
        }
    }
    return false;
}

bool trimalArgumentParser::check_outputs_coincidence()
{
    if((htmlOutFile != NULL) && (outfile != NULL) && (!appearErrors))
    {
        if(!strcmp(htmlOutFile, outfile))
        {
            cerr << endl << "ERROR: The output and html files should not be the same." << endl << endl;
            appearErrors = true;
            return true;
        }
    }
    return false;
}

bool trimalArgumentParser::check_col_numbering()
{
    if((columnNumbering) && (!appearErrors))
    {
        if((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1)
                && (gapThreshold == -1) && (conservationThreshold == -1) && (similarityThreshold == -1) &&  (consistencyThreshold == -1) && (!selectCols) && (!selectSeqs))
        {
            cerr << endl << "ERROR: This parameter can only be used with any trimming method." << endl << endl;
            appearErrors = true;
            return true;
        }
        else if(stats < 0)
        {
            cerr << endl << "ERROR: This parameter is not valid when statistics' parameters are defined." << endl << endl;
            appearErrors = true;
            return true;
        }
    }
    return false;
}

bool trimalArgumentParser::check_residue_and_sequence_overlap()
{
    if(((residuesOverlap != -1) || (sequenceOverlap != -1)) && (!appearErrors))
    {

        if((residuesOverlap != -1) && (sequenceOverlap == -1))
        {
            cerr << endl << "ERROR: The sequence overlap value should be defined." << endl << endl;
            appearErrors = true;
            return true;
        }

        else if((residuesOverlap == -1) && (sequenceOverlap != -1))
        {
            cerr << endl << "ERROR: The residue overlap value should be defined." << endl << endl;
            appearErrors = true;
            return true;
        }
    }
    return false;
}

bool trimalArgumentParser::check_html_output_interest()
{
    if((htmlOutFile != NULL) && (!appearErrors))
    {
        if((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1) &&
                (gapThreshold == -1) && (conservationThreshold == -1) && (similarityThreshold == -1) && (consistencyThreshold == -1) &&
                (!selectCols) && (!selectSeqs) && (residuesOverlap == -1) && (sequenceOverlap == -1) && (maxIdentity == -1) &&
                (clusters == -1))
        {
            cerr << endl << "ERROR: This parameter can only be used with any trimming method." << endl << endl;
            appearErrors = true;
            return true;
        }
    }
    return false;
}

bool trimalArgumentParser::check_output_file_with_statistics()
{
    if((stats < 0) && (!appearErrors))
    {
        stats--;

        if(((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)
                || (gapThreshold != -1) || (conservationThreshold != -1) || (similarityThreshold != -1)) && (outfile == NULL))
        {
            cerr << endl << "ERROR: An output file should be defined in order to get the alignment's statistics." << endl << endl;
            appearErrors = true;
            return true;
        }
    }
    return false;
}

bool trimalArgumentParser::check_combinations_among_thresholds()
{
    if((consistencyThreshold != -1) && (conservationThreshold != -1) && (!appearErrors))
    {

        if((gapThreshold != -1) || (similarityThreshold != -1))
        {
            cerr << endl << "ERROR: Combinations among thresholds are not allowed." << endl << endl;
            appearErrors = true;
            return true;
        }
    }
    return false;
}

bool trimalArgumentParser::check_automated_manual_incompatibilities()
{
    if((getComplementary) && (!appearErrors))
        if((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1)
                && (gapThreshold == -1) && (conservationThreshold == -1) && (similarityThreshold == -1) && (!selectCols) && (!selectSeqs)
                && (residuesOverlap == -1) && (sequenceOverlap == -1) && (maxIdentity == -1) && (clusters == -1))
        {
            cerr << endl << "ERROR: This parameter can only be used with either an automatic or a manual method." << endl << endl;
            appearErrors = true;
            return true;
        }
    /* ------------------------------------------------------------------------------------------------------ */

    /* ------------------------------------------------------------------------------------------------------ */
    if((terminalOnly) && (!appearErrors))
        if((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1)
                && (gapThreshold == -1) && (conservationThreshold == -1) && (similarityThreshold == -1) && (!selectCols) && (!selectSeqs)
                && (residuesOverlap == -1) && (sequenceOverlap == -1) && (maxIdentity == -1) && (clusters == -1))
        {
            cerr << endl << "ERROR: This parameter '-terminalonly' can only be used with either an automatic or a manual method." << endl << endl;
            appearErrors = true;
            return true;
        }
    return false;
}

bool trimalArgumentParser::check_multiple_files_comparison(char* argv[])
{
    bool hasError = false;
    if((compareset != -1) && (!appearErrors))
    {

        compareAlignmentsArray = new newAlignment*[numfiles];
        filesToCompare = new char*[numfiles];

        /* -------------------------------------------------------------------- */
        compare.open(argv[compareset], ifstream::in);

        for(i = 0; (i < numfiles)  && (!appearErrors); i++)
        {

            /* -------------------------------------------------------------------- */
            for(nline.clear(), compare.read(&c, 1); (c != '\n') && ((!compare.eof())); compare.read(&c, 1))
                nline += c;

            filesToCompare[i] = new char [nline.size() + 1];
            strcpy(filesToCompare[i], nline.c_str());
            /* -------------------------------------------------------------------- */

            /* -------------------------------------------------------------------- */
            compareAlignmentsArray[i] = new newAlignment();
            if(!compareAlignmentsArray[i] -> ReadWrite -> loadAlignment(filesToCompare[i]))
            {
                cerr << endl << "alignment not loaded: \"" << filesToCompare[i] << "\" Check the file's content." << endl << endl;
                hasError = true;
            }

            else
            {
                if(!compareAlignmentsArray[i] -> isFileAligned())
                {
                    cerr << endl << "ERROR: The sequences in the input alignment should be aligned in order to use this method." << endl << endl;
                    hasError = true;
                }
                else
                {
                    compareAlignmentsArray[i] -> SequencesMatrix = new sequencesMatrix(compareAlignmentsArray[i]);

                    if(compareAlignmentsArray[i] -> getNumAminos() > maxAminos)
                        maxAminos = compareAlignmentsArray[i] -> getNumAminos();

                    if((compareAlignmentsArray[i] -> getAlignmentType() != prevType) && (prevType != -1))
                    {
                        cerr << endl << "ERROR: The alignments' datatypes are different. Check your dataset." << endl << endl;
                        hasError = true;
                    }
                    else
                        prevType = compareAlignmentsArray[i] -> getAlignmentType();
                }
            }
        }
        /* -------------------------------------------------------------------- */

        /* -------------------------------------------------------------------- */
        if (!appearErrors)
        {
            if(forceFile == NULL)
            {
                compareVect = new float[maxAminos];
                if((stats >= 0) && (outfile != NULL))
                    referFile = compareFiles::algorithm(compareAlignmentsArray, filesToCompare, compareVect, numfiles, true);
                else
                    referFile = compareFiles::algorithm(compareAlignmentsArray, filesToCompare, compareVect, numfiles, false);

                if(windowSize != -1)
                    compareFiles::applyWindow(compareAlignmentsArray[referFile] -> getNumAminos(), windowSize, compareVect);
                else if(consistencyWindow != -1)
                    compareFiles::applyWindow(compareAlignmentsArray[referFile] -> getNumAminos(), consistencyWindow, compareVect);

                origAlig -> ReadWrite -> loadAlignment (filesToCompare[referFile]);
            }
            else
            {
                compareVect = new float[origAlig -> getNumAminos()];
                appearErrors = !compareFiles::forceComparison(compareAlignmentsArray, numfiles, origAlig, compareVect);

                if((windowSize != -1) && (!appearErrors))
                    compareFiles::applyWindow(origAlig -> getNumAminos(), windowSize, compareVect);
                else if((consistencyWindow != -1) && (!appearErrors))
                    compareFiles::applyWindow(origAlig -> getNumAminos(), consistencyWindow, compareVect);
            }
        }
        
        /* -------------------------------------------------------------------- */

        /* -------------------------------------------------------------------- */
        for(i = 0; i < numfiles; i++)
        {
            delete compareAlignmentsArray[i];
            delete filesToCompare[i];
        }
        /* -------------------------------------------------------------------- */
    }

    return hasError;

}
bool trimalArgumentParser::check_block_size()
{
    if((!appearErrors) && (origAlig -> getNumAminos() < (blockSize/4)))
    {
        cerr << endl << "ERROR: The block size value is too big. Please, choose another one smaller than residues number / 4." << endl << endl;
        appearErrors = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::check_backtranslations()
{
    if (!appearErrors)
    {
        if (backtransFile == NULL)
        {
            if(splitByStopCodon)
            {
                cerr << endl << "ERROR: The -splitbystopcodon parameter can be only set up with backtranslation functionality." << endl << endl;
                appearErrors = true;
                return true;
            }
            if(ignoreStopCodon)
            {
                cerr << endl << "ERROR: The -ignorestopcodon parameter can be only set up with backtranslation functionality." << endl << endl;
                appearErrors = true;
                return true;
            }
        }
        else if(!origAlig -> isFileAligned())
        {
            cerr << endl << "ERROR: The input protein file has to be aligned to carry out the backtranslation process" << endl << endl;
            appearErrors = true;
            return true;
        }
    }
    return false;
}

bool trimalArgumentParser::check_coding_sequences_type()
{
    if((!appearErrors) && (backtransFile != NULL) && (backtranslationAlig -> getAlignmentType() != DNAType && backtranslationAlig -> getAlignmentType() != DNADeg))
    {
        cerr << endl << "ERROR: Check your Coding sequences file. It has been detected other kind of biological sequences." << endl << endl;
        appearErrors = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::check_ignore_or_splitby_stop_codon()
{
    if((!appearErrors) && (ignoreStopCodon) && (splitByStopCodon))
    {
        cerr << endl << "ERROR: Incompatibility of -ignorestopcodon & -splitbystopcodon parameters. Choose one." << endl << endl;
        appearErrors = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::check_and_prepare_coding_sequence()
{
    if((!appearErrors)  && (backtransFile != NULL) && 
        (!backtranslationAlig -> prepareCodingSequence(splitByStopCodon, ignoreStopCodon, origAlig)))
    {
        // Error informing is made by prepareCodingSequence function.
        appearErrors = true;
        return true;
    }
    return false;
}

bool trimalArgumentParser::check_correspondence()
{
    if((!appearErrors) && (backtransFile != NULL))
    {

        sequencesNames = new string[backtranslationAlig -> getNumSpecies()];
        sequencesLengths = new int[backtranslationAlig -> getNumSpecies()];
        backtranslationAlig -> getSequences(sequencesNames, sequencesLengths);

        if(origAlig -> checkCorrespondence(sequencesNames, sequencesLengths, backtranslationAlig -> getNumSpecies(), 3) != true)
        {
            appearErrors = true;
            return true;
        }
    }
    return false;
}

void trimalArgumentParser::check_cw_argument()
{
    if((!appearErrors) && (windowSize != -1) && (compareset != -1))
        cerr << "INFO: Try with specific comparison file window value. parameter -cw." << endl << endl;
}

int trimalArgumentParser::perform()
{
    if (appearErrors) return -1;

    if(conservationThreshold == -1)
        conservationThreshold  = 0;
    /* -------------------------------------------------------------------- */

    origAlig -> Cleaning -> setTrimTerminalGapsFlag(terminalOnly);
    origAlig -> setKeepSequencesFlag(keepSeqs);
    origAlig -> setKeepSeqsHeaderFlag(keepHeader);

    /* -------------------------------------------------------------------- */
    set_window_size();

    if(blockSize != -1)
        origAlig -> setBlockSize(blockSize);

    if(outformat != -1)
        origAlig -> setOutputFormat(outformat, shortNames);

    if (create_or_use_similarity_matrix())
        return -1;

    print_statistics();

    if(backtransFile != NULL)
        seqMatrix = origAlig -> SequencesMatrix;
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    clean_alignment();
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if(singleAlig == NULL)
    {
        singleAlig = origAlig;
        origAlig = NULL;
    }
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if((htmlOutFile != NULL) && (!appearErrors))
        if(!origAlig -> ReadWrite ->
           alignmentSummaryHTML(htmlOutFile, singleAlig -> getNumAminos(), singleAlig -> getNumSpecies(),
                                singleAlig -> getCorrespResidues(), singleAlig -> getCorrespSequences(), compareVect))
        {
            cerr << endl << "ERROR: It's imposible to generate the HTML output file." << endl << endl;
            appearErrors = true;
        }
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if(backtransFile != NULL)
    {

        if(sequencesNames != NULL) delete [] sequencesNames;
        sequencesNames = new string[singleAlig -> getNumSpecies()];

        singleAlig -> getSequences(sequencesNames);

        singleAlig = backtranslationAlig -> getTranslationCDS(singleAlig -> getNumAminos(), singleAlig -> getNumSpecies(),
                     singleAlig -> getCorrespResidues(), sequencesNames, seqMatrix, singleAlig);
    }
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if((outfile != NULL) && (!appearErrors))
    {
        if(!singleAlig -> ReadWrite -> saveAlignment(outfile))
        {
            cerr << endl << "ERROR: It's imposible to generate the output file." << endl << endl;
            appearErrors = true;
        }
    }
    else if((stats >= 0) && (!appearErrors))
        singleAlig -> ReadWrite -> printAlignment();
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if((columnNumbering) && (!appearErrors))
        singleAlig -> Statistics -> printCorrespondence();
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    delete singleAlig;
    delete origAlig;
    delete[] compareAlignmentsArray;

    delete similMatrix;
    delete []delColumns;

    delete[] filesToCompare;
    delete[] compareVect;

    delete[] outfile;
    delete[] htmlOutFile;

    delete[] infile;
    delete[] matrixFile;
    /* -------------------------------------------------------------------- */

    return 0;
}

void trimalArgumentParser::print_statistics()
{
    if(sgc)
    {
        origAlig -> Statistics -> printStatisticsGapsColumns();
        stats++;
        if(stats < -1)
            cout << endl;
    }

    if(sgt)
    {
        origAlig -> Statistics -> printStatisticsGapsTotal();
        stats++;
        if(stats < -1)
            cout << endl;
    }

    if(scc)
    {
        origAlig -> Statistics -> printStatisticsConservationColumns();
        stats++;
        if(stats < -1)
            cout << endl;
    }

    if(sct)
    {
        origAlig -> Statistics -> printStatisticsConservationTotal();
        stats++;
        if(stats < -1)
            cout << endl;
    }

    if(sident)
    {
        origAlig -> printSeqIdentity();
        stats++;
        if(stats < -1)
            cout << endl;
    }

    if(compareset != -1)
    {
        if(sfc)
            compareFiles::printStatisticsFileColumns(origAlig -> getNumAminos(), compareVect);
        if(sft)
            compareFiles::printStatisticsFileAcl(origAlig -> getNumAminos(), compareVect);
    }
}

bool trimalArgumentParser::create_or_use_similarity_matrix()
{
        if((strict) || (strictplus) || (automated1) || (similarityThreshold != -1.0) || (scc == 1) || (sct == 1))
    {
        similMatrix = new similarityMatrix();

        if(matrixFile != NULL)
            similMatrix -> loadSimMatrix(matrixFile);

        else
        {
            if((origAlig -> getAlignmentType()) == AAType)
                similMatrix -> defaultAASimMatrix();
            else
                similMatrix -> defaultNTSimMatrix();
        }

        if(!origAlig -> Statistics -> setSimilarityMatrix(similMatrix))
        {
            cerr << endl << "ERROR: It's impossible to process the Similarity Matrix." << endl << endl;
            return true;
        }
    }
    return false;
}

void trimalArgumentParser::clean_alignment()
{
    if(nogaps)
        singleAlig = origAlig -> Cleaning -> cleanGaps(0, 0, getComplementary);

    else if(noallgaps)
        singleAlig = origAlig -> Cleaning -> cleanNoAllGaps(getComplementary);

    else if(gappyout)
        singleAlig = origAlig -> Cleaning -> clean2ndSlope(getComplementary);

    else if(strict)
        singleAlig = origAlig -> Cleaning -> cleanCombMethods(getComplementary, false);

    else if(strictplus)
        singleAlig = origAlig -> Cleaning -> cleanCombMethods(getComplementary, true);

    else if(automated1)
    {
        if(origAlig -> Cleaning -> selectMethod() == GAPPYOUT)
            singleAlig = origAlig -> Cleaning -> clean2ndSlope(getComplementary);
        else
            singleAlig = origAlig -> Cleaning -> cleanCombMethods(getComplementary, false);
    }
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if(consistencyThreshold != -1)
        singleAlig = origAlig -> Cleaning -> cleanCompareFile(consistencyThreshold, conservationThreshold, compareVect, getComplementary);
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if((residuesOverlap != -1) && (sequenceOverlap != -1))
    {
        singleAlig = origAlig   -> Cleaning -> cleanSpuriousSeq(residuesOverlap, (sequenceOverlap/100), getComplementary);
        singleAlig = singleAlig -> Cleaning -> cleanNoAllGaps(false);
    }
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if(similarityThreshold != -1.0)
    {
        if(gapThreshold != -1.0)
            singleAlig = origAlig -> Cleaning -> clean(conservationThreshold, gapThreshold, similarityThreshold, getComplementary);
        else
            singleAlig = origAlig -> Cleaning -> cleanConservation(conservationThreshold, similarityThreshold, getComplementary);
    }
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    else if(gapThreshold != -1.0)
        singleAlig = origAlig -> Cleaning -> cleanGaps(conservationThreshold, gapThreshold, getComplementary);
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if((selectCols) || (selectSeqs))
    {

        /* -------------------------------------------------------------------- */
        if(delColumns != NULL)
        {
            num = delColumns[0];
            if(delColumns[num] >= origAlig -> getNumAminos())
            {
                cerr << endl << "ERROR: This option only accepts integer numbers between 0 and the number of columns - 1." << endl << endl;
                appearErrors = true;
            }
            else
                singleAlig = origAlig -> Cleaning -> removeColumns(delColumns, 1, num, getComplementary);
        }
        /* -------------------------------------------------------------------- */

        /* -------------------------------------------------------------------- */
        if(delSequences != NULL)
        {
            num = delSequences[0];
            if(delSequences[num] >= origAlig -> getNumSpecies())
            {
                cerr << endl << "ERROR: This option only accepts integer numbers between 0 and the number of sequences - 1." << endl << endl;
                appearErrors = true;
            }
            else
            {
                singleAlig = origAlig -> Cleaning -> removeSequences(delSequences, 1, num, getComplementary);
                singleAlig = singleAlig -> Cleaning -> cleanNoAllGaps(false);
            }
        }
        /* -------------------------------------------------------------------- */
    }
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if(maxIdentity != -1)
    {
        singleAlig = origAlig -> Cleaning -> getClustering (maxIdentity);
        singleAlig = singleAlig -> Cleaning -> cleanNoAllGaps(false);
    }
    else if(clusters != -1)
    {
        if(clusters > origAlig -> getNumSpecies())
        {
            cerr << endl << "ERROR:The number of clusters from the newAlignment  can not be larger than the number of sequences from that alignment." << endl << endl;
            appearErrors = true;
        }
        else
        {
            singleAlig = origAlig -> Cleaning -> getClustering(origAlig -> Cleaning -> getCutPointClusters(clusters));
            singleAlig = singleAlig -> Cleaning -> cleanNoAllGaps(false);
        }
    }
}

void trimalArgumentParser::set_window_size()
{
    if(windowSize != -1)
    {
        gapWindow = windowSize;
        similarityWindow = windowSize;
    }
    else
    {
        if(gapWindow == -1)
            gapWindow = 0;
        if(similarityWindow == -1)
            similarityWindow = 0;
    }
    origAlig -> setWindowsSize(gapWindow, similarityWindow);
}

void trimalArgumentParser::menu(void)
{

    cout << endl;
    cout << "trimAl v" << VERSION << ".rev" << REVISION  << " build[" << BUILD
         << "]. " << AUTHORS << endl << endl;

    cout << "trimAl webpage: http://trimal.cgenomics.org" << endl << endl;

    cout << "This program is free software: you can redistribute it and/or modify " << endl
         << "it under the terms of the GNU General Public License as published by " << endl
         << "the Free Software Foundation, the last available version." << endl << endl;

    cout << "Please cite:" << endl
         << "\t\ttrimAl: a tool for automated newAlignment  trimming in large-scale phylogenetic analyses."
         << "\n\t\tSalvador Capella-Gutierrez; Jose M. Silla-Martinez; Toni Gabaldon."
         << "\n\t\tBioinformatics 2009, 25:1972-1973." << endl << endl;

    cout << "Basic usage" << endl
         << "\ttrimal -in <inputfile> -out <outputfile> -(other options)." << endl << endl;

    cout << "Common options (for a complete list please see the User Guide or visit http://trimal.cgenomics.org):" << endl << endl;
    cout << "    -h                       " << "Print this information and show some examples." << endl;
    cout << "    --version                " << "Print the trimAl version." << endl << endl;

    cout << "    -in <inputfile>          " << "Input file in several formats (clustal, fasta, NBRF/PIR, nexus, phylip3.2, phylip)." << endl << endl;

    cout << "    -compareset <inputfile>  " << "Input list of paths for the files containing the alignments to compare." << endl;
    cout << "    -forceselect <inputfile> " << "Force selection of the given input file in the files comparison method." << endl << endl;

    cout << "    -backtrans <inputfile>   " << "Use a Coding Sequences file to get a backtranslation for a given AA alignment" << endl;
    cout << "    -ignorestopcodon         " << "Ignore stop codons in the input coding sequences" << endl;
    cout << "    -splitbystopcodon        " << "Split input coding sequences up to first stop codon appearance" << endl << endl;


    cout << "    -matrix <inpufile>       " << "Input file for user-defined similarity matrix (default is Blosum62)." << endl << endl;

    cout << "    -out <outputfile>        " << "Output newAlignment  in the same input format (default stdout). (default input format)" << endl;
    cout << "    -htmlout <outputfile>    " << "Get a summary of trimal's work in an HTML file." << endl << endl;

    cout << "    -keepheader              " << "Keep original sequence header including non-alphanumeric characters." << endl;
    cout << "                             " << "Only available for input FASTA format files. (future versions will extend this feature)" << endl << endl;

    cout << "    -nbrf                    " << "Output file in NBRF/PIR format" << endl;
    cout << "    -mega                    " << "Output file in MEGA format" << endl;
    cout << "    -nexus                   " << "Output file in NEXUS format" << endl;
    cout << "    -clustal                 " << "Output file in CLUSTAL format" << endl << endl;

    cout << "    -fasta                   " << "Output file in FASTA format" << endl;
    cout << "    -fasta_m10               " << "Output file in FASTA format. Sequences name length up to 10 characters." << endl << endl;

    cout << "    -phylip                  " << "Output file in PHYLIP/PHYLIP4 format" << endl;
    cout << "    -phylip_m10              " << "Output file in PHYLIP/PHYLIP4 format. Sequences name length up to 10 characters." << endl;
    cout << "    -phylip_paml             " << "Output file in PHYLIP format compatible with PAML" << endl;
    cout << "    -phylip_paml_m10         " << "Output file in PHYLIP format compatible with PAML. Sequences name length up to 10 characters." << endl;
    cout << "    -phylip3.2               " << "Output file in PHYLIP3.2 format" << endl;
    cout << "    -phylip3.2_m10           " << "Output file in PHYLIP3.2 format. Sequences name length up to 10 characters." << endl << endl;

    cout << "    -complementary           " << "Get the complementary alignment." << endl;
    cout << "    -colnumbering            " << "Get the relationship between the columns in the old and new alignment." << endl << endl;

    cout << "    -selectcols { n,l,m-k }  " << "Selection of columns to be removed from the alignment. Range: [0 - (Number of Columns - 1)]. (see User Guide)." << endl;
    cout << "    -selectseqs { n,l,m-k }  " << "Selection of sequences to be removed from the alignment. Range: [0 - (Number of Sequences - 1)]. (see User Guide)." << endl << endl;

    cout << "    -gt -gapthreshold <n>    " << "1 - (fraction of sequences with a gap allowed). Range: [0 - 1]" << endl;
    cout << "    -st -simthreshold <n>    " << "Minimum average similarity allowed. Range: [0 - 1]" << endl;
    cout << "    -ct -conthreshold <n>    " << "Minimum consistency value allowed.Range: [0 - 1]" << endl;
    cout << "    -cons <n>                " << "Minimum percentage of the positions in the original newAlignment  to conserve. Range: [0 - 100]" << endl << endl;

    cout << "    -nogaps                  " << "Remove all positions with gaps in the alignment." << endl;
    cout << "    -noallgaps               " << "Remove columns composed only by gaps." << endl;
    cout << "    -keepseqs                " << "Keep sequences even if they are composed only by gaps." << endl << endl;

    cout << "    -gappyout                " << "Use automated selection on \"gappyout\" mode. This method only uses "
         << "information based on gaps' distribution. (see User Guide)." << endl;
    cout << "    -strict                  " << "Use automated selection on \"strict\" mode. (see User Guide)." << endl;
    cout << "    -strictplus              " << "Use automated selection on \"strictplus\" mode. (see User Guide)."  << endl;
    cout << "                             " << "(Optimized for Neighbour Joining phylogenetic tree reconstruction)."<< endl << endl;

    cout << "    -automated1              " << "Use a heuristic selection of the automatic method based on similarity statistics. "
         << "(see User Guide). (Optimized for Maximum Likelihood phylogenetic tree reconstruction)."
         << endl << endl;

    cout << "    -terminalonly            " << "Only columns out of internal boundaries (first and last column without gaps) are " << endl;
    cout << "                             " << "candidated to be trimmed depending on the applied method" << endl;

    cout << "    -block <n>               " << "Minimum column block size to be kept in the trimmed alignment. Available with manual"
         << " and automatic (gappyout) methods" << endl << endl;


    cout << "    -resoverlap              " << "Minimum overlap of a positions with other positions in the column to be considered a "
         << "\"good position\". Range: [0 - 1]. (see User Guide)." << endl;
    cout << "    -seqoverlap              " << "Minimum percentage of \"good positions\" that a sequence must have in order to be conserved. Range: [0 - 100]"
         << "(see User Guide)." << endl << endl;

    cout << "    -clusters <n>            " << "Get the most Nth representatives sequences from a given alignment. Range: [1 - (Number of sequences)]" << endl;
    cout << "    -maxidentity <n>         " << "Get the representatives sequences for a given identity threshold. Range: [0 - 1]." << endl << endl;

    cout << "    -w <n>                   " << "(half) Window size, score of position i is the average of the window (i - n) to (i + n)."
         << endl;
    cout << "    -gw <n>                  " << "(half) Window size only applies to statistics/methods based on Gaps." << endl;
    cout << "    -sw <n>                  " << "(half) Window size only applies to statistics/methods based on Similarity." << endl;
    cout << "    -cw <n>                  " << "(half) Window size only applies to statistics/methods based on Consistency." << endl << endl;

    cout << "    -sgc                     " << "Print gap scores for each column in the input alignment." << endl;
    cout << "    -sgt                     " << "Print accumulated gap scores for the input alignment." << endl;
    cout << "    -ssc                     " << "Print similarity scores for each column in the input alignment." << endl;
    cout << "    -sst                     " << "Print accumulated similarity scores for the input alignment." << endl;
    cout << "    -sfc                     " << "Print sum-of-pairs scores for each column from the selected alignment"
         << endl;
    cout << "    -sft                     " << "Print accumulated sum-of-pairs scores for the selected alignment"
         << endl;
    cout << "    -sident                  " << "Print identity scores for all sequences in the input alignment. (see User Guide)."
         << endl << endl;
}

void trimalArgumentParser::examples(void)
{

    cout << "Some Examples:" << endl << endl;

    cout << "1) Removes all positions in the newAlignment  with gaps in 10% or more of" << endl
         << "   the sequences, unless this leaves less than 60% of original alignment. " << endl
         << "   In such case, print the 60% best (with less gaps) positions." << endl << endl;

    cout << "   trimal -in <inputfile> -out <outputfile> -gt 0.9 -cons 60" << endl << endl;

    cout << "2) As above but, the gap score is averaged over a window starting" << endl
         << "   3 positions before and ending 3 positions after each column." << endl << endl;

    cout << "   trimal -in <inputfile> -out <outputfile> -gt 0.9 -cons 60 -w 3" << endl << endl;

    cout << "3) Use an automatic method to decide optimal thresholds, based in the gap scores" << endl
         << "   from input alignment. (see User Guide for details)." << endl << endl;

    cout << "   trimal -in <inputfile> -out <outputfile> -gappyout" << endl << endl;

    cout << "4) Use automatic methods to decide optimal thresholds, based on the combination " << endl
         << "   of gap and similarity scores. (see User Guide for details)." << endl << endl;

    cout << "   trimal -in <inputfile> -out <outputfile> -strictplus" << endl << endl;

    cout << "5) Use an heuristic to decide the optimal method for trimming the alignment. " << endl
         << "   (see User Guide for details)." << endl << endl;

    cout << "   trimal -in <inputfile> -out <outputfile> -automated1" << endl << endl;

    cout << "6) Use residues and sequences overlap thresholds to delete some sequences from the " << endl
         << "   alignemnt. (see User Guide for details)." << endl << endl;

    cout << "   trimal -in <inputfile> -out <outputfile> -resoverlap 0.8 -seqoverlap 75" << endl << endl;

    cout << "7) Selection of columns to be deleted from the alignment. The selection can " << endl
         << "   be a column number or a column number interval. Start from 0" << endl << endl;

    cout << "   trimal -in <inputfile> -out <outputfile> -selectcols { 0,2,3,10,45-60,68,70-78 }" << endl << endl;

    cout << "8) Get the complementary newAlignment  from the newAlignment  previously trimmed." << endl << endl;

    cout << "   trimal -in <inputfile> -out <outputfile> -selectcols { 0,2,3,10,45-60,68,70-78 } -complementary" << endl << endl;

    cout << "9) Selection of sequences to be deleted from the alignment. Start in 0" << endl << endl;

    cout << "   trimal -in <inputfile> -out <outputfile> -selectseqs { 2,4,8-12 } " << endl << endl;

    cout << "10) Select the 5 most representative sequences from the alignment" << endl << endl;

    cout << "   trimal -in <inputfile> -out <outputfile> -clusters 5 " << endl << endl;
}


