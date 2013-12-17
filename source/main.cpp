/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    2009-2013 Capella-Gutierrez S. and Gabaldon, T.
              [scapella, tgabaldon]@crg.es

    This file is part of trimAl.

    trimAl is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl. If not, see <http://www.gnu.org/licenses/>.

***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#include <fstream>
#include <iostream>
#include <iomanip>

#include <stdlib.h>
#include <string.h>

#include "compareFiles.h"
#include "compareFiles.h"
#include "alignment.h"
#include "defines.h"
#include "utils.h"

void menu(void);
void examples(void);

int main(int argc, char *argv[]){

  /* Parameters Control */
  bool appearErrors = false, complementary = false, colnumbering = false, nogaps = false, noallgaps = false, gappyout = false,
       strict = false, strictplus = false, automated1 = false, sgc = false, sgt = false, scc = false, sct = false, sfc = false,
       sft = false, sident = false, selectSeqs = false, selectCols = false, shortNames = false, splitbystop = false, terminal = false,
       keepSeqs = false, keepHeader = false, ignorestop = false;

  float conserve = -1, gapThreshold = -1, simThreshold = -1, comThreshold = -1, resOverlap = -1, seqOverlap = -1, maxIdentity = -1;
  int outformat = -1, prevType = -1, compareset = -1, stats = 0, windowSize = -1, gapWindow = -1, simWindow = -1, conWindow = -1,
      blockSize = -1, clusters = -1;

  /* Others varibles */
  ifstream compare;
  float *compareVect = NULL;
  alignment **compAlig  = NULL;
  string nline, *seqNames = NULL;
  sequencesMatrix *seqMatrix = NULL;
  similarityMatrix *similMatrix = NULL;
  alignment *origAlig = NULL, *singleAlig = NULL, *backtranslation = NULL;

  int i = 1, lng, num = 0, maxAminos = 0, numfiles = 0, referFile = 0, *delColumns = NULL, *delSequences = NULL, *seqLengths = NULL;
  char c, *forceFile = NULL, *infile = NULL, *backtransFile = NULL, *outfile = NULL, *outhtml = NULL, *matrix = NULL,
       **filesToCompare = NULL, line[256];

  /* ------------------------------------------------------------------------------------------------------ */

  /* Exec: TrimAl - Shows the menu. */

  /* ------------------------------------------------------------------------------------------------------ */
  if(argc == 1) {
    menu();
    return 0;
  }

  /* ------------------------------------------------------------------------------------------------------ */

  /*                                        Help and Version Menu                                           */

  /* ------------------------------------------------------------------------------------------------------ */
  if(!strcmp(argv[i], "-h") && (i+1 == argc)) {
    menu(); examples();
    return 0;
  }

  if(!strcmp(argv[i], "--version") && (i+1 == argc)) {
    cout << endl << "trimAl v" << VERSION << ".rev" << REVISION
         << " build[" << BUILD << "]" << endl << endl;
    return 0;
  }

  /***** ***** ***** ***** ***** ***** ***** Parameters Processing ***** ***** ***** ***** ***** ***** *****/
  origAlig = new alignment;

  while(i < argc) {

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                Input and Output files and format output                                */

   /* Option -in ------------------------------------------------------------------------------------------- */
    if(!strcmp(argv[i], "-in") && (i+1 != argc) && (infile == NULL)) {

      if((sfc) || (sft) || (comThreshold != -1)) {
        cerr << endl << "ERROR: Not allowed in combination of file comparision." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((compareset == -1) || (forceFile != NULL)) {
        lng = strlen(argv[++i]);
        infile = new char[lng + 1];
        strcpy(infile, argv[i]);

        if(!origAlig -> loadAlignment(infile)) {
          cerr << endl << "ERROR: Alignment not loaded: \"" << infile << "\" Check the file's content." << endl << endl;
          appearErrors = true;
        }
      }

      else {
        if(compareset != -1)
          cerr << endl << "ERROR: Option \"" << argv[i] << "\" not valid. A reference file exists with alignments to compare." << endl << endl;
        if(forceFile != NULL)
          cerr << endl << "ERROR: Option \"" << argv[i] << "\" not valid. A alignment file has been setting up to be compare with a set of alignmets." << endl << endl;
        appearErrors = true;
        i++;
      }
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -out ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-out")) && (i+1 != argc) && (outfile == NULL)) {
      lng = strlen(argv[++i]);
      outfile = new char[lng + 1];
      strcpy(outfile, argv[i]);
    }

   /* Option -htmlout -------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-htmlout")) && (i+1 != argc) && (outhtml == NULL)) {
      lng = strlen(argv[++i]);
      outhtml = new char[lng + 1];
      strcpy(outhtml, argv[i]);
    }

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                           Output File format                                           */

   /* Option -clustal -------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-clustal") && (outformat == -1))
      outformat = 1;

   /* Option -fasta -------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-fasta") && (outformat == -1))
      outformat = 8;

   /* Option -fasta-m10 -------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-fasta_m10") && (outformat == -1)) {
      outformat = 8; shortNames = true;
   }

   /* Option -nbrf ------------------------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-nbrf") && (outformat == -1))
      outformat = 3;

   /* Option -nexus ------------------------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-nexus") && (outformat == -1))
      outformat = 17;

   /* Option -mega ------------------------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-mega") && (outformat == -1))
      outformat = 21;

   /* Option -phylip3.2 --------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-phylip3.2") && (outformat == -1))
      outformat = 11;

   /* Option -phylip3.2-m10 ----------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-phylip3.2_m10") && (outformat == -1)) {
      outformat = 11; shortNames = true;
    }

   /* Option -phylip --------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-phylip") && (outformat == -1))
      outformat = 12;

   /* Option -phylip-m10 ----------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-phylip_m10") && (outformat == -1)) {
      outformat = 12; shortNames = true;
    }

   /* Option -phylip_paml ---------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-phylip_paml") && (outformat == -1))
      outformat = 13;

   /* Option -phylip_paml-m10 ------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-phylip_paml_m10") && (outformat == -1)) {
      outformat = 13; shortNames = true;
    }

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                         Similarity Matrix File                                         */

   /* Option -matrix --------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-matrix") && (i+1 != argc) && (matrix == NULL)) {
      lng = strlen(argv[++i]);
      matrix = new char[lng + 1];
      strcpy(matrix, argv[i]);
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /*                                   File with a alignments' set to compare                               */

   /* Option -compareset ----------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-compareset") && (i+1 != argc) && (compareset == -1)) {

      if(infile == NULL) {
        compare.open(argv[++i], ifstream::in);
        if(!compare) {
          cerr << endl << "ERROR: Check the reference file with the alignments to compare." << endl << endl;
          appearErrors = true;
        }

        while(compare.getline(line, 256)) numfiles++;
        compare.close();

        compareset = i;
      }

      else {
        cerr << endl << "ERROR: Option \"" << argv[i] << "\" not valid. A single alignment file has been set by the user." << endl << endl;
        appearErrors = true;
        i++;
      }
    }

    /* Option -forceselect ----------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-forceselect") && (i+1 != argc) && (forceFile == NULL)) {

      if(infile == NULL) {
        lng = strlen(argv[++i]);
        forceFile = new char[lng + 1];
        strcpy(forceFile, argv[i]);
        if(!origAlig -> loadAlignment(forceFile)) {
          cerr << endl << "ERROR: Alignment not loaded: \"" << forceFile << "\" Check the file's content." << endl << endl;
          appearErrors = true;
        }
      }

      else {
        cerr << endl << "ERROR: Option \"" << argv[i] << "\" not valid. A single alignment file has been setting it up" << endl << endl;
        appearErrors = true;
        i++;
      }
    }

    /* Option -backtrans -------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-backtrans") && (i+1 != argc) && (backtransFile == NULL)) {

      lng = strlen(argv[++i]);
      backtransFile = new char[lng + 1];
      strcpy(backtransFile, argv[i]);

      backtranslation = new alignment;
      if(!backtranslation -> loadAlignment(backtransFile)) {
        cerr << endl << "ERROR: Alignment not loaded: \"" << backtransFile << "\" Check the file's content." << endl << endl;
        appearErrors = true;
      }
    }

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                  Manual Method Values. Deleting columns                                */

   /* Option -gt, gapthreshold ----------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-gapthreshold") || !strcmp(argv[i], "-gt")) && (i+1 != argc) && (gapThreshold == -1)) {

      if((selectCols) || (selectSeqs)) {
        cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns." << endl << endl;
        appearErrors = true;
      }

      else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
        cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        appearErrors = true;
      }

      else {
        if(utils::isNumber(argv[++i])) {
          gapThreshold = 1 - atof(argv[i]);
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
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -st -simthreshold ----------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-simthreshold") || !strcmp(argv[i], "-st")) && (i+1 != argc) && (simThreshold == -1)) {

      if((selectCols) || (selectSeqs)) {
        cerr << endl << "ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns." << endl << endl;
        appearErrors = true;
      }

      else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
        cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        appearErrors = true;
      }

      else {
        if(utils::isNumber(argv[++i])) {
          simThreshold = atof(argv[i]);
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
   /* ------------------------------------------------------------------------------------------------------ */


   /* Option -ct -conthreshold ----------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-conthreshold") || !strcmp(argv[i], "-ct")) && (i+1 != argc) && (comThreshold == -1)) {

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
        if(utils::isNumber(argv[++i])) {
          comThreshold = atof(argv[i]);
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

   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -cons ----------------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-cons")) && (i+1 != argc) && (conserve == -1)) {

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
        if(utils::isNumber(argv[++i])) {
          conserve = atof(argv[i]);
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
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -selectcols -------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-selectcols")) && (selectCols == false) && ((i+3) < argc) && (!strcmp(argv[++i], "{")) && (!strcmp(argv[i+2], "}"))) {

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

      else if((delColumns = utils::readNumbers(argv[++i])) == NULL) {
        cerr << endl << "ERROR: Impossible to parser the sequences number" << endl << endl;
        appearErrors = true;
      }

      else selectCols = true;
      i++;
    }

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                   Automated Methods. Deleting Columns                                  */

   /* Option -nogaps --------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-nogaps") && (!nogaps)) {

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
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -noallgaps --------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-noallgaps") && (!noallgaps)) {

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

   /* Option -keepseqs --------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-keepseqs") && (!keepSeqs)) {
      keepSeqs = true;
    }

   /* Option -keepseqs --------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-keepheader") && (!keepHeader)) {
      keepHeader = true;
    }

   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -gappyout ------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-gappyout") && (!strict)) {

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
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -strict --------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-strict") && (!strict)) {

      if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
        cerr << endl << "ERROR: Not allowed in combination of window values." << endl << endl;
        appearErrors = true;
      }

      //~ else if(blockSize != -1) {
        //~ cerr << endl << "ERROR: Not allowed in combination of column block size value." << endl << endl;
        //~ appearErrors = true;
      //~ }

      else if((nogaps) || (noallgaps) || (gappyout) || (strictplus) || (automated1)) {
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
        strict = true;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -strictplus ----------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-strictplus")) && (!strictplus)) {

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
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -automated1 ----------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-automated1")) && (!automated1)) {

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
//~
      //~ else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || (comThreshold != -1) || (delColumns != NULL)) {
        //~ cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        //~ appearErrors = true;
      //~ }

      else
        automated1 = true;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                 Manual Method Values. Deleting sequences                                */

   /* Option -coloverlap ------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-resoverlap")) && (i+1 != argc) && (resOverlap == -1)) {

      if((selectCols) || (selectSeqs)) {
        cerr << endl << "ERROR: Not allowed in combination of methods such as manual selection of sequences/columns." << endl << endl;
        appearErrors = true;
      }

      else {
        if(utils::isNumber(argv[++i])) {
          resOverlap = atof(argv[i]);
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
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -seqoverlap ----------------------------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-seqoverlap")) && (i+1 != argc) && (seqOverlap == -1)) {

      if((selectCols) || (selectSeqs)) {
        cerr << endl << "ERROR: Not allowed in combination of methods such as manual selection of sequences/columns." << endl << endl;
        appearErrors = true;
      }

      else {
        if(utils::isNumber(argv[++i])) {
          seqOverlap = atof(argv[i]);
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

   /* Option -selectseqs -------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-selectseqs")) && (selectSeqs == false) && ((i+3) < argc) && (!strcmp(argv[++i], "{")) && (!strcmp(argv[i+2], "}"))) {

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

      else if((delSequences = utils::readNumbers(argv[++i])) == NULL) {
        cerr << endl << "ERROR: Impossible to parser the sequences number" << endl << endl;
        appearErrors = true;
      }

      else selectSeqs = true;
      i++;
    }

   /* Option -maxidentity ----------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-maxidentity")) && (i+1 != argc) && (maxIdentity == -1)) {

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
        if(utils::isNumber(argv[++i])) {
          maxIdentity = atof(argv[i]);
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
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -clusters ----------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-clusters")) && (i+1 != argc) && (clusters == -1)) {

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
        if(utils::isNumber(argv[++i])) {
          clusters = atoi(argv[i]);
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
   /* ------------------------------------------------------------------------------------------------------ */

   /* ------------------------------------------------------------------------------------------------------ */

   /* Other methods: Just remove the terminal gaps from an alignment keeping the columns that are in the middle
    * of the sequences independently of the trimming method used */

   /* ------------------------------------------------------------------------------------------------------ */
   /* Option -terminalonly --------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-terminalonly")) && (!terminal)) {
      terminal = true;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                           Windows Size Values                                           */

   /* Option -w -------------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-w") && (i+1 != argc) && (windowSize == -1)){

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
        if(utils::isNumber(argv[i+1])) {
          windowSize = atoi(argv[++i]);
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
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -gw -------------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-gw") && (i+1 != argc) && (gapWindow == -1)){

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
        if(utils::isNumber(argv[i+1])) {
          gapWindow = atoi(argv[++i]);
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
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -sw -------------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-sw") && (i+1 != argc) && (simWindow == -1)){

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
        if(utils::isNumber(argv[i+1])) {
          simWindow = atoi(argv[++i]);
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
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -cw -------------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-cw") && (i+1 != argc) && (conWindow == -1)){

      if(windowSize != -1) {
        cerr << endl << "ERROR: Not allowed in combination of general window value." << endl << endl;
        appearErrors = true;
      }

      if((selectCols) || (selectSeqs)) {
        cerr << endl << "ERROR: It's imposible to use this windows size in combination of manual selection method." << endl << endl;
        appearErrors = true;
      }

      else {
        if(utils::isNumber(argv[i+1])) {
          conWindow = atoi(argv[++i]);
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
   /* ------------------------------------------------------------------------------------------------------ */

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                             Block Size Value                                           */

   /* Option -block -------------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-block") && (i+1 != argc) && (blockSize == -1)){

      if(selectCols) {
        cerr << endl << "ERROR: It's imposible to set a block size value in combination with a column manual selection" << endl << endl;
        appearErrors = true;
      }

      else if(conserve != -1) {
        cerr << endl << "ERROR: It's imposible to ask for a minimum percentage of the input alignment in combination with column block size" << endl << endl;
        appearErrors = true;
      }

      //~ else if((nogaps) || (noallgaps) || (strict) || (strictplus) || (automated1)) {
      else if((nogaps) || (noallgaps)) {
        cerr << endl << "ERROR: Not allowed in combination of automatic methods." << endl << endl;
        appearErrors = true;
      }

      else {
        if(utils::isNumber(argv[i+1])) {
          blockSize = atoi(argv[++i]);
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
   /* ------------------------------------------------------------------------------------------------------ */

   /*                                               Statistics                                               */

   /* Option -sgc ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-sgc")) && (!sgc)) {
      sgc = true;
      stats--;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -sgt ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-sgt")) && (!sgt)) {
      sgt = true;
      stats--;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -scc ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-ssc")) && (!scc)) {
      scc = true;
      stats--;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -sct ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-sst")) && (!sct)) {
      sct = true;
      stats--;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -sident --------------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-sident")) && (!sident)) {
      sident = true;
      stats--;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -sfc ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-sfc")) && (!sfc)) {

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
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -sft ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-sft")) && (!sft)) {

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
   /* ------------------------------------------------------------------------------------------------------ */

   /*                                            Others parameters                                           */

   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -complementary -------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-complementary")) && (complementary == false)) {
      complementary = true;
    }

   /* Option -colnumbering ------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-colnumbering")) && (colnumbering == false)) {
      colnumbering = true;
    }

   /* Option -splitbystopcodon ------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-splitbystopcodon")) && (splitbystop == false)) {
      splitbystop = true;
    }

   /* Option -ignorestopcodon ------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-ignorestopcodon")) && (ignorestop == false)) {
      ignorestop = true;
    }

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                          Not Valids Parameters                                         */

   /* ------------------------------------------------------------------------------------------------------ */
    else {
      cerr << endl << "ERROR: Parameter \"" << argv[i] << "\" not valid." << endl << endl;
      appearErrors = true;
    }
   /* ------------------------------------------------------------------------------------------------------ */
    i++;

    if(appearErrors)
      break;

  }
   /* ------------------------------------------------------------------------------------------------------ */

   /*                                       Postprocessing Parameters                                        */


  /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors) && (infile != NULL) && (forceFile != NULL)) {
     cerr << endl << "ERROR: You can not use a single alignmet at the same "
        << "time that you force the alignment selection." << endl << endl;
     appearErrors = true;
  }
  /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors) && (compareset == -1) && (forceFile != NULL)) {
     cerr << endl << "ERROR: You can not force the alignment selection without set"
        << " an alignment dataset against to compare it." << endl << endl;
     appearErrors = true;
  }
  /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors) && (infile == NULL) && (compareset == -1) && (forceFile == NULL) && (backtransFile != NULL)) {
     cerr << endl << "ERROR: It is impossible to use a Coding Sequences file to apply the back translation method"
              << " without define an input alignment." << endl << endl;
     appearErrors = true;
  }
  /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors) && (infile != NULL)) {

    if(((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1) ||
      (gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || (selectCols) || (selectSeqs) ||
      (resOverlap != -1) || (seqOverlap != -1) || (stats < 0)) &&
      (!origAlig -> isFileAligned())) {
        cerr << endl << "ERROR: The sequences in the input alignment should be aligned in order to use trimming method." << endl << endl;
        appearErrors = true;
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors) && (windowSize != -1) && (compareset != -1))
    cerr << "INFO: Try with specific comparison file window value. parameter -cw." << endl << endl;
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */
  if((matrix != NULL) && (!appearErrors)) {
    if((!strict) && (!strictplus) && (!automated1) && (simThreshold == -1.0) && (!scc) && (!sct)) {
      cerr << endl << "ERROR: The Similarity Matrix can only be used with methods that use this matrix." << endl << endl;
      appearErrors = true;
    }

    if((gapWindow != -1) ||((compareset == -1) && (conWindow != -1))) {
      cerr << endl << "ERROR: The Similarity Matrix can only be used with general/similarity windows size." << endl << endl;
      appearErrors = true;
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */
  if((complementary) && (!appearErrors))
    if((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1)
      && (gapThreshold == -1) && (conserve == -1) && (simThreshold == -1) && (!selectCols) && (!selectSeqs)
    && (resOverlap == -1) && (seqOverlap == -1) && (maxIdentity == -1) && (clusters == -1)) {
      cerr << endl << "ERROR: This parameter can only be used with either an automatic or a manual method." << endl << endl;
      appearErrors = true;
    }
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */
  if((terminal) && (!appearErrors))
    if((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1)
      && (gapThreshold == -1) && (conserve == -1) && (simThreshold == -1) && (!selectCols) && (!selectSeqs)
    && (resOverlap == -1) && (seqOverlap == -1) && (maxIdentity == -1) && (clusters == -1)) {
      cerr << endl << "ERROR: This parameter '-terminalonly' can only be used with either an automatic or a manual method." << endl << endl;
      appearErrors = true;
    }
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */
  if((colnumbering) && (!appearErrors)) {
    if((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1)
      && (gapThreshold == -1) && (conserve == -1) && (simThreshold == -1) &&  (comThreshold == -1) && (!selectCols) && (!selectSeqs)) {
      cerr << endl << "ERROR: This parameter can only be used with any trimming method." << endl << endl;
      appearErrors = true;
    }
    else if(stats < 0) {
      cerr << endl << "ERROR: This parameter is not valid when statistics' parameters are defined." << endl << endl;
      appearErrors = true;
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */
  if((outhtml != NULL) && (outfile != NULL) && (!appearErrors)) {
    if(!strcmp(outhtml, outfile)) {
      cerr << endl << "ERROR: The output and html files should not be the same." << endl << endl;
      appearErrors = true;
    }
  }

  /* ------------------------------------------------------------------------------------------------------ */

  if((outhtml != NULL) && (!appearErrors)) {
   if((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1) &&
      (gapThreshold == -1) && (conserve == -1) && (simThreshold == -1) && (comThreshold == -1) &&
      (!selectCols) && (!selectSeqs) && (resOverlap == -1) && (seqOverlap == -1) && (maxIdentity == -1) &&
    (clusters == -1)) {
      cerr << endl << "ERROR: This parameter can only be used with any trimming method." << endl << endl;
      appearErrors = true;
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */

  //~ if((outhtml != NULL) && (!appearErrors)) {
   //~ if(((gapThreshold != -1) || (simThreshold != -1)) && (comThreshold != -1)) {
      //~ cerr << endl << "ERROR: Impossible to generate the HTML file using two consecutive trimming methods." << endl << endl;
      //~ appearErrors = true;
    //~ }
  //~ }
  /* ------------------------------------------------------------------------------------------------------ */


  /* ------------------------------------------------------------------------------------------------------ */
  if(((resOverlap != -1) || (seqOverlap != -1)) && (!appearErrors)) {

    if((resOverlap != -1) && (seqOverlap == -1)) {
      cerr << endl << "ERROR: The sequence overlap value should be defined." << endl << endl;
      appearErrors = true;
    }

    else if((resOverlap == -1) && (seqOverlap != -1)) {
      cerr << endl << "ERROR: The residue overlap value should be defined." << endl << endl;
      appearErrors = true;
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */
  if((stats < 0) && (!appearErrors)) {
    stats--;

    if(((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)
      || (gapThreshold != -1) || (conserve != -1) || (simThreshold != -1)) && (outfile == NULL)) {
      cerr << endl << "ERROR: An output file should be defined in order to get the alignment's statistics." << endl << endl;
      appearErrors = true;
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */
  if((comThreshold != -1) && (conserve != -1) && (!appearErrors)) {

    if((gapThreshold != -1) || (simThreshold != -1)) {
      cerr << endl << "ERROR: Combinations among thresholds are not allowed." << endl << endl;
      appearErrors = true;
    }
  }
  /* **** ***** ***** ***** ***** ***** **** **************************** **** ***** ***** ***** ***** ***** **** */

  /* **** ***** ***** ***** ***** ***** ***** Files Comparison Methods ***** ***** ***** ***** ***** ***** **** */
  if((compareset != -1) && (!appearErrors)) {

    compAlig = new alignment*[numfiles];
    filesToCompare = new char*[numfiles];

    /* -------------------------------------------------------------------- */
    compare.open(argv[compareset], ifstream::in);

    for(i = 0; (i < numfiles)  && (!appearErrors); i++) {

      /* -------------------------------------------------------------------- */
      for(nline.clear(), compare.read(&c, 1); (c != '\n') && ((!compare.eof())); compare.read(&c, 1))
        nline += c;

      filesToCompare[i] = new char [nline.size() + 1];
      strcpy(filesToCompare[i], nline.c_str());
      /* -------------------------------------------------------------------- */

      /* -------------------------------------------------------------------- */
      compAlig[i] = new alignment;
      if(!compAlig[i] -> loadAlignment(filesToCompare[i])) {
        cerr << endl << "Alignment not loaded: \"" << filesToCompare[i] << "\" Check the file's content." << endl << endl;
        appearErrors = true;
      }

      else {
        if(!compAlig[i] -> isFileAligned()) {
          cerr << endl << "ERROR: The sequences in the input alignment should be aligned in order to use this method." << endl << endl;
          appearErrors = true;
        } else {
          compAlig[i] -> sequenMatrix();

          if(compAlig[i] -> getNumAminos() > maxAminos)
            maxAminos = compAlig[i] -> getNumAminos();

          if((compAlig[i] -> getTypeAlignment() != prevType) && (prevType != -1)) {
            cerr << endl << "ERROR: The alignments' datatypes are different. Check your dataset." << endl << endl;
            appearErrors = true;
          } else
            prevType = compAlig[i] -> getTypeAlignment();
        }
      }
    }
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if((!appearErrors) && (forceFile == NULL)) {

      compareVect = new float[maxAminos];
      if((stats >= 0) && (outfile != NULL))
        referFile = compareFiles::algorithm(compAlig, filesToCompare, compareVect, numfiles, true);
      else
        referFile = compareFiles::algorithm(compAlig, filesToCompare, compareVect, numfiles, false);

      if(windowSize != -1)
        compareFiles::applyWindow(compAlig[referFile] -> getNumAminos(), windowSize, compareVect);
      else if(conWindow != -1)
        compareFiles::applyWindow(compAlig[referFile] -> getNumAminos(), conWindow, compareVect);

      origAlig -> loadAlignment(filesToCompare[referFile]);

    } else if((!appearErrors) && (forceFile != NULL)) {

      compareVect = new float[origAlig -> getNumAminos()];
      appearErrors = !(compareFiles::forceComparison(compAlig, numfiles, origAlig, compareVect));

      if((windowSize != -1) && (!appearErrors))
        compareFiles::applyWindow(origAlig -> getNumAminos(), windowSize, compareVect);
      else if((conWindow != -1) && (!appearErrors))
        compareFiles::applyWindow(origAlig -> getNumAminos(), conWindow, compareVect);
    }
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    for(i = 0; i < numfiles; i++) {
      delete compAlig[i];
      delete filesToCompare[i];
    }
    /* -------------------------------------------------------------------- */
  }

  /* **** ***** ***** ***** ***** ***** **** **************************** **** ***** ***** ***** ***** ***** **** */

  /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors) && (origAlig -> getNumAminos() < (blockSize/4))) {
     cerr << endl << "ERROR: The block size value is too big. Please, choose another one smaller than residues number / 4." << endl << endl;
     appearErrors = true;
  }

  if((!appearErrors) && (backtransFile != NULL) && (backtranslation -> getTypeAlignment() != DNAType && backtranslation -> getTypeAlignment() != DNADeg)) {
     cerr << endl << "ERROR: Check your Coding sequences file. It has been detected other kind of biological sequences." << endl << endl;
     appearErrors = true;
  }
  /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors) && (origAlig -> isFileAligned() != true) && (backtransFile != NULL)) {
     cerr << endl << "ERROR: The input protein file has to be aligned to carry out the backtranslation process" << endl << endl;
     appearErrors = true;
  }
  /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors) && (backtransFile == NULL) && (splitbystop)) {
     cerr << endl << "ERROR: The -splitbystopcodon parameter can be only set up with backtranslation functionality." << endl << endl;
     appearErrors = true;
  }
  /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors) && (backtransFile == NULL) && (ignorestop)) {
     cerr << endl << "ERROR: The -ignorestopcodon parameter can be only set up with backtranslation functionality." << endl << endl;
     appearErrors = true;
  }
  /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors) && (ignorestop) && (splitbystop)) {
     cerr << endl << "ERROR: Incompatibility of -ignorestopcodon & -splitbystopcodon parameters. Choose one." << endl << endl;
     appearErrors = true;
  }
  /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors)  && (backtransFile != NULL) && (backtranslation -> prepareCodingSequence(splitbystop, ignorestop, origAlig) != true))
    appearErrors = true;

  /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors) && (backtransFile != NULL)) {

    seqNames = new string[backtranslation -> getNumSpecies()];
    seqLengths = new int[backtranslation -> getNumSpecies()];
    backtranslation -> getSequences(seqNames, seqLengths);

    if(origAlig -> checkCorrespondence(seqNames, seqLengths, backtranslation -> getNumSpecies(), 3) != true)
      appearErrors = true;
  }
  /* ------------------------------------------------------------------------------------------------------ */

  /* **** ***** ***** ***** ***** ***** **** End of Parameters Processing **** ***** ***** ***** ***** ***** **** */


  /* **** ***** ***** ***** ***** ***** **** Errors Control **** ***** ***** ***** ***** ***** **** */
  if(appearErrors) {

    delete singleAlig;
    delete origAlig;
    delete[] compAlig;

    delete similMatrix;
    delete []delColumns;

    delete[] filesToCompare;
    delete[] compareVect;

    delete[] outfile;
    delete[] outhtml;

    delete[] infile;
    delete[] matrix;

    if(forceFile != NULL) delete forceFile;
    if(backtransFile != NULL) delete backtransFile;
    if(backtranslation != NULL) delete backtranslation;

    return -1;
  }
  /* **** ***** ***** ***** ***** ***** ** End Errors Control ** ***** ***** ***** ***** ***** **** */

  /* -------------------------------------------------------------------- */
  if(conserve == -1)
    conserve  = 0;
  /* -------------------------------------------------------------------- */

  origAlig -> trimTerminalGaps(terminal);
  origAlig -> setKeepSequencesFlag(keepSeqs);
  origAlig -> setKeepSeqsHeaderFlag(keepHeader);

  /* -------------------------------------------------------------------- */
  if(windowSize != -1) {
    gapWindow = windowSize;
    simWindow = windowSize;
  }
  else {
    if(gapWindow == -1)
      gapWindow = 0;
    if(simWindow == -1)
      simWindow = 0;
  }
  origAlig -> setWindowsSize(gapWindow, simWindow);

  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(blockSize != -1)
    origAlig -> setBlockSize(blockSize);

  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(outformat != -1)
    origAlig -> setOutputFormat(outformat, shortNames);
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if((strict) || (strictplus) || (automated1) || (simThreshold != -1.0) || (scc == 1) || (sct == 1)) {
    similMatrix = new similarityMatrix();

    if(matrix != NULL)
      similMatrix -> loadSimMatrix(matrix);

    else {
      if((origAlig -> getTypeAlignment()) == AAType)
        similMatrix -> defaultAASimMatrix();
      else
        similMatrix -> defaultNTSimMatrix();
    }

    if(!origAlig -> setSimilarityMatrix(similMatrix)) {
      cerr << endl << "ERROR: It's imposible to proccess the Similarity Matrix." << endl << endl;
      return -1;
    }
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(sgc) {
    origAlig -> printStatisticsGapsColumns();
    stats++;
    if(stats < -1)
      cout << endl;
  }

  if(sgt) {
    origAlig -> printStatisticsGapsTotal();
    stats++;
    if(stats < -1)
      cout << endl;
  }

  if(scc) {
    origAlig -> printStatisticsConservationColumns();
    stats++;
    if(stats < -1)
      cout << endl;
  }

  if(sct) {
    origAlig -> printStatisticsConservationTotal();
    stats++;
    if(stats < -1)
      cout << endl;
  }

  if(sident) {
    origAlig -> printSeqIdentity();
    stats++;
    if(stats < -1)
      cout << endl;
  }

  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(compareset != -1) {
    if(sfc)
      compareFiles::printStatisticsFileColumns(origAlig -> getNumAminos(), compareVect);
    if(sft)
      compareFiles::printStatisticsFileAcl(origAlig -> getNumAminos(), compareVect);
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(backtransFile != NULL)
    seqMatrix = origAlig -> getSeqMatrix();
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(nogaps)
    singleAlig = origAlig -> cleanGaps(0, 0, complementary);

  else if(noallgaps)
    singleAlig = origAlig -> cleanNoAllGaps(complementary);

  else if(gappyout)
    singleAlig = origAlig -> clean2ndSlope(complementary);

  else if(strict)
    singleAlig = origAlig -> cleanCombMethods(complementary, false);

  else if(strictplus)
    singleAlig = origAlig -> cleanCombMethods(complementary, true);

  else if(automated1) {
    if(origAlig -> selectMethod() == GAPPYOUT)
      singleAlig = origAlig -> clean2ndSlope(complementary);
    else
      singleAlig = origAlig -> cleanCombMethods(complementary, false);
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(comThreshold != -1)
    singleAlig = origAlig -> cleanCompareFile(comThreshold, conserve, compareVect, complementary);
 /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if((resOverlap != -1) && (seqOverlap != -1)) {
    singleAlig = origAlig -> cleanSpuriousSeq(resOverlap, (seqOverlap/100), complementary);
    singleAlig = singleAlig -> cleanNoAllGaps(false);
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(simThreshold != -1.0) {
    if(gapThreshold != -1.0)
      singleAlig = origAlig -> clean(conserve, gapThreshold, simThreshold, complementary);
    else
      singleAlig = origAlig -> cleanConservation(conserve, simThreshold, complementary);
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  else if(gapThreshold != -1.0)
    singleAlig = origAlig -> cleanGaps(conserve, gapThreshold, complementary);
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if((selectCols) || (selectSeqs)) {

    /* -------------------------------------------------------------------- */
    if(delColumns != NULL) {
    num = delColumns[0];
      if(delColumns[num] >= origAlig -> getNumAminos()) {
        cerr << endl << "ERROR: This option only accepts integer numbers between 0 and the number of columns - 1." << endl << endl;
        appearErrors = true;
      }
      else
        singleAlig = origAlig -> removeColumns(delColumns, 1, num, complementary);
    }
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if(delSequences != NULL) {
      num = delSequences[0];
      if(delSequences[num] >= origAlig -> getNumSpecies()) {
        cerr << endl << "ERROR: This option only accepts integer numbers between 0 and the number of sequences - 1." << endl << endl;
        appearErrors = true;
      }
      else {
        singleAlig = origAlig -> removeSequences(delSequences, 1, num, complementary);
        singleAlig = singleAlig -> cleanNoAllGaps(false);
      }
    }
    /* -------------------------------------------------------------------- */
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(maxIdentity != -1) {
    singleAlig = origAlig -> getClustering(maxIdentity);
  singleAlig = singleAlig -> cleanNoAllGaps(false);
  }
  else if(clusters != -1) {
  if(clusters > origAlig -> getNumSpecies()) {
        cerr << endl << "ERROR:The number of clusters from the alignment can not be larger than the number of sequences from that alignment." << endl << endl;
        appearErrors = true;
    } else {
    singleAlig = origAlig -> getClustering(origAlig -> getCutPointClusters(clusters));
    singleAlig = singleAlig -> cleanNoAllGaps(false);
  }
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(singleAlig == NULL) {
    singleAlig = origAlig;
    origAlig = NULL;
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if((outhtml != NULL) && (!appearErrors))
    if(!origAlig -> alignmentSummaryHTML(outhtml, singleAlig -> getNumAminos(), singleAlig -> getNumSpecies(),
                                     singleAlig -> getCorrespResidues(), singleAlig -> getCorrespSequences(), compareVect)) {
      cerr << endl << "ERROR: It's imposible to generate the HTML output file." << endl << endl;
      appearErrors = true;
    }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(backtransFile != NULL) {

  if(seqNames != NULL) delete [] seqNames;
    seqNames = new string[singleAlig -> getNumSpecies()];

  singleAlig -> getSequences(seqNames);

  singleAlig = backtranslation -> getTranslationCDS(singleAlig -> getNumAminos(), singleAlig -> getNumSpecies(),
                                                      singleAlig -> getCorrespResidues(), seqNames, seqMatrix, singleAlig);
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if((outfile != NULL) && (!appearErrors)) {
    if(!singleAlig -> saveAlignment(outfile)) {
      cerr << endl << "ERROR: It's imposible to generate the output file." << endl << endl;
      appearErrors = true;
    }
  }
  else if((stats >= 0) && (!appearErrors))
    singleAlig -> printAlignment();
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if((colnumbering) && (!appearErrors))
    singleAlig -> printCorrespondence();
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  delete singleAlig;
  delete origAlig;
  delete[] compAlig;

  delete similMatrix;
  delete []delColumns;

  delete[] filesToCompare;
  delete[] compareVect;

  delete[] outfile;
  delete[] outhtml;

  delete[] infile;
  delete[] matrix;
  /* -------------------------------------------------------------------- */

  return 0;
}

void menu(void) {

  cout << endl;
  cout << "trimAl v" << VERSION << ".rev" << REVISION  << " build[" << BUILD
       << "]. " << AUTHORS << endl << endl;

  cout << "trimAl webpage: http://trimal.cgenomics.org" << endl << endl;

  cout << "This program is free software: you can redistribute it and/or modify " << endl
       << "it under the terms of the GNU General Public License as published by " << endl
       << "the Free Software Foundation, the last available version." << endl << endl;

  cout << "Please cite:" << endl
       << "\t\ttrimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses."
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

  cout << "    -out <outputfile>        " << "Output alignment in the same input format (default stdout). (default input format)" << endl;
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
  cout << "    -cons <n>                " << "Minimum percentage of the positions in the original alignment to conserve. Range: [0 - 100]" << endl << endl;

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

void examples(void) {

  cout << "Some Examples:" << endl << endl;

  cout << "1) Removes all positions in the alignment with gaps in 10% or more of" << endl
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

  cout << "8) Get the complementary alignment from the alignment previously trimmed." << endl << endl;

  cout << "   trimal -in <inputfile> -out <outputfile> -selectcols { 0,2,3,10,45-60,68,70-78 } -complementary" << endl << endl;

  cout << "9) Selection of sequences to be deleted from the alignment. Start in 0" << endl << endl;

  cout << "   trimal -in <inputfile> -out <outputfile> -selectseqs { 2,4,8-12 } " << endl << endl;

  cout << "10) Select the 5 most representative sequences from the alignment" << endl << endl;

  cout << "   trimal -in <inputfile> -out <outputfile> -clusters 5 " << endl << endl;
}

