/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    readAl v1.3: a tool for automated alignment conversion among different
                 formats.

    Copyright (C) 2009 Capella-Gutierrez S. and Gabaldon, T.
                       [scapella, tgabaldon]@crg.es

    This file is part of readAl.

    readAl is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    readAl is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with readAl. If not, see <http://www.gnu.org/licenses/>.

 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#include <stdlib.h>
#include <string.h>

#include "compareFiles.h"
#include "alignment.h"
#include "utils.h"

#define VERSION 1.3
#define REVISION 20100902

#define DNAType 1
#define RNAType 2
#define AAType  3

void menu(void);

int main(int argc, char *argv[]){

  /* Parameters Control */
  alignment compAlig;
  int outformat = -1, i = 1;
  char *infile = NULL, *outfile = NULL;
  bool appearErrors = false, format = false, type = false, reverse = false, shortNames = false;

  /* Exec: readAl - Shows the menu. */
  if(argc == 1) {
    menu();
    return 0;
  }

  /***** ***** ***** ***** ***** ***** ***** Parameters Processing ***** ***** ***** ***** ***** ***** *****/
  /* ------------------------------------------------------------------------------------------------------ */

  /*                                        Help and Version Menu                                           */

  /* ------------------------------------------------------------------------------------------------------ */
  if(!strcmp(argv[i], "-h") && (i+1 == argc)) {
    menu();
    return 0;
  }

  if(!strcmp(argv[i], "--version") && (i+1 == argc)) {
    cout << endl << "readAl " << VERSION << "rev" << REVISION << endl << endl;
    return 0;
  }

  /* ------------------------------------------------------------------------------------------------------ */
  while(i < argc) {

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                Input and Output files and format output                                */

   /* Option -in ------------------------------------------------------------------------------------------- */
    if(!strcmp(argv[i], "-in") && (i+1 != argc) && (infile == NULL)) {

      infile = new char[strlen(argv[++i]) + 1];
      strcpy(infile, argv[i]);

      if(!compAlig.loadAlignment(infile)) {
        cerr << endl << "ERROR: Alignment not loaded: \"" << infile << "\" Check the file's content." << endl << endl;
        appearErrors = true;
      }
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -out ------------------------------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-out") && (i+1 != argc) && (outfile == NULL)) {
      outfile = new char[strlen(argv[++i]) + 1];
      strcpy(outfile, argv[i]);
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /*                                         Information File format                                        */

   /* Option -format ------------------------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-reverse") && (!reverse))
      reverse = true;
   /* ------------------------------------------------------------------------------------------------------ */

   /*                                         Information File format                                        */

   /* Option -format ------------------------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-format") && (!format))
      format = true;

   /* Option -type ------------------------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-type") && (!type))
      type = true;

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

   /* Option -html ---------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-html") && (outformat == -1))
      outformat = 100;

   /* Option -onlyseqs ------------------------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-onlyseqs") && (outformat == -1))
      outformat = 99;
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

  /***** ***** ***** ***** ***** ***** ***** Postprocessing parameters  ***** ***** ***** ***** ***** ***** *****/

 /* ------------------------------------------------------------------------------------------------------ */

  if((infile == NULL) && (!appearErrors)) {
    cerr << endl << "ERROR: An input file has to be defined." << endl << endl;
    appearErrors = true;
  }

 /* ------------------------------------------------------------------------------------------------------ */

  if((outformat == -1) && (!format) && (!type) && (!reverse) && (!appearErrors)) {
    cerr << endl << "ERROR: An option has to be defined." << endl << endl;
    appearErrors = true;
  }

 /* ------------------------------------------------------------------------------------------------------ */

  if(((outformat != -1) || (reverse)) && ((format) || (type)) && (!appearErrors)) {
    cerr << endl << "ERROR: Only one option has to be choosen: either an output format or the input file information "
         << "(Format and/or type options)." << endl << endl;
    appearErrors = true;
  }

 /* ------------------------------------------------------------------------------------------------------ */

 /* ------------------------------------------------------------------------------------------------------ */

  if(((format) || (type)) && (outfile != NULL) && (!appearErrors)) {
    cerr << endl << "ERROR: Not should defined an output file in order to get input file information." << endl << endl;
    appearErrors = true;
  }

 /* ------------------------------------------------------------------------------------------------------ */

  if(appearErrors) {
    delete[] outfile;
    delete[] infile;

    return -1;
  }

  /* ------------------------------------------------------------------------------------------------------ */

  /***** ***** ***** ***** ***** ***** **** End of Parameters Processing **** ***** ***** ***** ***** ***** *****/

  if((format) || (type)) {
    cout << endl << infile;

    if(format) {
      switch(compAlig.formatInputFile()) {
        case 1:  cout << "\t" << "clustal";           break;
        case 3:  cout << "\t" << "nbrf/pir";          break;
        case 8:  cout << "\t" << "fasta";             break;
        case 11: cout << "\t" << "phylip3.2";         break;
        case 12: cout << "\t" << "phylip";            break;
        case 17: cout << "\t" << "nexus";             break;
        case 21: cout << "\t" << "mega_interl";       break;
        case 22: cout << "\t" << "mega_noninterl";    break;
        default: cout << "\t" << "unknown";
      }
      switch(compAlig.typeInputFile()) {
        case SINGLE: cout << "\t" << "single";        break;
        case MULTI:  cout << "\t" << "multi";         break;
        default:     cout << "\t" << "unknown";
      }
    }

    if(type) {
      switch(compAlig.getTypeAlignment()) {
        case DNAType:
        case RNAType: cout << "\t" << "nucleotides";  break;
        case AAType:  cout << "\t" << "aminoacids";   break;
        default: cout << "\t" << "unknown";
      }
    }
    cout << endl << endl;
  }

  else {

    /* -------------------------------------------------------------------- */
    if(outformat != -1) compAlig.setOutputFormat(outformat, shortNames);
    if(reverse) compAlig.setReverse();
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if(outfile != NULL) {
      if(!compAlig.saveAlignment(outfile)) {
        cerr << endl << "ERROR: It's impossible to generate the output file." << endl << endl;
        return -1;
      }
    }
    else {
      compAlig.printAlignment();
    }
    /* -------------------------------------------------------------------- */
  }

  /* -------------------------------------------------------------------- */
  delete[] outfile;
  delete[] infile;
  /* -------------------------------------------------------------------- */

  return 0;
}

void menu(void) {

  cout << endl;
  cout << "readAl " << VERSION << "rev" << REVISION << ". Copyright (C) 2009. Salvador Capella-Gutierrez and "
       << "Toni Gabaldón." << endl << endl;

  cout << "readAl webpage: http://trimal.cgenomics.org" << endl << endl;

  cout << "This program is free software: you can redistribute it and/or modify " << endl
       << "it under the terms of the GNU General Public License as published by " << endl
       << "the Free Software Foundation, the last available version." << endl << endl;

  cout << "Basic usage" << endl
       << "\treadal -in <inputfile> -out <outputfile> -(format options)." << endl << endl;

  cout << "    -h                       " << "Prints this information and show some examples." << endl;
  cout << "    --version                " << "Prints the readAl version." << endl << endl;

  cout << "    -in <inputfile>          " << "Input file in several formats (clustal, fasta, NBRF/PIR, nexus, phylip3.2, phylip, mega)."
                                          << endl;
  cout << "    -out <outputfile>        " << "Output alignment in the same input format (default stdout)." << endl << endl;

  cout << "    -reverse                 " << "Output the reverse input alignment." << endl << endl;

  cout << "    -format                  " << "Prints the the input file format (fasta, clustal, ...) and kind (single or multi)." << endl;
  cout << "    -type                    " << "Prints the input file type (Nucleotides or Amino Acids)." << endl << endl;

  cout << "    -clustal                 " << "Output file in CLUSTAL format" << endl;
  cout << "    -fasta                   " << "Output file in FASTA format" << endl;
  cout << "    -fasta_m10               " << "Output file in FASTA format. Sequences name length up to 10 characters." << endl;
  cout << "    -nbrf                    " << "Output file in NBRF/PIR format" << endl;
  cout << "    -nexus                   " << "Output file in NEXUS format" << endl;
  cout << "    -mega                    " << "Output file in MEGA format" << endl;
  cout << "    -phylip_paml             " << "Output file in PHYLIP format compatible with PAML" << endl;
  cout << "    -phylip_paml_m10         " << "Output file in PHYLIP format compatible with PAML. Sequences name length up to 10 characters." << endl;
  cout << "    -phylip3.2               " << "Output file in PHYLIP3.2 format" << endl;
  cout << "    -phylip3.2_m10           " << "Output file in PHYLIP3.2 format. Sequences name length up to 10 characters." << endl;
  cout << "    -phylip                  " << "Output file in PHYLIP/PHYLIP4 format" << endl;
  cout << "    -phylip_m10              " << "Output file in PHYLIP/PHYLIP4 format. Sequences name length up to 10 characters." << endl << endl;

  cout << "    -html                    " << "Output file in HTML with residues colored using the CLUSTAL scheme" << endl << endl;

  cout << "    -onlyseqs                " << "Output file with the sequences without gaps in FASTA format" << endl << endl;
}

