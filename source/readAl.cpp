/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    readAl v1.4: a tool for automated alignment conversion among different
                 formats.

    2009-2013 Capella-Gutierrez S. and Gabaldon, T.
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

***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#include <stdlib.h>
#include <string.h>

#include "alignment.h"
#include "defines.h"
#include "utils.h"

void menu(void);

int main(int argc, char *argv[]) {

  /* Input alignment */
  alignment inAlig;

  /* Local variables */
  string align_format;
  int i, outformat = -1;
  char *infile = NULL, *outfile = NULL;
  bool errors = false, reverse = false, shortNames = false, format = false, \
    type = false, info = false;

  /* If there is no parameters: Inform about readAl options and finish */
  if(argc == 1) {
    menu();
    return 0;
  }

  i = 1;
  /* If option -h has been used, inform about readAl options and finish */
  if(!strcmp(argv[i], "-h") && (i+1 == argc)) {
    menu();
    return 0;
  }

  /* Inform about current readAl version/revision/build and finish */
  if(!strcmp(argv[i], "--version") && (i+1 == argc)) {
      cout << endl << "readAl v" << VERSION << ".rev" << REVISION << " build["
        << BUILD << "]" << endl << endl;
    return 0;
  }

  /* Catch different input options and then check whether there is a valid
   * combination of parameters */
  while(i < argc) {

    /* Input alignment option: -in */
    if(!strcmp(argv[i], "-in") && (i+1 != argc) && (infile == NULL)) {
      /* Allocate memory for storing input alignment filename */
      infile = new char[strlen(argv[++i]) + 1];
      strcpy(infile, argv[i]);

      /* Load input alignment and inform about it if something is wrong */
      if(!inAlig.loadAlignment(infile)) {
        cerr << endl << "ERROR: Alignment not loaded: \"" << infile
          << "\" Check the file's content." << endl << endl;
        errors = true;
      }
    }

    /* Output filename option: -out */
    else if(!strcmp(argv[i], "-out") && (i+1 != argc) && (outfile == NULL)) {
      /* Allocate memory for storing output alignment filename */
      outfile = new char[strlen(argv[++i]) + 1];
      strcpy(outfile, argv[i]);
    }

    /* Get information about input file format */
    else if(!strcmp(argv[i], "-format") && (!format))
      format = true;

    /* Get information about input file residues type */
    else if(!strcmp(argv[i], "-type") && (!type))
      type = true;

    /* Get general information about input file: seqs number, average seq length,
     * etc */
    else if(!strcmp(argv[i], "-info") && (!info))
      info = true;

    /* Get input sequences reverse option: -reverse */
    else if(!strcmp(argv[i], "-reverse") && (!reverse))
      reverse = true;

    /* For all output format options is checked if more
     * than one output format has been required */

    /* Set output alignment format to CLUSTAL: -clustal */
    else if(!strcmp(argv[i], "-clustal") && (outformat == -1))
      outformat = 1;

    /* Set output alignment format to FASTA: -fasta */
    else if(!strcmp(argv[i], "-fasta") && (outformat == -1))
      outformat = 8;

   /* Set output alignment format to FASTA and ask for using only
    * up to 10 characters for sequences name: -fasta_m10 */
    else if(!strcmp(argv[i], "-fasta_m10") && (outformat == -1)) {
      outformat = 8;
      shortNames = true;
    }

    /* Set output alignment format to NBRF/PIR: -nbrf */
    else if(!strcmp(argv[i], "-nbrf") && (outformat == -1))
      outformat = 3;

    /* Set output alignment format to NEXUS: -nexus */
    else if(!strcmp(argv[i], "-nexus") && (outformat == -1))
      outformat = 17;

    /* Set output alignment format to MEGA: -mega */
    else if(!strcmp(argv[i], "-mega") && (outformat == -1))
      outformat = 21;

    /* Set output alignment format to PHYLIP3.2 (sequential): -phylip3.2 */
    else if(!strcmp(argv[i], "-phylip3.2") && (outformat == -1))
      outformat = 11;

    /* Set output alignment format to PHYLIP3.2 (sequential) and ask for
     * using only up to 10 characters for sequences name: -phylip3.2_m10 */
    else if(!strcmp(argv[i], "-phylip3.2_m10") && (outformat == -1)) {
      outformat = 11;
      shortNames = true;
    }

    /* Set output alignment format to PHYLIP (interleaved): -phylip */
    else if(!strcmp(argv[i], "-phylip") && (outformat == -1))
      outformat = 12;

    /* Set output alignment format to PHYLIP (interleaved) and ask for
     * using only up to 10 characters for sequences name: -phylip_m10 */
    else if(!strcmp(argv[i], "-phylip_m10") && (outformat == -1)) {
      outformat = 12; shortNames = true;
    }

    /* Set output alignment format to PHYLIP compatible with programs
     * such as PAML: -phylip_paml */
    else if(!strcmp(argv[i], "-phylip_paml") && (outformat == -1))
      outformat = 13;

    /* Set output alignment format to PHYLIP compatible with programs such as
     * PAML and ask for using only up to 10 characters for sequences name:
     * -phylip_paml_m10 */
    else if(!strcmp(argv[i], "-phylip_paml_m10") && (outformat == -1)) {
      outformat = 13;
      shortNames = true;
    }

    /* Set output alignment format to HTML, that means residues will be colored
     * according to its physic-chemical properties using CLUSTAL color scheme:
     * -html */
    else if(!strcmp(argv[i], "-html") && (outformat == -1))
      outformat = 100;

    /* Get unaligned sequences from input file: -onlyseqs */
    else if(!strcmp(argv[i], "-onlyseqs") && (outformat == -1))
      outformat = 99;

    /* Inform about no valid options */
    else {
      cerr << endl << "ERROR: Parameter \"" << argv[i] << "\" not valid."
        << endl << endl;
      errors = true;
    }
    i++;

    /* If any error has been detected, break input options loop
     * and then process detected error */
    if(errors)
      break;
  }

  /* Final verifications to detect any possible mistake in the input options */
  /* It is mandatory to provide an input file. Otherwise, inform about it */
  if((infile == NULL) && (!errors)) {
    cerr << endl << "ERROR: An input file has to be defined." << endl << endl;
    errors = true;
  }

  /* It is mandatory to choose an option for processing input alignment */
  if((outformat == -1) && (!reverse) && (!format) && (!type) && (!info)
    && (!errors)) {
    cerr << endl << "ERROR: An option has to be chosen." << endl << endl;
    errors = true;
  }

  /* Only one option can be selected when an output file is not defined */
  if((outfile == NULL) && ((outformat != -1) || reverse) && (format || type \
    || info) && (!errors)) {
    cerr << endl << "ERROR: Only one option can be selected: either an output "
      << "format or get information about input file when an output file is "
      << "not defined" << endl << endl;
    errors = true;
  }

  /* Does not make any sense to define any output file when
   * only information about input alignment is requested */
  if(((outfile != NULL) && outformat == -1 && !reverse) && (format || type \
    || info) && (!errors)) {
    cerr << endl << "ERROR: An output file should not be provided when only "
      << "information about input alignment is requested" << endl << endl;
    errors = true;
  }

  /* If no error has been detected, process input file */
  if(!errors) {

    /* Print information about input alignment */
    if((format) || (type) || (info)) {
      cout << "## Input filename\t'" << infile << "'" << endl;

      if(format) {
        /* Input file format */
        if (inAlig.formatInputFile() == 1)
          align_format = "clustal";
        else if (inAlig.formatInputFile() == 3)
          align_format = "nbrf/pir";
        else if (inAlig.formatInputFile() == 8)
          align_format = "fasta";
        else if (inAlig.formatInputFile() == 11)
          align_format = "phylip3.2";
        else if (inAlig.formatInputFile() == 12)
          align_format = "phylip";
        else if (inAlig.formatInputFile() == 17)
          align_format = "nexus";
        else if (inAlig.formatInputFile() == 21)
          align_format = "mega_interleaved";
        else if (inAlig.formatInputFile() == 22)
          align_format = "mega_sequential";
        else
          align_format = "unknown";

        /* Inform about if sequences are aligned or not */
        cout << "## Input file format\t" << align_format << endl
          << "## Input file aligned\t" << (inAlig.isFileAligned() ? "YES":"NO")
          << endl;
      }

      if(type) {
        /* Inform about biological datatype */
        if (inAlig.getTypeAlignment() == DNAType)
          cout << "## Input file datatype\tnucleotides:dna" << endl;
        else if (inAlig.getTypeAlignment() == DNADeg)
          cout << "## Input file datatype\tnucleotides:dna_degenerate_codes"
            << endl;
        else if (inAlig.getTypeAlignment() == RNAType)
          cout << "## Input file datatype\tnucleotides:rna" << endl;
        else if (inAlig.getTypeAlignment() == RNADeg)
          cout << "## Input file datatype\tnucleotides:rna_degenerate_codes"
            << endl;
        else if (inAlig.getTypeAlignment() == AAType)
          cout << "## Input file datatype\tamino-acids" << endl;
        else
          cout << "## Input file datatype\tunknown" << endl;
      }

      if(info)
        inAlig.printAlignmentInfo(cout);
    }

    if((outfile != NULL) || (outformat != -1) || reverse || shortNames) {
      /* Set output format */
      if(outformat != -1 || shortNames)
        inAlig.setOutputFormat(outformat, shortNames);
      /* Ask for getting the reverse of input file */
      if(reverse)
        inAlig.setReverse();

      /* If a outfile has been provided, try to generate output file */
      if(outfile != NULL) {
        if(!inAlig.saveAlignment(outfile)) {
          cerr << endl << "ERROR: Impossible to generate OUTPUT file." << endl
            << endl;
          return -1;
        }
      /* ... otherwise dump outfile content to standard output */
      } else {
        inAlig.printAlignment();
      }
    }
  }

  /* Deallocate local memory */
  delete [] infile;
  delete [] outfile;

  /* Inform about readAl execution */
  return (errors == true ? -1 : 0);
}

void menu(void) {

  cout << endl
    << "readAl v" << VERSION << ".rev" << REVISION << " build[" << BUILD
    << "]. " << AUTHORS << endl << endl

    << "readAl webpage: http://trimal.cgenomics.org" << endl << endl

    << "This program is free software: you can redistribute it and/or modify "
    << endl
    << "it under the terms of the GNU General Public License as published by "
    << endl
    << "the Free Software Foundation, the last available version." << endl
    << endl

    << "Basic usage" << endl
    << "\treadal -in <inputfile> -out <outputfile> [options]." << endl << endl

    << "\t-h                   " << "Show this information." << endl
    << "\t--version            " << "Show readAl version." << endl << endl

    << "\t-in <inputfile>      " << "Input file in several formats." << endl
    << "\t-out <outputfile>    " << "Output file name (default STDOUT)." << endl
    << endl

    << "\t-format              " << "Print information about input file format "
    << "and if sequences are aligned or not." << endl

    << "\t-type                " << "Print information about biological "
    << "sequences datatype (e.g. nucleotides:dna, nucleotides:rna, aminoacids, etc)"
    << endl

    << "\t-info                " << "Print information about sequences number, "
    << "average sequence length, max & min sequence length"
    << endl << endl

    << "\t-onlyseqs            " << "Generate output with only residues from "
    << "input file" << endl << endl

    << "\t-html                " << "Output residues colored according their "
    << "physicochemical properties. HTML file." << endl << endl

    << "\t-reverse             " << "Output the reverse of sequences in "
    << "input file." << endl << endl

    << "\t-nbrf                " << "Output file in NBRF/PIR format" << endl
    << "\t-mega                " << "Output file in MEGA format" << endl

    << "\t-nexus               " << "Output file in NEXUS format" << endl
    << "\t-clustal             " << "Output file in CLUSTAL format" << endl
    << endl

    << "\t-fasta               " << "Output file in FASTA format" << endl
    << "\t-fasta_m10           " << "Output file in FASTA format. Sequences "
    << "name up to 10 characters." << endl << endl

    << "\t-phylip              " << "Output file in PHYLIP/PHYLIP4 format"
    << endl
    << "\t-phylip_m10          " << "Output file in PHYLIP/PHYLIP4 format. "
    << "Sequences name up to 10 characters." << endl
    << "\t-phylip_paml         " << "Output file in PHYLIP format compatible "
    << "with PAML" << endl
    << "\t-phylip_paml_m10     " << "Output file in PHYLIP format compatible "
    << "with PAML. Sequences name up to 10 characters." << endl
    << "\t-phylip3.2           " << "Output file in PHYLIP3.2 format" << endl
    << "\t-phylip3.2_m10       " << "Output file in PHYLIP3.2 format. Sequences"
    << " name up to 10 characters." << endl << endl;
}
