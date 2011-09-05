/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    2009-2011 Capella-Gutierrez S. and Gabaldon, T.
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

#include "compareFiles.h"
#include "alignment.h"

#define LONG 80

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method compares a set of alignment in order to select the most
 * consistent one respect of the other ones. To compute the consistency
 * values we use the proportion of residue pairs per column in the aligs
 * to compare */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int compareFiles::algorithm(alignment **vectAlignments, char **fileNames, float *columnsValue, int numAlignments, bool verbosity) {

  int *numResiduesAlig, *correspNames, *columnSeqMatrix, *columnSeqMatrixAux;
  int i, j, k, l, m, numSeqs, pairRes, hits, alig = 0;
  float max = 0, value = 0, **vectHits;
  bool appearErrors = false;
  string *names;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Get some parameters from the alignment that has
   * been selected */
  numSeqs = vectAlignments[0] -> getNumSpecies();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate dinamic local memory */
  names = new string[numSeqs];
  correspNames = new int[numSeqs];
  numResiduesAlig = new int[numSeqs];
  columnSeqMatrix = new int[numSeqs];
  vectHits = new float*[numAlignments];
  columnSeqMatrixAux = new int[numSeqs];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check that all of alignment has the same number of
   * sequence as well as there exists a correspondence
   * between the names for each pars of aligs. */
  for(i = 1; i < numAlignments; i++) {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if(numSeqs != vectAlignments[i] -> getNumSpecies()) {
      cerr << endl << "ERROR: The files to compare do not have "
           << "the same number of sequences" << endl << endl;
      appearErrors = true;
      break;
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    vectAlignments[i] -> getSequences(names);
    if(!vectAlignments[0] -> getSeqNameOrder(names, correspNames)) {
      cerr << endl << "ERROR: The files to compare do not"
           << " have the sequence names" << endl << endl;
      appearErrors = true;
      break;
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Changes the order in sequences number matrix
   * according to the order in the selected alignment */
  for(i = 1; ((i < numAlignments) && (!appearErrors)); i++) {
    vectAlignments[i] -> getSequences(names);
    vectAlignments[0] -> getSeqNameOrder(names, correspNames);
    vectAlignments[i] -> setSeqMatrixOrder(correspNames);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Get back the residues number for each alignment */
  for(i = 0; ((i < numAlignments) && (!appearErrors)); i++)
    numResiduesAlig[i] =  vectAlignments[i] -> getNumAminos();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Start the comparison among the alignments */
  for(i = 0; ((i < numAlignments) && (!appearErrors)); i++, value = 0) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If it's necessary, we print some information */
    if(verbosity)
      cout << endl;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Initialize the hits vector for each alignment */
    vectHits[i] = new float[numResiduesAlig[i]];
    utils::initlVect(vectHits[i], numResiduesAlig[i], 0);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    for(j = 0, pairRes = 0, hits = 0; j < numResiduesAlig[i]; j++, pairRes = 0, hits = 0) {

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* Get back each column from the current selected
       * alignment */
      vectAlignments[i] -> getColumnSeqMatrix(j, columnSeqMatrix);
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* For each position from the previous recovered
       * columns, we carry out the analysis to detect the
       * same residues pair in the rest of the alignmetns */
      for(k = 0; k < numSeqs; k++) {

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* If there is a valid residue, we go ahead with
         * the analysis */
        if(columnSeqMatrix[k] != 0) {
          /* ***** ***** ***** ***** ***** ***** ***** ***** */
          for(l = 0; l < i; l++) {
			/* Recover the residue pairs from the others aligs */
            vectAlignments[l] -> getColumnSeqMatrix(columnSeqMatrix[k], k, columnSeqMatrixAux);
            /* ***** ***** ***** ***** ***** ***** ***** ***** */
			/* and count the similar residue pairs */
            for(m = k + 1; m < numSeqs; m++)
              if(columnSeqMatrix[m] != 0) {
                if(columnSeqMatrix[m] == columnSeqMatrixAux[m])
                  hits++;
                pairRes++;
              }
            }
          /* ***** ***** ***** ***** ***** ***** ***** ***** */

          /* ***** ***** ***** ***** ***** ***** ***** ***** */
          for(l = i + 1; l < numAlignments; l++) {
			/* Recover the residue pairs from the others aligs */
            vectAlignments[l] -> getColumnSeqMatrix(columnSeqMatrix[k], k, columnSeqMatrixAux);
            /* ***** ***** ***** ***** ***** ***** ***** ***** */
			/* and count the similar residue pairs */
            for(m = k + 1; m < numSeqs; m++)
              if(columnSeqMatrix[m] != 0) {
                if(columnSeqMatrix[m] == columnSeqMatrixAux[m])
                  hits++;
                pairRes++;
              }
            }
          /* ***** ***** ***** ***** ***** ***** ***** ***** */
        }
      }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
	  /* For each column, compute the hits proportion for
	   * every residue pair against the rest of alignments */
      if(pairRes != 0) {
         vectHits[i][j] += ((1.0 * hits)/pairRes);
         value += vectHits[i][j];
      }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
	/* The method can offer some information about the
	 * comparison progression */
    if(verbosity)
      cout << "File:\t\t" << fileNames[i] << endl << "Values:\t\tSequences: " << numSeqs
           << "\tResidues: " << numResiduesAlig[i] <<  "\tPond. Hits: " << setw(8)
           << value << "\t%Consistency: " << value/numResiduesAlig[i] << endl;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
	/* Keep the alignment with higher consistency value */
    if((value/numResiduesAlig[i]) > max) {
      alig = i;
      max = value/numResiduesAlig[i];
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Prints the alignment that have been selected */
  if((verbosity) && (!appearErrors)) {
    cout << "\t\t\t\t\t--------------" << endl;
    cout << endl << "File Selected:\t" << fileNames[alig] << endl << "Value:\t\t" << max << endl << endl;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* The method returns a vector with the consistency
   * value for each column in the selected alignment */
  if((columnsValue != NULL) && (!appearErrors)) {
    utils::initlVect(columnsValue, numResiduesAlig[alig], -1);
    for(i = 0; i < numResiduesAlig[alig]; i++)
      columnsValue[i] = vectHits[alig][i];
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate memmory */
  for(i = 0; ((i < numAlignments) && (!appearErrors)); i++)
    delete [] vectHits[i];
  delete [] vectHits;

  delete [] names;
  delete [] correspNames;
  delete [] numResiduesAlig;
  delete [] columnSeqMatrix;
  delete [] columnSeqMatrixAux;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Return the selected alignment index or an error
   * flag otherwise */
  if(appearErrors) return -1;
  else return alig;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method returns the consistency value vector for a given alignment
 * against a set of alignments with the same sequences */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool compareFiles::forceComparison(alignment **vectAlignments, int numAlignments, alignment *selected, float *columnsValue) {

  int *correspNames, *columnSeqMatrix, *columnSeqMatrixAux;
  int i, j, k, ll, numResidues, numSeqs, pairRes, hit;
  bool appearErrors = false;
  string *names;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Get some parameters from the alignment that has
   * been selected */
  numResidues  = selected -> getNumAminos();
  numSeqs      = selected -> getNumSpecies();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Initialize the vector where we are going to store
   * the proportion of hits for each column in the
   * selected alignment */
  utils::initlVect(columnsValue, numResidues, 0);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate dinamic local memory */
  names = new string[numSeqs];
  correspNames = new int[numSeqs];
  columnSeqMatrix = new int[numSeqs];
  columnSeqMatrixAux = new int[numSeqs];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check that all of alignment has the same number of
   * sequence as well as there exists a correspondence
   * between the names for each pars of aligs. */
  for(i = 0; i < numAlignments; i++) {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if(numSeqs != vectAlignments[i] -> getNumSpecies()) {
      cerr << endl << "ERROR: The files to compare do not have "
           << "the same number of sequences" << endl << endl;
      appearErrors = true;
      break;
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    vectAlignments[i] -> getSequences(names);
    if(!selected -> getSeqNameOrder(names, correspNames)) {
      cerr << endl << "ERROR: The files to compare do not"
           << " have the sequence names" << endl << endl;
      appearErrors = true;
      break;
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Changes the order in sequences number matrix
   * according to the order in the selected alignment */
  for(i = 0; i < numAlignments; i++) {
    vectAlignments[i] -> getSequences(names);
    selected -> getSeqNameOrder(names, correspNames);
    vectAlignments[i] -> setSeqMatrixOrder(correspNames);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Do the same analysis for each column */
  for(i = 0, pairRes = 0, hit = 0; ((i < numResidues) && (!appearErrors)); i++, pairRes = 0, hit = 0) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* We get back the sequence position for each residue
     * from every column in the selected alignment */
    utils::initlVect(columnSeqMatrix, numSeqs, 0);
    selected -> getColumnSeqMatrix(i, columnSeqMatrix);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* For each residue pairs, we look for it in the rest
     * of alignments */
    for(j = 0; j < numSeqs; j++) {

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* If there is a residue, not a gap, we carry out the
       * assesment */
      if(columnSeqMatrix[j] != 0) {
        for(k = 0; k < numAlignments; k++) {

          /* ***** ***** ***** ***** ***** ***** ***** ***** */
          /* We look for the same residue in the same row in
           * the rest of alignments */
          utils::initlVect(columnSeqMatrixAux, numSeqs, 0);
          vectAlignments[k] -> getColumnSeqMatrix(columnSeqMatrix[j], j, columnSeqMatrixAux);
          /* ***** ***** ***** ***** ***** ***** ***** ***** */

          /* ***** ***** ***** ***** ***** ***** ***** ***** */
          /* We count when we get the same residue pairs in the
           * rest of alignments */
          for(ll = j + 1; ll < numSeqs; ll++)
            if(columnSeqMatrix[ll] != 0) {
              if(columnSeqMatrix[ll] == columnSeqMatrixAux[ll])
                hit++;
              pairRes++;
            }
          /* ***** ***** ***** ***** ***** ***** ***** ***** */
        }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Store the hits proportion for each column */
    if(pairRes != 0) columnsValue[i] += ((1.0 * hit)/pairRes);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate dinamic memory */
  delete [] names;
  delete [] correspNames;
  delete [] columnSeqMatrix;
  delete [] columnSeqMatrixAux;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If it was fine, return true. Otherwise, return
   * false */
  if(appearErrors) return false;
  else return true;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method applies a specific windows size to a selected alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool compareFiles::applyWindow(int columns, int halfWindow, float *columnsValue) {

  int i, j, window;
  float *vectAux;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If windows size is greater than 1/4 of alignment
   *length, trimAl rejects this windows size */
  if(halfWindow > columns/4) return false;
  else window = (2 * halfWindow + 1);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate local memory. Copy the array values to
   * auxiliar memory. */
  vectAux = new float[columns];
  utils::copyVect(columnsValue, vectAux, columns);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* For each column from the selected alignment,
   * compute the average for its consistency values */
  for(i = 0; i < columns; i++) {
	/* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* This average is computed from halfWindow positions
	 * before to halfWindow positions after */
    for(j = i - halfWindow, columnsValue[i] = 0; j <= i + halfWindow; j++) {
      if(j < 0) columnsValue[i] += vectAux[-j];
      else if(j >= columns)
        columnsValue[i] += vectAux[((2 * columns - j) - 2)];
      else columnsValue[i] += vectAux[j];
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Finally, the column value is divided by the window
	 * size in order to compute the average score. */
    columnsValue[i] /= window;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate dinamic memory */
  delete [] vectAux;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If everything is OK, return true */
  return true;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Print the consistency value for each column from the selected alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void compareFiles::printStatisticsFileColumns(int numAminos, float *compareVect) {

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Prepare the header information */
  cout << "| Residue\tConsistency |" << endl;
  cout << "| Number \t   Value    |" << endl;
  cout << "+---------------------------+" << endl;
  cout.precision(10);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Print the consistency values for each column from
   * the selected alignment */
  for(int i = 0; i < numAminos; i++)
    cout << "  " << setw(5) << i + 1 << "\t"
	     << "\t" << compareVect[i] << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Print the consistency values accumulative distribution for the selected
 * alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void compareFiles::printStatisticsFileAcl(int numAminos, float *compareVect) {

  float refer, *vectAux;
  int i, num;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate dinamic memory to copy the input vector
   * and sort it */
  vectAux = new float[numAminos];
  utils::copyVect(compareVect, vectAux, numAminos);
  utils::quicksort(vectAux, 0, numAminos-1);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Set the output precision and print the header */
  cout << "| Number of\t        \t|\t Cumulative \t% "
       << "Cumulative\t|   Consistency   |" << endl;
  cout << "| Residues \t% Length\t|\tNumberResid.\t   "
       << "Length   \t|      Value      |" << endl;
  cout << "+-------------------------------+------------"
       << "---------------------------+-----------------+"
	   << endl;
  cout.precision(10);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Fix the initial values to count how many columns
   * has the same consistency value */
  refer = vectAux[0];
  num = 1;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Print the accumulative distribution */
  for(i = 1; i < numAminos; i++) {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
	/* When the method detects a new consistency value
	 * print the previous value as well as its frequency
	 * and starts to count how many columns are for this
	 * new value */
    if(refer != vectAux[i]) {
      cout << "  " << num << "\t\t" << setw(10) << ((float) num/numAminos * 100.0)
           << "\t\t" << i << "\t\t" << setw(10) << ((float) i/numAminos * 100.0)
            << "\t" << setw(15) << refer << endl;
      refer = vectAux[i];
      num = 1;
    }
    else num++;
	/* ***** ***** ***** ***** ***** ***** ***** ***** */
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Print the last consistency value as well as its
   * frequency */
  cout << "  " << num << "\t\t" << setw(10) << ((float) num/numAminos * 100.0)
       << "\t\t" << i << "\t\t" << setw(10) << ((float) i/numAminos * 100.0)
	   << "\t"  << setw(15) << refer << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate dinamic memory */
  delete [] vectAux;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}
