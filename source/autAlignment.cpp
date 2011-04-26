/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.3: a tool for automated alignment trimming in large-scale
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
#include "alignment.h"

#define GAPPYOUT 1
#define STRICT 2

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function computes the identities values between the sequences from
 * the alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::calculateSeqIdentity(void) {

  char indet;
  int i, j, k, hit, dst;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Depending on the alignment type, the indetermination
   * symbol will be one or other */
  if(getTypeAlignment() == AAType) indet = 'X';
  else indet = 'N';
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Create the identities matrix to store the identity
   * score among the different sequences in the alig */
  identities = new float*[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* For each sequence, we have to compute its identity
   * value respect of the otherones in the alignment */
  for(i = 0; i < sequenNumber; i++) {
    identities[i] = new float[sequenNumber];

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
	/* It's a symmetric matrix, we copy the values that
	 * have been already computed */
    for(j = 0; j < i; j++)
      identities[i][j] = identities[j][i];
    identities[i][i] = 0;
	/* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
	/* For the rest of the sequences, we compute the
	 * identity value between the selected sequence and
	 * these ones */
    for(j = i + 1; j < sequenNumber; j++) {
      for(k = 0, hit = 0, dst = 0; k < residNumber; k++) {
		/* ***** ***** ***** ***** ***** ***** ***** ***** */
		/* If one of the two positions is a valid residue,
		 * we take it into account for the distance */
        if(((sequences[i][k] != indet) && (sequences[i][k] != '-')) ||
           ((sequences[j][k] != indet) && (sequences[j][k] != '-'))) {
          dst++;
		  /* If both position have the same residue, we
		   * count a hit */
          if(sequences[i][k] == sequences[j][k])
            hit++;
        }
        /* ***** ***** ***** ***** ***** ***** ***** ***** */
      }
	  /* ***** ***** ***** ***** ***** ***** ***** ***** */
	  /* The identity value between two given sequence is
	   * the proportion between number of hits and number
	   * of residue, commons and non-commons, from these
	   * sequences */
      identities[i][j] = (float) hit/dst;
	  /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function computes some parameters from the input alignment such as
 * identity average, identity average from each sequence and its most similar
 * one, etc, to select which one is the best automated method to trim this
 * alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int alignment::selectMethod(void) {

  float mx, avg, maxSeq = 0, avgSeq = 0;
  int i, j;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Ask for the sequence identities assesment */
  if(identities == NULL)
    calculateSeqIdentity();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Once we have the identities among all possible
   * combinations between each pair of sequence. We
   * compute the average identity as well as the
   * average identity for each sequence with its most
   * similar one */
  for(i = 0; i < sequenNumber; i++) {
    for(j = 0, mx = 0, avg = 0; j < sequenNumber; j++) {
      if(i != j) {
        mx  = mx < identities[i][j] ? identities[i][j] : mx;
        avg += identities[i][j];
      }
    }
    avgSeq += avg/(sequenNumber - 1);
    maxSeq += mx;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  avgSeq = avgSeq/sequenNumber;
  maxSeq = maxSeq/sequenNumber;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* With the different parameters, we decide wich one
   * is the best automated method, based on a previous
   * simulated data benchmarks, to trim the alig */
  if(avgSeq >= 0.55)      return GAPPYOUT;
  else if(avgSeq <= 0.38) return STRICT;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Sometimes we need to use more parameters to select
   * the best automated method, always based on our
   * benchmarks, to trim the input alignment */
  else {
    if(sequenNumber <= 20) return GAPPYOUT;
    else {
      if((maxSeq >= 0.5) && (maxSeq <= 0.65)) return GAPPYOUT;
      else return STRICT;
    }
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method prints different identity values computed from the alignment.
 * In this method, we asses the identity values matrix as well as diferent
 * average values. Moreover, the method computes which one is the most
 * similar sequence, in term of identity values, for each one in this alig */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::printSeqIdentity(void) {

  int i, j, k, pos, maxLongName;
  float mx, avg, maxSeq = 0, avgSeq = 0, **maxs;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Ask for the sequence identities assesment */
  if(identities == NULL)
    calculateSeqIdentity();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* For each sequence, we look for its most similar
   * one */
  maxs = new float*[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {
    maxs[i] = new float[2];

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    for(k = 0, mx = 0, avg = 0, pos = i; k < sequenNumber; k++) {
      if(i != k) {
        if(mx < identities[i][k]) {
          mx  = identities[i][k];
          pos = k;
        }
        avg += identities[i][k];
      }
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    avgSeq += avg/(sequenNumber - 1);
    maxSeq += mx;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

	/* ***** ***** ***** ***** ***** ***** ***** ***** */
    maxs[i][0] = mx;
    maxs[i][1] = pos;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  avgSeq = avgSeq/sequenNumber;
  maxSeq = maxSeq/sequenNumber;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Once the method has computed all of different
   * values, it prints it */
  for(i = 0, maxLongName = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  cout.precision(6);
  cout << endl << "#Mean Percentage of identity:"
       << "                             " << avgSeq;
  cout << endl << "#Mean Percentage of identity "
       << "with most similar sequence: " << maxSeq;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  cout << endl << endl << "#Percentage of identity matrix:";
  for(i = 0; i < sequenNumber; i++) {
    cout << endl << setw(maxLongName + 2) << left << seqsName[i] << "  ";
    for(j = 0; j < sequenNumber; j++)
      cout << setiosflags(ios::left) << setw(10) << identities[i][j] * 100 << "  ";
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  cout << endl << endl << "#Percentage of identity with most similar sequence:" << endl;
  for(i = 0; i < sequenNumber; i++)
    cout << setw(maxLongName + 2) << left << seqsName[i] << "  " << setiosflags(ios::left)
	     << setw(5) << maxs[i][0] * 100   << "    " << seqsName[(int) maxs[i][1]] << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  cout << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}
