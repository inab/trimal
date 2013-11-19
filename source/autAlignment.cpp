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
#include "alignment.h"
#include "defines.h"

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function computes the identities values between the sequences from
 * the alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::calculateSeqIdentity(void) {

  int i, j, k, hit, dst;
  char indet;

  /* Depending on alignment type, indetermination symbol will be one or other */
  indet = getTypeAlignment() == AAType ? 'X' : 'N';

  /* Create identities matrix to store identities scores */
  identities = new float*[sequenNumber];

  /* For each seq, compute its identity score against the others in the MSA */
  for(i = 0; i < sequenNumber; i++) {
    identities[i] = new float[sequenNumber];

    /* It's a symmetric matrix, copy values that have been already computed */
    for(j = 0; j < i; j++)
      identities[i][j] = identities[j][i];
    identities[i][i] = 0;

    /* Compute identity scores for the current sequence against the rest */
    for(j = i + 1; j < sequenNumber; j++) {
      for(k = 0, hit = 0, dst = 0; k < residNumber; k++) {
      /* If one of the two positions is a valid residue,
       * count it for the common length */
        if(((sequences[i][k] != indet) && (sequences[i][k] != '-')) ||
           ((sequences[j][k] != indet) && (sequences[j][k] != '-'))) {
          dst++;
          /* If both positions are the same, count a hit */
          if(sequences[i][k] == sequences[j][k])
            hit++;
        }
      }

      /* Identity score between two sequences is the ratio of identical residues
       * by the total length (common and no-common residues) among them */
      identities[i][j] = (float) hit/dst;
    }
  }
}

void alignment::calculateRelaxedSeqIdentity(void) {
  /* Raw approximation of sequence identity computation designed for reducing
   * comparisons for huge alignemnts */

  int i, j, k, hit;

  /* Create identities matrix to store identities scores */
  identities = new float*[sequenNumber];

  /* For each seq, compute its identity score against the others in the MSA */
  for(i = 0; i < sequenNumber; i++) {
    identities[i] = new float[sequenNumber];

    /* It's a symmetric matrix, copy values that have been already computed */
    for(j = 0; j < i; j++)
      identities[i][j] = identities[j][i];
    identities[i][i] = 0;

    /* Compute identity score between the selected sequence and the others */
    for(j = i + 1; j < sequenNumber; j++) {
      for(k = 0, hit = 0; k < residNumber; k++) {
        /* If both positions are the same, count a hit */
        if(sequences[i][k] == sequences[j][k])
          hit++;
      }
    /* Raw identity score is computed as the ratio of identical residues between
     * alignment length */
      identities[i][j] = (float) hit/residNumber;
    }
  }
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
  float mx, avg, maxAvgSeq = 0, maxSeq = 0, avgSeq = 0, **maxs;

  /* Ask for the sequence identities assesment */
  if(identities == NULL)
    calculateSeqIdentity();

  /* For each sequence, we look for its most similar one */
  maxs = new float*[sequenNumber];

  for(i = 0; i < sequenNumber; i++) {
    maxs[i] = new float[2];

    /* Get the most similar sequence to the current one in term of identity */
    for(k = 0, mx = 0, avg = 0, pos = i; k < sequenNumber; k++) {
      if(i != k) {
        avg += identities[i][k];
        if(mx < identities[i][k]) {
          mx = identities[i][k];
          pos = k;
        }
      }
    }
    /* Update global average variables*/
    avgSeq += avg/(sequenNumber - 1);
    maxAvgSeq += mx;

    /* Save the maximum average identity value for each sequence */
    maxs[i][0] = mx;
    maxs[i][1] = pos;
  }

  /* Compute general averages */
  avgSeq = avgSeq/sequenNumber;
  maxAvgSeq = maxAvgSeq/sequenNumber;

  /* Compute longest sequences name */
  for(i = 0, maxLongName = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());

  /* Once the method has computed all of different values, it prints it */
  cout.precision(4);
  cout << fixed;

  for(i = 0, maxSeq = 0; i < sequenNumber; i++)
    if(maxs[i][0] > maxSeq)
      maxSeq = maxs[i][0];

  cout << endl << "## MaxIdentity\t" << maxSeq;
  cout << endl << "#> MaxIdentity\tGet the maximum identity value for any pair "
    << "of sequences in the alignment" << endl;

  cout << endl << "## AverageIdentity\t" << avgSeq;
  cout << endl << "#> AverageIdentity\tAverage identity between all sequences";

  cout << endl << endl << "## Identity sequences matrix";
  for(i = 0; i < sequenNumber; i++) {
    cout << endl << setw(maxLongName + 2) << left << seqsName[i] << "\t";
    for(j = 0; j < i; j++)
      cout << setiosflags(ios::left) << setw(10) << identities[i][j] << "\t";
    cout << setiosflags(ios::left) << setw(10) << 1.00 << "\t";
    for(j = i + 1; j < sequenNumber; j++)
      cout << setiosflags(ios::left) << setw(10) << identities[i][j] << "\t";
  }
  cout << endl;

  cout << endl << "## AverageMostSimilarIdentity\t" << maxAvgSeq;
  cout << endl << "#> AverageMostSimilarIdentity\t Average identity between "
    << "most similar pair-wise sequences";

  cout << endl << endl << "## Identity for most similar pair-wise sequences "
    << "matrix" << endl;
  for(i = 0; i < sequenNumber; i++)
    cout << setw(maxLongName + 2) << left << seqsName[i]
      << "\t" << setiosflags(ios::left) << setw(5)
      << maxs[i][0] << "\t" << seqsName[(int) maxs[i][1]] << endl;
  cout << endl;
}

/* *** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *** */
/*                                                                           */
/*                             NEW CODE: feb/2012                            */
/*                                                                           */
/* *** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *** */
void alignment::calculateColIdentity(float *ColumnIdentities) {

  int i, j, counter, pos, max, columnLen;
  char letter, indet, gapSymbol;
  string column;

  /* Initialize some data for make computation more precise */
  indet = getTypeAlignment() == AAType ? 'X' : 'N';
  gapSymbol = '-';

  /* Compute identity score for the most frequent residue, it can be as well
   * gaps and indeterminations, for each column */
  for(i = 0, max = 0; i < residNumber; i++, max = 0, column.clear()) {

    /* Get residues from each column in capital letters */
    for(j = 0; j < sequenNumber; j++)
      /* Discard gaps and indeterminations from calculations */
      if((toupper(sequences[j][i]) != indet) && (sequences[j][i] != gapSymbol))
        column += toupper(sequences[j][i]);
    columnLen = column.size();

    /* Count letter frequency. It only matter the frequency. Use some shorcuts
     * to speed-up the process */
    while (!column.empty()) {
      letter = column[0];
      counter = 0;
      pos = 0;

      do {
        counter += 1;
        column.erase(pos, 1);
        pos = column.find(letter, pos);
      } while(pos != (int) string::npos);

      /* Keep only the most frequent residue */
      if(counter > max)
        max = counter;
      /* If column size is smaller than the current max, stop the count */
      if((int) column.size() < max)
        break;
    }

    /* Store column identity values */
    if(columnLen != 0)
      ColumnIdentities[i] = float(max)/columnLen;
  }
}

void alignment::printColumnsIdentity_DescriptiveStats(void) {

  float *colIdentities, avg, std, max, min;
  int i, positions;

  /* Allocate local memory for the computation */
  colIdentities = new float[residNumber];

  utils::initlVect(colIdentities, residNumber, -1);
  calculateColIdentity(colIdentities);

  for(i = 0, max = 0, min = 1, avg = 0, positions = 0; i < residNumber; i++) {
    if(colIdentities[i] != -1) {
      /* Compute on-the-fly max and min scores. Store accumulative score */
      avg += colIdentities[i];
      max = (colIdentities[i] > max) ? colIdentities[i] : max;
      min = (colIdentities[i] < min) ? colIdentities[i] : min;
      /* Count how many columns have a value score */
      positions += 1;
    }
  }
  /* Compute average identity column score */
  avg /= positions;

  /* Compute standard desviation */
  for(i = 0, std = 0; i < residNumber; i++)
    if(colIdentities[i] != -1)
      std += pow((colIdentities[i] - avg), 2);
  std = sqrt(std/positions);

  /* Print general descriptive stats */
  cout << "#maxColIdentity\t" << max << endl;
  cout << "#minColIdentity\t" << min << endl;
  cout << "#avgColIdentity\t" << avg << endl;
  cout << "#stdColIdentity\t" << std << endl;
}




