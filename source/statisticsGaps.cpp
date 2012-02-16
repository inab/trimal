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

#include "statisticsGaps.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  statisticsGaps::statisticsGaps(char **, int, int)                                                                   |
|                                                                                                                      |
|       Class constructor. This method uses the inputs parameters to put the information in the new object that        |
|       has been created.                                                                                              |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

statisticsGaps::statisticsGaps(string *alignmentMatrix, int species, int aminos, int dataType_) {

  int i, j;
  char indet;

  columnLength = species;
  columns =      aminos;
  maxGaps =      0;
  halfWindow =   0;
  dataType = dataType_;

  if(dataType == AAType)
    indet = 'X';
  else
    indet = 'N';

  /* Memory allocation for the vectors and its initialization */
  gapsInColumn =       new int[columns];
  utils::initlVect(gapsInColumn, columns, 0);

  aminosXInColumn =    new int[columns];
  utils::initlVect(aminosXInColumn, aminos, 0);

  gapsWindow =         new int[columns];
  utils::initlVect(gapsWindow, columns, 0);

  numColumnsWithGaps = new int[species+1];
  utils::initlVect(numColumnsWithGaps, columnLength+1, 0);

  /* Count the gaps and indeterminations of each columns */
  for(i = 0; i < columns; i++) {
    for(j = 0; j < columnLength; j++) {
      if(alignmentMatrix[j][i] == '-')
        gapsInColumn[i]++;
      else if(alignmentMatrix[j][i] == indet)
        aminosXInColumn[i]++;
    }

    /* Increase the number of colums with the number of gaps of the last processed column */
    numColumnsWithGaps[gapsInColumn[i]]++;
    gapsWindow[i] = gapsInColumn[i];
    if(gapsWindow[i] > maxGaps) maxGaps = gapsWindow[i];
  }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  statisticsGaps::statisticsGaps(void)                                                                                |
|                                                                                                                      |
|       Class constructor.                                                                                             |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

statisticsGaps::statisticsGaps(void) {

  /* Initializate all values to NULL or 0 */
  gapsInColumn = NULL;
  numColumnsWithGaps = NULL;
  aminosXInColumn = NULL;
  gapsWindow = NULL;

  columns = 	  0;
  columnLength  = 0;
  maxGaps = 	  0;
  halfWindow = 	  0;
  dataType =      0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  statisticsGaps::~statisticsGaps(void)                                                                               |
|                                                                                                                      |
|       Class destroyer.                                                                                               |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

statisticsGaps::~statisticsGaps(void) {

  /* Only free memory if there is previous memory allocation */
  if(gapsInColumn != NULL){
    delete[] gapsInColumn;
    delete[] numColumnsWithGaps;
    delete[] aminosXInColumn;
    delete[] gapsWindow;
  }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool statisticsGaps::applyWindow(int)                                                                               |
|                                                                                                                      |
|       This method computes for each column's alignment its gapwindows' value. For this purpose, the method uses the  |
|       values that previously has been calculated and the window's size value.                                        |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool statisticsGaps::applyWindow(int _halfWindow) {

  int i, j, window;

  /* If one of this conditions is true, we return FALSE:                         */
  /*    .- If already exists a previously calculated vector for this window size */
  /*    .- If halfWinSize value is greater than 1/4 of alignment length        */
  if((halfWindow == _halfWindow) || (_halfWindow > columns/4))
     return false;

  /* Initializate to 0 the vector that will store the number of gaps of each column */
  /* and the vector that will store window processing results                       */
  utils::initlVect(numColumnsWithGaps, columnLength+1, 0);
  utils::initlVect(gapsWindow, columns, 0);

  /* Initializate maximum gaps' number per column value and store the mediumWinSize value in the object. */
  maxGaps = 0;
  halfWindow = _halfWindow;
  window = (2 * halfWindow + 1);

  /* We calculate some statistics for every column in the alignment,and the maximum gaps' number value */
  for(i = 0; i < columns; i++) {
    /* Sum the total number of gaps for the considered window */
    for(j = i - halfWindow, gapsWindow[i] = 0; j <= i + halfWindow; j++) {
      if(j < 0)
        gapsWindow[i] += gapsInColumn[-j];
      else if(j >= columns)
        gapsWindow[i] += gapsInColumn[((2 * columns - j) - 2)];
      else
        gapsWindow[i] += gapsInColumn[j];
    }

    /* Calculate, and round to the nearest integer, the number of gaps for the i column */
    gapsWindow[i] = utils::roundInt(((double) gapsWindow[i]/window));
    /* Increase in 1 the number of colums with the same number of gaps than column i */
    numColumnsWithGaps[gapsWindow[i]]++;

    /* Update the max number of gaps in the alignment, if neccesary */
    if(gapsWindow[i] > maxGaps)
      maxGaps = gapsWindow[i];
  }
  return true;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  int *statisticsGaps::getGapsWindow(void)                                                                            |
|                                                                                                                      |
|       This method returns a pointer to gaps window's vector                                                          |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int *statisticsGaps::getGapsWindow(void) {

  return gapsWindow;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  double statisticsGaps::calcCutPoint(float gapBaseLine, float gapThreshold)                                          |
|                                                                                                                      |
|       This method computes and selects the cut point between two inputs parameters. To knows, the gapBaseline that   |
|       is the minimum columns' percentage to conserve in the new alignment and the gapThreshold that is the maximum   |
|       gaps' percentage that allow us in each alignment's column.                                                     |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double statisticsGaps::calcCutPoint(float gapBaseLine, float gapThreshold) {

  double cutThr, cutBL;
  int i, acum;

  /* We calculate the number of gaps represented by the gaps' percentage. This gaps' number will be the maximum gaps'
     number permitted in the clean alignment. */
  cutThr = (double) columnLength * gapThreshold;

  /* Now, we calculate the number of columns that allow us conserve the columns' percentage fixed by gapBaseline. */
  cutBL = utils::roundInt(((double) (columns * gapBaseLine) / 100.0));

  /* We look the number of gaps that allow us conserve the columns' percentage given by gapBaseLine */
  for(i = 0, acum = 0; i < columnLength; i++) {
    acum += numColumnsWithGaps[i];
    if(acum >= cutBL) break;
  }

  /* If we don't have a exact value for the next of gaps, we compute this value as the gaps' number minus the
     proportion necessary to achieve the columns' number fixed. */
  if(numColumnsWithGaps[i])
    cutBL = (double) (i - ((float) (acum - cutBL)/numColumnsWithGaps[i]));
  else
    cutBL = 0;

  /* Return the maximum value of the two calculated cuts (gapThreshold cut and gapBaseLine cut) */
  return (utils::max(cutThr, cutBL));
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  int statisticsGaps::calcCutPointMixSlope(void)                                                                      |
|                                                                                                                      |
|       This method computes and selects the cut point based on the maximum rate between the first slope ratio between |
|       gaps' percentage in the columns and alignments' length and the "second" slope (slope between three consecutive |
|       points using only the first and third of them) ratio between gaps' percentage in the columns and alignments'   |
|       length.                                                                                                        |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int statisticsGaps::calcCutPointMixSlope(void) {

  float delta = 0, maxSlope = -1, *firstSlopeVector, *secondSlopeVector;
  int prev, pprev, maxIter, row = 1, act = 0, max = 0;

  /* We build two slope vectors, one vector for the first slope and another one for the second. */
  firstSlopeVector = new float[maxGaps+1];
  secondSlopeVector = new float[maxGaps+1];

  /* We initialize them with -1.0 value and fix the maximum iteractions' number as maximun gaps' number plus 1. */
  utils::initlVect(firstSlopeVector, maxGaps, -1.0);
  utils::initlVect(secondSlopeVector, maxGaps, -1.0);
  maxIter = maxGaps + 1;

  /* Until to achieve the maximum iteractions' number. */
  while(act < maxIter) {

    /* We look for a first point to second slope. */
    while((numColumnsWithGaps[act]) == 0) act++;
    pprev = act; if((act+1) >= maxIter) break;

    /* We look for a first point to first slope. */
    do { act++; } while((numColumnsWithGaps[act]) == 0);
    prev = act; if((act+1) >= maxIter) break;

    /* We look for a second point to first and second slope. */
    do { act++; } while((numColumnsWithGaps[act]) == 0);
    if(act >= maxIter) break;

    /* Calculate the first slope between the earlier previus and previus points. */
    firstSlopeVector[prev] =  ((float) (prev - pprev) / columnLength);
    firstSlopeVector[prev] /= ((float) numColumnsWithGaps[prev] / columns);

    /* Calculate the first slope between the previus and current points. */
    firstSlopeVector[act] =  ((float) (act - prev) / columnLength);
    firstSlopeVector[act] /= ((float) numColumnsWithGaps[act] / columns);

    /* Calculate the second slope between the earlier previus and current points. */
    secondSlopeVector[act] = ((float) (act - pprev) / columnLength);
    secondSlopeVector[act] /= ((float) (numColumnsWithGaps[act] + numColumnsWithGaps[prev]) / columns);

    /* If the ratio between ... */
    if((secondSlopeVector[pprev] != -1.0) || (firstSlopeVector[pprev] != -1.0)) {

      /* .- first slope previus and first slope earlier previus points. */
      if(firstSlopeVector[pprev] != -1.0) {
        delta = firstSlopeVector[prev]/firstSlopeVector[pprev];
        row = pprev;
      }

      /* .- first slope current and first slope previus points. */
      if(delta < (firstSlopeVector[act]/firstSlopeVector[prev])) {
        delta = firstSlopeVector[act]/firstSlopeVector[prev];
        row = prev;
      }

      /* .- second slope current and second slope earlier previus points. */
      if(secondSlopeVector[pprev] != -1.0) {
        if(delta < (secondSlopeVector[act]/secondSlopeVector[pprev])) {
          delta = secondSlopeVector[act]/secondSlopeVector[pprev];
          row = pprev;
        }
      }

      /* ... is better that current maxSlope then we modify the maxSlope with the best ratio. */
      if(delta > maxSlope) {
        maxSlope = delta;
        max = row;
      }
    }
    act = prev;
  }

  /* We delete the local memory. */
  delete[] firstSlopeVector;
  delete[] secondSlopeVector;

  /* and, finally, we return the cut point. */
  return max;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  int statisticsGaps::calcCutPoint2ndSlope(void)                                                                      |
|                                                                                                                      |
|       This method computes and selects the cut point based on the maximum "second" slope (slope between three        |
|       consecutive points using only the first and third of them) ratio between gaps' percentage in the columns and   |
|       alignments' length.                                                                                            |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int statisticsGaps::calcCutPoint2ndSlope(void) {

  float maxSlope = -1, *secondSlopeVector;
  int prev, pprev, maxIter, act = 0, max = 0;

  /* We build one slope vector and fix the maximum iteractions' number as the gaps'number plus 1. */
  secondSlopeVector = new float[maxGaps+1];
  utils::initlVect(secondSlopeVector, maxGaps, -1.0);
  maxIter = maxGaps + 1;

  /* Find the lowest number of gaps into the input alignment. If there are few
   * points, it is possible that lowest number of gaps is returned as the thres
   * hold. It could happen input alignment does not have columns with no-gaps */
  for(act = 0, max = 0; numColumnsWithGaps[act] == 0; act++)
    max = act + 1;

  act = 0;
  while(act < maxIter) {

    /* We look for a first point to second slope. */
    while((numColumnsWithGaps[act]) == 0)
      act++;
    pprev = act;
    if((act+1) >= maxIter)
      break;

    /* We look for a first point to first slope. */
    do {
      act++;
    } while((numColumnsWithGaps[act]) == 0);
    prev = act;
    if((act+1) >= maxIter)
      break;

    /* We look for a second point to first and second slope. */
    do {
      act++;
    } while((numColumnsWithGaps[act]) == 0);
    if(act >= maxIter)
      break;

    /* Calculate the second slope between the earlier previous and current points. */
    secondSlopeVector[act] = ((float) (act - pprev) / columnLength);
    secondSlopeVector[act] /= ((float) (numColumnsWithGaps[act] + numColumnsWithGaps[prev]) / columns);

    /* If the ratio between second slope current and second slope earlier previous points. */
    if(secondSlopeVector[pprev] != -1.0) {
      if((secondSlopeVector[act]/secondSlopeVector[pprev]) > maxSlope) {
        maxSlope = (secondSlopeVector[act]/secondSlopeVector[pprev]);
        max = pprev;
      }
    } else if(secondSlopeVector[prev] != -1.0) {
      if((secondSlopeVector[act]/secondSlopeVector[prev]) > maxSlope) {
        maxSlope = (secondSlopeVector[act]/secondSlopeVector[prev]);
        max = pprev;
      }
    }
    act = prev;
  }

  /* We deallocate local memory. */
  delete[] secondSlopeVector;

  /* Finally, we return the selected cut point. */
  return max;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void statisticsGaps::printGapsColumns(void)                                                                         |
|                                                                                                                      |
|       This method shows the gaps' percentage per each column in the alignment.                                       |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void statisticsGaps::printGapsColumns(void) {

  int *vectAux;

  /* We allocate a local vector to recovery information on it */
  vectAux = new int[columns];

  /* We decide about the information's source then we get the information. */
  if(halfWindow == 0)
    utils::copyVect(gapsInColumn, vectAux, columns);
  else
   utils::copyVect(gapsWindow, vectAux, columns);

  /* Fix the precision of output */
  /* We set the output precision and print the header. */
  cout << "| Residue\t  % Gaps \t   Gap Score   |" << endl;
  cout << "| Number \t         \t               |" << endl;
  cout << "+----------------------------------------------+" << endl;
  cout.precision(10);

  /* Show the information that have been requered */
  for(int i = 0; i < columns; i++)
    cout << "  " << setw(5) << i << "\t\t" << setw(10) << (vectAux[i] * 100.0)/columnLength
         << "\t" << setw(7) << 1 -((vectAux[i] * 1.0)/columnLength) << endl;

  /* Finally, we deallocate the local memory */
  delete[] vectAux;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void statisticsGaps::printGapsAcl(void)                                                                             |
|                                                                                                                      |
|       This method shows the gaps' statistics for the alignment.                                                      |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void statisticsGaps::printGapsAcl(void) {

  int acm, i;

  /* Fix the precision of output */
  cout << "| Number of\t        \t|\t Cumulative \t% Cumulative\t|\tNumber of Gaps\t  % Gaps  \tGap Score  |"  << endl;
  cout << "| Residues \t% Length\t|\tNumberResid.\t   Length   \t|\t  per Column  \tper Column\tper Column |" << endl;
  cout << "+-------------------------------+-----------------------------"
       << "----------+--------------------------------------------------+" << endl;
  cout.precision(10);

  /* Count for each gaps' number the columns' number with that gaps' number. */
  for(i = 0, acm = 0; i <= maxGaps; i++) {

    /* If the columns' number with this gaps' number is not equal to zero, we will count the columns. */
    if(numColumnsWithGaps[i] != 0) {

      /* Compute and prints the accumulative values for the gaps in the alignment. */
      acm += numColumnsWithGaps[i];
      cout << "  " << setiosflags(ios::left) << numColumnsWithGaps[i] << "\t\t" << setw(10) << (numColumnsWithGaps[i] * 100.0)/columns
           << "\t\t" << acm << "\t\t" << setw(10) << (acm * 100.0)/columns
           << "\t\t" << i << "\t\t" << setw(10) << (i * 1.0)/columnLength << "\t"<< setw(10) << 1 - ((i * 1.0)/columnLength) << endl;
    }
  }
}
