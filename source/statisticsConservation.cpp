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

#include "statisticsConservation.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  statisticsConservation::statisticsConservation(char **, int, int)                                                   |
|                                                                                                                      |
|       Class constructor. This method uses the inputs parameters to put the information in the new object that        |
|       has been created.                                                                                              |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

statisticsConservation::statisticsConservation(string *alignmentMatrix, int species, int aminos, int dataType_) {

  /* Initializate values to its corresponds values */
  columns = aminos;
  sequences = species;
  dataType = dataType_;
  halfWindow = -1;

  /* Allocate memory to the structures and initializates it */
  Q = new float[columns];
  utils::initlVect(Q, columns, 0);

  MDK = new float[columns];
  utils::initlVect(MDK, columns, 0);

  MDK_Window = new float[columns];
  utils::initlVect(MDK_Window, columns, 0);

  matrixIdentity = new float*[sequences];
  for(int i = 0; i < sequences; i++){
    matrixIdentity[i] = new float[sequences];
    utils::initlVect(matrixIdentity[i], sequences, 0);
  }

  /* Initializate the similarity matrix to NULL. */
  simMatrix = NULL;

  /* Calculation methods call */
  calculateMatrixIdentity(alignmentMatrix);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  statisticsConservation::statisticsConservation(void)                                                                |
|                                                                                                                      |
|       Class constructor.                                                                                             |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

statisticsConservation::statisticsConservation(void) {

  /* Initializate all values to 0 */
  columns = 0;
  sequences = 0;
  halfWindow = 0;

  /* and the pointers to NULL */
  Q = NULL;
  MDK = NULL;
  MDK_Window = NULL;

  matrixIdentity = NULL;
  simMatrix = NULL;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  statisticsConservation::~statisticsConservation(void)                                                               |
|                                                                                                                      |
|       Class destroyer.                                                                                               |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

statisticsConservation::~statisticsConservation(void) {

  int i;

  /* Deallocate memory, if it have been allocated previously. */
  if(Q != NULL) {
    delete[] Q;
    delete[] MDK;
    delete[] MDK_Window;

    for(i = 0; i < sequences; i++)
      delete[] matrixIdentity[i];
    delete[] matrixIdentity;
  }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void statisticsConservation::calculateMatrixIdentity(char **, int, int)                                             |
|                                                                                                                      |
|       This method computes the matrix identity between all the sequences in the alignment.                           |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void statisticsConservation::calculateMatrixIdentity(string *alignmentMatrix) {

  char indet;
  int i, j, k, sum, length;

  if(dataType == AAType)
    indet = 'X';
  else
    indet = 'N';

  /* For each sequences' pair */
  for(i = 0; i < sequences; i++) {
    for(j = i + 1; j < sequences; j++) {

      /* For each position in the alignment of that pair than we are processing */
      for(k = 0, sum = 0, length = 0; k < columns; k++) {

        /* If we find a element that is not a gap or an X aminoacid in the first sequence of the pair */
        if((alignmentMatrix[i][k] != '-') && (alignmentMatrix[i][k] != indet)) {

          /* If we also find a valid element in the second sequence  */
          if((alignmentMatrix[j][k] != '-') && (alignmentMatrix[j][k] != indet))

            /* If the two valid elements are the same increase the sum */
            if(alignmentMatrix[j][k] ==  alignmentMatrix[i][k])
              sum++;

          /* Increase the length of the sequence free of gaps and X elements */
          length++;
        }

        /* If the first processed element is invalid and in the second we find a valid element increase the length of
           the sequence free of gaps and X elements */
        else if((alignmentMatrix[j][k] != '-') && (alignmentMatrix[j][k] != indet))
          length++;
      }

      /* Calculate the value of matrixidn for columns j and i */
      matrixIdentity[j][i] = (100.0 - ((float) sum/ length) * 100.0);
      matrixIdentity[i][j] = matrixIdentity[j][i];
    }
  }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool statisticsConservation::calculateVectors(char **, int *)                                                       |
|                                                                                                                      |
|       This method computes the distance between pairs for each column in the alignment.                              |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool statisticsConservation::calculateVectors(string *alignmentMatrix, int *gaps) {

  char indet;
  int i, j, k;
  float num, den;

  if(dataType == AAType)
    indet = 'X';
  else
    indet = 'N';

  /* A conservation matrix must be defined. If not, return false */
  if(simMatrix == NULL)
    return false;

  /* For each column calculate the Q value and the MD value using an equation */
  for(i = 0; i < columns; i++) {
    /* For each AAs/Nucleotides' pair in the column we compute its distance */
    for(j = 0, num = 0, den = 0; j < sequences; j++) {
      /* We don't compute the distant if the first element is a indeterminate (X) or a gap (-) element. */
      if((alignmentMatrix[j][i] != '-') && (alignmentMatrix[j][i] != indet))
        for(k = j + 1; k < sequences; k++)
          /* We don't compute the distant between the pair if the second element is a indeterminate or a gap element */
          if((alignmentMatrix[k][i] != '-') && (alignmentMatrix[k][i] != indet)) {
            /* We use the identity value for the two pairs and its distance based on similarity matrix's value. */
            num += matrixIdentity[j][k] * simMatrix -> getDistance(alignmentMatrix[j][i], alignmentMatrix[k][i]);;
            den += matrixIdentity[j][k];
          }
    }
    /* If we are procesing a column with only one AA/nucleotide, the denominator is 0 and we don't execute the division
       and we set the Q[i] value to 0. */
    Q[i] = (den == 0) ? 0 : num / den;
    MDK[i] = (float) exp(-Q[i]);

    /* If the column has 80% or more gaps then we set its conservation value to 0 */
    if(gaps != NULL)
      if(((float) gaps[i] / sequences) >= 0.8) MDK[i] = 0;

    /* If the MDK value is more than 1, we normalized this value to 1. */
    if(MDK[i] > 1) MDK[i] = 1;
  }

  return true;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool statisticsConservation::applyWindow(int)                                                                       |
|                                                                                                                      |
|       This method computes for each column's alignment its conservationwindows' value. For this purpose, the method  |
|       uses the values that previously has been calculated and the window's size value.                               |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool statisticsConservation::applyWindow(int _halfWindow) {

  int i, j, window;

  /* If one of this conditions is true, we return FALSE:                         */
  /*    .- If already exists a previously calculated vector for this window size */
  /*    .- If mediumWinSize value is greater than 1/4 of alignment length        */
  if((halfWindow == _halfWindow) || (_halfWindow > columns/4))
     return false;

  halfWindow = _halfWindow;
  window = 2 * halfWindow + 1;

  /* Do the average window calculations */
  for(i = 0; i < columns; i++) {
    for(j = i - halfWindow; j <= i + halfWindow; j++) {
      if(j < 0) MDK_Window[i] += MDK[-j];
      else if(j >= columns) MDK_Window[i] += MDK[((2 * columns - j) - 2)];
      else MDK_Window[i] += MDK[j];
    }

    /* Calculate the similiraty value for the i column */
    MDK_Window[i] = MDK_Window[i] / (float) window;
  }
  return true;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool statisticsConservation::isDefinedWindow(void)                                                                  |
|                                                                                                                      |
|       This method returns true if a similarity matrix has been defined and false in others cases.                    |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool statisticsConservation::isDefinedWindow(void) {

  return (halfWindow != -1);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool statisticsConservation::getMdkwVector(void)                                                                    |
|                                                                                                                      |
|       This method returns a pointer to conservation values' vector.                                                  |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

float *statisticsConservation::getMdkwVector(void) {

  return MDK_Window;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool statisticsConservation::setSimilarityMatrix(similarityMatrix *)                                                |
|                                                                                                                      |
|       This method associated a pointer to similarity matrix gives as input parameter. If a conservation matrix is    |
|       being used the methods return false and doesn't do anything.                                                   |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool statisticsConservation::setSimilarityMatrix(similarityMatrix *sm) {

  /* Checks if a similarity matrix is being used. */
  if(sm == NULL)
    return false;

  /* if a similarity matrix isn't being used, we associate a pointer gives as input parameter to object simMatrix's
     pointer and return true. */
  simMatrix = sm;
  return true;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool statisticsConservation::isSimMatrixDef(void)                                                                   |
|                                                                                                                      |
|       This method returns true if a similarity matrix is being used and false in others cases.                       |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool statisticsConservation::isSimMatrixDef(void) {

  return (simMatrix != NULL);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  float statisticsConservation::calcCutPoint(float baseLine, float conservationPct)                                   |
|                                                                                                                      |
|       This method computes and selects the cut point value based on alignment's conservation. For this purpose, this |
|       select the minimum conservation value bewteen the conservationPct -a conservation value- and the conservation  |
|       value allows us conserve the columns' number associated to baseline.                                           |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

float statisticsConservation::calcCutPoint(float baseLine, float conservationPct) {

  float cutBL, cutCons, *vectAux;
  int i;

  /* Allocate memory */
  vectAux = new float[columns];

  /* Sort a copy of the MDK_Window vector, and take the value of the column that marks the % baseline */
  utils::copyVect(MDK_Window, vectAux, columns);
  utils::quicksort(vectAux, 0, columns-1);

  for(i = columns - 1; i >= 0; i--)
    if(vectAux[i] < conservationPct) break;
  cutCons = vectAux[i];

  cutBL = vectAux[(int) ((float)(columns - 1) * (100.0 - baseLine)/100.0)];

  /* Deallocate memory */
  delete[] vectAux;

  /* Return the minimum of the baseLine cut and conservationPct value */
  return (cutBL < cutCons ? cutBL : cutCons);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void statisticsConservation::printConservationColumns(void)                                                         |
|                                                                                                                      |
|       This method prints the conservation's value for each column in the alignment.                                  |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void statisticsConservation::printConservationColumns(void) {

  int i;

  /* We set the output precision and print the header. */
  cout << "| Residue\t Similarity  |" << endl;
  cout << "| Number \t    Value    |" << endl;
  cout << "+----------------------------+" << endl;
  cout.precision(10);

  /* If MDK_Window vector is defined, we use it to print the conservation's values. */
  if(MDK_Window != NULL)
    for(i = 0; i < columns; i++)
      cout << "  " << setw(5) << i << "\t\t" << setw(7) << MDK_Window[i] << endl;

  /* In others cases, we uses the MDK vector to print the conservation's vlaues. */
  else
    for(i = 0; i < columns; i++)
      cout << "  " << setw(5) << i << "\t\t" << setw(7) << MDK[i] << endl;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void statisticsConservation::printConservationAcl(void)                                                             |
|                                                                                                                      |
|       This method prints the accumulative statistics related to conservation in the alignment.                       |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void statisticsConservation::printConservationAcl(void) {

  float refer, *vectAux;
  int i, num, acm;

  /* Allocate memory */
  vectAux = new float[columns];

  /* Select the conservation's value source and copy that vector in a auxiliar vector */
  if(MDK_Window != NULL) utils::copyVect(MDK_Window, vectAux, columns);
  else utils::copyVect(MDK, vectAux, columns);

  /* Sort the auxiliar vector. */
  utils::quicksort(vectAux, 0, columns-1);

  /* We set the output precision and print the header. */
  cout << "| Number of\t        \t|\t Cumulative \t% Cumulative\t|   Similarity   |" << endl;
  cout << "| Residues \t% Length\t|\tNumberResid.\t   Length   \t|     Value      |" << endl;
  cout << "+-------------------------------+---------------------------------------+----------------+" << endl;
  cout.precision(10);


  /* Initializate some values */
  refer = vectAux[columns-1];
  acm = 0; num = 1;

  /* Count the columns with the same conservation's value and compute this information to shows the accunulative
     statistics in the alignment. */
  for(i = columns-2; i >= 0; i--) {
    acm++;

    if(refer != vectAux[i]) {
      cout << "  " << num << "\t\t" << setw(10) << ((float) num/columns * 100.0)
           << "\t\t" << acm << "\t\t" << setw(10) << ((float) acm/columns * 100.0) << "\t"
           << setw(15) << refer << endl;
      refer = vectAux[i];
      num = 1;
    }
    else num++;
  }
  acm++;
  cout << "  " << num << "\t\t" << setw(10) << ((float) num/columns * 100.0)
       << "\t\t" << acm << "\t\t" << setw(10) << ((float) acm/columns * 100.0) << "\t"
       << setw(15) << refer << endl;

  /* Deallocate the reserved memory. */
  delete [] vectAux;
}
