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

#ifndef SIMMatrix
#define SIMMatrix

#define NUMAMINOS 20
#define TAMABC 28
#define LINE_LENGTH 256
#define REFER 65

#include "values.h"
#endif

#include "similarityMatrix.h"
#include "utils.h"

#include <iostream>

#include <string.h>
#include <stdlib.h>

using namespace std;

extern char listSym[21];
extern float defaultMatrix[20][20];

/*+++++++++++++++++++++++++++++++++++++++++++++
| similarityMatrix::similarityMatrix()        |
|      Class constructor.                     |
+++++++++++++++++++++++++++++++++++++++++++++*/

similarityMatrix::similarityMatrix(){
  numPositions = 0;
  vhash        = NULL;
  simMat       = NULL;
  distMat      = NULL;
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| void similarityMatrix::memoryAllocation(int)                 |
|      This method allocates memory for some class attributes  |
|      with a number of positios given as the method parameter |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void similarityMatrix::memoryAllocation(int nPos){
  int i, j;

  /* Initializate square table dimension to store the distances */
  /* and to store the similarity matrix.                        */
  if(numPositions != 0) memoryDeletion();
  numPositions = nPos;

  /* Reserve memory for all structures */
  vhash = new int[TAMABC];

  simMat = new float *[nPos];
  distMat = new float *[nPos];

  for(i = 0; i < nPos; i++) {
    simMat[i]  = new float[nPos];
    distMat[i] = new float[nPos];

    for(j = 0; j < nPos; j++) {
      distMat[i][j] = 0.0;
      simMat[i][j] = 0.0;
    }
  }
}


/*++++++++++++++++++++++++++++++++++++++++++++++
| similarityMatrix::~similarityMatrix()        |
|      Class destructor  .                     |
++++++++++++++++++++++++++++++++++++++++++++++*/

similarityMatrix::~similarityMatrix(){

  if(numPositions != 0) memoryDeletion();

}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| void similarityMatrix::memoryDeletion()                 |
|      This method deletes all previously reserved memory |
|      for the object attributes                          |
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void similarityMatrix::memoryDeletion(){
  int i;

  for(i = 0; i < numPositions; i++){
    delete[] simMat[i]; delete[] distMat[i];
  }

  delete[] simMat;
  delete[] distMat;
  delete[] vhash;

  numPositions = 0;
  vhash        = NULL;
  simMat       = NULL;
  distMat      = NULL;
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| bool similarityMatrix::loadSimMatrix(char *)                       |
|      This method loads a similarity matrix from a file             |
|      and checks if the file format is correct. In that case        |
|      the method return true. If the format is incorrect the method |
|      returns false.                                                |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool similarityMatrix::loadSimMatrix(char *fn){
  char aux[LINE_LENGTH+1], first[LINE_LENGTH], listSym[LINE_LENGTH+1];
  int i,j, k; float sum; bool firstColumn = true; ifstream file;

  /* We try to open the file, if we can't open the file */
  /* we return false.                                   */
  file.open(fn);
  if(file.fail()) return false;

  /* Read the first line of the file and, depending on the   */
  /* line length (free of spaces an tabulators), we allocate */
  /* memory for the object structures                        */
  file.getline(aux, LINE_LENGTH);
  utils::removeSpaces(aux, listSym);
  memoryAllocation(strlen(listSym));

  for(i = 0; i < TAMABC; i++) vhash[i] = -1;

  /* We create the hashing vector */
  for(i = 0; i < numPositions; i++) {
    listSym[i] = (char)toupper((int)listSym[i]);

    if((listSym[i] >= 'A') && (listSym[i] <= 'Z')) {
      if((vhash[listSym[i] - 'A']) != -1) {
	memoryDeletion(); return false;
      }
      vhash[listSym[i] - 'A'] = i;

    } else {
      memoryDeletion(); return false;
    }
  }

  for(i = 0; i < numPositions; i++) {
    /* Read the first symbol of the line */
    j = 0; file >> first;

    /* If the format includes the first aminoacid in the line */
    if(firstColumn) {
      first[0] = (char)toupper((int)first[0]);

      /* Format checking. The first token must be a valid number */
      if(((first[0] >= '0' && first[0] <= '9') || (first[0] == '-' && (first[1] >= '0' && first[1] <= '9'))) && i > 0) {
        memoryDeletion(); return false;
      }

      /* If in the token is a character, there is "first column" */
      /* in the format of the alignment                          */
      if((first[0] >= 'A') && (first[0] <= 'Z')) {
	firstColumn = true;

	if((vhash[first[0] - 'A']) == -1){
	  memoryDeletion(); return false;
	}
      }

      /* If we have read a number there is no "first column" */
      /* in the format of the alignment                      */
      else if((first[0] >= '0' && first[0] <= '9') || (first[0] == '-' && (first[1] >= '0' && first[1] <= '9'))){
	firstColumn = false; j = 1;

	simMat[i][0] = atof(first);
	first[0] = listSym[i];
      }

    } else {
      j = 1;

      /* Do some checkings */
      if((first[0] >= 'A') && (first[0] <= 'Z') && (i > 0)) {
        memoryDeletion(); return false;
      }

      simMat[i][0] = atof(first);
      first[0] = listSym[i];
    }

    /* Read the corresponding number row */
    for(; j < numPositions; j++)
      file >> simMat[vhash[first[0] - 'A']][j];
  }

  /* Calculate the average between two simmetric positions            */
  /* respect to the diagonal of the matrix (between [i][j] and [j][i] */
  /* If the input is a non-symmetric matrix, the output will be a     */
  /* symmetric matrix                                                 */

  for(i = 0; i < numPositions; i++) {
    for(j = i+1; j < numPositions; j++) {
      if(simMat[i][j] != simMat[j][i]) {
        float value = (simMat[i][j] + simMat[j][i]) / 2.0;
        simMat[i][j] = value;
        simMat[j][i] = value;
      }
    }
  }

  /* Calculate the distances between aminoacids */
  /* based on Euclidean distance                */

  for(j = 0; j < numPositions; j++){
    for(i = 0; i < numPositions; i++) {
      if((i != j) && (distMat[i][j] == 0.0)){
        for(k = 0, sum = 0; k < numPositions; k++)
          sum += ((simMat[k][j] - simMat[k][i]) * (simMat[k][j] - simMat[k][i]));
        sum = (float) sqrt(sum);
        distMat[i][j] = sum;
        distMat[j][i] = sum;
      }
    }
  }

  file.close();
  return true;
}

void similarityMatrix::defaultAASimMatrix(void) {

  int i,j, k;
  float sum;

  memoryAllocation(20);
  for(i = 0; i < TAMABC; i++)
    vhash[i] = -1;

  /* We create the hashing vector */
  for(i = 0; i < numPositions; i++)
    vhash[listAASym[i] - 'A'] = i;

  for(i = 0; i < numPositions; i++)
    for(j = 0; j < numPositions; j++)
      simMat[i][j] = defaultAAMatrix[i][j];

  /* Calculate the distances between aminoacids */
  /* based on Euclidean distance                */
  for(j = 0; j < numPositions; j++) {
    for(i = 0; i < numPositions; i++) {
      if((i != j) && (distMat[i][j] == 0.0)) {
        for(k = 0, sum = 0; k < numPositions; k++)
          sum += ((simMat[k][j] - simMat[k][i]) * (simMat[k][j] - simMat[k][i]));
        sum = (float) sqrt(sum);
        distMat[i][j] = sum;
        distMat[j][i] = sum;
      }
    }
  }
}

void similarityMatrix::defaultNTSimMatrix(void) {
  int i,j, k;
  float sum;

  memoryAllocation(5);
  for(i = 0; i < TAMABC; i++)
    vhash[i] = -1;

  /* We create the hashing vector */
  for(i = 0; i < numPositions; i++)
    vhash[listNTSym[i] - 'A'] = i;

  for(i = 0; i < numPositions; i++)
    for(j = 0; j < numPositions; j++)
      simMat[i][j] = defaultNTMatrix[i][j];

  /* Calculate the distances between aminoacids */
  /* based on Euclidean distance                */
  for(j = 0; j < numPositions; j++) {
    for(i = 0; i < numPositions; i++) {
      if((i != j) && (distMat[i][j] == 0.0)) {
        for(k = 0, sum = 0; k < numPositions; k++)
          sum += ((simMat[k][j] - simMat[k][i]) * (simMat[k][j] - simMat[k][i]));
        sum = (float) sqrt(sum);
        distMat[i][j] = sum;
        distMat[j][i] = sum;
      }
    }
  }
}



/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| void similarityMatrix::printMatrix()                                 |
|      This method prints the similarity matrix to the standard output |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void similarityMatrix::printMatrix(){

  for(int i = 0; i < numPositions; i++){
    for(int j = 0; j < numPositions; j++)
      cerr << setw(5) << right << simMat[i][j];
    cerr << endl;
  }
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| float similarityMatrix::getDistance(char,char)                            |
|      This method returns the distance between the two characters given    |
|      The two character can be aminoacid characters, nucleotide characters |
|      or any kind of characters. This depends on the defined characters in |
|      the similarity matrix file.                                          |
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
float similarityMatrix::getDistance(char a, char b){
  int numa, numb; char chA, chB;

  chA = (char)toupper((int) a);
  chB = (char)toupper((int) b);

  /* Search the first character position */
  if((chA >= 'A') && (chA <= 'Z')) numa = vhash[chA - 'A'];
  else { cerr << "Error: the symbol '" << a << "' is incorrect" << endl; return -1; }

  /* Search the second character position */
  if((chB >= 'A') && (chB <= 'Z')) numb = vhash[chB - 'A'];
  else { cerr << "Error: the symbol '" << b << "' is incorrect" << endl; return -1; }

  /* We check if the two character postions are valid positions */
  if(numa == -1) {
    cerr << "Error: the symbol '" << a << "' accesing the matrix is not defined in this object" << endl;
    return -1;
  }

  if(numb == -1) {
    cerr << "Error: the symbol '" << b << "' accesing the matrix is not defined in this object" << endl;
    return -1;
  }

  /* Return the distance value between a and b */
  return distMat[numa][numb];
}
