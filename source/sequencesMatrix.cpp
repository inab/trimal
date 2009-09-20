/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 

    trimAl v1.3: a tool for automated alignment trimming in large-scale 
                 phylogenetics analyses 

    Copyright (C) 2009 Capella-Gutierrez S. and Gabaldon, T.
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

 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#include "sequencesMatrix.h"

sequencesMatrix::sequencesMatrix(void) {

  columns = 0;
  columnLength = 0;
  matrix = NULL;

}

sequencesMatrix::sequencesMatrix(string *alignmentMatrix, int species, int aminos) {
  int i, j, k;

  columnLength = aminos;
  columns =      species;

  matrix = new int*[columns];
  for(i = 0; i < columns; i++) { 
    matrix[i] = new int[columnLength];
    utils::initlVect(matrix[i], columnLength, 0);
  }

  /* Determinate the sequence for each alignment specie */
  for(i = 0, k = 1; i < columns; i++, k = 1) {
    for(j = 0; j < columnLength; j++) {
      if(alignmentMatrix[i][j] != '-') {
        matrix[i][j] = k;
        k++;
      }
    }
  }
}

sequencesMatrix::~sequencesMatrix(void) {
  int i;

  if(matrix != NULL) {
    for(i = 0; i < columns; i++) 
      delete [] matrix[i];
    delete [] matrix;
  }

  matrix = NULL;
  columnLength = 0;
  columns = 0;
}

void sequencesMatrix::printMatrix(void) {
  int i, j, k;

  for(i = 0; i < columnLength; i += 20) {
    for(j = 0; j < columns; j++) {		
      for(k = i; k < (20 + i) && k < columnLength; k++) {
        cout << setw(4) << matrix[j][k] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
}

void sequencesMatrix::getColumn(int column, int *columnSeqMatrix) {

  int i;

  if(column < columnLength)
    for(i = 0; i < columns; i++)
      columnSeqMatrix[i] = matrix[i][column];

  else
    for(i = 0; i < columns; i++)
      columnSeqMatrix[i] = 0;

}

void sequencesMatrix::getColumn(int value, int row, int *columnSeqMatrix) {
  int i, j;

  for(i = 0; i < columnLength; i++)
    if(matrix[row][i] == value) break;

  if(i < columnLength)
    for(j = 0; j < columns; j++)
      columnSeqMatrix[j] = matrix[j][i];

  else
    for(j = 0; j < columns; j++)
      columnSeqMatrix[j] = -1;
}

void sequencesMatrix::setOrder(int *order) {
  int i, j, **resg;
    
  resg = new int*[columns];
  for(i = 0; i < columns; i++)
    resg[i] = new int[columnLength];
  
  for(i = 0; i < columns; i++)
    for(j = 0; j < columnLength; j++)
	  resg[i][j] = matrix[order[i]][j];

  for(i = 0; i < columns; i++) { 
    for(j = 0; j < columnLength; j++)
	  matrix[i][j] = resg[i][j];
	delete [] resg[i];
  }  
  delete [] resg;
}
