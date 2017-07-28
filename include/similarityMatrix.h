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

#ifndef SIMILARITYMATRIX_H
#define SIMILARITYMATRIX_H

#include <math.h>
#include <ctype.h>
#include <stdlib.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

class similarityMatrix{
  int *vhash;
  float **simMat;
  float **distMat;
  int numPositions;

 private:
  void memoryAllocation(int);
  void memoryDeletion();

 public:
  similarityMatrix();

  ~similarityMatrix();

  bool loadSimMatrix(char *);

  void defaultAASimMatrix();

  void defaultNTSimMatrix();

  float getDistance(char, char);

  void printMatrix();
};
#endif
