/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.5.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    2009-2020
        Fernandez-Rodriguez V.  (victor.fernandez@bsc.es)
        Capella-Gutierrez S.    (salvador.capella@bsc.es)
        Gabaldon, T.            (tgabaldon@bsc.es)

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

#include "defines.h"

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

  void defaultNTDegeneratedSimMatrix();

  void alternativeSimilarityMatrices(int, int);

  float getDistance(char, char);

  void printMatrix();
};
#endif
