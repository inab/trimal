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

#ifndef STATISTICS_CONSERVATION_H
#define STATISTICS_CONSERVATION_H

#include <math.h>
#include <iostream>
#include <iomanip>

#include "similarityMatrix.h"
#include "statisticsGaps.h"
#include "utils.h"

#define DNAType 1
#define RNAType 2
#define AAType  3

using namespace std;

/* ***************************************************************************************************************** */
/*                                       Header Class File: StatisticsConservation.                                  */
/* ***************************************************************************************************************** */

class statisticsConservation{
 private:

  /* Number of columns and sequences of the alignment */
  int columns;
  int sequences;

  /* Sequence's Datatype: DNA, RNA or Amino Acids. */
  int dataType;

  /* Half window size */
  int halfWindow;

  /* Conservation vectors */
  float *Q;
  float *MDK;
  float *MDK_Window;

  /* Identity weight matrix between alignment rows */
  float **matrixIdentity;

  /* Similarity matrix used to conservation calculations */
  similarityMatrix *simMatrix;

  /* Private methods */
  /* Computes the matrix identity between alignment's columns. */
  void calculateMatrixIdentity(string *alignmentMatrix);

 public:

  /* Constructors without any parameters */
  statisticsConservation(void);

  /* Constructors using parameters */
  statisticsConservation(string *, int, int, int);

  /* Destroyer */
  ~statisticsConservation(void);

  /* This methods allows us compute the alignment's conservation's values. */
  bool calculateVectors(string *, int *);

  /* Allows us compute the conservationWindow's values. */
  bool applyWindow(int);

  /* Returns if a windows size value has been defined or not. */
  bool isDefinedWindow(void);

  /* This methods returns a pointer to conservationWindow's vector */
  float *getMdkwVector(void);

  /* Associates a pointer to similarity matrix. This matrix is needed to compute the conservation's values. */
  bool setSimilarityMatrix(similarityMatrix *);

  /* Returns if a similarity matrix is being used or not. */
  bool isSimMatrixDef(void);

  /* Computes and selects the cut point values based on conservation's values. */
  float calcCutPoint(float, float);

  /* Prints the conservation's value for each alignment's column. */
  void printConservationColumns(void);

  /* Computes and prints the accumulative statistics associated to the alignment. */
  void printConservationAcl(void);

};
#endif
