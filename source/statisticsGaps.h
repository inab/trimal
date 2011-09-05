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
#ifndef STATISTICSGAPS_H
#define STATISTICSGAPS_H

#include <iostream>
#include <iomanip>

#include "utils.h"

#define DNAType 1
#define RNAType 2
#define AAType  3

using namespace std;

/* ***************************************************************************************************************** */
/*                                           Header Class File: StatisticsGaps.                                      */
/* ***************************************************************************************************************** */

class statisticsGaps {

  int columns;
  int columnLength;
  int maxGaps;
  int halfWindow;
  int dataType;

  int *gapsInColumn;
  int *numColumnsWithGaps;
  int *aminosXInColumn;
  int *gapsWindow;

 public:

  /* Class constructor without parameters. */
  statisticsGaps(void);

  /* Class destroyer. */
  ~statisticsGaps(void);

  /* Class constructor with parameters. */
  statisticsGaps(string *, int, int, int);

  /* Methods allows us compute the gapWindows' values. */
  bool applyWindow(int);

  /* This methods returns a gaps' vector reference. */
  int *getGapsWindow(void);

  /* Allows compute and select the cut point value. */
  double calcCutPoint(float, float);

  /* Automatic method to find a cut point value using the first and the second slopes. */
  int calcCutPointMixSlope(void);

  /* Automatic method to compute a cut point valur using the second slope approach. */
  int calcCutPoint2ndSlope(void);

  /* This methods print the gaps' percentage of each column in the alignment. */
  void printGapsColumns(void);

  /* This methods prints the statistics for the alignment relates to gaps. */
  void printGapsAcl(void);

};
#endif
