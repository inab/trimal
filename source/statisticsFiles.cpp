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

#include "statisticsFiles.h"

statisticsFiles::statisticsFiles() {
  columns = 0;
  columnLength = 0;
  sequencesMatrix = NULL;
}

statisticsFiles::statisticsFiles(char **alignmentMatrix, int species, int aminos) {
  int i;

  columnLength = species;
  columns =      aminos;

  sequencesMatrix = new int*[columnLength];
  for(i = 0; i < columnLength; i++)
    sequencesMatrix[i] = new int[columns];
}

statisticsFiles::~statisticsFiles() {
  int i;

  for(i = 0; i < columnLength; i++)
    delete[] sequencesMatrix[i];
  delete[] sequencesMatrix;

  sequencesMatrix = NULL;
  columnLength = 0;
  columns = 0;
}
