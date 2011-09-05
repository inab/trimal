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

#ifndef STATISTICSFILES_H
#define STATISTICSFILES_H

#include <iostream>
#include <iomanip>

#include "utils.h"

using namespace std;

/** \brief Class containing a sequences matrix
 *
 * This class stores the alignment sequences matrix. It provides
 * methods to \b build the sequences matrix and print the matrix.
 * It also provides methods for look to a column in the matrix and
 * for look to value at the position (row, column) in the matrix.
 */

class sequencesMatrix {
  int resNumber;
  int seqsNumber;

  /* Sequences Matrix */
  int **matrix;

  /* Sequences Name */
  string *seqsName;

  public:

  /* Constructors */

  /** \brief Null constructor.
   *
   * This construction method initializates all attributes
   * of the new object with 0 or NULL value.
   */
  sequencesMatrix(void);

  /* Copy constructor */
  sequencesMatrix(string *, string *, int, int);

  sequencesMatrix &operator=(const sequencesMatrix &);

  /* Destructor */

  /** \brief Destructor.
   *
   * Destruction method that frees, if exists, previously allocated memory.
   */
  ~sequencesMatrix();

  /* Basics Operations. */

  /** \brief Sequences Matrix printing method.
   *
   * Method that prints the alignment sequences matrix.
   */
  void printMatrix();

  /** \brief Column for looking to method.
   * \param column Column number at sequences matrix.
   * \param numResidueseqMatrix Vector where storage a column's sequences matrix.
   *
   * Method that storages a column's sequences matrix in a vector.
   */
  void getColumn(int, int *);

  /** \brief Column for looking to method.
   * \param value to look in a row's sequences matrix.
   * \param row where to look for a value.
   * \param numResidueseqMatrix Vector where storage a column's sequences matrix.
   *
   * Method that looks to value in a row and storages a column's, corresponding to row,
   * sequences matrix in a vector.
   */
  void getColumn(int, int, int *);

  void setOrder(int *);

  void removeColumns(int, int, int *, int *);

  bool getSequence(string, int *);

  int getSeqNumber(void);

  int getResidNumber(void);

};

#endif
