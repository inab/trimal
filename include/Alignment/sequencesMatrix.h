/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2009-2019
        Fernandez-Rodriguez V.  (victor.fernandez@bsc.es)
        Capella-Gutierrez S.    (salvador.capella@bsc.es)
        Gabaldon, T.            (tgabaldon@crg.es)

    This file is part of trimAl/readAl.

    trimAl/readAl are free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl/readAl are distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl/readAl. If not, see <http://www.gnu.org/licenses/>.

***************************************************************************** */

#ifndef STATISTICSFILES_H
#define STATISTICSFILES_H

#include <iostream>
#include <iomanip>

#include "Alignment/Alignment.h"

using namespace std;

// class Alignment;

/** \brief Internal Alignment Class that represents a sequences matrix
 *
 * A Sequence matrix is a 2D matrix that represents
 *  MSA without any gaps included.\n
 *
 * This class stores the alignment sequences matrix.\n It provides
 * methods to \b build the sequences matrix and print the matrix.\n
 * It also provides methods for look to a column in the matrix and
 * for look to value at the position (row, column) in the matrix.\n
 * See: \n
 * - Statistics::Consistency::compareAndChoose
 * - Statistics::Consistency::forceComparison
 */
class Alignment::sequencesMatrix {
    /**
     * \brief Number of residues per sequence.
     */
    int resNumber;
    /**
     * \brief Number of sequences in the matrix.
     */
    int seqsNumber;

    /**
     * \brief Matrix container. \n
     * It contains the sequences and residues of each sequence.
     */
    int **matrix;

    /**
     * \brief Sequences names container.
     */
    string *seqsName;

public:

    /** \brief Null constructor.
     *
     * This construction method initializes all attributes
     * of the new object with 0 or nullptr value.
     */
    sequencesMatrix();

    /** \brief Manual constructor.
     *
     * This construction method initializes all attributes
     * using the information passed as arguments.
     *
     * \param alignmentMatrix
     *  Sequences of the alignment to convert to sequenceMatrix
     * \param alignmentSeqsName
     *  Sequences names
     * \param sequences
     *  Number of sequences
     * \param residues
     *  Number of residues
     */
    sequencesMatrix(string *alignmentMatrix, string *alignmentSeqsName,
                    int sequences, int residues);

    /** \brief Automatic constructor.
     *
     * This construction method initializes all attributes
     * using the information present in the alignment pointer passed.
     *
     * \param parent Alignment to associate the sequences matrix to.
     */
    explicit sequencesMatrix(Alignment *parent);

    sequencesMatrix &operator=(const sequencesMatrix &);

    /**
     * \brief Destructor.
     *
     * Destruction method that frees, if exists, previously allocated memory.
     */
    ~sequencesMatrix();

    /** \brief Sequences Matrix printing method.
     *
     * Method that prints the alignment sequences matrix.
     *
     * \warning Not In Use
     */
    void printMatrix();

    /** \brief Method to get a column out of the matrix
     *  \param column
     *   Column number at sequences matrix.
     *  \param [out] numResidueseqMatrix
     *   Vector where storage a column's sequences matrix.
     *
     *  Method that storages a column's sequences matrix in a vector.
     */
    void getColumn(int column, int *numResidueseqMatrix);

    /** \brief
     * Method that looks to value in a row and stores a column's,
     * corresponding to row, sequences matrix in a vector.
     * \param value to look in a row's sequences matrix.
     * \param row where to look for a value.
     * \param columnSeqMatrix Vector where storage a column's sequences matrix.
     *
     */
    void getColumn(int value, int row, int *columnSeqMatrix);

    /**
     * \brief Method that reorders the stored sequences with a given order list.
     * \param order Pointer that contains the new order we want to give the matrix.
     */
    void setOrder(int *order);

    /**
     * \brief Method to obtain a sequence based on its name.
     * \param seqName Name of the sequence to find.
     * \param[out] sequence Where to store the sequence.
     * \return \b True if found. \b False if not.
     */
    bool getSequence(string seqName, int *sequence);

    /**
     * \brief Number of sequences getter.
     * \return Number of sequences.
     */
    int getSeqNumber();

    /**
     * \brief Number of residues getter.
     * \return Number of residues.
     */
    int getResidNumber();

private:
    /**
     * \brief Pointer to the alignment this object is related to.
     */
    Alignment *alig;
};

#endif
