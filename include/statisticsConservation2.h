/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    2009-2015 Capella-Gutierrez S. and Gabaldon, T.

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

#ifndef STATISTICS_CONSERVATION2_H
#define STATISTICS_CONSERVATION2_H

#include <math.h>
#include <iostream>
#include <iomanip>

#include "similarityMatrix.h"
#include "statisticsGaps.h"
#include "../include/utils.h"

// #define DNAType 1
// #define RNAType 2
// #define AAType  3
class newAlignment;
using namespace std;

/* ***************************************************************************************************************** */
/*                                       Header Class File: StatisticsConservation.                                  */
/* ***************************************************************************************************************** */

class statisticsConservation2 {
public:

    newAlignment * _alignment;

    int
    residues,
    sequences,
    halfWindow = -1;

    /* Conservation vectors */
    float *Q;
    float *MDK;
    float *MDK_Window;

    /** \brief Identity weight matrix between alignment rows */
    float **matrixIdentity;

    /** \brief Similarity matrix used to conservation calculations */
    similarityMatrix *simMatrix;


public:
    /** \brief Computes the matrix identity between alignment's columns. */
    void calculateMatrixIdentity();

    /** \brief Constructor without any parameters */
    statisticsConservation2(newAlignment * parentAlignment);

    /** \brief Destructor */
    ~statisticsConservation2(void);

    /**
        \brief Method to calculate the conservation values of a alignment matrix.
        \param gaps Vector containing the gaps info.
    */
    bool calculateVectors(int *gaps);

    /**
        \brief Allows us compute the conservationWindow's values.
        \param _halfWindow Half window size to apply.
        \return \b False if there is a previously computed vector for this window size or half window size is greater than 1/4 of the alignment length.
    */
    bool applyWindow(int _halfWindow);
    
    /** 
        \brief Returns if a windows size value has been defined or not.
        \return \b True if a windows size has been defined.\b False otherwise.
    */
    bool isDefinedWindow(void);

    /**
        \brief This methods returns a pointer to conservationWindow's vector
        \return Conservation window vector.
    */
    float *getMdkwVector(void);

    /**
        \brief Stores a valid similarity matrix point to use.
        \param sm Similarity matrix pointer to associate.
        \return \b True if sm is valid, \b False if it's null
    */
    bool setSimilarityMatrix(similarityMatrix * sm);

    /**
        \brief Returns if a similarity matrix is being used or not.
        \return \b True if there is a similarity matrix set, \b False otherwise.
    */
    bool isSimMatrixDef(void);

    /**
        \brief Computes and selects the cut point values based on conservation's values.
        \todo This description seems a little vague.
        \param baseLine Percentage of columns desired.
        \param conservationPct Percentage of conservation desired.
    */
    double calcCutPoint(float baseLine, float conservationPct);

    /** \brief Prints the conservation's value for each alignment's column. */
    void printConservationColumns(void);

    /** \brief Computes and prints the accumulative statistics associated to the alignment. */
    void printConservationAcl(void);

};
#endif
