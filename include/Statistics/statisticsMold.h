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

#ifndef STATISTICS_MOLD_H
#define STATISTICS_MOLD_H

#include <sstream>
#include <iomanip>
#include <string>

// Forward declaration
class newAlignment;

/// \brief Class to handle Conservation / Similarity statistics
class statisticsMold {
public:

    newAlignment * _alignment;

    int
    residues                = -1,
    sequences               = -1,
    halfWindowApplied       = -1;

    /* Conservation vectors */
    float *values           = nullptr;
    float *valuesWindow     = nullptr;

    int * refCounter        = nullptr;

    const static std::string statName;


public:
    statisticsMold(newAlignment *parentAlignment, statisticsMold *mold);

    explicit statisticsMold(newAlignment * parentAlignment);

    ~statisticsMold();

    virtual bool calculate() { return false; };

    virtual bool applyWindow(int _halfWindow);
    
    virtual bool isDefinedWindow();

    virtual float * getVector();

    virtual void printByColumn(bool calculateRelative = false); // NOLINT ->
    // Default values for virtual methods is 'prohibited'

    virtual void printAccumulative(bool calculateRelative = false); // NOLINT ->
    // Default values for virtual methods is 'prohibited'

};
#endif
