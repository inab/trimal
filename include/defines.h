#ifndef DEFINES
#define DEFINES

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v1.4: a tool for automated alignment conversion among different
                 formats.

    statAl v1.4: a tool for getting descriptive alignment features/scores.

    2009-2015 Capella-Gutierrez S. and Gabaldon, T.
              [scapella, tgabaldon]@crg.es

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

***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#define BUILD "2015-05-21"
#define VERSION 1.4
#define REVISION 22
#define AUTHORS "2009-2015. Salvador Capella-Gutierrez and Toni Gabaldón."

// #define DNAType 1
// #define RNAType 2
// #define AAType  3
// #define DNADeg  4
// #define RNADeg  5

#define GAPPYOUT 1
#define STRICT   2

#define DELIMITERS     "   \t\n"
#define OTHDELIMITERS  "   \t\n,:"
#define OTH2DELIMITERS "   \n,:;"

#define HTMLBLOCKS 120
#define PHYLIPDISTANCE 10

 /**
\file defines.h
\enum SequenceTypes
\details 
<b> Binary Tag Enum.</b> \n
<a href="http://www.alanzucconi.com/2015/07/26/enum-flags-and-bitwise-operators/" > Alan Zucconi </a> has a great explanation of these.\n\n
 
This enum contains the tags needed to obtain sequence types by composition of its main elements.\n
 
To obtain a composed alignment type you can bitwise OR the tags: ComposedTag = SequenceTypes::DNA | SequenceTypes::DEG\n
You can Fuzzy and Exact check both simple and composed tags\n\n
 
To \b FUZZY check an alignment type you can: \n
    _alignment -> getAlignmentType() & SequenceTypes::DNA\n
        <i> This will return true if the _alignment type is DNA, ignoring the rest of tags.</i>\n\n
    _alignment -> getAlignmentType() & (SequenceTypes::DNA | SequenceTypes::DEG)\n
        <i> This will return true if the _alignment type is DNA Deg, ignoring the rest of tags.</i>\n\n
    
To \b EXACT check an alignment type you can: \n
    _alignment -> getAlignmentType() == SequenceTypes::DNA \n
        <i> This will return true if the _alignment type is only DNA. DNA Deg would result in false</i>.\n\n
    _alignment -> getAlignmentType() == (SequenceTypes::DNA | SequenceTypes::DEG) \n
        <i> This will return true if the _alignment is type DNA Deg, additional tags (like SequenceTypes::DNA | SequenceTypes::RNA | SequenceTypes::DEG) would resut in false.</i>\n\n
*/
enum SequenceTypes
{
    NotDefined = 0,
    
    DNA = 1 << 1, ///< DNA Tag = 2
    RNA = 1 << 2, ///< RNA Tag = 4
    AA  = 1 << 3, ///< AA Tag = 8
    
    DEG = 1 << 4 ///< Degraded Tag = 16
};

#endif
