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

#ifndef DEFINES
#define DEFINES

#define BUILD "2019-05-08"
#define VERSION "2"
#define REVISION "RC"
#define AUTHORS "2009-2019. Victor Fernández-Rodríguez, Toni Gabaldón, and Salvador Capella-Gutierrez"

#define GAPPYOUT 1
#define STRICT   2

#define DELIMITERS     "   \t\n"
#define OTHDELIMITERS  "   \t\n,:"
#define OTH2DELIMITERS "   \n,:;"

#define HTMLBLOCKS 120
#define PHYLIPDISTANCE 10

#define MINPARALLELSIZE 50
#define NUMTHREADS 8

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
    alig -> getAlignmentType() & SequenceTypes::DNA\n
        <i> This will return true if the alig type is DNA, ignoring the rest of tags.</i>\n\n
    alig -> getAlignmentType() & (SequenceTypes::DNA | SequenceTypes::DEG)\n
        <i> This will return true if the alig type is DNA Deg, ignoring the rest of tags.</i>\n\n
    
To \b EXACT check an alignment type you can: \n
    alig -> getAlignmentType() == SequenceTypes::DNA \n
        <i> This will return true if the alig type is only DNA. DNA Deg would result in false</i>.\n\n
    alig -> getAlignmentType() == (SequenceTypes::DNA | SequenceTypes::DEG) \n
        <i> This will return true if the alig is type DNA Deg, additional tags (like SequenceTypes::DNA | SequenceTypes::RNA | SequenceTypes::DEG) would resut in false.</i>\n\n
*/
enum SequenceTypes
{
    // Not Defined Tag = 0
    /// = 0
    NotDefined = 0u,
    
    // DNA Tag = 2
    /// 1 << 1 = 2
    DNA = 1u << 1u,

    // RNA Tag = 4
    /// 1 << 2 = 4
    RNA = 1u << 2u,
    
    // AA Tag = 8
    /// 1 << 3 = 8
    AA  = 1u << 3u,
    
    // Degraded Tag = 16
    /// 1 << 4 = 16
    DEG = 1u << 4u
};

#endif
