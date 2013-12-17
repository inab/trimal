/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v1.4: a tool for automated alignment conversion among different
                 formats.

    statAl v1.4: a tool for getting descriptive alignment features/scores.

    2009-2013 Capella-Gutierrez S. and Gabaldon, T.
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

#define BUILD "2013-12-17"
#define VERSION 1.4
#define REVISION 15
#define AUTHORS "2009-2013. Salvador Capella-Gutierrez and Toni Gabaldón."

#define DNAType 1
#define RNAType 2
#define AAType  3
#define DNADeg  4
#define RNADeg  5

#define GAPPYOUT 1
#define STRICT   2

#define DELIMITERS     "   \t\n"
#define OTHDELIMITERS  "   \t\n,:"
#define OTH2DELIMITERS "   \n,:;"

#define HTMLBLOCKS 120
#define PHYLIPDISTANCE 10
