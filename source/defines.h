/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.5.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v1.5.0: a tool for automated alignment conversion among different
                 formats.

    statAl v1.5.0: a tool for getting descriptive alignment features/scores.

    2009-2020
        Fernandez-Rodriguez V.  (victor.fernandez@bsc.es)
        Capella-Gutierrez S.    (salvador.capella@bsc.es)
        Gabaldon, T.            (tgabaldon@bsc.es)

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

#define BUILD "2024-05-27"
#define VERSION 1.5
#define REVISION 0
#define AUTHORS "2009-2020. Victor Fernández-Rodríguez, Toni Gabaldón, and Salvador Capella-Gutiérrez"

#define DNAType 1
#define RNAType 2
#define AAType  3
#define DNADeg  4
#define RNADeg  5

#define SINGLE  1
#define MULTI   2

#define GAPPYOUT 1
#define STRICT   2

#define DELIMITERS     "   \t\n"
#define OTHDELIMITERS  "   \t\n,:"
#define OTH2DELIMITERS "   \n,:;"

#define HTMLBLOCKS 120
#define PHYLIPDISTANCE 10
