/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated newAlignment  trimming in large-scale
                 phylogenetics analyses.

    2009-2013 Capella-Gutierrez S. and Gabaldon, T.
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

#include "../include/trimalArgumentParser.h"
#include "../include/reportsystem.h"

int main(int argc, char *argv[])
{
    debug.IsDebug = true;
    
//     debug.Level = VerboseLevel::ERROR;
    
    debug.log(VerboseLevel::ERROR) << "This is the origin branch" << endl;

    trimAlManager Trimal = trimAlManager();

    Trimal.parseArguments(argc, argv);
    
    Trimal.processArguments(argv);

    return Trimal.perform();
    
}
