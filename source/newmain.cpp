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

int main(int argc, char *argv[])
{
//     cout << "Starting" << endl;
    trimAlManager Trimal = trimAlManager();
//     cout << "Parsing" << endl;
    Trimal.parseArguments(argc, argv);
//     cout << "Processing" << endl;
    Trimal.processArguments(argv);
//     cout << "Performing" << endl;
    
//     std::string title = "New graphic";
//     std::string filename = "filename.svg";
//     
//     utils::streamSVG(NULL, NULL, 0, false, & title, & filename,);
    
    return Trimal.perform();
}
