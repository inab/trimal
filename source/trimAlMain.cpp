/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated Alignment  trimming in large-scale
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

#include "InternalBenchmarker.h"
#include "trimalManager.h"
#include "reportsystem.h"

int main(int argc, char *argv[]) {

    // We have to check the output parameter of the timerFactory
    //  before it's called the first time
    #if TimingReport
        timerFactory.checkOutputParameter(argc, argv);
    #endif

    // Debug tag to prevent Debug Messages to be printed on user side
    debug.IsDebug = true;

    // To allow timing the whole program, we must encapsulate
    //  inside it's own scope. Thus, we store the returnValue outside
    //  the scope.
    int returnValue;
    {
        // Create a BenchmarkSnapshot that will report times and/or memory usage
        // upon its destruction
        StartTiming("int main(int argc, char *argv[]) ");

        // Create trimAl Manager
        trimAlManager trimAl;

        trimAl.parseArguments(argc, argv);

        // Process them: Incompatibilities and co-dependencies
        trimAl.processArguments(argv);

        // Perform the required actions
        returnValue = trimAl.perform();
    }

    // Return the exit code obtained in the perform method.
    return returnValue;
}
