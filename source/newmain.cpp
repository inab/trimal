#include <TimerFactory.h>
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

#include <../include/ReadWriteMS/vcf_statish.h>


int main(int argc, char *argv[]) {

#if TimingReport
    // We have to check the output parameter of the timerFactory
    //  before it's called the first time, AKA main function timer
    timerFactory.checkOutputParameter(argc, argv);
#endif

    // This is to prevent debug calls.
    // They will be silenced to the final user
    //  although they appear in code.
    debug.IsDebug = true;

    // We store the return value before returning if as program exit code
    //  to allow the main function timer to be deleted.
    int returnValue;
    {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("int main(int argc, char *argv[]) ");

        // We create the trimAl Manager
        trimAlManager trimAl;

        // Feed the arguments
        trimAl.parseArguments(argc, argv);

        // Process them: Incompatibilities and co-dependencies
        trimAl.processArguments(argv);

        // Perform the required actions
        returnValue = trimAl.perform();
    }

    // Return the exit code
    return returnValue;
}
