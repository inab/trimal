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

#ifndef MTRIMAL_TIMERFACTORY_H
#define MTRIMAL_TIMERFACTORY_H

//#define TimingReport true

#define BenchmarkTimes false
#define BenchmarkMemory false
/// Macro only valid if compiling with GCC
#define BenchmarkingExtendedNames true

#define TimingReport BenchmarkTimes || BenchmarkMemory

#if TimingReport

#include <algorithm>
#include <iostream>
#include <cstring>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <utility>
#include <chrono>
#include <vector>
#include <ctime>
#include <stack>
#include <map>

typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_point;

#define GCC_COMPILER (defined(__GNUC__) && !defined(__clang__))

// Detect if we are using GCC
#if GCC_COMPILER
// If using GCC, check if the developer wants extended names
#if BenchmarkingExtendedNames
// Helping macros to allow the use of location-dependent macros within macros
#define S1(x) #x
#define S2(x) S1(x)
// Macro that summarizes the extended names
#define LOCATION std::string(__PRETTY_FUNCTION__) + " :: " + std::string(__FILE__) + ":" + std::string(S2(__LINE__))

// In case we use GCC, we ignore the name provided.
#define StartTiming(name) auto macroCustomTimer = timerFactory.getTimer(LOCATION)
#else
// In case we don't want extended names, we use the pretty_function name
#define StartTiming(name) auto macroCustomTimer = timerFactory.getTimer(__PRETTY_FUNCTION__)
#endif

#else
#define StartTiming(name) auto macroCustomTimer = timerFactory.getTimer(name)
#endif

#define SetTimingOfstream(_filename) timerFactory.SetOutput(_filename)



// Forward Declaration
class InternalBenchmarker;

class BenchmarkSnapshot {
    friend class TimerFactory;

#if BenchmarkTimes
    /// Creation timestamp
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
#endif
    /// Name of the timer
    std::string name;

#if BenchmarkMemory
    /// Method to chek current memory usage
    static int currentMemoryUsage();
#endif

public:

    static std::string separator, spacer;
    /// Constructor
    explicit BenchmarkSnapshot(std::string name, InternalBenchmarker *timerFactory);
    /// TimeFactory pointer
    InternalBenchmarker* timerFactory;
    /// Destructor
    ~BenchmarkSnapshot();

};



class InternalBenchmarker {
    friend class BenchmarkSnapshot;
public:
    /// Class Destructor
    ~InternalBenchmarker();
    /// Method to start timing an scope
    BenchmarkSnapshot getTimer(std::string name);

    /// Class constructor
    explicit InternalBenchmarker(const std::string &outFilename) {
        out = new std::ofstream(outFilename);
        startPoint = std::chrono::high_resolution_clock::now();
    };

    void SetOutput(const std::string &_filename)
    {
        delete out;
        out = new std::ofstream(_filename);
    }
    /// Method to stop tracking and report
    void reportTotal();

    bool checkOutputParameter(int argc, char **argv);

private:
    int lastLevel = 0;
    std::vector<int> levels;
    /// Timer Pool
    std::vector<BenchmarkSnapshot *> _pool;
    /// Ostream where to output the report
    std::ostream * out;
    /// Cursor position of the output file.
    std::streampos cursorPos{0};

    std::map<std::string, long> totalDuration;
    std::map<std::string, long> relativeDuration;

    std::string lastName = "";
    int timesStepRepeated = 1;
    time_point firstStepStart{};

    time_point startPoint;
};

extern InternalBenchmarker timerFactory;
#else

#define StartTiming(name) {}
#define SetTimingOfstream(_stream) {}

#endif

#endif //MTRIMAL_TIMERFACTORY_H
