//
// Created by bioinfo on 19/01/18.
//

#ifndef MTRIMAL_TIMERFACTORY_H
#define MTRIMAL_TIMERFACTORY_H

//#define TimingReport true

#define BenchmarkTimes false
#define BenchmarkMemory false

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

#define StartTiming(name) auto macroCustomTimer = timerFactory.getTimer(name)
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
