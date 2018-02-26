//
// Created by bioinfo on 19/01/18.
//

#ifndef MTRIMAL_TIMERFACTORY_H
#define MTRIMAL_TIMERFACTORY_H

//#define TimingReport true

#define BenchmarkTimes true
#define BenchmarkMemory false

#define TimingReport BenchmarkTimes || BenchmarkMemory

#if TimingReport

#define StartTiming(name) auto macroCustomTimer = timerFactory.getTimer(name)
#define SetTimingOfstream(_filename) timerFactory.SetOutput(_filename)
#include <chrono>
#include <string>
#include <stack>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <cstring>

class TimerFactory;

class Timer {
    friend class TimerFactory;

#if BenchmarkTimes
    /// Creation timestamp
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    /// Siblings duration, it serves to remove this time from the own calculations
    std::chrono::milliseconds siblings_duration;
#endif
    /// Name of the timer
    std::string name;
    /// How many childs does this timer have
    int childCount = 0;

#if BenchmarkMemory
    /// Method to chek current memory usage
    static int currentMemoryUsage();
#endif

public:
    static std::string separator, spacer;

    /// Constructor
    explicit Timer(std::string name, TimerFactory *timerFactory);
    /// TimeFactory pointer
    TimerFactory* timerFactory;
    /// Destructor
    ~Timer();

};



class TimerFactory {
    friend class Timer;
public:
    /// Class Destructor
    ~TimerFactory();
    /// Method to start timing an scope
    Timer getTimer(std::string name);

#if BenchmarkTimes
    /// Class constructor
    explicit TimerFactory(const std::string &outFilename)
            : lastTimer(""), accTimer(0),
              wholeStats(std::map<std::string, std::chrono::milliseconds>()),
              wholeUniqueStats(std::map<std::string, std::chrono::milliseconds>()),
              total(std::chrono::milliseconds(0)) {
        out = new std::ofstream(outFilename);
    };
#else
    TimerFactory(std::string outFilename)
            : lastTimer("") {
        out = new std::ofstream(outFilename);
    };
#endif

    void SetOutput(const std::string &_filename)
    {
        delete out;
        out = new std::ofstream(_filename);
    }
    /// Method to stop tracking and report
    void reportTotal();


    bool checkOutputParameter(int argc, char **argv);

private:
    /// Timer Level
    int timer = 1;
    /// Timer Pool
    std::vector<Timer *> _pool;
    /// Last Timer name. It serves to check if we are repeating the same step over and over.
    std::string lastTimer;
#if BenchmarkTimes
    /// Accumulative timer
    std::chrono::milliseconds accTimer;
    /// Map that contains the whole timing stats
    std::map<std::string, std::chrono::milliseconds> wholeStats;
    /// Map that contains the timing stats of timers, without including it's childs.
    std::map<std::string, std::chrono::milliseconds> wholeUniqueStats;
    /// Total time since the start of the first timer
    std::chrono::milliseconds total;
#endif
    /// How many times the same counter has been triggered
    int accTimerCount{0};

    /// Ostream where to output the report
    std::ostream * out;
    /// Cursor position of the output file.
    std::streampos cursorPos{0};
};

extern TimerFactory timerFactory;
#else

#define StartTiming(name) {}
#define SetTimingOfstream(_stream) {}

#endif

#endif //MTRIMAL_TIMERFACTORY_H
