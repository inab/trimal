#include "TimerFactory.h"
#include <utility>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <fstream>
#include <reportsystem.h>

#if TimingReport
TimerFactory timerFactory = TimerFactory(""); // NOLINT

//region TimerFactory
Timer TimerFactory::getTimer(std::string name) {
    // Create a new timer
    Timer timer = Timer(std::move(name), this);
    // Store its pointer in the current queue
    _pool.push_back(&timer);
    // Return the timer
    return timer; // NOLINT
}

void TimerFactory::reportTotal() {

    *out << "\n"
         << Timer::spacer
         << "'~ "
         #if BenchmarkTimes
         << " Duration: " << total.count() << " ms."
         #endif
         #if BenchmarkMemory
         << " VmRSS: " << Timer::currentMemoryUsage() << " kb"
#endif
            ;


#if BenchmarkTimes
    // Write a separator line to split the tree and the chart.
    *out << "\n";
    // Calculate the max size of the timer names
    int maxSize = -1;
    for (auto &&value : wholeStats) {
        maxSize = std::max(maxSize, static_cast<const int &>(value.first.substr(0, value.first.find("(")).length()));
        maxSize = std::max(maxSize, static_cast<const int &>(value.first.substr(value.first.find("(")).length()));
    }

    maxSize = std::max(/*len of "Avg Total Duration   "*/ 21, maxSize);
    maxSize += 5;

    // Titles
    *out
            << "\n" << "Final Report: " << lastTimer
            << "\n" << "Total Time: " << total.count() << "ms."

            // Titles separators
            << "\n"
            << std::setfill('-')
            << " " << std::setw(maxSize - 1) << " "
            << " " << std::setw(25) << " "
            << " " << std::setw(25) << " "
            << " " << std::setw(25) << " "
            << " " << std::setw(25) << " "
            << std::setfill(' ')

            // Titles
            << "\n"
            << " │" << std::setw(maxSize - 5) << "Key " << " │ "
            << " │" << std::setw(21) << "Duration " << " │ "
            << " │" << std::setw(21) << "Rel Duration " << " │ "
            << " │" << std::setw(21) << "Total Duration " << " │ "
            << " │" << std::setw(21) << "Rel Total Duration " << " │ "
            << "\n"

            // Titles separators
            << std::setfill('-')
            << " " << std::setw(maxSize - 1) << " "
            << " " << std::setw(25) << " "
            << " " << std::setw(25) << " "
            << " " << std::setw(25) << " "
            << " " << std::setw(25) << " "
            << std::setfill(' ')
            << "\n";

    // Foreach Timer triggered
    for (auto &&value : wholeStats) {
        // output the name
        *out << std::setw(maxSize - 2) << value.first.substr(0, value.first.find("(")) << "  "
             // Output the unique timing
             << std::setw(21) << wholeUniqueStats[value.first].count() << " ms "
             << std::fixed
             // Output the relative unique timing
             << std::setw(21) << wholeUniqueStats[value.first].count() * 100.F / total.count() << " %  "
             << std::defaultfloat
             // Output the total timing
             << std::setw(21) << value.second.count() << " ms "
             << std::fixed
             // Output the relative total timing
             << std::setw(21) << value.second.count() * 100.F / total.count() << " %  "
             << std::defaultfloat
             // Output the parameters of the method
             << "\n"
             << std::setw(maxSize - 2) << value.first.substr(value.first.find("(")) << "\n";
    }
#endif
}

TimerFactory::~TimerFactory() {

    // Output the chart
    reportTotal();
    // Delete the out pointer in case it's not a global one (cout, cerr, clog)
    if (out->rdbuf() != std::cout.rdbuf() &&
        out->rdbuf() != std::cerr.rdbuf() &&
        out->rdbuf() != std::clog.rdbuf())
        delete out;
}

bool TimerFactory::checkOutputParameter(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        if ((!strcmp(argv[i], "-timetrackerout")) && ((i) + 1 != argc)) {
            size_t argumentLength = strlen(argv[++i]);
            char * out = new char[argumentLength + 1];
            strcpy(out, argv[i]);
            this->SetOutput(out);
            std::cout << "Time Tracker Output changed to " << out << std::endl;
            delete [] out;
            return true;
        }
    }
    this->SetOutput("/dev/null");
    return false;
}
//endregion

//region Timer

std::string Timer::separator = std::string("│"); // NOLINT
std::string Timer::spacer = std::string("  "); // NOLINT


Timer::Timer(std::string name, TimerFactory *timerFactory) :
        name(std::move(name)),
        timerFactory(timerFactory)
#if BenchmarkTimes
        , siblings_duration(0)
#endif
{

#if BenchmarkTimes
    // Save the current timestamp
    start = std::chrono::high_resolution_clock::now();
#endif
    // If the Timer Manager has a current timer
    if (!timerFactory->_pool.empty()) {
        // We add a child to the last timer
        // This means that hierarchically we are a child of the previous timer
        timerFactory->_pool.back()->childCount += 1;
    }

    // We check if the last finished timer is not the same as the current one
    if (timerFactory->lastTimer != this->name) {
        // We add a new line
        *timerFactory->out << "\n";
        // And output the tree form
        int i = 0;
        for (i = 0; i < timerFactory->timer - 1; i++) {
            *timerFactory->out << spacer << separator;
        }

        if (timerFactory->timer == 1)
            *timerFactory->out << spacer << "┌──┬─> ";
        else
            *timerFactory->out << spacer << "├──┬─> ";
        // We output the numbered ID
        for (Timer *timer : timerFactory->_pool) {
            *timerFactory->out << timer->childCount << ".";
        }
        // And output the current memory located.
        *timerFactory->out << "  " << this->name;
    }

    // Increase the number of timers we have.
    timerFactory->timer++;

}

Timer::~Timer() {
    // Reduce the number of timers out Manager has
    timerFactory->timer--;
#if BenchmarkTimes
    // Calculate the diff of timestamps since creation
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);
#endif
    // Check if the last timer finished is the same as current.
    if (timerFactory->lastTimer == name && childCount == 0) {
        // Increase the counter of times we've called the same timer in different calls
        timerFactory->accTimerCount += 1;
#if BenchmarkTimes
        // Increase the counter of time elapsed within the same timer in different calls
        timerFactory->accTimer += duration;
#endif
        // Move the position to the cursorPos, to remove the last call.
        // The cursorPos will be the same on each repetition
        timerFactory->out->seekp(timerFactory->cursorPos);
        // Add the tree separators
        for (int i = 0; i < timerFactory->timer; i++) *timerFactory->out << spacer << separator;
        // Output the timing info: Name, Repetitions, Average and Accumulative.
        *timerFactory->out << spacer
                           << "└─>" //<< name
                           << " Repetitions: " << timerFactory->accTimerCount
                           #if BenchmarkTimes
                           << " Average: " << (timerFactory->accTimer.count() / timerFactory->accTimerCount) << "ms."
                           << " Acc " << timerFactory->accTimer.count() << "ms."
                           #endif
                           #if BenchmarkMemory
                           << " VmRSS: " << currentMemoryUsage() << " kb";
#else
        ;
#endif

    } else {
        // Set the counter of times we've called the same timer in different calls as 1
        timerFactory->accTimerCount = 1;
#if BenchmarkTimes
        // Set the counter of time elapsed within the same timer in different calls to the current timer duration
        timerFactory->accTimer = duration;
#endif
        // Set the last finished timer name as the current timer name
        timerFactory->lastTimer = name;

        // Output will be on the nex line.
        *timerFactory->out << "\n";
        // We store the current cursor position
        timerFactory->cursorPos = timerFactory->out->tellp();
        // Output the separator
        for (int i = 0; i < timerFactory->timer; i++) *timerFactory->out << spacer << separator;
        // Output the timing info: Name, Duration, Accumulative
        *timerFactory->out << spacer
                << "└─>" //<< name
               #if BenchmarkTimes
               << " Duration: " << (duration - siblings_duration).count() << " ms."
               << " Acc: " << duration.count() << " ms."
               #endif
               #if BenchmarkMemory
               << " VmRSS: " << currentMemoryUsage() << " kb"
#endif
                                                              ;
        // Output a new line
        *timerFactory->out << "\n";
        // Output the tree separators
        for (int i = 0; i < timerFactory->timer; i++) *timerFactory->out << spacer << separator;
    }
    // Remove the current pointer from que queue.
    timerFactory->_pool.pop_back();
#if BenchmarkTimes
    // If we are inside another timer, add the current duration to it's sibilings duration
    if (!timerFactory->_pool.empty()) {
        timerFactory->_pool.back()->siblings_duration += duration;
    }
    if (timerFactory->wholeStats.find(name) != timerFactory->wholeStats.end()) {
        timerFactory->wholeStats[name] += duration;
        timerFactory->wholeUniqueStats[name] += (duration - siblings_duration);
        // If it's the first time we trigger this timer, add it to the timing maps.
    }
        // If not, increase its duration
    else {
        timerFactory->wholeStats[name] = duration;
        timerFactory->wholeUniqueStats[name] = (duration - siblings_duration);
    }
    // Increase the total duration of the program.
    timerFactory->total += (duration - siblings_duration);
#endif
}

#if BenchmarkMemory
int Timer::currentMemoryUsage() {
    std::ifstream proc_status_fhandle("/proc/self/status");
    std::string s;
    int line = 0;
    while (std::getline(proc_status_fhandle, s)) {
        ++line;
        if (!s.compare(0, 6, "VmRSS:")) {
            int value = atoi(&(s.substr(7, s.npos))[0]);
            return value;
        }
    }
    return -1;
}
#endif
//endregion

#endif
