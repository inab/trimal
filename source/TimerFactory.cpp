#include "TimerFactory.h"
#include <utility>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <fstream>
#include <reportsystem.h>

#if TimingReport
TimerFactory timerFactory = TimerFactory("");

Timer TimerFactory::getTimer(std::string name) {
    Timer timer = Timer(std::move(name), this);
    _pool.push_back(&timer);
    return timer; // NOLINT
}

void TimerFactory::reportTotal() {
    *out << "\n";
    int maxSize = -1;
    for (auto &&value : wholeStats)
    {
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
            << " |" << std::setw(maxSize - 5) << "Key " << " | "
            << " |" << std::setw(21) << "Duration " << " | "
            << " |" << std::setw(21) << "Rel Duration " << " | "
            << " |" << std::setw(21) << "Total Duration " << " | "
            << " |" << std::setw(21) << "Rel Total Duration " << " | "
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

    for (auto &&value : wholeStats) {
        *out << std::setw(maxSize - 2) << value.first.substr(0, value.first.find("(")) << "  "
             << std::setw(21) << wholeUniqueStats[value.first].count() << " ms "
             << std::fixed
             << std::setw(21) << wholeUniqueStats[value.first].count() * 100.F / total.count() << " %  "
             << std::defaultfloat
             << std::setw(21) << value.second.count() << " ms "
             << std::fixed
             << std::setw(21) << value.second.count() * 100.F / total.count() << " %  "
             << std::defaultfloat
             << "\n"
             << std::setw(maxSize - 2) << value.first.substr(value.first.find("(")) << "\n"
                ;
    }
}

TimerFactory::~TimerFactory() {
    reportTotal();
    if (out->rdbuf() != std::cout.rdbuf() &&
        out->rdbuf() != std::cerr.rdbuf() &&
        out->rdbuf() != std::clog.rdbuf())
        delete out;
}

Timer::Timer(std::string name, TimerFactory *timerFactory) :
        name(std::move(name)), timerFactory(timerFactory), siblings_duration(0) {

    start = std::chrono::high_resolution_clock::now();

    if (!timerFactory->_pool.empty())
    {
        timerFactory->_pool.back()->childCount += 1;
    }

    if (timerFactory->lastTimer != this->name) {
        *timerFactory->out << "\n";
        for (int i = 0; i < timerFactory->timer - 1; i++) {
            *timerFactory->out << "|";
        }
        *timerFactory->out << "|----> 1.";
        for (Timer* timer : timerFactory->_pool)
        {
            *timerFactory->out << timer->childCount << ".";
        }

        *timerFactory->out << "  " << this->name << " " << currentMemoryUsage() << " kb";
    }
    timerFactory->timer++;
}

int Timer::currentMemoryUsage() {
    std::ifstream proc_status_fhandle("/proc/self/status");
    std::string s;
    int line=0;
    while(std::getline(proc_status_fhandle, s)){
        ++line;
        if (!s.compare(0, 6, "VmRSS:"))
        {
            int value = atoi(&(s.substr(7, s.npos))[0]);
//            debug << s << "\n";
            return value;
        }
    }
}

Timer::~Timer() {
    timerFactory->timer--;
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);
    if (timerFactory->lastTimer == name) {
        timerFactory->accTimerCount += 1;
        timerFactory->accTimer += duration;
        timerFactory->out->seekp(timerFactory->cursorPos);
        for (int i = 0; i < timerFactory->timer; i++) *timerFactory->out << "|";
        *timerFactory->out << "`·~ " << name
                           << " Repetitions: " << timerFactory->accTimerCount
                           << " Average: " << (timerFactory->accTimer.count() / timerFactory->accTimerCount) << "ms."
                           << " Acc " << timerFactory->accTimer.count() << "ms. "
                            << currentMemoryUsage() << " kb";

    } else {
        timerFactory->accTimerCount = 1;
        timerFactory->accTimer = duration;
        timerFactory->lastTimer = name;

        *timerFactory->out << "\n";
        timerFactory->cursorPos = timerFactory->out->tellp();
        for (int i = 0; i < timerFactory->timer; i++) *timerFactory->out << "|";
        *timerFactory->out << "`·~ " << name
                           << " Duration: " << (duration - siblings_duration).count() << "ms."
                           << " Acc " << duration.count() << "ms. "
                           << currentMemoryUsage() << " kb";
        *timerFactory->out << "\n";
        for (int i = 0; i < timerFactory->timer; i++) *timerFactory->out << "|";
    }
    timerFactory->_pool.pop_back();
    if (!timerFactory->_pool.empty()) {
        timerFactory->_pool.back()->siblings_duration += duration;
    }
    if (timerFactory->wholeStats.find(name) != timerFactory->wholeStats.end()) {
        timerFactory->wholeStats[name] += duration;
        timerFactory->wholeUniqueStats[name] += (duration - siblings_duration);
    } else {
        timerFactory->wholeStats[name] = duration;
        timerFactory->wholeUniqueStats[name] = (duration - siblings_duration);
    }
    timerFactory->total += (duration - siblings_duration);
}

#endif
