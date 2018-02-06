//
// Created by bioinfo on 19/01/18.
//

#ifndef MTRIMAL_TIMERFACTORY_H
#define MTRIMAL_TIMERFACTORY_H

#define TimingReport true

#if TimingReport

#define StartTiming(name) auto macroCustomTimer = timerFactory.getTimer(name)
#define SetTimingOstream(_stream) timerFactory.out = _stream;
#include <chrono>
#include <string>
#include <stack>
#include <iostream>
#include <map>

class TimerFactory;

class Timer {
    friend TimerFactory;
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::milliseconds siblings_duration;
    std::string name;


public:
    explicit Timer(std::string name, TimerFactory *timerFactory);;
    TimerFactory* timerFactory;
    ~Timer();;
};

class TimerFactory {
    friend Timer;
public:
    ~TimerFactory();
    Timer getTimer(std::string name);;

    explicit TimerFactory(std::ostream * out)
            : out(out), lastTimer(""), accTimer(0),
              wholeStats(std::map<std::string, std::chrono::milliseconds>()),
              wholeUniqueStats(std::map<std::string, std::chrono::milliseconds>()),
              total(std::chrono::milliseconds(0)) {};

    void reportTotal();

    std::ostream * out;
private:

    int timer = 0;
    std::stack<Timer *> _pool;
    std::string lastTimer;
    std::chrono::milliseconds accTimer;
    int accTimerCount{0};

    std::map<std::string, std::chrono::milliseconds> wholeStats;
    std::map<std::string, std::chrono::milliseconds> wholeUniqueStats;
    std::chrono::milliseconds total;
    std::streampos cursorPos{0};
};

extern TimerFactory timerFactory;
#else
#define StartTiming(name) {}
#define SetTimingOstream(_stream) {}
#endif
#endif //MTRIMAL_TIMERFACTORY_H
