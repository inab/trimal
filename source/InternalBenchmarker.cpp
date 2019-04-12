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

#include "InternalBenchmarker.h"

#if TimingReport
    // Global variable timerFactory.
    InternalBenchmarker timerFactory = InternalBenchmarker(""); // NOLINT

    //region InternalBenchmarker
    BenchmarkSnapshot InternalBenchmarker::getTimer(std::string name) {
        // Create a new timer
        BenchmarkSnapshot timer = BenchmarkSnapshot(std::move(name), this);
        // Store its pointer in the current queue
        _pool.emplace_back(&timer);
        // Return the timer
        return timer; // NOLINT
    }

    void InternalBenchmarker::reportTotal() {
        double duration = std::chrono::duration_cast<std::chrono::milliseconds>
                (std::chrono::high_resolution_clock::now() - startPoint).count();

        *timerFactory.out
                << "\n\nTotal Duration: " << duration << "ms\n"
                "Number of different counters: " << totalDuration.size();

        duration /= 100.0F;

        for (auto &total: totalDuration) {
            *timerFactory.out << "\n\n"
                            << total.first << "\n"
                            << std::setw(25) << std::left << "\tTotal Time Elapsed:"
                            << total.second << "("
                            << total.second / duration << "%)\n"
                            << std::setw(25) << std::left << "\tRelative Time Elapsed:"
                            << (total.second - relativeDuration[total.first]) << "("
                            << (total.second - relativeDuration[total.first]) / duration << "%)";
        }
        *timerFactory.out << "\n";
    }

    InternalBenchmarker::~InternalBenchmarker() {

        // Output the chart
        reportTotal();
        // Delete the out pointer in case it's not a global one (cout, cerr, clog)
        if (out->rdbuf() != std::cout.rdbuf() &&
            out->rdbuf() != std::cerr.rdbuf() &&
            out->rdbuf() != std::clog.rdbuf())
            delete out;
    }

    bool InternalBenchmarker::checkOutputParameter(int argc, char **argv) {
        for (int i = 1; i < argc; i++) {
            if ((!strcmp(argv[i], "-timetrackerout")) && ((i) + 1 != argc)) {
                size_t argumentLength = strlen(argv[++i]);
                auto *out = new char[argumentLength + 1];
                strcpy(out, argv[i]);
                this->SetOutput(out);
                std::cout << "Time Tracker Output changed to " << out << std::endl;
                delete[] out;
                return true;
            }
        }
        this->SetOutput("/dev/null");
        return false;
    }
    //endregion

    //region BenchmarkSnapshot

    std::string BenchmarkSnapshot::separator = std::string("│ "); // NOLINT
    std::string BenchmarkSnapshot::spacer = std::string("  "); // NOLINT


    BenchmarkSnapshot::BenchmarkSnapshot(std::string name, InternalBenchmarker *timerFactory) :
            name(std::move(name)),
            timerFactory(timerFactory) {
        start = std::chrono::high_resolution_clock::now();
        timerFactory->levels.emplace_back(timerFactory->lastLevel + 1);
        timerFactory->lastLevel = 0;
        if (!timerFactory->_pool.empty() && timerFactory->lastName == this->name) {
            timerFactory->out->seekp(timerFactory->cursorPos);
            timerFactory->timesStepRepeated++;
        } else
            timerFactory->timesStepRepeated = 1;

        timerFactory->cursorPos = timerFactory->out->tellp();
        for (int i = 1; i < timerFactory->levels.size(); i++)
            *timerFactory->out << separator;
        *timerFactory->out << "[STT] ";
        for (int &level : timerFactory->levels)
            *timerFactory->out << level << ".";
        *timerFactory->out << " " << this->name
                        #if BenchmarkMemory
                        << currentMemoryUsage() << " kb"
                        #endif
                        << "\n";
    }

    BenchmarkSnapshot::~BenchmarkSnapshot() {
        timerFactory->lastLevel = timerFactory->levels.back();
        timerFactory->_pool.pop_back();

        for (int i = 1; i < timerFactory->levels.size(); i++)
            *timerFactory->out << separator;
        *timerFactory->out << "[END] ";
        for (int &level : timerFactory->levels) {
            *timerFactory->out << level << ".";
        }

        long duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
        if (!timerFactory->_pool.empty())
            timerFactory->relativeDuration[timerFactory->_pool.back()->name] += duration;
        timerFactory->totalDuration[name] += duration;
    
    #if BenchmarkTimes
            *timerFactory->out /*<< " " << this->name*/ << " duration: "
                                                    << duration
    #endif

    #if BenchmarkMemory
        << " " << currentMemoryUsage() << " kb"
    #endif
                ;
        if (timerFactory->timesStepRepeated > 1)
            *timerFactory->out << ". " << timerFactory->timesStepRepeated << " repetitions";
        *timerFactory->out << "\n";

        timerFactory->levels.pop_back();
        timerFactory->lastName = this->name;
        timerFactory->firstStepStart = start;
    }

    #if BenchmarkMemory

        int BenchmarkSnapshot::currentMemoryUsage() {
            static std::ifstream proc_status_fhandle("/proc/self/status");
            proc_status_fhandle.seekg(0);
            std::string s;
            while (std::getline(proc_status_fhandle, s)) {
                if (!s.compare(0, 6, "VmRSS:")) {
                    int value = atoi(&(s.substr(7, std::string::npos))[0]);
                    return value;
                }
            }
            return -1;
        }

    #endif
    //endregion

#endif
