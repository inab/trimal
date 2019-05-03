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
#include "reportsystem.h"

reporting::reportManager debug = reporting::reportManager();

void reporting::reportManager::PrintCodesAndMessages() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void reporting::reportManager::PrintCodesAndMessages() ");
    switch (Level) {
        case VerboseLevel::NONE:
            std::cout << "[VerboseLevel] None" << std::endl;
            break;
        case VerboseLevel::INFO:
            std::cout << "[VerboseLevel] Info" << std::endl;
            break;
        case VerboseLevel::WARNING:
            std::cout << "[VerboseLevel] Warning" << std::endl;
            break;
        case VerboseLevel::ERROR:
            std::cout << "[VerboseLevel] Error" << std::endl;
            break;
    }

    for (int i = 1; i < InfoCode::__MAXINFO; i++) {
        report((InfoCode) i);
    }

    for (int i = 1; i < WarningCode::__MAXWARNING; i++) {
        report((WarningCode) i);
    }

    for (int i = 1; i < ErrorCode::__MAXERROR; i++) {
        report((ErrorCode) i);
    }
}

void reporting::reportManager::report(ErrorCode message, std::string *vars) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void reporting::reportManager::report(ErrorCode message, std::string *vars) ");
    if (Level > VerboseLevel::ERROR) {
            delete[] vars;
    } else {
        if (vars == nullptr) {
            std::cerr << "[ERROR " << std::setw(3) << std::setfill('0') << message << "] " << ErrorMessages.at(message) << std::endl << std::setfill(' ');
            return;
        }

        std::string s(ErrorMessages.at(message));

        std::string FindWord = "[tag]";

        int counter = 0;

        std::size_t index;
        while ((index = s.find(FindWord)) != std::string::npos)
            s.replace(index, FindWord.length(), vars[counter++]);

        std::cerr << "[ERROR " << std::setw(3) << std::setfill('0') << message << "] " << s << std::endl << std::setfill(' ');

        delete[] vars;
    }
}

void reporting::reportManager::report(ErrorCode message, const char *vars) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void reporting::reportManager::report(ErrorCode message, char *vars) ");
    if (Level > VerboseLevel::ERROR) return;

    if (vars == nullptr) {
        std::cerr << "[ERROR " << std::setw(3) << std::setfill('0') << message << "] " << ErrorMessages.at(message) << std::endl << std::setfill(' ');
        return;
    }

    std::string s(ErrorMessages.at(message));

    std::string FindWord = "[tag]";

    std::string Vars = vars;

    std::size_t index;
    while ((index = s.find(FindWord)) != std::string::npos)
        s.replace(index, FindWord.length(), Vars);

    std::cerr << "[ERROR " << std::setw(3) << std::setfill('0') << message << "] " << s << std::endl << std::setfill(' ');

}

void reporting::reportManager::report(WarningCode message, std::string *vars) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void reporting::reportManager::report(WarningCode message, std::string *vars) ");
    if (Level > VerboseLevel::WARNING) {
            delete[] vars;
    } else {
        if (vars == nullptr) {
            std::cout << "[WARNING " << std::setw(3) << std::setfill('0') << message << "] " << WarningMessages.at(message) << std::endl << std::setfill(' ');
            return;
        }

        std::string s(WarningMessages.at(message));

        std::string FindWord = "[tag]";

        int counter = 0;

        std::size_t index;
        while ((index = s.find(FindWord)) != std::string::npos)
            s.replace(index, FindWord.length(), vars[counter++]);

        std::cout << "[WARNING " << std::setw(3) << std::setfill('0') << message << "] " << s << std::endl << std::setfill(' ');

        delete[] vars;
    }
}

void reporting::reportManager::report(WarningCode message, const char *vars) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void reporting::reportManager::report(WarningCode message, char *vars) ");
    if (Level > VerboseLevel::WARNING) return;

    if (vars == nullptr) {
        std::cout << "[WARNING " << std::setw(3) << std::setfill('0') << message << "] " << WarningMessages.at(message) << std::endl << std::setfill(' ');
        return;
    }

    std::string s(WarningMessages.at(message));

    std::string FindWord = "[tag]";

    std::string Vars = vars;

    std::size_t index;
    while ((index = s.find(FindWord)) != std::string::npos)
        s.replace(index, FindWord.length(), Vars);

    std::cout << "[WARNING " << std::setw(3) << std::setfill('0') << message << "] " << s << std::endl << std::setfill(' ');

}

void reporting::reportManager::report(InfoCode message, std::string *vars) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void reporting::reportManager::report(InfoCode message, std::string *vars) ");
    if (Level > VerboseLevel::INFO) {
            delete[] vars;
    } else {
        if (vars == nullptr) {
            std::cout << "[INFO " << std::setw(3) << std::setfill('0') << message << "] " << InfoMessages.at(message) << std::endl << std::setfill(' ');
            return;
        }

        std::string s(InfoMessages.at(message));

        std::string FindWord = "[tag]";

        int counter = 0;

        std::size_t index;
        while ((index = s.find(FindWord)) != std::string::npos)
            s.replace(index, FindWord.length(), vars[counter++]);

        std::cout << "[INFO " << std::setw(3) << std::setfill('0') << message << "] " << s << std::endl << std::setfill(' ');

        delete[] vars;
    }
}

void reporting::reportManager::report(InfoCode message, const char *vars) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void reporting::reportManager::report(InfoCode message, char *vars) ");
    if (Level > VerboseLevel::INFO) return;

    if (vars == nullptr) {
        std::cout << "[INFO " << std::setw(3) << std::setfill('0') << message << "] " << InfoMessages.at(message) << std::endl << std::setfill(' ');
        return;
    }

    std::string s(InfoMessages.at(message));

    std::string FindWord = "[tag]";

    std::string Vars = vars;

    std::size_t index;
    while ((index = s.find(FindWord)) != std::string::npos)
        s.replace(index, FindWord.length(), Vars);

    std::cout << "[INFO " << std::setw(3) << std::setfill('0') << message << "] " << s << std::endl << std::setfill(' ');

}
