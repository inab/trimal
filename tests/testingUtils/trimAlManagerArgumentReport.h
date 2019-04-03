//
// Created by vfernandez on 08/03/19.
//

#ifndef TRIMAL_TRIMALMANAGERARGUMENTREPORT_H
#define TRIMAL_TRIMALMANAGERARGUMENTREPORT_H

#include "map"
#include "../../include/trimalManager.h"
#include "functional"

std::ostream &operator<<(std::ostream &os, const trimAlManager::argumentReport &item);

namespace catch2Utils
{
    namespace trimAlManager
    {
        enum parseArgumentsEnum
        {
            parseArgumentsErrored = 0,
            parseArgumentsCorrect = 1,
            parseArgumentsFinal   = 2,

            processArgumentsErrored = 3,
            processArgumentsNonErrored = 4
        };

        const std::map<parseArgumentsEnum , const std::string> parseArgumentsMessages
        {
            {parseArgumentsErrored,             "ParseArguments returns Errored"},
            {parseArgumentsCorrect,             "ParseArguments returns Correct"},
            {parseArgumentsFinal,               "ParseArguments returns Final"},

            {processArgumentsErrored,           "ProcessArguments returns errored"},
            {processArgumentsNonErrored,        "ProcessArguments returns non-errored"}
        };

        const std::string getMessageFor(parseArgumentsEnum val);
        const std::string getMessageFor(::trimAlManager::argumentReport argument);

        typedef std::tuple<bool,::trimAlManager::argumentReport, bool, bool> testArgumentsConfig;
        typedef std::function<void(void)> testArgumentsCallback;

        void test_arguments(
                std::vector<const char *> args,
                ::trimAlManager &manager,
                bool parseArguments,
                ::trimAlManager::argumentReport expectedParseArgumentsResult,
                bool processArguments,
                bool expectingProcessArgumentsError,
                const testArgumentsCallback& callback = {});

        void test_arguments(
                std::vector<const char *> args,
                ::trimAlManager &manager,
                testArgumentsConfig testConfig,
                const testArgumentsCallback& callback = {});
    }
}

#endif //TRIMAL_TRIMALMANAGERARGUMENTREPORT_H
