//
// Created by vfernandez on 08/03/19.
//

#include "ostream"
#include "trimalManager.h"
#include "trimAlManagerArgumentReport.h"
#include "../catch.hpp"
#include "ScopedRedirect.h"

std::ostream &operator<<(std::ostream &os, const trimAlManager::argumentReport &item) {

    switch (item) {
        case trimAlManager::NotRecognized:
            os << "NotRecognized";
            break;
        case trimAlManager::Recognized:
            os << "Recognized";
            break;
        case trimAlManager::Errored:
            os << "Errored";
            break;
        case trimAlManager::Final:
            os << "Final";
            break;
        default:
            os << static_cast<std::underlying_type<trimAlManager::argumentReport>::type>(item);
            break;
    }

    return os;
}

const std::string catch2Utils::trimAlManager::getMessageFor(catch2Utils::trimAlManager::parseArgumentsEnum val) {
    return parseArgumentsMessages.at(val);
}

const std::string catch2Utils::trimAlManager::getMessageFor(::trimAlManager::argumentReport val) {

    switch(val)
    {
        case ::trimAlManager::argumentReport::Errored:
            return parseArgumentsMessages.at(parseArgumentsErrored);
        case ::trimAlManager::argumentReport::Recognized:
            return parseArgumentsMessages.at(parseArgumentsCorrect);
        case ::trimAlManager::argumentReport::Final:
            return parseArgumentsMessages.at(parseArgumentsFinal);
        default:
            return getMessageFor((catch2Utils::trimAlManager::parseArgumentsEnum)val);
    }
}

void catch2Utils::trimAlManager::test_arguments(
        std::vector<const char *> args,
        ::trimAlManager &manager,
        bool parseArguments,
        ::trimAlManager::argumentReport expectedParseArgumentsResult,
        bool processArguments,
        bool expectingProcessArgumentsError,
        const testArgumentsCallback& callback)
{

    // Start the section for the expected report.
    //  Get the correct Section Name
    THEN(getMessageFor(expectedParseArgumentsResult)) {

        // If the first check should be done
        if (parseArguments) {

            // Get the result from the call
            ::trimAlManager::argumentReport report;

            // Capture the output of the execution with the capturingLambda
            auto output = catch2Utils::capturingLambda([&]() -> void {
                // Perform the processArguments call
                report = (::trimAlManager::argumentReport)
                        manager.parseArguments(
                                // With the correct number on argc
                                args.size(),
                                // And the values on argv
                                const_cast<char **>(args.data()));
            });

            // "\u202F" => NARROW NO-BREAK SPACE -> Multiple \n are ignored on CATCH.
            // To trick the system and allow to multiple line breaks,
            //  we use an 'invisible' character to allow separate blocks with
            //  several line breaks
            output += "\n\u202F\n";

            // Report the captured output
            INFO(output);
            CAPTURE(args);
            // If the result of parseArguments is not the expected
            CHECK(report == expectedParseArgumentsResult);

            output.clear();

            // If the second check should be done
            if (processArguments) {
                // Get the correct Section Name
                AND_THEN(getMessageFor( expectingProcessArgumentsError ?
                        processArgumentsErrored : processArgumentsNonErrored)) {
                    bool hasErrors;
                    output = catch2Utils::capturingLambda([&]() -> void {
                        hasErrors = manager.processArguments(
                                const_cast<char **>(args.data()));
                    });
                    // "\u202F" => NARROW NO-BREAK SPACE -> Multiple \n are ignored on CATCH.
                    // To trick the system and allow to multiple line breaks,
                    //  we use an 'invisible' character to allow separate blocks with
                    //  several line breaks
                    output = "\n\u202F\n  " + output + "\n\u202F\n";

                    INFO(output);
                    CHECK(hasErrors == expectingProcessArgumentsError);
                }
            } else
                WARN("Second test skipped");
        } else
            FAIL("First check skipped");

        if (callback)
            callback();
    }
}

void catch2Utils::trimAlManager::test_arguments(
        std::vector<const char *> args,
        ::trimAlManager &manager,
        testArgumentsConfig testConfig,
        const testArgumentsCallback& callback)
{
    test_arguments(
            std::move(args), manager,
            std::get<0>(testConfig), std::get<1>(testConfig),
            std::get<2>(testConfig), std::get<3>(testConfig),
                    std::move(callback));
}