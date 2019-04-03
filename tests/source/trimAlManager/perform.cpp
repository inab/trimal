#include "trimAlManagerIncludes.h"

using namespace catch2Utils::trimAlManager;
using namespace trimAlManagerTestIncludes;

void testPerform(
        trimAlManager & manager,
        int expectedResult,
        const std::function<void(void)>& callback = {})
{
    int performResult;
    std::string capturedOutput = catch2Utils::capturingLambda(
            [&performResult, &manager](){
                performResult = manager.perform();
            });
    INFO(capturedOutput);
    CHECK(performResult == expectedResult);

    if (callback) callback();
}

void applyAutomatedMethod(
        trimAlManager & manager,
        int expectedReturn,
        const std::function<void(void)>& callback = {})
{
    GIVEN("Cleaning method")
    {
        for(auto & method : getAutomatedMethods(manager))
        {
            method.reference = true;
            WHEN("Applying " << method.argument)
            {
                THEN("Performs returns non-errored")
                {
                    testPerform(manager, expectedReturn);
                }
                if (callback) callback();
            }
        }
    }
}

void applyStat(
        trimAlManager & manager,
        int expectedReturn,
        const std::function<void(void)>& callback = {})
{
    GIVEN("Requested Stat Output -sXX | -sident | -soverlap")
    {
        for (auto& arg : getStatsArguments(manager, true, false))
        {
            arg.reference = true;
            WHEN("Requesting " << &arg.argument[1])
            {
                THEN("Performs returns non-errored")
                {
                    testPerform(manager, expectedReturn);
                }
                if (callback) callback();
            }
        }
    }
}


//region Providing compareset arguments
TEST_CASE("Perform", "[manager][perform]") {
    catch2Utils::init();

    trimAlManager manager;

    GIVEN("Input alignment")
    {
        manager.origAlig = manager.getFormatManager().loadAlignment(alignmentPath);
        testPerform(manager, 0);

        WHEN("Previous steps have errored")
        {
            manager.appearErrors = true;
            THEN("Perform returns errored")
            {
                testPerform(manager, 1);
            }
        }

        applyAutomatedMethod(manager, 0, [&manager](){ applyStat(manager, 0); });
        applyStat(manager, 0);
    }
}