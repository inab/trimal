#include "../trimAlManagerIncludes.h"

using namespace catch2Utils::trimAlManager;
using namespace trimAlManagerTestIncludes;

// Valgrinded

TEST_CASE("Automated method argument parse", "[manager][arguments][automated]") {
    catch2Utils::init();

    trimAlManager manager;

    auto automatedMethods = getAutomatedMethods(manager);

    for (auto & method : automatedMethods) {
        GIVEN(method.argument) {
            AND_GIVEN("No Input") {
                test_arguments(
                        {"", method.argument},
                        manager,
                        true, trimAlManager::argumentReport::Errored,
                        false, true);
            }
            AND_GIVEN("Correct input") {
                test_arguments(
                        {"",  method.argument,
                         "-in", alignmentPath},
                        manager,
                        true, trimAlManager::argumentReport::Recognized,
                        true, false,
                        [&method, &automatedMethods](){
                            THEN("The variable in the manager is set correctly")
                            {
                                REQUIRE(method.reference);
                            }
                            THEN("Other methods stay false")
                            {
                                for(auto & otherMethod : automatedMethods)
                                {
                                    if (otherMethod == method) continue;
                                    CHECK(!otherMethod.reference);
                                }
                            }
                        });
            }
        }
    }

    GIVEN("More than one automated method") {
        for (auto &method1 : automatedMethods) {
            for (auto &method2 : automatedMethods) {
                if (method1 == method2) continue;
                GIVEN(method1.argument << " and " << method2.argument) {
                    AND_GIVEN("Input") {
                        test_arguments(
                                {"", "-in", alignmentPath,
                                 method1.argument, method2.argument},
                                manager,
                                true, trimAlManager::argumentReport::Recognized,
                                true, true);
                    }
                    AND_GIVEN("No Input") {
                        test_arguments(
                                {"", method1.argument, method2.argument},
                                manager,
                                true, trimAlManager::argumentReport::Errored,
                                true, true);
                    }
                }
            }
        }
    }
}
