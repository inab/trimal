#include "../trimAlManagerIncludes.h"

using namespace catch2Utils::trimAlManager;
using namespace trimAlManagerTestIncludes;

// Valgrinded

//region Providing report output
TEST_CASE("Reporting argument parse", "[manager][arguments]") {

    catch2Utils::init();

    trimAlManager manager;

    auto reportArguments = getReportArguments(manager);

    for (auto& method : reportArguments) {
        GIVEN(method.argument+1) {
            WHEN("No input is provided") {
                WHEN("Output path is provided") {
                    WHEN("Trimming method is not provided") {
                        test_arguments(
                                {"", method.argument, "./tmp"},
                                manager,
                                trimAlManager::argumentReport::Wrong,
                                true);
                    }
                    WHEN("Trimming method is provided") {
                        test_arguments(
                                {"",
                                 "-strict",
                                 method.argument, "./tmp"},
                                manager,
                                trimAlManager::argumentReport::Wrong,
                                true);
                    }
                }

                WHEN("Output path is not provided") {
                    WHEN("Trimming method is not provided") {
                        test_arguments(
                                {"",
                                 method.argument},
                                manager,
                                trimAlManager::argumentReport::Wrong,
                                true);
                    }
                    WHEN("Trimming method is provided") {
                        test_arguments(
                                {"",
                                 "-strict",
                                 method.argument},
                                manager,
                                trimAlManager::argumentReport::Wrong,
                                true);
                    }
                }
            }

            WHEN("Input is provided") {
                WHEN("Output path is provided") {
                    WHEN("Trimming method is not provided") {
                        test_arguments(
                                {"",
                                 "-in", alignmentPath,
                                 method.argument, "./tmp"},
                                manager,
                                trimAlManager::argumentReport::Recognized,
                                true);
                    }
                    WHEN("Trimming method is provided") {
                        test_arguments(
                                {"",
                                 "-in", alignmentPath,
                                 "-strict",
                                 method.argument, "./tmp"},
                                manager,
                                trimAlManager::argumentReport::Recognized,
                                false,
                                [&method](){
                                    THEN("Values are set correctly for " << method.argument+1)
                                    {
                                        CHECK(strcmp(*method.pointerToPointerVar, "./tmp") == 0);
                                    }
                                });
                    }
                }

                WHEN("Output path is not provided") {
                    WHEN("Trimming method is not provided") {
                        test_arguments(
                                {"",
                                 "-in", alignmentPath,
                                 method.argument},
                                manager,
                                trimAlManager::argumentReport::Wrong,
                                true);
                    }
                    WHEN("Trimming method is provided") {
                        test_arguments(
                                {"",
                                 "-in", alignmentPath,
                                 "-strict",
                                 method.argument},
                                manager,
                                trimAlManager::argumentReport::Wrong,
                                true);
                    }
                }
            }
        }
    }

    WHEN("More than one report argument has been provided") {
        for (auto & method1 : reportArguments) {
            for (auto & method2 : reportArguments) {
                if (method1 == method2) continue;

                GIVEN((method1.argument+1) << " and " << (method2.argument+1)) {
                    WHEN("Both arguments point to different outputs") {
                        test_arguments(
                                {"",
                                 "-in", alignmentPath,
                                 "-strict",
                                 method1.argument, "./tmp1", method2.argument, "./tmp2"},
                                manager,
                                trimAlManager::argumentReport::Recognized,
                                false,
                                [&method1, &method2](){
                                    THEN("Values are set correctly for " << method1.argument+1)
                                    {
                                        CHECK(strcmp(*method1.pointerToPointerVar, "./tmp1") == 0);
                                    }
                                    THEN("Values are set correctly for " << method2.argument+1)
                                    {
                                        CHECK(strcmp(*method2.pointerToPointerVar, "./tmp2") == 0);
                                    }
                                });
                    }

                    WHEN("Both arguments point to same outputs") {
                        test_arguments(
                                {"",
                                 "-in", alignmentPath,
                                 "-strict",
                                 method1.argument, "./tmp1", method2.argument, "./tmp1"},
                                manager,
                                trimAlManager::argumentReport::Recognized,
                                true);
                    }
                }
            }
        }
    }
}
//endregion