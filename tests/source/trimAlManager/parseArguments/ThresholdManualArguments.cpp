#include "../trimAlManagerIncludes.h"

using namespace catch2Utils::trimAlManager;
using namespace trimAlManagerTestIncludes::thresholds;
using namespace trimAlManagerTestIncludes;

inline void applyCons(
        std::vector<const char *> &args,
        trimAlManager &manager,
        const trimAlManager::argumentReport expectedReport,
        bool expectedHasErrors,
        const testArgumentsCallback &callback = {})
{
    const auto comparator = [](const char *& valueB) {
        return [&valueB](const char *& valueA) {
            return !strcmp(valueA, valueB);
        };
    };

    const auto contains = [&comparator](std::vector<const char*>& container, const char *&& value) -> bool {
        return std::find_if(container.begin(), container.end(), comparator(value)) != container.end();
    };

    auto itct = contains(args, "-ct");
    auto itgt = contains(args, "-gt");
    auto itst = contains(args, "-st");

    // Due to incompatibilities among -ct and -cons providing also -gt or -st
    if ((itct ) && (itgt || itst ))
    {
        expectedHasErrors = true;
    }

    GIVEN("Cons argument") {
        args.push_back("-cons");
        WHEN("No cons value is provided") {
            test_arguments(args, manager,
                           true, trimAlManager::argumentReport::Errored,
                           true, true);
        }
        WHEN("Cons value provided") {
            GIVEN("-10 as cons value") {
                args.push_back("-10");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Errored,
                               true, true);
            }

            for (const char *i : {"0", "50", "100"}) {
                GIVEN(i << " as cons value") {
                    args.push_back(i);
                    test_arguments(args, manager,
                                   true, expectedReport,
                                   true, expectedHasErrors,
                                   [&]() {
                                       if (expectedHasErrors || expectedReport != trimAlManager::argumentReport::Recognized)
                                           return;
                                       THEN("Conservation threshold is set correctly")
                                       {
                                           CHECK(manager.conservationThreshold == std::stoi(i));
                                       }
                                       if (callback) callback();
                                   });
                }
            }

            GIVEN("110 as cons value") {
                args.push_back("110");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Errored,
                               true, true);
            }
        }
    }
}


/**
 * @brief Method to apply a window argument, embedded on a Case (currentCase)
 *          to a current argument call (args) and test it's combinations
 * @param currentCase thresholdCase containing the window argument.
 * Tuple of <CaseName, CaseArgument, WindowArgument>
 * @param args Current argument call to test agains the manager
 * @param manager trimAlManager to perform the calls to parseArguments and processArguments
 * @param expectedReport Expected result of parseArguments ~in the correct call~
 * @param expectedHasErrors Expected result of processArguments ~in the correct call~
 */
inline void applyWindow(
        const thresholdCase &currentCase,
        std::vector<const char *> &args,
        trimAlManager &manager,
        const trimAlManager::argumentReport expectedReport,
        const bool expectedHasErrors,
        const testArgumentsCallback &callback = {}) {

    float tmpFloat = 0;
    const thresholdCase generalCase{"General", "", "-w", tmpFloat, manager.windowSize};

    GIVEN(currentCase.name << " window") {
        args.push_back(currentCase.windowArgument);
        WHEN(currentCase.name << " window value is not provided") {
            test_arguments(args, manager,
                           true, trimAlManager::argumentReport::Errored, true, true);

            applyCons(args, manager, trimAlManager::argumentReport::Errored, true);

            if (currentCase != generalCase) {
                applyWindow(
                        generalCase, args, manager,
                        trimAlManager::argumentReport::Errored, true);
            }
        }
        WHEN(currentCase.name << " window value is provided") {
            for (auto &value : std::vector<reportChecker>
                    {
                            {"-1",   trimAlManager::argumentReport::Errored, true},
                            {"-0.1", trimAlManager::argumentReport::Errored, true},
                            {"1",    expectedReport,                         expectedHasErrors},
                            {"5",    expectedReport,                         expectedHasErrors},
                    }) {
                WHEN(currentCase.name << " window value is " << value.argumentValue) {
                    args.push_back(value.argumentValue);

                    test_arguments(
                            args, manager,
                            true, value.expectedReport,
                            true, value.expectedReturn,
                            [&]() {
                                if (value.expectedReturn || value.expectedReport != trimAlManager::argumentReport::Recognized)
                                    return;
                                THEN(currentCase.name << " window is set correctly")
                                {
                                    CHECK(currentCase.windowRef == std::stoi(value.argumentValue));
                                }
                                if (callback) callback();
                            });

                    applyCons(args, manager,
                              value.expectedReport, value.expectedReturn,
                              [&]() {
                                  if (value.expectedReturn || value.expectedReport != trimAlManager::argumentReport::Recognized)
                                      return;
                                  THEN(currentCase.name << " window is set correctly")
                                  {
                                      CHECK(currentCase.windowRef == std::stoi(value.argumentValue));
                                  }
                                  if (callback) callback();
                              });

                    if (currentCase != generalCase) {
                        applyWindow(
                                generalCase, args, manager,
                                value.expectedReport, true);
                    }
                }
            }
        }
    }
}

/**
 * @brief Method to test a Case:
 * Against passing a threshold argument, not passing a threshold argument.
 *  - Against passing a threshold argument with and without providing values (In ranges).
 * Against passing a window argument, not passing a window argument.
 *  - Against passing a window argument with and without providing values (In ranges).
 * Against passing ~other cases~ windows
 *  - Against passing ~other cases window arguments~ with and without providing values (In ranges).
 * @param currentCase Current thresholdCase to be tested.
 * Tuple of <CaseName, CaseArgument, WindowArgument>
 * @param allCases Vector containing all cases to be tested (for ~other window arguments~)
 * In case allCases contains currentCase, it will be skipped from ~other window arguments~
 * @param args Current call to be tested
 * @param manager trimAlManager to perform the calls to parseArguments and processArguments
 * @param expectedReportOnOK Expected result of parseArguments ~in the correct call~
 * @param expectedErrorsOnOk Expected result of processArguments ~in the correct call~
 */
inline void performSimpleThreshold(
        const thresholdCase &currentCase,
        const std::vector<thresholdCase> &allCases,
        std::vector<const char *> &args,
        trimAlManager &manager,
        const trimAlManager::argumentReport expectedReportOnOK = trimAlManager::argumentReport::Recognized,
        const bool expectedErrorsOnOk = false,
        const testArgumentsCallback &callback = {}) {
    GIVEN(currentCase.name << " Threshold argument") {
        args.push_back(currentCase.thresholdArgument);

        GIVEN("No " << currentCase.name << " Threshold Value") {
            test_arguments(args, manager,
                           true, trimAlManager::argumentReport::Errored, true, true);

            applyWindow(
                    currentCase, args, manager,
                    trimAlManager::argumentReport::Errored, true);

            applyCons(args, manager,
                      trimAlManager::argumentReport::Errored, true);

            for (const thresholdCase &otherCase: allCases) {
                if (otherCase == currentCase) continue;
                applyWindow(
                        otherCase, args, manager,
                        trimAlManager::argumentReport::Errored, true);
            }
        }
        GIVEN(currentCase.name << " Threshold Value") {
            for (auto &value : std::vector<reportChecker>
                    {
                            {"-1.0", trimAlManager::argumentReport::Errored, true},
                            {"-0.1", trimAlManager::argumentReport::Errored, true},
                            {"0.0",  expectedReportOnOK,                     expectedErrorsOnOk},
                            {"0.5",  expectedReportOnOK,                     expectedErrorsOnOk},
                            {"1.0",  expectedReportOnOK,                     expectedErrorsOnOk},
                            {"1.1",  trimAlManager::argumentReport::Errored, true},
                    }) {
                WHEN(currentCase.name << " threshold value is " << value.argumentValue) {
                    args.push_back(value.argumentValue);
                    test_arguments(args, manager,
                                   true, value.expectedReport, true, value.expectedReturn,
                                   [&]() {
                                       if (value.expectedReturn || value.expectedReport != trimAlManager::argumentReport::Recognized)
                                           return;
                                       THEN(currentCase.name << " threshold is set correctly")
                                       {
                                           if (strcmp(currentCase.thresholdArgument, "-gt") == 0)
                                           {
//                                               WARN("Gap threshold is stored differently. This breaks consistency, as all the other stats thresholds are stored straightforward, and no reversed.");
                                               CHECK(currentCase.thresholdRef == 1.F - std::stof(value.argumentValue));
                                           }
                                           else CHECK(currentCase.thresholdRef == std::stof(value.argumentValue));
                                       }
                                       if (callback) callback();
                                   });

                    applyWindow(
                            currentCase, args, manager,
                            value.expectedReport, value.expectedReturn,
                            [&]() {
                                if (value.expectedReturn || value.expectedReport != trimAlManager::argumentReport::Recognized)
                                    return;
                                THEN(currentCase.name << " threshold is set correctly")
                                {
                                    if (strcmp(currentCase.thresholdArgument, "-gt") == 0)
                                    {
//                                        WARN("Gap threshold is stored differently. This breaks consistency, as all the other stats thresholds are stored straightforward, and no reversed.");
                                        CHECK(currentCase.thresholdRef == 1.F - std::stof(value.argumentValue));
                                    }
                                    else CHECK(currentCase.thresholdRef == std::stof(value.argumentValue));
                                }
                                if (callback) callback();
                            });

                    applyCons(args, manager, value.expectedReport, value.expectedReturn,
                              [&]() {
                                  if (value.expectedReturn || value.expectedReport != trimAlManager::argumentReport::Recognized)
                                      return;
                                  THEN(currentCase.name << " threshold is set correctly")
                                  {
                                      if (strcmp(currentCase.thresholdArgument, "-gt") == 0)
                                      {
//                                          WARN("Gap threshold is stored differently. This breaks consistency, as all the other stats thresholds are stored straightforward, and no reversed.");
                                          CHECK(currentCase.thresholdRef == 1.F - std::stof(value.argumentValue));
                                      }
                                      else CHECK(currentCase.thresholdRef == std::stof(value.argumentValue));
                                  }
                                  if (callback) callback();
                              });

                    for (const thresholdCase &otherCase: allCases) {
                        if (otherCase == currentCase) continue;
                        applyWindow(
                                otherCase, args, manager,
                                value.expectedReport, true);
                    }
                }
            }
        }
    }
}

/**
 * @brief Method that performs tests against combinations of two tests:
 * This method contains the same sections as performSimpleThreshold,
 * but on each of those sections, it performs the preparation of the section
 * for the currentCase. Instead of calling to test_arguments or applyWindow,
 * it calls on performSimpleThreshold ~on otherCase~
 *
 * This means the call to performCombination(a,b,...) yields different results
 * than the call to performCombination(b,a,...), and thus, both should be called
 * to achieve the whole combination.
 * @param currentCase Current Case to be tested
 * @param otherCase Case to test combinations against.
 * @param args Arguments to pass to the manager
 * @param manager trimAlManager to perform the calls to parseArguments and processArguments
 */
inline void performCombination(
        const thresholdCase &currentCase,
        const thresholdCase &otherCase,
        std::vector<const char *> &args,
        trimAlManager &manager,
        const testArgumentsCallback &callback = {}) {
    GIVEN(currentCase.name << " along with " << otherCase.name) {
        GIVEN(currentCase.name << " Threshold argument") {
            args.push_back(currentCase.thresholdArgument);

            GIVEN("No " << currentCase.name << " Threshold Value") {
                // Passing {} as allCases prevents it from applying 'other cases windows'
                // This is part of the combination test and thus, must be skipped
                // Otherwise it will expect Errored on cases
                //  where it should expect Recognized
                performSimpleThreshold(
                        otherCase, {},
                        args, manager,
                        trimAlManager::argumentReport::Errored, true);
            }
            GIVEN(currentCase.name << " Threshold Value") {
                for (auto &value : std::vector<reportChecker>
                        {
                                {"-1.0", trimAlManager::argumentReport::Errored,    true},
                                {"-0.1", trimAlManager::argumentReport::Errored,    true},
                                {"0.0",  trimAlManager::argumentReport::Recognized, false},
                                {"0.5",  trimAlManager::argumentReport::Recognized, false},
                                {"1.0",  trimAlManager::argumentReport::Recognized, false},
                                {"1.1",  trimAlManager::argumentReport::Errored,    true},
                        }) {
                    WHEN(currentCase.name << " threshold value is " << value.argumentValue) {
                        args.push_back(value.argumentValue);
                        // Passing {} as allCases prevents it from applying 'other cases windows'
                        // This is part of the combination test and thus, must be skipped
                        // Otherwise it will expect Errored on cases where it should be Recognized
                        performSimpleThreshold(
                                otherCase, {},
                                args, manager,
                                value.expectedReport, value.expectedReturn,
                                [&]() {
                                    if (value.expectedReturn || value.expectedReport != trimAlManager::argumentReport::Recognized)
                                        return;
                                    THEN(currentCase.name << " threshold is set correctly")
                                    {
                                        if (strcmp(currentCase.thresholdArgument, "-gt") == 0)
                                        {
//                                            WARN("Gap threshold is stored differently. This breaks consistency, as all the other stats thresholds are stored straightforward, and no reversed.");
                                            CHECK(currentCase.thresholdRef == 1.F - std::stof(value.argumentValue));
                                        }
                                        else CHECK(currentCase.thresholdRef == std::stof(value.argumentValue));
                                    }
                                    if (callback) callback();
                                });
                    }
                }
            }
        }
    }
}

/**
 * @brief Method to obtain all the combinations against all cases.
 * It calls to performCombination(a,b,...) for each pair possible in cases
 * @param cases Set of Cases to test their combinations
 * @param args Args to be passed to manager
 * @param manager trimAlManager to perform the calls to parseArguments and processArguments
 */
inline void performCombinations(
        const std::vector<thresholdCase> &cases,
        std::vector<const char *> &args,
        trimAlManager &manager) {
    GIVEN("Combinations among thresholds") {
        for (const thresholdCase &currentCase : cases) {
            for (const thresholdCase &otherCase: cases) {
                if (currentCase == otherCase) continue;
                performCombination(currentCase, otherCase, args, manager);
            }
        }
    }
}

// Valgrinded

TEST_CASE("Thresholds argument parse", "[manager][arguments][thresholds][window]") {

    catch2Utils::init();

    trimAlManager manager;

    float tmpFloat = 0;
    const thresholdCase generalCase{"General", "", "-w", tmpFloat, manager.windowSize};

    WHEN("Storing the Gap Threshold")
    {
        WARN(
                "Gap threshold is stored differently to the rest.\n"
                "While all the stats thresholds are stored straightforward: (CMD-Input:\"1.0\" => StoredValue:1.0F)\n"
                "Gap threshold is stored reversed: (CMD-Input:\"1.0\" => StoredValue:0.0F)\n"
                "This breaks the consistency among the code, as it can be seen on this example:\n"
                "\t- While all the tests could be generic, in the case of checking threshold values,\n"
                "\t\twe have to check splicitly if it is the gaps one.\n"
                "This is a code smell.");
    }

    GIVEN("Alignment as input") {
        auto thresholdCases = getThresholdCases(manager, false);
        std::vector<const char *> args =
                {"", "-in", "../dataset/example.091.AA.strNOG.ENOG411BWBU.fasta"};

        for (auto & x : thresholdCases)
            performSimpleThreshold(x, thresholdCases, args, manager);

        GIVEN("No threshold") {
            for (auto & x : thresholdCases)
                applyWindow(x, args, manager, trimAlManager::argumentReport::Recognized, true);
            applyWindow(generalCase, args, manager, trimAlManager::argumentReport::Recognized, true);
        }

        performCombinations(thresholdCases, args, manager);
    }

    GIVEN("Compareset as input") {
        std::vector<const char *> args =
                {"", "-compareset", "../dataset/alignments_comparison.1"};

        auto thresholdCases = getThresholdCases(manager, true);

        for (auto & x : thresholdCases)
            performSimpleThreshold(x, thresholdCases, args, manager);

        GIVEN("No threshold") {
            for (auto & x : thresholdCases)
                applyWindow(x, args, manager, trimAlManager::argumentReport::Recognized, true);
            applyWindow(generalCase, args, manager, trimAlManager::argumentReport::Recognized, true);
        }

        performCombinations(thresholdCases, args, manager);
    }

}