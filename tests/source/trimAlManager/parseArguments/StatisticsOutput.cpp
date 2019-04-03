#include "../trimAlManagerIncludes.h"

using namespace catch2Utils::trimAlManager;
using namespace trimAlManagerTestIncludes;

/**
 * @brief Method to call to test_arguments with all the combinations of arguments
 * present in a vector of arguments (argumentsToCombine),
 * appending the combination to a 'prefix argument call list' (argumentPrefix).
 *
 * The method will:
 *
 *  - If explicitly requested, (performCallback) perform the call
 *  to test_arguments, using argumentPrefix as argv.
 *
 *  - If iterator (i) points to a valid position it will call this method two times
 *  with an increased value of the iterator (i):
 *
 *      - First, with the new argument added to argumentPrefix, requesting execution
 *      of test_callback. (As is the first time this combination has been seen).
 *
 *      - Second, without the new argument added, NOT requesting execution of test_arguments
 *      (As it would do the same as the current execution).
 *      This skips the current argument, to obtain combinations ~without~ that argument.
 *
 * @param argumentsToCombine Vector of arguments to obtain the combination of them.
 * @param argumentPrefix Initial argument list to pass to the manager.
 * On the first call, it contains the prefix common to all combinations.
 * On subsequent calls, it contains the prefix, plus the added arguments from argumentsToCombine
 * @param manager trimAlManager to perform the call to test_arguments
 * @param performCallback Whether this call to the method should trigger the call to test_arguments
 * @param i argumentsToCombine iterator variable
 */
void perform(
        const std::vector<argumentReference> &argumentsToCombine,
        const std::vector<const char *> &argumentPrefix,
        trimAlManager &manager,
        const bool performCallback = false,
        const ulong i = 0UL,
        std::vector<argumentReference> argumentsToCheck = {})
{
    // Only perform the test when explicitly requested
    if (performCallback)
        test_arguments(
            argumentPrefix, manager,
            true, trimAlManager::argumentReport::Recognized, true, false,
            [&argumentsToCheck]()
            {
                THEN("Variables are set correctly in the manager")
                {
                    for(auto & argumentPair : argumentsToCheck)
                    {
                        THEN(&argumentPair.argument[1] << " is set correctly.")
                            CHECK(argumentPair.reference);
                    }
                }
            });

    // If we haven't arrived to the end of the vector, there are 2 more cases
    if (i < argumentsToCombine.size())
    {
        GIVEN(argumentsToCombine[i].argument)
        {
            // A case with a new argument added
            std::vector<const char *> newCall = argumentPrefix;
            newCall.push_back(argumentsToCombine[i].argument);
            argumentsToCheck.push_back(argumentsToCombine[i]);
            perform(argumentsToCombine, newCall, manager, true, i + 1, argumentsToCheck);
        }

        // A case with no new argument added, but with a new iterator value
        perform(argumentsToCombine, argumentPrefix, manager, false, i + 1, argumentsToCheck);
    }
}

// Valgrinded

TEST_CASE("Statistics report argument parse", "[manager][arguments]") {

    using namespace trimAlManagerTestIncludes;
    // Tests initialization
    trimAlManager manager;

    catch2Utils::init();

    GIVEN("Input")
    {
        perform(
                getStatsArguments(manager, true, false),
                {"", "-in", "../dataset/example.091.AA.strNOG.ENOG411BWBU.fasta"},
                manager);
    }

    GIVEN("Compareset")
    {
        perform(
                getStatsArguments(manager, true, true),
                {"", "-compareset", "../dataset/alignments_comparison.1"},
                manager);
    }
}