#include "../trimAlManagerIncludes.h"

using namespace catch2Utils::trimAlManager;
using namespace trimAlManagerTestIncludes;

// Valgrinded

TEST_CASE("Input argument parse", "[manager][arguments][in][input]") {

    catch2Utils::init();

    trimAlManager manager;

    WHEN("Input is correct")
    {
        test_arguments(
                {"", "-in", alignmentPath},
                manager,
                true, trimAlManager::argumentReport::Recognized,
                true, false,
                [&manager](){
                    THEN("The variables in the manager are set correctly")
                    {
                        REQUIRE(manager.infile != nullptr);
                        CHECK(strcmp(manager.infile, alignmentPath) == 0);
                    }
                });

        AND_WHEN("Input argument is passed two times")
        {
            test_arguments(
                    {"",
                     "-in", alignmentPath,
                     "-in", alignmentPath},
                    manager,
                    true, trimAlManager::argumentReport::Errored,
                    true, true);
        }
    }

    WHEN("Input is incorrect")
    {
        test_arguments(
                {"", "-in", "../dataset/example..010.AA.fasta"},
                manager,
                true, trimAlManager::argumentReport::Errored,
                true, true);

        AND_WHEN("Input argument is passed two times")
        {
            test_arguments(
                    {"",
                     "-in", "../dataset/example..010.AA.fasta",
                     "-in", "../dataset/example..010.AA.fasta"},
                    manager,
                    true, trimAlManager::argumentReport::Errored,
                    true, true);
        }
    }


    WHEN("Input is a folder")
    {
        test_arguments(
                {"", "-in", "../dataset/"},
                manager,
                true, trimAlManager::argumentReport::Errored,
                true, true);
        AND_WHEN("Input argument is passed two times")
        {
            test_arguments(
                    {"",
                     "-in", "../dataset/",
                     "-in", "../dataset/"},
                    manager,
                    true, trimAlManager::argumentReport::Errored,
                    true, true);
        }
    }

}
