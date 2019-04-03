#include "../trimAlManagerIncludes.h"
#include <stdlib.h>

using namespace catch2Utils::trimAlManager;
using namespace trimAlManagerTestIncludes;

// Valgrinded

TEST_CASE("Block", "[manager][arguments][block]") {

    catch2Utils::init();

    trimAlManager manager;

    std::vector<const char *> args {""};

    GIVEN("Block argument")
    {
        args.push_back("-block");
        test_arguments(args, manager,
                       true, trimAlManager::argumentReport::Errored,
                       true, true);
        AND_GIVEN("Block argument value")
        {
            std::unique_ptr<Alignment> alig {
                    manager.getFormatManager().loadAlignment(alignmentPath)
            };

            std::vector<std::tuple<int, trimAlManager::argumentReport ,bool>> values {
                    {0, trimAlManager::argumentReport::Errored, true},
                    {1, trimAlManager::argumentReport::Recognized, false},
                    {alig->numberOfResidues * 2, trimAlManager::argumentReport::Recognized, true},
                    {alig->numberOfResidues / 4, trimAlManager::argumentReport::Recognized, false}
            };

            for(auto & value : values)
            {
                WHEN("Given value is " << std::get<0>(value))
                {
                    std::unique_ptr<char[]> arr{new char[10]};
                    std::sprintf(arr.get(),"%d",std::get<0>(value));
                    args.push_back(arr.get());

                    test_arguments(args, manager,
                                   true, trimAlManager::argumentReport::Errored,
                                   true, true);
                    AND_GIVEN("Input")
                    {
                        args.push_back("-in");
                        args.push_back(alignmentPath);

                        if (std::get<2>(value) &&
                                std::get<1>(value) == trimAlManager::argumentReport::Recognized) {
                            test_arguments(
                                    args, manager,
                                    true, std::get<1>(value),
                                    true, std::get<2>(value),
                                    [&manager, &value]() {
                                        THEN("Blocksize is stored correctly") {
                                            CHECK(manager.blockSize == std::get<0>(value));
                                        }
                                    });
                        } else {
                            test_arguments(
                                    args, manager,
                                    true, std::get<1>(value),
                                    true, std::get<2>(value));
                        }
                    }
                }
            }

            WHEN("Given value is ATCG")
            {
                args.push_back("ATCG");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Errored,
                               true, true);
                AND_GIVEN("Input")
                {
                    args.push_back("-in");
                    args.push_back(alignmentPath);
                    test_arguments(
                            args, manager,
                            true, trimAlManager::argumentReport::Errored,
                            true, true);
                }
            }
        }
    }

}
