#include "../trimAlManagerIncludes.h"
#include "../../../../include/FormatHandling/FormatManager.h"
#include <random>
#include <sstream>

using namespace catch2Utils::trimAlManager;
using namespace trimAlManagerTestIncludes;

void getRangeCombinations(
        std::stringstream& mixture,
        std::vector<int>& valuesToBeAdded,
        int includedMax )
{
    static std::random_device rd;
    static std::mt19937 mt(rd());
    static std::uniform_real_distribution<float> dist(0.0, 1.0);

    float pos = dist(mt);

    if (pos < 0.4F)
    {
        if (mixture.tellp() != 0)
            mixture << ",";
        int a = (int) (dist(mt) * includedMax);
        mixture << a;
        valuesToBeAdded.push_back(a);
        valuesToBeAdded.push_back(a);
    }
    else if (pos < 0.8F)
    {
        if (mixture.tellp() != 0)
            mixture << ",";
        int a = (int) (dist(mt) * (includedMax - 1));
        int b = (int) (dist(mt) * (includedMax - a) + a);
        mixture << a << "-" << b;
        valuesToBeAdded.push_back(a);
        valuesToBeAdded.push_back(b);
    }
}

TEST_CASE("Select (cols|seqs) argument parse", "[manager][arguments][manual][selectX]") {
    // Tests initialization
    catch2Utils::init();

    trimAlManager manager;

    // Useful for selectCols and selectSeqs arguments.
    //  Provides a stringstream filled with ranges and points
    //  in the format X,X,Y-Y,X,Y-Y,X,X
    // The order of ranges and points is random, to provide better insight

    std::unique_ptr<Alignment> alig =
            std::unique_ptr<Alignment>(
                    manager.getFormatManager()
                            .loadAlignment(alignmentPath));



    for (auto & entry : getSelectX(manager, alig.get()))
    {
        GIVEN(entry.name << " argument")
        {
            AND_GIVEN("Input Alignment")
            {
                AND_GIVEN(entry.name << " in a correct range")
                {
                    WHEN(entry.name << " points to a single position")
                    {
                        for (int i = 0, ix = entry.maxValue; i < ix; i++)
                        {
                            WHEN(entry.name << " is " << std::to_string(i))
                            {
                                test_arguments(
                                        {"",
                                         "-in", alignmentPath,
                                         entry.argument, "{", std::to_string(i).c_str(), "}"},
                                        manager,
                                        true, trimAlManager::argumentReport::Recognized,
                                        true, false,
                                        [&entry, i](){
                                            THEN("Values are set correctly")
                                            {
                                                int * del = *entry.pointerToPointerVar;
                                                REQUIRE(del[0] == 2);
                                                CHECK(del[1] == i);
                                                CHECK(del[2] == i);
                                            }
                                        });
                            }
                        }
                    }

                    WHEN(entry.name << " points to several positions")
                    {
                        for (int i = 0, ix = entry.maxValue; i < ix; i++)
                        {
                            // GIVEN statement needed to reduce the time used in the loops
                            //      as catch2 tends to go over to
                            //      check if the test has been performed
                            // With the GIVEN statement, we skip the inner loop.
                            GIVEN("Start " << i)
                            {
                                for (int x = i; x < ix; x++)
                                {
                                    WHEN(entry.name << " is " << std::to_string(i) << "," << std::to_string(x))
                                    {
                                        test_arguments(
                                                {"",
                                                 "-in", alignmentPath,
                                                 entry.argument, "{",
                                                 (std::to_string(i) + "," + std::to_string(x)).c_str(),
                                                 "}"},
                                                manager,
                                                true, trimAlManager::argumentReport::Recognized,
                                                true, false,
                                                [&entry, i, x](){
                                                    THEN("Values are set correctly")
                                                    {
                                                        int * del = *entry.pointerToPointerVar;
                                                        REQUIRE(del[0] == 4);
                                                        CHECK(del[1] == i);
                                                        CHECK(del[2] == i);
                                                        CHECK(del[3] == x);
                                                        CHECK(del[4] == x);
                                                    }
                                                });
                                    }
                                }
                            }
                        }
                    }

                    WHEN(entry.name << " points to a range")
                    {
                        for (int i = 0, ix = entry.maxValue; i < ix; i++)
                        {
                            // GIVEN statement needed to reduce the time used in the loops
                            //      as catch2 tends to go over to
                            //      check if the test has been performed
                            // With the GIVEN statement, we skip the inner loop.
                            GIVEN("Start " << i) {
                                for (int x = i; x < ix; x++) {
                                    WHEN(entry.name << " is " << std::to_string(i) << "-" << std::to_string(x)) {
                                        test_arguments(
                                                {"",
                                                 "-in", alignmentPath,
                                                 entry.argument, "{",
                                                 (std::to_string(i) + "-" + std::to_string(x)).c_str(),
                                                 "}"},
                                                manager,
                                                true, trimAlManager::argumentReport::Recognized,
                                                true, false,
                                                [&entry, i, x](){
                                                    THEN("Values are set correctly")
                                                    {
                                                        int * del = *entry.pointerToPointerVar;
                                                        REQUIRE(del[0] == 2);
                                                        CHECK(del[1] == i);
                                                        CHECK(del[2] == x);
                                                    }
                                                });
                                    }
                                }
                            }
                        }
                    }

                    WHEN(entry.name << " points to a mixture of ranges and positions")
                    {
                        for (int i = 0; i < 25; i++)
                        {
                            GIVEN("Replicate " << i)
                            {
                                std::stringstream ss;
                                std::vector<int> positions;
                                for (int x = 0; x < 100; x++)
                                    getRangeCombinations(ss, positions, entry.maxValue - 1);

                                std::string str = ss.str();

                                CAPTURE(str);
                                test_arguments(
                                        {"",
                                         "-in", alignmentPath,
                                         entry.argument, "{", str.c_str() ,"}"},
                                        manager,
                                        true, trimAlManager::argumentReport::Recognized,
                                        true, false,
                                        [&entry, &positions](){
                                            THEN("Values are set correctly")
                                            {
                                                // Pointer to del(Sequences|Residues)
                                                int * del = *entry.pointerToPointerVar;
                                                REQUIRE(*del == positions.size());

                                                // Move the del(S|R) to the next position (first valid)
                                                del++;

                                                for (int i = 0; i < positions.size(); i++)
                                                {
                                                    CAPTURE(i);
                                                    CHECK(del[i] == positions[i]);
                                                }

                                            }
                                        });
                            }
                        }
                    }
                }

                AND_GIVEN(entry.name << " with an out of bounds")
                {
                    test_arguments(
                            {"",
                             "-in", alignmentPath,
                             entry.argument, "{", "20000000", "}"},
                            manager,
                            true, trimAlManager::argumentReport::Recognized,
                            true, true);
                }
            }
            AND_GIVEN("No Input")
            {
                test_arguments(
                        {"",entry.argument, "{", "0", "}"},
                        manager,
                        true, trimAlManager::argumentReport::Errored,
                        false, true);
            }
        }
    }
}

TEST_CASE("Cluster method argument parse", "[manager][arguments][manual][clusters]") {
    // Tests initialization
    catch2Utils::init();

    trimAlManager manager;

    // Useful for selectCols and selectSeqs arguments.
    //  Provides a stringstream filled with ranges and points
    //  in the format X,X,Y-Y,X,Y-Y,X,X
    // The order of ranges and points is random, to provide better insight

    std::unique_ptr<Alignment> alig =
            std::unique_ptr<Alignment>(
                    manager.getFormatManager()
                            .loadAlignment(alignmentPath));

    GIVEN("No Input")
    {
        GIVEN("No number of clusters")
            test_arguments(
                    {"", "-clusters"},
                    manager,
                    true, trimAlManager::argumentReport::Errored,
                    true, true);

        GIVEN("-1 clusters")
            test_arguments(
                    {"", "-clusters", std::to_string(-1).c_str()},
                    manager,
                    true, trimAlManager::argumentReport::Errored,
                    true, true);

        for (int i = 0; i < 2; i++)
            GIVEN(std::to_string(i) << " clusters")
                test_arguments(
                        {"", "-clusters", std::to_string(i).c_str()},
                        manager,
                        true, trimAlManager::argumentReport::Errored,
                        false, true);
    }

    GIVEN("Alignment as Input")
    {
        GIVEN("No number of clusters")
            test_arguments(
                    {"",
                     "-in", alignmentPath,
                     "-clusters"},
                    manager,
                    true, trimAlManager::argumentReport::Errored,
                    true, true);

        GIVEN("-1/" << alig->numberOfSequences << " clusters")
            test_arguments(
                    {"",
                     "-in", alignmentPath,
                     "-clusters", "-1"},
                    manager,
                    true, trimAlManager::argumentReport::Errored,
                    true, true);

        GIVEN("0/" << alig->numberOfSequences << " clusters")
            test_arguments(
                    {"",
                     "-in", alignmentPath,
                     "-clusters", "0"},
                    manager,
                    true, trimAlManager::argumentReport::Errored,
                    true, true);

        for (int i = 1; i <= alig->numberOfSequences; i++)
            GIVEN(std::to_string(i) << "/" << alig->numberOfSequences << " clusters")
                test_arguments(
                        {"",
                         "-in", alignmentPath,
                         "-clusters", std::to_string(i).c_str()},
                        manager,
                        true, trimAlManager::argumentReport::Recognized,
                        true, false,
                        [&manager, i](){
                          THEN("Clusters is set correctly in the manager")
                          {
                              CHECK(manager.clusters == i);
                          }
                        });

        GIVEN(alig->numberOfSequences + 1 << "/" << alig->numberOfSequences << " clusters")
            test_arguments(
                    {"", "-in", alignmentPath,
                     "-clusters", std::to_string(alig->numberOfSequences + 1).c_str()},
                    manager,
                    true, trimAlManager::argumentReport::Recognized,
                    true, true);
    }
}

TEST_CASE("MaxIdentity", "[manager][arguments][manual][maxidentity]")
{
    catch2Utils::init();

    trimAlManager manager;
    GIVEN("Max Identity argument")
    {
        GIVEN("No value for max identity")
        {
            test_arguments(
                    {"", "-in", alignmentPath,
                     "-maxidentity"},
                    manager,
                    true, trimAlManager::argumentReport::Errored,
                    true, true);
        }

        GIVEN("Value for max identity")
        {
            for (auto & value : std::vector<std::tuple<const char *, trimAlManager::argumentReport, bool>>
            {
                    {"-1",      trimAlManager::argumentReport::Errored,     true},
                    {"-0.1",    trimAlManager::argumentReport::Errored,     true},
                    {"0",       trimAlManager::argumentReport::Recognized,  false},
                    {"0.5",     trimAlManager::argumentReport::Recognized,  false},
                    {"1",       trimAlManager::argumentReport::Recognized,  false},
            })
            {
                GIVEN(std::get<0>(value) << " as value for max identity")
                {
                    test_arguments(
                            {"", "-in", alignmentPath,
                             "-maxidentity", std::get<0>(value)},
                            manager,
                            true, std::get<1>(value),
                            true, std::get<2>(value),
                                    [&manager, &value](){
                                if (std::get<2>(value) || std::get<1>(value) != trimAlManager::argumentReport::Recognized)
                                    return;
                                THEN("MaxIdentity is set correctly")
                                {
                                    CHECK(manager.maxIdentity == std::stof(std::get<0>(value)));
                                }
                            });
                }
            }
        }
    }
}

void iterateCombinations(
        const std::vector<argumentReference>& booleanModifiers,
        std::vector<const char*> args,
        trimAlManager& manager,
        const trimAlManager::argumentReport expectedReport,
        const bool expectedError,
        std::vector<argumentReference> booleanModifiersApplied = {},
        const uint32_t i = 0,
        const bool performRequested = false)
{
    if (performRequested)
    {
        test_arguments(
                args, manager,
                        true, expectedReport,
                        true, expectedError,
            // Check the value has been set correctly
            [&booleanModifiersApplied](){
                THEN("The variable in the manager is set correctly")
                {
                    for(auto & booleanModifier: booleanModifiersApplied)
                    {
                        THEN(booleanModifier.argument + 1 << " variable value is true")
                        {
                            CHECK(booleanModifier.reference);
                        }
                    }
                }
            });
    }

    if (i < booleanModifiers.size())
    {
        iterateCombinations(
                booleanModifiers, args, manager,
                expectedReport, expectedError, booleanModifiersApplied,
                i + 1, false);


        GIVEN(booleanModifiers[i].argument)
        {
            args.push_back(booleanModifiers[i].argument);
            booleanModifiersApplied.push_back(booleanModifiers[i]);
            iterateCombinations(
                    booleanModifiers, args, manager,
                    expectedReport, expectedError, booleanModifiersApplied,
                    i + 1, true);
        }

    }
}

// Valgrinded

TEST_CASE("Boolean modifiers", "[manager][arguments][manual][terminalonly][automated]")
{
    catch2Utils::init();

    trimAlManager manager;



    std::vector<const char *> args {"", "-in", alignmentPath};

    auto automatedMethods = getAutomatedMethods(manager);

    for(auto& automatedMethod : automatedMethods)
    {
        GIVEN(&automatedMethod.argument[1] << " and input")
        {
            args.push_back(automatedMethod.argument);
            iterateCombinations(
                    getBooleanModifiers(manager), args, manager,
                    trimAlManager::argumentReport::Recognized, false);
        }
    }

}