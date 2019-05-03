#include "../trimAlManagerIncludes.h"

using namespace catch2Utils::trimAlManager;

/**
 * Method to check test-case adding input information
 * @param args arguments to pass to test_arguments
 * @param manager manager reference to current manager
 * @param expectedReport expected report of trimAlManager::parseArguments
 * @param callback Callback to be applied by test_arguments when input applied is correct
 */
void addInput(
        std::vector<const char *>& args,
        trimAlManager& manager,
        const trimAlManager::argumentReport expectedReport,
        const testArgumentsCallback& callback = {})
{
    GIVEN("Input argument")
    {
        args.push_back("-in");
        WHEN("When input is an alignment")
        {
            args.push_back("../dataset/example.004.AA.fasta");
            test_arguments(args, manager,
                    true, expectedReport,
                    true, true, callback);
        }
//        WHEN("When input is a folder")
//        {
//            args.push_back("../dataset/");
//            test_arguments(args, manager,
//                    true, trimAlManager::argumentReport::Errored,
//                    true, true);
//        }
//        WHEN("When input is not provided")
//        {
//            test_arguments(args, manager,
//                    true, trimAlManager::argumentReport::Errored,
//                    true, true);
//        }
    }
}

// Valgrinded FAIL

TEST_CASE("Compareset argument parse", "[manager][arguments][compareset]") {

    catch2Utils::init();

    trimAlManager manager;

    std::vector<const char*> args{"", "-compareset"};

    WHEN("Compareset is an alignment") {

        args.push_back("../dataset/example.004.AA.fasta");

        GIVEN("No forceselect") {
            test_arguments(args, manager,
                    true, trimAlManager::argumentReport::Errored,
                    true, true);
            addInput(args, manager, trimAlManager::argumentReport::Errored);
        }

        GIVEN("Forceselect")
        {
            args.push_back("-forceselect");
            GIVEN("Forceselect as alignment") {
                args.push_back("../dataset/example.004.AA.fasta");
                test_arguments(args, manager,
                        true, trimAlManager::argumentReport::Errored,
                        true, true);
                addInput(args, manager, trimAlManager::argumentReport::Errored);
            }

            GIVEN("Forceselect as folder") {
                args.push_back("../dataset/");
                test_arguments(args, manager,
                        true, trimAlManager::argumentReport::Errored,
                        true, true);
                addInput(args, manager, trimAlManager::argumentReport::Errored);
            }

            GIVEN("Forceselect as empty") {
                test_arguments(args, manager,
                        true, trimAlManager::argumentReport::Errored,
                        true, true);
                addInput(args, manager, trimAlManager::argumentReport::Errored);
            }
        }

    }

    WHEN("Compareset is correct") {
        args.push_back("../dataset/alignments_comparison.1");


        test_arguments(args, manager,
                true, trimAlManager::argumentReport::Recognized,
                true, false,
                [&](){
                    THEN("The variables in the manager are set correctly")
                    {
                        REQUIRE(manager.compareset != nullptr);
                        CHECK(strcmp(manager.compareset, "../dataset/alignments_comparison.1") == 0);
                        CHECK(manager.infile == nullptr);
                    }
                });
        addInput(args, manager, trimAlManager::argumentReport::Recognized);

        GIVEN("Forceselect argument")
        {
            args.push_back("-forceselect");

            GIVEN("Forceselect as alignment") {
                WHEN("ForceSelect is an alignment of the set") {
                    args.push_back("../dataset/example.001.AA.phy");
                    test_arguments(args, manager,
                                   true, trimAlManager::argumentReport::Recognized,
                                   true, false,
                                   [&](){
                                       THEN("The variables in the manager are set correctly")
                                       {
                                           REQUIRE(manager.compareset != nullptr);
                                           REQUIRE(manager.forceFile != nullptr);

                                           CHECK(strcmp(manager.compareset, "../dataset/alignments_comparison.1") == 0);
                                           CHECK(strcmp(manager.forceFile, "../dataset/example.001.AA.phy") == 0);
                                           CHECK(manager.infile == nullptr);
                                       }
                                   });
                    addInput(args, manager, trimAlManager::argumentReport::Recognized);
                }

                WHEN("ForceSelect is NOT an alignment of the set") {
                    INFO("\u202F\n"
                         "------------------------------------------------------------\n"
                         "Explanation:\n"
                         "  The suite should not allow to input an \n"
                         "    external alignment (not present in the compareset file).\n"
                         "  This is due to uncertainity produced to the final user:\n"
                         "  Ex: If the alignment is not on the set,\n "
                         "    is it used to calculate CT ?\n"
                         "------------------------------------------------------------\n"
                         "\u202F")
                    // "\u202F" => NARROW NO-BREAK SPACE -> Multiple \n are ignored on CATCH.
                    // To trick the system and allow to multiple line breaks,
                    //  we use an 'invisible' character to allow separate blocks with
                    //  several line breaks
                    args.push_back("../dataset/example.004.AA.fasta");
                    test_arguments(args, manager,
                                   true, trimAlManager::argumentReport::Errored,
                                   true, true);
                    addInput(args, manager, trimAlManager::argumentReport::Errored);
                }
            }

            GIVEN("Forceselect as folder") {
                args.push_back("../dataset/");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Errored,
                               true, true);
                addInput(args, manager, trimAlManager::argumentReport::Errored);
            }

            GIVEN("Forceselect as empty") {
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Errored,
                               true, true);
                addInput(args, manager, trimAlManager::argumentReport::Errored);
            }
        }

    }

    WHEN("Compareset is a folder") {
        args.push_back("../dataset");
        test_arguments(args, manager,
                false, trimAlManager::argumentReport::Errored,
                true, true);

        GIVEN("Forceselect")
        {
            args.push_back("-forceselect");
            GIVEN("Forceselect as alignment") {
                args.push_back("../dataset/example.004.AA.fasta");
                test_arguments(args, manager,
                        true, trimAlManager::argumentReport::Errored,
                        true, true);
                addInput(args, manager, trimAlManager::argumentReport::Errored);
            }

            GIVEN("Forceselect as folder") {
                args.push_back("../dataset/");
                test_arguments(args, manager,
                        true, trimAlManager::argumentReport::Errored,
                        true, true);
                addInput(args, manager, trimAlManager::argumentReport::Errored);
            }
            GIVEN("Forceselect as empty") {
                test_arguments(args, manager,
                        true, trimAlManager::argumentReport::Errored,
                        true, true);
                addInput(args, manager, trimAlManager::argumentReport::Errored);
            }
        }
    }
}
