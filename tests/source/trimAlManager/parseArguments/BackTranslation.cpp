#include "../trimAlManagerIncludes.h"

using namespace catch2Utils::trimAlManager;
//using namespace trimAlManagerTestIncludes;
using namespace trimAlManagerTestIncludes::backtranslation;


// Valgrinded

TEST_CASE("Backtranslation argument parse", "[manager][arguments][backtranslation]") {

    // Tests initialization
    trimAlManager manager;
    catch2Utils::init();

    GIVEN("Backtranslation argument") {
        WHEN("Correct codon file as backtranslation file")
        {
            std::vector<const char *> args {
                "",
                "-in", backtranslationAlignmentPath,
                "-backtrans", backtranslationCodonPath};

            test_arguments(args, manager,
                        true, trimAlManager::argumentReport::Recognized,
                        true, false,
                           [&](){
                               THEN("The variables in the manager are set correctly")
                               {
                                   AND_THEN("Infile is set correctly")
                                   {
                                       REQUIRE(manager.infile != nullptr);
                                       CHECK(strcmp(manager.backtransFile, backtranslationCodonPath) == 0);
                                   }
                                   AND_THEN("Backtransfile is set correctly")
                                    {
                                        REQUIRE(manager.backtransFile != nullptr);
                                        CHECK(strcmp(manager.infile, backtranslationAlignmentPath) == 0);
                                    }
                                   AND_THEN("IgnoreStopCodon is set corretly")
                                       CHECK_FALSE(manager.ignoreStopCodon);

                                   AND_THEN("SplitByStopCodon is set corretly")
                                       CHECK_FALSE(manager.splitByStopCodon);
                               }
                           });

            GIVEN("IgnoreStopCodon argument")
            {
                args.push_back("-ignorestopcodon");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Recognized,
                               true, false,
                               [&](){
                                   THEN("The variables in the manager are set correctly")
                                   {
                                       AND_THEN("Infile is set correctly")
                                       {
                                           REQUIRE(manager.infile != nullptr);
                                           CHECK(strcmp(manager.infile, backtranslationAlignmentPath) == 0);
                                       }
                                       AND_THEN("Backtransfile is set correctly")
                                       {
                                           REQUIRE(manager.backtransFile != nullptr);
                                           CHECK(strcmp(manager.backtransFile, backtranslationCodonPath) == 0);
                                       }

                                       AND_THEN("IgnoreStopCodon is set corretly")
                                        CHECK(manager.ignoreStopCodon);

                                       AND_THEN("SplitByStopCodon is set corretly")
                                           CHECK_FALSE(manager.splitByStopCodon);
                                   }
                               });
            }

            GIVEN("Splitbystopcodon argument")
            {
                args.push_back("-splitbystopcodon");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Recognized,
                               true, false,
                               [&](){
                                   THEN("The variables in the manager are set correctly")
                                   {
                                       AND_THEN("Infile is set correctly")
                                       {
                                           REQUIRE(manager.infile != nullptr);
                                           CHECK(strcmp(manager.infile, backtranslationAlignmentPath) == 0);
                                       }
                                       AND_THEN("Backtransfile is set correctly")
                                       {
                                           REQUIRE(manager.backtransFile != nullptr);
                                           CHECK(strcmp(manager.backtransFile, backtranslationCodonPath) == 0);
                                       }

                                       AND_THEN("SplitByStopCodon is set corretly")
                                           CHECK(manager.splitByStopCodon);

                                       AND_THEN("IgnoreStopCodon is set corretly")
                                           CHECK_FALSE(manager.ignoreStopCodon);
                                   }
                               });
            }

            GIVEN("Splitbystopcodon and ignorestopcodon argument")
            {
                args.push_back("-ignorestopcodon");
                args.push_back("-splitbystopcodon");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Recognized,
                               true, true);
            }

        }
        WHEN("Alignment file as backtranslation file")
        {
            std::vector<const char *> args{
                    "",
                    "-in", backtranslationAlignmentPath,
                    "-backtrans", backtranslationAlignmentPath};

            test_arguments(args, manager,
                           true, trimAlManager::argumentReport::Recognized,
                           true, true);

            GIVEN("IgnoreStopCodon argument")
            {
                args.push_back("-ignorestopcodon");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Recognized,
                               true, true);
            }

            GIVEN("Splitbystopcodon argument")
            {
                args.push_back("-splitbystopcodon");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Recognized,
                               true, true);
            }

            GIVEN("Splitbystopcodon and ignorestopcodon argument")
            {
                args.push_back("-ignorestopcodon");
                args.push_back("-splitbystopcodon");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Recognized,
                               true, true);
            }
        }
        WHEN("Folder as backtranslation file")
        {
            std::vector<const char *> args{
                    "",
                    "-in", backtranslationAlignmentPath,
                    "-backtrans", "../dataset/"};

            test_arguments(args, manager,
                           true, trimAlManager::argumentReport::Errored,
                           true, true);

            GIVEN("IgnoreStopCodon argument")
            {
                args.push_back("-ignorestopcodon");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Errored,
                               true, true);
            }

            GIVEN("Splitbystopcodon argument")
            {
                args.push_back("-splitbystopcodon");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Errored,
                               true, true);
            }

            GIVEN("Splitbystopcodon and ignorestopcodon argument")
            {
                args.push_back("-ignorestopcodon");
                args.push_back("-splitbystopcodon");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Errored,
                               true, true);
            }
        }
        WHEN("Backtranslation file is not provided")
        {
            std::vector<const char *> args{
                    "",
                    "-in", backtranslationAlignmentPath,
                    "-backtrans"};

            test_arguments(args, manager,
                           true, trimAlManager::argumentReport::Errored,
                           true, true);

            GIVEN("IgnoreStopCodon argument")
            {
                args.push_back("-ignorestopcodon");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Errored,
                               true, true);
            }

            GIVEN("Splitbystopcodon argument")
            {
                args.push_back("-splitbystopcodon");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Errored,
                               true, true);
            }

            GIVEN("Splitbystopcodon and ignorestopcodon argument")
            {
                args.push_back("-ignorestopcodon");
                args.push_back("-splitbystopcodon");
                test_arguments(args, manager,
                               true, trimAlManager::argumentReport::Errored,
                               true, true);
            }
        }
    }

    GIVEN("No Backtranslation argument")
    {
        std::vector<const char *> args{
                "",
                "-in", backtranslationAlignmentPath};

        GIVEN("IgnoreStopCodon argument")
        {
            args.push_back("-ignorestopcodon");
            test_arguments(args, manager,
                           true, trimAlManager::argumentReport::Recognized,
                           true, true);
        }

        GIVEN("Splitbystopcodon argument")
        {
            args.push_back("-splitbystopcodon");
            test_arguments(args, manager,
                           true, trimAlManager::argumentReport::Recognized,
                           true, true);
        }

        GIVEN("Splitbystopcodon and ignorestopcodon argument")
        {
            args.push_back("-ignorestopcodon");
            args.push_back("-splitbystopcodon");
            test_arguments(args, manager,
                           true, trimAlManager::argumentReport::Recognized,
                           true, true);
        }
    }

}