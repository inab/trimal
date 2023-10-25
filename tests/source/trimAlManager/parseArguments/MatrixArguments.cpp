#include "../trimAlManagerIncludes.h"
#include "../../../../include/FormatHandling/FormatManager.h"
#include <random>
#include <sstream>
#include "../../../testingUtils/hasher.h"
#include "trimalManager.h"

using namespace catch2Utils::trimAlManager;
using namespace trimAlManagerTestIncludes;

trimAlManagerPrivateExposerCreator("similarity_matrix") {

    static void create_or_use_similarity_matrix(
            trimAlManager &manager,
            const std::function<void(bool)> &callback = {})
    {
        if (callback) callback(manager.create_or_use_similarity_matrix());
        else manager.create_or_use_similarity_matrix();
    };

    static void getSimilarityMatrix(
            trimAlManager &manager,
            const std::function<void(statistics::similarityMatrix*)> &callback = {})
    {
        if (callback) callback(manager.similMatrix);
    }
};

// Valgrinded

TEST_CASE("Matrix", "[manager][arguments][matrix]") {
    // Tests initialization
    catch2Utils::init();

    trimAlManager manager;

    std::vector<const char *> args{""};

    GIVEN("Matrix argument")
    {
        args.push_back("-matrix");
        test_arguments(
                args,
                manager,
                true, trimAlManager::argumentReport::Errored,
                true, true);

        AND_GIVEN("Correct value input")
        {
            args.push_back("../dataset/matrix.BLOSUM62");
            test_arguments(
                    args,
                    manager,
                    true, trimAlManager::argumentReport::Errored,
                    true, true);

            GIVEN("Input file")
            {
                args.push_back("-in");
                args.push_back(alignmentPath);
                test_arguments(
                        args,
                        manager,
                        true, trimAlManager::argumentReport::Recognized,
                        true, true);

                GIVEN("Similarity trimming method")
                {
                    args.push_back("-strict");
                    test_arguments(
                            args,
                            manager,
                            true, trimAlManager::argumentReport::Recognized,
                            true, false,
                            [&manager](){
                                THEN("Variable matrixFile is set correctly")
                                {
                                    REQUIRE(manager.matrixFile != nullptr);
                                    CHECK(strcmp(manager.matrixFile, "../dataset/matrix.BLOSUM62") == 0);
                                }

                                trimAlManagerPrivateExposerCaller("similarity_matrix")
                                ::create_or_use_similarity_matrix(manager, [](bool result) -> void {
                                    THEN("When calling create_or_use_similarity_matrix") {
                                        CHECK(result);
                                    }
                                });

                                trimAlManagerPrivateExposerCaller("similarity_matrix")
                                ::getSimilarityMatrix(manager, [](statistics::similarityMatrix *matrix) {
                                    THEN("Matrix is not null") {
                                        CHECK(matrix != nullptr);
                                    }
                                });
                            });

                }
                GIVEN("Non-similarity trimming method")
                {
                    args.push_back("-gappyout");
                    test_arguments(
                            args,
                            manager,
                            true, trimAlManager::argumentReport::Recognized,
                            true, true,
                            [&manager](){
                                THEN("Variable matrixFile is set correctly")
                                {
                                    REQUIRE(manager.matrixFile != nullptr);
                                    CHECK(strcmp(manager.matrixFile, "../dataset/matrix.BLOSUM62") == 0);
                                }

                                trimAlManagerPrivateExposerCaller("similarity_matrix")
                                ::create_or_use_similarity_matrix(manager, [](bool result) -> void {
                                    THEN("When calling create_or_use_similarity_matrix") {
                                        CHECK(result);
                                    }
                                });

                                trimAlManagerPrivateExposerCaller("similarity_matrix")
                                ::getSimilarityMatrix(manager, [](statistics::similarityMatrix *matrix) {
                                    THEN("Matrix is null") {
                                        CHECK(matrix == nullptr);
                                    }
                                });
                            });
                }
            }
        }
    }

    GIVEN("Alternative matrix argument")
    {
        args.push_back("--degenerated_nt_identity");
        test_arguments(
                args,
                manager,
                true, trimAlManager::argumentReport::Errored,
                true, true);

        AND_GIVEN("Correct value input")
        {
            GIVEN("Input file")
            {
                args.push_back("-in");
                args.push_back(alignmentPath);
                test_arguments(
                        args,
                        manager,
                        true, trimAlManager::argumentReport::Recognized,
                        true, true);

                GIVEN("Similarity trimming method")
                {
                    args.push_back("-strict");
                    test_arguments(
                            args,
                            manager,
                            true, trimAlManager::argumentReport::Recognized,
                            true, false,
                            [&manager](){
                                THEN("Variable matrixFile is null")
                                {
                                    REQUIRE(manager.matrixFile == nullptr);
                                }

                                trimAlManagerPrivateExposerCaller("similarity_matrix")
                                ::create_or_use_similarity_matrix(manager, [](bool result) -> void {
                                    THEN("Calling create_or_use_similarity_matrix") {
                                        CHECK(result);
                                    }
                                });

                                trimAlManagerPrivateExposerCaller("similarity_matrix")
                                ::getSimilarityMatrix(manager, [](statistics::similarityMatrix *matrix) {
                                    THEN("Matrix is not null") {
                                        CHECK(matrix != nullptr);
                                    }
                                });
                            });

                }
                GIVEN("Non-similarity trimming method")
                {
                    args.push_back("-gappyout");
                    test_arguments(
                            args,
                            manager,
                            true, trimAlManager::argumentReport::Recognized,
                            true, true,
                            [&manager](){
                                THEN("Variable matrixFile is null")
                                {
                                    REQUIRE(manager.matrixFile == nullptr);
                                }

                                trimAlManagerPrivateExposerCaller("similarity_matrix")
                                ::create_or_use_similarity_matrix(manager, [](bool result) -> void {
                                    THEN("Calling create_or_use_similarity_matrix") {
                                        CHECK(result);
                                    }
                                });

                                trimAlManagerPrivateExposerCaller("similarity_matrix")
                                ::getSimilarityMatrix(manager, [](statistics::similarityMatrix *matrix) {
                                    THEN("Matrix is null") {
                                        CHECK(matrix == nullptr);
                                    }
                                });
                            });
                }
            }
        }
    }








}