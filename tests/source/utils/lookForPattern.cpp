#include "../../catch.hpp"
#include "utils.h"

TEST_CASE("lookForPattern", "[utils]") {

    WHEN("looking for pattern")
    {
        std::string data, pattern;
        GIVEN("AAATTT - A")
        {
            data = "AAATTT";
            pattern = "A";
            CAPTURE(data, pattern);
            REQUIRE(utils::lookForPattern(data, pattern, 0.5F));

            REQUIRE(data == "AAATTT");
            REQUIRE(pattern == "A");
        }

        GIVEN("AAATTT - a")
        {
            data = "AAATTT";
            pattern = "a";
            CAPTURE(data, pattern);
            REQUIRE(utils::lookForPattern(data, pattern, 0.5F));

            REQUIRE(data == "AAATTT");
            REQUIRE(pattern == "a");}

        GIVEN("AATTTT - A")
        {
            data = "AATTTT";
            pattern = "A";
            CAPTURE(data, pattern);
            REQUIRE_FALSE(utils::lookForPattern(data, pattern, 0.5F));

            REQUIRE(data == "AATTTT");
            REQUIRE(pattern == "A");
        }

        GIVEN("AATTTT - a")
        {
            data = "AATTTT";
            pattern = "a";
            CAPTURE(data, pattern);
            REQUIRE_FALSE(utils::lookForPattern(data, pattern, 0.5F));

            REQUIRE(data == "AATTTT");
            REQUIRE(pattern == "a");
        }

        GIVEN("ATATCGCG - AT")
        {
            data = "ATATCGCG";
            pattern = "AT";
            CAPTURE(data, pattern);
            REQUIRE(utils::lookForPattern(data, pattern, 0.5F));

            REQUIRE(data == "ATATCGCG");
            REQUIRE(pattern == "AT");
        }

        GIVEN("ATATCGCG - at")
        {
            data = "ATATCGCG";
            pattern = "at";
            CAPTURE(data, pattern);
            REQUIRE(utils::lookForPattern(data, pattern, 0.5F));

            REQUIRE(data == "ATATCGCG");
            REQUIRE(pattern == "at");
        }

        GIVEN("aTaTCGCG - At")
        {
            data = "ATATCGCG";
            pattern = "At";
            CAPTURE(data, pattern);
            REQUIRE(utils::lookForPattern(data, pattern, 0.5F));

            REQUIRE(data == "ATATCGCG");
            REQUIRE(pattern == "At");
        }
    }

    WHEN("Determining Colors")
    {
        WARN("This test needs profound design by the original author or reference to compare against.");
    }
}

TEST_CASE("[B] lookForPattern", "[!benchmark][utils]") {

    WHEN("looking for pattern")
    {
        std::string data, pattern;

        for (int patternPow = 1; patternPow < 6; patternPow++)
        {
            auto patternSize = (ulong) pow(2, patternPow);

            pattern.reserve(patternSize);

            for (unsigned long d = pattern.size(); d < patternSize; d++)
                pattern += (char)(rand() % ('Z' - 'A') + 'A');

            for (int dataPow = 1; dataPow < 14; dataPow++)
            {
                auto dataSize = (ulong) pow(2, dataPow);

                data.reserve(dataSize);

                for (unsigned long d = data.size(); d < dataSize; d++)
                    data += (char)(rand() % ('Z' - 'A') + 'A');

                BENCHMARK("Pattern: " + std::to_string(patternSize) + "     Data: " + std::to_string(dataSize))
                    utils::lookForPattern(data, pattern, 0.5F);
            }
            std::cout << std::endl;
         }
    }
}