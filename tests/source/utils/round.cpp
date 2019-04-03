#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-flp30-c"

#include "../../catch.hpp"
#include "utils.h"

TEST_CASE("round", "[utils]") {
    SECTION("roundToInf") {
        GIVEN("doubles") {
            double x;
            WHEN("In range [0,1) == 0")for (x = 0; x < 1; x += 0.1)
                    CHECK(utils::roundToInf(x) == 0);

            x = 1;
            CHECK(utils::roundToInf(x) == 1);

            WHEN("In range (1,2] == 1")for (; x < 2; x += 0.1)
                    CHECK(utils::roundToInf(x) == 1);
        }

        GIVEN("floats") {
            float x;
            WHEN("In range [0,1) == 0")for (x = 0; x < 1; x += 0.1)
                    CHECK(utils::roundToInf(x) == 0);

            x = 1;
            CHECK(utils::roundToInf(x) == 1);

            WHEN("In range (1,2] == 1")for (; x < 2; x += 0.1)
                    CHECK(utils::roundToInf(x) == 1);
        }

        GIVEN("ints [0-1000)")THEN("Input == Output")for (int i = 0; i < 1000; i++)
                    CHECK(utils::roundToInf(i) == i);

        GIVEN("uints [0-1000)")THEN("Input == Output")for (unsigned int i = 0; i < 1000; i++)
                    CHECK(utils::roundToInf(i) == i);
    }

    SECTION("roundToSup") {
        GIVEN("doubles") {
            double x;
            WHEN("In range [0,1) == 1")for (x = 0; x <= 0.9; x += 0.1)
                    GIVEN(std::to_string(x)) CHECK(utils::roundToSup(x) == 1);

            x = 1;
            CHECK(utils::roundToInf(x) == 1);

            WHEN("In range (1,2] == 2")for (; x < 2; x += 0.1)
                    GIVEN(x) CHECK(utils::roundToSup(x) == 2);
        }

        GIVEN("floats") {
            float x;
            WHEN("In range [0,1) == 1")for (x = 0; x <= 0.9F; x += 0.1)
                    GIVEN(x) CHECK(utils::roundToSup(x) == 1);

            x = 1;
            CHECK(utils::roundToSup(x) == 2);

            WHEN("In range (1,2] == 2")for (; x < 2; x += 0.1)
                    GIVEN(x) CHECK(utils::roundToSup(x) == 2);
        }

        GIVEN("ints [0-10)")THEN("Input == Output + 1")for (int i = 0; i < 10; i++)
                    GIVEN(i) CHECK(utils::roundToSup(i) == i + 1);

        GIVEN("uints [0-10)")THEN("Input == Output + 1")for (unsigned int i = 0; i < 10; i++)
                    GIVEN(i) CHECK(utils::roundToSup(i) == i + 1);
    }

    SECTION("roundToInt") {
        GIVEN("doubles") {
            double x;
            WHEN("In range [0,0.5) == 1")for (x = 0; x <= 0.4; x += 0.1)
                    GIVEN(std::to_string(x))CHECK(utils::roundInt(x) == 0);
            WHEN("In range [0.5,1) == 1")for (x = 0.5; x <= 0.9; x += 0.1)
                    GIVEN(std::to_string(x))CHECK(utils::roundInt(x) == 1);
            x = 1;
            CHECK(utils::roundInt(x) == 1);
        }

        GIVEN("floats") {
            float x;
            WHEN("In range [0,0.5) == 1")for (x = 0; x <= 0.4; x += 0.1)
                    GIVEN(std::to_string(x))CHECK(utils::roundInt(x) == 0);
            WHEN("In range [0.5,1) == 1")for (x = 0.5; x <= 0.9; x += 0.1)
                    GIVEN(std::to_string(x))CHECK(utils::roundInt(x) == 1);
            x = 1;
            CHECK(utils::roundInt(x) == 1);
        }


        GIVEN("ints [0-10)")THEN("Input == Output + 1")for (int i = 0; i < 10; i++)
                    GIVEN(i) CHECK(utils::roundInt(i) == i);

        GIVEN("uints [0-10)")THEN("Input == Output + 1")for (unsigned int i = 0; i < 10; i++)
                    GIVEN(i) CHECK(utils::roundInt(i) == i);
    }
}

TEST_CASE("[B] round", "[!benchmark][utils]") {
    SECTION("roundToInf") {
        GIVEN("doubles")for (double x = 0; x < 1; x += 0.1)
                BENCHMARK("roundToInf double " + std::to_string(x))utils::roundToInf(x);

        GIVEN("floats")for (float x = 0; x < 1; x += 0.1)
                BENCHMARK("roundToInf float " + std::to_string(x))utils::roundToInf(x);
    }

    SECTION("roundToSup") {
        GIVEN("doubles")for (double x = 0; x < 1; x += 0.1)
                BENCHMARK("roundToSup double " + std::to_string(x))utils::roundToSup(x);

        GIVEN("floats")for (float x = 0; x < 1; x += 0.1)
                BENCHMARK("roundToSup float " + std::to_string(x))utils::roundToSup(x);
    }

    SECTION("roundToInt") {
        GIVEN("doubles")for (double x = 0; x < 1; x += 0.1)
                BENCHMARK("roundInt double " + std::to_string(x))utils::roundInt(x);

        GIVEN("floats")for (float x = 0; x < 1; x += 0.1)
                BENCHMARK("roundInt float " + std::to_string(x))utils::roundInt(x);
    }
}

#pragma clang diagnostic pop