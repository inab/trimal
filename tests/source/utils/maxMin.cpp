#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-flp30-c"

#include "../../catch.hpp"
#include "utils.h"

TEST_CASE("maxmin", "[utils]") {
    SECTION("Max") {
        GIVEN("ints") {
            CHECK(utils::max(0, 0) == 0);
            CHECK(utils::max(0, 1) == 1);
            CHECK(utils::max(1, 0) == 1);
            CHECK(utils::max(0, -1) == 0);
            CHECK(utils::max(-1, 0) == 0);
        }
        GIVEN("floats") {
            CHECK(utils::max(0.0F, 0.0F) == 0.0F);
            CHECK(utils::max(0.0F, 1.0F) == 1.0F);
            CHECK(utils::max(1.0F, 0.0F) == 1.0F);
            CHECK(utils::max(0.0F, -1.0F) == 0.0F);
            CHECK(utils::max(-1.0F, 0.0F) == 0.0F);
        }
        GIVEN("doubles") {
            CHECK(utils::max(0.0, 0.0) == 0.0);
            CHECK(utils::max(0.0, 1.0) == 1.0);
            CHECK(utils::max(1.0, 0.0) == 1.0);
            CHECK(utils::max(0.0, -1.0) == 0.0);
            CHECK(utils::max(-1.0, 0.0) == 0.0);
        }
    }

    SECTION("Min") {
        GIVEN("ints") {
            CHECK(utils::min(0, 0) == 0);
            CHECK(utils::min(0, 1) == 0);
            CHECK(utils::min(1, 0) == 0);
            CHECK(utils::min(0, -1) == -1);
            CHECK(utils::min(-1, 0) == -1);
        }
        GIVEN("floats") {
            CHECK(utils::min(0.0F, 0.0F) == 0.0F);
            CHECK(utils::min(0.0F, 1.0F) == 0.0F);
            CHECK(utils::min(1.0F, 0.0F) == 0.0F);
            CHECK(utils::min(0.0F, -1.0F) == -1.0F);
            CHECK(utils::min(-1.0F, 0.0F) == -1.0F);
        }
        GIVEN("doubles") {
            CHECK(utils::min(0.0, 0.0) == 0.0);
            CHECK(utils::min(0.0, 1.0) == 0.0);
            CHECK(utils::min(1.0, 0.0) == 0.0);
            CHECK(utils::min(0.0, -1.0) == -1.0);
            CHECK(utils::min(-1.0, 0.0) == -1.0);
        }
    }
}

TEST_CASE("[B] maxmin", "[!benchmark][utils]") {
    SECTION("Max") {
        GIVEN("ints") {
            WHEN("Benchmarking") {
                BENCHMARK("Max int (0,0)")utils::max(0, 0);
                BENCHMARK("Max int (0,1)")utils::max(0, 1);
                BENCHMARK("Max int (1,0)")utils::max(1, 0);
                BENCHMARK("Max int (0,-1)")utils::max(0, -1);
                BENCHMARK("Max int (-1,0)")utils::max(-1, 0);
            }
        }
        GIVEN("floats") {
            WHEN("Benchmarking") {
                BENCHMARK("Max float (0,0)")utils::max(0.0F, 0.0F);
                BENCHMARK("Max float (0,1)")utils::max(0.0F, 1.0F);
                BENCHMARK("Max float (1,0)")utils::max(1.0F, 0.0F);
                BENCHMARK("Max float (0,-1)")utils::max(0.0F, -1.0F);
                BENCHMARK("Max float (-1,0)")utils::max(-1.0F, 0.0F);
            }
        }
        GIVEN("doubles") {
            WHEN("Benchmarking") {
                BENCHMARK("Max float (0,0)")utils::max(0.0, 0.0);
                BENCHMARK("Max float (0,1)")utils::max(0.0, 1.0);
                BENCHMARK("Max float (1,0)")utils::max(1.0, 0.0);
                BENCHMARK("Max float (0,-1)")utils::max(0.0, -1.0);
                BENCHMARK("Max float (-1,0)")utils::max(-1.0, 0.0);
            }
        }
    }

    SECTION("Min") {
        GIVEN("ints") {
            WHEN("Benchmarking") {
                BENCHMARK("Min int (0,0)")utils::min(0, 0);
                BENCHMARK("Min int (0,1)")utils::min(0, 1);
                BENCHMARK("Min int (1,0)")utils::min(1, 0);
                BENCHMARK("Min int (0,-1)")utils::min(0, -1);
                BENCHMARK("Min int (-1,0)")utils::min(-1, 0);
            }
        }
        GIVEN("floats") {
            WHEN("Benchmarking") {
                BENCHMARK("Min float (0,0)")utils::min(0.0F, 0.0F);
                BENCHMARK("Min float (0,1)")utils::min(0.0F, 1.0F);
                BENCHMARK("Min float (1,0)")utils::min(1.0F, 0.0F);
                BENCHMARK("Min float (0,-1)")utils::min(0.0F, -1.0F);
                BENCHMARK("Min float (-1,0)")utils::min(-1.0F, 0.0F);
            }
        }
        GIVEN("doubles") {
            WHEN("Benchmarking") {
                BENCHMARK("Min float (0,0)")utils::min(0.0, 0.0);
                BENCHMARK("Min float (0,1)")utils::min(0.0, 1.0);
                BENCHMARK("Min float (1,0)")utils::min(1.0, 0.0);
                BENCHMARK("Min float (0,-1)")utils::min(0.0, -1.0);
                BENCHMARK("Min float (-1,0)")utils::min(-1.0, 0.0);
            }
        }
    }
}

#pragma clang diagnostic pop