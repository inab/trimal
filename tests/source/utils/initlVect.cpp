#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-flp30-c"
#include "../../catch.hpp"
#include "utils.h"

TEST_CASE("initlVect", "[utils]")
{
    GIVEN("An int vector")
    {
        std::unique_ptr<int[]> a(nullptr);
        WHEN("Vector is null")
        {
            WHEN("Provided size 0")
            {
                CHECK_NOTHROW(utils::initlVect(a.get(), 0, 0));
            }
            WHEN("Provided size 0")
            {
                FAIL("Fatal Error");
                CHECK_NOTHROW(utils::initlVect(a.get(), 0, 0));
            }
        }
        WHEN("Vector is not null")
        {
             a.reset(new int[2]);
            utils::initlVect(a.get(), 2, 1);
            THEN("Filling it fully")
            {
                CHECK(a[0] == 1);
                CHECK(a[1] == 1);
            }
            THEN("Filling it partially")
            {
                utils::initlVect(a.get(), 1, 2);
                CHECK(a[0] == 2);
                CHECK(a[1] == 1);
            }

        }
    }

    GIVEN("A float vector")
    {
        std::unique_ptr<float[]> a(nullptr);
        WHEN("Vector is null")
        {
            WHEN("Provided size 0")
            {
                CHECK_NOTHROW(utils::initlVect(a.get(), 0, 0));
            }
        }
        WHEN("Vector is not null")
        {
            a.reset(new float[2]);
            utils::initlVect(a.get(), 2, 1);
            THEN("Filling it fully")
            {
                CHECK(a[0] == 1);
                CHECK(a[1] == 1);
            }
            THEN("Filling it partially")
            {
                utils::initlVect(a.get(), 1, 2);
                CHECK(a[0] == 2);
                CHECK(a[1] == 1);
            }

        }
    }
}

TEST_CASE("[B] initlVect", "[!benchmark][utils]")
{
    GIVEN("An int vector")
    {
        std::unique_ptr<int[]> a(nullptr);

        for (int i = 1; i < 10; i++)
        {
            int size = (int)pow(2, i);
            a.reset(new int[size]);
            BENCHMARK("Filling int vector with " + std::to_string(size) + " elements")
            {
                utils::initlVect(a.get(), size, 0);
            }
        }
    }

    GIVEN("A float vector")
    {
        std::unique_ptr<float[]> a(nullptr);

        for (int i = 1; i < 10; i++)
        {
            int size = (int)pow(2, i);
            a.reset(new float[size]);
            BENCHMARK("Filling float vector with " + std::to_string(size) + " elements")
            {
                utils::initlVect(a.get(), size, 0);
            }
        }

    }
}





#pragma clang diagnostic pop