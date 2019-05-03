#include "../../catch.hpp"
#include "utils.h"

TEST_CASE("removeSpaces", "[utils]") {

    std::unique_ptr<char[]> a(nullptr);
    std::unique_ptr<char[]> b(nullptr);

    WHEN("Provided two nullptr")
    {
        FAIL("Fatal Error");
        utils::removeSpaces(a.get(), b.get());
    }

    WHEN("Provided b as nullptr")
    {
        a.reset(new char [3]);
        FAIL("Fatal Error");
        utils::removeSpaces(a.get(), nullptr);
    }

    WHEN("Provided a as nullptr")
    {
        b.reset(new char [3]);
        FAIL("Fatal Error");
        utils::removeSpaces(a.get(), nullptr);
    }

    WHEN("Provided only spaces")
    {
        a.reset(new char[2]{' ', '\0'});
        b.reset(new char[2]{'-', '\0'});
        utils::removeSpaces(a.get(), b.get());
        REQUIRE(b[0] == '\0');
    }

    WHEN("Provided only tabs")
    {
        a.reset(new char[2]{'\t', '\0'});
        b.reset(new char[2]{'-', '\0'});
        utils::removeSpaces(a.get(), b.get());
        REQUIRE(b[0] == '\0');
    }


    WHEN("Provided non-spaces string")
    {
        a.reset(new char[2]{'-', '\0'});
        b.reset(new char[2]{' ', '\0'});
        utils::removeSpaces(a.get(), b.get());
        REQUIRE(b[0] == '-');
        REQUIRE(b[1] == '\0');
    }

    WHEN("Provided mixed string letters and spaces")
    {
        a.reset(new char[5]{'-', ' ', 'a', ' ', '\0'});
        b.reset(new char[5]{' ', '\0'});
        utils::removeSpaces(a.get(), b.get());
        REQUIRE(b[0] == '-');
        REQUIRE(b[1] == 'a');
        REQUIRE(b[2] == '\0');
    }

    WHEN("Provided mixed string letters and tabs")
    {
        a.reset(new char[5]{'-', '\t', 'a', '\t', '\0'});
        b.reset(new char[5]{' ', '\0'});
        utils::removeSpaces(a.get(), b.get());
        REQUIRE(b[0] == '-');
        REQUIRE(b[1] == 'a');
        REQUIRE(b[2] == '\0');
    }

    WHEN("Provided mixed string letters, spaces and tabs")
    {
        a.reset(new char[5]{'-', '\t', 'a', ' ', '\0'});
        b.reset(new char[5]{' ', '\0'});
        utils::removeSpaces(a.get(), b.get());
        REQUIRE(b[0] == '-');
        REQUIRE(b[1] == 'a');
        REQUIRE(b[2] == '\0');
    }
}

TEST_CASE("[B] removeSpaces", "[!benchmark][utils]") {

    std::unique_ptr<char[]> a(nullptr);
    std::unique_ptr<char[]> b(nullptr);

    WHEN("Benchmarking")
    {
        GIVEN("Random elements")
        {
            for (int i = 1; i < 10; i++)
            {
                int size = (int)pow(2, i);
                a.reset(new char[size]);
                b.reset(new char[size]);
                BENCHMARK("Remove spaces with " + std::to_string(size) + " random elements")
                {
                    utils::removeSpaces(a.get(), b.get());
                }
            }
        }

        GIVEN("All-spaces elements")
        {
            for (int i = 1; i < 10; i++)
            {
                int size = (int)pow(2, i);
                a.reset(new char[size]);
                std::fill(a.get(), a.get() + size - 1, ' ');
                a[size - 1] = '\0';
                b.reset(new char[size]);
                BENCHMARK("Remove spaces with " + std::to_string(size) + " all-spaces elements")
                {
                    utils::removeSpaces(a.get(), b.get());
                }
            }
        }

        GIVEN("No-spaces elements")
        {
            for (int i = 1; i < 10; i++)
            {
                int size = (int)pow(2, i);
                a.reset(new char[size]);
                std::fill(a.get(), a.get() + size - 1, '-');
                a[size - 1] = '\0';
                b.reset(new char[size]);
                BENCHMARK("Remove spaces with " + std::to_string(size) + " no-spaces elements")
                {
                    utils::removeSpaces(a.get(), b.get());
                }
            }
        }
    }
}