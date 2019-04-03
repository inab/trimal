#include "../../catch.hpp"
#include "utils.h"

TEST_CASE("compareCharArray", "[utils]") {

    std::unique_ptr<char[]> a(nullptr), b(nullptr);

    GIVEN("Two nullptr arrays")
    {
        CAPTURE(a.get(), b.get());
        FAIL("Fatal error");
        CHECK(utils::compare(a.get(), b.get()));
    }

    GIVEN("One nullptr array and one non nullptr")
    {
        a.reset(new char[2] {'a', '\0'});
        CAPTURE(a.get(), b.get());
        FAIL("Fatal error");
        CHECK_FALSE(utils::compare(a.get(), b.get()));
    }

    GIVEN("Two non nullptr arrays")
    {
        WHEN("Both arrays are the same")
        {
            a.reset(new char[2] {'a', '\0'});
            b.reset(new char[2] {'a', '\0'});
            CAPTURE(a.get(), b.get());
            CHECK(utils::compare(a.get(), b.get()));
        }

        WHEN("A starts as B, but is shorter")
        {
            a.reset(new char[2] {'a', '\0'});
            b.reset(new char[4] {'a', 'n', 'b', '\0'});
            CAPTURE(a.get(), b.get());
            CHECK_FALSE(utils::compare(a.get(), b.get()));
        }

        WHEN("B starts as A, but is shorter")
        {
            a.reset(new char[4] {'a', 'n', 'b', '\0'});
            b.reset(new char[2] {'a', '\0'});
            CAPTURE(a.get(), b.get());
            CHECK_FALSE(utils::compare(a.get(), b.get()));
        }
    }
}