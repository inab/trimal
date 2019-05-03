#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-flp30-c"
#include "../../catch.hpp"
#include "utils.h"

TEST_CASE("isNumber", "[utils]")
{
    std::unique_ptr<char[]> tentativeNumber(nullptr);
    GIVEN("nullptr")
    {
        CAPTURE(tentativeNumber.get());
        FAIL("Fatal Error");
        CHECK_FALSE(utils::isNumber(tentativeNumber.get()));
    }
    GIVEN("'0\\0' -> true")
    {
        tentativeNumber.reset(
                new char[2] { '0', '\0'});
        CAPTURE(tentativeNumber.get());
        REQUIRE(utils::isNumber(tentativeNumber.get()));
    }
    GIVEN("'\\0' -> false")
    {
        tentativeNumber.reset(
                new char[1] { '\0'});
        CAPTURE(tentativeNumber.get());
        REQUIRE_FALSE(utils::isNumber(tentativeNumber.get()));
    }
    GIVEN("'10e10\\0' -> true")
    {
        tentativeNumber.reset(
                new char[6] { '1', '0', 'e', '1', '0', '\0'});
        CAPTURE(tentativeNumber.get());
        REQUIRE(utils::isNumber(tentativeNumber.get()));
    }
    GIVEN("'10E10\\0' -> true")
    {
        tentativeNumber.reset(
                new char[6] { '1', '0', 'E', '1', '0', '\0'});
        CAPTURE(tentativeNumber.get());
        REQUIRE(utils::isNumber(tentativeNumber.get()));
    }
    GIVEN("'10--0\\0' -> false")
    {
        tentativeNumber.reset(
                new char[6] { '1', '0', '-', '-', '0', '\0'});
        CAPTURE(tentativeNumber.get());
        REQUIRE_FALSE(utils::isNumber(tentativeNumber.get()));
    }
    GIVEN("'10EE0\\0' -> false")
    {
        tentativeNumber.reset(
                new char[6] { '1', '0', 'E', 'E', '0', '\0'});
        CAPTURE(tentativeNumber.get());
        REQUIRE_FALSE(utils::isNumber(tentativeNumber.get()));
    }
    GIVEN("'E\\0' -> false")
    {
        tentativeNumber.reset(
                new char[2] { 'E', '\0'});
        CAPTURE(tentativeNumber.get());
        REQUIRE_FALSE(utils::isNumber(tentativeNumber.get()));
    }
    GIVEN("'E1\\0' -> false")
    {
        tentativeNumber.reset(
                new char[3] { 'E', '1', '\0'});
        CAPTURE(tentativeNumber.get());
        REQUIRE_FALSE(utils::isNumber(tentativeNumber.get()));
    }
    GIVEN("'EeE-E\\0' -> false")
    {
        tentativeNumber.reset(
                new char[6] { 'E', 'e', 'E', '-', 'E', '\0'});
        CAPTURE(tentativeNumber.get());
        REQUIRE_FALSE(utils::isNumber(tentativeNumber.get()));
    }
    GIVEN("'1.1\\0' -> true")
    {
        tentativeNumber.reset(
                new char[4] { '1', '.', '1', '\0'});
        CAPTURE(tentativeNumber.get());
        REQUIRE(utils::isNumber(tentativeNumber.get()));
    }

}
#pragma clang diagnostic pop