#include "../../catch.hpp"
#include "utils.h"

TEST_CASE("readLine", "[utils]") {

    std::stringstream ss;
    std::unique_ptr<char[]> a(nullptr);

    GIVEN("Simple line")
    {
        ss << "HELLO";
        a.reset(utils::readLine(ss));
        std::string str(a.get());

        REQUIRE(str == "HELLO");
    }

    GIVEN("Two lines - Unix")
    {
        ss << "HELLO\nBYE";

        a.reset(utils::readLine(ss));
        std::string str(a.get());
        REQUIRE(str == "HELLO");

        a.reset(utils::readLine(ss));
        str = std::string(a.get());
        REQUIRE(str == "BYE");
    }

    GIVEN("Two lines with spaces on the left - Unix")
    {
        ss << "   HELLO\n   BYE";

        a.reset(utils::readLine(ss));
        std::string str(a.get());
        REQUIRE(str == "HELLO");

        a.reset(utils::readLine(ss));
        str = std::string(a.get());
        REQUIRE(str == "BYE");
    }

    GIVEN("Two lines - Windows")
    {
        ss << "HELLO\r\nBYE";

        a.reset(utils::readLine(ss));
        std::string str(a.get());
        REQUIRE(str == "HELLO");

        a.reset(utils::readLine(ss));
        str = std::string(a.get());
        REQUIRE(str == "BYE");
    }

    GIVEN("Two lines with spaces on the left - Windows")
    {
        ss << "   HELLO\r\n   BYE";

        a.reset(utils::readLine(ss));
        std::string str(a.get());
        REQUIRE(str == "HELLO");

        a.reset(utils::readLine(ss));
        str = std::string(a.get());
        REQUIRE(str == "BYE");
    }

}