#include "../../catch.hpp"
#include "utils.h"

TEST_CASE("readLine", "[utils]") {

    std::stringstream ss;
    std::unique_ptr<char[]> a(nullptr);

    GIVEN("Simple line")
    {
        ss << "HELLO";
        std::string buffer;
        a.reset(utils::readLine(ss, buffer));
        std::string str(a.get());

        REQUIRE(str == "HELLO");
    }

    GIVEN("Two lines - Unix")
    {
        ss << "HELLO\nBYE";

        std::string buffer;
        a.reset(utils::readLine(ss, buffer));
        std::string str(a.get());
        REQUIRE(str == "HELLO");

        a.reset(utils::readLine(ss, buffer));
        str = std::string(a.get());
        REQUIRE(str == "BYE");
    }

    GIVEN("Two lines with spaces on the left - Unix")
    {
        ss << "   HELLO\n   BYE";

        std::string buffer;
        a.reset(utils::readLine(ss, buffer));
        std::string str(a.get());
        REQUIRE(str == "HELLO");

        a.reset(utils::readLine(ss, buffer));
        str = std::string(a.get());
        REQUIRE(str == "BYE");
    }

    GIVEN("Two lines - Windows")
    {
        ss << "HELLO\r\nBYE";

        std::string buffer;
        a.reset(utils::readLine(ss, buffer));
        std::string str(a.get());
        REQUIRE(str == "HELLO");

        a.reset(utils::readLine(ss, buffer));
        str = std::string(a.get());
        REQUIRE(str == "BYE");
    }

    GIVEN("Two lines with spaces on the left - Windows")
    {
        ss << "   HELLO\r\n   BYE";

        std::string buffer;
        a.reset(utils::readLine(ss, buffer));
        std::string str(a.get());
        REQUIRE(str == "HELLO");

        a.reset(utils::readLine(ss, buffer));
        str = std::string(a.get());
        REQUIRE(str == "BYE");
    }

}