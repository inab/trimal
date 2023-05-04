#include "../../catch.hpp"
#include "utils.h"

TEST_CASE("stringManipulation", "[utils]") {

    WHEN("Trimming line") {
        std::unique_ptr<char[]> ptr;

        GIVEN("Line without comments") {
            std::string str("Hello");
            ptr.reset(utils::trimLine(str));
            CAPTURE(str);
            REQUIRE(ptr.get() == std::string("Hello"));
        }

        GIVEN("Line with a [] comment") {
            std::string str("Hello[my] darling");
            ptr.reset(utils::trimLine(str));
            CAPTURE(str);
            REQUIRE(ptr.get() == std::string("Hello darling"));
        }

        GIVEN("Line with two [] comments") {
            std::string str("Hello[my] darling[.]");
            ptr.reset(utils::trimLine(str));
            CAPTURE(str);
            REQUIRE(ptr.get() == std::string("Hello darling"));
        }

        GIVEN("Line with a \"\" comment") {
            std::string str("Hello\"my\" darling\".\"");
            ptr.reset(utils::trimLine(str));
            CAPTURE(str);
            REQUIRE(ptr.get() == std::string("Hello darling"));
        }

        GIVEN("Line with two \"\" comments and two [] comments") {
            std::string str("Hello\"my\"[my] darling\".\"[.]");
            ptr.reset(utils::trimLine(str));
            CAPTURE(str);
            REQUIRE(ptr.get() == std::string("Hello darling"));
        }

        GIVEN("Line with double comments [\"\"] ") {
            std::string str("Hello[\"my\"]");
            ptr.reset(utils::trimLine(str));
            CAPTURE(str);
            REQUIRE(ptr.get() == std::string("Hello"));
        }

        GIVEN("Line with double comments \"[]\" ") {
            std::string str("Hello\"[my]\"");
            ptr.reset(utils::trimLine(str));
            CAPTURE(str);
            REQUIRE(ptr.get() == std::string("Hello"));
        }

        GIVEN("Line with (malformed) double comments \"[\"] ") {
            std::string str("Hello\"[my\"]");
            ptr.reset(utils::trimLine(str));
            CAPTURE(str);
            REQUIRE(ptr == nullptr);
        }

        GIVEN("Line with (malformed) double comments [\"]\" ") {
            std::string str("Hello[\"my]\"");
            ptr.reset(utils::trimLine(str));
            CAPTURE(str);
            REQUIRE(ptr == nullptr);
        }
    }

    WHEN("Reversing line") {
        GIVEN("ABDCEF") {
            std::string str("ABCDEF");
            CAPTURE(str);
            CHECK(utils::getReverse(str) == "FEDCBA");
        }

        GIVEN("One letter string") {
            std::string str("A");
            CAPTURE(str);
            CHECK(utils::getReverse(str) == "A");
        }

        GIVEN("Empty string") {
            std::string str("");
            CAPTURE(str);
            CHECK(utils::getReverse(str).empty());
        }
    }

    WHEN("Removing character from line") {
        GIVEN("ABDCEF - 'A'") {
            std::string str("ABCDEF");
            CAPTURE(str);
            CHECK(utils::removeCharacter('A', str) == "BCDEF");
        }

        GIVEN("ABDCEF - 'F'") {
            std::string str("ABCDEF");
            CAPTURE(str);
            CHECK(utils::removeCharacter('F', str) == "ABCDE");
        }

        GIVEN("ABCDEF - '-'") {
            std::string str("ABCDEF");
            CAPTURE(str);
            CHECK(utils::removeCharacter('-', str) == "ABCDEF");
        }

        GIVEN(" \"\" - '-'") {
            std::string str;
            REQUIRE(str.empty());
            CAPTURE(str);
            CHECK(utils::removeCharacter('-', str).empty());
        }
    }

    WHEN("Counting character in a line") {
        GIVEN("ABCDEF - 'A'") {
            std::string str("ABCDEF");
            CAPTURE(str);
            CHECK(utils::countCharacter('A', str) == 1);
        }

        GIVEN("ABCDEF - 'H'") {
            std::string str("ABCDEF");
            CAPTURE(str);
            CHECK(utils::countCharacter('H', str) == 0);
        }

        GIVEN("ABABBABB - 'A'") {
            std::string str("ABABBABB");
            CAPTURE(str);
            CHECK(utils::countCharacter('A', str) == 3);
        }

        GIVEN("ABABBABB - 'B'") {
            std::string str("ABABBABB");
            CAPTURE(str);
            CHECK(utils::countCharacter('B', str) == 5);
        }
    }

    WHEN("Reading numbers") {
        std::unique_ptr<int[]> ptr;
        GIVEN("Empty string") {
            std::string str;
            CAPTURE(str);
            REQUIRE(str.empty());
            ptr.reset(utils::readNumbers(str));
            {
                CAPTURE(ptr.get()[0]);
                REQUIRE(ptr.get()[0] == 0);
            }
        }

        GIVEN("String with number") {
            std::string str("1");
            CAPTURE(str);
            REQUIRE(!str.empty());
            ptr.reset(utils::readNumbers(str));
            {
                REQUIRE(ptr.get()[0] == 2);
                CHECK(ptr.get()[1] == 1);
                CHECK(ptr.get()[2] == 1);
            }
        }

        GIVEN("String with number list") {
            std::string str("1, 2, 3, 4");
            CAPTURE(str);
            REQUIRE(!str.empty());
            ptr.reset(utils::readNumbers(str));

            REQUIRE(ptr.get()[0] == 8);
            CHECK(ptr.get()[1] == 1);
            CHECK(ptr.get()[2] == 1);
            CHECK(ptr.get()[3] == 2);
            CHECK(ptr.get()[4] == 2);
            CHECK(ptr.get()[5] == 3);
            CHECK(ptr.get()[6] == 3);
            CHECK(ptr.get()[7] == 4);
            CHECK(ptr.get()[8] == 4);
        }

        GIVEN("String with broken number list") {
            std::string str("1, 2, 3, ");
            CAPTURE(str);
            REQUIRE(!str.empty());
            ptr.reset(utils::readNumbers(str));

            CHECK(ptr.get()[0] == 6);
            REQUIRE(ptr.get()[0] > 5);
            CHECK(ptr.get()[1] == 1);
            CHECK(ptr.get()[2] == 1);
            CHECK(ptr.get()[3] == 2);
            CHECK(ptr.get()[4] == 2);
            CHECK(ptr.get()[5] == 3);
            CHECK(ptr.get()[6] == 3);
            if (ptr.get()[0] == 8)
            {
                CAPTURE(ptr.get()[7]);
                CAPTURE(ptr.get()[8]);
                FAIL("There are more number positions than numbers");
            }
        }

        GIVEN("String with range") {
            std::string str("1-4");
            CAPTURE(str);
            REQUIRE(!str.empty());
            ptr.reset(utils::readNumbers(str));

            REQUIRE(ptr.get()[0] == 2);
            CHECK(ptr.get()[1] == 1);
            CHECK(ptr.get()[2] == 4);
        }

        GIVEN("String with broken range") {
            std::string str("1-");
            CAPTURE(str);
            REQUIRE(!str.empty());
            ptr.reset(utils::readNumbers(str));

            REQUIRE(ptr.get() == nullptr);
        }

        GIVEN("String with ranges and numbers") {
            std::string str("1-4, 10, 20-30");
            CAPTURE(str);
            REQUIRE(!str.empty());
            ptr.reset(utils::readNumbers(str));

            REQUIRE(ptr.get()[0] == 6);
            CHECK(ptr.get()[1] == 1);
            CHECK(ptr.get()[2] == 4);

            CHECK(ptr.get()[3] == 10);
            CHECK(ptr.get()[4] == 10);

            CHECK(ptr.get()[5] == 20);
            CHECK(ptr.get()[6] == 30);
        }
    }

    WHEN("Replacing in strings")
    {
        WHEN("In Place")
        {
            GIVEN("A string with a [tag]")
            {
                std::string str("A string with a [tag]");
                utils::ReplaceStringInPlace(str, "[tag]", "removed tag");
                REQUIRE(str == "A string with a removed tag");
            }

            GIVEN("A string with no tag")
            {
                std::string str("A string with no tag");
                utils::ReplaceStringInPlace(str, "[tag]", "removed tag");
                REQUIRE(str == "A string with no tag");
            }
        }

        WHEN("Keeping original")
        {
            GIVEN("A string with a [tag]")
            {
                std::string str("A string with a [tag]");
                std::string res =
                        utils::ReplaceString(str, "[tag]", "removed tag");
                REQUIRE(res == "A string with a removed tag");
                REQUIRE(str == "A string with a [tag]");
            }

            GIVEN("A string with no tag")
            {
                std::string str("A string with no tag");
                std::string res =
                        utils::ReplaceString(str, "[tag]", "removed tag");
                REQUIRE(str == "A string with no tag");
                REQUIRE(&str != &res);
            }
        }
    }

}