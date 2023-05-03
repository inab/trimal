#include "../trimAlManagerIncludes.h"

using namespace catch2Utils::trimAlManager;
using namespace trimAlManagerTestIncludes::formats;

// Valgrinded

//region Providing compareset arguments
TEST_CASE("Format argument parse", "[manager][arguments][formats]") {

    catch2Utils::init();

    trimAlManager manager;

    GIVEN("Legacy Options")
    {
        for(const char* legacy_format : legacyFormats)
        {
            GIVEN(legacy_format) test_arguments(
                {"", "-in", "../dataset/example.093.DNA.fasta", legacy_format},
                manager,
                trimAlManager::argumentReport::Recognized,
                false);

            for(const char* legacy_format2 : legacyFormats)
            {
                if (legacy_format == legacy_format2) continue;

                GIVEN(legacy_format << " and " << legacy_format2)
                    test_arguments(
                        {"", "-in", "../dataset/example.093.DNA.fasta",
                         legacy_format, legacy_format2},
                        manager,
                        trimAlManager::argumentReport::Recognized,
                        false);
            }
        }
    }

    GIVEN("Non legacy Options")
    {
        for(const char* format : nonLegacyFormats)
        {
            GIVEN(format) test_arguments(
                        {"", "-in", "../dataset/example.093.DNA.fasta",
                         "-formats", format},
                        manager,
                        trimAlManager::argumentReport::Recognized,
                        false);

            for(const char* format2 : nonLegacyFormats)
            {
                if (format == format2) continue;

                GIVEN(format << " and " << format2)
                    test_arguments(
                            {"", "-in", "../dataset/example.093.DNA.fasta",
                             "-formats", format, format2},
                            manager,
                            trimAlManager::argumentReport::Recognized,
                            false);
            }
        }

        GIVEN("No Format")
        {
            test_arguments(
                    {"", "-in", "../dataset/example.093.DNA.fasta",
                     "-formats"},
                    manager,
                    trimAlManager::argumentReport::Wrong,
                    true);
        }

        GIVEN("Incorrect Format")
        {
            test_arguments(
                    {"", "-in", "../dataset/example.093.DNA.fasta",
                     "-formats", "thisisnoformat"},
                    manager,
                    trimAlManager::argumentReport::Wrong,
                    true);
        }
    }

    GIVEN("Legacy and Non-Legacy argument combinations")
    {
        for(const char* legacyFormat : legacyFormats)
        {
            for(const char* nonLegacyFormat : nonLegacyFormats)
            {
                GIVEN(legacyFormat << " and " << nonLegacyFormat)
                    test_arguments(
                            {"", "-in", "../dataset/example.093.DNA.fasta",
                             "-formats", legacyFormat, "-formats", nonLegacyFormat},
                            manager,
                            trimAlManager::argumentReport::Recognized,
                            false);
            }
        }

        for(const char* nonLegacyFormat : nonLegacyFormats)
        {
            for(const char* legacyFormat : legacyFormats)
            {
                GIVEN(nonLegacyFormat << " and " << legacyFormat)
                    test_arguments(
                            {"", "-in", "../dataset/example.093.DNA.fasta",
                             "-formats", nonLegacyFormat, legacyFormat},
                            manager,
                            trimAlManager::argumentReport::Recognized,
                            false);
            }
        }
    }

}
