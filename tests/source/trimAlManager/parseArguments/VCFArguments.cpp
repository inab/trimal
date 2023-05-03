#include "../trimAlManagerIncludes.h"

using namespace catch2Utils::trimAlManager;
using namespace trimAlManagerTestIncludes;

typedef std::function<void(
        std::vector<const char *>& args,
        trimAlManager & manager,
        const trimAlManager::argumentReport expectedReport,
        const bool expectedErrors,
        const std::function<void(void)>& callback)> model;

void AddQuality(
        std::vector<const char *>& args,
        trimAlManager & manager,
        const trimAlManager::argumentReport expectedReport,
        const bool expectedErrors,
        const std::function<void(void)>& callback = {})
{
    GIVEN("Quality argument")
    {
        args.push_back("-minquality");
        test_arguments(
                args, manager,
                trimAlManager::argumentReport::Wrong,
                true);
        GIVEN("Quality value")
        {
            float qualityValue = 0.5F;
            args.push_back("0.5");
            test_arguments(
                    args, manager,
                    expectedReport,
                    expectedErrors,
                    [&manager, qualityValue, expectedErrors, expectedReport](){
                        if (!expectedErrors && expectedReport != trimAlManager::argumentReport::Recognized)
                        THEN("Value is set correctly in the manager")
                            CHECK(manager.minQuality == qualityValue);
                    });

            if (callback) callback();
        }
    }
}

void AddCoverage(
        std::vector<const char *>& args,
        trimAlManager & manager,
        const trimAlManager::argumentReport expectedReport,
        const bool expectedErrors,
        const std::function<void(void)>& callback = {})
{
    GIVEN("Coverage argument")
    {
        args.push_back("-mincoverage");
        test_arguments(
                args, manager,
                trimAlManager::argumentReport::Wrong,
                true);
        GIVEN("Coverage value")
        {
            float coverageValue = 0.5F;
            args.push_back("0.5");
            test_arguments(
                    args, manager,
                    expectedReport,
                    expectedErrors,
                    [&manager, coverageValue, expectedErrors, expectedReport](){
                        if (!expectedErrors && expectedReport != trimAlManager::argumentReport::Recognized)
                        THEN("Value is set correctly in the manager")
                            CHECK(manager.minCoverage == coverageValue);
                    });

            if (callback) callback();
        }
    }
}

void AddIgnoreFilter(
        std::vector<const char *>& args,
        trimAlManager & manager,
        const trimAlManager::argumentReport expectedReport,
        const bool expectedErrors,
        const std::function<void(void)>& callback = {})
{
    GIVEN("IgnoreFilter")
    {
        args.push_back("-ignorefilter");
        test_arguments(
                args, manager,
                expectedReport,
                expectedErrors,
                [&manager, expectedErrors, expectedReport](){
                    if (!expectedErrors && expectedReport != trimAlManager::argumentReport::Recognized)
                        THEN("Value is set correctly in the manager")
                            CHECK(manager.ignoreFilter);
                });
        if (callback) callback();
    }

    GIVEN("No Ignorefilter")
    {
        THEN("Value is set correctly in the manager")
            CHECK_FALSE(manager.ignoreFilter);
        if (callback) callback();
    }
}

void AddVCF(
        std::vector<const char *>& args,
        trimAlManager & manager,
        const trimAlManager::argumentReport expectedReport,
        const bool expectedErrors,
        const std::function<void(void)>& callback = {})
{
    GIVEN("VCF argument")
    {
        args.push_back("-vcf");
        test_arguments(
                args, manager,
                trimAlManager::argumentReport::Wrong,
                true);
        GIVEN("VCF file")
        {
            const char * vcfFile = "../dataset/GCA_000247815.2_current_ids.vcf";
            args.push_back(vcfFile);
            test_arguments(
                    args, manager,
                    expectedReport,
                    expectedErrors,
                    [&manager, &vcfFile, expectedErrors, expectedReport](){
                        if (!expectedErrors && expectedReport != trimAlManager::argumentReport::Recognized)
                            THEN("Value is set correctly in the manager")
                            {
                                REQUIRE(manager.vcfs->size() == 1);
                                CHECK(strcmp(manager.vcfs->at(0).c_str(), vcfFile) == 0);
                            }
                    });

            if (callback) callback();
        }
    }
}

void applyCalls(
        std::vector<model>& models,
        std::vector<const char *>& args,
        trimAlManager & manager,
        const trimAlManager::argumentReport expectedReport,
        const bool expectedErrors, uint iterator = 0)
{
    if (iterator < models.size())
    {
        models.at(iterator)(args, manager, expectedReport, expectedErrors,
                [&](){ applyCalls(models, args, manager, expectedReport, expectedErrors, iterator + 1); });
        applyCalls(models, args, manager, expectedReport, expectedErrors, iterator + 1);
    }
}


TEST_CASE("VCF Arguments", "[manager][arguments][vcf][ngs]") {
    catch2Utils::init();

    trimAlManager manager;

    std::vector<const char *> args {""};

    std::vector<model> models {AddIgnoreFilter, AddCoverage, AddQuality};

    GIVEN("No Input")
    {
        GIVEN("No VCF")
        {
            applyCalls(models, args, manager, trimAlManager::argumentReport::Wrong, true);
        }

        AddVCF(args, manager, trimAlManager::argumentReport::Wrong, true,
               [&](){ applyCalls(models, args, manager, trimAlManager::argumentReport::Wrong, true); });
    }

    GIVEN("Input")
    {
        args.push_back("-in");
        args.push_back("../dataset/example.091.AA.strNOG.ENOG411BWBU.fasta");

        GIVEN("No VCF")
            applyCalls(models, args, manager, trimAlManager::argumentReport::Recognized, true);

        AddVCF(args, manager, trimAlManager::argumentReport::Recognized, false,
                [&](){ applyCalls(models, args, manager, trimAlManager::argumentReport::Recognized, false); });
    }


}
