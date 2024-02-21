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
                true, trimAlManager::argumentReport::Errored,
                true, true);
        GIVEN("Quality value")
        {
            float qualityValue = 0.5F;
            args.push_back("0.5");
            test_arguments(
                    args, manager,
                    true, expectedReport,
                    true, expectedErrors,
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
                true, trimAlManager::argumentReport::Errored,
                true, true);
        GIVEN("Coverage value")
        {
            float coverageValue = 0.5F;
            args.push_back("0.5");
            test_arguments(
                    args, manager,
                    true, expectedReport,
                    true, expectedErrors,
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
                true, expectedReport,
                true, expectedErrors,
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
                true, trimAlManager::argumentReport::Errored,
                true, true);
        GIVEN("VCF file")
        {
            const char * vcfFile = "/home/vfernandez/Dropbox/VCF/B1012M.bam.flt.vcf";
            args.push_back(vcfFile);
            test_arguments(
                    args, manager,
                    true, expectedReport,
                    true, expectedErrors,
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
        const bool expectedErrors, uint32_t iterator = 0)
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
            applyCalls(models, args, manager, trimAlManager::argumentReport::Errored, true);
        }

        AddVCF(args, manager, trimAlManager::argumentReport::Errored, true,
               [&](){ applyCalls(models, args, manager, trimAlManager::argumentReport::Errored, true); });
    }

    GIVEN("Input")
    {
        args.push_back("-in");
        args.push_back("/home/vfernandez/git/PreTFM/NGS_trimAl/dataset/ngs/CANGA_REF.fasta");

        GIVEN("No VCF")
            applyCalls(models, args, manager, trimAlManager::argumentReport::Recognized, true);

        AddVCF(args, manager, trimAlManager::argumentReport::Recognized, false,
                [&](){ applyCalls(models, args, manager, trimAlManager::argumentReport::Recognized, false); });
    }


}
