#ifndef VCF_STATISH_H
#define VCF_STATISH_H

#include <string>
#include <vector>
#include <newAlignment.h>

class vcf_statish
{
public:
    
    vcf_statish();
    
    ~vcf_statish();

    void readVCF(std::vector< newAlignment*> sources,
                 std::vector<std::string> filenames,
                 float minQuality = 0,
                 float minCoverage = 0,
                 bool ignoreFilter = false,
                 const char * const replacementChar = NULL);
private:
    
    void printApeek(std::vector<newAlignment *> & sources);
    
    void applyVariantCallingFiles(
        std::vector<newAlignment*>      & sources,
        std::vector<std::string>        & filenames,
        std::vector<std::string>        & contigs,
        std::vector<std::vector<int>>   & donorsPositions,
        float minQuality,
        float minCoverage,
        bool ignoreFilter,
        const char * const replacementChar);

    void extendAlignments(
        std::vector<newAlignment*>  & sources,
        std::vector<std::string>    & contigs,
        std::vector<std::string>    & donors);

    void obtainContigsAndDonors(
        std::vector<std::string>        & filenames,
        std::vector<std::string>        & donors,
        std::vector<std::string>        & contigs,
        std::vector<std::vector<int>>   & donorsPositions);

    int getTagFromCharArray(char * array, char separator);

    char getCharFromTag(int tag);
};

#endif // VCF_STATISH_H
