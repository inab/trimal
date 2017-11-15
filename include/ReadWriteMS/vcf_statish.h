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
    
    void readVCF(std::vector< newAlignment*> sources, std::vector<std::string> filenames, float minQuality = 0, float minCoverage = 0, bool ignoreFilter = false);
};

class vcf_file
{
public:
    
    class snp_entry
    {
    public:
        char *  chromosome;
        int     position;
        char *  ref;
        char *  alt;
        float   quality;
        bool    filter;
        std::vector<int> present_in = std::vector<int>();
    };
    
    
    std::vector<char *>     filenames   = std::vector<char *>();
    std::vector<char *>     entries     = std::vector<char *>();
    std::vector<snp_entry>  snps        = std::vector<snp_entry>();
};

#endif // VCF_STATISH_H
