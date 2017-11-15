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
};

#endif // VCF_STATISH_H
