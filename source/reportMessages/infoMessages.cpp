#include "reportsystem.h"

const std::map<InfoCode, const char *> reporting::reportManager::InfoMessages =
{
    {InfoCode::CuttingSequence,
            "Cutting sequence \"[tag]\" at first appearance of stop codon \"[tag]\""
                    " (residue \"[tag]\") at position [tag] (length: [tag] \")"},

    {InfoCode::WindowSizeCompareset,
            "Window size (-w) is provided. "
                    "It's recommended to use specific consistency window size (-cw)"
                    " when using -compareset option"},
    {InfoCode ::AddingSNP,
            "Applying SNP to \"[tag]\":\"[tag]\" at position [tag] \"[tag]\"->\"[tag]\""
    }
};