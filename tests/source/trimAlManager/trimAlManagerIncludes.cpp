#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err58-cpp"
#include "trimAlManagerIncludes.h"

namespace trimAlManagerTestIncludes
{
    const char * const alignmentPath =
            "../dataset/example.004.AA.fasta";

    const char * const aminoAcidAlignmentPath =
            "../dataset/example.005.AA.fasta";

    const char * const dnaAlignmentPath =
            "../dataset/example.092.DNA.fasta";

    const char * const rnaAlignmentPath =
            "../dataset/example.096.RNA.fasta";

    const char * const ambiguousAminoAcidAlignmentPath =
            "../dataset/example.097.ambiguous.AA.fasta";

    const char * const degenerateDnaAlignmentPath =
            "../dataset/example.098.deg.DNA.fasta";

    const char * const degenerateRnaAlignmentPath =
            "../dataset/example.099.deg.RNA.fasta";

    const char * const alternativeAminoAcidAlignmentPath =
            "../dataset/example.100.alt.AA.fasta";

    namespace formats
    {
        std::vector<const char *> legacyFormats {
                "-nbrf", "-mega", "-nexus", "-clustal", "-fasta", "-fasta_m10", "-phylip",
                "-phylip_m10", "-phylip_paml", "-phylip_paml_m10", "-phylip3.2", "-phylip3.2_m10"
        };

        std::vector<const char *> nonLegacyFormats {
                "clustal", "fasta_m10", "fasta", "html", "mega_sequential",
                "nexus_m10", "nexus", "phylip32_m10", "phylip32",
                "phylip40_m10", "phylip40", "phylip_paml_m10", "phylip_paml", "pir"
        };
    }

    namespace backtranslation
    {
        const char * backtranslationAlignmentPath =
                "../dataset/example.091.AA.strNOG.ENOG411BWBU.fasta";
        const char * backtranslationCodonPath =
                "../dataset/example.091.AA.strNOG.ENOG411BWBU.codon.fa";
    }

    const std::vector<argumentReference>
        getAutomatedMethods(trimAlManager & manager)
    {
        return
        {
                {"-strict",     manager.strict},
                {"-strictplus", manager.strictplus},
                {"-gappyout",   manager.gappyout},
                {"-automated1", manager.automated1},
                {"-nogaps",     manager.nogaps},
                {"-noallgaps",  manager.noallgaps},
        };
    }

    const std::vector<argumentReference>
        getStatsArguments(trimAlManager& manager, bool includeGeneral, bool includeSpecial)
    {
        assert(includeGeneral || includeSpecial);

        if (includeGeneral)
        {
            if (!includeSpecial) return
            {
                {"-sgc",        manager.sgc},     // Gap per column
                {"-sgt",        manager.sgt},     // Gap accumulated
                {"-ssc",        manager.ssc},     // Similarity per column
                {"-sst",        manager.sst},     // Similarity accumulated
                {"-sident",     manager.sident},  // Identity scores
                {"-soverlap",   manager.soverlap} // Overlap score matrix
            };
            else return
            {
                {"-sgc",        manager.sgc},      // Gap per column
                {"-sgt",        manager.sgt},      // Gap accumulated
                {"-ssc",        manager.ssc},      // Similarity per column
                {"-sst",        manager.sst},      // Similarity accumulated
                {"-sident",     manager.sident},   // Identity scores
                {"-soverlap",   manager.soverlap}, // Overlap score matrix
                {"-sfc",        manager.sfc},      // Consistency per column
                {"-sft",        manager.sft}       // Consistency accumulated
            };
        } else {
            if (includeSpecial) return
            {
                {"-sfc", manager.sfc},            // Consistency per column
                {"-sft", manager.sft}             // Consistency accumulated
            };
        }
    }

    const std::vector<selectXCase>
        getSelectX(trimAlManager &manager, Alignment *alignment)
    {
        if (alignment != nullptr) return
        {
                {"Select Seqs", "-selectseqs", alignment->numberOfSequences, &manager.delSequences},
                {"Select Cols", "-selectcols", alignment->numberOfResidues,  &manager.delColumns},
        };
        else return
        {
                {"Select Seqs", "-selectseqs", -1,  &manager.delSequences},
                {"Select Cols", "-selectcols", -1,  &manager.delColumns},
        };
    }

    const std::vector<argumentReference>
        getBooleanModifiers(trimAlManager & manager)
    {
        return
        {
            {"-terminalonly",       manager.terminalOnly},
            {"-keepseqs",           manager.keepSeqs},
            {"-keepheader",         manager.getFormatManager().keepHeader},
            {"-complementary",      manager.getComplementary},
            {"-colnumbering",       manager.columnNumbering}
        };
    }

    const std::vector<reportCase>
        getReportArguments(trimAlManager & manager)
    {
        return
        {
            {"-htmlout",    &manager.htmlOutFile},
            {"-svgout",     &manager.svgOutFile},
            {"-svgstats",   &manager.svgStatsOutFile},
        };
    }

    const std::vector<thresholds::thresholdCase>
        getThresholdCases(trimAlManager & manager, bool includeConsistency)
    {
        if (includeConsistency) return
        {
            {"Consistency", "-ct", "-cw", manager.consistencyThreshold, manager.consistencyWindow},
            {"Gaps",        "-gt", "-gw", manager.gapThreshold,         manager.gapWindow},
            {"Similarity",  "-st", "-sw", manager.similarityThreshold,  manager.similarityWindow},
        };
        else return
        {
            {"Gaps",        "-gt", "-gw", manager.gapThreshold,         manager.gapWindow},
            {"Similarity",  "-st", "-sw", manager.similarityThreshold,  manager.similarityWindow},
        };
    }
}

#pragma clang diagnostic pop