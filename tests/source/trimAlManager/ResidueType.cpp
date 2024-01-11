#include "trimAlManagerIncludes.h"

using namespace catch2Utils::trimAlManager;
using namespace trimAlManagerTestIncludes;

TEST_CASE("Check Alignment Type") {

    catch2Utils::init();
    trimAlManager manager;
    
    GIVEN("Amino acid MSA")
    {
        manager.origAlig = manager.getFormatManager().loadAlignment(aminoAcidAlignmentPath);
        CHECK(manager.origAlig->getAlignmentType() == SequenceTypes::AA);
    }

    GIVEN("DNA MSA")
    {
        manager.origAlig = manager.getFormatManager().loadAlignment(dnaAlignmentPath);
        CHECK(manager.origAlig->getAlignmentType() == SequenceTypes::DNA);
    }

    GIVEN("RNA MSA")
    {
        manager.origAlig = manager.getFormatManager().loadAlignment(rnaAlignmentPath);
        CHECK(manager.origAlig->getAlignmentType() == SequenceTypes::RNA);
    }

    GIVEN("Ambiguous amino acid MSA")
    {
        manager.origAlig = manager.getFormatManager().loadAlignment(ambiguousAminoAcidAlignmentPath);
        CHECK(manager.origAlig->getAlignmentType() == (SequenceTypes::AA | SequenceTypes::DEG));
    }

    GIVEN("Degenerate DNA MSA")
    {
        manager.origAlig = manager.getFormatManager().loadAlignment(degenerateDnaAlignmentPath);
        CHECK(manager.origAlig->getAlignmentType() == (SequenceTypes::DNA | SequenceTypes::DEG));
    }

    GIVEN("Degenerate RNA MSA")
    {
        manager.origAlig = manager.getFormatManager().loadAlignment(degenerateRnaAlignmentPath);
        CHECK(manager.origAlig->getAlignmentType() == (SequenceTypes::RNA | SequenceTypes::DEG));
    }

    GIVEN("Alternative amino acid MSA")
    {
        manager.origAlig = manager.getFormatManager().loadAlignment(alternativeAminoAcidAlignmentPath);
        CHECK(manager.origAlig->getAlignmentType() == SequenceTypes::AA);
    }
}
