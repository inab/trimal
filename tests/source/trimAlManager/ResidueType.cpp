#include "trimAlManagerIncludes.h"

using namespace catch2Utils::trimAlManager;
using namespace trimAlManagerTestIncludes;

TEST_CASE("Check Alignment Type") {

    catch2Utils::init();
    trimAlManager manager;
    
    GIVEN("Amino acid residues")
    {
        manager.origAlig = manager.getFormatManager().loadAlignment(aminoAcidAlignmentPath);
        CHECK(manager.origAlig->getAlignmentType() == SequenceTypes::AA);
    }

    GIVEN("DNA residues")
    {
        manager.origAlig = manager.getFormatManager().loadAlignment(dnaAlignmentPath);
        CHECK(manager.origAlig->getAlignmentType() == SequenceTypes::AA);
        // add warning check
    }

    GIVEN("RNA residues")
    {
        manager.origAlig = manager.getFormatManager().loadAlignment(rnaAlignmentPath);
        CHECK(manager.origAlig->getAlignmentType() == SequenceTypes::RNA);
    }
}
