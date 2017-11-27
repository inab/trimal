#include "../../include/Catch/catch.hpp"

#include "../../include/trimalArgumentParser.h"
#include "../../include/reportsystem.h"
#include "../../source/Catch/Matchers/ArrayMatcher.cpp"

SCENARIO ( "Alignment methods work correctly", "[alignment][aligMethods]" ) {
    newAlignment alig;

    debug.Level = VerboseLevel::NONE;

    alig.sequenNumber = 6;
    alig.originalSequenNumber = 6;

    alig.residNumber = 60;
    alig.originalResidNumber = 60;

    alig.sequences = new string [6] {
        "--------FAYTAPD---LLLIGFLLKTVA-T-FG--DTWF-----QLWQGLDLNKMPVF",
        "----------DPAVL----FV--IMLGTIT-K-FS--SEWF-----FAWLGLEINMMVII",
        "----------GLGKV---IVY-GIVLGTKS-DQFSNWVVWL-----FPWNGLQIHMMGII",
        "-----------PTIL---NIA-GLHMETDI-N-FS--LAWF-----QAWGGLEINKQAIL",
        "----------ASGAI---LTL-GIYLFTLC-AVIS--VSWY-----LAWLGLEINMMAII",
        "AAAAAAAA----ALL---TYL-GLFLGTDY-----EN---FAAAAANAWLGLEINMMAQI"
    };

    alig.seqsName = new string [6] {
        "Sequence 0",
        "Sequence 1",
        "Sequence 2",
        "Sequence 3",
        "Sequence 4",
        "Sequence 5"
    };

    alig.fillMatrices ( true );

    int * saveResidues = new int[60];
    for ( int i = 0; i < 60; i++ )
        saveResidues[i] = i;

    REQUIRE_THAT ( alig.saveResidues, ArrayContentsEqual ( saveResidues, alig.residNumber ) );

    delete [] saveResidues;

    WHEN ( "Testing Alignment::getTranslationCDS" ) {
        WARN ( "Not implemented" );
    }

    WHEN ( "Getting number of species" ) {
        REQUIRE ( alig.getNumSpecies() == 6 );
    }

    WHEN ( "Getting number of residues" ) {
        REQUIRE ( alig.getNumAminos() == 60 );
    }

    WHEN ( "Getting number of residues" ) {
        REQUIRE ( alig.getNumAminos() == 60 );
    }

    WHEN ( "Setting Windows Size" ) {
        alig.setWindowsSize ( 10, 20 );
        REQUIRE ( alig.Statistics->ghWindow == 10 );
        REQUIRE ( alig.Statistics->shWindow == 20 );
    }

    WHEN ( "Setting Block Size" ) {
        alig.setBlockSize ( 40 );
        REQUIRE ( alig.Cleaning-> blockSize == 40 );

        AND_WHEN ( "Getting Block Size" ) {
            REQUIRE ( alig.getBlockSize() == 40 );
        }
    }

    WHEN ( "Setting Keep Sequences Flag" ) {
        alig.setKeepSequencesFlag ( true );
        REQUIRE ( alig.Cleaning->keepSequences == true );
    }

    WHEN ( "Calculating Sequences Identity" ) {
        FAIL ( "Not implemented" );
        alig.calculateSeqIdentity();

        for ( int i = 0; i < alig.sequenNumber; i++ ) {
            CAPTURE ( std::vector<float> ( alig.identities[i], alig.identities[i] + alig.sequenNumber ) );
            CHECK ( false );
        }
    }

    WHEN ( "Calculating Sequences Overlap" ) {
        FAIL ( "Not implemented" );
        alig.calculateSeqOverlap();

        for ( int i = 0; i < alig.sequenNumber; i++ ) {
            CAPTURE ( std::vector<float> ( alig.overlaps[i], alig.overlaps[i] + alig.sequenNumber ) );
            CHECK ( false );
        }
    }

    WHEN ( "Getting Sequences" ) {
        string * names = new string[alig.sequenNumber];
        string * expectedNames = new string[alig.sequenNumber] {
            "Sequence 0",
            "Sequence 1",
            "Sequence 2",
            "Sequence 3",
            "Sequence 4",
            "Sequence 5"
        };

        GIVEN ( "Out string array" ) {
            alig.getSequences ( names );
            CAPTURE ( std::vector<string> ( names, names + alig.sequenNumber ) );
            CAPTURE ( std::vector<string> ( expectedNames, expectedNames + alig.sequenNumber ) );
            REQUIRE_THAT ( names,
                           ArrayContentsEqual ( expectedNames, alig.sequenNumber ) );
        }

        GIVEN ( "Out string array and out lengths" ) {
            int * sizePerSequence = new int[alig.sequenNumber],
            * expectedSizePerSequence = new int[alig.sequenNumber] { 40, 35, 40, 36, 38, 44 };
            alig.getSequences ( names, sizePerSequence );

            // Isolated blocks are necessary to link CAPTURE events with REQUIRE_THAT events, otherwise CAPTURE events get stacked.
            {
                CAPTURE ( std::vector<string> ( names, names + alig.sequenNumber ) );
                CAPTURE ( std::vector<string> ( expectedNames, expectedNames + alig.sequenNumber ) );
                REQUIRE_THAT ( names,
                               ArrayContentsEqual ( expectedNames, alig.sequenNumber ) );
            }

            {
                CAPTURE ( std::vector<int> ( sizePerSequence, sizePerSequence + alig.sequenNumber ) );
                CAPTURE ( std::vector<int> ( expectedSizePerSequence, expectedSizePerSequence + alig.sequenNumber ) );
                REQUIRE_THAT ( sizePerSequence,
                               ArrayContentsEqual ( expectedSizePerSequence, alig.sequenNumber ) );
            }

            delete [] sizePerSequence;
            delete [] expectedSizePerSequence;
        }

        GIVEN ( "Out string array, out sequences, out lengths" ) {
            int * sizePerSequence = new int[alig.sequenNumber],
            * expectedSizePerSequence = new int[alig.sequenNumber] { 40, 35, 40, 36, 38, 44 };
            string * sequences = new string[alig.sequenNumber],
            * expectedSequences = new string[alig.sequenNumber] {
                "FAYTAPDLLLIGFLLKTVATFGDTWFQLWQGLDLNKMPVF",
                "DPAVLFVIMLGTITKFSSEWFFAWLGLEINMMVII",
                "GLGKVIVYGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII",
                "PTILNIAGLHMETDINFSLAWFQAWGGLEINKQAIL",
                "ASGAILTLGIYLFTLCAVISVSWYLAWLGLEINMMAII",
                "AAAAAAAAALLTYLGLFLGTDYENFAAAAANAWLGLEINMMAQI"
            };
            alig.getSequences ( names, sequences, sizePerSequence );

            // Isolated blocks are necessary to link CAPTURE events with REQUIRE_THAT events, otherwise CAPTURE events get stacked.
            {
                CAPTURE ( std::vector<string> ( names, names + alig.sequenNumber ) );
                CAPTURE ( std::vector<string> ( expectedNames, expectedNames + alig.sequenNumber ) );
                REQUIRE_THAT ( names,
                               ArrayContentsEqual ( expectedNames, alig.sequenNumber ) );
            }

            {
                CAPTURE ( std::vector<int> ( sizePerSequence, sizePerSequence + alig.sequenNumber ) );
                CAPTURE ( std::vector<int> ( expectedSizePerSequence, expectedSizePerSequence + alig.sequenNumber ) );
                REQUIRE_THAT ( sizePerSequence,
                               ArrayContentsEqual ( expectedSizePerSequence, alig.sequenNumber ) );
            }

            {
                CAPTURE ( std::vector<string> ( sequences, sequences + alig.sequenNumber ) );
                CAPTURE ( std::vector<string> ( expectedSequences, expectedSequences + alig.sequenNumber ) );
                REQUIRE_THAT ( sequences,
                               ArrayContentsEqual ( expectedSequences, alig.sequenNumber ) );
            }
            delete [] sizePerSequence;
            delete [] expectedSizePerSequence;

            delete [] sequences;
            delete [] expectedSequences;
        }

        delete [] names;
        delete [] expectedNames;
    }

    WHEN ( "Getting sequences order by name" ) {
        string * Sequences = new string[alig.sequenNumber] {
            "Sequence 0",
            "Sequence 1",
            "Sequence 2",
            "Sequence 3",
            "Sequence 4",
            "Sequence 5"
        };

        int * expectedSequenceOrder = new int[alig.sequenNumber] { 0, 1, 2, 3, 4, 5 },
        * outSequenOrder = new int[alig.sequenNumber];


        GIVEN ( "Transformation 1" ) {
            REQUIRE ( alig.getSequenceNameOrder ( Sequences, outSequenOrder ) );

            CAPTURE ( std::vector<int> ( outSequenOrder, outSequenOrder + alig.sequenNumber ) );
            CAPTURE ( std::vector<int> ( expectedSequenceOrder, expectedSequenceOrder + alig.sequenNumber ) );
            REQUIRE_THAT ( expectedSequenceOrder, ArrayContentsEqual ( outSequenOrder, alig.sequenNumber ) );
        }

        GIVEN ( "Transformation 2" ) {
            std::swap ( Sequences[0], Sequences[3] );
            std::swap ( expectedSequenceOrder[0], expectedSequenceOrder[3] );

            REQUIRE ( alig.getSequenceNameOrder ( Sequences, outSequenOrder ) );

            CAPTURE ( std::vector<int> ( outSequenOrder, outSequenOrder + alig.sequenNumber ) );
            CAPTURE ( std::vector<int> ( expectedSequenceOrder, expectedSequenceOrder + alig.sequenNumber ) );
            REQUIRE_THAT ( expectedSequenceOrder, ArrayContentsEqual ( outSequenOrder, alig.sequenNumber ) );
        }

        GIVEN ( "Transformation 3" ) {
            std::swap ( Sequences[4], Sequences[1] );
            std::swap ( expectedSequenceOrder[4], expectedSequenceOrder[1] );

            REQUIRE ( alig.getSequenceNameOrder ( Sequences, outSequenOrder ) );

            CAPTURE ( std::vector<int> ( outSequenOrder, outSequenOrder + alig.sequenNumber ) );
            CAPTURE ( std::vector<int> ( expectedSequenceOrder, expectedSequenceOrder + alig.sequenNumber ) );
            REQUIRE_THAT ( expectedSequenceOrder, ArrayContentsEqual ( outSequenOrder, alig.sequenNumber ) );
        }

        GIVEN ( "More Sequences in test than alignment" ) {
            delete [] Sequences;
            Sequences = new string[alig.sequenNumber + 1] {
                "Sequence 0",
                "Sequence 1",
                "Sequence 2",
                "Sequence 3",
                "Sequence 4",
                "Sequence 5",
                "Sequence 6",
            };
            // This is a difficult assertion as we don't know the size of the initial array.
            // To check this, we should pass the names array size or use vectors as they contain the sizes.
//             REQUIRE_FALSE(alig.getSequenceNameOrder(Sequences, outSequenOrder));
        }

        GIVEN ( "Less Sequences in test than alignment" ) {
            delete [] Sequences;
            Sequences = new string[alig.sequenNumber - 1] {
                "Sequence 0",
                "Sequence 1",
                "Sequence 2",
                "Sequence 3",
                "Sequence 4"
            };

            REQUIRE_FALSE ( alig.getSequenceNameOrder ( Sequences, outSequenOrder ) );
        }

        GIVEN ( "Sequence not present" ) {
            Sequences[2] = "Troll Sequence";

            REQUIRE_FALSE ( alig.getSequenceNameOrder ( Sequences, outSequenOrder ) );
        }


        delete [] Sequences;
        delete [] outSequenOrder;
        delete [] expectedSequenceOrder;
    }

    WHEN ( "Getting alignment type, the return value is the same as utils::checkAlignmentType" )
    {
        REQUIRE ( alig.getAlignmentType() == utils::checkAlignmentType ( alig.sequenNumber, alig.residNumber, alig.sequences ) );
    }
    
    WHEN("Checking if aligned")
    {
        REQUIRE(alig.isAligned);
    }
    
    WHEN("Preparing Coding Sequence in a AA alignment")
    {
        REQUIRE_FALSE(alig.prepareCodingSequence(true, true, NULL));
    }
}



