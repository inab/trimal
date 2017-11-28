#include "../../include/Catch/catch.hpp"
#include "../../include/Catch/picojson.h"

#include "../../include/trimalArgumentParser.h"
#include "../../source/Catch/Matchers/ArrayMatcher.cpp"

#include <fstream>

#ifndef ReportMissingFiles
#define ReportMissingFiles true
#endif

void populate ( int *& s, picojson::value value ) {
    picojson::array pArray = value.get<picojson::array>();
    s = new int[pArray.size()];
    int i = 0;
    for ( auto it = pArray.begin(); it != pArray.end(); it++, i++ ) {
        s[i] = it->get<double>();
    }
}

void populate ( float *& s, picojson::value value ) {
    picojson::array pArray = value.get<picojson::array>();
    s = new float[pArray.size()];
    int i = 0;
    for ( auto it = pArray.begin(); it != pArray.end(); it++, i++ ) {
        s[i] = it->get<double>();
    }
}

void populate ( double *& s, picojson::value value ) {
    picojson::array pArray = value.get<picojson::array>();
    s = new double[pArray.size()];
    int i = 0;
    for ( auto it = pArray.begin(); it != pArray.end(); it++, i++ ) {
        s[i] = it->get<double>();
    }
}

void populate ( string *& s, picojson::value value ) {
    picojson::array pArray = value.get<picojson::array>();
    s = new string[pArray.size()];
    int i = 0;
    for ( auto it = pArray.begin(); it != pArray.end(); it++, i++ ) {
        s[i] = it->to_str();
    }
}

bool loadJSON ( const std::string & filename, picojson::value & alignmentTestData ) {
    ifstream testData;
    testData.open ( "./dataset/testingFiles/alignmentIntermediates/" + filename + ".json" );
    if ( testData.is_open() ) {
        std::string content ( ( std::istreambuf_iterator<char> ( testData ) ),
                              ( std::istreambuf_iterator<char>() ) );

        testData.close();
        // parse the input
        picojson::parse ( alignmentTestData, content );
        std::string err = picojson::get_last_error();
        if ( ! err.empty() ) {
            WARN ( err );
            return false;
        }

        // check if the type of the value is "object"
        if ( ! alignmentTestData.is<picojson::object>() ) {
            WARN ( "JSON is not an object" );
            return false;
        }
        return true;
    }
    return false;
}

SCENARIO ( "Alignment methods work correctly", "[alignment][aligMethods]" ) {

    std::vector<string> testingFilesVector = {
        "example.001.AA.clw",
        "example.001.AA.msl",
        "example.001.AA.phy",
        "example.002.AA.clw",
        "example.002.AA.phy",
        "example.003.AA.clw",
        "example.004.AA.fasta",
        "example.005.AA.fasta",
        "example.006.AA.pir",
        "example.007.AA.fasta",
        "example.007.AA.only_seqs",
        "example.009.AA.fasta",
        "example.010.AA.fasta",
        "example.011.AA.YKL197C.clw",
        "example.011.AA.YKL197C.fasta",
        "example.011.AA.YKL197C.phy",
        "example.012.AA.SuperAlignment.phy",
        "example.013.AA.SuperAlignment.phy",
        "example.014.AA.EggNOG.COG0591.fasta",
        "example.015.AA.bctoNOG.ENOG41099F3.fasta",
        "example.016.AA.bctoNOG.ENOG41099FB.fasta",
        "example.017.AA.bctoNOG.ENOG41099FJ.fasta",
        "example.018.AA.bctoNOG.ENOG41099FV.fasta",
        "example.019.AA.bctoNOG.ENOG41099HI.fasta",
        "example.020.AA.bctoNOG.ENOG41099HN.fasta",
        "example.021.AA.bctoNOG.ENOG41099I5.fasta",
        "example.022.AA.bctoNOG.ENOG41099IZ.fasta",
        "example.023.AA.bctoNOG.ENOG41099K3.fasta",
        "example.024.AA.bctoNOG.ENOG41099KM.fasta",
        "example.025.AA.bctoNOG.ENOG41099KP.fasta",
        "example.026.AA.bctoNOG.ENOG41099MV.fasta",
        "example.027.AA.bctoNOG.ENOG41099NY.fasta",
        "example.028.AA.bctoNOG.ENOG41099PA.fasta",
        "example.029.AA.bctoNOG.ENOG41099Q3.fasta",
        "example.030.AA.bctoNOG.ENOG41099RG.fasta",
        "example.031.AA.bctoNOG.ENOG41099UK.fasta",
        "example.032.AA.bctoNOG.ENOG41099UW.fasta",
        "example.033.AA.bctoNOG.ENOG41099VK.fasta",
        "example.034.AA.bctoNOG.ENOG41099WA.fasta",
        "example.035.AA.bctoNOG.ENOG41099WF.fasta",
        "example.036.AA.bctoNOG.ENOG41099XJ.fasta",
        "example.037.AA.bctoNOG.ENOG41099XP.fasta",
        "example.038.AA.bctoNOG.ENOG41099Y4.fasta",
        "example.039.AA.bctoNOG.ENOG41099YD.fasta",
        "example.040.AA.bctoNOG.ENOG4109A32.fasta",
        "example.041.AA.bctoNOG.ENOG4109A5T.fasta",
        "example.042.AA.bctoNOG.ENOG4109A9M.fasta",
        "example.043.AA.bctoNOG.ENOG4109ADN.fasta",
        "example.044.AA.bctoNOG.ENOG4109AED.fasta",
        "example.045.AA.bctoNOG.ENOG4109AGT.fasta",
        "example.046.AA.bctoNOG.ENOG4109AGW.fasta",
        "example.047.AA.bctoNOG.ENOG4109AIC.fasta",
        "example.048.AA.bctoNOG.ENOG4109AJ3.fasta",
        "example.049.AA.bctoNOG.ENOG4109AY5.fasta",
        "example.050.AA.bctoNOG.ENOG4109B8Z.fasta",
        "example.051.AA.bctoNOG.ENOG4109BCJ.fasta",
        "example.052.AA.bctoNOG.ENOG4109CTU.fasta",
        "example.053.AA.bctoNOG.ENOG4109CVC.fasta",
        "example.054.AA.bctoNOG.ENOG4109FIT.fasta",
        "example.055.AA.bctoNOG.ENOG4109GY9.fasta",
        "example.056.AA.bctoNOG.ENOG4109IPJ.fasta",
        "example.057.AA.bctoNOG.ENOG4109SZ2.fasta",
        "example.058.AA.strNOG.ENOG411BBR6.fasta",
        "example.059.AA.strNOG.ENOG411BBRR.fasta",
        "example.060.AA.strNOG.ENOG411BBWK.fasta",
        "example.061.AA.strNOG.ENOG411BCDZ.fasta",
        "example.062.AA.strNOG.ENOG411BCX3.fasta",
        "example.063.AA.strNOG.ENOG411BDBU.fasta",
        "example.064.AA.strNOG.ENOG411BDKC.fasta",
        "example.065.AA.strNOG.ENOG411BDSZ.fasta",
        "example.066.AA.strNOG.ENOG411BDUE.fasta",
        "example.067.AA.strNOG.ENOG411BDX3.fasta",
        "example.068.AA.strNOG.ENOG411BE45.fasta",
        "example.069.AA.strNOG.ENOG411BE8B.fasta",
        "example.070.AA.strNOG.ENOG411BEUV.fasta",
        "example.071.AA.strNOG.ENOG411BEZ0.fasta",
        "example.072.AA.strNOG.ENOG411BF1S.fasta",
        "example.073.AA.strNOG.ENOG411BFCW.fasta",
        "example.074.AA.strNOG.ENOG411BFPF.fasta",
        "example.075.AA.strNOG.ENOG411BFQS.fasta",
        "example.075.AA.strNOG.ENOG411BFQS.nxs",
        "example.076.AA.strNOG.ENOG411BH75.fasta",
        "example.077.AA.strNOG.ENOG411BH79.fasta",
        "example.078.AA.strNOG.ENOG411BH99.fasta",
        "example.079.AA.strNOG.ENOG411BJDC.fasta",
        "example.080.AA.strNOG.ENOG411BJIF.fasta",
        "example.081.AA.strNOG.ENOG411BK9X.fasta",
        "example.082.AA.strNOG.ENOG411BKC5.fasta",
        "example.083.AA.strNOG.ENOG411BMKC.fasta",
        "example.084.AA.strNOG.ENOG411BNP9.fasta",
        "example.085.AA.strNOG.ENOG411BQTJ.fasta",
        "example.086.AA.strNOG.ENOG411BR1D.fasta",
        "example.087.AA.strNOG.ENOG411BRCH.fasta",
        "example.088.AA.strNOG.ENOG411BSXF.fasta",
        "example.089.AA.strNOG.ENOG411BV9B.fasta",
        "example.090.AA.strNOG.ENOG411BVKR.fasta",
        "example.091.AA.strNOG.ENOG411BWBU.codon.fa",
        "example.091.AA.strNOG.ENOG411BWBU.fasta",
        "example.092.DNA.fasta",
        "example.093.DNA.fasta",
        "example.094.DNADeg.sequential_phy",
    };

    for ( string & filename : testingFilesVector ) {

        picojson::value alignmentTestData;
        if ( !loadJSON ( filename, alignmentTestData ) ) {
#if ReportMissingFiles
            WARN ( "Failed to load JSON containing tests results\nSkipping all tests for: " + filename );
#endif
            continue;
        }

        if ( alignmentTestData.contains ( "filename" ) ) {
            GIVEN ( alignmentTestData.get ( "filename" ).to_str() ) {
                newAlignment alig;
                int i = 0;
                {
                    INFO ( "Impossible to test current alignment as it doesn't contain all the information required to make the testing" );
                    REQUIRE ( alignmentTestData.contains ( "sequenNumber" ) );
                    REQUIRE ( alignmentTestData.contains ( "residNumber" ) );
                    REQUIRE ( alignmentTestData.contains ( "sequences" ) );
                    REQUIRE ( alignmentTestData.contains ( "aligned" ) );
                    REQUIRE ( alignmentTestData.contains ( "names" ) );
                }

                alig.sequenNumber = alignmentTestData.get ( "sequenNumber" ).get<double>();
                alig.originalSequenNumber = alig.sequenNumber;
                alig.residNumber = alignmentTestData.get ( "residNumber" ).get<double>();
                alig.originalResidNumber = alig.residNumber;
                alig.isAligned = alignmentTestData.get ( "aligned" ).get<bool>();

                populate ( alig.sequences, alignmentTestData.get ( "sequences" ) );
                populate ( alig.seqsName, alignmentTestData.get ( "names" ) );

                alig.fillMatrices ( alig.isAligned );

                int * saveResidues = new int[alig.residNumber];
                for ( i = 0; i < alig.residNumber; i++ ) {
                    saveResidues[i] = i;
                }

                {
                    REQUIRE ( saveResidues != NULL );
                    REQUIRE ( alig.saveResidues != NULL );

                    INFO ( "Save residues array in alignment doesn't match with expected saveResidues" );
                    auto expected = std::vector<int> ( saveResidues, saveResidues + alig.residNumber );
                    CAPTURE ( expected );
                    auto obtained = std::vector<int> ( alig.saveResidues, alig.saveResidues + alig.residNumber );
                    CAPTURE ( obtained );
                    REQUIRE_THAT ( alig.saveResidues, ArrayContentsEqual ( saveResidues, alig.residNumber ) );
                }

                WHEN ( "Testing Alignment::getTranslationCDS" ) {
                    WARN ( "Not implemented" );
                }

                WHEN ( "Getting number of species" ) {
                    REQUIRE ( alig.getNumSpecies() == alignmentTestData.get ( "sequenNumber" ).get<double>() );
                }

                WHEN ( "Getting number of residues" ) {
                    REQUIRE ( alig.getNumAminos() == alignmentTestData.get ( "residNumber" ).get<double>() );
                }

                WHEN ( "Getting Sequences" ) {
                    string * names = new string[alig.sequenNumber];
                    string * expectedNames;

                    populate ( expectedNames, alignmentTestData.get ( "names" ) );

                    if ( alignmentTestData.contains ( "names" ) ) {
                        GIVEN ( "Out string array" ) {
                            alig.getSequences ( names );
                            CAPTURE ( std::vector<string> ( names, names + alig.sequenNumber ) );
                            CAPTURE ( std::vector<string> ( expectedNames, expectedNames + alig.sequenNumber ) );
                            REQUIRE_THAT ( names,
                                           ArrayContentsEqual ( expectedNames, alig.sequenNumber ) );
                        }
                    } else WARN ( "Alignment file does not contain 'names' variable.\nSkipping test 'Out string array'" );

                    if ( alignmentTestData.contains ( "names" ) ) {
                        if ( alignmentTestData.contains ( "noGapsSize" ) ) {
                            GIVEN ( "Out string array and out lengths" ) {
                                int * sizePerSequence = new int[alig.sequenNumber], * expectedSizePerSequence;
                                populate ( expectedSizePerSequence, alignmentTestData.get ( "noGapsSize" ) );

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
                        } else WARN ( "Alignment file does not contain 'noGapsSize' variable.\nSkipping test 'Out string array and out lengths'" );
                    } else WARN ( "Alignment file does not contain 'names' variable.\nSkipping test 'Out string array and out lengths'" );

                    if ( alignmentTestData.contains ( "names" ) ) {
                        if ( alignmentTestData.contains ( "noGapsSize" ) ) {
                            if ( alignmentTestData.contains ( "sequences" ) ) {
                                GIVEN ( "Out string array, out sequences, out lengths" ) {
                                    int * sizePerSequence = new int[alig.sequenNumber], * expectedSizePerSequence;
                                    populate ( expectedSizePerSequence, alignmentTestData.get ( "noGapsSize" ) );

                                    string * sequences = new string[alig.sequenNumber], * expectedSequences;
                                    populate ( expectedSequences, alignmentTestData.get ( "sequences" ) );

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
                            } else WARN ( "Alignment file does not contain 'sequences' variable.\nSkipping test 'Out string array, out sequences, out lengths'" );
                        } else WARN ( "Alignment file does not contain 'noGapsSize' variable.\nSkipping test 'Out string array, out sequences, out lengths'" );
                    } else WARN ( "Alignment file does not contain 'names' variable.\nSkipping test 'Out string array, out sequences, out lengths'" );
                }
            }
        }

//         else WARN ( "file " << filename << ".json not found.\nSkipping all test for alignment" );
    }
    /*



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

        WHEN ( "Getting alignment type, the return value is the same as utils::checkAlignmentType" ) {
            REQUIRE ( alig.getAlignmentType() == utils::checkAlignmentType ( alig.sequenNumber, alig.residNumber, alig.sequences ) );
        }

        WHEN ( "Checking if aligned" ) {
            REQUIRE ( alig.isAligned );
        }

        WHEN ( "Preparing Coding Sequence in a AA alignment" ) {
            REQUIRE_FALSE ( alig.prepareCodingSequence ( true, true, NULL ) );
        }*/
}





