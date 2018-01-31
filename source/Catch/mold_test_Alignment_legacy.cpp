#include "Catch/catch.hpp"
#include "Catch/picojson.h"
#include "Catch/catchhelperfunctions.h"

#include "trimalArgumentParser.h"
#include "Matchers/ArrayMatcher.cpp"

#include <fstream>

#ifndef ReportMissingFiles
#define ReportMissingFiles false
#endif



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
//         "example.011.AA.YKL197C.clw",
//         "example.011.AA.YKL197C.fasta",
//         "example.011.AA.YKL197C.phy",
//         "example.012.AA.SuperAlignment.phy",
//         "example.013.AA.SuperAlignment.phy",
//         "example.014.AA.EggNOG.COG0591.fasta",
//         "example.015.AA.bctoNOG.ENOG41099F3.fasta",
//         "example.016.AA.bctoNOG.ENOG41099FB.fasta",
//         "example.017.AA.bctoNOG.ENOG41099FJ.fasta",
//         "example.018.AA.bctoNOG.ENOG41099FV.fasta",
//         "example.019.AA.bctoNOG.ENOG41099HI.fasta",
//         "example.020.AA.bctoNOG.ENOG41099HN.fasta",
//         "example.021.AA.bctoNOG.ENOG41099I5.fasta",
//         "example.022.AA.bctoNOG.ENOG41099IZ.fasta",
//         "example.023.AA.bctoNOG.ENOG41099K3.fasta",
//         "example.024.AA.bctoNOG.ENOG41099KM.fasta",
//         "example.025.AA.bctoNOG.ENOG41099KP.fasta",
//         "example.026.AA.bctoNOG.ENOG41099MV.fasta",
//         "example.027.AA.bctoNOG.ENOG41099NY.fasta",
//         "example.028.AA.bctoNOG.ENOG41099PA.fasta",
//         "example.029.AA.bctoNOG.ENOG41099Q3.fasta",
//         "example.030.AA.bctoNOG.ENOG41099RG.fasta",
//         "example.031.AA.bctoNOG.ENOG41099UK.fasta",
//         "example.032.AA.bctoNOG.ENOG41099UW.fasta",
//         "example.033.AA.bctoNOG.ENOG41099VK.fasta",
//         "example.034.AA.bctoNOG.ENOG41099WA.fasta",
//         "example.035.AA.bctoNOG.ENOG41099WF.fasta",
//         "example.036.AA.bctoNOG.ENOG41099XJ.fasta",
//         "example.037.AA.bctoNOG.ENOG41099XP.fasta",
//         "example.038.AA.bctoNOG.ENOG41099Y4.fasta",
//         "example.039.AA.bctoNOG.ENOG41099YD.fasta",
//         "example.040.AA.bctoNOG.ENOG4109A32.fasta",
//         "example.041.AA.bctoNOG.ENOG4109A5T.fasta",
//         "example.042.AA.bctoNOG.ENOG4109A9M.fasta",
//         "example.043.AA.bctoNOG.ENOG4109ADN.fasta",
//         "example.044.AA.bctoNOG.ENOG4109AED.fasta",
//         "example.045.AA.bctoNOG.ENOG4109AGT.fasta",
//         "example.046.AA.bctoNOG.ENOG4109AGW.fasta",
//         "example.047.AA.bctoNOG.ENOG4109AIC.fasta",
//         "example.048.AA.bctoNOG.ENOG4109AJ3.fasta",
//         "example.049.AA.bctoNOG.ENOG4109AY5.fasta",
//         "example.050.AA.bctoNOG.ENOG4109B8Z.fasta",
//         "example.051.AA.bctoNOG.ENOG4109BCJ.fasta",
//         "example.052.AA.bctoNOG.ENOG4109CTU.fasta",
//         "example.053.AA.bctoNOG.ENOG4109CVC.fasta",
//         "example.054.AA.bctoNOG.ENOG4109FIT.fasta",
//         "example.055.AA.bctoNOG.ENOG4109GY9.fasta",
//         "example.056.AA.bctoNOG.ENOG4109IPJ.fasta",
//         "example.057.AA.bctoNOG.ENOG4109SZ2.fasta",
//         "example.058.AA.strNOG.ENOG411BBR6.fasta",
//         "example.059.AA.strNOG.ENOG411BBRR.fasta",
//         "example.060.AA.strNOG.ENOG411BBWK.fasta",
//         "example.061.AA.strNOG.ENOG411BCDZ.fasta",
//         "example.062.AA.strNOG.ENOG411BCX3.fasta",
//         "example.063.AA.strNOG.ENOG411BDBU.fasta",
//         "example.064.AA.strNOG.ENOG411BDKC.fasta",
//         "example.065.AA.strNOG.ENOG411BDSZ.fasta",
//         "example.066.AA.strNOG.ENOG411BDUE.fasta",
//         "example.067.AA.strNOG.ENOG411BDX3.fasta",
//         "example.068.AA.strNOG.ENOG411BE45.fasta",
//         "example.069.AA.strNOG.ENOG411BE8B.fasta",
//         "example.070.AA.strNOG.ENOG411BEUV.fasta",
//         "example.071.AA.strNOG.ENOG411BEZ0.fasta",
//         "example.072.AA.strNOG.ENOG411BF1S.fasta",
//         "example.073.AA.strNOG.ENOG411BFCW.fasta",
//         "example.074.AA.strNOG.ENOG411BFPF.fasta",
//         "example.075.AA.strNOG.ENOG411BFQS.fasta",
//         "example.075.AA.strNOG.ENOG411BFQS.nxs",
//         "example.076.AA.strNOG.ENOG411BH75.fasta",
//         "example.077.AA.strNOG.ENOG411BH79.fasta",
//         "example.078.AA.strNOG.ENOG411BH99.fasta",
//         "example.079.AA.strNOG.ENOG411BJDC.fasta",
//         "example.080.AA.strNOG.ENOG411BJIF.fasta",
//         "example.081.AA.strNOG.ENOG411BK9X.fasta",
//         "example.082.AA.strNOG.ENOG411BKC5.fasta",
//         "example.083.AA.strNOG.ENOG411BMKC.fasta",
//         "example.084.AA.strNOG.ENOG411BNP9.fasta",
//         "example.085.AA.strNOG.ENOG411BQTJ.fasta",
//         "example.086.AA.strNOG.ENOG411BR1D.fasta",
//         "example.087.AA.strNOG.ENOG411BRCH.fasta",
//         "example.088.AA.strNOG.ENOG411BSXF.fasta",
//         "example.089.AA.strNOG.ENOG411BV9B.fasta",
//         "example.090.AA.strNOG.ENOG411BVKR.fasta",
//         "example.091.AA.strNOG.ENOG411BWBU.codon.fa",
//         "example.091.AA.strNOG.ENOG411BWBU.fasta",
//         "example.092.DNA.fasta",
//         "example.093.DNA.fasta",
//         "example.094.DNADeg.sequential_phy",
    };

    for ( string & filename : testingFilesVector ) {

        GIVEN ( filename ) {
            static picojson::value testData;
            if ( !loadJSON ( filename, testData ) ) {
#if ReportMissingFiles
                WARN ( "Failed to load JSON containing tests results\nSkipping all tests for: " + filename );
#endif
                continue;
            }

            if (
                testData.contains ( "sequenNumber" ) &&
                testData.contains ( "residNumber" ) &&
                testData.contains ( "sequences" ) &&
                testData.contains ( "filename" ) &&
                testData.contains ( "aligned" ) &&
                testData.contains ( "names" )
            ) {

                int i;
                newAlignment alig;

                alig.sequenNumber = testData.get ( "sequenNumber" ).get<double>();
                alig.originalSequenNumber = alig.sequenNumber;
                alig.residNumber = testData.get ( "residNumber" ).get<double>();
                alig.originalResidNumber = alig.residNumber;
                alig.isAligned = testData.get ( "aligned" ).get<bool>();

                populate ( alig.sequences, testData.get ( "sequences" ) );
                populate ( alig.seqsName, testData.get ( "names" ) );

                alig.fillMatrices ( alig.isAligned );

                if ( testData.contains ( "aligned" ) && testData.get ( "aligned" ).get<bool>() ) {
                    alig.sgaps = new statisticsGaps ( &alig );
                    alig.scons = new statisticsConservation2 ( &alig );

                    int * saveResidues = new int[alig.residNumber];
                    for ( i = 0; i < alig.residNumber; i++ ) {
                        saveResidues[i] = i;
                    }

                    WHEN ( "Checking save residues is initialized correctly" ) {
                        REQUIRE ( saveResidues != NULL );
                        REQUIRE ( alig.saveResidues != NULL );

                        INFO ( "Save residues array in alignment doesn't match with expected saveResidues" );
                        auto expected = std::vector<int> ( saveResidues, saveResidues + alig.residNumber );
                        CAPTURE ( expected );
                        auto obtained = std::vector<int> ( alig.saveResidues, alig.saveResidues + alig.residNumber );
                        CAPTURE ( obtained );
                        REQUIRE_THAT ( alig.saveResidues, ArrayContentsEqual ( saveResidues, alig.residNumber ) );

                        if ( saveResidues != NULL )
                            delete [] saveResidues;
                        saveResidues = NULL;
                    }

                    if ( saveResidues != NULL )
                        delete [] saveResidues;
                    saveResidues = NULL;
                }

// NOTE NOT IMPLEMENTED
//                 WHEN ( "Testing Alignment::getTranslationCDS" ) {
//                     WARN ( "Not implemented" );
//                 }

                WHEN ( "Getting number of species" ) {
                    REQUIRE ( alig.getNumSpecies() == testData.get ( "sequenNumber" ).get<double>() );
                }

                WHEN ( "Getting number of residues" ) {
                    REQUIRE ( alig.getNumAminos() == testData.get ( "residNumber" ).get<double>() );
                }

                WHEN ( "Getting Sequences" ) {
                    string * names = new string[alig.sequenNumber];
                    string * expectedNames;

                    populate ( expectedNames, testData.get ( "names" ) );

                    GIVEN ( "Out string array" ) {
                        alig.getSequences ( names );
                        CAPTURE ( std::vector<string> ( names, names + alig.sequenNumber ) );
                        CAPTURE ( std::vector<string> ( expectedNames, expectedNames + alig.sequenNumber ) );
                        REQUIRE_THAT ( names,
                                       ArrayContentsEqual ( expectedNames, alig.sequenNumber ) );
                    }

                    GIVEN ( "Out string array and out lengths" ) {
                        if ( testData.contains ( "names" ) ) {
                            if ( testData.contains ( "noGapsSize" ) ) {
                                int * sizePerSequence = new int[alig.sequenNumber], * expectedSizePerSequence;
                                populate ( expectedSizePerSequence, testData.get ( "noGapsSize" ) );

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
                            } else WARN ( "Alignment testfile does not contain 'noGapsSize' variable.\nSkipping test 'Out string array and out lengths'" );
                        } else WARN ( "Alignment testfile does not contain 'names' variable.\nSkipping test 'Out string array and out lengths'" );
                    }

                    GIVEN ( "Out string array, out sequences, out lengths" ) {
                        if ( testData.contains ( "names" ) ) {
                            if ( testData.contains ( "noGapsSize" ) ) {
                                if ( testData.contains ( "noGapsSequences" ) ) {
                                    int * sizePerSequence = new int[alig.sequenNumber], * expectedSizePerSequence;
                                    populate ( expectedSizePerSequence, testData.get ( "noGapsSize" ) );

                                    string * sequences = new string[alig.sequenNumber], * expectedSequences;
                                    populate ( expectedSequences, testData.get ( "noGapsSequences" ) );

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




                                } else WARN ( "Alignment testfile does not contain 'sequences' variable.\nSkipping test 'Out string array, out sequences, out lengths'" );
                            } else WARN ( "Alignment testfile does not contain 'noGapsSize' variable.\nSkipping test 'Out string array, out sequences, out lengths'" );
                        } else WARN ( "Alignment testfile does not contain 'names' variable.\nSkipping test 'Out string array, out sequences, out lengths'" );
                    }
                    delete [] names;
                    delete [] expectedNames;
                }

                WHEN ( "Calculate Gaps" ) {
                    if ( testData.contains ( "gapsPerColumn" ) ) {

                        if ( testData.contains ( "aligned" ) ) {
                            if ( testData.get ( "aligned" ).get<bool>() ) {
                                int * GapsPerColumn = alig.sgaps->getGapsWindow(), * expectedGapsPerColumn;
                                populate ( expectedGapsPerColumn, testData.get ( "gapsPerColumn" ) );

                                REQUIRE_THAT ( GapsPerColumn, ArrayContentsEqual ( expectedGapsPerColumn, alig.residNumber ) );

                                delete [] expectedGapsPerColumn;
                            } //else WARN ( "Alignment testfile is not aligned.\nSkipping test 'Calculate Gaps'" );
                        } else WARN ( "Alignment testfile does not contain 'aligned' variable\nSkipping test 'Calculate Gaps'" );
                    } else WARN ( "Alignment testfile does not contain 'gapsPerColumn' variable\nSkipping test 'Calculate Gaps'" );
                }





            } // end GIVEN ( testData.get ( "filename" ).to_str() )
        } // end if ( testData.contains ( all ) )
    } // end for ( string & filename : testingFilesVector )
}






