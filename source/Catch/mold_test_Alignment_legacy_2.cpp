#ifndef mold_test_includes
#define mold_test_includes
#include "Catch/catch.hpp"
#include "Catch/picojson.h"
#include "Catch/catchhelperfunctions.h"

#include "trimalArgumentParser.h"
#include "Matchers/ArrayMatcher.cpp"

#include <fstream>

#endif

static picojson::value testData;
// EndOfHeader // <- Do not remove this comment, as is used by TestMaker.


SCENARIO ("Alignment methods work correctly "
                  "in input_filename",
          "[alignment][aligMethods]") {


        if (!loadJSON("input_filename", testData)) {

            WARN ( "Failed to load JSON containing tests results"
                           "\nSkipping all tests for: input_filename");
        } else
        {
            if (
                    testData.contains("sequenNumber") &&
                    testData.contains("residNumber") &&
                    testData.contains("sequences") &&
                    testData.contains("filename") &&
                    testData.contains("aligned") &&
                    testData.contains("names")
                    ) {

                int i;
                newAlignment alig;

                alig.sequenNumber = static_cast<int>(testData.get("sequenNumber").get<double>());
                alig.originalSequenNumber = alig.sequenNumber;
                alig.residNumber = static_cast<int>(testData.get("residNumber").get<double>());
                alig.originalResidNumber = alig.residNumber;
                alig.isAligned = testData.get("aligned").get<bool>();

                populate(alig.sequences, testData.get("sequences"));
                populate(alig.seqsName, testData.get("names"));

                alig.fillMatrices(alig.isAligned);

                if (testData.contains("aligned") && testData.get("aligned").get<bool>()) {
                    alig.sgaps = new statisticsGaps(&alig);
                    alig.scons = new statisticsConservation2(&alig);

                    int *saveResidues = new int[alig.residNumber];
                    for (i = 0; i < alig.residNumber; i++) {
                        saveResidues[i] = i;
                    }

                    WHEN ("Checking save residues is initialized correctly") {
                        REQUIRE (saveResidues != NULL);
                        REQUIRE (alig.saveResidues != NULL);

                        INFO ("Save residues array in alignment doesn't match with expected saveResidues");
                        auto expected = std::vector<int>(saveResidues, saveResidues + alig.residNumber);
                        CAPTURE (expected);
                        auto obtained = std::vector<int>(alig.saveResidues, alig.saveResidues + alig.residNumber);
                        CAPTURE (obtained);
                        REQUIRE_THAT (alig.saveResidues, ArrayContentsEqual(saveResidues, alig.residNumber));

                        delete[] saveResidues;
                        saveResidues = NULL;
                    }

                    if (saveResidues != NULL)
                        delete[] saveResidues;
                    saveResidues = NULL;
                }

                WHEN ("Getting number of species") {
                    REQUIRE (alig.getNumSpecies() == testData.get("sequenNumber").get<double>());
                }

                WHEN ("Getting number of residues") {
                    REQUIRE (alig.getNumAminos() == testData.get("residNumber").get<double>());
                }

                WHEN ("Getting Sequences") {
                    string *names = new string[alig.sequenNumber];
                    string *expectedNames;

                    populate(expectedNames, testData.get("names"));

                    GIVEN ("Out string array") {
                        alig.getSequences(names);
                        CAPTURE (std::vector<string>(names, names + alig.sequenNumber));
                        CAPTURE (std::vector<string>(expectedNames, expectedNames + alig.sequenNumber));
                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber));
                    }

                    GIVEN ("Out string array and out lengths") {
                        if (testData.contains("names")) {
                            if (testData.contains("noGapsSize")) {
                                int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                                populate(expectedSizePerSequence, testData.get("noGapsSize"));

                                alig.getSequences(names, sizePerSequence);

                                // Isolated blocks are necessary to link CAPTURE events with REQUIRE_THAT events,
                                // otherwise CAPTURE events get stacked.
                                {
                                    CAPTURE (std::vector<string>(names, names + alig.sequenNumber));
                                    CAPTURE (std::vector<string>(expectedNames, expectedNames + alig.sequenNumber));
                                    REQUIRE_THAT (names,
                                                  ArrayContentsEqual(expectedNames, alig.sequenNumber));
                                }

                                {
                                    CAPTURE (std::vector<int>(sizePerSequence, sizePerSequence + alig.sequenNumber));
                                    CAPTURE (std::vector<int>(expectedSizePerSequence,
                                                              expectedSizePerSequence + alig.sequenNumber));
                                    REQUIRE_THAT (sizePerSequence,
                                                  ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber));
                                }

                                delete[] sizePerSequence;
                                delete[] expectedSizePerSequence;
                            } else
                                WARN ("Alignment testfile does not contain 'noGapsSize' variable."
                                              "\nSkipping test 'Out string array and out lengths'");
                        } else
                            WARN ("Alignment testfile does not contain 'names' variable."
                                          "\nSkipping test 'Out string array and out lengths'");
                    }

                    GIVEN ("Out string array, out sequences, out lengths") {
                        if (testData.contains("names")) {
                            if (testData.contains("noGapsSize")) {
                                if (testData.contains("noGapsSequences")) {
                                    int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                                    populate(expectedSizePerSequence, testData.get("noGapsSize"));

                                    string *sequences = new string[alig.sequenNumber], *expectedSequences;
                                    populate(expectedSequences, testData.get("noGapsSequences"));

                                    alig.getSequences(names, sequences, sizePerSequence);

                                    // Isolated blocks are necessary to link CAPTURE events with REQUIRE_THAT events,
                                    // otherwise CAPTURE events get stacked.
                                    {
                                        CAPTURE (std::vector<string>(names,
                                                                     names + alig.sequenNumber));
                                        CAPTURE (std::vector<string>(expectedNames,
                                                                     expectedNames + alig.sequenNumber));
                                        REQUIRE_THAT (names,
                                                      ArrayContentsEqual(expectedNames, alig.sequenNumber));
                                    }

                                    {
                                        CAPTURE (std::vector<int>(sizePerSequence,
                                                                  sizePerSequence + alig.sequenNumber));
                                        CAPTURE (std::vector<int>(expectedSizePerSequence,
                                                                  expectedSizePerSequence + alig.sequenNumber));
                                        REQUIRE_THAT (sizePerSequence,
                                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber));
                                    }

                                    {
                                        CAPTURE (std::vector<string>(sequences, sequences + alig.sequenNumber));
                                        CAPTURE (std::vector<string>(expectedSequences,
                                                                     expectedSequences + alig.sequenNumber));
                                        REQUIRE_THAT (sequences,
                                                      ArrayContentsEqual(expectedSequences, alig.sequenNumber));
                                    }
                                    delete[] sizePerSequence;
                                    delete[] expectedSizePerSequence;

                                    delete[] sequences;
                                    delete[] expectedSequences;


                                } else
                                    WARN ("Alignment testfile does not contain 'sequences' variable."
                                                  "\nSkipping test 'Out string array, out sequences, out lengths'");
                            } else
                                WARN ("Alignment testfile does not contain 'noGapsSize' variable."
                                              "\nSkipping test 'Out string array, out sequences, out lengths'");
                        } else
                            WARN ("Alignment testfile does not contain 'names' variable."
                                          "\nSkipping test 'Out string array, out sequences, out lengths'");
                    }
                    delete[] names;
                    delete[] expectedNames;
                }

                WHEN ("Calculate Gaps") {
                    if (testData.contains("gapsPerColumn")) {

                        if (testData.contains("aligned")) {
                            if (testData.get("aligned").get<bool>()) {
                                int *GapsPerColumn = alig.sgaps->getGapsWindow(), *expectedGapsPerColumn;
                                populate(expectedGapsPerColumn, testData.get("gapsPerColumn"));

                                REQUIRE_THAT (GapsPerColumn, ArrayContentsEqual(expectedGapsPerColumn, alig.residNumber));

                                delete[] expectedGapsPerColumn;
                            } //else WARN ( "Alignment testfile is not aligned.\nSkipping test 'Calculate Gaps'" );
                        } else
                            WARN ("Alignment testfile does not contain 'aligned' variable\nSkipping test 'Calculate Gaps'");
                    } else
                        WARN ("Alignment testfile does not contain 'gapsPerColumn' variable\nSkipping test 'Calculate Gaps'");
                }


            } // end GIVEN ( testData.get ( "filename" ).to_str() )
        } // end if ( testData.contains ( all ) )
//    }


//    } // end for ( string & filename : testingFilesVector )
}