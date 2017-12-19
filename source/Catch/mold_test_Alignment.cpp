#include "Catch/catch.hpp"
#include "Catch/picojson.h"
#include "Catch/catchhelperfunctions.h"

#include "trimalArgumentParser.h"
#include "Matchers/ArrayMatcher.cpp"
#include <vector>

#include <fstream>

int i;

// tag input_filename will be replaced for each name in ./database
// EndOfHeader // <- Do not remove this comment, as is used by TestMaker.
SCENARIO ("Alignment methods work correctly in input_filename", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("input_filename", testData)) {

        alig.sequenNumber = testData.get("sequenNumber").get<double>();
        alig.originalSequenNumber = alig.sequenNumber;
        alig.residNumber = testData.get("residNumber").get<double>();
        alig.originalResidNumber = alig.residNumber;
        alig.isAligned = testData.get("aligned").get<bool>();

        populate(alig.sequences, testData.get("sequences"));
        populate(alig.seqsName, testData.get("names"));

        alig.fillMatrices(alig.isAligned);

        if (testData.contains("aligned") && testData.get("aligned").get<bool>()) {
            alig.sgaps = new statisticsGaps(&alig);
            alig.scons = new statisticsConservation2(&alig);
            // Init Similarity Matrix
            {
                int SequenceType = testData.get("SequenceType").get<double>();
                similarityMatrix * sm = new similarityMatrix();
                switch (SequenceType)
                {
                    case 1: case 2:
                        sm->defaultNTSimMatrix();
                        break;
                    case 3:
                        sm->defaultAASimMatrix();
                        break;
                    case 4: case 5:
                        sm->defaultNTDegeneratedSimMatrix();
                        break;
                    default:
                        sm->defaultNTSimMatrix();
                }
                alig.scons->setSimilarityMatrix(sm);
            }

        }

        WHEN ("Checking save residues is initialized correctly") {

            if (testData.contains("aligned") && testData.get("aligned").get<bool>()) {

                int *saveResidues = new int[alig.residNumber];
                for (i = 0; i < alig.residNumber; i++) {
                    saveResidues[i] = i;
                }

                REQUIRE (saveResidues != NULL);
                REQUIRE (alig.saveResidues != NULL);

                INFO ("Save residues array in alignment doesn't match with expected saveResidues");
                auto expected =
                        std::vector<int>(saveResidues,
                                         saveResidues + alig.residNumber);

                CAPTURE (expected);

                auto obtained =
                        std::vector<int>(alig.saveResidues,
                                         alig.saveResidues + alig.residNumber);

                CAPTURE (obtained);

                REQUIRE_THAT (alig.saveResidues,
                              ArrayContentsEqual(saveResidues, alig.residNumber,
                                      "./dataset/testingFiles/alignmentErrors/input_filename.saveResidues.error.tsv"));

                delete[] saveResidues;
            }
        }

        WHEN ("Getting number of species") {
            REQUIRE (alig.getNumSpecies() ==
                             testData.get("sequenNumber").get<double>());
        }

        WHEN ("Getting number of residues") {
            REQUIRE (alig.getNumAminos() ==
                             testData.get("residNumber").get<double>());
        }

        WHEN ("Getting Sequences") {
            string *names = new string[alig.sequenNumber];
            string *expectedNames;

            populate(expectedNames, testData.get("names"));

            GIVEN ("Out string array") {
                alig.getSequences(names);

                REQUIRE_THAT (names,
                              ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                 "./dataset/testingFiles/alignmentErrors/input_filename.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/input_filename.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/input_filename.sizePerSequence1.error.tsv"));

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
                            int *sizePerSequence = new int[alig.sequenNumber],
                                    *expectedSizePerSequence;

                            populate(expectedSizePerSequence,
                                     testData.get("noGapsSize"));

                            string *sequences = new string[alig.sequenNumber],
                                    *expectedSequences;
                            populate(expectedSequences,
                                     testData.get("noGapsSequences"));

                            alig.getSequences(names,
                                              sequences,
                                              sizePerSequence);

                            REQUIRE_THAT (names,
                                          ArrayContentsEqual(expectedNames,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/input_filename.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/input_filename.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/input_filename.sequences.error.tsv"));

                            delete[] sizePerSequence;
                            delete[] expectedSizePerSequence;

                            delete[] sequences;
                            delete[] expectedSequences;
                        }


                    } else
                        WARN ("Alignment testfile does not contain 'sequences' variable."
                                      "\nSkipping test 'Out string array, out sequences, out lengths'");
                } else
                    WARN ("Alignment testfile does not contain 'noGapsSize' variable."
                                  "\nSkipping test 'Out string array, out sequences, out lengths'");

            } else WARN ("Alignment testfile does not contain 'names' variable."
                                 "\nSkipping test 'Out string array, out sequences, out lengths'");

            delete[] names;
            delete[] expectedNames;
        }

        WHEN ("Calculate Gaps per column") {
            if (testData.contains("gapsPerColumn")) {
                if (testData.contains("aligned")) {
                    if (testData.get("aligned").get<bool>()) {
                        alig.sgaps = new statisticsGaps(&alig);
                        int *GapsPerColumn = alig.sgaps->getGapsWindow(), *expectedGapsPerColumn;
                        populate(expectedGapsPerColumn, testData.get("gapsPerColumn"));

                        REQUIRE_THAT (GapsPerColumn, ArrayContentsEqual(expectedGapsPerColumn, alig.residNumber,
                                                                        "./dataset/testingFiles/alignmentErrors/input_filename.gapsPerColumn.error.tsv"));

                        delete[] expectedGapsPerColumn;
                    } //else WARN ( "Alignment testfile is not aligned.\nSkipping test 'Calculate Gaps'" );
                } else
                    WARN ("Alignment testfile does not contain 'aligned' variable"
                                  "\nSkipping test 'Calculate Gaps'");
            } else
                WARN ("Alignment testfile does not contain 'gapsPerColumn' variable"
                              "\nSkipping test 'Calculate Gaps'");
        }

        WHEN ("Calculate No All Gaps") {
            if (testData.contains("cleanNoAllGaps")) {
                if (testData.contains("aligned")) {
                    if (testData.get("aligned").get<bool>()) {
                        int * cleanNoAllGaps;
                        populate(cleanNoAllGaps, testData.get("cleanNoAllGaps"));
                        alig.sgaps = new statisticsGaps(&alig);

                        newAlignment * newAl = alig.Cleaning->cleanNoAllGaps(false);

                        REQUIRE_THAT(newAl->saveResidues, ArrayContentsEqual(cleanNoAllGaps, alig.residNumber,
                                                                           "./dataset/testingFiles/alignmentErrors/input_filename.noAllGaps.error.tsv"));
                        delete newAl;
                        delete [] cleanNoAllGaps;
                    }
                } else
                    WARN ("Alignment testfile does not contain 'aligned' variable"
                                  "\nSkipping test 'Calculate Gaps'");
            } else
                WARN ("Alignment testfile does not contain 'gapsPerColumn' variable"
                              "\nSkipping test 'Calculate Gaps'");
        }

        WHEN ("Calculate Clean 2nd Slope") {
            if (testData.contains("clean2ndSlope")) {
                if (testData.contains("aligned")) {
                    if (testData.get("aligned").get<bool>()) {
                        int * cleanNoAllGaps;
                        populate(cleanNoAllGaps, testData.get("clean2ndSlope"));
                        alig.sgaps = new statisticsGaps(&alig);

                        newAlignment * newAl = alig.Cleaning->clean2ndSlope(false);

                        REQUIRE_THAT(newAl->saveResidues, ArrayContentsEqual(cleanNoAllGaps, alig.residNumber,
                                                                           "./dataset/testingFiles/alignmentErrors/input_filename.2ndSlope.error.tsv"));
                        delete newAl;
                        delete [] cleanNoAllGaps;
                    }
                } else
                    WARN ("Alignment testfile does not contain 'aligned' variable"
                                  "\nSkipping test 'Calculate Gaps'");
            } else
                WARN ("Alignment testfile does not contain 'gapsPerColumn' variable"
                              "\nSkipping test 'Calculate Gaps'");
        }

        WHEN ("Calculate Clean Comb Methods Lax") {
            if (testData.contains("cleanCombMethods lax")) {
                if (testData.contains("aligned")) {
                    if (testData.get("aligned").get<bool>()) {
                        int * cleanNoAllGaps;
                        populate(cleanNoAllGaps, testData.get("clean2ndSlope"));
                        alig.sgaps = new statisticsGaps(&alig);

                        newAlignment * newAl = alig.Cleaning->cleanCombMethods(false, true);

                        REQUIRE_THAT(newAl->saveResidues, ArrayContentsEqual(cleanNoAllGaps, alig.residNumber,
                                                                           "./dataset/testingFiles/alignmentErrors/input_filename.combMethodLax.error.tsv"));
                        delete newAl;
                        delete [] cleanNoAllGaps;
                    }
                } else
                    WARN ("Alignment testfile does not contain 'aligned' variable"
                                  "\nSkipping test 'Calculate Gaps'");
            } else
                WARN ("Alignment testfile does not contain 'gapsPerColumn' variable"
                              "\nSkipping test 'Calculate Gaps'");
        }

        WHEN ("Calculate Clean Comb Methods Strict") {
            if (testData.contains("cleanCombMethods strict")) {
                if (testData.contains("aligned")) {
                    if (testData.get("aligned").get<bool>()) {
                        int * cleanNoAllGaps;
                        populate(cleanNoAllGaps, testData.get("clean2ndSlope"));
                        alig.sgaps = new statisticsGaps(&alig);

                        newAlignment * newAl = alig.Cleaning->cleanCombMethods(false, false);

                        REQUIRE_THAT(newAl->saveResidues, ArrayContentsEqual(cleanNoAllGaps, alig.residNumber,
                                                                           "./dataset/testingFiles/alignmentErrors/input_filename.combMethodStrict.error.tsv"));
                        delete newAl;
                        delete [] cleanNoAllGaps;
                    }
                } else
                    WARN ("Alignment testfile does not contain 'aligned' variable"
                                  "\nSkipping test 'Calculate Gaps'");
            } else
                WARN ("Alignment testfile does not contain 'gapsPerColumn' variable"
                              "\nSkipping test 'Calculate Gaps'");
        }


    } else
        WARN("No Expected Test Results from input_filename");

}
