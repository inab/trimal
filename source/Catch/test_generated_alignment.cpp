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
SCENARIO ("Alignment methods work correctly in alignments_comparison.1", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("alignments_comparison.1", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/alignments_comparison.1.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/alignments_comparison.1.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/alignments_comparison.1.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/alignments_comparison.1.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/alignments_comparison.1.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/alignments_comparison.1.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/alignments_comparison.1.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/alignments_comparison.1.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/alignments_comparison.1.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/alignments_comparison.1.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/alignments_comparison.1.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/alignments_comparison.1.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from alignments_comparison.1");

}

SCENARIO ("Alignment methods work correctly in alignments_comparison.2", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("alignments_comparison.2", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/alignments_comparison.2.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/alignments_comparison.2.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/alignments_comparison.2.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/alignments_comparison.2.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/alignments_comparison.2.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/alignments_comparison.2.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/alignments_comparison.2.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/alignments_comparison.2.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/alignments_comparison.2.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/alignments_comparison.2.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/alignments_comparison.2.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/alignments_comparison.2.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from alignments_comparison.2");

}

SCENARIO ("Alignment methods work correctly in alignments_comparison.3", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("alignments_comparison.3", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/alignments_comparison.3.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/alignments_comparison.3.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/alignments_comparison.3.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/alignments_comparison.3.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/alignments_comparison.3.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/alignments_comparison.3.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/alignments_comparison.3.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/alignments_comparison.3.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/alignments_comparison.3.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/alignments_comparison.3.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/alignments_comparison.3.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/alignments_comparison.3.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from alignments_comparison.3");

}

SCENARIO ("Alignment methods work correctly in example.001.AA.clw", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.001.AA.clw", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.001.AA.clw.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.001.AA.clw.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.001.AA.clw.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.001.AA.clw.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.001.AA.clw.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.001.AA.clw.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.001.AA.clw.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.001.AA.clw.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.001.AA.clw.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.001.AA.clw.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.001.AA.clw.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.001.AA.clw.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.001.AA.clw");

}

SCENARIO ("Alignment methods work correctly in example.001.AA.msl", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.001.AA.msl", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.001.AA.msl.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.001.AA.msl.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.001.AA.msl.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.001.AA.msl.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.001.AA.msl.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.001.AA.msl.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.001.AA.msl.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.001.AA.msl.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.001.AA.msl.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.001.AA.msl.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.001.AA.msl.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.001.AA.msl.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.001.AA.msl");

}

SCENARIO ("Alignment methods work correctly in example.001.AA.phy", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.001.AA.phy", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.001.AA.phy.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.001.AA.phy.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.001.AA.phy.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.001.AA.phy.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.001.AA.phy.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.001.AA.phy.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.001.AA.phy.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.001.AA.phy.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.001.AA.phy.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.001.AA.phy.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.001.AA.phy.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.001.AA.phy.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.001.AA.phy");

}

SCENARIO ("Alignment methods work correctly in example.002.AA.clw", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.002.AA.clw", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.002.AA.clw.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.002.AA.clw.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.002.AA.clw.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.002.AA.clw.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.002.AA.clw.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.002.AA.clw.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.002.AA.clw.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.002.AA.clw.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.002.AA.clw.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.002.AA.clw.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.002.AA.clw.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.002.AA.clw.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.002.AA.clw");

}

SCENARIO ("Alignment methods work correctly in example.002.AA.phy", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.002.AA.phy", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.002.AA.phy.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.002.AA.phy.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.002.AA.phy.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.002.AA.phy.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.002.AA.phy.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.002.AA.phy.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.002.AA.phy.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.002.AA.phy.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.002.AA.phy.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.002.AA.phy.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.002.AA.phy.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.002.AA.phy.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.002.AA.phy");

}

SCENARIO ("Alignment methods work correctly in example.003.AA.clw", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.003.AA.clw", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.003.AA.clw.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.003.AA.clw.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.003.AA.clw.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.003.AA.clw.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.003.AA.clw.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.003.AA.clw.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.003.AA.clw.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.003.AA.clw.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.003.AA.clw.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.003.AA.clw.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.003.AA.clw.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.003.AA.clw.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.003.AA.clw");

}

SCENARIO ("Alignment methods work correctly in example.004.AA.fasta", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.004.AA.fasta", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.004.AA.fasta.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.004.AA.fasta.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.004.AA.fasta.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.004.AA.fasta.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.004.AA.fasta.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.004.AA.fasta.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.004.AA.fasta.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.004.AA.fasta.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.004.AA.fasta.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.004.AA.fasta.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.004.AA.fasta.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.004.AA.fasta.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.004.AA.fasta");

}

SCENARIO ("Alignment methods work correctly in example.005.AA.fasta", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.005.AA.fasta", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.005.AA.fasta.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.005.AA.fasta.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.005.AA.fasta.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.005.AA.fasta.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.005.AA.fasta.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.005.AA.fasta.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.005.AA.fasta.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.005.AA.fasta.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.005.AA.fasta.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.005.AA.fasta.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.005.AA.fasta.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.005.AA.fasta.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.005.AA.fasta");

}

SCENARIO ("Alignment methods work correctly in example.006.AA.pir", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.006.AA.pir", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.006.AA.pir.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.006.AA.pir.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.006.AA.pir.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.006.AA.pir.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.006.AA.pir.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.006.AA.pir.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.006.AA.pir.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.006.AA.pir.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.006.AA.pir.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.006.AA.pir.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.006.AA.pir.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.006.AA.pir.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.006.AA.pir");

}

SCENARIO ("Alignment methods work correctly in example.007.AA.fasta", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.007.AA.fasta", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.007.AA.fasta.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.007.AA.fasta.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.007.AA.fasta.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.007.AA.fasta.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.007.AA.fasta.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.007.AA.fasta.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.007.AA.fasta.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.007.AA.fasta.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.007.AA.fasta.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.007.AA.fasta.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.007.AA.fasta.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.007.AA.fasta.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.007.AA.fasta");

}

SCENARIO ("Alignment methods work correctly in example.007.AA.only_seqs", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.007.AA.only_seqs", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.007.AA.only_seqs.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.007.AA.only_seqs.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.007.AA.only_seqs.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.007.AA.only_seqs.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.007.AA.only_seqs.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.007.AA.only_seqs.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.007.AA.only_seqs.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.007.AA.only_seqs.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.007.AA.only_seqs.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.007.AA.only_seqs.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.007.AA.only_seqs.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.007.AA.only_seqs.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.007.AA.only_seqs");

}

SCENARIO ("Alignment methods work correctly in example.009.AA.fasta", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.009.AA.fasta", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.009.AA.fasta.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.009.AA.fasta.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.009.AA.fasta.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.009.AA.fasta.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.009.AA.fasta.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.009.AA.fasta.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.009.AA.fasta.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.009.AA.fasta.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.009.AA.fasta.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.009.AA.fasta.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.009.AA.fasta.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.009.AA.fasta.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.009.AA.fasta");

}

SCENARIO ("Alignment methods work correctly in example.010.AA.fasta", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.010.AA.fasta", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.010.AA.fasta.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.010.AA.fasta.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.010.AA.fasta.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.010.AA.fasta.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.010.AA.fasta.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.010.AA.fasta.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.010.AA.fasta.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.010.AA.fasta.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.010.AA.fasta.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.010.AA.fasta.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.010.AA.fasta.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.010.AA.fasta.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.010.AA.fasta");

}

SCENARIO ("Alignment methods work correctly in example.011.AA.YKL197C.clw", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.011.AA.YKL197C.clw", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.clw.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.clw.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.clw.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.clw.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.clw.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.clw.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.clw.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.clw.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.clw.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.clw.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.clw.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.clw.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.011.AA.YKL197C.clw");

}

SCENARIO ("Alignment methods work correctly in example.011.AA.YKL197C.fasta", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.011.AA.YKL197C.fasta", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.fasta.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.fasta.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.fasta.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.fasta.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.fasta.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.fasta.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.fasta.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.fasta.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.fasta.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.fasta.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.fasta.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.fasta.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.011.AA.YKL197C.fasta");

}

SCENARIO ("Alignment methods work correctly in example.011.AA.YKL197C.phy", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.011.AA.YKL197C.phy", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.phy.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.phy.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.phy.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.phy.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.phy.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.phy.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.phy.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.phy.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.phy.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.phy.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.phy.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.011.AA.YKL197C.phy.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.011.AA.YKL197C.phy");

}

SCENARIO ("Alignment methods work correctly in example.012.AA.SuperAlignment.phy", "[alignment][aligMethods]") {
    static picojson::value testData;

    newAlignment alig;

    if (loadJSON("example.012.AA.SuperAlignment.phy", testData)) {

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
                                      "./dataset/testingFiles/alignmentErrors/example.012.AA.SuperAlignment.phy.saveResidues.error.tsv"));

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
                                                 "./dataset/testingFiles/alignmentErrors/example.012.AA.SuperAlignment.phy.names1.error.tsv"));
            }

            GIVEN ("Out string array and out lengths") {
                if (testData.contains("names")) {
                    if (testData.contains("noGapsSize")) {
                        int *sizePerSequence = new int[alig.sequenNumber], *expectedSizePerSequence;
                        populate(expectedSizePerSequence, testData.get("noGapsSize"));

                        alig.getSequences(names, sizePerSequence);

                        REQUIRE_THAT (names,
                                      ArrayContentsEqual(expectedNames, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.012.AA.SuperAlignment.phy.names2.error.tsv"));

                        REQUIRE_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence, alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/example.012.AA.SuperAlignment.phy.sizePerSequence1.error.tsv"));

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
                                                             "./dataset/testingFiles/alignmentErrors/example.012.AA.SuperAlignment.phy.names4.error.tsv"));

                            REQUIRE_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.012.AA.SuperAlignment.phy.sizePerSequence2.error.tsv"));

                            REQUIRE_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/example.012.AA.SuperAlignment.phy.sequences.error.tsv"));

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
                                                                        "./dataset/testingFiles/alignmentErrors/example.012.AA.SuperAlignment.phy.gapsPerColumn.error.tsv"));

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
                                                                           "./dataset/testingFiles/alignmentErrors/example.012.AA.SuperAlignment.phy.noAllGaps.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.012.AA.SuperAlignment.phy.2ndSlope.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.012.AA.SuperAlignment.phy.combMethodLax.error.tsv"));
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
                                                                           "./dataset/testingFiles/alignmentErrors/example.012.AA.SuperAlignment.phy.combMethodStrict.error.tsv"));
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
        WARN("No Expected Test Results from example.012.AA.SuperAlignment.phy");

}

