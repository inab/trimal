#include "Catch/catch.hpp"
#include "Catch/picojson.h"
#include "Catch/catchhelperfunctions.h"

#include "trimalArgumentParser.h"
#include "Matchers/ArrayMatcher.cpp"
#include <vector>

#include <fstream>

int i;

// Map of old and new ID's of sequence types
std::map<int,int> AlignmentTypeCorrespondence = {
        {1, SequenceTypes::DNA},
        {2, SequenceTypes::RNA},
        {3, SequenceTypes::AA},
        {4, SequenceTypes::DNA | SequenceTypes::DEG},
        {5, SequenceTypes::RNA | SequenceTypes::DEG}
};

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

        WHEN ("Checking save residues is initialized correctly") {

            if (testData.contains("aligned") && testData.get("aligned").get<bool>()) {

                int *saveResidues = new int[alig.residNumber];
                for (i = 0; i < alig.residNumber; i++) {
                    saveResidues[i] = i;
                }

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
                              ArrayContentsEqual(saveResidues,
                                                 alig.residNumber,
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

        WHEN("Choosing an automated method")
        {
            if (testData.contains("aligned")) {
                if (testData.get("aligned").get<bool>()) {
                    // TODO Clean this. We shouldn't be casting. The original value should be in int, not string.
                    REQUIRE(std::to_string(alig.Cleaning->selectMethod()) == testData.get("AutoMethod").get<string>());
                }
            }
        }

        WHEN("Calculating Cons Vectors")
        {
            if (testData.contains("aligned")) {
                if (testData.get("aligned").get<bool>()) {
                    alig.scons = new statisticsConservation2(&alig);
                    alig.sgaps = new statisticsGaps(&alig);
                    int SequenceType = testData.get("SequenceType").get<double>();
                    similarityMatrix *sm = new similarityMatrix();
                    switch (SequenceType) {
                        case 1: case 2: default:
                            sm->defaultNTSimMatrix(); break;
                        case 3:
                            sm->defaultAASimMatrix(); break;
                        case 4: case 5:
                            sm->defaultNTDegeneratedSimMatrix(); break;
                    }
                    alig.scons->setSimilarityMatrix(sm);
                    alig.scons->calculateVectors(alig.sgaps->getGapsWindow());

                    WHEN("Calculating MDK")
                    {
                        float * mdk_expected;
                        populate(mdk_expected, testData.get("mdk"));

                        float * mdk_obtained = alig.scons->getMdkwVector();

                        REQUIRE(mdk_obtained != NULL);

                        CHECK_THAT(mdk_obtained,
                                     ArrayContentsEqual(mdk_expected,
                                                        alig.residNumber,
                                                        "./dataset/testingFiles/alignmentErrors/input_filename.MDK.error.tsv"));
                        delete [] mdk_expected;
//                        delete [] mdk_obtained;
                    }

                    WHEN("Calculating Q")
                    {
                        float * Q_expected;
                        populate(Q_expected, testData.get("Q"));

                        float * Q_obtained = alig.scons->Q;

                        REQUIRE(Q_obtained != NULL);

                        CHECK_THAT(Q_expected,
                                     ArrayContentsEqual(Q_obtained,
                                                        alig.residNumber,
                                                        "./dataset/testingFiles/alignmentErrors/input_filename.Q.error.tsv"));
                        delete [] Q_expected;
//                        delete [] Q_obtained;
                    }

                    delete sm;
                }
            }
        }

        WHEN("Getting alignment type")
        {
            // TODO Result JSON should include the new mapping ids of SequenceTypes instead of the original.
            REQUIRE(alig.getAlignmentType() == AlignmentTypeCorrespondence[static_cast<int>(testData.get("SequenceType").get<double>())]);
        }

        WHEN ("Getting Sequences") {

            if (testData.contains("names")) {
                string *names = new string[alig.sequenNumber],
                        *expectedNames;
                populate(expectedNames, testData.get("names"));

                GIVEN ("Out names") {

                    alig.getSequences(names);

                    CHECK_THAT (names,
                                  ArrayContentsEqual(expectedNames,
                                                     alig.sequenNumber,
                                                     "./dataset/testingFiles/alignmentErrors/input_filename.names1.error.tsv"));
                }

                if (testData.contains("noGapsSize")) {
                    int *sizePerSequence = new int[alig.sequenNumber],
                            *expectedSizePerSequence;
                    populate(expectedSizePerSequence, testData.get("noGapsSize"));

                    GIVEN ("Out names and out lengths") {

                        alig.getSequences(names,
                                          sizePerSequence);

                        CHECK_THAT (names,
                                      ArrayContentsEqual(expectedNames,
                                                         alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/input_filename.names2.error.tsv"));

                        CHECK_THAT (sizePerSequence,
                                      ArrayContentsEqual(expectedSizePerSequence,
                                                         alig.sequenNumber,
                                                         "./dataset/testingFiles/alignmentErrors/input_filename.sizePerSequence1.error.tsv"));
                    }

                    if (testData.contains("noGapsSequences")) {
                        string *sequences = new string[alig.sequenNumber],
                                *expectedSequences;
                        populate(expectedSequences,
                                 testData.get("noGapsSequences"));

                        GIVEN ("Out names, out sequences, out lengths") {
                            alig.getSequences(names,
                                              sequences,
                                              sizePerSequence);

                            CHECK_THAT (names,
                                          ArrayContentsEqual(expectedNames,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/input_filename.names4.error.tsv"));

                            CHECK_THAT (sizePerSequence,
                                          ArrayContentsEqual(expectedSizePerSequence,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/input_filename.sizePerSequence2.error.tsv"));

                            CHECK_THAT (sequences,
                                          ArrayContentsEqual(expectedSequences,
                                                             alig.sequenNumber,
                                                             "./dataset/testingFiles/alignmentErrors/input_filename.sequences.error.tsv"));
                        }

                        delete[] sequences;
                        delete[] expectedSequences;
                    } else
                        WARN ("Alignment testfile does not contain 'sequences' variable."
                                      "\nSkipping test 'Out string array, out sequences, out lengths'");
                    delete[] sizePerSequence;
                    delete[] expectedSizePerSequence;
                } else
                    WARN ("Alignment testfile does not contain 'noGapsSize' variable."
                                  "\nSkipping test 'Out string array and out lengths'");
                delete[] names;
                delete[] expectedNames;
            } else
                WARN ("Alignment testfile does not contain 'names' variable."
                              "\nSkipping test 'Out string array'");


        }

        WHEN ("Calculating Gaps per column") {
            if (testData.contains("gapsPerColumn")) {
                if (testData.contains("aligned")) {
                    if (testData.get("aligned").get<bool>()) {
                        alig.sgaps = new statisticsGaps(&alig);
                        int *GapsPerColumn = alig.sgaps->getGapsWindow(), *expectedGapsPerColumn;
                        populate(expectedGapsPerColumn, testData.get("gapsPerColumn"));

                        CHECK_THAT (GapsPerColumn,
                                      ArrayContentsEqual(expectedGapsPerColumn,
                                                         alig.residNumber,
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

        WHEN ("Trimming by No All Gaps") {
            if (testData.contains("cleanNoAllGaps")) {
                if (testData.contains("aligned")) {
                    if (testData.get("aligned").get<bool>()) {
                        alig.sgaps = new statisticsGaps(&alig);
                        int *cleanNoAllGaps;
                        populate(cleanNoAllGaps, testData.get("cleanNoAllGaps"));

                        newAlignment *newAl = alig.Cleaning->cleanNoAllGaps(false);

                        CHECK_THAT(newAl->saveResidues,
                                     ArrayContentsEqual(cleanNoAllGaps,
                                                        alig.residNumber,
                                                        "./dataset/testingFiles/alignmentErrors/input_filename.noAllGaps.error.tsv"));
                        delete newAl;
                        delete[] cleanNoAllGaps;
                    }
                } else
                    WARN ("Alignment testfile does not contain 'aligned' variable"
                                  "\nSkipping test 'Calculate Gaps'");
            } else
                WARN ("Alignment testfile does not contain 'gapsPerColumn' variable"
                              "\nSkipping test 'Calculate Gaps'");
        }

        WHEN ("Trimming by Clean 2nd Slope") {
            if (testData.contains("clean2ndSlope")) {
                if (testData.contains("aligned")) {
                    if (testData.get("aligned").get<bool>()) {
                        alig.sgaps = new statisticsGaps(&alig);

                        WHEN("Calculating 2nd Slope")
                            CHECK(alig.sgaps->calcCutPoint2ndSlope() == testData.get("2ndSlopeCutPoint").get<double>());

                        WHEN("Applying 2nd Slope")
                        {
                            int *cleanNoAllGaps;
                            populate(cleanNoAllGaps, testData.get("clean2ndSlope"));

                            newAlignment *newAl = alig.Cleaning->clean2ndSlope(false);

                            CHECK_THAT(newAl->saveResidues,
                                         ArrayContentsEqual(cleanNoAllGaps,
                                                            alig.residNumber,
                                                            "./dataset/testingFiles/alignmentErrors/input_filename.2ndSlope.error.tsv"));
                            delete newAl;
                            delete[] cleanNoAllGaps;
                        }
                    }
                } else
                    WARN ("Alignment testfile does not contain 'aligned' variable"
                                  "\nSkipping test 'Calculate Gaps'");
            } else
                WARN ("Alignment testfile does not contain 'gapsPerColumn' variable"
                              "\nSkipping test 'Calculate Gaps'");
        }

        WHEN ("Trimming by Clean Comb Methods Lax") {
            if (testData.contains("cleanCombMethods lax")) {
                if (testData.contains("aligned")) {
                    if (testData.get("aligned").get<bool>()) {

                        alig.sgaps = new statisticsGaps(&alig);
                        alig.scons = new statisticsConservation2(&alig);
                        int SequenceType = testData.get("SequenceType").get<double>();
                        similarityMatrix *sm = new similarityMatrix();
                        switch (SequenceType) {
                            case 1: case 2: default:
                                sm->defaultNTSimMatrix(); break;
                            case 3:
                                sm->defaultAASimMatrix(); break;
                            case 4: case 5:
                                sm->defaultNTDegeneratedSimMatrix(); break;
                        }
                        alig.scons->setSimilarityMatrix(sm);

                        int *cleanNoAllGaps;
                        populate(cleanNoAllGaps, testData.get("clean2ndSlope"));

                        newAlignment *newAl = alig.Cleaning->cleanCombMethods(false, true);

                        CHECK_THAT(newAl->saveResidues,
                                     ArrayContentsEqual(cleanNoAllGaps,
                                                        alig.residNumber,
                                                        "./dataset/testingFiles/alignmentErrors/input_filename.combMethodLax.error.tsv"));
                        delete newAl;
                        delete sm;
                        delete[] cleanNoAllGaps;
                    }
                } else
                    WARN ("Alignment testfile does not contain 'aligned' variable"
                                  "\nSkipping test 'Calculate Gaps'");
            } else
                WARN ("Alignment testfile does not contain 'gapsPerColumn' variable"
                              "\nSkipping test 'Calculate Gaps'");
        }

        WHEN ("Trimming by Clean Comb Methods Strict") {
            if (testData.contains("cleanCombMethods strict")) {
                if (testData.contains("aligned")) {
                    if (testData.get("aligned").get<bool>()) {

                        alig.sgaps = new statisticsGaps(&alig);
                        alig.scons = new statisticsConservation2(&alig);
                        int SequenceType = testData.get("SequenceType").get<double>();
                        similarityMatrix *sm = new similarityMatrix();
                        switch (SequenceType) {
                            case 1: case 2: default:
                                sm->defaultNTSimMatrix(); break;
                            case 3:
                                sm->defaultAASimMatrix(); break;
                            case 4: case 5:
                                sm->defaultNTDegeneratedSimMatrix(); break;
                        }
                        alig.scons->setSimilarityMatrix(sm);

                        int *cleanNoAllGaps;
                        populate(cleanNoAllGaps, testData.get("clean2ndSlope"));

                        newAlignment *newAl = alig.Cleaning->cleanCombMethods(false, false);

                        CHECK_THAT(newAl->saveResidues,
                                     ArrayContentsEqual(cleanNoAllGaps,
                                                        alig.residNumber,
                                                        "./dataset/testingFiles/alignmentErrors/input_filename.combMethodStrict.error.tsv"));
                        delete newAl;
                        delete sm;
                        delete[] cleanNoAllGaps;
                    }
                } else
                    WARN ("Alignment testfile does not contain 'aligned' variable"
                                  "\nSkipping test 'Calculate Gaps'");
            } else
                WARN ("Alignment testfile does not contain 'gapsPerColumn' variable"
                              "\nSkipping test 'Calculate Gaps'");
        }

        WHEN("Trimming by Clean Gaps") {
            if (testData.contains("cleanGaps")) {
                if (testData.contains("aligned") && testData.get("aligned").get<bool>()) {
                    alig.sgaps = new statisticsGaps(&alig);
                    auto cleanGapsData = testData.get("cleanGaps").get<picojson::object>();
                    for (auto baselineIT = cleanGapsData.begin(), end = cleanGapsData.end(); baselineIT != end; baselineIT++) {
                        GIVEN (baselineIT->first) {
                            auto cleanGapsDataSecondLevel = baselineIT->second.get<picojson::object>();
                            for (auto gapPctIT = cleanGapsDataSecondLevel.begin(), end2 = cleanGapsDataSecondLevel.end();
                                 gapPctIT != end2;
                                 gapPctIT++)
                            {
                                GIVEN(gapPctIT->first) {
                                    float baseline = std::stof(baselineIT->first.substr(baselineIT->first.find('=') + 1));
                                    float gapsPct = std::stof(gapPctIT->first.substr(gapPctIT->first.find('=') + 1));
                                    int *saveResidues;
                                    populate(saveResidues, gapPctIT->second);

                                    newAlignment *newAl = alig.Cleaning->cleanGaps(baseline, gapsPct, false);

                                    std::stringstream ss;
                                    ss << "./dataset/testingFiles/alignmentErrors/input_filename.CleanGaps"
                                       << ".BL" << std::to_string(baseline)
                                       << ".GPCT" << std::to_string(gapsPct)
                                       << ".error.tsv";

                                    CHECK_THAT(newAl->saveResidues,
                                                 ArrayContentsEqual(saveResidues,
                                                                    newAl->residNumber,
                                                                    ss.rdbuf()->str()));

                                    delete newAl;

                                    delete[] saveResidues;
                                }
                            }
                        }
                    }
                }
            }
        }

        WHEN("Calculating consCutPoint") {
            if (testData.contains("consCutPoint")) {
                if (testData.contains("aligned") && testData.get("aligned").get<bool>()) {
                    alig.sgaps = new statisticsGaps(&alig);
                    alig.scons = new statisticsConservation2(&alig);
                    int SequenceType = static_cast<int>(testData.get("SequenceType").get<double>());
                    similarityMatrix *sm = new similarityMatrix();
                    switch (SequenceType) {
                        case 1: case 2: default:
                            sm->defaultNTSimMatrix(); break;
                        case 3:
                            sm->defaultAASimMatrix(); break;
                        case 4: case 5:
                            sm->defaultNTDegeneratedSimMatrix(); break;
                    }
                    alig.scons->setSimilarityMatrix(sm);

                    auto cleanGapsData = testData.get("consCutPoint").get<picojson::object>();

                    for (auto baselineIT = cleanGapsData.begin(), end = cleanGapsData.end(); baselineIT != end; baselineIT++) {
                        GIVEN (baselineIT->first) {
                            auto cleanGapsDataSecondLevel = baselineIT->second.get<picojson::object>();
                            for (auto conPctIT = cleanGapsDataSecondLevel.begin(), end2 = cleanGapsDataSecondLevel.end();
                                 conPctIT != end2;
                                 conPctIT++) {
                                GIVEN(conPctIT->first) {
                                    float baseline = std::stof(baselineIT->first.substr(baselineIT->first.find('=') + 1));
                                    float consPct = std::stof(conPctIT->first.substr(conPctIT->first.find('=') + 1));
                                    alig.scons->calculateVectors(alig.sgaps->getGapsWindow());
                                    CHECK(alig.scons->calcCutPoint(baseline, consPct) == conPctIT->second.get<double>());
                                }
                            }
                        }
                    }
                    delete sm;
                }
            }
        }
    } else
        WARN("No Expected Test Results from input_filename");

}
