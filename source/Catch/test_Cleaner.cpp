#include "../../include/Catch/catch.hpp"

#include "../../include/trimalArgumentParser.h"
#include "../../include/reportsystem.h"
#include "../../source/Catch/Matchers/ArrayMatcher.cpp"

SCENARIO("Cleaner module can apply transformations to alignments", "[cleaner][alignment]")
{
    newAlignment alig;

    debug.Level = VerboseLevel::NONE;
    
    alig.sequences = new string [6];
    
    alig.sequenNumber = 6;
    alig.originalSequenNumber = 6;
    
    alig.residNumber = 60;
    alig.originalResidNumber = 60;
    
    alig.sequences[0] = "--------FAYTAPD---LLLIGFLLKTVA-T-FG--DTWF-----QLWQGLDLNKMPVF";
    alig.sequences[1] = "----------DPAVL----FV--IMLGTIT-K-FS--SEWF-----FAWLGLEINMMVII";
    alig.sequences[2] = "----------GLGKV---IVY-GIVLGTKS-DQFSNWVVWL-----FPWNGLQIHMMGII";
    alig.sequences[3] = "-----------PTIL---NIA-GLHMETDI-N-FS--LAWF-----QAWGGLEINKQAIL";
    alig.sequences[4] = "----------ASGAI---LTL-GIYLFTLC-AVIS--VSWY-----LAWLGLEINMMAII";
    alig.sequences[5] = "AAAAAAAA----ALL---TYL-GLFLGTDY-----EN---FAAAAANAWLGLEINMMAQI";
    
    alig.fillMatrices(true);
    
    int * saveResidues = new int[60];
    for (int i = 0; i < 60; i++)
        saveResidues[i] = i;
    
    REQUIRE_THAT(alig.saveResidues, ArrayContentsEqual(saveResidues, alig.residNumber));
    
    delete [] saveResidues;
    
    WHEN("Selecting a method")
    {
        THEN("Strict is selected")
            REQUIRE(alig.Cleaning->selectMethod() == STRICT);
    }
    
    WHEN("Calculate Gap Stats")
    {
        if (alig.sgaps == NULL)
            alig.sgaps = new statisticsGaps(&alig);
        THEN("Window 0")
        {
            alig.sgaps->applyWindow(0);
            int * vals = alig.sgaps->getGapsWindow();
            int * knownVals = new int[60] { 5,5,5,5,5,5,5,5,5,5,2,1,0,0,0,6,6,6,1,0,0,5,1,0,0,0,0,0,0,0,6,1,4,1,1,4,4,1,1,1,0,5,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
            
            CAPTURE(std::vector<int>(vals, vals + 60));
            CAPTURE(std::vector<int>(knownVals, knownVals + 60));
            REQUIRE_THAT(vals, ArrayContentsEqual(knownVals, 60));
            
            delete [] knownVals;
        }
    }
    
    WHEN("CleanByCutValueOverpass")
    {
        if (alig.sgaps == NULL)
            alig.sgaps = new statisticsGaps(&alig);
        alig.sgaps->applyWindow(0);
        
        newAlignment * newAl = alig.Cleaning->cleanByCutValueOverpass(0.5, 0.5, alig.sgaps->getGapsWindow(), false);
        
        int * knownVals = new int [60]
        {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,12,13,14,-1,-1,-1,-1,19,20,-1,-1,23,24,25,26,27,28,29,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,40,-1,-1,-1,-1,-1,46,47,48,49,50,51,52,53,54,55,56,57,58,59};
        
        CAPTURE(std::vector<int>(newAl->saveResidues, newAl->saveResidues + 60));
        CAPTURE(std::vector<int>(knownVals, knownVals + 60));
        REQUIRE_THAT(newAl->saveResidues, ArrayContentsEqual(knownVals, 60));
        
        delete [] knownVals;
        delete newAl;
    }
    
    WHEN("Calculate Sim Stats")
    {
        if (alig.scons == NULL)
        {
            similarityMatrix * sm = new similarityMatrix();
            sm->defaultNTSimMatrix();
            alig.scons = new statisticsConservation2(&alig);
            alig.scons->setSimilarityMatrix(sm);
            alig.Statistics->calculateConservationStats();
            alig.scons->applyWindow(0);
            
        }
        
        THEN("MDK_Vector")
        {
            float * vals = alig.scons->getMdkwVector();
            float * knownCons = new float[60]{0,0,0,0,0,0,0,0,0,0,1,1,0.35546314716339111328,1,1,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
            
            CAPTURE(std::vector<float>(vals, vals + 60));
            CAPTURE(std::vector<float>(knownCons, knownCons + 60));
            REQUIRE_THAT(vals, ArrayContentsEqual(knownCons, 60));
            
            delete [] knownCons;
        }
        
        THEN("Sim Values")
        {
            newAlignment * newAl = alig.Cleaning->cleanByCutValueFallBehind(0.5,0.5, alig.scons->getMdkwVector(), false);
            
            int * knownVals = new int[60]{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,10,11,-1,13,14,-1,-1,-1,18,19,20,-1,22,23,24,25,26,27,28,29,-1,31,32,33,34,35,36,37,38,39,40,-1,-1,-1,-1,-1,46,47,48,49,50,51,52,53,54,55,56,57,58,59};
            int * vals = alig.saveResidues;
            
            CAPTURE(std::vector<int>(newAl->saveResidues, newAl->saveResidues + 60));
            CAPTURE(std::vector<int>(knownVals, knownVals + 60));
            REQUIRE_THAT(newAl->saveResidues, ArrayContentsEqual(knownVals, 60));
            
            delete [] knownVals;
            delete newAl;
        }
    }
    
    WHEN("cleanByCutValue Composed")
    {
        if (alig.scons == NULL)
        {
            similarityMatrix * sm = new similarityMatrix();
            sm->defaultNTSimMatrix();
            alig.scons = new statisticsConservation2(&alig);
            alig.scons->setSimilarityMatrix(sm);
            alig.Statistics->calculateConservationStats();
            alig.scons->applyWindow(0);
        }
        
        newAlignment * newAl = alig.Cleaning->cleanByCutValueOverpassOrEquals(0.5, alig.sgaps->getGapsWindow(), 0.5, 0.5, alig.scons->getMdkwVector(), false);
        
        int * knownVals = new int[60]{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,13,14,-1,-1,-1,-1,19,20,-1,-1,23,24,25,26,27,28,29,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,40,-1,-1,-1,-1,-1,46,47,48,49,50,51,52,53,54,55,56,57,58,59};
        
        CAPTURE(std::vector<int>(newAl->saveResidues, newAl->saveResidues + 60));
        CAPTURE(std::vector<int>(knownVals, knownVals + 60));
        REQUIRE_THAT(newAl->saveResidues, ArrayContentsEqual(knownVals, 60));
        
        delete newAl;
        delete [] knownVals;
    }

}
    


