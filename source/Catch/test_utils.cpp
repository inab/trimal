#include "../../include/Catch/catch.hpp"

#include "../../include/trimalArgumentParser.h"
#include "../../include/reportsystem.h"

#include "../include/utils.h"
#include "Matchers/ArrayMatcher.cpp"

SCENARIO("Multiple Utilities", "[Utils]")
{
    WHEN("Initialize a vector with an int value")
    {
        int * vector = new int[20];
        
        THEN("Array is filled with 0")
        {
            int * vector_to = new int[20] { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
            utils::initlVect(vector, 20, 0);
            REQUIRE_THAT(vector_to, ArrayContentsEqual( vector, 20));
            delete [] vector_to;
        }
        
        THEN("Array is filled with 1")
        {
            int * vector_to = new int[20] { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
            utils::initlVect(vector, 20, 1);
            REQUIRE_THAT(vector_to, ArrayContentsEqual( vector, 20));
            delete [] vector_to;
        }
        
        delete [] vector;
    }
    
    WHEN("Initialize a vector with a float value")
    {
        float * vector = new float[20];
        
        THEN("Array is filled with 0F")
        {
            float * vector_to = new float[20] { 0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F };
            utils::initlVect(vector, 20, 0.F);
            REQUIRE_THAT(vector_to, ArrayContentsEqual( vector, 20));
            delete [] vector_to;
        }
        
        THEN("Array is filled with 1F")
        {
            float * vector_to = new float[20] { 1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F };
            utils::initlVect(vector, 20, 1.F);
            REQUIRE_THAT(vector_to, ArrayContentsEqual( vector, 20));
            delete [] vector_to;
        }
        
        delete [] vector;
    }
    
    WHEN("Copy a vector")
    {
        
        THEN("Float arrays equal")
        {
            float * vector = new float[20];
            
            float * vector_to = new float[20] { 0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F };
            
            utils::initlVect(vector, 20, 0.F);
            REQUIRE_THAT(vector_to, ArrayContentsEqual( vector, 20));
            
            delete [] vector_to;
            delete [] vector;
        }
        
        THEN("Int arrays equal")
        {
            int * vector = new int[20];
            
            int * vector_to = new int[20] { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
            
            utils::copyVect(vector_to, vector, 20);
            
            REQUIRE_THAT(vector_to, ArrayContentsEqual( vector, 20));
            
            delete [] vector_to;
            delete [] vector;
        }
        
    }
    
    WHEN("Rounding")
    {
        THEN("Value -1.9")
        {
            double value = -1.9;
            CAPTURE(value);
            CHECK(utils::roundInt(value) == -2);
            CHECK(utils::roundToInf(value) == -2);
            CHECK(utils::roundToSup(value) == -1);
        }
        
        THEN("Value -1.5")
        {
            double value = -1.5;
            CAPTURE(value);
            CHECK(utils::roundInt(value) == -2);
            CHECK(utils::roundToInf(value) == -2);
            CHECK(utils::roundToSup(value) == -1);
        }
        
        THEN("Value -1.1")
        {
            double value = -1.1;
            CAPTURE(value);
            CHECK(utils::roundInt(value) == -1);
            CHECK(utils::roundToInf(value) == -2);
            CHECK(utils::roundToSup(value) == -1);
        }
        
        THEN("Value -1")
        {
            double value = -1;
            CAPTURE(value);
            CHECK(utils::roundInt(value) == -1);
            CHECK(utils::roundToInf(value) == -1);
            CHECK(utils::roundToSup(value) == -1);
        }
        
        THEN("Value 0")
        {
            double value = -0;
            CAPTURE(value);
            CHECK(utils::roundInt(value) == 0);
            CHECK(utils::roundToInf(value) == 0);
            CHECK(utils::roundToSup(value) == 0);
        }
        
        THEN("Value 1.0")
        {
            double value = 1.0;
            CAPTURE(value);
            CHECK(utils::roundInt(value) == 1);
            CHECK(utils::roundToInf(value) == 1);
            CHECK(utils::roundToSup(value) == 1);
        }
        
        THEN("Value 1.1")
        {
            double value = 1.1;
            CAPTURE(value);
            CHECK(utils::roundInt(value) == 1);
            CHECK(utils::roundToInf(value) == 1);
            CHECK(utils::roundToSup(value) == 2);
        }
        
        THEN("Value 1.5")
        {
            double value = 1.5;
            CAPTURE(value);
            CHECK(utils::roundInt(value) == 2);
            CHECK(utils::roundToInf(value) == 1);
            CHECK(utils::roundToSup(value) == 2);
        }
        
        THEN("Value 1.9")
        {
            double value = 1.9;
            CAPTURE(value);
            CHECK(utils::roundInt(value) == 2);
            CHECK(utils::roundToInf(value) == 1);
            CHECK(utils::roundToSup(value) == 2);
        }
    }
    
    WHEN("Comparing numbers")
    {
        THEN("Int values 0 : 1")
        {
            CHECK(utils::max(0, 1) == 1);
            CHECK(utils::min(0, 1) == 0);
        }
        
        THEN("Float values 0.F : 1.F")
        {
            CHECK(utils::max(0.F, 1.F) == 1.F);
            CHECK(utils::min(0.F, 1.F) == 0.F);
        }
                
        THEN("Double values 0. : 1.")
        {
            CHECK(utils::max(0., 1.) == 1.);
            CHECK(utils::min(0., 1.) == 0.);
        }
    }
}
    


