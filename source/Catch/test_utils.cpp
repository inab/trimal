#include "../../include/Catch/catch.hpp"

#include "../../include/trimalArgumentParser.h"
#include "../../include/reportsystem.h"

#include "../include/utils.h"
#include "Matchers/ArrayMatcher.cpp"

SCENARIO ( "Array utils", "[utils][array]" ) {
    WHEN ( "Initialize a vector with an int value" ) {
        int * vector = new int[20];

        THEN ( "Array is filled with 0" ) {
            int * vector_to = new int[20] { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
            utils::initlVect ( vector, 20, 0 );

            CAPTURE ( std::vector<int> ( vector, vector + 20 ) );
            CAPTURE ( std::vector<int> ( vector_to, vector_to + 20 ) );
            REQUIRE_THAT ( vector_to, ArrayContentsEqual ( vector, 20 ) );
            delete [] vector_to;
        }

        AND_THEN ( "Array is filled with 1" ) {
            int * vector_to = new int[20] { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
            utils::initlVect ( vector, 20, 1 );

            CAPTURE ( std::vector<int> ( vector, vector + 20 ) );
            CAPTURE ( std::vector<int> ( vector_to, vector_to + 20 ) );
            REQUIRE_THAT ( vector_to, ArrayContentsEqual ( vector, 20 ) );
            delete [] vector_to;
        }

        delete [] vector;
    }

    WHEN ( "Initialize a vector with a float value" ) {
        float * vector = new float[20];

        THEN ( "Array is filled with 0F" ) {
            float * vector_to = new float[20] { 0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F };
            utils::initlVect ( vector, 20, 0.F );

            CAPTURE ( std::vector<int> ( vector, vector + 20 ) );
            CAPTURE ( std::vector<int> ( vector_to, vector_to + 20 ) );
            REQUIRE_THAT ( vector_to, ArrayContentsEqual ( vector, 20 ) );
            delete [] vector_to;
        }

        AND_THEN ( "Array is filled with 1F" ) {
            float * vector_to = new float[20] { 1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F,1.F };
            utils::initlVect ( vector, 20, 1.F );

            CAPTURE ( std::vector<int> ( vector, vector + 20 ) );
            CAPTURE ( std::vector<int> ( vector_to, vector_to + 20 ) );
            REQUIRE_THAT ( vector_to, ArrayContentsEqual ( vector, 20 ) );
            delete [] vector_to;
        }

        delete [] vector;
    }

    WHEN ( "Copy a vector" ) {

        GIVEN ( "Float arrays equal" ) {
            float * vector = new float[20];

            float * vector_to = new float[20] { 0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F };

            utils::initlVect ( vector, 20, 0.F );

            CAPTURE ( std::vector<int> ( vector, vector + 20 ) );
            CAPTURE ( std::vector<int> ( vector_to, vector_to + 20 ) );
            REQUIRE_THAT ( vector_to, ArrayContentsEqual ( vector, 20 ) );

            delete [] vector_to;
            delete [] vector;
        }

        GIVEN ( "Int arrays equal" ) {
            int * vector = new int[20];

            int * vector_to = new int[20] { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };

            utils::copyVect ( vector_to, vector, 20 );

            CAPTURE ( std::vector<int> ( vector, vector + 20 ) );
            CAPTURE ( std::vector<int> ( vector_to, vector_to + 20 ) );
            REQUIRE_THAT ( vector_to, ArrayContentsEqual ( vector, 20 ) );

            delete [] vector_to;
            delete [] vector;
        }
        
        WHEN ( "Compare char arrays" ) {
            char * A, * B;
            GIVEN ( "Unitary arrays" ) {
                THEN ( "Are equal" ) {
                    A = new char[2] {'A'};
                    B = new char[2] {'A'};
                    REQUIRE ( utils::compare ( A,B ) );
                    delete [] A;
                    delete [] B;
                }
                THEN ( "Are not equal" ) {
                    A = new char[2] {'A'};
                    B = new char[2] {'B'};
                    REQUIRE_FALSE ( utils::compare ( A,B ) );
                    delete [] A;
                    delete [] B;
                }
            }

            GIVEN ( "Non unitary arrays" ) {
                THEN ( "Are equal" ) {
                    A = new char[3] {'A', 'B'};
                    B = new char[3] {'A', 'B'};
                    REQUIRE ( utils::compare ( A,B ) );
                    delete [] A;
                    delete [] B;
                }
                THEN ( "Are not equal" ) {
                    A = new char[3] {'A', 'C'};
                    B = new char[3] {'B', 'D'};
                    REQUIRE_FALSE ( utils::compare ( A,B ) );
                    delete [] A;
                    delete [] B;
                }
            }
            GIVEN ( "Different sized arrays" ) {
                THEN ( "Are not equal" ) {
                    A = new char[3] {'A', 'B'};
                    B = new char[2] {'B'};
                    REQUIRE_FALSE ( utils::compare ( A,B ) );
                    delete [] A;
                    delete [] B;
                }
            }
        }
    }
}

SCENARIO ( "Number utils", "[utils][number]" ) {
    WHEN ( "Rounding" ) {
        
        struct testStruct
        {
            float valueToTest;
            int roundExpected;
            int floorExpected;
            int ceilExpected;
            
            testStruct(float value, int round, int floor, int ceil) : valueToTest(value), roundExpected(round), floorExpected(floor), ceilExpected(ceil) { } 
        };
        
        std::vector<testStruct> testValues = {
            testStruct(-1.9, -2, -2, -1),
            testStruct(-1.5, -2, -2, -1),
            testStruct(-1.1, -1, -2, -1),
            testStruct(-1.0, -1, -1, -1),
            testStruct(-0.0, +0, +0, +0),
            testStruct(+1.0, +1, +1, +1),
            testStruct(+1.1, +1, +1, +2),
            testStruct(+1.5, +2, +1, +2),
            testStruct(+1.9, +2, +1, +2),
        };
        
        for (testStruct & testValue : testValues)
        {
            WHEN("Number" + std::to_string(testValue.valueToTest))
            {
                CAPTURE ( testValue.valueToTest );
                CHECK ( utils::roundInt ( testValue.valueToTest ) == testValue.roundExpected );
                CHECK ( utils::roundToInf ( testValue.valueToTest ) == testValue.floorExpected );
                CHECK ( utils::roundToSup ( testValue.valueToTest ) == testValue.ceilExpected );
            }
        }
    }

    WHEN ( "Comparing numbers" ) {
        THEN ( "Int values 0 : 1" ) {
            CHECK ( utils::max ( 0, 1 ) == 1 );
            CHECK ( utils::min ( 0, 1 ) == 0 );
        }

        THEN ( "Float values 0.F : 1.F" ) {
            CHECK ( utils::max ( 0.F, 1.F ) == 1.F );
            CHECK ( utils::min ( 0.F, 1.F ) == 0.F );
        }

        THEN ( "Double values 0. : 1." ) {
            CHECK ( utils::max ( 0., 1. ) == 1. );
            CHECK ( utils::min ( 0., 1. ) == 0. );
        }
    }

    WHEN ( "Checking if string is number" ) {
        char * characters;
        GIVEN ( "'2'" ) {
            THEN ( "is a number" ) {
                characters = new char[2] {'2'};
                REQUIRE ( utils::isNumber ( characters ) );
                delete [] characters;
            }
        }

        GIVEN ( "'0'" ) {
            THEN ( "is a number" ) {
                characters = new char[2] {'0'};
                REQUIRE ( utils::isNumber ( characters ) );
                delete [] characters;
            }
        }

        GIVEN ( "Empty array" ) {
            THEN ( "is not a number" ) {
                characters = new char[1] {};
                REQUIRE_FALSE ( utils::isNumber ( characters ) );
                delete [] characters;
            }
        }

        GIVEN ( "'2.2'" ) {
            THEN ( "is a number" ) {
                characters = new char[4] {'2', '.', '2'};
                REQUIRE ( utils::isNumber ( characters ) );
                delete [] characters;
            }
        }

        GIVEN ( "'-2'" ) {
            THEN ( "is a number" ) {
                characters = new char[3] {'-', '2'};
                REQUIRE ( utils::isNumber ( characters ) );
                delete [] characters;
            }
        }

        GIVEN ( "'-2e-3'" ) {
            THEN ( "is a number" ) {
                characters = new char[6] {'-', '2', 'e', '-', '3'};
                REQUIRE ( utils::isNumber ( characters ) );
                delete [] characters;
            }
        }

        GIVEN ( "'-2E3'" ) {
            THEN ( "is a number" ) {
                characters = new char[5] {'-', '2', 'E', '3'};
                REQUIRE ( utils::isNumber ( characters ) );
                delete [] characters;
            }
        }

        GIVEN ( "String" ) {
            THEN ( "is not a number" ) {
                characters = new char[7] {'S', 't', 'r', 'i', 'n', 'g'};
                REQUIRE_FALSE ( utils::isNumber ( characters ) );
                delete [] characters;
            }
        }
    }
    WHEN ( "Reading numbers" ) {
        GIVEN ( "A number '2'" ) {
            int * numbers = new int [2] { 2, 2 };
            int * result = utils::readNumbers ( "2" );
            CAPTURE ( std::vector<int> ( numbers, numbers + 2 ) );
            REQUIRE ( result != NULL );
            CAPTURE ( std::vector<int> ( result, result + result[0] ) );
            REQUIRE_THAT ( result, ArrayContentsEqual ( numbers, 2 ) );
            delete [] numbers;
            delete [] result;
        }

        GIVEN ( "A number '22'" ) {
            int * numbers = new int [2] { 2, 22 };
            int * result = utils::readNumbers ( "22" );
            CAPTURE ( std::vector<int> ( numbers, numbers + 2 ) );
            REQUIRE ( result != NULL );
            CAPTURE ( std::vector<int> ( result, result + result[0] ) );
            REQUIRE_THAT ( result, ArrayContentsEqual ( numbers, 2 ) );
            delete [] numbers;
            delete [] result;
        }

        GIVEN ( "A number sequence '2, 3'" ) {
            int * numbers = new int [4] { 4, 2, 2, 3 };
            int * result = utils::readNumbers ( "2, 3" );
            CAPTURE ( std::vector<int> ( numbers, numbers + 4 ) );
            REQUIRE ( result != NULL );
            CAPTURE ( std::vector<int> ( result, result + result[0] ) );
            REQUIRE_THAT ( result, ArrayContentsEqual ( numbers, result[0] ) );
            delete [] numbers;
            delete [] result;
        }

        GIVEN ( "A number sequence '2,3,4'" ) {
            int * numbers = new int [6] { 6, 2, 2, 3, 3, 4 };
            int * result = utils::readNumbers ( "2,3,4" );
            CAPTURE ( std::vector<int> ( numbers, numbers + 6 ) );
            REQUIRE ( result != NULL );
            CAPTURE ( std::vector<int> ( result, result + result[0] ) );
            REQUIRE_THAT ( result, ArrayContentsEqual ( numbers, result[0] ) );
            delete [] numbers;
            delete [] result;
        }

        GIVEN ( "A number sequence '2,3-4,6'" ) {
            int * numbers = new int [6] { 6, 2, 2, 3, 4, 6 };
            int * result = utils::readNumbers ( "2,3-4,6" );
            CAPTURE ( std::vector<int> ( numbers, numbers + 6 ) );
            REQUIRE ( result != NULL );
            CAPTURE ( std::vector<int> ( result, result + result[0] ) );
            REQUIRE_THAT ( result, ArrayContentsEqual ( numbers, result[0] ) );
            delete [] numbers;
            delete [] result;
        }

        GIVEN ( "A number sequence '2,3-4,6'" ) {
            int * numbers = new int [6] { 6, 2, 2, 3, 4, 6 };
            int * result = utils::readNumbers ( "2,3-4,6" );
            CAPTURE ( std::vector<int> ( numbers, numbers + 6 ) );
            REQUIRE ( result != NULL );
            CAPTURE ( std::vector<int> ( result, result + result[0] ) );
            REQUIRE_THAT ( result, ArrayContentsEqual ( numbers, result[0] ) );
            delete [] numbers;
            delete [] result;
        }

        GIVEN ( "A number sequence with foreign symbol '2,3;4,6'" ) {
            int * result = utils::readNumbers ( "2,3;4,6" );
            REQUIRE ( result == NULL );
            if ( result != NULL )
                delete [] result;
        }

    }
}

SCENARIO ( "String utils", "[utils][string]" ) {
    WHEN ( "Removing spaces" ) {
        char * A, * B, * C;
        GIVEN ( "An spaced array" ) {
            A = new char [10] { 'H',' ', 'e', ' ', 'l', ' ', 'l', ' ', 'o' };
            B = new char [10];
            C = new char [10] { 'H','e', 'l', 'l', 'o' };
            utils::removeSpaces ( A, B );
            REQUIRE_THAT ( B, ArrayContentsEqual ( C, 6 ) );
            delete [] A;
            delete [] B;
            delete [] C;
        }

        GIVEN ( "An non spaced array" ) {
            A = new char [6] { 'H','e', 'l', 'l', 'o' };
            B = new char [6];
            utils::removeSpaces ( A, B );
            REQUIRE_THAT ( B, ArrayContentsEqual ( A, 6 ) );
            delete [] A;
            delete [] B;
        }

        GIVEN ( "An spaced array with numbers" ) {
            A = new char [10] { 'H',' ', '3', ' ', 'l', ' ', 'l', ' ', '0' };
            B = new char [10];
            C = new char [10] { 'H','3', 'l', 'l', '0' };
            utils::removeSpaces ( A, B );
            REQUIRE_THAT ( B, ArrayContentsEqual ( C, 6 ) );
            delete [] A;
            delete [] B;
            delete [] C;
        }
    }

    WHEN ( "Trimming lines" ) {
        char * lineToTrim = new char [11] { 'A', 'B', 'C', 'D', 'F', 'G', 'H', 'I', 'J', 'K' };
        GIVEN ( "A string with single paired brackets '[]'" ) {
            THEN ( "String is trimmed" ) {
                lineToTrim[1] = '[';
                lineToTrim[3] = ']';
                char * result = utils::trimLine ( lineToTrim );
                char * expectedResult = new char [8] { 'A', 'F', 'G', 'H', 'I', 'J', 'K' };
                CAPTURE ( std::vector<char> ( lineToTrim, lineToTrim + 10 ) );
                CAPTURE ( std::vector<char> ( result, result + 10 ) );
                CAPTURE ( std::vector<char> ( expectedResult, expectedResult + 10 ) );
                REQUIRE_THAT ( result, ArrayContentsEqual ( expectedResult, strlen ( expectedResult ) ) );
                delete [] result;
                delete [] expectedResult;
            }
        }

        GIVEN ( "A string with single paired quotation marks '\"\"'" ) {
            THEN ( "String is trimmed" ) {
                lineToTrim[1] = '"';
                lineToTrim[3] = '"';
                char * result = utils::trimLine ( lineToTrim );
                char * expectedResult = new char [8] { 'A', 'F', 'G', 'H', 'I', 'J', 'K' };
                CAPTURE ( std::vector<char> ( lineToTrim, lineToTrim + 10 ) );
                CAPTURE ( std::vector<char> ( result, result + 10 ) );
                CAPTURE ( std::vector<char> ( expectedResult, expectedResult + 10 ) );
                REQUIRE_THAT ( result, ArrayContentsEqual ( expectedResult, strlen ( expectedResult ) ) );
                delete [] result;
                delete [] expectedResult;
            }
        }

        GIVEN ( "A string with single paired brackets '{}'" ) {
            THEN ( "String remains unchanged" ) {
                lineToTrim[1] = '{';
                lineToTrim[3] = '}';
                char * result = utils::trimLine ( lineToTrim );
                CAPTURE ( std::vector<char> ( lineToTrim, lineToTrim + 10 ) );
                CAPTURE ( std::vector<char> ( result, result + 10 ) );
                REQUIRE_THAT ( result, ArrayContentsEqual ( lineToTrim, strlen ( lineToTrim ) ) );
                delete [] result;
            }
        }

        GIVEN ( "A string with single unpaired brackets '[]'" ) {
            THEN ( "String is not trimmed and return value is NULL" ) {
                lineToTrim[1] = '[';
                char * result = utils::trimLine ( lineToTrim );
                CAPTURE ( std::vector<char> ( lineToTrim, lineToTrim + 10 ) );
                REQUIRE ( result == NULL );
                delete [] result;
            }
        }

        GIVEN ( "A string with single unpaired quotation marks '\"\"'" ) {
            THEN ( "String is not trimmed and return value is NULL" ) {
                lineToTrim[1] = '"';
                char * result = utils::trimLine ( lineToTrim );
                CAPTURE ( std::vector<char> ( lineToTrim, lineToTrim + 10 ) );
                REQUIRE ( result == NULL );
                delete [] result;
            }
        }

        GIVEN ( "A string with single paired brackets '[]' and quotation marks '\"\"'" ) {
            THEN ( "String is trimmed" ) {
                lineToTrim[1] = '[';
                lineToTrim[3] = ']';
                lineToTrim[4] = '"';
                lineToTrim[6] = '"';
                char * result = utils::trimLine ( lineToTrim );
                char * expectedResult = new char [8] { 'A', 'I', 'J', 'K' };
                CAPTURE ( std::vector<char> ( lineToTrim, lineToTrim + 10 ) );
                CAPTURE ( std::vector<char> ( result, result + 10 ) );
                CAPTURE ( std::vector<char> ( expectedResult, expectedResult + 10 ) );
                REQUIRE_THAT ( result, ArrayContentsEqual ( expectedResult, strlen ( expectedResult ) ) );
                delete [] result;
                delete [] expectedResult;
            }
        }

        GIVEN ( "A string with single paired brackets '[]' with quotation marks '\"\"' inside the paired brackets" ) {
            THEN ( "String is not trimmed and return value is NULL" ) {
                lineToTrim[1] = '[';
                lineToTrim[4] = ']';
                lineToTrim[2] = '"';
                lineToTrim[3] = '"';
                char * result = utils::trimLine ( lineToTrim );
                char * expectedResult = new char [7] { 'A', 'G', 'H', 'I', 'J', 'K' };
                CAPTURE ( std::vector<char> ( lineToTrim, lineToTrim + 10 ) );
                CAPTURE ( std::vector<char> ( result, result + 10 ) );
                CAPTURE ( std::vector<char> ( expectedResult, expectedResult + 10 ) );
                REQUIRE_THAT ( result, ArrayContentsEqual ( expectedResult, strlen ( expectedResult ) ) );
                delete [] result;
                delete [] expectedResult;
            }
        }

        delete [] lineToTrim;
    }

    WHEN ( "Reversing a string" ) {
        GIVEN ( "Empty string" )
        THEN ( "Empty string is returned" )
        REQUIRE ( utils::getReverse ( "" ) == "" );

        GIVEN ( "ABC" )
        THEN ( "CBA is returned" )
        REQUIRE ( utils::getReverse ( "ABC" ) == "CBA" );

        GIVEN ( "Single character" )
        THEN ( "Same character is returned" )
        REQUIRE ( utils::getReverse ( "A" ) == "A" );
    }

    WHEN ( "Removing character" ) {
        GIVEN ( "Empty string" ) {
            string in_ = "";
            string out_ = utils::removeCharacter ( '-', in_ );
            REQUIRE ( out_ == "" );
        }

        GIVEN ( "A string with no removable character" ) {
            string in_ = "This string does not contain gap symbol";
            string out_ = utils::removeCharacter ( '-', in_ );
            REQUIRE ( out_ == "This string does not contain gap symbol" );
        }

        GIVEN ( "A string with removable character" ) {
            string in_ = "-This- -string- -does- -contain- -gap- -symbols-";
            string out_ = utils::removeCharacter ( '-', in_ );
            REQUIRE ( out_ == "This string does contain gap symbols" );
        }
    }

    WHEN ( "Replacing string" ) {
        string origin = "This is a test [string] and thus it {should} contain some kind of *tags*. Even [missmatched} {ones]";

        WHEN ( "In place" ) {
            string result = "This is a test variable and thus it {should} contain some kind of *tags*. Even [missmatched} {ones]";
            utils::ReplaceStringInPlace ( origin, "[string]", "variable" );
            REQUIRE ( origin == result );
        }

        WHEN ( "New string" ) {
            string expectedResult = "This is a test variable and thus it {should} contain some kind of *tags*. Even [missmatched} {ones]";
            string result = utils::ReplaceString ( origin, "[string]", "variable" );
            REQUIRE ( origin != expectedResult );
            REQUIRE ( expectedResult == result );
        }

        WHEN ( "In place multiple calls" ) {
            string result = "This is a test variable and thus it must contain some kind of TAGS. Even with some spaces";
            utils::ReplaceStringInPlace ( origin, "[string]", "variable" );
            utils::ReplaceStringInPlace ( origin, "{should}", "must" );
            utils::ReplaceStringInPlace ( origin, "*tags*", "TAGS" );
            utils::ReplaceStringInPlace ( origin, "[missmatched}", "with some" );
            utils::ReplaceStringInPlace ( origin, "{ones]", "spaces" );
            REQUIRE ( origin == result );
        }
    }
}

SCENARIO ( "Memory utils", "[utils][memory]" ) {
    WHEN ( "Swapping variables" ) {

        GIVEN ( "A pair of ints" ) {
            int A = 5,  B = 10;
            utils::swap ( &A, &B );
            REQUIRE ( A == 10 );
            REQUIRE ( B == 5 );
        }

        GIVEN ( "A pair of int arrays" ) {
            THEN ( "Only the first item is be swapped" ) {
                int * A  = new int[4] { 1, 2, 3, 4 },
                * B  = new int[4] { 5, 6, 7, 8 },

                * AA = new int[4] { 5, 2, 3, 4 },
                * BB = new int[4] { 1, 6, 7, 8 };

                utils::swap ( A, B );
                REQUIRE_THAT ( A, ArrayContentsEqual ( AA, 4 ) );
                REQUIRE_THAT ( B, ArrayContentsEqual ( BB, 4 ) );

                delete [] A;
                delete [] B;
                delete [] AA;
                delete [] BB;
            }

            THEN ( "The whole array is swapped" ) {
                int * A  = new int[4] { 1, 2, 3, 4 },
                * B  = new int[4] { 5, 6, 7, 8 },

                * AA = new int[4] { 1, 2, 3, 4 },
                * BB = new int[4] { 5, 6, 7, 8 };

                utils::swap ( &A, &B );
                REQUIRE_THAT ( A, ArrayContentsEqual ( BB, 4 ) );
                REQUIRE_THAT ( B, ArrayContentsEqual ( AA, 4 ) );

                delete [] A;
                delete [] B;
                delete [] AA;
                delete [] BB;
            }
        }

        GIVEN ( "A pair of floats" ) {
            float A = 5,  B = 10;
            utils::swap ( &A, &B );
            REQUIRE ( A == 10 );
            REQUIRE ( B == 5 );
        }

        GIVEN ( "A pair of float arrays" ) {
            THEN ( "Only the first item is be swapped" ) {
                float   * A  = new float[4] { 1, 2, 3, 4 },
                * B  = new float[4] { 5, 6, 7, 8 },

                * AA = new float[4] { 5, 2, 3, 4 },
                * BB = new float[4] { 1, 6, 7, 8 };

                utils::swap ( A, B );
                REQUIRE_THAT ( A, ArrayContentsEqual ( AA, 4 ) );
                REQUIRE_THAT ( B, ArrayContentsEqual ( BB, 4 ) );

                delete [] A;
                delete [] B;
                delete [] AA;
                delete [] BB;
            }
        }
    }

    WHEN ( "Quicksorting" ) {
        GIVEN ( "An int array" ) {
            int * A = new int [5] { 4, 3, 5, 2, 1 };
            THEN ( "We order it fully" ) {
                int * B = new int [5] { 1, 2, 3, 4, 5 };
                utils::quicksort ( A, 0, 4 );
                CAPTURE ( std::vector<int> ( A, A + 5 ) );
                CAPTURE ( std::vector<int> ( B, B + 5 ) );
                REQUIRE_THAT ( A, ArrayContentsEqual ( B, 5 ) );

                delete [] B;
            }

            THEN ( "We order it partially" ) {
                int * B = new int [5] { 4, 2, 3, 5, 1 };
                utils::quicksort ( A, 1, 3 );
                CAPTURE ( std::vector<int> ( A, A + 5 ) );
                CAPTURE ( std::vector<int> ( B, B + 5 ) );
                REQUIRE_THAT ( A, ArrayContentsEqual ( B, 5 ) );

                delete [] B;
            }

            THEN ( "We order only two positions" ) {
                int * B = new int [5] { 4, 3, 2, 5, 1 };
                utils::quicksort ( A, 2, 3 );
                CAPTURE ( std::vector<int> ( A, A + 5 ) );
                CAPTURE ( std::vector<int> ( B, B + 5 ) );
                REQUIRE_THAT ( A, ArrayContentsEqual ( B, 5 ) );
            }

            delete [] A;
        }

        GIVEN ( "A float array" ) {
            float * A = new float [5] { 4.f, 3.f, 5.f, 6.f, 2.f };

            THEN ( "We order it fully" ) {
                float * B = new float [5] { 2.f, 3.f, 4.f, 5.f, 6.f };
                utils::quicksort ( A, 0, 4 );
                CAPTURE ( std::vector<float> ( A, A + 5 ) );
                CAPTURE ( std::vector<float> ( B, B + 5 ) );
                REQUIRE_THAT ( A, ArrayContentsEqual ( B, 5 ) );
            }

            THEN ( "We order it partially" ) {
                float * B = new float [5] { 4.f, 3.f, 5.f, 6.f, 2.f };
                utils::quicksort ( A, 1, 3 );
                CAPTURE ( std::vector<float> ( A, A + 5 ) );
                CAPTURE ( std::vector<float> ( B, B + 5 ) );
                REQUIRE_THAT ( A, ArrayContentsEqual ( B, 5 ) );
            }

            THEN ( "We order only two positions" ) {
                float * B = new float [5] { 3.f, 4.f, 5.f, 6.f, 2.f };
                utils::quicksort ( A, 0, 1 );
                CAPTURE ( std::vector<float> ( A, A + 5 ) );
                CAPTURE ( std::vector<float> ( B, B + 5 ) );
                REQUIRE_THAT ( A, ArrayContentsEqual ( B, 5 ) );
            }

            delete [] A;
        }
    }
}

SCENARIO ( "File utils", "[utils][file]" ) {
    WHEN ( "Checking files functions" ) {
        GIVEN ( "Empty file" ) {
            ifstream F_;
            F_.open ( "./dataset/testingFiles/efile.txt" );

            THEN ( "We detect is empty" ) {
                REQUIRE_FALSE ( utils::checkFile ( F_ ) );
            }

            THEN ( "Reading line returns NULL" ) {
                REQUIRE ( utils::readLine ( F_ ) == NULL );
            }

            F_.close();
        }
        GIVEN ( "A file containing a line with 'Hello\\n'" ) {
            ifstream F_;
            F_.open ( "./dataset/testingFiles/testfile01.txt" );
            WHEN ( "Reading file" ) {
                char * result = utils::readLine ( F_ );
                WHEN ( "First line contains the correct data" ) {
                    char * expectedResult = new char [6] { 'H', 'e', 'l', 'l', 'o' };
                    CAPTURE ( result );
                    REQUIRE_FALSE ( result == NULL );
                    CHECK_THAT ( expectedResult, ArrayContentsEqual ( result, 5 ) );
                    delete [] expectedResult;
                }

                WHEN ( "Reading second line returns NULL" ) {
                    delete [] result;
                    result = utils::readLine ( F_ );
                    REQUIRE ( result == NULL );
                }
                if ( result != NULL )
                    delete [] result;
            }
            F_.close();
        }

        GIVEN ( "A file containing a line with 'Hello\\r'" ) {
            ifstream F_;
            F_.open ( "./dataset/testingFiles/testfile02.txt" );
            WHEN ( "Reading file" ) {
                char * result = utils::readLine ( F_ );

                WHEN ( "First line contains the correct data" ) {
                    char * expectedResult = new char [6] { 'H', 'e', 'l', 'l', 'o' };
                    CAPTURE ( result );
                    REQUIRE_FALSE ( result == NULL );
                    CHECK_THAT ( expectedResult, ArrayContentsEqual ( result, 5 ) );
                    delete [] expectedResult;
                }

                WHEN ( "Reading second line returns NULL" ) {
                    delete [] result;
                    result = utils::readLine ( F_ );
                    REQUIRE ( result == NULL );
                }
                if ( result != NULL )
                    delete [] result;
            }
            F_.close();
        }

        GIVEN ( "A file containing a line with some spaces before 'Hello\\n'" ) {
            ifstream F_;
            F_.open ( "./dataset/testingFiles/testfile03.txt" );
            WHEN ( "Reading file" ) {
                char * result = utils::readLine ( F_ );

                WHEN ( "First line contains the correct data" ) {
                    char * expectedResult = new char [6] { 'H', 'e', 'l', 'l', 'o' };
                    CAPTURE ( result );
                    REQUIRE_FALSE ( result == NULL );
                    CHECK_THAT ( expectedResult, ArrayContentsEqual ( result, 5 ) );
                    delete [] expectedResult;
                }

                WHEN ( "Reading second line returns NULL" ) {
                    delete [] result;
                    result = utils::readLine ( F_ );
                    REQUIRE ( result == NULL );
                }
                if ( result != NULL )
                    delete [] result;
            }
            F_.close();
        }

        GIVEN ( "A file containing a line with some tabs before 'Hello\\n'" ) {
            ifstream F_;
            F_.open ( "./dataset/testingFiles/testfile04.txt" );
            WHEN ( "Reading file" ) {
                char * result = utils::readLine ( F_ );

                WHEN ( "First line contains the correct data" ) {
                    char * expectedResult = new char [6] { 'H', 'e', 'l', 'l', 'o' };
                    CAPTURE ( result );
                    REQUIRE_FALSE ( result == NULL );
                    CHECK_THAT ( expectedResult, ArrayContentsEqual ( result, 5 ) );
                    delete [] expectedResult;
                }

                WHEN ( "Reading second line returns NULL" ) {
                    delete [] result;
                    result = utils::readLine ( F_ );
                    REQUIRE ( result == NULL );
                }
                if ( result != NULL )
                    delete [] result;
            }
            F_.close();
        }

        GIVEN ( "A file containing a line with some tabs and spaces before 'Hello\\n'" ) {
            ifstream F_;
            F_.open ( "./dataset/testingFiles/testfile05.txt" );
            WHEN ( "Reading file" ) {
                char * result = utils::readLine ( F_ );

                WHEN ( "First line contains the correct data" ) {
                    char * expectedResult = new char [6] { 'H', 'e', 'l', 'l', 'o' };
                    CAPTURE ( result );
                    REQUIRE_FALSE ( result == NULL );
                    CHECK_THAT ( expectedResult, ArrayContentsEqual ( result, 5 ) );
                    delete [] expectedResult;
                }

                WHEN ( "Reading second line returns NULL" ) {
                    delete [] result;
                    result = utils::readLine ( F_ );
                    REQUIRE ( result == NULL );
                }
                if ( result != NULL )
                    delete [] result;
            }
            F_.close();
        }

        GIVEN ( "A file containing spaces and tabs before a '\\n'" ) {
            ifstream F_;
            F_.open ( "./dataset/testingFiles/testfile06.txt" );
            WHEN ( "Reading file" ) {
                WHEN ( "Reading a line returns null" ) {
                    char * result = utils::readLine ( F_ );
                    REQUIRE ( result == NULL );
                    delete [] result;
                }
            }
            F_.close();
        }

        GIVEN ( "A file containing a line with 'Hello\\r\\n'" ) {
            ifstream F_;
            F_.open ( "./dataset/testingFiles/testfile07.txt" );
            WHEN ( "Reading file" ) {
                char * result = utils::readLine ( F_ );

                WHEN ( "First line contains the correct data" ) {
                    char * expectedResult = new char [6] { 'H', 'e', 'l', 'l', 'o' };
                    CAPTURE ( result );
                    REQUIRE_FALSE ( result == NULL );
                    CHECK_THAT ( expectedResult, ArrayContentsEqual ( result, 5 ) );
                    delete [] expectedResult;
                }

                WHEN ( "Reading second line returns NULL" ) {
                    delete [] result;
                    result = utils::readLine ( F_ );
                    REQUIRE ( result == NULL );
                }

                WHEN ( "Reading third line returns NULL" ) {
                    delete [] result;
                    utils::readLine ( F_ );
                    result = utils::readLine ( F_ );
                    REQUIRE ( result == NULL );
                }
                if ( result != NULL )
                    delete [] result;
            }
            F_.close();
        }
    }
}

SCENARIO ( "Alignment utils", "[utils][alignment]" ) {
    WHEN ( "Checking sequences types" ) {
        int sequenNumber = 6, residNumber = 60;
        string * sequences = new string[6];

        GIVEN ( "AA Sequences" ) {
            GIVEN ( "With gaps" ) {
                sequences[0] = "--------FAYTAPD---LLLIGFLLKTVA-T-FG--DTWF-----QLWQGLDLNKMPVF";
                sequences[1] = "----------DPAVL----FV--IMLGTIT-K-FS--SEWF-----FAWLGLEINMMVII";
                sequences[2] = "----------GLGKV---IVY-GIVLGTKS-DQFSNWVVWL-----FPWNGLQIHMMGII";
                sequences[3] = "-----------PTIL---NIA-GLHMETDI-N-FS--LAWF-----QAWGGLEINKQAIL";
                sequences[4] = "----------ASGAI---LTL-GIYLFTLC-AVIS--VSWY-----LAWLGLEINMMAII";
                sequences[5] = "AAAAAAAA----ALL---TYL-GLFLGTDY-----EN---FAAAAANAWLGLEINMMAQI";

                REQUIRE ( utils::checkAlignmentType ( sequenNumber, residNumber, sequences ) == SequenceTypes::AA );
            }
        }

        GIVEN ( "DNA Sequences" ) {
            WHEN ( "Pure" ) {
                GIVEN ( "With gaps" ) {
                    sequences[0] = "ATGCTGTGTGTCTGTA--------------------------TGTTGAAAGTGTCGTCAA";
                    sequences[1] = "ATGCTGTGTGTCTGTA--------------------------TGTTGAAAGTGTCGTCAA";
                    sequences[2] = "ATGCTGTGTGTCTGTA--------------------------TGTTGAAAGTGTCGTCAA";
                    sequences[3] = "ATGCTGTGTGTCTGTA--------------------------TGTTGAAAGTGTCGTCAA";
                    sequences[4] = "ATGCTGTGTGTCTGTA--------------------------TGTTGAAAGTGTCGTCAA";
                    sequences[5] = "ATGCTGTGTGTCTGTA--------------------------TGTTGAAAGTGTCGTCAA";

                    REQUIRE ( utils::checkAlignmentType ( sequenNumber, residNumber, sequences ) == SequenceTypes::DNA );
                }

                GIVEN ( "Without gaps" ) {
                    sequences[0] = "ATGCTGTGTGTCTGTAATGCTGTGTGTCTGTAATGCTGTGTGTGTTGAAAGTGTCGTCAA";
                    sequences[1] = "ATGCTGTGTGTCTGTAATGCTGTGTGTCTGTAATGCTGTGTGTGTTGAAAGTGTCGTCAA";
                    sequences[2] = "ATGCTGTGTGTCTGTAATGCTGTGTGTCTGTAATGCTGTGTGTGTTGAAAGTGTCGTCAA";
                    sequences[3] = "ATGCTGTGTGTCTGTAATGCTGTGTGTCTGTAATGCTGTGTGTGTTGAAAGTGTCGTCAA";
                    sequences[4] = "ATGCTGTGTGTCTGTAATGCTGTGTGTCTGTAATGCTGTGTGTGTTGAAAGTGTCGTCAA";
                    sequences[5] = "ATGCTGTGTGTCTGTAATGCTGTGTGTCTGTAATGCTGTGTGTGTTGAAAGTGTCGTCAA";

                    REQUIRE ( utils::checkAlignmentType ( sequenNumber, residNumber, sequences ) == SequenceTypes::DNA );
                }
            }

            WHEN ( "Degenerated" ) {
                GIVEN ( "With gaps" ) {
                    sequences[0] = "ATGCTGTGTGTCTGTA--------------------------TGTTGAAAGTGTCGTCAA";
                    sequences[1] = "ATGCTGTGTGTCTGTA--------------------------TGTTGAAAGTGTCGTCAA";
                    sequences[2] = "ATGCTGTGTGTCTGTA--------------------------TGTTGAAAGTGTCGTCAA";
                    sequences[3] = "YKVHKVKVKVKHKVKY--------------------------KVKKVYYYVKVKHVKHYY";
                    sequences[4] = "YKVHKVKVKVKHKVKY--------------------------KVKKVYYYVKVKHVKHYY";
                    sequences[5] = "YKVHKVKVKVKHKVKY--------------------------KVKKVYYYVKVKHVKHYY";

                    REQUIRE ( utils::checkAlignmentType ( sequenNumber, residNumber, sequences ) == ( SequenceTypes::DNA | SequenceTypes::DEG ) );
                }

                GIVEN ( "Without gaps" ) {
                    sequences[0] = "ATGCTGTGTGTCTGTAATGCTGTGTGTCTGTAATGCTGTGTGTGTTGAAAGTGTCGTCAA";
                    sequences[1] = "ATGCTGTGTGTCTGTAATGCTGTGTGTCTGTAATGCTGTGTGTGTTGAAAGTGTCGTCAA";
                    sequences[2] = "ATGCTGTGTGTCTGTAATGCTGTGTGTCTGTAATGCTGTGTGTGTTGAAAGTGTCGTCAA";
                    sequences[3] = "YKVHKVKVKVKHKVKYYKVHKVKVKVKHKVKYYKVHKVKVKVKVKKVYYYVKVKHVKHYY";
                    sequences[4] = "YKVHKVKVKVKHKVKYYKVHKVKVKVKHKVKYYKVHKVKVKVKVKKVYYYVKVKHVKHYY";
                    sequences[5] = "YKVHKVKVKVKHKVKYYKVHKVKVKVKHKVKYYKVHKVKVKVKVKKVYYYVKVKHVKHYY";

                    REQUIRE ( utils::checkAlignmentType ( sequenNumber, residNumber, sequences ) == ( SequenceTypes::DNA | SequenceTypes::DEG ) );
                }
            }
        }

        GIVEN ( "RNA Sequences" ) {
            WHEN ( "Pure" ) {
                GIVEN ( "With gaps" ) {
                    sequences[0] = "AUGCUGUGUGUCUGUA--------------------------UGUUGAAAGUGUCGUCAA";
                    sequences[1] = "AUGCUGUGUGUCUGUA--------------------------UGUUGAAAGUGUCGUCAA";
                    sequences[2] = "AUGCUGUGUGUCUGUA--------------------------UGUUGAAAGUGUCGUCAA";
                    sequences[3] = "AUGCUGUGUGUCUGUA--------------------------UGUUGAAAGUGUCGUCAA";
                    sequences[4] = "AUGCUGUGUGUCUGUA--------------------------UGUUGAAAGUGUCGUCAA";
                    sequences[5] = "AUGCUGUGUGUCUGUA--------------------------UGUUGAAAGUGUCGUCAA";

                    REQUIRE ( utils::checkAlignmentType ( sequenNumber, residNumber, sequences ) == SequenceTypes::RNA );
                }

                GIVEN ( "Without gaps" ) {
                    sequences[0] = "AUGCUGUGUGUCUGUAAUGCUGUGUGUCUGUAAUGCUGUGUGUGUUGAAAGUGUCGUCAA";
                    sequences[1] = "AUGCUGUGUGUCUGUAAUGCUGUGUGUCUGUAAUGCUGUGUGUGUUGAAAGUGUCGUCAA";
                    sequences[2] = "AUGCUGUGUGUCUGUAAUGCUGUGUGUCUGUAAUGCUGUGUGUGUUGAAAGUGUCGUCAA";
                    sequences[3] = "AUGCUGUGUGUCUGUAAUGCUGUGUGUCUGUAAUGCUGUGUGUGUUGAAAGUGUCGUCAA";
                    sequences[4] = "AUGCUGUGUGUCUGUAAUGCUGUGUGUCUGUAAUGCUGUGUGUGUUGAAAGUGUCGUCAA";
                    sequences[5] = "AUGCUGUGUGUCUGUAAUGCUGUGUGUCUGUAAUGCUGUGUGUGUUGAAAGUGUCGUCAA";

                    REQUIRE ( utils::checkAlignmentType ( sequenNumber, residNumber, sequences ) == SequenceTypes::RNA );
                }
            }

            WHEN ( "Degenerated" ) {
                GIVEN ( "With gaps" ) {
                    sequences[0] = "AUGCUGUGUGUCUGUA--------------------------UGUUGAAAGUGUCGUCAA";
                    sequences[1] = "AUGCUGUGUGUCUGUA--------------------------UGUUGAAAGUGUCGUCAA";
                    sequences[2] = "AUGCUGUGUGUCUGUA--------------------------UGUUGAAAGUGUCGUCAA";
                    sequences[3] = "YKVHKVKVKVKHKVKY--------------------------KVKKVYYYVKVKHVKHYY";
                    sequences[4] = "YKVHKVKVKVKHKVKY--------------------------KVKKVYYYVKVKHVKHYY";
                    sequences[5] = "YKVHKVKVKVKHKVKY--------------------------KVKKVYYYVKVKHVKHYY";

                    REQUIRE ( utils::checkAlignmentType ( sequenNumber, residNumber, sequences ) == ( SequenceTypes::RNA | SequenceTypes::DEG ) );
                }

                GIVEN ( "Without gaps" ) {
                    sequences[0] = "AUGCUGUGUGUCUGUAAUGCUGUGUGUCUGUAAUGCUGUGUGUGUUGAAAGUGUCGUCAA";
                    sequences[1] = "AUGCUGUGUGUCUGUAAUGCUGUGUGUCUGUAAUGCUGUGUGUGUUGAAAGUGUCGUCAA";
                    sequences[2] = "AUGCUGUGUGUCUGUAAUGCUGUGUGUCUGUAAUGCUGUGUGUGUUGAAAGUGUCGUCAA";
                    sequences[3] = "YKVHKVKVKVKHKVKYYKVHKVKVKVKHKVKYYKVHKVKVKVKVKKVYYYVKVKHVKHYY";
                    sequences[4] = "YKVHKVKVKVKHKVKYYKVHKVKVKVKHKVKYYKVHKVKVKVKVKKVYYYVKVKHVKHYY";
                    sequences[5] = "YKVHKVKVKVKHKVKYYKVHKVKVKVKHKVKYYKVHKVKVKVKVKKVYYYVKVKHVKHYY";

                    REQUIRE ( utils::checkAlignmentType ( sequenNumber, residNumber, sequences ) == ( SequenceTypes::RNA | SequenceTypes::DEG ) );
                }
            }
        }

        GIVEN ( "Fully degenerated sequences" ) {
            GIVEN ( "With gaps" ) {
                sequences[0] = "YKVHKVKVKVKHKVKY--------------------------KVKKVYYYVKVKHVKHYY";
                sequences[1] = "YKVHKVKVKVKHKVKY--------------------------KVKKVYYYVKVKHVKHYY";
                sequences[2] = "YKVHKVKVKVKHKVKY--------------------------KVKKVYYYVKVKHVKHYY";
                sequences[3] = "YKVHKVKVKVKHKVKY--------------------------KVKKVYYYVKVKHVKHYY";
                sequences[4] = "YKVHKVKVKVKHKVKY--------------------------KVKKVYYYVKVKHVKHYY";
                sequences[5] = "YKVHKVKVKVKHKVKY--------------------------KVKKVYYYVKVKHVKHYY";

                REQUIRE ( utils::checkAlignmentType ( sequenNumber, residNumber, sequences ) == ( SequenceTypes::DNA | SequenceTypes::DEG ) );
            }
            GIVEN ( "Without gaps" ) {
                sequences[0] = "YKVHKVKVKVKHKVKYYKVHKVKVKVKHKVKYYKVHKVKVKVKVKKVYYYVKVKHVKHYY";
                sequences[1] = "YKVHKVKVKVKHKVKYYKVHKVKVKVKHKVKYYKVHKVKVKVKVKKVYYYVKVKHVKHYY";
                sequences[2] = "YKVHKVKVKVKHKVKYYKVHKVKVKVKHKVKYYKVHKVKVKVKVKKVYYYVKVKHVKHYY";
                sequences[3] = "YKVHKVKVKVKHKVKYYKVHKVKVKVKHKVKYYKVHKVKVKVKVKKVYYYVKVKHVKHYY";
                sequences[4] = "YKVHKVKVKVKHKVKYYKVHKVKVKVKHKVKYYKVHKVKVKVKVKKVYYYVKVKHVKHYY";
                sequences[5] = "YKVHKVKVKVKHKVKYYKVHKVKVKVKHKVKYYKVHKVKVKVKVKKVYYYVKVKHVKHYY";

                REQUIRE ( utils::checkAlignmentType ( sequenNumber, residNumber, sequences ) == ( SequenceTypes::DNA | SequenceTypes::DEG ) );
            }
        }

        delete [] sequences;
    }

    WHEN ( "CHECK COLORS" ) {
        WARN ( "Test not implemented" );
    }
}

SCENARIO ( "Gap Ranges", "[utils][stats][gap]" ) {

    struct testStruct {
        
        int sequenNumber;
        int gapNumber;
        int expectedResult;

        testStruct ( int sequences, int gaps, int result ) : sequenNumber ( sequences ), gapNumber ( gaps ), expectedResult ( result ) { }
    };

    std::vector<testStruct> testValues = {
        testStruct ( 100, 0, 11 ),
        testStruct ( 100, 25, 10 ),
        testStruct ( 100, 50, 9 ),
        testStruct ( 100, 65, 8 ),
        testStruct ( 100, 75, 7 ),
        testStruct ( 100, 80, 6 ),
        testStruct ( 100, 85, 5 ),
        testStruct ( 100, 90, 4 ),
        testStruct ( 100, 95, 3 ),
        testStruct ( 100, 99, 2 ),

        testStruct ( 100, 100, 0 ),

        testStruct ( 1000, 0, 11 ),
        testStruct ( 1000, 250, 10 ),
        testStruct ( 1000, 500, 9 ),
        testStruct ( 1000, 650, 8 ),
        testStruct ( 1000, 750, 7 ),
        testStruct ( 1000, 800, 6 ),
        testStruct ( 1000, 850, 5 ),
        testStruct ( 1000, 900, 4 ),
        testStruct ( 1000, 950, 3 ),
        testStruct ( 1000, 990, 2 ),
        testStruct ( 1000, 999, 1 ),
        testStruct ( 1000, 1000, 0 ),
    };

    for ( testStruct & testValue : testValues ) {
        GIVEN ( std::to_string ( testValue.sequenNumber ) + " sequences; " + std::to_string ( testValue.gapNumber ) + " gaps" ) {

            WHEN ( "Direct Method" ) {
                REQUIRE ( utils::GetGapStep ( &testValue.gapNumber, testValue.sequenNumber ) == testValue.expectedResult );
            }
            WHEN ( "Precalculated Method" ) {
                REQUIRE ( utils::GetGapStep ( &testValue.gapNumber, 1.F / testValue.sequenNumber ) == testValue.expectedResult );
            }
        }
    }

}

SCENARIO ( "Sim Ranges", "[utils][stats][sim]" ) {

    std::vector<std::pair<float, int>> testValues = {
        std::pair<float, int> ( 0.F, 11 ),
        std::pair<float, int> ( 1.F, 0 ),
        std::pair<float, int> ( .750F, 10 ),
        std::pair<float, int> ( .500F, 9 ),
        std::pair<float, int> ( .250F, 8 ),
        std::pair<float, int> ( .100F, 7 ),
        std::pair<float, int> ( .010F, 6 ),
        std::pair<float, int> ( .001F, 5 ),
        std::pair<float, int> ( 1e-4F, 4 ),
        std::pair<float, int> ( 1e-5F, 3 ),
        std::pair<float, int> ( 1e-6F, 2 ),
        std::pair<float, int> ( 1e-7F, 1 )
    };

    for ( std::pair<float, int> & pair : testValues ) {
        GIVEN ( std::to_string ( pair.first ) ) {
            REQUIRE ( utils::GetSimStep ( &pair.first ) == pair.second );
        }
    }
}

SCENARIO ( "Cons Ranges", "[utils][stats][cons]" ) {

    std::vector<std::pair<float, int>> testValues = {
        std::pair<float, int> ( 0.F, 0 ),
        std::pair<float, int> ( 1.F, 11 ),
        std::pair<float, int> ( .750F, 10 ),
        std::pair<float, int> ( .500F, 9 ),
        std::pair<float, int> ( .350F, 8 ),
        std::pair<float, int> ( .250F, 7 ),
        std::pair<float, int> ( .200F, 6 ),
        std::pair<float, int> ( .150F, 5 ),
        std::pair<float, int> ( .100F, 4 ),
        std::pair<float, int> ( .050F, 3 ),
        std::pair<float, int> ( .001F, 2 ),
        std::pair<float, int> ( .0001F, 1 )
    };

    for ( std::pair<float, int> & pair : testValues ) {
        GIVEN ( std::to_string ( pair.first ) ) {
            INFO( std::setprecision(20) << pair.first );
            REQUIRE ( utils::GetConsStep ( &pair.first ) == pair.second );
        }
    }
}
