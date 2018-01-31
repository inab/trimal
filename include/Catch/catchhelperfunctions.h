#ifndef catch_helper
#define catch_helper
#include <string.h>
#include <fstream>
#include "picojson.h"

void populate ( int *& s, picojson::value& value ) {
    picojson::array pArray = value.get<picojson::array>();
    s = new int[pArray.size()];
    for (int i = 0; i < pArray.size(); i++ ) {
        s[i] = static_cast<int>(pArray[i].get<double>());
    }
}

void populate ( float *& s, picojson::value& value ) {
    picojson::array pArray = value.get<picojson::array>();
    s = new float[pArray.size()];
    for (int i = 0; i < pArray.size(); i++ ) {
        s[i] = static_cast<float>(pArray[i].get<double>());
    }
}

void populate ( double *& s, picojson::value& value ) {
    picojson::array pArray = value.get<picojson::array>();
    s = new double[pArray.size()];
    for (int i = 0; i < pArray.size(); i++ ) {
        s[i] = pArray[i].get<double>();
    }
}

void populate ( std::string *& s, picojson::value& value ) {
    picojson::array pArray = value.get<picojson::array>();
    s = new std::string[pArray.size()];
    for (int i = 0; i < pArray.size(); i++ ) {
        s[i] = pArray[i].to_str();
    }
}

bool loadJSON ( const std::string & filename, picojson::value & alignmentTestData ) {

    // region Avoid reloading the same JSON again.
    static std::string lastFilename = "";
    static bool lastValue = false;
    if ( lastFilename == filename ) {
        return lastValue;
    }
    // endregion

    std::ifstream testData;
    testData.open ( "./dataset/testingFiles/alignmentIntermediates/" + filename + ".json" );
    if ( testData.is_open() ) {
        std::string content ( ( std::istreambuf_iterator<char> ( testData ) ),
                              ( std::istreambuf_iterator<char>() ) );

        testData.close();
        // parse the input
        picojson::parse ( alignmentTestData, content );
        const std::string &err = picojson::get_last_error();

        if ( ! err.empty() ) {
            WARN ( err );
            lastValue = false;
        } else if ( ! alignmentTestData.is<picojson::object>() ) {
            WARN ( "JSON is not an object" );
            lastValue = false;
        } else
            lastValue = true;
    } else
        lastValue = false;

    lastFilename = filename;
    return lastValue;
}

#endif
