#ifndef catch_helper
#define catch_helper
#include <string.h>

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

void populate ( std::string *& s, picojson::value value ) {
    picojson::array pArray = value.get<picojson::array>();
    s = new std::string[pArray.size()];
    int i = 0;
    for ( auto it = pArray.begin(); it != pArray.end(); it++, i++ ) {
        s[i] = it->to_str();
    }
}

#endif
