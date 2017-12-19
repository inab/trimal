#include "../include/Catch/catch.hpp"
#include <algorithm>
#include <cmath>
#include <iterator>

#include <limits>
#include <cmath>
#include <fstream>

// // The matcher class
template <class T>
class ArrayComparer : public Catch::MatcherBase<T *> {
    T * object_2;
    int _arraySize;
    std::string _filename;
public:
    ArrayComparer ( T * begin, int size, std::string& filename)
            : object_2( begin ), _arraySize( size ), _filename(filename)
    {}

    // Performs the test for this matcher
    bool match( T * object_1 ) const override {
        bool result = std::equal(object_1, object_1 + _arraySize, object_2);
        if (!result && _filename != "")
        {
            std::ofstream ofile(_filename);
            for (int i = 0; i < _arraySize; i++)
            {
                ofile << object_1[i] << "\t" << object_2[i] << std::endl;
            }
        }
        return result;
    }

    // Produces a string describing what this matcher does. It should
    // include any provided data (the begin/ end in this case) and
    // be written as if it were stating a fact (in the output it will be
    // preceded by the value under test).
    std::string describe() const override {
        std::ostringstream ss;
        ss << "content's is equal to " << object_2 << std::endl ;
        return ss.str();
    }
};

template<> 
inline bool ArrayComparer<double>::match(double * i) const
{
    for (int x = 0; x < _arraySize; x++)
    {
        double I = i[x], M = object_2[x];
        if (fabs(I - M) > std::numeric_limits<double>::epsilon())
            return false;
    }
    return true;
}

template<> 
inline bool ArrayComparer<float>::match(float * i) const
{
    for (int x = 0; x < _arraySize; x++)
    {
        float I = i[x], M = object_2[x];
        if (std::fabs(I - M) > std::numeric_limits<float>::epsilon())
            return false;
    }
    return true;
}

// The builder function
template <class T>
inline ArrayComparer<T> ArrayContentsEqual ( T * begin, int size, std::string filename = "") {
    return ArrayComparer<T> ( begin, size, filename);
}
