#include "../include/Catch/catch.hpp"
#include <algorithm>
#include <iterator>

#include <limits>
#include <math.h>

// // The matcher class
template <class T>
class ArrayComparer : public Catch::MatcherBase<T *> {
    T * m_begin;
    int size;
public:
    ArrayComparer ( T * begin, int size ) : m_begin( begin ), size( size ) {}

    // Performs the test for this matcher
    virtual bool match( T * i ) const override {
        
        return std::equal(i, i + size, m_begin);
    }

    // Produces a string describing what this matcher does. It should
    // include any provided data (the begin/ end in this case) and
    // be written as if it were stating a fact (in the output it will be
    // preceded by the value under test).
    virtual std::string describe() const override {
        std::ostringstream ss;
        ss << "content's is equal to " << m_begin << std::endl ;
        return ss.str();
    }
};

template<> 
inline bool ArrayComparer<double>::match(double * i) const
{
    for (int x = 0; x < size; x++)
    {
        double I = i[x], M = m_begin[x];
        if (fabs(i[x] - m_begin[x]) > std::numeric_limits<double>::epsilon()) 
            return false;
    }
    return true;
}

template<> 
inline bool ArrayComparer<float>::match(float * i) const
{
    for (int x = 0; x < size; x++)
    {
        float I = i[x], M = m_begin[x];
        if (fabs(i[x] - m_begin[x]) > std::numeric_limits<float>::epsilon()) 
            return false;
    }
    return true;
}

// The builder function
template <class T>
inline ArrayComparer<T> ArrayContentsEqual ( T * begin, int size ) {
    return ArrayComparer<T> ( begin, size );
}
