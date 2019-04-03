//
// Created by vfernandez on 06/03/19.
//

#include <functional>
#include <iostream>
#include <sstream>
#include <cstring>

#ifndef capturingL
#define capturingL
namespace catch2Utils {

    /**
     * @brief https://stackoverflow.com/a/4753399
     */
    class ScopedRedirect {
    public:
        ScopedRedirect(std::ostream &inOriginal, std::ostream &inRedirect) :
                mOriginal(inOriginal),
                mOldBuffer(inOriginal.rdbuf(inRedirect.rdbuf())) {}

        ~ScopedRedirect() { mOriginal.rdbuf(mOldBuffer); }

    private:
        ScopedRedirect(const ScopedRedirect &);

        ScopedRedirect &operator=(const ScopedRedirect &);

        std::ostream &mOriginal;
        std::streambuf *mOldBuffer;
    };

    std::string capturingLambda(std::function<void()> && function);

// Useful to keep coherence with Catch2 framework, but ignoring defines
//#define CAPTURE_LAMBDA(x, y) std::string x = catch2Utils::capturingLambda(y);
//#undef CAPTURE_LAMBDA
}

#endif
