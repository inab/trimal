//
// Created by vfernandez on 06/03/19.
//

#include <functional>
#include <iostream>
#include <sstream>
#include <cstring>
#include "ScopedRedirect.h"

namespace catch2Utils {

    std::string capturingLambda(std::function<void()> && function)
    {
        std::stringstream ss("");

        ScopedRedirect redirect1(std::cout , ss);
        ScopedRedirect redirect2(std::cerr , ss);

        function();

        if (ss.str().empty())
            return "-- No output captured --";
        else
        {
            std::string str = ss.str();
            if (!str.empty())
            {
                if (str[str.length()-1] == '\n')
                    str.erase(str.length()-1);
            }
            return  "-----  Captured output -----\n" +
                    str +
                    "\n  -- End of captured output --";
        }
    }
}

