//
// Created by vfernandez on 08/03/19.
//

#ifndef TRIMAL_CHANGEPATHTOROOT_H
#define TRIMAL_CHANGEPATHTOROOT_H

#include "experimental/filesystem"
#include <iostream>

namespace catch2Utils {

    void init() {
        static bool todo = true;
        if (todo) {
            auto curPath = std::experimental::filesystem::current_path();

            while (curPath.string().find("cmake-build-") != std::string::npos)
                curPath = curPath.parent_path();

            if (std::experimental::filesystem::exists(curPath.append("bin"))) {
                std::experimental::filesystem::current_path(curPath);
                std::cout << "Current path: " << curPath << std::endl;
                std::flush(std::cout);
            } else {
                std::cout << "Failed to stablish a correct path for testing the module";
                std::flush(std::cout);
                exit(1);
            }
            todo = false;
        }
    }

}

#endif //TRIMAL_CHANGEPATHTOROOT_H
