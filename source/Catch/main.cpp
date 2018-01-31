// /*
//#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#define CATCH_CONFIG_RUNNER

#include "../../include/Catch/catch.hpp"

#include <glob.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <cstring>
#include <utils.h>

std::vector<std::string> globVector(const std::string& pattern){
    glob_t glob_result{};
    glob((pattern + "*").c_str(), GLOB_TILDE, nullptr, &glob_result);
    std::vector<std::string> files;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        files.emplace_back(std::string(glob_result.gl_pathv[i] + std::strlen(&pattern[0])));
    }
    globfree(&glob_result);
    return files;
}

void cleanResults()
{

    // We obtain the path names of all files and directories of the data set folder
    std::vector<std::string> filenames = globVector("./dataset/testingFiles/alignmentErrors/");

    int counter = 0;
    // We check if the pathname is a file, directory, or something else
    // https://stackoverflow.com/a/146938
    for (const std::string &filename : filenames)
    {
        // Get the full path length
        size_t filenameSize = std::strlen("./dataset/testingFiles/alignmentErrors/") + std::strlen(&filename[0]);
        // region Store it in a new variable
        auto pathname = std::string();
        pathname += "./dataset/testingFiles/alignmentErrors/";
        pathname += filename;
        // endregion
        struct stat s{};
        if(stat(&pathname[0], &s) == 0 )
        {
            if( s.st_mode & S_IFREG ) // Is file?
            {
                remove(&pathname[0]);
            }
        }
    }
}


int main(int argc, char *argv[]) {

    cleanResults();

    int result = Catch::Session().run(argc, argv);

    // global clean-up...

    return result;
}

