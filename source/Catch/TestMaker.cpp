// Target to make tests for each file in dataset.
// This allows to have a full integration of the tests,
//      as they will be hardcoded in files, while maintaining dynamism

#include <glob.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <cstring>
#include <utils.h>

// Method to make a ls on a directory
// https://stackoverflow.com/a/24703135
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

int main(int argc, char *argv[]) {

    std::ifstream moldFile("./source/Catch/mold_test_Alignment.cpp");
    std::ofstream testFile("./source/Catch/test_generated_alignment.cpp");

    std::string mold = std::string(std::istreambuf_iterator<char>(moldFile), std::istreambuf_iterator<char>());

    unsigned long scenarioPos = mold.find("// EndOfHeader");
    scenarioPos = mold.find('.', scenarioPos) + 2 /* The point itself and the newline */;
    testFile << mold.substr(0, scenarioPos);
    mold = mold.substr(scenarioPos);

    // We obtain the path names of all files and directories of the data set folder
    std::vector<std::string> filenames = globVector("./dataset/");

    // Output File

    int counter = 0;
    // We check if the pathname is a file, directory, or something else
    // https://stackoverflow.com/a/146938
    for (const std::string &filename : filenames)
    {
        if (counter++ == 20) break;
        // Get the full path length
        size_t filenameSize = std::strlen("./dataset/") + std::strlen(&filename[0]);
        // region Store it in a new variable
        auto pathname = std::string();
        pathname += "./dataset/";
        pathname += filename;
        // endregion
        struct stat s{};
        if(stat(&pathname[0], &s) == 0 )
        {
            if( s.st_mode & S_IFREG ) // Is file?
            {
                std::cout << pathname << std::endl;
                testFile << utils::ReplaceString(mold, "input_filename", filename);
                testFile << "\n";
            }
        }
    }
}