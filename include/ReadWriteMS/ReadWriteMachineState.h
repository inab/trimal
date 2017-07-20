#ifndef READWRITEMS_H
#define READWRITEMS_H

#include <vector>
#include <string>
#include "../../include/newAlignment.h"

class ReadWriteBaseState;

/**
 \brief Class that handles the load and save functions of an alignment.
 */
class ReadWriteMS
{
public:
    ReadWriteMS();
    ~ReadWriteMS();
    
private:
    /**
     \brief Vector that contains the available formats to load/save from.
            They are loaded into the format in the constructor function.
     */
    std::vector<ReadWriteBaseState*> available_states;
    /**
     \brief Function that adds a newState to the available_states vector.
            This should be called on the constructor function foreach format existent.
     \param newState Pointer to the newState we want to instantiate.
     */
    void addState(ReadWriteBaseState* newState);
    
public:
    
    bool hasOutputFile  = true;
    bool shortNames     = false;
    bool keepHeader     = false;
    bool reverse        = false;
    
    // LEGACY PARAMETERS
    bool format         = false;
    bool type           = false;
    bool info           = false;
    
    /**
     \brief Function that creates an alignment from a filepath.
            It automatically detects the format of the file.
    \param inFile Filepath of the alignment to load.
    \return Pointer to the alignment if it could be loaded. Null if not.
     */
    newAlignment* loadAlignment(std::string inFile);
    
    /**
     \brief Function to save an alignment to a file.
            It searched on the available_states, the one that can write the alignment in the correct format.
     \param outFile Filepath to save the alignment.
     \param outFormat Format in which save the alignment.
    \param alignment Alignment
            */
    void saveAlignment(std::string outFile, std::string outFormat, newAlignment* alignment);
    
    /**
     \brief Function that takes multile files, loads them and saves in a cumulus of formats, using an outPattern.
     \param inFile Vector of files to load, reformat and save.
     \param outPattern Path and name of the new files. The function changes some tokens on the original string to obtain multiple versions:
                [in]        Token that is changed with the original filename without extension.
                [format]    Token that is changed with the new format name.
                [extension] Token that is changed with the format file extensions.
     \param outFormats Output formats that original files should reformat to. 
     */
    void processFile(std::vector< std::string >* inFile, std::string* outPattern, std::vector< std::string >* outFormats);
    
    /**
     \brief Function to obtain the format name of a given file.
     \param inFile Filepath of the file which we want to obtain it's format.
     */
    std::string getInputStateName(std::string inFile);
    
    /**
     \brief Function to obtain all format names available to this object.
     \todo This function doesn't take in mind half-states: States that can only load or save an alignment, as, for example, HTML format.
            This could be implemented with some tags on the classes, and this function may be replaced with three new: getAllFormatsAvailable, getInputFormats, getOutputFormats.
     */
    std::string getFormatsAvailable();
};


#endif // READWRITEMS_H
