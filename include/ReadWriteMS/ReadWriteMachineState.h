#ifndef READWRITEMS_H
#define READWRITEMS_H

#include <vector>
#include <string>
#include "../../include/newAlignment.h"

class ReadWriteBaseState;

/**
 \brief Class to handle \link ReadWriteBaseState ReadWriteBaseStates \endlink,
  which represents formats.\n
 It serves as a proxy to the formats, so the code outside the ReadWriteMachine
  is format-agnostic.
 */
class ReadWriteMS
{
public:
    ReadWriteMS();
    ~ReadWriteMS();
    
private:
    /**
     \brief Vector that contains the available formats to load/save from.\n
            They are loaded into the format in the constructor function.
     \note The formats (\link ReadWriteBaseState \endlink) are loaded statically,
      no dynamically, so when a new format is implemented,
      it has to be added <i>by hand</i>
     */
    std::vector<ReadWriteBaseState*> available_states;
    /**
     \brief Function that adds a newState to the
             \link ReadWriteMS::available_states \endlink vector.\n
            This should be called on the constructor function foreach format existent.
     \param newState Pointer to the newState we want to instantiate.
     */
    void addState(ReadWriteBaseState* newState);
    
public:
    
    /** \brief Tag to know if the machine
        has an output file or it has to output to console. */
    bool hasOutputFile  = true;
    /** \brief Tag to know if the machine should keep original headers.*/
    bool keepHeader     = false;
    /** \brief Tag to know if sequences should be reversed before saving them.*/
    bool reverse        = false;
    
    // LEGACY PARAMETERS
    /** \brief Tag to know if the machine should output
        the format information about the alignment.*/
    bool format         = false;
    /** \brief Tag to know if the machine should output
        the type of the alignment.*/
    bool type           = false;
    /** \brief Tag to know if the machine should output
        the information of the alignment.*/
    bool info           = false;
    
    /**
    \brief Function that creates an alignment given a file path.
        It automatically detects the format of the file.
    \param inFile File path of the alignment to load.
    \return <b>Pointer to the alignment</b> if it could be loaded.\n
            <b>Null</b> if not.
        */
    newAlignment* loadAlignment(std::string inFile);

    /**
    \brief Function to save an alignment to a file.
        It searches among the available_states one that can write the alignment
         in the specified format.
    \param outPattern File path to save the alignment.
    \param outFormats Format in which save the alignment.
    \param alignment Alignment
            */
    bool saveAlignment(std::string outPattern, std::vector< std::string >* outFormats, newAlignment* alignment);
    
    /**
     \brief Function that takes multiple files,
             loads them and saves in a cumulus of formats, using an outPattern.
     \param inFile Vector of files to load, reformat and save.
     \param outPattern Path and name of the new files.\n
               The function changes some optional tokens on the original string to obtain multiple versions: \n
                - <b> [in] </b>        Token that is changed with the original filename without extension.\n
                - <b> [format] </b>    Token that is changed with the new format name.\n
                - <b> [extension] </b> Token that is changed with the format file extensions.
     \note Without using tags, when outputting several files
     \param outFormats Output formats that original files should reformat to. 
     */
    void loadAndSaveMultipleAlignments(std::vector< std::string >* inFile, std::string* outPattern, std::vector< std::string >* outFormats);
    
    /**
     \brief Function to obtain the format name of a given file.
     \param inFile File path of the file which we want to obtain it's format.
     */
    std::string getFileFormatName(std::string inFile);
    
    /**
     \brief Function to obtain all format names available
      by this object that can \b load an alignment.
     */
    std::string getInputFormatsAvailable();
    
    /**
     \brief Function to obtain all format names available
      by this object that can \b save an alignment.
     */
    std::string getOutputFormatsAvailable();
    
    /**
     \brief Function to divide an alignment into different alignments,
             each one with a sequence from the original.\n
            This function does a deep copy of each sequence,
             so the original alignment can be deleted after being splitted*/
    std::vector<newAlignment*> splitAlignmentKeeping(newAlignment& alignment);
};


#endif // READWRITEMS_H
