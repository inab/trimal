#ifndef READWRITES_H
#define READWRITES_H

#include "../../include/newAlignment.h"
#include <fstream>


class ReadWriteMS;
/**
 \brief Abstract class to implement.\n
        This class implementations are instantiated by the ReadWrite Machine State.\n
        Each implementation is related to an alignment format, and should know how to:\n
        - Recognize if a file is in its format\n
        - Load a file in its format and return an alignment object\n
        - Save an alignment to a file in this format\n
        - Recognize its acronyms.
 */
class ReadWriteBaseState
{
public:
    
    /** \brief Bool tag that indicates wether or not this state can load an alignment in the corresponding format.*/
    bool canLoad = false;
    
    /** \brief Bool tag that indicates wether or not this state can save an alignment in the corresponding format.*/
    bool canSave = false;
    
    /**
     \brief String that contains a well known acronym to the format
     */
    string name;
    
    /**
     \brief String that contains the main extension of the format
     */
    string extension;
    
    /**
     \brief Function to check if a file is in the implemented format.
     \param origin File Handler to the file.
     \return A number that represents the confidence recognizing this format.
     \attention Take in mind incompatibilities with other formats when you implement this function:\n\n
                
                Example:
                
                Fasta may recognize a file if it starts with '>', returning a 1.
                Pir files also begin with '>', but they also recognize the ';' after the '>' and 2 letter symbol.
                When the Pir recognizes its format, it returns a 2.
                So, when the Pir format is recognized, it's preferred over the Fasta Format.
                            
     */
    virtual int CheckAlignment(istream* origin) = 0;
    
    /**
     \brief Function to load a file in the current format and return an alignment object.
     \param filename Filename of the file to load.
     \return <b> Alignment Object</b> with the information of the file. <b> nullptr</b> if there was any error. 
     */
    virtual newAlignment* LoadAlignment(string filename) = 0;
    /**
     \brief Function to save to a file in the current format an alignment object.
     \param alignment Alignment to save.
     \param output File Handler where to save the formatted alignment;
     \param FileName File Name of the File Handler.
     \return <b> True </b>if file could be saved. <b> False </b> otherwise.
     \todo Remove the File Name Argument as it's stored on the alignment itself.
     */
    virtual bool SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName) = 0;
    
    /**
     \brief Function that recognizes acronyms of the format.
     \param FormatName acronym to recognize
     \return <b>True</b> if recognized, <b>False</b> otherwise.
     */
    virtual bool RecognizeOutputFormat(std::string FormatName) {
        return (name == FormatName);
    }
    
    /**
     \brief Class Destructor
     */
    virtual ~ReadWriteBaseState() {};
    
protected:
    /**
     * \brief Pointer to the object that contains and manages this state\n
     */
    ReadWriteMS* Machine;
};

#endif // READWRITES_H
