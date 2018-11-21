#ifndef READWRITES_H
#define READWRITES_H

#include "../../include/newAlignment.h"
#include "../../include/reportsystem.h"

#include <iomanip>

class ReadWriteMS;

/**
 \brief Abstract class that represent a format.
        Formats must inherit from this class and
         will be handled by \link ReadWriteMS \endlink

        Each class should know how to:
        - Recognize if a file is in its format\n
        - Load a file in its format and return an alignment object\n
        - Save an alignment to a file in this format\n
        - Recognize its format acronyms.
 */
class ReadWriteBaseState {
public:

    /** \brief Bool tag that indicates whether or not this state
     * can load an alignment in the corresponding format.*/
    bool canLoad = false;

    /** \brief Bool tag that indicates whether or not this state
     * can save an alignment in the corresponding format.*/
    bool canSave = false;

    /**
     \brief String that contains a well known acronym to the format
     */
    std::string name;

    /**
     \brief String that contains the main extension of the format
     */
    std::string extension;

    /**
     \brief Function to check if a file is in the implemented format.
     \param origin File Handler to the file.
     \return A number that represents the confidence recognizing this format.
     \attention Take in mind incompatibilities with
                other formats when you implement this function:\n\n
                <i>Fasta may recognize a file if it starts with '>', returning a 1.\n
                Pir files also begin with '>', but they also recognize the ';' after the '>' and 2 letter symbol.\n
                When the Pir recognizes its format, it returns a 2.\n
                So, when the Pir format is recognized, it's preferred over the Fasta Format,
                 when both of them can load the alignment bur pir extracts more information.</i>\n
     */
    virtual int CheckAlignment(std::istream *origin) = 0;

    /**
     \brief Function to load a file in the current format and return an alignment object.
     \param filename Filename of the file to load.
     \return    <b> Alignment</b> loaded with the information of the file. \n
                <b> nullptr</b> if there was any error.
     */
    virtual newAlignment *LoadAlignment(std::string& filename) = 0;

    /**
     \brief Function to save a \link newAlignment \endlink to a file.
     \param alignment Alignment to save.
     \param output File Handler where to save the formatted alignment;
     \param FileName File Name of the File Handler.
     \return    <b> True </b>if file could be saved.\n
                <b> False </b> otherwise.
     \todo Remove the File Name Argument as it's stored on the alignment itself.
     */
    virtual bool SaveAlignment(newAlignment *alignment, std::ostream *output, std::string *FileName) = 0;

    /**
     \brief Function that recognizes acronyms of the format.
     \param FormatName acronym to recognize
     \return <b>True</b> if recognized\n
             <b>False</b> otherwise.
     */
    virtual bool RecognizeOutputFormat(std::string& FormatName) {
        return (name == FormatName);
    }

    /**
     \brief Class Destructor
     */
    virtual ~ReadWriteBaseState() {};

protected:
    /**
     \brief Pointer to the \link ReadWriteMS \endlink
     that contains and manages this state
     */
    ReadWriteMS *Machine;
};

#endif // READWRITES_H
