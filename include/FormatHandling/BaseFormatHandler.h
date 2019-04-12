/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2009-2019
        Fernandez-Rodriguez V.  (victor.fernandez@bsc.es)
        Capella-Gutierrez S.    (salvador.capella@bsc.es)
        Gabaldon, T.            (tgabaldon@crg.es)

    This file is part of trimAl/readAl.

    trimAl/readAl are free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl/readAl are distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl/readAl. If not, see <http://www.gnu.org/licenses/>.

***************************************************************************** */

#ifndef READWRITES_H
#define READWRITES_H

#include "Alignment/Alignment.h"
#include "reportsystem.h"

#include <iomanip>

namespace FormatHandling {
class FormatManager;

/**
 \brief Abstract class that serves as template for Format Handlers.
        Formats must inherit from this class and
         will be handled by \link FormatManager \endlink

        Each class should know how to:
        - Recognize if a file is in its format\n
        - Load a file in its format and return an alignment object\n
        - Save an alignment to a file in this format\n
        - Recognize its format acronyms.
 */
class BaseFormatHandler {
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
    virtual Alignment *LoadAlignment(const std::string &filename) = 0;

    /**
     \brief Function to save a \link Alignment \endlink to a file.
     \param alignment Alignment to save.
     \param output File Handler where to save the formatted alignment;
     \return    <b> True </b>if file could be saved.\n
                <b> False </b> otherwise.
     */
    virtual bool SaveAlignment(const Alignment &alignment, std::ostream *output) = 0;

    /**
     \brief Function that recognizes acronyms of the format.
     \param FormatName acronym to recognize
     \return <b>True</b> if recognized\n
             <b>False</b> otherwise.
     */
    virtual bool RecognizeOutputFormat(const std::string &FormatName) {
        return (name == FormatName);
    }

    /**
     \brief Class Destructor
     */
    virtual ~BaseFormatHandler() = default;

protected:
    /**
     \brief Pointer to the \link FormatManager \endlink
     that contains and manages this state
     */
    FormatHandling::FormatManager *Machine;
};
}

#endif // READWRITES_H
