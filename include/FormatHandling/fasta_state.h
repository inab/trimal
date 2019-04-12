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

#ifndef FASTASTATE_H
#define FASTASTATE_H

#include "BaseFormatHandler.h"

namespace FormatHandling {
class fasta_state : public BaseFormatHandler {
public:

    explicit fasta_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "fasta";
        extension = "fasta";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    Alignment *LoadAlignment(const std::string &filename) override;

    bool SaveAlignment(const Alignment &alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(const std::string &FormatName) override;

};
}
#endif // FASTASTATE_H
