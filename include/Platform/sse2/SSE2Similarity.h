/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2022
        Larralde M.             (martin.larralde@embl.de)

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

#ifndef TRIMAL_SSE2SIMILARITY_H
#define TRIMAL_SSE2SIMILARITY_H

#ifdef HAVE_SSE2

#include "Statistics/Similarity.h"

namespace statistics {
    class SSE2Similarity: public Similarity {
    // private:
    //     std::vector<char> ascii_vhash;
    //     std::vector<char> colgap;
    //     std::string column;
    public:
        SSE2Similarity(Alignment* parentAlignment);
        void calculateMatrixIdentity() override;
        bool calculateVectors(bool cutByGap) override;
    };
}

#endif
#endif
