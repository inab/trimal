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

#include "reportsystem.h"

const std::map<InfoCode, const char *> reporting::reportManager::InfoMessages =
{
    {InfoCode::CuttingSequence,
            "Cutting sequence \"[tag]\" at first appearance of stop codon \"[tag]\""
                    " (residue \"[tag]\") at position [tag] (length: [tag] \")"},

    {InfoCode::WindowSizeCompareset,
            "Window size (-w) is provided. "
                    "It's recommended to use specific consistency window size (-cw)"
                    " when using -compareset option"},
    {InfoCode ::AddingSNP,
            "Applying SNP to \"[tag]\":\"[tag]\" at position [tag] \"[tag]\"->\"[tag]\""
    },

    {InfoCode::RemovingDuplicateSequences,
                "Removing sequence \"[tag]\" as it is a duplicate of \"[tag]\"."}
};