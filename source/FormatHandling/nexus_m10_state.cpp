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

#include "FormatHandling/nexus_m10_state.h"

#include "FormatHandling/FormatManager.h"
#include "utils.h"

namespace FormatHandling {
int nexus_m10_state::CheckAlignment(std::istream *origin) {
    return 0;
}

Alignment *nexus_m10_state::LoadAlignment(const std::string &filename) {
    return nullptr;
}

bool nexus_m10_state::SaveAlignment(const Alignment &alignment, std::ostream *output) {
    /* Generate output alignment in NEXUS format setting only alignment block */

    int i, j, k, l, maxLongName = 0;
    std::string *tmpMatrix;

    // Check whether sequences in the alignment are aligned or not.
    // Warn about it if there are not aligned.
    if (!alignment.isAligned) {
        debug.report(ErrorCode::UnalignedAlignmentToAlignedFormat, new std::string[1]{this->name});
        return false;
    }

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    if (Machine->reverse)
    {
        /* Allocate local memory for generating output alignment */
        tmpMatrix = new std::string[alignment.originalNumberOfSequences];
        for(i = 0; i < alignment.originalNumberOfSequences; i++)
            tmpMatrix[i] = utils::getReverse(alignment.sequences[i]);
    }
    else tmpMatrix = alignment.sequences;

    // Compute maximum sequences name length
    for (i = 0; (i < alignment.originalNumberOfSequences); i++)
        if (alignment.saveSequences[i] != -1)
            maxLongName = utils::max(maxLongName, alignment.seqsName[i].size());

    if (maxLongName > PHYLIPDISTANCE) {
        maxLongName = PHYLIPDISTANCE;
        debug.report(WarningCode::HeaderWillBeCut, new std::string[1]{std::string(name)});
    }

    // Compute output file datatype
    alignment.getAlignmentType();

    std::string alignmentInfo {alignment.alignmentInfo};
    // Remove characters like ";" from input alignment information line
    while ((int) alignment.alignmentInfo.find(';') != (int) std::string::npos)
        alignmentInfo.erase(alignment.alignmentInfo.find(';'), 1);

    // Print Alignment header
    *output << "#NEXUS\nBEGIN DATA;\n DIMENSIONS NTAX="
            << alignment.numberOfSequences << " NCHAR=" << alignment.numberOfResidues << ";\n";

    // Print alignment datatype
    if (alignment.getAlignmentType() & SequenceTypes::DNA)
        *output << "FORMAT DATATYPE=DNA INTERLEAVE=yes GAP=-";
    else if (alignment.getAlignmentType() & SequenceTypes::RNA)
        *output << "FORMAT DATATYPE=RNA INTERLEAVE=yes GAP=-";
    else if (alignment.getAlignmentType() & SequenceTypes::AA)
        *output << "FORMAT DATATYPE=PROTEIN INTERLEAVE=yes GAP=-";

    i = 0;
    // Using information from input alignment. Use only some tags.
    while ((j = alignmentInfo.find(' ', i)) != (int) std::string::npos) {

        if ((alignmentInfo.substr(i, j - i)).compare(0, 7, "MISSING") == 0 ||
            (alignmentInfo.substr(i, j)).compare(0, 7, "missing") == 0)
            *output << " " << (alignmentInfo.substr(i, j - i));

        else if ((alignmentInfo.substr(i, j)).compare(0, 9, "MATCHCHAR") == 0 ||
                 (alignmentInfo.substr(i, j)).compare(0, 9, "matchchar") == 0)
            *output << " " << (alignmentInfo.substr(i, j - i));

        i = j + 1;
    }
    *output << ";\n";

    // Add a header indicating the number of residues of each sequence.
    for (i = 0; i < alignment.originalNumberOfSequences; i++) {
        if (alignment.saveSequences[i] == -1)
            continue;

        *output << "[Name: "
                << std::setw(maxLongName + 4) << std::left << alignment.seqsName[i].substr(0, maxLongName)
                << "Len: " << alignment.numberOfResidues << "]\n";
    }
    *output << "\nMATRIX";

    // Start filling the file with sequence names and sequences.
    for (j = 0, k = 0; j < alignment.originalNumberOfResidues;) {
        // Move until next not rejected residue
        if (alignment.saveResidues[j] == -1) {
            j++;
            continue;
        }

        // Iterate over the sequences
        for (i = 0; i < alignment.originalNumberOfSequences; i++) {
            // Skip rejected sequences
            if (alignment.saveSequences[i] == -1) continue;
            // Output sequence name
            *output << "\n" << std::setw(maxLongName + 5) << std::left << alignment.seqsName[i].substr(0, maxLongName);
            // Add residues per block.
            // k = residue position;
            // l = residues added on current line
            for (k = j, l = 0; k < alignment.originalNumberOfResidues && l < 50;) {
                // Don't save residues that have been rejected
                if (alignment.saveResidues[k] == -1) {
                    k++;
                    continue;
                }
                // Add the residue at position K
                *output << alignment.sequences[i][k];
                k++;
                l++;
                // If a block of 10 residues have been added, append a space to sequence
                if (l % 10 == 0 && l != 50) *output << " ";
            }
        }
        // Add the line to split blocks
        *output << "\n";
        j = k;
    }

    // End of MATRIX
    *output << "\n;\nEND;\n";

    /* Deallocate local memory */
    if (Machine->reverse)
        delete [] tmpMatrix;

    return true;
}

bool nexus_m10_state::RecognizeOutputFormat(const std::string &FormatName) {
    if (BaseFormatHandler::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "nexus_m10";
}

}