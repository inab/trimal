#include "../../include/ReadWriteMS/nexus_m10_state.h"

#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/defines.h"
#include "../../include/utils.h"

int nexus_m10_state::CheckAlignment(std::istream *origin) {
    return 0;
}

newAlignment *nexus_m10_state::LoadAlignment(std::string& filename) {
    return nullptr;
}

bool nexus_m10_state::SaveAlignment(newAlignment *alignment, std::ostream *output, std::string *FileName) {
    /* Generate output alignment in NEXUS format setting only alignment block */

    int i, j, k, l, maxLongName = 0;
    std::string *tmpMatrix;

    // Check whether sequences in the alignment are aligned or not.
    // Warn about it if there are not aligned.
    if (!alignment->isAligned) {
        debug.report(ErrorCode::UnalignedAlignmentToAlignedFormat, new std::string[1]{this->name});
        return false;
    }

    // Allocate local memory for generating output alignment
    tmpMatrix = new std::string[alignment->originalSequenNumber];

    // Depending on alignment orientation: forward or reverse. Copy directly
    // sequence information or get firstly the reversed sequences and then
    // copy it into local memory
    for (i = 0; i < alignment->originalSequenNumber; i++)
        tmpMatrix[i] = (!Machine->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    // Compute maximum sequences name length
    for (i = 0; (i < alignment->originalSequenNumber); i++)
        if (alignment->saveSequences[i] != -1)
            maxLongName = utils::max(maxLongName, alignment->seqsName[i].size());

    maxLongName = std::min(maxLongName, PHYLIPDISTANCE);
    // Compute output file datatype
    alignment->getAlignmentType();

    // Remove characters like ";" from input alignment information line
    while ((int) alignment->aligInfo.find(';') != (int) std::string::npos)
        alignment->aligInfo.erase(alignment->aligInfo.find(';'), 1);

    // Print Alignment header
    *output << "#NEXUS\nBEGIN DATA;\n DIMENSIONS NTAX="
            << alignment->sequenNumber << " NCHAR=" << alignment->residNumber << ";\n";

    // Print alignment datatype
    if (alignment->getAlignmentType() & SequenceTypes::DNA)
        *output << "FORMAT DATATYPE=DNA INTERLEAVE=yes GAP=-";
    else if (alignment->getAlignmentType() & SequenceTypes::RNA)
        *output << "FORMAT DATATYPE=RNA INTERLEAVE=yes GAP=-";
    else if (alignment->getAlignmentType() & SequenceTypes::AA)
        *output << "FORMAT DATATYPE=PROTEIN INTERLEAVE=yes GAP=-";

    i = 0;
    // Using information from input alignment. Use only some tags.
    while ((j = alignment->aligInfo.find(' ', i)) != (int) std::string::npos) {

        if ((alignment->aligInfo.substr(i, j - i)).compare(0, 7, "MISSING") == 0 ||
            (alignment->aligInfo.substr(i, j)).compare(0, 7, "missing") == 0)
            *output << " " << (alignment->aligInfo.substr(i, j - i));

        else if ((alignment->aligInfo.substr(i, j)).compare(0, 9, "MATCHCHAR") == 0 ||
                 (alignment->aligInfo.substr(i, j)).compare(0, 9, "matchchar") == 0)
            *output << " " << (alignment->aligInfo.substr(i, j - i));

        i = j + 1;
    }
    *output << ";\n";

    // Add a header indicating the number of residues of each sequence.
    for (i = 0; i < alignment->originalSequenNumber; i++) {
        if (alignment->saveSequences[i] == -1)
            continue;

        *output << "[Name: "
                << std::setw(maxLongName + 4) << std::left << alignment->seqsName[i].substr(0, maxLongName)
                << "Len: " << alignment->residNumber << "]\n";
    }
    *output << "\nMATRIX";

    // Start filling the file with sequence names and sequences.
    for (j = 0, k = 0; j < alignment->originalResidNumber;) {
        // Move until next not rejected residue
        if (alignment->saveResidues[j] == -1) {
            j++;
            continue;
        }

        // Iterate over the sequences
        for (i = 0; i < alignment->originalSequenNumber; i++) {
            // Skip rejected sequences
            if (alignment->saveSequences[i] == -1) continue;
            // Output sequence name
            *output << "\n" << std::setw(maxLongName + 5) << std::left << alignment->seqsName[i].substr(0, maxLongName);
            // Add residues per block.
            // k = residue position;
            // l = residues added on current line
            for (k = j, l = 0; k < alignment->originalResidNumber && l < 50;) {
                // Don't save residues that have been rejected
                if (alignment->saveResidues[k] == -1) {
                    k++;
                    continue;
                }
                // Add the residue at position K
                *output << alignment->sequences[i][k];
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
    delete[] tmpMatrix;

    return true;
}

bool nexus_m10_state::RecognizeOutputFormat(std::string& FormatName) {
    if (ReadWriteBaseState::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "nexus_m10";
}

