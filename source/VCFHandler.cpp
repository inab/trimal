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

#include "VCFHandler.h"

namespace ngs {

    namespace IUPAC {

        int getTagFromCharArray(
                const char *const array,
                const char separator)
        {
            if (array[1] != separator)
                return -1;
            // Loop variables
            size_t c, maxlen;
            // Current tag value
            int curval = 0;
            // Iterate over the sequence, bit-adding the tag
            for (c = 0, maxlen = strlen(array); c < maxlen; c++) {
                switch (array[c]) {
                    case 'A':
                        curval |= 1 << 0;
                        break; // = 1
                    case 'C':
                        curval |= 1 << 1;
                        break; // = 2
                    case 'T':
                        curval |= 1 << 2;
                        break; // = 4
                    case 'G':
                        curval |= 1 << 3;
                        break; // = 8
                    default:
                        break;
                }

                // Check if we haven't arrived to end of sequence
                if (++c < maxlen)
                {
                    // check the next element in array is the separator char
                    if (array[c] == separator)
                        continue;
                    // Otherwise, it's not a well-formed sequence to convert
                    else
                        return -1;
                }
            }
            // If we have arrived to the end of the sequence
            //      it was well formed and parsed.
            if (c == maxlen)
                return curval;
            // Otherwise, we have early stopped.
            return -1;
        }

        char getCharFromTag(const int tag) {
            switch (tag) {

                // One Base
                // case 1://A = 1
                // case 2://C = 2
                // case 4://T = 4
                // case 8://G = 8

                // Two Bases
                case 3:
                    return 'M';// aMine AC
                case 5:
                    return 'W'; // Weak AT
                case 9:
                    return 'R'; // puRine AG
                case 6:
                    return 'Y';// pYrimidine CT
                case 10:
                    return 'S';// Strong CG
                case 12:
                    return 'K';// Keto TG

                    // Three Bases
                case 14:
                    return 'B';// CTG not A
                case 13:
                    return 'D';// ATG not C
                case 11:
                    return 'V';// ACG not T
                case 7:
                    return 'H';// ACT not G

                    // All Four Bases
                case 15:
                    return 'N';// ACTG

                default:
                    return '-';
            }
        }


    }

    namespace __internal {

         void printApeek(AlignmentVector &sources) {
            for (Alignment *A : sources) {
                std::cout << A->seqsName[0] << std::endl;

                for (int X = 0; X < A->numberOfSequences; X++) {
                    std::cout << "\t>" << A->seqsName[X] << std::endl;
                    std::cout << "\t" << A->sequences[X].substr(0, 50) << std::endl;
                }
            }
        }

        inline bool checkContigsWithReference(
                const AlignmentVector &sources,
                const StringVector &contigs)
        {
            // Store the result value, as we want to report all errors
            //      and thus, not skipping on first error
            bool checkIn = true;
            // Iterate over the list of contigs
            for (const std::string & contig : contigs) {
                int i;
                // Iterate all over the sources
                for (i = 0; i < sources.size(); i++) {
                    // And early skip on first coincidence
                    if (!strcmp(&contig[0], &sources[i]->seqsName[0][0]))
                        break;
                }
                // If we haven't early skipped, then, contig doesn't have
                //      an alignment assigned, and thus, there is no
                //      sequence on the original alignment that correspond
                //      to the contig.
                if (i == sources.size()) {
                    debug.report(ErrorCode::NoReferenceSequenceForContig, &contig[0]);
                    checkIn = false;
                }
            }
            // Return the result value, which is false in case any contig
            //      doesn't have a reference sequence on the original alignment
            return checkIn;
        }

        inline  void increaseSequencesInAlignment(
                const AlignmentVector &sources,
                const StringVector &donors)
        {
            for (Alignment *nA : sources) {

                // Store the original information: Reference sequence and name
                std::string seq = nA->sequences[0];
                std::string nam = nA->seqsName[0];

                // Delete the original sequences
                delete[] nA->sequences;
                // Create a new sequences array
                nA->sequences = new std::string[donors.size() + 1];
                // Add the reference sequence to the new sequences array
                nA->sequences[0] = seq;

                // Delete the original sequence names
                delete[] nA->seqsName;
                // Create a new sequences names array
                nA->seqsName = new std::string[donors.size() + 1];
                // Add the reference sequence name to the new names array
                nA->seqsName[0] = nam;

                // For each donor, add a copy of
                //      the reference sequence and the donor name
                for (int i = 1; i < donors.size() + 1; i++) {
                    nA->sequences[i] = seq;
                    nA->seqsName[i] = donors[i - 1];
                }

                // Set the number of sequences to the number of donors + reference
                nA -> originalNumberOfSequences =
                nA -> numberOfSequences         = donors.size() + 1;

                // Initialize the alignment
                nA->fillMatrices(true, false);
            }
        }

        inline  bool extractFeature(
                const char *const line,
                vcfFeature &feature)
        {
            if (line[0] == '#')
                return false;

            std::string tmpLine{line};
            char *tmp;

            size_t strlenVal;

            feature.reset();

            // Contig
            tmp = std::strtok(&tmpLine[0], "\t");
            strlenVal = strlen(tmp);
            feature.contig = new char[strlenVal + 1];
            std::memmove(feature.contig, tmp, strlenVal);
            feature.contig[strlenVal] = '\0';

            // Position
            tmp = std::strtok(nullptr, "\t");
            feature.position = atoi(tmp) - 1;

            // ID
            std::strtok(nullptr, "\t");

            // Ref
            tmp = std::strtok(nullptr, "\t");
            strlenVal = strlen(tmp);
            feature.ref = new char[strlenVal + 1];
            std::memmove(feature.ref, tmp, strlenVal);
            feature.ref[strlenVal] = '\0';

            // Alt
            tmp = std::strtok(nullptr, "\t");
            strlenVal = strlen(tmp);
            feature.alt = new char[strlenVal + 1];
            std::memmove(feature.alt, tmp, strlenVal);
            feature.alt[strlenVal] = '\0';

            if (strlen(feature.alt) > 1)
            {
                int alt = ngs::IUPAC::getTagFromCharArray(feature.alt, ':');
                if (alt != -1)
                {
                    delete [] feature.alt;
                    feature.alt = new char[2] {ngs::IUPAC::getCharFromTag(alt), '\0'};
                }
            }

            // Quality
            tmp = std::strtok(nullptr, "\t");
            feature.quality = (float) atof(tmp);

            // Filter
            tmp = std::strtok(nullptr, "\t");
            feature.filter = std::strcmp(tmp, "PASS") == 0;

            // Info
            std::strtok(nullptr, "\t");

            // Format
            tmp = std::strtok(nullptr, "\t");
            char *format = new char[strlen(tmp) + 1];
            std::memmove(format, tmp, strlen(tmp) + 1);

            // Donors
            int counter = 0;
            tmp = std::strtok(nullptr, "\t");
            while (tmp != nullptr) {
                feature.donorsInfo.emplace_back(tmp);
                tmp = std::strtok(nullptr, "\t");
            }

            // Format -> Read Depth Index
            counter = 0;
            tmp = std::strtok(format, ":");
            feature.readDepthIndex = -1;
            while (tmp != nullptr) {
                if (strlen(tmp) > 1 && tmp[0] == 'D' && tmp[1] == 'P') {
                    feature.readDepthIndex = counter;
                    break;
                }
                tmp = std::strtok(nullptr, ":");
                counter++;
            }

            delete[] format;

            return true;
        }

        inline bool filterFeature(
                const vcfFeature &tmpVariant,
                const float minQuality)
        {
            if (tmpVariant.quality < minQuality)
                return false;

            if (strlen(tmpVariant.ref) != 1)
                return false;

            if (strlen(tmpVariant.alt) != 1)
                return false;

            return true;
        }

        inline  void applyVariantCallingFiles(
                const AlignmentVector &sources,
                const StringVector &filenames,
                const StringVector &contigs,
                const std::vector<std::vector<int>> &donorsPositions,
                const float minQuality,
                const float minCoverage,
                const bool ignoreFilter,
                const char *const replacementChar)
        {
            char *line = new char[4096];

            // Create a feature object to reuse
            vcfFeature feature = vcfFeature();

            for (int vcfIndex = 0; vcfIndex < filenames.size(); vcfIndex++) {

                // Open the VCF file
                std::ifstream infile(filenames[vcfIndex]);
                if (!infile.good()) {
                    debug.report(ErrorCode::CantOpenFile, &filenames[vcfIndex][0]);
                    continue;
                }

                // Iterate over all VCF
                while (infile.getline(line, 4096, '\n')) {

                    // Extract feature line by line
                    if (!extractFeature(line, feature))
                        // Continue if line is not feature
                        continue;


                    // Check if feature-referred contig is found

                    // contigIndex is also sourcesIndex
                    int contigIndex;
                    {
                        for (contigIndex = 0; contigIndex < contigs.size(); contigIndex++)
                        {
                            if (contigs[contigIndex] == feature.contig)
                                break;
                        }
                        if (contigIndex == contigs.size()) {
                            std::cout << "Contig not Found" << std::endl;
                            continue;
                        }
                    }

                    // Get the alignment that refers to the contig
                    Alignment * alig = sources[contigIndex];

                    // Cannot apply the feature if the position to change
                    //      is bigger than the alignment
                    if (alig->originalNumberOfResidues <= feature.position)
                    {
                        debug.report(ErrorCode::SNPoutOfBounds,
                            new std::string[3]
                            {
                                std::to_string(feature.position),
                                filenames[vcfIndex],
                                std::to_string(alig->originalNumberOfResidues)
                            }
                        );
                        continue;
                    }

                    // Filter-out features with low quality or not snp
                    if (!filterFeature(feature, minQuality))
                        continue;

                    // Iterate over all donors
                    for (int donorID = 0; donorID < feature.donorsInfo.size(); donorID++)
                    {
                        // Sanity check to prevent OutOfBounds
                        if (donorID >= donorsPositions[vcfIndex].size())
                        {
                            debug.report(ErrorCode::MoreDonorsOnLineThanPresented,
                                    new std::string[1]{ filenames[vcfIndex]} );
                            continue;
                        }
                        // Get the sequence position.
                        const int seqPos = donorsPositions[vcfIndex][donorID];

                        // Get the character & reference for further use
                        char & donorSequenceChar =
                                alig->sequences[seqPos][feature.position];
                        const char & referenceSequenceChar =
                                alig->sequences[0][feature.position];

                        // In case we filter - out by VCF filter, apply the replacement char
                        if (!feature.filter && !ignoreFilter)
                        {
                            if (replacementChar) {
                                alig->sequences[seqPos][feature.position] = *replacementChar;
                            }
                            continue;
                        }

                        // Get the string that informs the user whether
                        //      a feature should be applied to a donor
                        const std::string & referenceString = feature.donorsInfo[donorID];

                        // Only apply if the string is valid - different to 0
                        if (referenceString != "0") {
                            std::string donorInfo{referenceString};
                            // Extract the read depth
                            char *tmp = std::strtok(&donorInfo[0], ":");
                            for (int V = 0; V < feature.readDepthIndex; V++)
                                tmp = std::strtok(nullptr, ":");
                            int readDepth = std::stoi(tmp);

                            // Filter out by read depth
                            if (readDepth < minCoverage) {
                                // If provided, change the character by the replacement-char
                                //      as we don't know if the character is equal to the
                                //      reference or the provided SNP
                                if (replacementChar)
                                    donorSequenceChar = *replacementChar;
                                break;
                            }

                            // Check if the reference character is the same as expected in the sequence to apply
                            if (donorSequenceChar == std::tolower(feature.ref[0]) ||
                                donorSequenceChar == std::toupper(feature.ref[0]))
                            {
                                debug.report(InfoCode::AddingSNP,
                                     new std::string[5]
                                     {
                                             contigs[contigIndex],
                                             alig->seqsName[seqPos],
                                             std::to_string(feature.position),
                                             std::string(1, feature.ref[0]),
                                             std::string(1, feature.alt[0])
                                     });
                                donorSequenceChar = feature.alt[0];
                            }
                            // Check if the alternative sequence is the same in the sequence to apply
                            else if (donorSequenceChar == feature.alt[0])
                            {
                                debug.report(WarningCode::SNPAlreadApplied,
                                     new std::string[5]
                                         {
                                                 contigs[contigIndex],
                                                 alig->seqsName[seqPos],
                                                 std::to_string(feature.position),
                                                 std::string(1, feature.ref[0]),
                                                 std::string(1, feature.alt[0])
                                         });
                            }
                            // Check if reference character is the same as expected in the reference sequence
                            else if (referenceSequenceChar == std::tolower(feature.ref[0]) ||
                                     referenceSequenceChar == std::toupper(feature.ref[0]))
                            {
                                debug.report(ErrorCode::OverwrittingSNP,
                                     new std::string[7]
                                         {
                                                 contigs[contigIndex],
                                                 alig->seqsName[seqPos],
                                                 std::to_string(feature.position),
                                                 std::string(1, referenceSequenceChar),
                                                 std::string(1, feature.alt[0]),
                                                 std::string(1, referenceSequenceChar),
                                                 std::string(1, donorSequenceChar),
                                         });
                                donorSequenceChar = feature.alt[0];
                            }
                            else {
                                debug.report(ErrorCode::ReferenceNucleotideNotCorresponding,
                                     new std::string[5]
                                     {
                                        contigs[contigIndex],
                                        alig->seqsName[donorID],
                                        std::to_string(feature.position),
                                        std::string(1, feature.ref[0]),
                                        std::string(1, referenceSequenceChar)
                                     });
                            }
                        }

                    }
                }
            }
            delete[] line;
        }


        inline  void extractContigsFromLine(
                char *const line,
                StringVector &contigs)
        {
            // Remove first two characters;
            memmove(line, line + 2, strlen(line + 2) + 2);
            // Print result
            char *field_name = std::strtok(line, "=");

            // We only want the contig fields
            if (!strcmp(field_name, "contig")) {
                std::strtok(nullptr, "=");
                char *field_info = std::strtok(nullptr, ",");
                char *fname = new char[std::strlen(field_info) + 1];
                memmove(fname, field_info, strlen(field_info));
                fname[std::strlen(field_info)] = '\0';

                int U;
                // Check if the contig has already been added.
                for (U = 0; U < contigs.size(); U++)
                    if (contigs[U] == fname)
                        break;

                // If not, add it to the vector.
                if (U == contigs.size())
                    contigs.emplace_back(fname);

                delete[] fname;
            }
        }

        inline void extractDonorsFromLine(
                char *const line,
                std::vector<int> &currentDonorPosition,
                StringVector &donors)
        {
            // We only want to parse the FORMAT line, which contains the donors.
            std::strtok(strstr(line, "FORMAT"), "\t");

            char *token = std::strtok(nullptr, "\t");
            while (token != nullptr) {
                char *fname = new char[std::strlen(token) + 1];
                memmove(fname, token, strlen(token));
                fname[std::strlen(token)] = '\0';

                int U;

                // Check every other donor already added.
                for (U = 0; U < donors.size(); U++) {
                    if (donors[U] == fname) {
                        break;
                    }
                }

                // If not present, we add it
                if (U == donors.size()) {
                    donors.emplace_back(fname);
                }

                // If already present, warn the user.
                else {
                    debug.report(WarningCode::DonorAlreadyAdded, &fname[0]);
                }

                // We sum 1 as the position 0 is for reference
                currentDonorPosition.emplace_back(U + 1);

                token = std::strtok(nullptr, "\t");
                delete[] fname;
            }
        }

        inline  void obtainContigsAndDonors(
                const StringVector &filenames,
                StringVector &donors,
                StringVector &contigs,
                std::vector<std::vector<int>> &donorsPositions)
        {
            const int bufSize = 4096;
            char *line = new char[bufSize];

            // Obtain contigs and donors from all VCF files.
            for (int C = 0; C < filenames.size(); C++) {

                // Add a new vector<int>, which holds the sequence index
                //      for all donors present on this VCF
                donorsPositions.emplace_back(std::vector<int>());

                // Open the VCF file
                std::ifstream infile(filenames[C]); // Deleted with RAII

                if (!infile.good()) {
                    debug.report(ErrorCode::CantOpenFile, &filenames[C][0]);
                    continue;
                }

                // Read file to get all the donors and contigs present on the VCF
                while (infile.getline(line, bufSize, '\n')) {
                    // Extract contigs
                    if (strncmp("##contig", line, 8) == 0)
                            extractContigsFromLine(line, contigs);
                    // Extract donors
                    else if (strncmp("#CHROM", line, 6) == 0)
                        extractDonorsFromLine(line, donorsPositions[C], donors);
                }
            }
            delete[] line;
        }
    }

     void readVCF(
             const AlignmentVector &sources,
             const StringVector &filenames,
             const float minQuality,
             const float minCoverage,
             const bool ignoreFilter,
             const char *const replacementChar)
    {
        // All donors present in the vcf files.
        StringVector donors = StringVector();
        // Contigs present on the VCF files. Each entry must be present in the sources vector.
        StringVector contigs = StringVector();
        // Position on the reference files of each donor.
        //  Donors positions is a vector, of vectors, of ints.
        //  The first level refers to the VCF index in filenames
        //  The second level refers to the donor index in the VCF
        //  The value obtained is the sequence index of that donor in an alignment.
        std::vector<std::vector<int>> donorsPositions = std::vector<std::vector<int>>();

        {
            using namespace ngs::__internal;

            // Get a list of donors, contigs and a map, called donorsPositions
            obtainContigsAndDonors(filenames, donors, contigs, donorsPositions);

            // Check if all contigs have a reference alignment.
            if (!checkContigsWithReference(sources, contigs))
                return;

            // Extend the reference files with new sequences from donors.
            increaseSequencesInAlignment(sources, donors);

            applyVariantCallingFiles(
                    sources, filenames, contigs, donorsPositions,
                    minQuality, minCoverage, ignoreFilter, replacementChar);
        }
    }
}
