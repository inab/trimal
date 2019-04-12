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

#ifndef VCF_STATISH_H
#define VCF_STATISH_H

#include "reportsystem.h"
#include "Alignment/Alignment.h"

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <vector>

namespace ngs {

    typedef std::vector<Alignment *> AlignmentVector;
    typedef std::vector<std::string> StringVector;

    namespace IUPAC {

        /**
         * Method to obtain a tag (see bottom of description) from a sequence of residues in the format\n
         * "CSCSCSCSC...C" where C is a character and S the separator.\n
         * It will error if more than one character is present between the separators. Ej: CSCCS \n\n
         *
         * <b> Binary Tag Enum.</b> \n
         * <a href="http://www.alanzucconi.com/2015/07/26/enum-flags-and-bitwise-operators/" > Alan Zucconi </a> has a great explanation of these.\n\n
         * <b> Values: </b>A = \i 1; C = \i 2; T = \i 4; G = \i 8;
         * @param [in] array Sequence of chars to convert
         * @param [in] separator Separator of the chars in the sequence
         * @return <b>tag</b> if correct\n<b>-1</b> if errored
         */
        inline int getTagFromCharArray(
                const char *const array,
                const char separator);

        /**
         * Method to obtain the char represented by a IUPAC tag.\n\n
         *
         * <b> Binary Tag Enum.</b> \n
         * <a href="http://www.alanzucconi.com/2015/07/26/enum-flags-and-bitwise-operators/" > Alan Zucconi </a> has a great explanation of these.\n\n
         * <b> Values: </b>A = \i 1; C = \i 2; T = \i 4; G = \i 8;
         * @param [in] tag Int representing the tag to convert.\n
         * Obtained by ::getTagFromCharArray
         * @return <b>Char</b> represented by the tag\n
         * <b>-</b> if tag not recognized
         */
        inline char getCharFromTag(
                const int tag);

    }

    namespace __internal {

        struct vcfFeature {
        public:
            bool filter;
            int position;
            float quality;
            float readDepthIndex;
            char
                * ref    = nullptr,
                * alt    = nullptr,
                * contig = nullptr;

            StringVector donorsInfo {};

            vcfFeature() : position(0), quality(0.0F), filter(false), readDepthIndex(0)
            {
                ref = nullptr;
                alt = nullptr;
                contig = nullptr;
            }

            void reset()
            {
                delete [] ref;
                ref = nullptr;

                delete [] alt;
                alt = nullptr;

                delete [] contig;
                contig = nullptr;

                donorsInfo.clear();
            }

            ~vcfFeature()
            {
                reset();
            }
        };

        void printApeek(AlignmentVector &sources);

        /**
         * Method to check if there is an alignment for each contig extracted
         * @param [in] sources Alignments obtained through splitting reference alignment
         * @param [in] contigs Contigs extracted from VCFs
         * @return \i True if all contigs have a reference /n\i False otherwise
         */
        inline bool checkContigsWithReference(
                const AlignmentVector &sources,
                const StringVector &contigs);

        /**
         * Method to increase the number of sequences on each splitted reference.
         * Each alignment should contain only one sequence -> the reference one.
         * The method will make a copy of the reference sequence and put the
         * donor name on it for each donor present on donors.
         * @param [in,out] sources Alignments to increase their sequences
         * @param [in] donors Donor names
         */
        inline void increaseSequencesInAlignment(
                const AlignmentVector &sources,
                const StringVector &donors);

        /**
         * Method to extract a variant from a char array
         * @param [in] line Char array containing information to form a variant
         * @param [out] feature Variant reference to store results
         * @return \i True if done \n\i False if line does starts by '#'
         */
        inline bool extractFeature(
                const char *const line,
                vcfFeature &feature);

        /**
         * Filter a variant based on multiple filters:
         * Size of reference if different than 1.
         * Size of alternative is different than 1.
         * Quality of variant is lesser than minQuality.ss
         * @param [in] tmpVariant Variant to filter
         * @param [in] minQuality Min quality of variant to be kept
         * @param [in] ignoreFilter Ignore filter field on VCF
         * @return \i True if variant is valid \n\i False otherwise
         */
        inline bool filterFeature(const vcfFeature &tmpVariant, const float minQuality);

        /**
         *
         * @param [in,out] sources Alignment vector containing all the result alignments
         * @param [in] filenames Filenames of the VCF files to be applied
         * @param [in] contigs Contigs found along all teh VCFs
         * @param [in] donorsPositions Structure to map donors of each VCF file to
         * it's corresponding sequence position
         * @param [in] minQuality Min quality to apply a feature
         * @param [in] minCoverage Min coverage to apply a feature
         * @param [in] ignoreFilter Ignore feature filter field on VCF?
         * @param [in] replacementChar Replacement char used when a SNP is not reliable
         */
        inline void applyVariantCallingFiles(
                const AlignmentVector &sources,
                const StringVector &filenames,
                const StringVector &contigs,
                const std::vector<std::vector<int>> &donorsPositions,
                const float minQuality,
                const float minCoverage,
                const bool ignoreFilter,
                const char *const replacementChar);

        /**
         * Method to extract contigs from a line of a VCF
         * @param [in] line Char array containing contigs information in VCF.
         * Format expected: "##contig=<ID=ContigName,length=x>"
         * @param [out] contigs Vector to store contigs found
         */
        inline void extractContigsFromLine(
                char *const line,
                StringVector &contigs);
        /**
         * Method to extract donors from a line of a VCF
         * @param [in] line Char array containing contigs information in VCF.
         * @param [in, out] currentDonorPosition Map of donors in VCF to sequence position in Alignment.
         * @param [out] donors List of donors found
         */
        inline void extractDonorsFromLine(
                char *const line,
                std::vector<int> &currentDonorPosition,
                StringVector &donors);
        /**
         *
         * @param [in] filenames Vector of VCF filenames to parse and apply.
         * @param [out] donors List of donors to be filled
         * @param [out] contigs List of contigs to be filled
         * @param [out] donorsPositions Map of donors in VCFs to sequence positions for each VCF file
         */
        inline void obtainContigsAndDonors(
                const StringVector &filenames,
                StringVector &donors,
                StringVector &contigs,
                std::vector<std::vector<int>> &donorsPositions);
    }

    /**
     *
     * @param [in,out] sources List of Alignments, each of them having one sequence, the reference one.\n
     * They fill be extended with one sequence for each donor found on the VCF dataset, and their SNP will be applied
     * to the correspondent sequence.
     * @param [in] filenames List of VCF files to parse and apply.
     * @param [in] minQuality Min quality a SNP should have to be applied.
     * @param [in] minCoverage Min coverage a SNP should have to be applied.
     * @param [in] ignoreFilter Whether to ignore filter field in VCF
     * @param [in] replacementChar Replacement char to apply when VCF-feature doesn't pass the filter.
     */
    void readVCF(
            const AlignmentVector &sources,
            const StringVector &filenames,
            const float minQuality,
            const float minCoverage,
            const bool ignoreFilter,
            const char *const replacementChar);
}


#endif // VCF_STATISH_H
