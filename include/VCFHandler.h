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

    namespace IUPAC {

        /**
         * Method to obtain a tag from a sequence of residues in the format
         * "CSCSCSCSC...C" where C is a character and S the separator.
         * It will error if more than one character is present between the separators. Ej: CSCCS
         * A = 1; C = 2; T = 3; G = 4;
         * @param array sequence of chars to convert
         * @param separator separator of the chars in the sequence
         * @return -1 if errored, tag if correct
         */
        inline int getTagFromCharArray(
                const char *const array,
                const char separator);

        /**
         * Method to obtain a letter from a IUPAC tag
         * @param tag int representing the tag to convert
         * @return the char represented by the tag
         */
        inline char getCharFromTag(
                const int tag);

    }

    namespace __internal {

        struct vcfFeature {
        public:
            int position;
            float quality;
            char * ref    = nullptr,
                    * alt    = nullptr,
                    * contig = nullptr;
            bool filter;
            float readDepthIndex;

            std::vector<std::string> donorsInfo {};

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

        void printApeek(std::vector<Alignment *> &sources);

        /**
         * Method to check if there is an alignment for each contig extracted
         * @param sources Alignments obtain through splitting reference alignment
         * @param contigs Contigs extracted from VCFS
         * @return \i True if all contigs have a reference /n\i False otherwise
         */
        inline bool checkContigsWithReference(
                const std::vector<Alignment *> &sources,
                const std::vector<std::string> &contigs);

        /**
         * Method to increase the number of sequences on each splitted reference
         * Each alignment should contain only one sequence -> the reference one
         * The method will make a copy of the reference sequence and put the
         * donor name on it for each donor present on donors
         * @param sources Alignments to increase their sequences
         * @param donors Donor names
         */
        inline void increaseSequencesInAlignment(
                const std::vector<Alignment *> &sources,
                const std::vector<std::string> &donors);

        /**
         * Method to check if there is an alignment for each contig
         * and make a copy of their reference sequence for each donor.
         * @param sources Alignments to check and increase in number of sequences
         * @param contigs Contigs found in VCF. Used to check if there is an
         * alignment for each of them
         * @param donors Donor names, each of them will have a sequence on each
         * of the alignments in sources
         */
        inline void extendAlignments(
                const std::vector<Alignment *> &sources,
                const std::vector<std::string> &contigs,
                const std::vector<std::string> &donors);

        /**
         * Method to extract a variant from a char array
         * @param line char array containing information to form a variant
         * @param feature variant reference to store results
         * @return \i True if done \n\i False if line does starts by '#'
         */
        inline bool extractFeature(
                const char *const line,
                vcfFeature &feature);

        /**
         * Filter a variant based on multiple filters:
         * Size of reference if different than 1
         * Size of alternative is different than 1
         * Quality of variant is lesser than minQuality
         * Filter field present on the VCF if ignoreFilter is false
         * @param tmpVariant Variant to filter
         * @param minQuality Min quality of variant to be kept
         * @param ignoreFilter Ignore filter field on VCF
         * @return \i True if variant is valid \n\i False otherwise
         */
        inline bool filterFeature(const vcfFeature &tmpVariant, const float minQuality);

//        /**
//         * Method to apply a non-filtered-out variant to a set of sources
//         * It will filter-out the variant for each donor, based on its coverage
//         * @param tmpVariant VCF Feature
//         * @param seqPos Sequence position
//         * @param minCoverage Coverage to apply a feature
//         * @param replacementChar Char to apply when a feature has been filtered out
//         * @param sources Alignments to apply the variants
//         * @param contigs Contigs
//         * @param filename
//         */
//        inline void applyUnfilteredFeature(
//                const vcfFeature &tmpVariant,
//                const int seqPos,
//                const float minCoverage,
//                const char *const replacementChar,
//                const std::vector<Alignment *> &sources,
//                const std::vector<std::string> &contigs,
//                const std::string &filename);
//
//        /**
//         *
//         * @param tmpVariant
//         * @param seqPos
//         * @param replacementChar
//         * @param sources
//         * @param contigs
//         */
//        inline void applyFilteredFeature(
//                const vcfFeature &tmpVariant,
//                const int seqPos,
//                const char *const replacementChar,
//                const std::vector<Alignment *> &sources,
//                const std::vector<std::string> &contigs);

        /**
         *
         * @param sources
         * @param filenames
         * @param contigs
         * @param donorsPositions
         * @param minQuality
         * @param minCoverage
         * @param ignoreFilter
         * @param replacementChar
         */
        inline void applyVariantCallingFiles(
                const std::vector<Alignment *> &sources,
                const std::vector<std::string> &filenames,
                const std::vector<std::string> &contigs,
                const std::vector<std::vector<int>> &donorsPositions,
                const float minQuality,
                const float minCoverage,
                const bool ignoreFilter,
                const char *const replacementChar);

        /**
         *
         * @param line
         * @param contigs
         */
        inline void extractContigsFromLine(
                char *const line,
                std::vector<std::string> &contigs);
        /**
         *
         * @param line
         * @param C
         * @param donors
         * @param donorsPositions
         */
        inline void extractDonorsFromLine(
                char *const line,
                int C,
                std::vector<std::string> &donors,
                std::vector<std::vector<int>> &donorsPositions);
        /**
         *
         * @param filenames
         * @param donors
         * @param contigs
         * @param donorsPositions
         */
        inline void obtainContigsAndDonors(
                const std::vector<std::string> &filenames,
                std::vector<std::string> &donors,
                std::vector<std::string> &contigs,
                std::vector<std::vector<int>> &donorsPositions);
    }

    /**
     *
     * @param sources
     * @param filenames
     * @param minQuality
     * @param minCoverage
     * @param ignoreFilter
     * @param replacementChar
     */
    void readVCF(
            const std::vector<Alignment *> &sources,
            const std::vector<std::string> &filenames,
            const float minQuality,
            const float minCoverage,
            const bool ignoreFilter,
            const char *const replacementChar);
}


#endif // VCF_STATISH_H
