#pragma once

#include "reportsystem.h"
#include "trimalManager.h"
#include "../../catch.hpp"
#include "../../testingUtils/ScopedRedirect.h"

#include "../../testingUtils/changePathToRoot.h"
#include "../../testingUtils/trimAlManagerArgumentReport.h"
#include "../../testingUtils/hasher.h"

namespace trimAlManagerTestIncludes
{
    extern const char * const alignmentPath;
    extern const char * const aminoAcidAlignmentPath;
    extern const char * const dnaAlignmentPath;
    extern const char * const rnaAlignmentPath;
    extern const char * const ambiguousAminoAcidAlignmentPath;
    extern const char * const degenerateDnaAlignmentPath;
    extern const char * const degenerateRnaAlignmentPath;
    extern const char * const alternativeAminoAcidAlignmentPath;

    namespace thresholds
    {
        /**
         * Structure to store a statistic
         */
        struct thresholdCase {
            
            const char * name; /**< Name of the case */
            
            const char * thresholdArgument; /**< Argument used to request the stat */
            
            const char * windowArgument; /**< Argument used to provide an argument for the stat */
            
            float & thresholdRef; /**< Reference variable that stores the threshold provided */
            
            int & windowRef; /**< Reference variable that stores the threshold provided */

            bool operator==(const thresholdCase& other) const
            {
                return other.name == this->name;
            }

            bool operator!=(const thresholdCase& other) const
            {
                return other.name != this->name;
            }
        };

        struct reportChecker
        {
            const char * argumentValue;
            trimAlManager::argumentReport expectedReport;
            bool expectedReturn;
        };
    }

    namespace formats
    {
        extern std::vector<const char *> legacyFormats;

        extern std::vector<const char *> nonLegacyFormats;
    }

    namespace backtranslation
    {
        extern const char * backtranslationAlignmentPath;
        extern const char * backtranslationCodonPath;
    }

    struct argumentReference
    {
        const char * argument;
        bool & reference;

        bool operator==(const argumentReference& other) const
        {
            return other.reference == this->reference;
        }
    };

    const std::vector<argumentReference>
            getAutomatedMethods(trimAlManager & manager);

    const std::vector<argumentReference>
            getStatsArguments(trimAlManager& manager, bool includeGeneral, bool includeSpecial);

    struct selectXCase {
        const char * name;
        const char * argument;
        const int maxValue;
        int ** pointerToPointerVar;
    };

    const std::vector<selectXCase>
            getSelectX(trimAlManager &manager, Alignment *alignment);

    const std::vector<argumentReference>
            getBooleanModifiers(trimAlManager & manager);

    struct reportCase
    {
        const char * argument;
        char ** pointerToPointerVar;

        bool operator==(const reportCase& other) const {
            return other.pointerToPointerVar == this->pointerToPointerVar;
        }
    };

    const std::vector<reportCase>
            getReportArguments(trimAlManager & manager);

    const std::vector<thresholds::thresholdCase>
            getThresholdCases(trimAlManager & manager, bool includeConsistency);
}

