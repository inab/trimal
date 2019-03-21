#ifndef NEXUSM10STATE_H
#define NEXUSM10STATE_H

#include "BaseFormatHandler.h"

namespace FormatHandling {
class nexus_m10_state : public BaseFormatHandler {
public:

    explicit nexus_m10_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "nexus_m10";
        extension = "nxs";
        canLoad = false;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    Alignment *LoadAlignment(const std::string &filename) override;

    bool SaveAlignment(const Alignment &alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(const std::string &FormatName) override;

};
}

#endif // NEXUSSTATE_H
