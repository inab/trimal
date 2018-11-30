#ifndef NEXUSM10STATE_H
#define NEXUSM10STATE_H

#include "BaseFormatHandler.h"

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

    Alignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(Alignment *alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};

#endif // NEXUSSTATE_H
