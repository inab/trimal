#ifndef NEXUSM10STATE_H
#define NEXUSM10STATE_H

#include "ReadWriteBaseState.h"

class nexus_m10_state : public ReadWriteBaseState {
public:

    explicit nexus_m10_state(ReadWriteMS *MachineState) {
        Machine = MachineState;
        name = "nexus_m10";
        extension = "nxs";
        canLoad = false;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    newAlignment *LoadAlignment(std::string filename) override;

    bool SaveAlignment(newAlignment *alignment, std::ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string FormatName) override;

};

#endif // NEXUSSTATE_H
