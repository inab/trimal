#ifndef NEXUSSTATE_H
#define NEXUSSTATE_H

#include "ReadWriteBaseState.h"

class nexus_state : public ReadWriteBaseState {
public:

    explicit nexus_state(ReadWriteMS *MachineState) {
        Machine = MachineState;
        name = "nexus";
        extension = "nxs";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(istream *origin) override;

    newAlignment *LoadAlignment(string filename) override;

    bool SaveAlignment(newAlignment *alignment, ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string FormatName) override;

};

#endif // NEXUSSTATE_H
