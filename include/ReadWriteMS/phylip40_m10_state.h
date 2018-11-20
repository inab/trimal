#ifndef PHYLIP40M10STATE_H
#define PHYLIP40M10STATE_H

#include "ReadWriteBaseState.h"

class phylip40_m10_state : public ReadWriteBaseState {
public:

    explicit phylip40_m10_state(ReadWriteMS *MachineState) {
        Machine = MachineState;
        name = "phylip40_m10";
        extension = "phy2";
        canLoad = false;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    newAlignment *LoadAlignment(std::string filename) override;

    bool SaveAlignment(newAlignment *alignment, std::ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string FormatName) override;

};

#endif // PHYLIP40STATE_H
