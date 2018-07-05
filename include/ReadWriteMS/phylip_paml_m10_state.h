#ifndef PHYLIPPAMLM10STATE_H
#define PHYLIPPAMLM10STATE_H

#include "ReadWriteBaseState.h"

class phylip_paml_m10_state : public ReadWriteBaseState {
public:

    explicit phylip_paml_m10_state(ReadWriteMS *MachineState) {
        Machine = MachineState;
        name = "phylip_paml_m10";
        extension = "phy";
        canLoad = false;
        canSave = true;
    };

    int CheckAlignment(istream *origin) override;

    newAlignment *LoadAlignment(string filename) override;

    bool SaveAlignment(newAlignment *alignment, ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string FormatName) override;

};

#endif // PHYLIPPAMLSTATE_H
