#ifndef PHYLIPPAMLSTATE_H
#define PHYLIPPAMLSTATE_H

#include "ReadWriteBaseState.h"

class phylip_paml_state : public ReadWriteBaseState {
public:

    explicit phylip_paml_state(ReadWriteMS *MachineState) {
        Machine = MachineState;
        name = "phylip_paml";
        extension = "phy";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(istream *origin) override;

    newAlignment *LoadAlignment(string filename) override;

    bool SaveAlignment(newAlignment *alignment, ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string FormatName) override;

};

#endif // PHYLIPPAMLSTATE_H
