#ifndef PHYLIP40STATE_H
#define PHYLIP40STATE_H

#include "ReadWriteBaseState.h"

class phylip40_state : public ReadWriteBaseState {
public:

    explicit phylip40_state(ReadWriteMS *MachineState) {
        Machine = MachineState;
        name = "phylip40";
        extension = "phy2";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(istream *origin) override;

    newAlignment *LoadAlignment(string filename) override;

    bool SaveAlignment(newAlignment *alignment, ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string FormatName) override;

};

#endif // PHYLIP40STATE_H
