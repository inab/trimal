#ifndef PIRSTATE_H
#define PIRSTATE_H

#include "ReadWriteBaseState.h"

class pir_state : public ReadWriteBaseState {
public:

    explicit pir_state(ReadWriteMS *MachineState) {
        Machine = MachineState;
        name = "pir";
        extension = "pir";
        canLoad = true;
        canSave = true;;
    };

    int CheckAlignment(istream *origin) override;

    newAlignment *LoadAlignment(string filename) override;

    bool SaveAlignment(newAlignment *alignment, ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string FormatName) override;

};

#endif // PIRSTATE_H
