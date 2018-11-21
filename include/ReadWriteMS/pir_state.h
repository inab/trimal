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
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    newAlignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(newAlignment *alignment, std::ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};

#endif // PIRSTATE_H
