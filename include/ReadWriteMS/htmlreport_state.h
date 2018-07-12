#ifndef HTMLSTATE_H
#define HTMLSTATE_H

#include "ReadWriteBaseState.h"

class htmlreport_state : public ReadWriteBaseState {
public:

    explicit htmlreport_state(ReadWriteMS *MachineState) {
        Machine = MachineState;
        name = "html";
        extension = "html";
        canSave = true;
    };

    int CheckAlignment(istream *origin) override;

    newAlignment *LoadAlignment(string filename) override;

    bool SaveAlignment(newAlignment *alignment, ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string FormatName) override;

};

#endif // HTMLSTATE_H
