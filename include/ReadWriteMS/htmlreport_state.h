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
        canLoad = false;
    };

    int CheckAlignment(std::istream *origin) override;

    newAlignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(newAlignment *alignment, std::ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};

#endif // HTMLSTATE_H
