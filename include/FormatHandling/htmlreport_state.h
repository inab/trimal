#ifndef HTMLSTATE_H
#define HTMLSTATE_H

#include "BaseFormatHandler.h"

class htmlreport_state : public BaseFormatHandler {
public:

    explicit htmlreport_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "html";
        extension = "html";
        canSave = true;
        canLoad = false;
    };

    int CheckAlignment(std::istream *origin) override;

    Alignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(Alignment *alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};

#endif // HTMLSTATE_H
