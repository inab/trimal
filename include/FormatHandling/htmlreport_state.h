#ifndef HTMLSTATE_H
#define HTMLSTATE_H

#include "BaseFormatHandler.h"

namespace FormatHandling {
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

    Alignment *LoadAlignment(const std::string &filename) override;

    bool SaveAlignment(const Alignment &alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(const std::string &FormatName) override;

};
}

#endif // HTMLSTATE_H
