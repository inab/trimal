#ifndef PHYLIP40M10STATE_H
#define PHYLIP40M10STATE_H

#include "BaseFormatHandler.h"

namespace FormatHandling {
class phylip40_m10_state : public BaseFormatHandler {
public:

    explicit phylip40_m10_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "phylip40_m10";
        extension = "phy2";
        canLoad = false;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    Alignment *LoadAlignment(const std::string &filename) override;

    bool SaveAlignment(const Alignment &alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(const std::string &FormatName) override;

};
}

#endif // PHYLIP40STATE_H
