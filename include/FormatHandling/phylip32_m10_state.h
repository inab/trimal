#ifndef PHYLIP32M10STATE_H
#define PHYLIP32M10STATE_H

#include "BaseFormatHandler.h"

class phylip32_m10_state : public BaseFormatHandler {
public:

    explicit phylip32_m10_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "phylip32_m10";
        extension = "phy";
        canLoad = false;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    Alignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(Alignment *alignment, std::ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};

#endif // PHYLIP32STATE_H
