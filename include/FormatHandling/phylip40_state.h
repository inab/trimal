#ifndef PHYLIP40STATE_H
#define PHYLIP40STATE_H

#include "BaseFormatHandler.h"

class phylip40_state : public BaseFormatHandler {
public:

    explicit phylip40_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "phylip40";
        extension = "phy2";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    Alignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(Alignment *alignment, std::ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};

#endif // PHYLIP40STATE_H
