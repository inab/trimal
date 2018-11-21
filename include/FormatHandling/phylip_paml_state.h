#ifndef PHYLIPPAMLSTATE_H
#define PHYLIPPAMLSTATE_H

#include "BaseFormatHandler.h"

class phylip_paml_state : public BaseFormatHandler {
public:

    explicit phylip_paml_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "phylip_paml";
        extension = "phy";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    newAlignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(newAlignment *alignment, std::ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};

#endif // PHYLIPPAMLSTATE_H
