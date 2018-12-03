#ifndef PHYLIPPAMLM10STATE_H
#define PHYLIPPAMLM10STATE_H

#include "BaseFormatHandler.h"

namespace FormatHandling {
class phylip_paml_m10_state : public BaseFormatHandler {
public:

    explicit phylip_paml_m10_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "phylip_paml_m10";
        extension = "phy";
        canLoad = false;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    Alignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(Alignment *alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};
}
#endif // PHYLIPPAMLSTATE_H