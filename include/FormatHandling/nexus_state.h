#ifndef NEXUSSTATE_H
#define NEXUSSTATE_H

#include "BaseFormatHandler.h"

namespace FormatHandling {
class nexus_state : public BaseFormatHandler {
public:

    explicit nexus_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "nexus";
        extension = "nxs";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    Alignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(Alignment *alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};
}

#endif // NEXUSSTATE_H
