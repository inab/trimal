#ifndef PIRSTATE_H
#define PIRSTATE_H

#include "BaseFormatHandler.h"

namespace FormatHandling {
class pir_state : public BaseFormatHandler {
public:

    explicit pir_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "pir";
        extension = "pir";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    Alignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(Alignment *alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};
}

#endif // PIRSTATE_H