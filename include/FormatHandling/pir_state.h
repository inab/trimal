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

    Alignment *LoadAlignment(const std::string &filename) override;

    bool SaveAlignment(const Alignment &alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(const std::string &FormatName) override;

};
}

#endif // PIRSTATE_H
