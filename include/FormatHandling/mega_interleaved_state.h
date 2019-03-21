#ifndef MEGAISTATE_H
#define MEGAISTATE_H

#include "BaseFormatHandler.h"

namespace FormatHandling {
class mega_interleaved_state : public BaseFormatHandler {
public:

    explicit mega_interleaved_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "mega_interleaved";
        extension = "mega";
        canLoad = true;
        canSave = false;
    };

    int CheckAlignment(std::istream *origin) override;

    Alignment *LoadAlignment(const std::string &filename) override;

    bool SaveAlignment(const Alignment &alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(const std::string &FormatName) override;

};
}

#endif // MEGAISTATE_H
