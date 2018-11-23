#ifndef MEGASSTATE_H
#define MEGASSTATE_H

#include "BaseFormatHandler.h"

class mega_sequential_state : public BaseFormatHandler {
public:

    explicit mega_sequential_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "mega_sequential";
        extension = "mega";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    Alignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(Alignment *alignment, std::ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};

#endif // MEGASSTATE_H
