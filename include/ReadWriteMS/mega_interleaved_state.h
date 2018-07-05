#ifndef MEGAISTATE_H
#define MEGAISTATE_H

#include "ReadWriteBaseState.h"

class mega_interleaved_state : public ReadWriteBaseState {
public:

    explicit mega_interleaved_state(ReadWriteMS *MachineState) {
        Machine = MachineState;
        name = "mega_interleaved";
        extension = "mega";
        canLoad = true;
        canSave = false;
    };

    int CheckAlignment(istream *origin) override;

    newAlignment *LoadAlignment(string filename) override;

    bool SaveAlignment(newAlignment *alignment, ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string FormatName) override;

};

#endif // MEGAISTATE_H
