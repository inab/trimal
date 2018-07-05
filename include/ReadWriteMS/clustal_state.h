#ifndef CLUSTALSTATE_H
#define CLUSTALSTATE_H

#include "ReadWriteBaseState.h"

class clustal_state : public ReadWriteBaseState {
public:

    explicit clustal_state(ReadWriteMS *MachineState) {
        Machine = MachineState;
        name = "clustal";
        extension = "clw";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(istream *origin) override;

    newAlignment *LoadAlignment(string filename) override;

    bool SaveAlignment(newAlignment *alignment, ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string FormatName) override;

};

#endif // CLUSTALSTATE_H
