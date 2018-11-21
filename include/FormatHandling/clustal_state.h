#ifndef CLUSTALSTATE_H
#define CLUSTALSTATE_H

#include "BaseFormatHandler.h"

class clustal_state : public BaseFormatHandler {
public:

    explicit clustal_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "clustal";
        extension = "clw";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    newAlignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(newAlignment *alignment, std::ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};

#endif // CLUSTALSTATE_H
