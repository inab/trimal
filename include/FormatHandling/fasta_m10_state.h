#ifndef FASTAM10STATE_H
#define FASTAM10STATE_H

#include "BaseFormatHandler.h"

class fasta_m10_state : public BaseFormatHandler {
public:

    explicit fasta_m10_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "fasta_m10";
        extension = "fasta";
        canLoad = false;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    newAlignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(newAlignment *alignment, std::ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};

#endif // FASTASTATE_H
