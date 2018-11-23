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

    Alignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(Alignment *alignment, std::ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};

#endif // FASTASTATE_H
