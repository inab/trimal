#ifndef FASTASTATE_H
#define FASTASTATE_H

#include "BaseFormatHandler.h"

class fasta_state : public BaseFormatHandler {
public:

    explicit fasta_state(FormatManager *MachineState) {
        Machine = MachineState;
        name = "fasta";
        extension = "fasta";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(std::istream *origin) override;

    newAlignment *LoadAlignment(std::string& filename) override;

    bool SaveAlignment(newAlignment *alignment, std::ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string& FormatName) override;

};

#endif // FASTASTATE_H
