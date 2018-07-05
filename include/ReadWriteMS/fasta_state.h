#ifndef FASTASTATE_H
#define FASTASTATE_H

#include "ReadWriteBaseState.h"

class fasta_state : public ReadWriteBaseState {
public:

    explicit fasta_state(ReadWriteMS *MachineState) {
        Machine = MachineState;
        name = "fasta";
        extension = "fasta";
        canLoad = true;
        canSave = true;
    };

    int CheckAlignment(istream *origin) override;

    newAlignment *LoadAlignment(string filename) override;

    bool SaveAlignment(newAlignment *alignment, ostream *output, std::string *FileName) override;

    bool RecognizeOutputFormat(std::string FormatName) override;

};

#endif // FASTASTATE_H
