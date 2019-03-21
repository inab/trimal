#ifndef FASTASTATE_H
#define FASTASTATE_H

#include "BaseFormatHandler.h"

namespace FormatHandling {
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

    Alignment *LoadAlignment(const std::string &filename) override;

    bool SaveAlignment(const Alignment &alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(const std::string &FormatName) override;

};
}
#endif // FASTASTATE_H
