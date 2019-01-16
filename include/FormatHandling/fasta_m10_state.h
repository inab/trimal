#ifndef FASTAM10STATE_H
#define FASTAM10STATE_H

#include "BaseFormatHandler.h"

namespace FormatHandling {
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

    Alignment *LoadAlignment(const std::string &filename) override;

    bool SaveAlignment(const Alignment &alignment, std::ostream *output) override;

    bool RecognizeOutputFormat(const std::string &FormatName) override;

};
}
#endif // FASTASTATE_H
