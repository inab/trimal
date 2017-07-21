#ifndef PHYLIPPAMLSTATE_H
#define PHYLIPPAMLSTATE_H

#include "ReadWriteBaseState.h"

class PhylipPamlState : public ReadWriteBaseState
{
public:
    
    PhylipPamlState(ReadWriteMS* MachineState) { Machine = MachineState; name="PHYLIP_PAML"; extension="phy"; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(string filename);
    virtual bool SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // PHYLIPPAMLSTATE_H
