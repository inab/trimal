#ifndef PHYLIP40STATE_H
#define PHYLIP40STATE_H

#include "ReadWriteBaseState.h"

class Phylip40State : public ReadWriteBaseState
{
public:
    
    Phylip40State(ReadWriteMS* MachineState) { Machine = MachineState; name="PHYLIP40"; extension="phy2"; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(string filename);
    virtual void SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // PHYLIP40STATE_H
