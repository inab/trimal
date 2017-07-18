#ifndef MEGAISTATE_H
#define MEGAISTATE_H

#include "ReadWriteBaseState.h"

class MegaInterleavedState : public ReadWriteBaseState
{
public:
    
    MegaInterleavedState(ReadWriteMS* MachineState) { Machine = MachineState; name="MEGAI"; extension="mega"; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(string filename);
    virtual void SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // MEGAISTATE_H
