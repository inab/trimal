#ifndef MEGAISTATE_H
#define MEGAISTATE_H

#include "readwrites.h"

class MegaInterleavedState : public readwrites
{
public:
    
    MegaInterleavedState(ReadWriteMS* MachineState) { Machine = MachineState; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(istream* origin);
    virtual void SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // MEGAISTATE_H
