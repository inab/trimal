#ifndef PHYLIPSTATE_H
#define PHYLIPSTATE_H

#include "readwrites.h"

class PhylipInterleavedState : public readwrites
{
public:
    
    PhylipInterleavedState(ReadWriteMS* MachineState) { Machine = MachineState; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(istream* origin);
    virtual void SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // PHYLIPSTATE_H
