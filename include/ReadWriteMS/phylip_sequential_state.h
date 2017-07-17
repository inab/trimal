#ifndef PHYLIPSEQUENTIALSTATE_H
#define PHYLIPSEQUENTIALSTATE_H

#include "readwrites.h"

class PhylipSequentialState : public readwrites
{
public:
    
    PhylipSequentialState(ReadWriteMS* MachineState) { Machine = MachineState; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(istream* origin);
    virtual void SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // PHYLIPSEQUENTIALSTATE_H
