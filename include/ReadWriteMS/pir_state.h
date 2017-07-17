#ifndef PIRSTATE_H
#define PIRSTATE_H

#include "readwrites.h"

class PirState : public readwrites
{
public:
    
    PirState(ReadWriteMS* MachineState) { Machine = MachineState; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(istream* origin);
    virtual void SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // PIRSTATE_H
