#ifndef FASTASTATE_H
#define FASTASTATE_H

#include "readwrites.h"

class FastaState : public readwrites
{
public:
    
    FastaState(ReadWriteMS* MachineState) { Machine = MachineState; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(istream* origin);
    virtual void SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // FASTASTATE_H
