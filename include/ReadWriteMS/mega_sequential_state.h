#ifndef MEGANISTATE_H
#define MEGANISTATE_H

#include "readwrites.h"

class MegaSequentialState : public readwrites
{
public:
    
    MegaSequentialState(ReadWriteMS* MachineState) { Machine = MachineState; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(istream* origin);
    virtual void SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // MEGANISTATE_H
