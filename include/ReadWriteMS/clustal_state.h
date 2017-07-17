#ifndef CLUSTALSTATE_H
#define CLUSTALSTATE_H

#include "readwrites.h"

class ClustalState : public readwrites
{
public:
    
    ClustalState(ReadWriteMS* MachineState) { Machine = MachineState; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(istream* origin);
    virtual void SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // CLUSTALSTATE_H
