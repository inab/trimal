#ifndef MEGASSTATE_H
#define MEGASSTATE_H

#include "ReadWriteBaseState.h"

class MegaSequentialState : public ReadWriteBaseState
{
public:
    
    MegaSequentialState(ReadWriteMS* MachineState) { Machine = MachineState; name="mega_sequential"; extension="mega"; canLoad=true ; canSave=true; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(string filename);
    virtual bool SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // MEGASSTATE_H
