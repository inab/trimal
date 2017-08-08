#ifndef PIRSTATE_H
#define PIRSTATE_H

#include "ReadWriteBaseState.h"

class PirState : public ReadWriteBaseState
{
public:
    
    PirState(ReadWriteMS* MachineState) { Machine = MachineState; name="pir"; extension="pir";canLoad=true ; canSave=true;; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(string filename);
    virtual bool SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // PIRSTATE_H
