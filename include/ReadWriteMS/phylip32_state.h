#ifndef PHYLIP32STATE_H
#define PHYLIP32STATE_H

#include "ReadWriteBaseState.h"

class Phylip32State : public ReadWriteBaseState
{
public:
    
    Phylip32State(ReadWriteMS* MachineState) { Machine = MachineState; name="phylip32"; extension="phy"; canLoad=true ; canSave=true; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(string filename);
    virtual bool SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // PHYLIP32STATE_H
