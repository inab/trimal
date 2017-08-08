#ifndef NEXUSSTATE_H
#define NEXUSSTATE_H

#include "ReadWriteBaseState.h"

class NexusState : public ReadWriteBaseState
{
public:
    
    NexusState(ReadWriteMS* MachineState) { Machine = MachineState; name="nexus"; extension="nxs"; canLoad=true ; canSave=true; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(string filename);
    virtual bool SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // NEXUSSTATE_H
