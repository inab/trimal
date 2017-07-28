#ifndef CLUSTALSTATE_H
#define CLUSTALSTATE_H

#include "ReadWriteBaseState.h"

class ClustalState : public ReadWriteBaseState
{
public:
    
    ClustalState(ReadWriteMS* MachineState) { Machine = MachineState; name="CLUSTAL"; extension="clw"; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(string filename);
    virtual bool SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // CLUSTALSTATE_H
