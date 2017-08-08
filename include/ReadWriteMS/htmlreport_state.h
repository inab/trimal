#ifndef HTMLSTATE_H
#define HTMLSTATE_H

#include "ReadWriteBaseState.h"

class HTMLState : public ReadWriteBaseState
{
public:
    
    HTMLState(ReadWriteMS* MachineState) { Machine = MachineState; name="html"; extension="html"; canSave=true; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(string filename);
    virtual bool SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // HTMLSTATE_H
