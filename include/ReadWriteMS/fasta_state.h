#ifndef FASTASTATE_H
#define FASTASTATE_H

#include "ReadWriteBaseState.h"

class FastaState : public ReadWriteBaseState
{
public:
    
    FastaState(ReadWriteMS* MachineState) { Machine = MachineState; name="FASTA"; extension="fasta"; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(string filename);
    virtual bool SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // FASTASTATE_H
