#ifndef READWRITES_H
#define READWRITES_H

#include "../../include/newAlignment.h"
#include <fstream>


class ReadWriteMS;
class ReadWriteBaseState
{
public:
    
    string name;
    
    virtual int CheckAlignment(istream* origin) = 0;
    virtual newAlignment* LoadAlignment(string filename) = 0;
    virtual void SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName) = 0;
    virtual bool RecognizeOutputFormat(std::string FormatName) = 0;
    
    virtual ~ReadWriteBaseState() {};
    
protected:
    
    ReadWriteMS* Machine;
};

#endif // READWRITES_H
