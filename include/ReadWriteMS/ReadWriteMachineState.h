#ifndef READWRITEMS_H
#define READWRITEMS_H

#include <vector>
#include <string>
#include "../../include/newAlignment.h"

class ReadWriteBaseState;
class ReadWriteMS
{
public:
    
    ReadWriteMS();
    ~ReadWriteMS();
    
private:
    std::vector<ReadWriteBaseState*> available_states;
    void addState(ReadWriteBaseState* newState);
    
public:
    
    bool hasOutputFile  = true;
    bool shortNames     = false;
    bool keepHeader     = false;
    bool reverse        = false;
    
    newAlignment* loadAlignment(std::string inFile);
    void saveAlignment(std::string outFile, std::string outFormat, newAlignment* alignment);

    void processFile(std::vector< std::string >* inFile, std::string* outPattern, std::vector< std::string >* outFormats);
};


#endif // READWRITEMS_H
