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
    
    bool shortNames = false;
    bool keepHeader = false;
    
    newAlignment* loadAlignment(std::string inFile);
    void saveAlignment(std::__cxx11::string outFile, std::__cxx11::string outFormat, newAlignment* alignment);
    
    void processFile(std::string inFile, std::string outFile, std::string outFormat);
    void processFile(std::string inFile, std::string outPattern, std::vector<std::string> outFormats[]);
     
//     void processFile(std::string inFiles[], std::string outFile, std::string outFormat);
//     void processFile(std::string inFiles[], std::string outPattern, std::string outFormats[]);
};


#endif // READWRITEMS_H
