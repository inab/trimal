#ifndef READWRITEMS_H
#define READWRITEMS_H

#include <vector>
#include <string>

class readwrites;
class ReadWriteMS
{
public:
    
    ReadWriteMS();
    ~ReadWriteMS();
    
private:
    std::vector<readwrites*> available_states;
    readwrites* outState;
    readwrites* inState;
    
    void addState(readwrites* newState);
public:
    
    bool shortNames = false;
    
    void processFile(std::string inFile, std::string outFile, std::string outFormat);
//     void processFile(std::string inFile, std::string outPattern, std::string outFormats[]);
//     
//     void processFile(std::string inFiles[], std::string outPattern, std::string outFormat);
//     void processFile(std::string inFiles[], std::string outPattern, std::string outFormats[]);
};


#endif // READWRITEMS_H
