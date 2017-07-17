#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include <string>


int main(int argc, char *argv[])
{
    ReadWriteMS MachineState = ReadWriteMS();
    MachineState.shortNames = true;
    
    MachineState.processFile(
        "./dataset/superAlignment1.phy2", 
        "./dataset/Testing/outputs.$.txt", 
        new std::vector<std::string>{"fasta", "clustal", "pir", "phylip32", "phylip40", "nexus", "mega", "html"});
    
//     MachineState.processFile("./dataset/superAlignment1.mega", "./dataset/Testing/meganon.txt", "meganon");
//     MachineState.processFile("./dataset/superAlignment1.mega", "./dataset/Testing/megai.txt", "megai");
}
