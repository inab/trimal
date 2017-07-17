#include "../../include/ReadWriteMS/ReadWriteMachineState.h"


int main(int argc, char *argv[])
{
    ReadWriteMS MachineState = ReadWriteMS();
    MachineState.shortNames = true;
    MachineState.processFile("./dataset/superAlignment1.mega", "./File.txt", "meganon");
    MachineState.processFile("./dataset/superAlignment1.mega", "./File2.txt", "megai");
}
