//
// Created by vfernandez on 27/02/19.
//

#include <iostream>
#include <limits>
#include <cxxabi.h>
#include <iomanip>
#include <typeinfo>

template<typename T>
std::string type_name()
{
    int status;
    std::string tname = typeid(T).name();
    char *demangled_name = abi::__cxa_demangle(tname.c_str(), NULL, NULL, &status);
    if(status == 0) {
        tname = demangled_name;
        std::free(demangled_name);
    }
    return tname;
}

template <typename T>
void printSizeOfTypes()
{
    std::cout   << std::setw(15) << std::right << type_name<T>()
                << std::setw(5)  << std::right << sizeof(T)
                << std::setw(25) << std::right << static_cast<unsigned long>(std::numeric_limits<T>::max())
                << std::setw(30) << std::right << static_cast<long>(std::numeric_limits<T>::min())

                << "\n";
}

typedef int internalType;

int main()
{
    printSizeOfTypes<char>();
    printSizeOfTypes<short>();
    printSizeOfTypes<int>();
    printSizeOfTypes<long>();

    printSizeOfTypes<unsigned char>();
    printSizeOfTypes<unsigned short>();
    printSizeOfTypes<unsigned int>();
    printSizeOfTypes<unsigned long>();

}