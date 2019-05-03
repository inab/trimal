/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2009-2019
        Fernandez-Rodriguez V.  (victor.fernandez@bsc.es)
        Capella-Gutierrez S.    (salvador.capella@bsc.es)
        Gabaldon, T.            (tgabaldon@crg.es)

    This file is part of trimAl/readAl.

    trimAl/readAl are free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl/readAl are distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl/readAl. If not, see <http://www.gnu.org/licenses/>.

***************************************************************************** */

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