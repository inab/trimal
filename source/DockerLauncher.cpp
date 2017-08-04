// # NOTE As there is no command processor this doesn't work at all. NOTE

// #include <iostream>
// #include <sstream>
// #include <cstring>
// #include <stdio.h>
// 
// 
// 
// int main(int argc, char *argv[])
// {
//     std::ostringstream s;
//     if (argc == 1) 
//     {
//         std::cout << "Running trimAl" << std::endl;
//         s << "/trimal" ;
//     }
//     else
//     {
// 
//         if (!std::strcmp(argv[1], "trimal"))
//         {
//             std::cout << "Hello there. You asked to run trimAl" << std::endl;
//             s << "/trimal" ;
//         }
//         else if (!std::strcmp(argv[1], "readal"))
//         {
//             std::cout << "Hello there. You asked to run readAl" << std::endl;
//             s << "/readal" ; 
//         }
//          else if (!std::strcmp(argv[1], "statal"))
//         {
//             std::cout << "Hello there. You asked to run statAl" << std::endl;
//             s << "/statal" ;
//         }
//         else
//         {
//             std::cerr << "Program " << argv[1] << " not recognized" << std::endl;
//             std::cout << "Exiting..." << std::endl;
//             return -1;
//         }
//         
//     }
//     
//     for (int i = 2; i < argc; i++)
//     {
//         s << " " << argv[i];
//     }
// 
//     char buff[512];
// 
//     FILE *program = popen(s.str().c_str(), "r");
//     if (!program) 
//         std::cout << "Ooopsy" << std::endl;
//     while(fgets(buff, sizeof(buff), program)!=NULL){
//         std::cout << buff;
//     }
//     pclose(program);
// }
