#include <ReadWriteMS/vcf_statish.h>
#include <fstream>

#include <reportsystem.h>
#include <cstring>
#include <vector>
#include <iomanip>

vcf_statish::vcf_statish()
{

}

vcf_statish::~vcf_statish()
{

}

void vcf_statish::readVCF(std::string filename)
{
    vcf_file vcff = vcf_file();
    
    std::cout << "~~> SNP's" << std::endl;
    
    std::ifstream infile;
    infile.open(filename);

    if (!infile.is_open())
    {
        debug.report(ErrorCode::CantOpenFile, new std::string[1] { filename });
    }
    
    char * line = new char [4096];
    int counter = 0;
    while (counter < 100 && infile.getline(line, 4096, '\n'))
    {
        
        if (line[0] == '#')
        {
            if (line[1] == '#')
            {
                // Remove first two characters;
                memmove (line, line+2, strlen (line+2) + 2);
                // Print result
                char * field_name = std::strtok(line, "=");
                
                if (!strcmp(field_name, "contig"))
                {           
                    std::strtok(NULL, "=");
                    char * field_info = std::strtok(NULL, ",");
                    char * fname = new char [std::strlen(field_info)];
                    memmove(fname, field_info, strlen(field_info));
                    vcff.filenames.push_back(fname);
                }
            }
            
            else
            {
                
                std::strtok(strstr(line, "FORMAT"), "\t");
                
                char * token = std::strtok(NULL, "\t");
                while (token != NULL)
                {
                    char * fname = new char [std::strlen(token)];
                    memmove (fname, token, strlen (token));
                    vcff.entries.push_back(fname);
                    token = std::strtok(NULL, "\t");
                }
                
            }
        }
        
        else
        {
            vcf_file::snp_entry snp = vcf_file::snp_entry();
            vcff.snps.push_back(snp);
            char * tmp;
            
//             std::cout << line << std::endl;
            
            tmp = std::strtok(line, "\t");
            snp.chromosome = new char[strlen(tmp)];
            std::memmove(snp.chromosome, tmp, strlen(tmp));
            
            tmp = std::strtok(NULL, "\t");
            snp.position = atoi(tmp);
            
            tmp = std::strtok(NULL, "\t");
            tmp = std::strtok(NULL, "\t");
            
            snp.ref = new char[strlen(tmp)];
            std::memmove(snp.ref, tmp, strlen(tmp));
            
            tmp = std::strtok(NULL, "\t");
            snp.alt = new char[strlen(tmp)];
            std::memmove(snp.alt, tmp, strlen(tmp));
            
            tmp = std::strtok(NULL, "\t");
            snp.quality = atof(tmp);
            
            tmp = std::strtok(NULL, "\t");
            snp.filter = std::strcmp(tmp, "PASS") ? false : true;
            
            
            // INFO
            std::strtok(NULL, "\t");
            // FORMAT
            std::strtok(NULL, "\t");
            // Individues with this SNP
            
            tmp = std::strtok(NULL, "\t");
            int counter = 0;
            while (tmp != NULL)
            {
                if (!strcmp(tmp, "0"))
                {
                    snp.present_in.push_back(counter);
                }
                std::cout   << vcff.entries[counter] << " "
                            << (std::strcmp(tmp, "0") ? "has the snp" : "doesn't have the snp") 
                            << " at " << snp.chromosome << " " << snp.position
                            << std::endl;
                counter ++;
                tmp = std::strtok(NULL, "\t");
            }
            
            std::cout   << std::left        << std::setw(15)
                        << snp.chromosome   << std::setw(15)
                        << snp.position     << std::setw(15)
                        << snp.ref          << std::setw(15) 
                        << "~>"             << std::setw(15)
                        << snp.alt          << std::setw(15)
                        << snp.quality      << std::setw(15)
                        << (snp.filter ? "PASS" : "FILTERED") 
                        
                        << std::endl;
        }
    }
//     std::cout << "~~>  Entries in genome" << std::endl;
//     for (char * fname : vcff.filenames)
//     {
//         std::cout << fname << std::endl;
//     }
//     
//     std::cout << "~~>  Individues" << std::endl;
//     for (char * fname : vcff.entries)
//     {
//         std::cout << fname << std::endl;
//     }
}
