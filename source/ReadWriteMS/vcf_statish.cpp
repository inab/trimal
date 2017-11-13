#include <ReadWriteMS/vcf_statish.h>
#include <fstream>

#include <reportsystem.h>
#include <cstring>
#include <vector>
#include <iomanip>

#include <algorithm>

vcf_statish::vcf_statish()
{

}

vcf_statish::~vcf_statish()
{

}

void vcf_statish::readVCF(std::vector<newAlignment*> sources, std::string filename)
{

    std::ifstream infile;
    infile.open(filename);

    if (!infile.is_open())
    {
        debug.report(ErrorCode::CantOpenFile, new std::string[1] { filename });
    }
    
    std::vector<std::string> donors = std::vector<std::string>();
    std::vector<std::string> sequences = std::vector<std::string>();
    
    char * line = new char [4096];
    while (infile.getline(line, 4096, '\n'))
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
                    char * fname = new char [std::strlen(field_info) + 1];
                    memmove(fname, field_info, strlen(field_info));
                    fname[std::strlen(field_info)] = '\0';
                    sequences.push_back(fname);
                    delete[] fname;
                }
            }
            
            else
            {
                
                std::strtok(strstr(line, "FORMAT"), "\t");
                
                char * token = std::strtok(NULL, "\t");
                while (token != NULL)
                {
                    char * fname = new char [std::strlen(token) + 1];
                    memmove (fname, token, strlen (token));
                    fname[std::strlen(token)] = '\0';
                    donors.push_back(fname);
                    token = std::strtok(NULL, "\t");
                    delete[] fname;
                }
                break;
            }
        }
    }
    
    {
        bool checkIn = true;
        for (int x = 0; x < sequences.size(); x++)
        {
            int i;
            for (i = 0; i < sources.size(); i++)
            {
                if (!strcmp(&sequences[x][0], &sources[i]->seqsName[0][0]) )
                    break;
            }
            
            if (i == sources.size())
            {
                std::cout << sequences[x] << " not found in sources" << endl;
                checkIn = false;
            }
        }
        
        std::cout << "Has all files?: " << (checkIn ? "True" : "False") << std::endl;
        
        if (checkIn)
        {
            for (newAlignment * nA : sources)
            {
                
                std::string * oldSequences  = nA->sequences;
                std::string * oldNames      = nA->seqsName;
                
                std::string seq = nA->sequences[0];
                std::string nam = nA->seqsName[0];
                
                nA -> sequences = new std::string[donors.size() + nA->originalSequenNumber];
                nA -> sequences[0] = seq;
                
                nA -> seqsName = new std::string[donors.size() + nA->originalSequenNumber];
                nA -> seqsName[0] = nam;
                
                for (int i = 1; i < donors.size() + 1; i++)
                {
                    nA->sequences[i] = seq.c_str();
                    nA->seqsName[i] = donors[i - 1];
                }
                
                for (int i = 0; i < nA->originalSequenNumber - 1; i++)
                {
                    int pos = i + donors.size() + 1;
                    nA->sequences[i + donors.size() + 1] = oldSequences[i + 1];
                    nA->seqsName[i + donors.size() + 1] = oldNames[i + 1];
                }
                
                nA->originalSequenNumber = donors.size() + nA->originalSequenNumber;
                nA->sequenNumber = nA->originalSequenNumber;
                
                delete [] oldSequences;
                delete [] oldNames;

            }
        }
        else 
        {
            exit(122);
        }
    }
    
    
    while (infile.getline(line, 4096, '\n'))
    {
        {
//             vcf_file::snp_entry snp = vcf_file::snp_entry();
            char * tmp;
            
            tmp = std::strtok(line, "\t");
            char * chromosome = new char[strlen(tmp) + 1];
            std::memmove(chromosome, tmp, strlen(tmp) + 1);
            
            tmp = std::strtok(NULL, "\t");
            int position = atoi(tmp);
            
            tmp = std::strtok(NULL, "\t");
            tmp = std::strtok(NULL, "\t");
            
            char * ref = new char[strlen(tmp) + 1];
            std::memmove(ref, tmp, strlen(tmp) + 1);
            
            tmp = std::strtok(NULL, "\t");
            char * alt = new char[strlen(tmp) + 1];
            std::memmove(alt, tmp, strlen(tmp) + 1);
            
            tmp = std::strtok(NULL, "\t");
            float quality = atof(tmp);
            
            tmp = std::strtok(NULL, "\t");
            bool filter = std::strcmp(tmp, "PASS") ? false : true;
            
            
            // INFO
            std::strtok(NULL, "\t");
            // FORMAT
            std::strtok(NULL, "\t");
            // Individues with this SNP
            
            tmp = std::strtok(NULL, "\t");
            
            std::cout   << std::left        << std::setw(15)
                        << chromosome   << std::setw(15)
                        << position     << std::setw(15)
                        << ref          << std::setw(15) 
                        << "~>"             << std::setw(15)
                        << alt          << std::setw(15)
                        << quality      << std::setw(15)
                        << (strlen(ref) == 1 && strlen(alt) == 1) << std::setw(15)
                        << (filter ? "PASS" : "FILTERED");
            
            int counter = 1;
            if (filter && strlen(ref) == 1 && strlen(alt) == 1)
            {
                int i;
                for (i = 0; i < sequences.size(); i++)
                {
                    if (sequences[i] == chromosome)
                        break;
                }
                if (i == sequences.size())
                {
                    std::cout << "Not Found" << std::endl;
                }
                else
                    while (tmp != NULL)
                    {
                        if (strlen(tmp) > 0)
                        {
                            std::cout << " " << donors[counter - 1];
                            
                            if (sources[i]->sequences[counter][position - 1] == ref[0])
                            {
                                sources[i]->sequences[counter][position - 1] = alt[0];
                                std::cout << "~";
                            }
                            else
                            {
                                std::cout << "¬";
                            }
                            
                        }
                        
                        counter ++;
                        tmp = std::strtok(NULL, "\t");
                    }
            }
            
            counter = 0;
            while (tmp != NULL)
            {
                if (!strcmp(tmp, "0"))
                {
                    std::cout << " " << donors[counter];
                }
                
                counter ++;
                tmp = std::strtok(NULL, "\t");
            }
            
            std::cout << std::endl;
            
            delete [] chromosome;
            delete [] ref;
            delete [] alt;
        }
    }
    
    std::cout << "~~> Sequences in FASTA" << std::endl;
    for (newAlignment * A : sources)
    {
        std::cout << A->seqsName[0] << std::endl;
        
        for (int X = 0; X < A->sequenNumber; X++)
        {
            std::cout << "\t" << A->seqsName[X] << std::endl;
            std::cout << "\t" << A->sequences[X].substr(0, 50) << std::endl;
        }
    }
    
//     std::cout << "~~>  Sequences in VCF" << std::endl;
//     for (std::string fname : sequences)
//     {
//         std::cout << fname << std::endl;
//     }
//     
//     std::cout << "~~>  Donors" << std::endl;
//     for (std::string fname : donors)
//     {
//         std::cout << fname << std::endl;
//     }
    
    delete [] line;
}
