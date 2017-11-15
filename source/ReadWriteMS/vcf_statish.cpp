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

void vcf_statish::readVCF(std::vector<newAlignment*> sources, std::vector<std::string> filenames, float minQuality, float minCoverage, bool ignoreFilter)
{
    // All donors present in the vcf files.
    std::vector<std::string> donors = std::vector<std::string>();
    // Contigs present on the VCF files. Each entry must be present in the sources vector.
    std::vector<std::string> contigs = std::vector<std::string>();
    // Position on the reference files of each donor.
    std::vector<std::vector<int>> donorsPositions = std::vector<std::vector<int>>();
    
    
    char * line = new char [4096];
    
    // Obtain contigs and donors from all VCF files.
    for (int C = 0; C < filenames.size(); C++)
    {
        std::string& filename = filenames[C];
        donorsPositions.push_back(std::vector<int>());
        
        std::ifstream infile;
        infile.open(filename);

        if (!infile.is_open())
        {
            debug.report(ErrorCode::CantOpenFile, &filename[0]);
        }
        
        // Read file to get all the donors and contigs present on the VCF
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
                    
                    // We only want the contig fields
                    if (!strcmp(field_name, "contig"))
                    {           
                        std::strtok(NULL, "=");
                        char * field_info = std::strtok(NULL, ",");
                        char * fname = new char [std::strlen(field_info) + 1];
                        memmove(fname, field_info, strlen(field_info));
                        fname[std::strlen(field_info)] = '\0';
                        
                        int U;
                        // Check if the contig has already been added.
                        for (U = 0; U < contigs.size(); U++)
                        {
                            if (contigs[U] == fname)
                            {
                                break;
                            }
                        }
                        
                        // If not, add it to the vector.
                        if (U == contigs.size())
                        {
                            contigs.push_back(fname);
                        }
                        
                        
                        delete[] fname;
                    }
                }
                
                else
                {
                    // We only want to parse the FORMAT line, which contains the donors.
                    std::strtok(strstr(line, "FORMAT"), "\t");
                    
                    char * token = std::strtok(NULL, "\t");
                    while (token != NULL)
                    {
                        char * fname = new char [std::strlen(token) + 1];
                        memmove (fname, token, strlen (token));
                        fname[std::strlen(token)] = '\0';
                        
                        int U;
                        
                        // Check every other donor already added.
                        for (U = 0; U < donors.size(); U++)
                        {
                            if (donors[U] == fname)
                            {
                                break;
                            }
                        }
                        
                        // If not present, we add it
                        if (U == donors.size())
                        {
                            donorsPositions[C].push_back(U);
                            donors.push_back(fname);
                        }
                        
                        // If already present, warn the user.
                        else
                        {
                            donorsPositions[C].push_back(U);
                            debug.report(WarningCode::DonorAlreadyAdded, &fname[0]);
                        }
                        
                        token = std::strtok(NULL, "\t");
                        delete[] fname;
                    }
                    break;
                }
            }
        }
        
        infile.close();
    }
    
    // Extend the files.
    {
        // Check if all contigs have a reference alignment.
        bool checkIn = true;
        for (int x = 0; x < contigs.size(); x++)
        {
            int i;
            for (i = 0; i < sources.size(); i++)
            {
                if (!strcmp(&contigs[x][0], &sources[i]->seqsName[0][0]) )
                    break;
            }
            
            if (i == sources.size())
            {
                debug.report(ErrorCode::NoReferenceSequenceForContig, &contigs[x][0]);
                checkIn = false;
            }
        }
        
        // Extend the reference files with new sequences from donors.
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
    }
    
    for (int C = 0; C < filenames.size(); C++)
    {
//         std::cout << filenames[C] << std::endl;
        
        std::ifstream infile;
        infile.open(filenames[C]);
        if (!infile.is_open())
        {
            debug.report(ErrorCode::CantOpenFile, &filenames[C][0]);
        }
        
        while (infile.getline(line, 4096, '\n'))
        {
            if (line[0] == '#') continue;
            {
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
                tmp = std::strtok(NULL, "\t");
                
                // FORMAT
                tmp = std::strtok(NULL, "\t");
                // Individues with this SNP
                char * format = new char[strlen(tmp) + 1];
                strcpy(format, tmp);
                std:vector<std::string> format_fields = std::vector<std::string>();
                
                char * tmp_format = strtok(format, ":");
                while (tmp_format != NULL)
                {
                    if (strlen(tmp_format) > 1 && tmp_format[0] == 'D' && tmp_format[1] == 'P')
                    {
//                         std::cout << tmp_format << std::endl;
                    }
                    tmp_format = strtok(NULL, ":");
                }
                
                delete [] format;
                
                tmp = std::strtok(line, "\t");
                
//                 std::cout   << std::left    << std::setw(15)
//                             << chromosome   << std::setw(15)
//                             << position     << std::setw(15)
//                             << ref          << std::setw(15) 
//                             << "~>"         << std::setw(15)
//                             << alt          << std::setw(15)
//                             << quality      << std::setw(15)
//                             << (strlen(ref) == 1 && strlen(alt) == 1) << std::setw(15)
//                             << (filter ? "PASS" : "FILTERED");
                
                int counter = 0;
                if (/*filter && quality > minQuality && strlen(ref) == 1 && strlen(alt) == 1*/ true)
                {
                    bool canPass = true;
                    if (!filter)
                    {
//                         std::cout << " Previously Filtered.";
                        canPass = false;
                    }
                    
                    if (quality < minQuality)
                    {
//                         std::cout << " Filtered by Quality.";
                        canPass = false;
                    }
                    
                    if (strlen(ref) != 1)
                    {
//                         std::cout << " Reference has more than 1 nucleotide.";
                        canPass = false;
                    }

                    if (canPass) 
                    {
                        
                        int i;
                        for (i = 0; i < contigs.size(); i++)
                        {
                            if (contigs[i] == chromosome)
                                break;
                        }
                        
                        if (i == contigs.size())
                        {
                            std::cout << "Not Found" << std::endl;
//                             canPass = false;
                        }
                        
                        else if (sources[i]->sequences[donorsPositions[C][counter] + 1].size() <= position - 1)
                        {
//                             if (donorsPositions[C].size() <= counter) exit(500);
                            debug.report(ErrorCode::SNPoutOfBounds, 
                                        new std::string[3] { 
                                                std::to_string(position), 
                                                filenames[C],
                                                std::to_string(sources[i]->sequences[donorsPositions[C][counter] + 1].size()) });
//                             canPass = false;
                        }
                        
                        else
                        {
                            counter = 0;
                            while (tmp != NULL)
                            {
                                if (strlen(tmp) > 0)
                                {
    //                                 std::cout << " " << donors[counter];
                                    
                                    if (sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] == ref[0] || 
                                        std::toupper(sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1]) == std::toupper(ref[0]))
                                    {
                                        int curval = 0;
                                        if (strlen(alt) > 1)
                                        {
                                            int c, maxlen;
                                            for (c = 0, maxlen = strlen(alt); c < maxlen; c++)
                                            {
                                                switch(alt[c])
                                                {
                                                    case 'A': curval |= 1 << 0; break; // = 1
                                                    case 'C': curval |= 1 << 1; break; // = 2
                                                    case 'T': curval |= 1 << 2; break; // = 4
                                                    case 'G': curval |= 1 << 3; break; // = 8
                                                    default: break;
                                                }
                                                if (++c < maxlen && c == ',')
                                                    continue;
                                            }
                                            if (c == maxlen + 1)
                                            {
                                                switch(curval)
                                                {
                                                    // One Base
    //                                              case 1://A = 1
    //                                              case 2://C = 2
    //                                              case 4://T = 4
    //                                              case 8://G = 8
                                                    
                                                    // Two Bases
                                                    case 3:  sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = 'M'; break; // aMine AC
                                                    case 5:  sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = 'W'; break; // Weak AT
                                                    case 9:  sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = 'R'; break; // puRine AG
                                                    case 6:  sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = 'Y'; break; // pYrimidine CT
                                                    case 10: sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = 'S'; break; // Strong CG
                                                    case 12: sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = 'K'; break; // Keto TG
                                                    
                                                    // Three Bases
                                                    case 14: sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = 'B'; break; // CTG not A
                                                    case 13: sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = 'D'; break; // ATG not C
                                                    case 11: sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = 'V'; break; // ACG not T
                                                    case 7:  sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = 'H'; break; // ACT not G
                                                    
                                                    // All Four Bases
                                                    case 15: sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = 'N'; break; // ACTG
                                                    default:break;
                                                }
    //                                             std::cout << "#";
                                            }
//                                             else
//                                             {
//                                                 std::cout << "¬";
//                                             }
                                        }
                                        else
                                        {
                                            sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = alt[0];
    //                                         std::cout << "~";
                                        }
    //                                     std::cout << "!";
                                    }
                                    else
                                    {
                                        debug.report(WarningCode::ReferenceNucleotideNotCorresponding, 
                                                    new std::string[5]{ 
                                                        sources[i]->seqsName[0], 
                                                        std::to_string(position),
                                                        filenames[C],
                                                        std::string(1, sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1]),
                                                        ref
                                                    } );
    //                                     std::cout << "¬";
                                    } 
                                }
                                counter ++;
                                tmp = std::strtok(NULL, "\t");
                            }
                            counter = 0;
                            while (tmp != NULL)
                            {
                                if (!strcmp(tmp, "0"))
                                {
    //                                 std::cout << " " << donors[counter];
                                }
                                counter ++;
                                tmp = std::strtok(NULL, "\t");
                            }
                    }
                        }
//                     std::cout << std::endl;
                }
                delete [] chromosome;
                delete [] ref;
                delete [] alt;
            }
        }
    
    }

    for (newAlignment * A : sources)
    {
        std::cout << A->seqsName[0] << std::endl;
        
        for (int X = 0; X < A->sequenNumber; X++)
        {
            std::cout << "\t" << A->seqsName[X] << std::endl;
            std::cout << "\t" << A->sequences[X].substr(0, 50) << std::endl;
        }
    }
    delete [] line;
}


