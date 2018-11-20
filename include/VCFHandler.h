#ifndef VCF_STATISH_H
#define VCF_STATISH_H

#include "../include/reportsystem.h"
#include "../include/newAlignment.h"

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <vector>

namespace ngs {
    
    namespace IUPAC
    {
        static int getTagFromCharArray ( char * array, char separator ) {
            int c, maxlen, curval = 0;
            for ( c = 0, maxlen = strlen ( array ); c < maxlen; c++ ) {
                switch ( array[c] ) {
                case 'A':
                    curval |= 1 << 0;
                    break; // = 1
                case 'C':
                    curval |= 1 << 1;
                    break; // = 2
                case 'T':
                    curval |= 1 << 2;
                    break; // = 4
                case 'G':
                    curval |= 1 << 3;
                    break; // = 8
                default:
                    break;
                }
                if ( ++c < maxlen && c == separator )
                    continue;
            }
            if ( c == maxlen + 1 )
                return curval;
            return -1;
        }

        static char getCharFromTag ( int tag ) {
            switch ( tag ) {

            // One Base
            // case 1://A = 1
            // case 2://C = 2
            // case 4://T = 4
            // case 8://G = 8

            // Two Bases
            case 3:
                return 'M';// aMine AC
            case 5:
                return 'W'; // Weak AT
            case 9:
                return 'R'; // puRine AG
            case 6:
                return 'Y';// pYrimidine CT
            case 10:
                return 'S';// Strong CG
            case 12:
                return 'K';// Keto TG

            // Three Bases
            case 14:
                return 'B';// CTG not A
            case 13:
                return 'D';// ATG not C
            case 11:
                return 'V';// ACG not T
            case 7:
                return 'H';// ACT not G

            // All Four Bases
            default:
            case 15:
                return 'N';// ACTG
            }
        }

    }
    
    namespace __internal
    {
        static void printApeek ( std::vector<newAlignment *> & sources ) {
            for ( newAlignment * A : sources ) {
                std::cout << A->seqsName[0] << std::endl;

                for ( int X = 0; X < A->sequenNumber; X++ ) {
                    std::cout << "\t>" << A->seqsName[X] << std::endl;
                    std::cout << "\t" << A->sequences[X].substr ( 0, 50 ) << std::endl;
                }
            }
        }

        static void extendAlignments (
            std::vector<newAlignment*>  & sources,
            std::vector<std::string>    & contigs,
            std::vector<std::string>    & donors ) {
            // Extend the files.
            {
                // Check if all contigs have a reference alignment.
                bool checkIn = true;
                for ( int x = 0; x < contigs.size(); x++ ) {
                    int i;
                    for ( i = 0; i < sources.size(); i++ ) {
                        if ( !strcmp ( &contigs[x][0], &sources[i]->seqsName[0][0] ) )
                            break;
                    }

                    if ( i == sources.size() ) {
                        debug.report ( ErrorCode::NoReferenceSequenceForContig, &contigs[x][0] );
                        checkIn = false;
                    }
                }

                // Extend the reference files with new sequences from donors.
                if ( checkIn ) {
                    for ( newAlignment * nA : sources ) {

                        std::string * oldSequences  = nA->sequences;
                        std::string * oldNames      = nA->seqsName;

                        std::string seq = nA->sequences[0];
                        std::string nam = nA->seqsName[0];

                        nA -> sequences = new std::string[donors.size() + nA->originalSequenNumber];
                        nA -> sequences[0] = seq;

                        nA -> seqsName = new std::string[donors.size() + nA->originalSequenNumber];
                        nA -> seqsName[0] = nam;

                        for ( int i = 1; i < donors.size() + 1; i++ ) {
                            nA->sequences[i] = seq.c_str();
                            nA->seqsName[i] = donors[i - 1];
                        }

                        for ( int i = 0; i < nA->originalSequenNumber - 1; i++ ) {
                            int pos = i + donors.size() + 1;
                            nA->sequences[pos] = oldSequences[i + 1];
                            nA->seqsName[pos] = oldNames[i + 1];
                        }

                        nA->originalSequenNumber = donors.size() + nA->originalSequenNumber;
                        nA->sequenNumber = nA->originalSequenNumber;

                        delete [] oldSequences;
                        delete [] oldNames;

                    }
                }
            }
        }

        static void applyVariantCallingFiles (
            std::vector<newAlignment*>      & sources ,
            std::vector<std::string>        & filenames,
            std::vector<std::string>        & contigs,
            std::vector<std::vector<int>>   & donorsPositions,
            float minQuality,
            float minCoverage,
            bool ignoreFilter,
            const char * const replacementChar ) {
            char * line = new char [4096];
            for ( int C = 0; C < filenames.size(); C++ ) {
                std::ifstream infile;
                infile.open ( filenames[C] );
                if ( !infile.is_open() ) {
                    debug.report ( ErrorCode::CantOpenFile, &filenames[C][0] );
                }

                std::vector<std::string> donorsInfo = std::vector<std::string> ( donorsPositions[C].size() );

                while ( infile.getline ( line, 4096, '\n' ) ) {
                    if ( line[0] == '#' ) continue;
                    {
                        char * tmp;

                        // Contig
                        tmp = std::strtok ( line, "\t" );
                        char * contig = new char[strlen ( tmp ) + 1];
                        std::memmove ( contig, tmp, strlen ( tmp ) + 1 );

                        // Position
                        tmp = std::strtok ( nullptr, "\t" );
                        int position = atoi ( tmp );

                        // ID
                        std::strtok ( nullptr, "\t" );
                        // Ref
                        tmp = std::strtok ( nullptr, "\t" );
                        char * ref = new char[strlen ( tmp ) + 1];
                        std::memmove ( ref, tmp, strlen ( tmp ) + 1 );

                        // Alt
                        tmp = std::strtok ( nullptr, "\t" );
                        char * alt = new char[strlen ( tmp ) + 1];
                        std::memmove ( alt, tmp, strlen ( tmp ) + 1 );

                        // Quality
                        tmp = std::strtok ( nullptr, "\t" );
                        float quality = atof ( tmp );

                        // Filter
                        tmp = std::strtok ( nullptr, "\t" );
                        bool filter = std::strcmp ( tmp, "PASS" ) ? false : true;

                        // Info
                        tmp = std::strtok ( nullptr, "\t" );

                        // Format
                        tmp = std::strtok ( nullptr, "\t" );
                        char * format = new char[strlen ( tmp ) + 1];
                        std::memmove ( format, tmp, strlen ( tmp ) + 1 );

                        // Donors
                        int counter = 0;
                        tmp = std::strtok ( nullptr, "\t" );
                        while ( tmp != nullptr ) {
                            donorsInfo[counter++] = tmp;
                            tmp = std::strtok ( nullptr, "\t" );
                        }

                        // Format -> Read Depth Index
                        counter = 0;
                        tmp = std::strtok ( format, ":" );
                        int readDepthIndex = -1;
                        while ( tmp != nullptr ) {
                            if ( strlen ( tmp ) > 1 && tmp[0] == 'D' && tmp[1] == 'P' ) {
                                readDepthIndex = counter;
                                break;
                            }
                            tmp = std::strtok ( nullptr, ":" );
                            counter++;
                        }

                        delete [] format;

                        counter = 0;
                        if ( /*filter && quality > minQuality && strlen(ref) == 1 && strlen(alt) == 1*/ true ) {
                            bool canPass = true;
                            if ( !ignoreFilter && !filter ) {
                                canPass = false;
                            }

                            if ( quality < minQuality ) {
                                canPass = false;
                            }

                            if ( strlen ( ref ) != 1 ) {
                                canPass = false;
                            }

                            if ( canPass ) {

                                int i;
                                for ( i = 0; i < contigs.size(); i++ ) {
                                    if ( contigs[i] == contig )
                                        break;
                                }

                                if ( i == contigs.size() ) {
                                    std::cout << "Contig not Found" << std::endl;
                                }

                                else if ( sources[i]->sequences[donorsPositions[C][counter] + 1].size() <= position - 1 ) {
                                    debug.report ( ErrorCode::SNPoutOfBounds, new std::string[3] {
                                                       std::to_string ( position ),
                                                       filenames[C],
                                                       std::to_string ( sources[i]->sequences[donorsPositions[C][counter] + 1].size() )
                                                   } );
                                }

                                else {
                                    counter = 0;
                                    for ( std::string& s : donorsInfo ) {
                                        if ( s != "0" ) {
                                            char * tmp = std::strtok ( &s[0], ":" );
                                            for ( int V = 0; V < readDepthIndex; V++ )
                                                tmp = std::strtok ( nullptr, ":" );
                                            int readDepth = std::stoi ( tmp );

                                            if ( readDepth < minCoverage ) {
                                                if ( replacementChar )
                                                    sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = *replacementChar;
                                                break;
                                            }

                                            if ( sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] == ref[0] ||
                                                    std::toupper ( sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] ) == std::toupper ( ref[0] ) ) {
                                                int curval = 0;
                                                if ( strlen ( alt ) > 1 ) {
                                                    int curval = ngs::IUPAC::getTagFromCharArray ( alt, ',' );
                                                    if ( curval != -1 ) {
                                                        sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = ngs::IUPAC::getCharFromTag ( curval );
                                                    }
                                                } else {
                                                    sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = alt[0];
                                                }
                                            } else {
                                                debug.report ( WarningCode::ReferenceNucleotideNotCorresponding,
                                                               new std::string[5] {
                                                                   sources[i]->seqsName[0],
                                                                   std::to_string ( position ),
                                                                   filenames[C],
                                                                   std::string ( 1, sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] ),
                                                                   ref
                                                               } );
                                            }
                                        }
                                        counter ++;
                                        tmp = std::strtok ( nullptr, "\t" );
                                    }
                                }
                            } else {
                                int i;
                                for ( i = 0; i < contigs.size(); i++ ) {
                                    if ( contigs[i] == contig )
                                        break;
                                }

                                if ( i == contigs.size() ) {
                                    std::cout << "Contig not Found" << std::endl;
                                }

                                if ( replacementChar )
                                    sources[i]->sequences[donorsPositions[C][counter] + 1][position - 1] = *replacementChar;
                            }
                        }
                        delete [] contig;
                        delete [] ref;
                        delete [] alt;
                    }
                }
            }
        }

        static void obtainContigsAndDonors (
            std::vector<std::string> & filenames,
            std::vector<std::string> & donors,
            std::vector<std::string> & contigs,
            std::vector<std::vector<int>> & donorsPositions ) {
            char * line = new char [4096];
            // Obtain contigs and donors from all VCF files.
            for ( int C = 0; C < filenames.size(); C++ ) {
                std::string& filename = filenames[C];
                donorsPositions.push_back ( std::vector<int>() );

                std::ifstream infile;
                infile.open ( filename );

                if ( !infile.is_open() ) {
                    debug.report ( ErrorCode::CantOpenFile, &filename[0] );
                }

                // Read file to get all the donors and contigs present on the VCF
                while ( infile.getline ( line, 4096, '\n' ) ) {
                    if ( line[0] == '#' ) {
                        if ( line[1] == '#' ) {
                            // Remove first two characters;
                            memmove ( line, line+2, strlen ( line+2 ) + 2 );
                            // Print result
                            char * field_name = std::strtok ( line, "=" );

                            // We only want the contig fields
                            if ( !strcmp ( field_name, "contig" ) ) {
                                std::strtok ( nullptr, "=" );
                                char * field_info = std::strtok ( nullptr, "," );
                                char * fname = new char [std::strlen ( field_info ) + 1];
                                memmove ( fname, field_info, strlen ( field_info ) );
                                fname[std::strlen ( field_info )] = '\0';

                                int U;
                                // Check if the contig has already been added.
                                for ( U = 0; U < contigs.size(); U++ ) {
                                    if ( contigs[U] == fname ) {
                                        break;
                                    }
                                }

                                // If not, add it to the vector.
                                if ( U == contigs.size() ) {
                                    contigs.push_back ( fname );
                                }


                                delete[] fname;
                            }
                        }

                        else {
                            // We only want to parse the FORMAT line, which contains the donors.
                            std::strtok ( strstr ( line, "FORMAT" ), "\t" );

                            char * token = std::strtok ( nullptr, "\t" );
                            while ( token != nullptr ) {
                                char * fname = new char [std::strlen ( token ) + 1];
                                memmove ( fname, token, strlen ( token ) );
                                fname[std::strlen ( token )] = '\0';

                                int U;

                                // Check every other donor already added.
                                for ( U = 0; U < donors.size(); U++ ) {
                                    if ( donors[U] == fname ) {
                                        break;
                                    }
                                }

                                // If not present, we add it
                                if ( U == donors.size() ) {
                                    donorsPositions[C].push_back ( U );
                                    donors.push_back ( fname );
                                }

                                // If already present, warn the user.
                                else {
                                    donorsPositions[C].push_back ( U );
                                    debug.report ( WarningCode::DonorAlreadyAdded, &fname[0] );
                                }

                                token = std::strtok ( nullptr, "\t" );
                                delete[] fname;
                            }
                            break;
                        }
                    }
                }

                infile.close();
            }

            delete [] line;
        }


    }

    static void readVCF(
            std::vector<newAlignment *> sources,
            std::vector<std::string> filenames,
            float minQuality,
            float minCoverage,
            bool ignoreFilter,
            const char *const replacementChar) {

        // All donors present in the vcf files.
        std::vector<std::string> donors = std::vector<std::string>();
        // Contigs present on the VCF files. Each entry must be present in the sources vector.
        std::vector<std::string> contigs = std::vector<std::string>();
        // Position on the reference files of each donor.
        std::vector<std::vector<int>> donorsPositions = std::vector<std::vector<int>>();

        ngs::__internal::obtainContigsAndDonors ( filenames, donors, contigs, donorsPositions );

        ngs::__internal::extendAlignments ( sources, contigs, donors );

        ngs::__internal::applyVariantCallingFiles (
                sources , filenames, contigs, donorsPositions,
                minQuality, minCoverage, ignoreFilter, replacementChar );
    }
}

#endif // VCF_STATISH_H
