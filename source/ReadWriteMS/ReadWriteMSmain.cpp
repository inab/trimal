#include "../../include/newAlignment.h"
#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include <string>
#ifndef debug
#define debug false
#endif

int parseArguments(int argc, char *argv[], ReadWriteMS* machine, std::vector<std::string>* inFiles, std::vector<std::string>* outFormats, std::string* outPattern)
{
    for(int i = 1; i < argc; i++ )
    {
        // Check if current argument is the '-in' argument.
        if (!strcmp(argv[i], "-in"))
        {
            // Check if this is the last argument.
            if (i >= argc -1)
            {
                cerr << "At least one file should be passed after the '-in' argument" << endl;
                return 1;
            }
            
            // Check if the next argument is a parameter or file argument.
            else if (argv[i + 1][0] == '-')
            {
                cerr << "At least one file should be passed after the '-in' argument and you passed argument " << argv[i + 1] << endl;
                return 1;
            }
            
            // Add every file that is not a parameter
            else while(++i != argc)
            {
                if (argv[i][0] == '-')
                {
                    i--;
                    break;
                }
                else
                {
                    inFiles->push_back(argv[i]);
                }
            }
        }
        
        // Check if current argument is the '-out' argument.
        else if (!strcmp(argv[i], "-out"))
        {
            // Check if this is the last argument.
            if (i >= argc -1)
            {
                cerr << "A file pattern should be passed after the '-out' argument" << endl;
                return 1;
            }
            else 
            {
                if (argv[i + 1][0] == '-')
                {
                    cerr << "A file pattern should be passed after the '-out' argument and you passed argument " << argv[i + 1] << endl;
                    return 1;
                }
                *outPattern = argv[++i];
                continue;
            }
        }
        
        // Check if current argument is the '-formats' argument.
        else if (!strcmp(argv[i], "-formats"))
        {
            if (i >= argc -1)
            {
                cerr << "A format should be passed after the '-formats' argument" << endl;
                return 1;
            }
            else while(++i != argc)
            {
                if (argv[i][0] == '-')
                {
                    i--;
                    break;
                }
                else
                {
                    outFormats->push_back(argv[i]);
                }
            }
        }
        
        else if (!strcmp(argv[i], "-reverse"))
        {
            machine->reverse = true;
        }
        
        //Compatibility with legacy options:
        else if (argv[i][0] == '-')
        {
            if (!strcmp(argv[i], "-html"))
                outFormats->push_back("html");
            
            else if (!strcmp(argv[i], "-nbrf"))
                outFormats->push_back("nbrf");
            
            else if (!strcmp(argv[i], "-mega"))
                outFormats->push_back("mega");
            
            else if (!strcmp(argv[i], "-nexus"))
                outFormats->push_back("nexus");
            
            else if (!strcmp(argv[i], "-clustal"))
                outFormats->push_back("clustal");
            
            else if (!strcmp(argv[i], "-fasta") || !strcmp(argv[i], "-onlyseqs"))
                outFormats->push_back("fasta");
            
            else if (!strcmp(argv[i], "-fasta_m10"))
            {
                outFormats->push_back("fasta");
                machine->shortNames = true;
            }
            
            else if (!strcmp(argv[i], "-phylip"))
                outFormats->push_back("phylip40");
            
            else if (!strcmp(argv[i], "-phylip_m10"))
            {
                outFormats->push_back("phylip40");
                machine->shortNames = true;
            }
            
            else if (!strcmp(argv[i], "-phylip_paml"))
                outFormats->push_back("phylippaml");
            
            else if (!strcmp(argv[i], "-phylip_paml_m10"))
            {
                outFormats->push_back("phylippaml");
                machine->shortNames = true;
            }
            
            else if (!strcmp(argv[i], "-phylip3.2"))
                outFormats->push_back("phylip32");
            
            else if (!strcmp(argv[i], "-phylip3.2_m10"))
            {
                outFormats->push_back("phylip32");
                machine->shortNames = true;
            }
            else if (!strcmp(argv[i], "-format"))
            {
                machine->format = true;
            }
            else if (!strcmp(argv[i], "-type"))
            {
                machine->type = true;
            }
            else if (!strcmp(argv[i], "-info"))
            {
                machine->info = true;
            }
        }
        
        // If a command is not recognized, give an error.
        else
        {
            cerr << argv[i] << " not recognized." << endl;
            return 1;
        }
    }

#if debug
    cout << "Input Files:" << endl;
    for (std::string ifile : *inFiles)
    {
        cout << "-> Input file: " << ifile << endl;
    }
    
    cout << "Out Formats:" << endl;
    for (std::string oformat : *outFormats)
    {
        cout << "-> Output format: " << oformat << endl;
    }
    
    cout << "Under the pattern\n-> " << *outPattern << endl;
#endif
    
    return 0;
}

int checkArguments(ReadWriteMS* machine, std::vector<std::string>* inFiles, std::vector<std::string>* outFormats, std::string* outPattern)
{
    int returnValue = 0;
    if (inFiles->size() == 0)
    {
        cerr        << "At least one input file must be provided" << endl;
        returnValue = 1;
    }
    // LEGACY OPTIONS
    if (machine->format || machine->info || machine->type)
    {
        if (inFiles->size() != 1)
        {
            cerr    << "Print alignment info from multiple files is not allowed" << endl;
            returnValue = 1;
        }
    }
    if (*outPattern == "")
    {
        if (inFiles->size() == 1 && outFormats->size() == 1 && !(machine->format || machine->info || machine->type))
            machine->hasOutputFile = false;
        else if (inFiles->size() > 1)
        {
            cerr    << "Terminal output option disabled when there are more than one input files." << endl;
            returnValue = 1;
        }
        else if (outFormats->size() > 1)
        {
            cerr    << "Terminal output option disabled when there are more than one output formats." << endl;
            returnValue = 1;
        }
        else if (outFormats->size() != 0)
        {
            cerr    << "Terminal output option not compatible with information printing (-info | -format | -type)." << endl 
                    << "Provide an output format or disable information printing." << endl;
            returnValue = 1;
        }
    }
    else if (outFormats->size() == 0)
    {
        cerr        << "At least one output format must be provided" << endl;
        returnValue = 1;
    }
    return returnValue;
}

int main(int argc, char *argv[])
{
    ReadWriteMS MachineState = ReadWriteMS();
    
    std::vector<std::string> outFormats = std::vector<std::string>();
    std::vector<std::string> inFiles  = std::vector<std::string>();
    std::string outPattern = "";
    
    int result = parseArguments(argc, argv, &MachineState, &inFiles, &outFormats, &outPattern);
    if (result != 0) return result;
    
    result = checkArguments(&MachineState, &inFiles, &outFormats, &outPattern);
    if (result != 0) return result;
    
    
//     MachineState.shortNames = false;
//     MachineState.keepHeader = false;
//     MachineState.reverse    = false;
    if(MachineState.format || MachineState.info || MachineState.type) {
        newAlignment* alignment = MachineState.loadAlignment(inFiles[0]);
        
        if (MachineState.format)
            /* Inform about if sequences are aligned or not */
            cout    << "## Input file format\t" << MachineState.getInputStateName(inFiles[0]) << endl
                    << "## Input file aligned\t" << (alignment->isAligned ? "YES":"NO")
                    << endl;

        if(MachineState.type) {
            /* Inform about biological datatype */
            if (alignment->getAlignmentType() == DNAType)
                cout << "## Input file datatype\tnucleotides:dna" << endl;
            else if (alignment->getAlignmentType() == DNADeg)
                cout << "## Input file datatype\tnucleotides:dna_degenerate_codes" << endl;
            else if (alignment->getAlignmentType() == RNAType)
                cout << "## Input file datatype\tnucleotides:rna" << endl;
            else if (alignment->getAlignmentType() == RNADeg)
                cout << "## Input file datatype\tnucleotides:rna_degenerate_codes" << endl;
            else if (alignment->getAlignmentType() == AAType)
                cout << "## Input file datatype\tamino-acids" << endl;
            else
                cout << "## Input file datatype\tunknown" << endl;
        }

        if(MachineState.info)
            alignment->printAlignmentInfo(cout);
        
        delete alignment;
    }
    else
        MachineState.processFile(&inFiles, &outPattern, &outFormats);
}
