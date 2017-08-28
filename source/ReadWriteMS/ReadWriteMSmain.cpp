#include "../../include/defines.h"
#include "../../include/newAlignment.h"
#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include <string>
#include "../../include/values.h"

int parseArguments(int argc, char *argv[], ReadWriteMS* machine, std::vector<std::string>* inFiles, std::vector<std::string>* outFormats, std::string* outPattern)
{
    if (argc == 1) return -1;
    for(int i = 1; i < argc; i++ )
    {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help"))
            return -1;
        // Check if current argument is the '-in' argument.
        if (!strcmp(argv[i], "-in"))
        {
            // Check if this is the last argument.
            if (i >= argc -1)
            {
                cerr << "ERROR: At least one file should be passed after the '-in' argument" << endl;
                return 1;
            }
            
            // Check if the next argument is a parameter or file argument.
            else if (argv[i + 1][0] == '-')
            {
                cerr << "ERROR: At least one file should be passed after the '-in' argument and you passed argument " << argv[i + 1] << endl;
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
        
        else if (!strcmp(argv[i], "-shortNames"))
        {
            machine->shortNames = true;
        }
        else if (!strcmp(argv[i], "-keepHeaders"))
        {
            machine->keepHeader = true;
        }
        
        //Compatibility with legacy options:

        else if (!strcmp(argv[i], "-html"))
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

        // If a command is not recognized, give an error.
        else
        {
            cerr << argv[i] << " not recognized or repeated." << endl;
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
        cerr        << "ERROR: At least one input file must be provided" << endl;
        returnValue = 1;
    }
    if (*outPattern == "")
    {
        if (inFiles->size() == 1 && outFormats->size() == 1 && !(machine->format || machine->info || machine->type))
            machine->hasOutputFile = false;
        else if (outFormats->size() != 0)
        {
            cerr    << "ERROR: Terminal output option not compatible with information printing (-info | -format | -type)" << endl 
                    << "Provide an output format or disable information printing." << endl;
            returnValue = 1;
        }
    }
    else if (outFormats->size() == 0)
    {
        cerr        << "ERROR: At least one output format must be provided" << endl;
        returnValue = 1;
    }
    return returnValue;
}

void menu()
{

    cout << endl
    << "readAl v" << VERSION << ".rev" << REVISION << " build[" << BUILD
    << "]. " << AUTHORS << endl << endl

    << "readAl webpage: http://trimal.cgenomics.org" << endl << endl

    << "This program is free software: you can redistribute it and/or modify "
    << endl
    << "it under the terms of the GNU General Public License as published by "
    << endl
    << "the Free Software Foundation, the last available version." << endl
    << endl

    << "BASIC USAGE" << endl << endl
    << "\treadalMS -in <inputfiles> -out <pattern> -format [formats] [options]." << endl << endl

    << "\t-h                    " << "Show this information." << endl
//     << "\t--version            " << "Show readAl version." << endl << endl

    << "\t-in <inputfiles>      " << "Input files in several formats. Separated by spaces." << endl
    << "\t                      " << "Available formats are: " << ReadWriteMS().getInputFormatsAvailable() << endl 
    << "\t-out <pattern>        " << "Output file name pattern (default STDOUT)." << endl
    << "\t                      " << "It will replace optional the tags [in]        -> Original filename without extension." << endl
    << "\t                      " << "                                  [format]    -> Output's format name" << endl
    << "\t                      " << "                                  [extension] -> Output's extension" << endl
    << endl
    
    << "\t-formats             " << "Formats you want the output to be converted to." << endl
    << "\t                     " << "Available formats are: " << ReadWriteMS().getOutputFormatsAvailable() << endl 
    << "\t                     " << "Being the HTML format not a format itself, but a colored report of the alignment files." << endl << endl
    << "\t-format              " << "Print information about input file format "
    << "and if sequences are aligned or not." << endl

    << "\t-type                " << "Print information about biological "
    << "sequences datatype (e.g. nucleotides:dna, nucleotides:rna, aminoacids, etc)"
    << endl

    << "\t-info                " << "Print information about sequences number, "
    << "average sequence length, max & min sequence length"
    << endl 
    
    << "\t-reverse             " << "Output the reverse of sequences in "
    << "input file." << endl << endl
    
    << "\t-shortNames          " << "Shortens the names so they fit on certain formats" << endl

    << "\t-keepHeaders         " << "Keeps the headers of the original format if it had any" << endl << endl

    
    << "LEGACY OPTIONS" << endl << "Take in mind that this arguments may be discontinued any time."<< endl << endl
    
    << "\t-onlyseqs            " << "Generate output with only residues from "
    << "input file" << endl << endl
    
    << "\t-html                " << "Output residues colored according their "
    << "physicochemical properties. HTML file." << endl << endl


    << "\t-nbrf                " << "Output file in NBRF/PIR format" << endl
    << "\t-mega                " << "Output file in MEGA format" << endl

    << "\t-nexus               " << "Output file in NEXUS format" << endl
    << "\t-clustal             " << "Output file in CLUSTAL format" << endl
    << endl

    << "\t-fasta               " << "Output file in FASTA format" << endl
    << "\t-fasta_m10           " << "Output file in FASTA format. Sequences "
    << "name up to 10 characters." << endl << endl

    << "\t-phylip              " << "Output file in PHYLIP/PHYLIP4 format"
    << endl
    << "\t-phylip_m10          " << "Output file in PHYLIP/PHYLIP4 format. "
    << "Sequences name up to 10 characters." << endl
    << "\t-phylip_paml         " << "Output file in PHYLIP format compatible "
    << "with PAML" << endl
    << "\t-phylip_paml_m10     " << "Output file in PHYLIP format compatible "
    << "with PAML. Sequences name up to 10 characters." << endl
    << "\t-phylip3.2           " << "Output file in PHYLIP3.2 format" << endl
    << "\t-phylip3.2_m10       " << "Output file in PHYLIP3.2 format. Sequences"
    << " name up to 10 characters." << endl << endl
    << "If you specify any m10 format, this will result in all formats having the sequences names shortened as this has the same effect as '-shortNames' argument" << endl << endl
    
    
    << "EXAMPLES OF USE" << endl << endl
    
    << "\treadalMS -in ./dataset/AA1.fas -out ./dataset/[in].output.[extension] -formats clustal" << endl
    << "\t  -> Will produce ./dataset/AA1.output.clw" << endl << endl
    
    << "\treadalMS -in ./dataset/example1.clw -out ./dataset/[in].[format].[extension] -formats fasta phylip32 phylip40" << endl
    << "\t  -> Will produce ./dataset/example1.FASTA.fasta ./dataset/example1.PHYLIP32.phy ./dataset/example1.PHYLIP40.phy" << endl << endl
    
    << "\treadalMS -in ./dataset/example1.clw -out ./dataset/[in]/[format].[extension] -formats fasta phylip32 phylip40" << endl
    << "\t  -> Will produce ./dataset/example1/FASTA.fasta ./dataset/example1/PHYLIP32.phy ./dataset/example1/PHYLIP40.phy" << endl
    << "\t     ONLY if ./dataset/example1/ already exists." << endl << endl 
    
    << "\treadalMS -in ./dataset/AA1.fas ./dataset/AA2.fas -out ./dataset/[in].output.[extension] -formats clustal pir" << endl
    << "\t  -> Will produce ./dataset/AA1.output.clw ./dataset/AA2.output.clw ./dataset/AA1.output.pir ./dataset/AA2.output.pir" << endl << endl
    
    << "\treadalMS -in ./dataset/AA1.fas -format -type -info" << endl
    << "\t  -> Will produce terminal output giving information about AA1.fas alignment file" << endl << endl
    
    << "\treadalMS -in ./dataset/AA1.fas ./dataset/AA2.fas -out ./dataset/[in].output.[extension] -formats html" << endl
    << "\t  -> Will produce ./dataset/AA1.output.html ./dataset/AA2.output.html" << endl
    << "\t     Those files are not indeed reformats of the original alignments, but an HTML colored report of the alignment file." << endl
    
    << endl;
}


int main(int argc, char *argv[])
{
    ReadWriteMS MachineState = ReadWriteMS();
    
    std::vector<std::string> outFormats = std::vector<std::string>();
    std::vector<std::string> inFiles  = std::vector<std::string>();
    std::string outPattern = "";
    
    int result = parseArguments(argc, argv, &MachineState, &inFiles, &outFormats, &outPattern);
    if (result == -1) 
    {
        menu();
        return 0;
    }
    else if (result != 0) return result;

    result = checkArguments(&MachineState, &inFiles, &outFormats, &outPattern);
    if (result != 0) return result;
    

    if(MachineState.format || MachineState.info || MachineState.type) {
        for (string str : inFiles)
        {
            newAlignment* alignment = MachineState.loadAlignment(str);
                cout << "Hello there, LOADED" << endl;
            if (alignment != nullptr)
            {
                cout << "## Alignment File: " << str << endl;
                
                if (MachineState.format)
                    /* Inform about if sequences are aligned or not */
                    cout    << "## Input file format\t" << MachineState.getFileFormatName(str) << endl
                            << "## Input file aligned\t" << (alignment->isAligned ? "YES":"NO")
                            << endl;

                if(MachineState.type) {
                    /* Inform about biological datatype */
                    if (alignment->getAlignmentType() == SequenceTypes::DNA)
                        cout << "## Input file datatype\tnucleotides:dna" << endl;
                    else if (alignment->getAlignmentType() == (SequenceTypes::DNA | SequenceTypes::DEG))
                        cout << "## Input file datatype\tnucleotides:dna_degenerate_codes" << endl;
                    else if (alignment->getAlignmentType() == SequenceTypes::RNA)
                        cout << "## Input file datatype\tnucleotides:rna" << endl;
                    else if (alignment->getAlignmentType() == (SequenceTypes::RNA | SequenceTypes::DEG))
                        cout << "## Input file datatype\tnucleotides:rna_degenerate_codes" << endl;
                    else if (alignment->getAlignmentType() == SequenceTypes::AA)
                        cout << "## Input file datatype\tamino-acids" << endl;
                    else
                        cout << "## Input file datatype\tunknown" << endl;
                }

                if(MachineState.info)
                    alignment->printAlignmentInfo(cout);
                
                cout << endl;
                
                if (outFormats.size() > 0)
                    MachineState.saveAlignment(outPattern, &outFormats, alignment);
                
                delete alignment;
            }
        }
    }
    else if (outFormats.size() != 0)
        MachineState.loadAndSaveMultipleAlignments(&inFiles, &outPattern, &outFormats);
    else
        cerr << "ERROR: An option has to be chosen" << endl;
}
