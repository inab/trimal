#include "../../include/ReadWriteMS/htmlreport_state.h"

#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include "../../include/defines.h"
#include "../../include/utils.h"

int htmlreport_state::CheckAlignment(std::istream* origin)
{
    return 0;
}

newAlignment* htmlreport_state::LoadAlignment(std::string filename)
{
    return nullptr;
}

bool htmlreport_state::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
    int i, j, kj, upper, k = 0, maxLongName = 0;
    std::string tmpColumn;
    char type;

    /* Allocate some local memory */
    tmpColumn.reserve(alignment->sequenNumber);

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!alignment->isAligned) {
        debug.report(ErrorCode::UnalignedAlignmentToAlignedFormat, new std::string[1] { this->name });
        return false;
    }

    /* Compute maximum sequences name length */
    maxLongName = 0;
    for(i = 0; i < alignment->sequenNumber; i++)
        maxLongName = utils::max(maxLongName, alignment->seqsName[i].size());


    /* Print HTML header into output file */
    *output << "<!DOCTYPE html>\n" 
            << "<html><head>\n" 
            << "    <meta http-equiv=\"Content-Type\" content=\"text/html;charset=ISO-8859-1\" />\n" 
            << "    <title>readAl v1.4</title>\n"
            << "    <style type=\"text/css\">\n"
            << "    #b  { background-color: #3366ff; }\n"
            << "    #r  { background-color: #cc0000; }\n"
            << "    #g  { background-color: #33cc00; }\n"
            << "    #p  { background-color: #ff6666; }\n"
            << "    #m  { background-color: #cc33cc; }\n"
            << "    #o  { background-color: #ff9900; }\n"
            << "    #c  { background-color: #46C7C7; }\n"
            << "    #y  { background-color: #FFFF00; }\n"
            << "    </style>\n  </head>\n\n  <body>\n  <pre>\n";

    /* Print sequences colored according to CLUSTAL scheme based on
     * physical-chemical properties */
    for(j = 0, upper = HTMLBLOCKS; j < alignment->residNumber; j += HTMLBLOCKS, upper += \
    HTMLBLOCKS) {

        *output << "\n";
        /* Print main columns number */
        *output << std::setw(maxLongName + 19) << std::right << (j + 10);
        for(i = j + 20; ((i <= alignment->residNumber) && (i <= upper)); i += 10)
            *output << std::setw(10) << std::right << i;

        /* Print special characters to delimit sequences blocks */
        *output << "\n" << std::setw(maxLongName + 10);
        for(i = j + 1; ((i <= alignment->residNumber) && (i <= upper)); i++)
            *output << (!(i % 10) ? "+" : "=");

        /* Print sequences themselves */
        for(i = 0; i < alignment->sequenNumber; i++) {

            /* Print sequences name */
            *output << "\n" << std::setw(maxLongName + 9) << std::left << alignment->seqsName[i];

            /* Print residues corresponding to current sequences block */
            for(k = j; ((k < alignment->residNumber) && (k < upper)); k++) {
                for(kj = 0, tmpColumn.clear(); kj < alignment->sequenNumber; kj++)
                    tmpColumn += alignment->sequences[kj][k];
                /* Determine residue color based on residues across the alig column */
                type = utils::determineColor(alignment->sequences[i][k], tmpColumn);
                if (type == 'w')
                    *output << alignment->sequences[i][k];
                else
                    *output << "<span id=" << type << ">" << alignment->sequences[i][k] << "</span>";
            }
        }
        *output << "\n";
    }

    /* Print HTML footer into output file */
    *output << "    </pre>\n  </body>\n</html>\n";
    
    return true;
}

bool htmlreport_state::RecognizeOutputFormat(std::string FormatName)
{
    if (ReadWriteBaseState::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "html" ||
           FormatName == "HTML" ||
           FormatName == "htmlreport";
}
