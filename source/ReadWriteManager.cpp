//
// Created by bioinfo on 5/06/17.
//

#include "../include/ReadWriteManager.h"
#include "../include/newAlignment.h"
#include "../include/defines.h"

using namespace std;

bool ReadWriteManager::alignmentSummaryHTML(char *destFile, int residues, int seqs,
  int *selectedRes, int *selectedSeq, float *consValues) {

    /* Generate an HTML file with a visual summary about which sequences/columns
     * have been selected and which have not */

    int i, j, k, kj, upper, minHTML, maxLongName, *gapsValues;
    string tmpColumn;
    float *simValues;
    bool *res, *seq;
    ofstream file;
    char type;

    /* Allocate some local memory */
    tmpColumn.reserve(_alignment->sequenNumber);

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!_alignment->isAligned) {
        cerr << endl << "ERROR: Sequences are not aligned." << endl << endl;
        return false;
    }

    /* Open output file and check that file pointer is valid */
    file.open(destFile);
    if(!file)
        return false;

    /* Compute maximum sequences name length. */
    maxLongName = 0;
    for(i = 0; i < _alignment->sequenNumber; i++)
        maxLongName = utils::max(maxLongName, _alignment->seqsName[i].size());

    /* Compute HTML blank spaces */
    minHTML = utils::max(25, maxLongName + 10);

    /* Initialize local variables to control which columns/sequences
     * will be kept in the output alignment */
    res = new bool[_alignment->residNumber];
    for(i = 0; i < _alignment->residNumber; i++)
        res[i] = false;

    seq = new bool[_alignment->sequenNumber];
    for(i = 0; i < _alignment->sequenNumber; i++)
        seq[i] = false;

    /* Record which columns/sequences from original alignment
     * have been kept in the final one */
    for(i = 0; i < residues; i++)
        res[selectedRes[i]] = true;
    for(i = 0; i < seqs; i++)
        seq[selectedSeq[i]] = true;

    /* Recover some stats about different scores from current alignment */
    gapsValues = NULL;
    if (_alignment->sgaps != NULL)
        gapsValues = _alignment->sgaps -> getGapsWindow();
    simValues = NULL;
    if (_alignment->scons != NULL)
        simValues = _alignment->scons -> getMdkwVector();

    /* Print HTML header into output file */
    file << "<!DOCTYPE html>" << endl << "<html><head>" << endl << "    <meta "
         << "http-equiv=\"Content-Type\" content=\"text/html;charset=ISO-8859-1\" />"
         << endl << "    <title>trimAl v1.4 Summary</title>" << endl
         << "    <style type=\"text/css\" media=\"all\">" << endl

         << "    #b  { background-color: #3366ff; }\n"
         << "    #r  { background-color: #cc0000; }\n"
         << "    #g  { background-color: #33cc00; }\n"
         << "    #p  { background-color: #ff6666; }\n"
         << "    #m  { background-color: #cc33cc; }\n"
         << "    #o  { background-color: #ff9900; }\n"
         << "    #c  { background-color: #46C7C7; }\n"
         << "    #y  { background-color: #FFFF00; }\n"

         << "    .sel  { background-color: #B9B9B9; }\n"
         << "    .nsel { background-color: #E9E9E9; }\n"

         /* Sets of colors for high-lighting scores intervals */
         << "    .c1   { background-color: #FFFBF2; }\n"
         << "    .c2   { background-color: #FFF8CC; }\n"
         << "    .c3   { background-color: #FAF0BE; }\n"
         << "    .c4   { background-color: #F0EAD6; }\n"
         << "    .c5   { background-color: #F3E5AB; }\n"
         << "    .c6   { background-color: #F4C430; }\n"
         << "    .c7   { background-color: #C2B280; color: white; }\n"
         << "    .c8   { background-color: #DAA520; color: white; }\n"
         << "    .c9   { background-color: #B8860B; color: white; }\n"
         << "    .c10  { background-color: #918151; color: white; }\n"
         << "    .c11  { background-color: #967117; color: white; }\n"
         << "    .c12  { background-color: #6E5411; color: white; }\n"

         /* Other HTML elements */
         << "    </style>\n  </head>\n\n" << "  <body>\n" << "  <pre>" << endl;

    /* Show information about how many sequences/residues have been selected */
    file << "    <span class=sel>Selected Sequences: " << setw(5) << right << seqs
         <<" /Selected Residues: " << setw(7) << right << residues << "</span>"
         << endl << "    <span class=nsel>Deleted Sequences:  " << setw(5) << right
         << _alignment->sequenNumber - seqs << " /Deleted Residues:  " << setw(7) << right
         << _alignment->residNumber - residues << "</span>" << endl;

    /* Print headers for different scores derived from input alignment/s */
    if (gapsValues != NULL)
        file << endl << setw(minHTML) << left << "    Gaps Scores:        "
             << "<span  class=c1>  =0=  </span><span  class=c2> <.001 </span>"
             << "<span  class=c3> <.050 </span><span  class=c4> <.100 </span>"
             << "<span  class=c5> <.150 </span><span  class=c6> <.200 </span>"
             << "<span  class=c7> <.250 </span><span  class=c8> <.350 </span>"
             << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
             << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";

    if (simValues != NULL)
        file << endl << setw(minHTML) << left << "    Similarity Scores:  "
             << "<span  class=c1>  =0=  </span><span  class=c2> <1e-6 </span>"
             << "<span  class=c3> <1e-5 </span><span  class=c4> <1e-4 </span>"
             << "<span  class=c5> <.001 </span><span  class=c6> <.010 </span>"
             << "<span  class=c7> <.100 </span><span  class=c8> <.250 </span>"
             << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
             << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";

    if (consValues != NULL)
        file << endl << setw(minHTML) << left << "    Consistency Scores: "
             << "<span  class=c1>  =0=  </span><span  class=c2> <.001 </span>"
             << "<span  class=c3> <.050 </span><span  class=c4> <.100 </span>"
             << "<span  class=c5> <.150 </span><span  class=c6> <.200 </span>"
             << "<span  class=c7> <.250 </span><span  class=c8> <.350 </span>"
             << "<span  class=c9> <.500 </span><span class=c10> <.750 </span>"
             << "<span class=c11> <1.00 </span><span class=c12>  =1=  </span>";

    if ((gapsValues != NULL) or (simValues == NULL) or (consValues == NULL))
        file << endl;

    /* Print Sequences in block of BLOCK_SIZE */
    for(j = 0, upper = HTMLBLOCKS; j < _alignment->residNumber; j += HTMLBLOCKS, upper += \
    HTMLBLOCKS) {

        /* Print main columns number */
        file << endl << setw(minHTML + 10) << right << (j + 10);
        for(i = j + 20; ((i <= _alignment->residNumber) && (i <= upper)); i += 10)
            file << setw(10) << right << (i);

        /* Print special characters to delimit sequences blocks */
        file << endl << setw(minHTML + 1) << right;
        for(i = j + 1; ((i <= _alignment->residNumber) && (i <= upper)); i++)
            file << (!(i % 10) ? "+" : "=");
        file << endl;

        /* Print sequences name */
        for(i = 0; i < _alignment->sequenNumber; i++) {
            file << "    <span class=" << ((seq[i]) ? "sel>" : "nsel>") << _alignment->seqsName[i]
                 << "</span>" << setw(minHTML - 4 - _alignment->seqsName[i].size()) << right << "";

            /* Print residues corresponding to current sequences block */
            for(k = j; ((k < _alignment->residNumber) && (k < upper)); k++) {
                for(kj = 0, tmpColumn.clear(); kj < _alignment->sequenNumber; kj++)
                    tmpColumn += _alignment->sequences[kj][k];
                /* Determine residue color based on residues across the alig column */
                type = utils::determineColor(_alignment->sequences[i][k], tmpColumn);
                if (type == 'w')
                    file << _alignment->sequences[i][k];
                else
                    file << "<span id=" << type << ">" << _alignment->sequences[i][k] << "</span>";
            }
            file << endl;
        }

        file << endl << setw(minHTML) << left << "    Selected Cols:      ";
        for(k = j; ((k < _alignment->residNumber) && (k < (j + HTMLBLOCKS))); k++)
            file << "<span class=" << (res[k] ? "sel" : "nsel") << "> </span>";
        file << endl;

        /* If there is not any score to print, skip this part of the function */
        if ((gapsValues == NULL) and (simValues == NULL) and (consValues == NULL))
            continue;

        /* Print score colors according to certain predefined thresholds */
        if (gapsValues != NULL) {
            file << endl << setw(minHTML) << left << "    Gaps Scores:        ";
            for(k = j; ((k < _alignment->residNumber) && (k < (j + HTMLBLOCKS))); k++)
                if(gapsValues[k] == 0)
                    file << "<span class=c12> </span>";
                else if(gapsValues[k] == _alignment->sequenNumber)
                    file << "<span class=c1> </span>";
                else if(1 - (float(gapsValues[k])/_alignment->sequenNumber) >= .750)
                    file << "<span class=c11> </span>";
                else if(1 - (float(gapsValues[k])/_alignment->sequenNumber) >= .500)
                    file << "<span class=c10> </span>";
                else if(1 - (float(gapsValues[k])/_alignment->sequenNumber) >= .350)
                    file << "<span  class=c9> </span>";
                else if(1 - (float(gapsValues[k])/_alignment->sequenNumber) >= .250)
                    file << "<span  class=c8> </span>";
                else if(1 - (float(gapsValues[k])/_alignment->sequenNumber) >= .200)
                    file << "<span  class=c7> </span>";
                else if(1 - (float(gapsValues[k])/_alignment->sequenNumber) >= .150)
                    file << "<span  class=c6> </span>";
                else if(1 - (float(gapsValues[k])/_alignment->sequenNumber) >= .100)
                    file << "<span  class=c5> </span>";
                else if(1 - (float(gapsValues[k])/_alignment->sequenNumber) >= .050)
                    file << "<span  class=c4> </span>";
                else if(1 - (float(gapsValues[k])/_alignment->sequenNumber) >= .001)
                    file << "<span  class=c3> </span>";
                else
                    file << "<span  class=c2> </span>";
        }
        if (simValues != NULL) {
            file << endl << setw(minHTML) << left << "    Similarity Scores:  ";
            for(k = j; ((k < _alignment->residNumber) && (k < (j + HTMLBLOCKS))); k++)
                if(simValues[k] == 1)
                    file << "<span class=c12> </span>";
                else if(simValues[k] == 0)
                    file << "<span class=c1> </span>";
                else if(simValues[k] >= .750)
                    file << "<span class=c11> </span>";
                else if(simValues[k] >= .500)
                    file << "<span class=c10> </span>";
                else if(simValues[k] >= .250)
                    file << "<span  class=c9> </span>";
                else if(simValues[k] >= .100)
                    file << "<span  class=c8> </span>";
                else if(simValues[k] >= .010)
                    file << "<span  class=c7> </span>";
                else if(simValues[k] >= .001)
                    file << "<span  class=c6> </span>";
                else if(simValues[k] >= 1e-4)
                    file << "<span  class=c5> </span>";
                else if(simValues[k] >= 1e-5)
                    file << "<span  class=c4> </span>";
                else if(simValues[k] >= 1e-6)
                    file << "<span  class=c3> </span>";
                else
                    file << "<span  class=c2> </span>";
        }
        if (consValues != NULL) {
            file << endl << setw(minHTML) << left << "    Consistency Scores: ";
            for(k = j; ((k < _alignment->residNumber) && (k < (j + HTMLBLOCKS))); k++)
                if(consValues[k] == 1)
                    file << "<span class=c12> </span>";
                else if(consValues[k] == 0)
                    file << "<span class=c1> </span>";
                else if(consValues[k] >= .750)
                    file << "<span class=c11> </span>";
                else if(consValues[k] >= .500)
                    file << "<span class=c10> </span>";
                else if(consValues[k] >= .350)
                    file << "<span  class=c9> </span>";
                else if(consValues[k] >= .250)
                    file << "<span  class=c8> </span>";
                else if(consValues[k] >= .200)
                    file << "<span  class=c7> </span>";
                else if(consValues[k] >= .150)
                    file << "<span  class=c6> </span>";
                else if(consValues[k] >= .100)
                    file << "<span  class=c5> </span>";
                else if(consValues[k] >= .050)
                    file << "<span  class=c4> </span>";
                else if(consValues[k] >= .001)
                    file << "<span  class=c3> </span>";
                else
                    file << "<span  class=c2> </span>";
        }
        file << endl;
    }

    /* Print HTML footer into output file */
    file << "    </pre>" << endl << "  </body>" << endl << "</html>" << endl;

    /* Close output file and deallocate local memory */
    file.close();
    delete [] seq;
    delete [] res;

    return true;
}

bool ReadWriteManager::alignmentColourHTML(ostream &file) {

    int i, j, kj, upper, k = 0, maxLongName = 0;
    string tmpColumn;
    char type;

    /* Allocate some local memory */
    tmpColumn.reserve(_alignment->sequenNumber);

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!_alignment->isAligned) {
        cerr << endl << "ERROR: Sequences are not aligned." << endl << endl;
        return false;
    }

    /* Compute maximum sequences name length */
    maxLongName = 0;
    for(i = 0; i < _alignment->sequenNumber; i++)
        maxLongName = utils::max(maxLongName, _alignment->seqsName[i].size());


    /* Print HTML header into output file */
    file << "<!DOCTYPE html>" << endl << "<html><head>" << endl << "    <meta "
         << "http-equiv=\"Content-Type\" content=\"text/html;charset=ISO-8859-1\" />"
         << endl << "    <title>readAl v1.4</title>" << endl
         << "    <style type=\"text/css\">" << endl
         << "    #b  { background-color: #3366ff; }\n"
         << "    #r  { background-color: #cc0000; }\n"
         << "    #g  { background-color: #33cc00; }\n"
         << "    #p  { background-color: #ff6666; }\n"
         << "    #m  { background-color: #cc33cc; }\n"
         << "    #o  { background-color: #ff9900; }\n"
         << "    #c  { background-color: #46C7C7; }\n"
         << "    #y  { background-color: #FFFF00; }\n"
         << "    </style>\n  </head>\n\n" << "  <body>\n  <pre>" << endl;

    /* Print sequences colored according to CLUSTAL scheme based on
     * physical-chemical properties */
    for(j = 0, upper = HTMLBLOCKS; j < _alignment->residNumber; j += HTMLBLOCKS, upper += \
    HTMLBLOCKS) {

        file << endl;
        /* Print main columns number */
        file << setw(maxLongName + 19) << right << (j + 10);
        for(i = j + 20; ((i <= _alignment->residNumber) && (i <= upper)); i += 10)
            file << setw(10) << right << i;

        /* Print special characters to delimit sequences blocks */
        file << endl << setw(maxLongName + 10);
        for(i = j + 1; ((i <= _alignment->residNumber) && (i <= upper)); i++)
            file << (!(i % 10) ? "+" : "=");

        /* Print sequences themselves */
        for(i = 0; i < _alignment->sequenNumber; i++) {

            /* Print sequences name */
            file << endl << setw(maxLongName + 9) << left << _alignment->seqsName[i];

            /* Print residues corresponding to current sequences block */
            for(k = j; ((k < _alignment->residNumber) && (k < upper)); k++) {
                for(kj = 0, tmpColumn.clear(); kj < _alignment->sequenNumber; kj++)
                    tmpColumn += _alignment->sequences[kj][k];
                /* Determine residue color based on residues across the alig column */
                type = utils::determineColor(_alignment->sequences[i][k], tmpColumn);
                if (type == 'w')
                    file << _alignment->sequences[i][k];
                else
                    file << "<span id=" << type << ">" << _alignment->sequences[i][k] << "</span>";
            }
        }
        file << endl;
    }

    /* Print HTML footer into output file */
    file << "    </pre>" << endl << "  </body>" << endl << "</html>" << endl;

    return true;
}

ReadWriteManager::ReadWriteManager(newAlignment* parent)
{
    _alignment = parent;
//     keepHeader =    false;
//     filename =      "";
//     aligInfo =      "";
    
}

ReadWriteManager::ReadWriteManager(newAlignment* parent, ReadWriteManager* mold)
{
    _alignment =    parent;
//     keepHeader =    mold->keepHeader;
//     filename =      mold->filename;
//     aligInfo =      mold->aligInfo;
}
