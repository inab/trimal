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

#include <utils.h>

#include "InternalBenchmarker.h"
#include "reportsystem.h"
#include "defines.h"
#include "values.h"
#include "utils.h"

namespace utils {

    void initlVect(int *vector, int tam, int valor) {

        for (int i = 0; i < tam; i++) vector[i] = valor;

    }

    void initlVect(float *vector, int tam, float valor) {

        for (int i = 0; i < tam; i++) vector[i] = valor;
    }

    void copyVect(int *vect1, int *vect2, int tam) {

        for (int i = 0; i < tam; i++) vect2[i] = vect1[i];

    }

    void copyVect(float *vect1, float *vect2, int tam) {

        for (int i = 0; i < tam; i++) vect2[i] = vect1[i];
    }

    int roundToInf(double number) {

        //    return std::floor(number);
        return ((int) number);
    }

    int roundInt(double number) {

        //    return std::round(number);
        return ((int) ((double) number + 0.5));
    }

    int roundToSup(double number) {

        //    return std::ceil(number);
        return ((int) ((double) number + 1.0));
    }

    int max(int x, int y) {

        if (x > y) return x;
        else return y;
    }

    float max(float x, float y) {

        if (x > y) return x;
        else return y;
    }

    double max(double x, double y) {

        if (x > y) return x;
        else return y;
    }

    int min(int x, int y) {

        if (x < y) return x;
        else return y;
    }

    float min(float x, float y) {

        if (x < y) return x;
        else return y;
    }

    double min(double x, double y) {

        if (x < y) return x;
        else return y;
    }

    bool isNumber(char *num) {

        int tam = strlen(num);
        if (tam == 0) return false;
        int i, flt = 1, expn = 1, sgn = 1;

        for (i = 0; i < tam; i++) {
            if (num[i] == '.' && flt)
                flt = 0;

            else if (((num[i] == 'e') || (num[i] == 'E')) && expn)
                expn = 0;

            else if (((num[i] == '+') || (num[i] == '-')) && sgn) {
                if (!expn) sgn = 0;
            } else if (num[i] > '9' || num[i] < '0')
                return false;
        }

        return true;

    }

    bool compare(char *a, char *b) {

        return (!strcmp(a, b));
    }

    void removeSpaces(char *in, char *out) {

        unsigned int i, j = 0;

        for (i = 0; i < strlen(in); i++) {

            if (in[i] != ' ' && in[i] != '\t') {
                out[j] = in[i];
                j++;
            }
        }
        out[j] = '\0';
    }


    void quicksort(float *vect, int ini, int fin) {
        // NOTE Introspective sort seems better - JMaria
        float elem_div;
        int i, j;

        if ((ini >= fin) || (fin < 0))
            return;

        elem_div = vect[fin];
        i = ini - 1;
        j = fin;

        while (true) {

            while (vect[++i] < elem_div)
                if (i == fin)
                    break;

            while (vect[--j] > elem_div)
                if (j == 0)
                    break;

            if (i < j)
                swap(&vect[i], &vect[j]);
            else
                break;
        }

        swap(&vect[i], &vect[fin]);

        quicksort(vect, ini, i - 1);
        quicksort(vect, i + 1, fin);
    }

    void swap(float *a, float *b) {

        float temp;

        temp = *a;
        *a = *b;
        *b = temp;

    }


    void quicksort(int *vect, int ini, int fin) {
        // NOTE Introspective sort seems better - JMaria
        int i, j, elem_div;

        if ((ini >= fin) || (fin < 0))
            return;

        elem_div = vect[fin];
        i = ini - 1;
        j = fin;

        while (true) {

            while (vect[++i] < elem_div)
                if (i == fin)
                    break;

            while (vect[--j] > elem_div)
                if (j == 0)
                    break;

            if (i < j)
                swap(&vect[i], &vect[j]);
            else
                break;
        }

        swap(&vect[i], &vect[fin]);

        quicksort(vect, ini, i - 1);
        quicksort(vect, i + 1, fin);
    }

    void swap(int *a, int *b) {
        int temp;

        temp = *a;
        *a = *b;
        *b = temp;

    }

    void quicksort(int **vect, int ini, int fin) {
        // NOTE Introspective sort seems better - JMaria
        float elem_div;
        int i, j;

        if ((ini >= fin) || (fin < 0))
            return;

        elem_div = vect[fin][0];
        i = ini - 1;
        j = fin;

        while (true) {

            while (vect[++i][0] < elem_div) if (i == fin) break;
            while (vect[--j][0] > elem_div) if (j == 0) break;

            if (i < j) swap(&vect[i], &vect[j]);
            else break;

        }

        swap(&vect[i], &vect[fin]);

        quicksort(vect, ini, i - 1);
        quicksort(vect, i + 1, fin);

    }

    void swap(int **a, int **b) {

        int *temp;

        temp = *a;
        *a = *b;
        *b = temp;
    }

    bool checkFile(std::ifstream &file) {
        // Check if a given file exists and its size is greater than 0 
        long begin, end;

        // Check whether input file exists or not 
        if (!file)
            return false;

        // Check input file sizes. A valid file should have a size grater than 0 
        begin = file.tellg();
        file.seekg(0, std::ios::end);
        end = file.tellg();
        file.seekg(0, std::ios::beg);
        // Compare difference between file start and end.
        // Depending on result, return True or False 
        if (!(end - begin))
            return false;
        return true;
    }

    char *readLine(std::ifstream &file) {
        StartTiming("char* readLine(std::ifstream &file) ");
        // Read a new line from current input stream. This function is better than
        // standard one since cares of operative system compability. It is useful
        // as well because remove tabs and blank spaces at lines beginning/ending 

        int state;
        std::string nline;
        char *line = nullptr;

        // Check it the end of the file has been reached or not
        if (file.eof())
            return nullptr;

        /* Store first line found. For -Windows & MacOS compatibility- carriage return
        * is considered as well as a new line character */
        //    for( ; (c != '\n') && (c != '\r') && ((!file.eof())); file.read(&c, 1))
        //        nline.resize(nline.size() + 1, c);
        getline(file, nline);

        // Remove blank spaces & tabs from the beginning of the line
        state = nline.find(' ', 0);
        while (state != (int) std::string::npos && state == 0) {
            nline.erase(state, 1);
            state = nline.find(' ', state);
        }

        state = nline.find('\t', 0);
        while (state != (int) std::string::npos && state == 0) {
            nline.erase(state, 1);
            state = nline.find('\t', state);
        }

        /* If there is nothing to return, give back a nullptr pointer ... */
        if (nline.empty())
            return nullptr;

        // Otherwise, initialize the appropiate data structure,
        // dump the data and return it 
        line = new char[nline.size() + 1];
        strcpy(line, nline.c_str());
        return line;
    }

    char *readLine(std::istream &file) {
        // Read a new line from current input stream. This function is better than
        // standard one since cares of operative system compability. It is useful
        // as well because remove tabs and blank spaces at lines beginning/ending

        int state;
        std::string nline;
        char *line = nullptr;

        // Check it the end of the file has been reached or not
        if (file.eof())
            return nullptr;


        /* Store first line found. For -Windows & MacOS compatibility- carriage return
        * is considered as well as a new line character */
        //    for( ; (c != '\n') && (c != '\r') && ((!file.eof())); file.read(&c, 1))
        //        nline.resize(nline.size() + 1, c);
        getline(file, nline);
        /* Remove blank spaces & tabs from the beginning of the line */

        state = nline.find(' ', 0);
        while (state != (int) std::string::npos && state == 0) {
            nline.erase(state, 1);
            state = nline.find(' ', state);
        }

        state = nline.find('\t', 0);
        while (state != (int) std::string::npos && state == 0) {
            nline.erase(state, 1);
            state = nline.find('\t', state);
        }

        // If there is nothing to return, give back a nullptr pointer ...
        if (nline.empty())
            return nullptr;

        // Otherwise, initialize the appropiate data structure,
        // dump the data and return it
        line = new char[nline.size() + 1];
        strcpy(line, &nline[0]);
        return line;
    }

    char *trimLine(std::string nline) {
        // This function is used to remove comments in between a biological sequence.
        // Removes all content surrounded by ("") or ([]).
        // It warns as well when a mismatch for these flags is found

        int pos, next;
        char *line;

        // Set-up lower and upper limit to look for comments inside of input std::string
        pos = -1;

        // Identify comments inside of input sequence and remove it
        while (true) {
            pos = nline.find('\"', (pos + 1));

            // When there is not any more a comment inside of sequence,
            // go out from this loop 
            if (pos == (int) std::string::npos)
                break;

            // Look for closing flag
            next = nline.rfind('\"', nline.size());

            // If a pair of comments flags '"' is found, remove everything inbetween
            if ((int) nline.find('\"', (pos + 1)) == next) {
                nline.erase(pos, (next - pos + 1));
                pos = -1;
            }

            // If there is only one flag '"' for comments inside of sequence,
            // user should be warned about that
            if (pos == next) {
                debug.report(ErrorCode::PossibleMissmatch);
                return nullptr;
            }
        }

        // Look for other kind of comments, in this case those with []
        while (true) {
            pos = -1;
            next = -1;

            // Search for last opened bracket. It is supposed to be the first one for
            // being close
            while ((pos = nline.find('[', (pos + 1))) != (int) std::string::npos)
                next = pos;

            // If no opening bracket has been found.
            // Check if there is any closing one
            if (next == -1) {
                // There are not any bracket in input std::string
                if ((int) nline.find(']', 0) == (int) std::string::npos)
                    break;
                // Otherwise, warn about the error
                debug.report(ErrorCode::BracketsMissmatchFound);
                return nullptr;
            }

            // Look for closest closing bracket to the opening one found
            pos = nline.find(']', (next + 1));

            //If no closing bracket has been found. Warn about the mismatch
            if (pos == (int) std::string::npos) {
                debug.report(ErrorCode::BracketsMissmatchFound);
                return nullptr;
            }

            // When both brackets have been found, remove comments inbetween them
            nline.erase(next, (pos - next + 1));
        }

        // Check if after removing all comments from input std::string there is still part
        // of sequences or not
        if (nline.empty())
            return nullptr;

        // Initialize and store resulting sequence into an appropiate structure
        line = new char[nline.size() + 1];
        strcpy(line, &nline[0]);

        return line;
    }

    std::string getReverse(const std::string &toReverse) {

        std::string line(toReverse.size(), ' ');
        long i, x;

        for (i = toReverse.size() - 1, x = 0; i >= 0; i--, x++)
            line[x] = toReverse[i];
        return line;
    }

    std::string removeCharacter(char c, std::string line) {

        size_t pos;

        pos = line.find(c, 0);
        while (pos != (int) std::string::npos) {
            line.erase(pos, 1);
            pos = line.find(c, pos);
        }

        return line;
    }


    int checkAlignmentType(int seqNumber, int residNumber, std::string *sequences) {
        // All codes from RNA - DNA - AA and DEG for both versions

        static const auto check =
                [](const std::string &pattern,
                   int &counter,
                   const char &character) -> void {
                    size_t pos(0);
                    while ((pos = pattern.find(character, pos)) != std::string::npos) {
                        pos++;
                        counter++;
                    }
                };

        // Inosinic Acid is a nucleotide, but it's not part of DNA or RNA
        static const std::string COMMON = "ACG";
        static const std::string RNA = "U"; // Removed ACG as is in common
        static const std::string DNA = "T"; // Removed ACG as is in common

        static const std::string DEG_NN = "RYSWKMBDHV"; // Indetermination N not added

        static const std::string AA = "PVLIMFYWHKRQEDSTUO"; // Indetermination N not added // Removed ACG as is in common
        static const std::string DEG_AA = "BZ"; // Indetermination X not added

        int rna = 0, dna = 0, deg_nn = 0, aa = 0, deg_aa = 0;

        for (int s = 0; s < seqNumber; s++) {
            for (char c : sequences[s]) {
                if (c == '-' || c == '?' || c == '.') continue;
                c = utils::toUpper(c);
                check(RNA, rna, c);
                check(DNA, dna, c);
                check(DEG_NN, deg_nn, c);
                check(AA, aa, c);
                check(DEG_AA, deg_aa, c);

                size_t pos(0);
                while ((pos = COMMON.find(c, pos)) != std::string::npos) {
                    pos++;
                    aa++;
                    rna++;
                    dna++;
                }
            }
        }

        if (aa > dna && aa > rna) {
            return SequenceTypes::AA | (deg_aa > 0 ? SequenceTypes::DEG : SequenceTypes::NotDefined);
        }
        if (dna > rna) {
            return SequenceTypes::DNA | (deg_nn > 0 ? SequenceTypes::DEG : SequenceTypes::NotDefined);
        }
        return SequenceTypes::RNA | (deg_nn > 0 ? SequenceTypes::DEG : SequenceTypes::NotDefined);


        int i, j, k, l, hitDNA, hitRNA, degenerate, gDNA, gRNA, extDNA, extRNA;
        float ratioDNA, ratioRNA;
        // Standard tables
        static char listRNA[11] = "AGCUNagcun";
        static char listDNA[11] = "AGCTNagctn";

        // Degenerate Nucleotides codes
        static char degeneratedCodes[21] = "MmRrWwSsYyKkVvHhDdBb";

        // For each sequences, this method locks at the 100 letters (excluding gaps).
        // The method is able to distinguish between pure DNA/RNA nucleotides or those
        // containing degenerate Nucleotide letters
        for (i = 0, gDNA = 0, gRNA = 0, extDNA = 0, extRNA = 0; i < seqNumber; i++) {

            // Looks at the 100 letters (excluding gaps) while doesn's get the sequence's end
            // When there are less than a 100 characters, break the loop before reaching that limit

            residNumber = (int) sequences[i].size();
            for (j = 0, k = 0, hitDNA = 0, hitRNA = 0, degenerate = 0; j < residNumber; j++)
                if (sequences[i][j] != '-' && sequences[i][j] != '.' && sequences[i][j] != '?') {
                    k++;

                    /* Recognizes between DNA and RNA. */
                    for (l = 0; l < (int) strlen(listDNA); l++)
                        if (listDNA[l] == sequences[i][j])
                            hitDNA++;

                    for (l = 0; l < (int) strlen(listRNA); l++)
                        if (listRNA[l] == sequences[i][j])
                            hitRNA++;

                    for (l = 0; l < (int) strlen(degeneratedCodes); l++)
                        if (degeneratedCodes[l] == sequences[i][j])
                            degenerate++;
                }

            // If input sequences have less than 95% of nucleotides, even when residues
            // are treated with degenerated codes, consider the input file as containing
            // amino-acidic sequences. 
            ratioDNA = float(degenerate + hitDNA) / k;
            ratioRNA = float(degenerate + hitRNA) / k;

            if (ratioDNA < 0.95 && ratioRNA < 0.95)
                return SequenceTypes::AA;

                // Identify precisely if nucleotides sequences are DNA/RNA strict or
                // any degenerate code has been used in the sequence
            else if (hitRNA > hitDNA && degenerate == 0)
                gRNA++;
            else if (hitRNA > hitDNA && degenerate != 0)
                extRNA++;
            else if (hitRNA < hitDNA && degenerate == 0)
                gDNA++;
            else if (hitRNA < hitDNA && degenerate != 0)
                extDNA++;
        }
        // Return the datatype with greater values, considering always degenerate
        // codes 
        if (extDNA != 0 && extDNA > extRNA)
            return (SequenceTypes::DNA | SequenceTypes::DEG);
        else if (extRNA != 0 && extDNA < extRNA)
            return (SequenceTypes::RNA | SequenceTypes::DEG);
        else if (gRNA > gDNA)
            return SequenceTypes::RNA;
        else if (gDNA >= gRNA)
            return SequenceTypes::DNA;
        else
            return SequenceTypes::NotDefined;
    }

    int *readNumbers(const std::string &line) {

        int i, nElems = 0;
        static int *numbers;

        size_t comma, separ, init;

        comma = size_t(-1);
        while ((comma = line.find(',', comma + 1)) != (int) std::string::npos)
            nElems += 2;
        nElems += 2;

        numbers = new int[nElems + 1];
        numbers[0] = nElems;

        init = 0;
        i = 1;

        do {
            comma = line.find(',', init);
            separ = line.find('-', init);

            if (((separ < comma) || (comma == (int) std::string::npos)) && (separ != (int) std::string::npos)) {
                numbers[i++] = atoi(line.substr(init, separ - init).c_str());
                numbers[i++] = atoi(line.substr(separ + 1, comma - separ - 1).c_str());
                init = comma + 1;
            } else if ((separ > comma) || (separ == (int) std::string::npos)) {
                numbers[i++] = atoi(line.substr(init, comma - init).c_str());
                numbers[i++] = atoi(line.substr(init, comma - init).c_str());
                init = comma + 1;
            }

            if (numbers[i - 2] < 0)
            {
                delete [] numbers;
                return nullptr;
            }
            if (numbers[i - 1] < numbers[i - 2])
            {
                delete [] numbers;
                return nullptr;
            }
            if (comma == (int) std::string::npos)
                break;

        } while (true);

        return numbers;
    }


    char determineColor(char res, const std::string &column) {

        if (toupper(res) == 'G')
            return 'o';

        else if (toupper(res) == 'P')
            return 'y';

        else if (res != '-') {
            switch (toupper(res)) {

                // (W, L, V, I, M, F): {50%, p}{60%, wlvimafcyhp}
                case 87:
                case 76:
                case 86:
                case 73:
                case 77:
                case 70:
                    if (lookForPattern(column, "p", 0.5)) return 'b';
                    else if (lookForPattern(column, "wlvimafcyhp", 0.6)) return 'b';
                    else return 'w';



                    // (A): {50%, p}{60%, wlvimafcyhp}{85% t,s,g}
                case 65:
                    if (lookForPattern(column, "p", 0.5)) return 'b';
                    else if (lookForPattern(column, "wlvimafcyhp", 0.6)) return 'b';
                    else if (lookForPattern(column, "t", 0.85)) return 'b';
                    else if (lookForPattern(column, "s", 0.85)) return 'b';
                    else if (lookForPattern(column, "g", 0.85)) return 'b';
                    else return 'w';



                    // BLUE: (C): {50%, p}{60%, wlvimafcyhp}{85% s}
                    // PINK: (C): {85%, c}
                case 67:
                    if (lookForPattern(column, "p", 0.5)) return 'b';
                    else if (lookForPattern(column, "wlvimafcyhp", 0.6)) return 'b';
                    else if (lookForPattern(column, "s", 0.85)) return 'b';
                    else if (lookForPattern(column, "c", 0.85)) return 'p';
                    else return 'w';



                    // (K, R): {60%, kr}{85%, q}
                case 75:
                case 82:
                    if (lookForPattern(column, "kr", 0.6)) return 'r';
                    else if (lookForPattern(column, "q", 0.85)) return 'r';
                    else return 'w';



                    // (T): {50%, ts}{60%, wlvimafcyhp }
                case 84:
                    if (lookForPattern(column, "ts", 0.5)) return 'g';
                    else if (lookForPattern(column, "wlvimafcyhp", 0.6)) return 'g';
                    else return 'w';



                    // (S): {50%, ts}{80%, wlvimafcyhp }
                case 83:
                    if (lookForPattern(column, "ts", 0.5)) return 'g';
                    else if (lookForPattern(column, "wlvimafcyhp", 0.8)) return 'g';
                    else return 'w';



                    // (N): {50%, n}{85%, d }
                case 78:
                    if (lookForPattern(column, "n", 0.5)) return 'g';
                    else if (lookForPattern(column, "d", 0.85)) return 'g';
                    else return 'w';



                    // (Q): {50%, qe}{60%, kr}
                case 81:
                    if (lookForPattern(column, "qe", 0.5)) return 'g';
                    else if (lookForPattern(column, "kr", 0.6)) return 'g';
                    else return 'w';



                    // (D): {50%, de, n}
                case 68:
                    if (lookForPattern(column, "de", 0.5)) return 'm';
                    else if (lookForPattern(column, "n", 0.5)) return 'm';
                    else return 'w';



                    // (E): {50%, de,qe}
                case 69:
                    if (lookForPattern(column, "de", 0.5)) return 'm';
                    else if (lookForPattern(column, "qe", 0.5)) return 'm';
                    else return 'w';



                    // (H,Y): {50%, p}{60%, wlvimafcyhp}
                case 72:
                case 89:
                    if (lookForPattern(column, "p", 0.5)) return 'c';
                    else if (lookForPattern(column, "wlvimafcyhp", 0.5)) return 'c';
                    else return 'w';

                default:
                    return 'w';
            }
        }
        return 'w';
    }


    bool lookForPattern(const std::string &data, const std::string &pattern, const float threshold) {

        float count = 0;
        int i, j;

        for (i = 0; i < (int) data.size(); i++) {
            for (j = 0; j < (int) pattern.size(); j++) {
                if (utils::toUpper(data[i]) == utils::toUpper(pattern[j])) {
                    count++;
                    break;
                }
            }
        }

        return (count / data.size()) >= threshold;
    }

    void ReplaceStringInPlace(std::string &subject,
                              const std::string &search,
                              const std::string &replace) {
        size_t pos = 0;
        while ((pos = subject.find(search, pos)) != std::string::npos) {
            subject.replace(pos, search.length(), replace);
            pos += replace.length();
        }
    }

    std::string ReplaceString(std::string subject,
                              const std::string &search,
                              const std::string &replace) {
        size_t pos = 0;
        while ((pos = subject.find(search, pos)) != std::string::npos) {
            subject.replace(pos, search.length(), replace);
            pos += replace.length();
        }
        return subject;
    }


    int GetGapStep(int *gapValue, int sequenNumber) {
        // Special cases. Upper and Lower limits.
        if (*gapValue == 0)
            return 11;

        if (*gapValue == sequenNumber)
            return 0;

        float relativeGap = 1.F - float(*gapValue) / sequenNumber;

        if (relativeGap >= .750)
            return 10;
        if (relativeGap >= .500)
            return 9;
        if (relativeGap >= .350)
            return 8;
        if (relativeGap >= .250)
            return 7;
        if (relativeGap >= .200)
            return 6;
        if (relativeGap >= .150)
            return 5;
        if (relativeGap >= .100)
            return 4;
        if (relativeGap >= .050)
            return 3;
        if (relativeGap >= .001)
            return 2;
        return 1;
    }

    int GetGapStep(int *gapValue, float inverseSequenNumber) {
        // Special cases. Upper and Lower limits.
        if (*gapValue == 0)
            return 11;

        float relativeGap = 1.F - float(*gapValue) * inverseSequenNumber;

        if (relativeGap == 0.F)
            return 0;
        if (relativeGap >= .750F)
            return 10;
        if (relativeGap >= .500F)
            return 9;
        if (relativeGap >= .350F)
            return 8;
        if (relativeGap >= .250F)
            return 7;
        if (relativeGap >= .200F)
            return 6;
        if (relativeGap >= .150F)
            return 5;
        if (relativeGap >= .100F)
            return 4;
        if (relativeGap >= .050F)
            return 3;
        if (relativeGap >= .001F)
            return 2;
        return 1;
    }

    int GetSimStep(float *simValue) {

        if (*simValue == 0.F)
            return 11;
        if (*simValue == 1.F)
            return 0;
        if (*simValue >= .750F)
            return 10;
        if (*simValue >= .500F)
            return 9;
        if (*simValue >= .250F)
            return 8;
        if (*simValue >= .100F)
            return 7;
        if (*simValue >= .010F)
            return 6;
        if (*simValue >= .001F)
            return 5;
        if (*simValue >= 1e-4F)
            return 4;
        if (*simValue >= 1e-5F)
            return 3;
        if (*simValue >= 1e-6F)
            return 2;
        return 1;
    }

    int GetConsStep(float *consValue) {
        // Special cases. Upper and Lower limits.
        if (*consValue == 1.F)
            return 11;
        if (*consValue == 0.F)
            return 0;

        if (*consValue >= .750)
            return 10;
        if (*consValue >= .500)
            return 9;
        if (*consValue >= .350)
            return 8;
        if (*consValue >= .250)
            return 7;
        if (*consValue >= .200)
            return 6;
        if (*consValue >= .150)
            return 5;
        if (*consValue >= .100)
            return 4;
        if (*consValue >= .050)
            return 3;
        if (*consValue >= .001)
            return 2;
        return 1;
    }

    bool fileExists(std::string &path) {
        struct stat buffer;
        return (stat(path.c_str(), &buffer) == 0);
    }

    bool fileExists(std::string &&path) {
        struct stat buffer;
        return (stat(path.c_str(), &buffer) == 0);
    }

    char toUpper(char c) {
        if (c >= 'a' and c <= 'z')
            return (char) (c & (~(1 << 5)));
        return c;
    }

    namespace TerminalColors {
        std::map<terminalColor, const std::string> colors = {
                {TerminalColors::RESET,     "\033[0m"},
                {TerminalColors::BLACK,     "\033[30m"},         /* Black      */
                {TerminalColors::RED,       "\033[31m"},         /* Red        */
                {TerminalColors::GREEN,     "\033[32m"},         /* Green      */
                {TerminalColors::YELLOW,    "\033[33m"},         /* Yellow     */
                {TerminalColors::BLUE,      "\033[34m"},         /* Blue       */
                {TerminalColors::MAGENTA,   "\033[35m"},         /* Magenta    */
                {TerminalColors::CYAN,      "\033[36m"},         /* Cyan       */
                {TerminalColors::WHITE,     "\033[37m"},         /* White      */
                {TerminalColors::BOLD,      "\033[1m"},          /* Bold Black */
                {TerminalColors::UNDERLINE, "\033[4m"}           /* Underline  */

        };
    }

}
