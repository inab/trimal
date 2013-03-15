/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v1.4: a tool for automated alignment conversion among different
                 formats.

    2009-2011 Capella-Gutierrez S. and Gabaldon, T.
              [scapella, tgabaldon]@crg.es

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

***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#include "utils.h"
#include "values.h"
#include "defines.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++
| void utils::initVect(int *, int, int)          |
|      This method is used to initializate all   |
|      positions of a vector with a given value. |
++++++++++++++++++++++++++++++++++++++++++++++++*/

void utils::initlVect(int *vector, int tam, int valor) {

  for(int i = 0; i < tam; i++) vector[i] = valor;

}

void utils::initlVect(float *vector, int tam, float valor) {

  for(int i = 0; i < tam; i++) vector[i] = valor;

}


/*+++++++++++++++++++++++++++++++++++++++++++++
| void utils::copyVect(int *, int *, int)     |
|      This method copies integer vector 1 to |
|      integer vector 2.                      |
+++++++++++++++++++++++++++++++++++++++++++++*/

void utils::copyVect(int *vect1, int *vect2, int tam) {

  for(int i = 0; i < tam; i++) vect2[i] = vect1[i];

}


/*+++++++++++++++++++++++++++++++++++++++++++++++
| void utils::copyVect(float *, float *, float) |
|      This method copies float vector 1 to     |
|      float vector 2.                          |
+++++++++++++++++++++++++++++++++++++++++++++++*/

void utils::copyVect(float *vect1, float *vect2, int tam) {

  for(int i = 0; i < tam; i++) vect2[i] = vect1[i];

}


/*+++++++++++++++++++++++++++++++++++++++++
| int utils::roundToInf(double)           |
|      This method rounds a double number |
|      to the inferior integer.           |
+++++++++++++++++++++++++++++++++++++++++*/

int utils::roundToInf(double number) {

  return ((int) number);
}


/*+++++++++++++++++++++++++++++++++++++++++
| int utils::roundInt(double)             |
|      This method rounds a double number |
|      to a integer.                      |
+++++++++++++++++++++++++++++++++++++++++*/

int utils::roundInt(double number) {

  return ((int) ((double) number + 0.5));
}


/*+++++++++++++++++++++++++++++++++++++++++
| int utils::roundToSup(double)           |
|      This method rounds a double number |
|      to the greater integer.            |
+++++++++++++++++++++++++++++++++++++++++*/

int utils::roundToSup(double number) {

  return ((int) ((double) number + 1.0));
}


/*+++++++++++++++++++++++++++++++++++++++++
| int utils::max(int, int)                |
|      This method returns the maximum    |
|      value of the two given arguments.  |
+++++++++++++++++++++++++++++++++++++++++*/

int utils::max(int x, int y) {

  if(x > y) return x;
  else      return y;
}

float utils::max(float x, float y) {

  if(x > y) return x;
  else      return y;
}

double utils::max(double x, double y) {

  if(x > y) return x;
  else      return y;
}

int utils::min(int x, int y) {

  if(x < y) return x;
  else      return y;
}

float utils::min(float x, float y) {

  if(x < y) return x;
  else      return y;
}

double utils::min(double x, double y) {

  if(x < y) return x;
  else      return y;
}

/*+++++++++++++++++++++++++++++++++++++++++++
| bool utils::isNumber(char *)              |
|      This method checks if the given      |
|      string is a float number.            |
+++++++++++++++++++++++++++++++++++++++++++*/

bool utils::isNumber(char *num){

  int tam = strlen(num);
  int i, flt = 1, expn = 1, sgn = 1;

  for(i = 0; i < tam; i++) {
    if(num[i] == '.' && flt)
      flt = 0;

    else if(((num[i] == 'e') ||(num[i] == 'E')) && expn)
      expn = 0;

    else if(((num[i] == '+') ||(num[i] == '-')) && sgn) {
      if(!expn) sgn = 0;
    }
    else if(num[i] > '9' || num[i] < '0')
      return false;
  }

  return true;

}


/*+++++++++++++++++++++++++++++++++++++++++++
| bool utils::compare(char *, char *)       |
|      This method compares the two strings |
|      given, and returns true if the two   |
|      strings are equal.                   |
+++++++++++++++++++++++++++++++++++++++++++*/

bool utils::compare(char *a, char *b){

  return(!strcmp(a,b));
}


/*++++++++++++++++++++++++++++++++++++++++++
| void utils::removeSpaces(char *, char *) |
|      This method removes spaces in the   |
|      input string and put the result in  |
|      the output string.                  |
++++++++++++++++++++++++++++++++++++++++++*/

void utils::removeSpaces(char *in, char *out){

  unsigned int i, j = 0;

  for(i = 0; i < strlen(in); i++){

    if(in[i] != ' ' && in[i] != '\t'){
      out[j] = in[i];
      j++;
    }
  }
  out[j] = '\0';
}


/*++++++++++++++++++++++++++++++++++++++++++
| void utils::quicksort(float *, int, int) |
|      This method sorts the vector using  |
|      the quicksort method.               |
++++++++++++++++++++++++++++++++++++++++++*/

void utils::quicksort(float *vect, int ini, int fin) {

  float elem_div;
  int i, j;

  if ((ini >= fin) || (fin < 0))
    return;

  elem_div = vect[fin];
  i = ini - 1;
  j = fin;

  while (1) {

    while (vect[++i] < elem_div)
      if(i == fin)
        break;

    while (vect[--j] > elem_div)
      if(j == 0)
        break;

    if(i < j)
      swap(&vect[i], &vect[j]);
    else
      break;
  }

  swap(&vect[i], &vect[fin]);

  quicksort(vect, ini, i - 1);
  quicksort(vect, i + 1, fin);
}


/*++++++++++++++++++++++++++++++++++++++++
| void utils::swap(float *, float *)     |
|      This method swaps the values in a |
|      and b.                            |
++++++++++++++++++++++++++++++++++++++++*/

void utils::swap(float *a, float *b){

  float temp;

  temp = *a;
  *a = *b;
  *b = temp;

}

/*++++++++++++++++++++++++++++++++++++++++++
| void utils::quicksort(float *, int, int) |
|      This method sorts the vector using  |
|      the quicksort method.               |
++++++++++++++++++++++++++++++++++++++++++*/

void utils::quicksort(int *vect, int ini, int fin) {

  int i, j, elem_div;

  if ((ini >= fin) || (fin < 0))
    return;

  elem_div = vect[fin];
  i = ini - 1;
  j = fin;

  while (1) {

    while (vect[++i] < elem_div)
      if(i == fin)
        break;

    while (vect[--j] > elem_div)
      if(j == 0)
        break;

    if(i < j)
      swap(&vect[i], &vect[j]);
    else
      break;
  }

  swap(&vect[i], &vect[fin]);

  quicksort(vect, ini, i - 1);
  quicksort(vect, i + 1, fin);
}


/*++++++++++++++++++++++++++++++++++++++++
| void utils::swap(float *, float *)     |
|      This method swaps the values in a |
|      and b.                            |
++++++++++++++++++++++++++++++++++++++++*/

void utils::swap(int *a, int *b) {

  int temp;

  temp = *a;
  *a = *b;
  *b = temp;

}


void utils::quicksort(int **vect, int ini, int fin) {

  float elem_div;
  int i, j;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if ((ini >= fin) || (fin < 0))
    return;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  elem_div = vect[fin][0];
  i = ini - 1;
  j = fin;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  while (true) {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    while (vect[++i][0] < elem_div) if(i == fin) break;
    while (vect[--j][0] > elem_div) if(j == 0)   break;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if(i < j) swap(&vect[i], &vect[j]);
    else break;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  swap(&vect[i], &vect[fin]);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  quicksort(vect, ini, i - 1);
  quicksort(vect, i + 1, fin);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

void utils::swap(int **a, int **b){

  int *temp;

  temp = *a;
  *a = *b;
  *b = temp;
}

/* *****************************************************************************
 *
 * START - Refactored code
 *
 * ************************************************************************** */

bool utils::checkFile(ifstream &file) {
  /* Check if a given file exists and its size is greater than 0 */
  long begin, end;

  /* Check whether input file exists or not */
  if(!file)
    return false;

  /* Check input file sizes. A valid file should have a size grater than 0 */
  begin = file.tellg();
  file.seekg(0, ios::end);
  end = file.tellg();
  file.seekg(0, ios::beg);
  /* Compare difference between file start and end.
   * Depending on result, return True or False */
  if(!(end - begin))
    return false;
  return true;
}

char* utils::readLine(ifstream &file) {
  /* Read a new line from current input stream. This function is better than
   * standard one since cares of operative system compability. It is useful
   * as well because remove tabs and blank spaces at lines beginning/ending */

  int state;
  char c = ' ';
  string nline;
  static char *line = NULL;

  /* Check it the end of the file has been reached or not */
  if(file.eof())
    return NULL;

  /* Store first line found. For -Windows & MacOS compatibility- carriage return
   * is considered as well as a new line character */
  for( ; (c != '\n') && (c != '\r') && ((!file.eof())); file.read(&c, 1))
    nline.resize(nline.size() + 1, c);

  /* Remove blank spaces & tabs from the beginning of the line */
  state = nline.find(" ", 0);
  while(state != (int) string::npos && state == 0) {
    nline.erase(state, 1);
    state = nline.find(" ", state);
  }

  state = nline.find("\t", 0);
  while(state != (int) string::npos && state == 0) {
    nline.erase(state, 1);
    state = nline.find("\t", state);
  }

  /* If there is nothing to return, give back a NULL pointer ... */
  if(nline.size() == 0)
    return NULL;

  /* Otherwise, initialize the appropiate data structure,
   * dump the data and return it */
  line = new char[nline.size() + 1];
  strcpy(line, nline.c_str());
  return line;
}

char* utils::trimLine(string nline) {
  /* This function is used to remove comments inbetween a biological sequence.
   * Remove all content surrounded by ("") or ([]). It wans as well when a
   * mismatch for these flags is found */

  int pos, next;
  static char *line;

  /* Set-up lower and upper limit to look for comments inside of input string */
  pos = -1;

  /* Identify comments inside of input sequence and remove it */
  while(true) {
    pos  = nline.find("\"", (pos + 1));

    /* When there is not any more a comment inside of sequence,
     * go out from this loop */
    if(pos == (int) string::npos)
      break;

    /* Look for closing flag */
    next = nline.rfind("\"", nline.size());

    /* If a pair of comments flags '"' is found, remove everything inbetween */
    if((int) nline.find("\"", (pos + 1)) == next) {
      nline.erase(pos, (next - pos + 1));
      pos = -1;
    }

    /* If there is only one flag '"' for comments inside of sequence,
     * user should be warned about that */
    if (pos == next) {
      cerr << endl << "ERROR: Possible (\") mismatch for comments" << endl;
      return NULL;
    }
  }

  /* Look for other kind of comments, in this case those with [] */
  while(true) {
    pos = -1;
    next = -1;

    /* Search for last opened bracket. It is supposed to be the first one for
     * being close */
    while((pos = nline.find("[", (pos + 1))) != (int) string::npos)
      next = pos;

    /* If no opening bracket has been found.
     * Check if there is any closing one */
    if (next == -1) {
      /* There are not any bracket in input string */
      if ((int) nline.find("]", 0) == (int) string::npos)
        break;
      /* Otherwise, warn about the error */
      cerr << endl << "ERROR: Brackets (]) mismatch found" << endl;
      return NULL;
    }

    /* Look for closest closing bracket to the opening one found */
    pos = nline.find("]", (next + 1));

    /* If no closing bracket has been found. Warn about the mismatch */
    if (pos == (int) string::npos) {
      cerr << endl << "ERROR: Brackets ([) mismatch found" << endl;
      return NULL;
    }

    /* When both brackets have been found, remove comments inbetween them */
    nline.erase(next, (pos - next + 1));
  }

  /* Check if after removing all comments from input string there is still part
   * of sequences or not */
  if(nline.size() == 0)
    return NULL;

  /* Initialize and store resulting sequence into an appropiate structure */
  line = new char[nline.size() + 1];
  strcpy(line, nline.c_str());

  return line;
}
/* *****************************************************************************
 *
 * END - Refactored code
 *
 * ************************************************************************** */
string utils::getReverse(string toReverse) {

  string line;
  int i;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(i = toReverse.size() - 1; i >= 0; i--)
    line += toReverse[i];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  return line;
}

string utils::removeCharacter(char c, string line) {

  int pos;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  pos = line.find(c, 0);
  while(pos != (int) string::npos) {
    line.erase(pos, 1);
    pos = line.find(c, pos);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  return line;
}

int utils::checkTypeAlignment(int seqNumber, int residNumber, string *sequences) {

  int i, j, k, l, hitDNA, hitRNA, degenerate, gDNA, gRNA, extDNA, extRNA;
  float ratioDNA, ratioRNA;
  /* Standard tables */
  char listRNA[11] = "AGCUNagcun";
  char listDNA[11] = "AGCTNagctn";

  /* Degenerate Nucleotides codes */
  char degeneratedCodes[21] = "MmRrWwSsYyKkVvHhDdBb";

  /* For each sequences, this method locks at the 100 letters (excluding gaps).
   * The method is able to distinguish between pure DNA/RNA nucleotides or those
   * containing degenerate Nucleotide letters */
  for(i = 0, gDNA = 0, gRNA = 0, extDNA = 0, extRNA = 0; i < seqNumber; i++) {

    /* Looks at the 100 letters (excluding gaps) while doesn's get the sequence's end */
    /* When there are less than a 100 characters, break the loop before reaching that limit */
    residNumber = (int) sequences[i].size();
    for(j = 0, k = 0, hitDNA = 0, hitRNA = 0, degenerate = 0; j < residNumber && k  < 100; j++)
      if(sequences[i][j] != '-' && sequences[i][j] != '.' && sequences[i][j] != '?') {
        k++;

        /* Recognizes between DNA and RNA. */
        for(l = 0; l < (int) strlen(listDNA); l++)
          if(listDNA[l] == sequences[i][j])
            hitDNA++;

        for(l = 0; l < (int) strlen(listRNA); l++)
          if(listRNA[l] == sequences[i][j])
            hitRNA++;

        for(l = 0; l < (int) strlen(degeneratedCodes); l++)
          if(degeneratedCodes[l] == sequences[i][j])
            degenerate++;
      }

    /* If input sequences have less than 95% of nucleotides, even when residues
     * are treated with degenerated codes, consider the input file as containing
     * amino-acidic sequences. */
    ratioDNA = float(degenerate + hitDNA)/k;
    ratioRNA = float(degenerate + hitRNA)/k;

    if(ratioDNA < 0.95 && ratioRNA < 0.95)
      return AAType;

    /* Identify precisely if nucleotides sequences are DNA/RNA strict or
     * any degenerate code has been used in the sequence */
    else if(hitRNA > hitDNA && degenerate == 0)
      gRNA++;
    else if(hitRNA > hitDNA && degenerate != 0)
      extRNA++;
    else if(hitRNA < hitDNA && degenerate == 0)
      gDNA++;
    else if(hitRNA < hitDNA && degenerate != 0)
      extDNA++;
  }
  /* Return the datatype with greater values, considering always degenerate
   * codes */
  if (extDNA != 0 && extDNA > extRNA)
    return DNADeg;
  else if (extRNA != 0 && extDNA < extRNA)
    return RNADeg;
  else if(gRNA > gDNA)
    return RNAType;
  else
    return DNAType;
}

int* utils::readNumbers(string line) {

  int i, comma, separ, init, nElems = 0;
  static int *numbers;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
 comma = -1;
  while((comma = line.find(",", comma + 1)) != (int) string::npos)
    nElems += 2;
  nElems += 2;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  numbers = new int[nElems + 1];
  numbers[0] = nElems;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  init = 0;
  i = 1;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  do {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    comma = line.find(",", init);
    separ = line.find("-", init);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if(((separ < comma) || (comma == (int) string::npos)) && (separ != (int) string::npos)) {
      numbers[i++] = atoi(line.substr(init, separ - init).c_str());
      numbers[i++] = atoi(line.substr(separ+1, comma - separ - 1).c_str());
      init = comma + 1;
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    else if((separ > comma) || (separ == (int) string::npos)) {
      numbers[i++] = atoi(line.substr(init, comma - init).c_str());
      numbers[i++] = atoi(line.substr(init, comma - init).c_str());
      init = comma + 1;
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if(numbers[i-2] < 0) return NULL;
    if(numbers[i-1] < numbers[i-2]) return NULL;
    if(comma == (int) string::npos) break;
   /* ***** ***** ***** ***** ***** ***** ***** ***** */
  } while(true);

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  return numbers;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}


char utils::determineColor(char res, string column) {

  /* ***** ***** ***** ***** ***** ***** ***** ***** ***** */
  if(toupper(res) == 'G')
    return 'o';
  /* ***** ***** ***** ***** ***** ***** ***** ***** ***** */
  else if(toupper(res) == 'P')
    return 'y';
  /* ***** ***** ***** ***** ***** ***** ***** ***** ***** */

  else if(res != '-') {
    switch(toupper(res)) {

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (W, L, V, I, M, F): {50%, p}{60%, wlvimafcyhp} */
      case 87: case 76:  case 86: case 73: case 77: case 70:
        if(lookForPattern(column, "p", 0.5))                return 'b';
        else if(lookForPattern(column, "wlvimafcyhp", 0.6)) return 'b';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (A): {50%, p}{60%, wlvimafcyhp}{85% t,s,g} */
      case 65:
        if(lookForPattern(column, "p", 0.5))                return 'b';
        else if(lookForPattern(column, "wlvimafcyhp", 0.6)) return 'b';
        else if(lookForPattern(column, "t", 0.85))          return 'b';
        else if(lookForPattern(column, "s", 0.85))          return 'b';
        else if(lookForPattern(column, "g", 0.85))          return 'b';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* BLUE: (C): {50%, p}{60%, wlvimafcyhp}{85% s}
       * PINK: (C): {85%, c}
      */
      case 67:
        if(lookForPattern(column, "p", 0.5))                return 'b';
        else if(lookForPattern(column, "wlvimafcyhp", 0.6)) return 'b';
        else if(lookForPattern(column, "s", 0.85))          return 'b';
        else if(lookForPattern(column, "c", 0.85))          return 'p';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (K, R): {60%, kr}{85%, q} */
      case 75: case 82:
        if(lookForPattern(column, "kr", 0.6))               return 'r';
        else if(lookForPattern(column, "q", 0.85))          return 'r';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (T): {50%, ts}{60%, wlvimafcyhp } */
      case 84:
        if(lookForPattern(column, "ts", 0.5))               return 'g';
        else if(lookForPattern(column, "wlvimafcyhp", 0.6)) return 'g';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (S): {50%, ts}{80%, wlvimafcyhp } */
      case 83:
        if(lookForPattern(column, "ts", 0.5))               return 'g';
        else if(lookForPattern(column, "wlvimafcyhp", 0.8)) return 'g';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (N): {50%, n}{85%, d } */
      case 78:
        if(lookForPattern(column, "n", 0.5))                return 'g';
        else if(lookForPattern(column, "d", 0.85))          return 'g';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (Q): {50%, qe}{60%, kr} */
      case 81:
        if(lookForPattern(column, "qe", 0.5))               return 'g';
        else if(lookForPattern(column, "kr", 0.6))          return 'g';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (D): {50%, de, n} */
      case 68:
        if(lookForPattern(column, "de", 0.5))               return 'm';
        else if(lookForPattern(column, "n", 0.5))           return 'm';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (E): {50%, de,qe} */
      case 69:
        if(lookForPattern(column, "de", 0.5))               return 'm';
        else if(lookForPattern(column, "qe", 0.5))          return 'm';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (H,Y): {50%, p}{60%, wlvimafcyhp} */
      case 72: case 89:
        if(lookForPattern(column, "p", 0.5))                return 'c';
        else if(lookForPattern(column, "wlvimafcyhp", 0.5)) return 'c';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
    }
  }
  return 'w';
}


bool utils::lookForPattern(string column, string dataset, float level) {

  float count = 0;
  int i, j;

  for(i = 0; i < (int) column.size(); i++) {
    for(j = 0; j < (int) dataset.size(); j++) {
      if(toupper(column[i]) == toupper(dataset[j])) {
        count++; break;
      }
    }
  }

  if((count/column.size()) >= level)
    return true;
  else return false;
}

