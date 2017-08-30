/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v1.4: a tool for automated alignment conversion among different
                 formats.

    2009-2015 Capella-Gutierrez S. and Gabaldon, T.
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

#include "../include/utils.h"
#include "../include/values.h"
#include "../include/defines.h"


void utils::initlVect(int *vector, int tam, int valor) {

  for(int i = 0; i < tam; i++) vector[i] = valor;

}

void utils::initlVect(float *vector, int tam, float valor) {

  for(int i = 0; i < tam; i++) vector[i] = valor;
}

void utils::copyVect(int *vect1, int *vect2, int tam) {

  for(int i = 0; i < tam; i++) vect2[i] = vect1[i];

}

void utils::copyVect(float *vect1, float *vect2, int tam) {

  for(int i = 0; i < tam; i++) vect2[i] = vect1[i];
}

int utils::roundToInf(double number) {

  return ((int) number);
}

int utils::roundInt(double number) {

  return ((int) ((double) number + 0.5));
}

int utils::roundToSup(double number) {

  return ((int) ((double) number + 1.0));
}

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

bool utils::compare(char *a, char *b){

  return(!strcmp(a,b));
}

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

void utils::swap(float *a, float *b){

  float temp;

  temp = *a;
  *a = *b;
  *b = temp;

}

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
  char *line = NULL;

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

char* utils::readLine(std::istream& file)
{
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
  strcpy(line, &nline[0]);
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
  strcpy(line, &nline[0]);

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

int utils::checkAlignmentType(int seqNumber, int residNumber, string *sequences) {

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
    //~ for(j = 0, k = 0, hitDNA = 0, hitRNA = 0, degenerate = 0; j < residNumber && k  < 100; j++)
    for(j = 0, k = 0, hitDNA = 0, hitRNA = 0, degenerate = 0; j < residNumber; j++)
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
      return SequenceTypes::AA;

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
    return SequenceTypes::DNA | SequenceTypes::DEG;
  else if (extRNA != 0 && extDNA < extRNA)
    return SequenceTypes::RNA | SequenceTypes::DEG;
  else if(gRNA > gDNA)
    return SequenceTypes::RNA;
  else
    return SequenceTypes::DNA;
}

int* utils::readNumbers_StartEnd(string line) {

  int comma, nElems = 0;
  static int *numbers;

 comma = -1;
  while((comma = line.find(",", comma + 1)) != (int) string::npos)
    nElems += 2;

  //~ If there is more than two numbers separated by a comma, return NULL
  if(nElems != 2)
    return NULL;

  numbers = new int[2];
  comma = line.find(",", 0);
  numbers[0] = atoi(line.substr(0, comma).c_str());
  numbers[1] = atoi(line.substr(comma+1).c_str());

  return numbers;
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
    if(numbers[i-2] < 0)
      return NULL;
    if(numbers[i-1] < numbers[i-2])
      return NULL;
    if(comma == (int) string::npos)
      break;
   /* ***** ***** ***** ***** ***** ***** ***** ***** */
  } while(true);

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  return numbers;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}


char utils::determineColor(char res, string column) {

    char up = toupper(res);
  /* ***** ***** ***** ***** ***** ***** ***** ***** ***** */
  if(up == 'G')
    return 'o';
  /* ***** ***** ***** ***** ***** ***** ***** ***** ***** */
  else if(up == 'P')
    return 'y';
  /* ***** ***** ***** ***** ***** ***** ***** ***** ***** */

  else if(res != '-') {
    switch(up) {

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (W, L, V, I, M, F): {50%, p}{60%, wlvimafcyhp} */
      case 87: case 76:  case 86: case 73: case 77: case 70:
        if(lookForPattern(column, "P", 0.5))                return 'b';
        else if(lookForPattern(column, "WLVIMAFCYHP", 0.6)) return 'b';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (A): {50%, p}{60%, wlvimafcyhp}{85% t,s,g} */
      case 65:
        if(lookForPattern(column, "P", 0.5))                return 'b';
        else if(lookForPattern(column, "WLVIMAFCYHP", 0.6)) return 'b';
        else if(lookForPattern(column, "T", 0.85))          return 'b';
        else if(lookForPattern(column, "S", 0.85))          return 'b';
        else if(lookForPattern(column, "u", 0.85))          return 'b';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* BLUE: (C): {50%, p}{60%, wlvimafcyhp}{85% s}
       * PINK: (C): {85%, c}
      */
      case 67:
        if(lookForPattern(column, "P", 0.5))                return 'b';
        else if(lookForPattern(column, "WLVIMAFCYHP", 0.6)) return 'b';
        else if(lookForPattern(column, "S", 0.85))          return 'b';
        else if(lookForPattern(column, "C", 0.85))          return 'p';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (K, R): {60%, kr}{85%, q} */
      case 75: case 82:
        if(lookForPattern(column, "KR", 0.6))               return 'r';
        else if(lookForPattern(column, "Q", 0.85))          return 'r';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (T): {50%, ts}{60%, wlvimafcyhp } */
      case 84:
        if(lookForPattern(column, "TS", 0.5))               return 'g';
        else if(lookForPattern(column, "WLVIMAFCYHP", 0.6)) return 'g';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (S): {50%, ts}{80%, wlvimafcyhp } */
      case 83:
        if(lookForPattern(column, "TS", 0.5))               return 'g';
        else if(lookForPattern(column, "WLVIMAFCYHP", 0.8)) return 'g';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (N): {50%, n}{85%, d } */
      case 78:
        if(lookForPattern(column, "N", 0.5))                return 'g';
        else if(lookForPattern(column, "D", 0.85))          return 'g';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (Q): {50%, qe}{60%, kr} */
      case 81:
        if(lookForPattern(column, "QE", 0.5))               return 'g';
        else if(lookForPattern(column, "KR", 0.6))          return 'g';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (D): {50%, de, n} */
      case 68:
        if(lookForPattern(column, "DE", 0.5))               return 'm';
        else if(lookForPattern(column, "N", 0.5))           return 'm';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (E): {50%, de,qe} */
      case 69:
        if(lookForPattern(column, "DE", 0.5))               return 'm';
        else if(lookForPattern(column, "QE", 0.5))          return 'm';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      /* (H,Y): {50%, p}{60%, wlvimafcyhp} */
      case 72: case 89:
        if(lookForPattern(column, "P", 0.5))                return 'c';
        else if(lookForPattern(column, "WLVIMAFCYHP", 0.5)) return 'c';
        else                                                return 'w';
      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
    }
  }
  return 'w';
}


bool utils::lookForPattern(string const& column, string dataset, float level) {

  float count = 0;
  int i, j;

  for(i = 0; i < (int) column.size(); i++) {
    for(j = 0; j < (int) dataset.size(); j++) {
      if(toupper(column[i]) == dataset[j]) {
        count++; break;
      }
    }
  }

  if((count/column.size()) >= level)
    return true;
  else return false;
}

void utils::ReplaceStringInPlace(std::string& subject, const std::string& search, const std::string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
}

std::string utils::ReplaceString(std::string subject, const std::string& search,
                          const std::string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
    return subject;
}

int utils::GetGapStep(int * gapValue, int sequenNumber)
{
    // Special cases. Upper and Lower limits.
    if(*gapValue == 0)
        return 11;
    
    if(*gapValue == sequenNumber)
        return 0;
    
    float relativeGap = 1.F - float(*gapValue) / sequenNumber;
    
    if(relativeGap >= .750)
        return 10;
    if(relativeGap >= .500)
        return 9;
    if(relativeGap >= .350)
        return 8;
    if(relativeGap >= .250)
        return 7;
    if(relativeGap >= .200)
        return 6;
    if(relativeGap >= .150)
        return 5;
    if(relativeGap >= .100)
        return 4;
    if(relativeGap >= .050)
        return 3;
    if(relativeGap >= .001)
        return 2;
    return 1;
}

int utils::GetGapStep(int * gapValue, float inverseSequenNumber)
{
    // Special cases. Upper and Lower limits.
    if(*gapValue == 0)
        return 11;

    float relativeGap = 1.F - float(*gapValue) * inverseSequenNumber;
    
    if(relativeGap == 1.F)
        return 0;
    if(relativeGap >= .750)
        return 10;
    if(relativeGap >= .500)
        return 9;
    if(relativeGap >= .350)
        return 8;
    if(relativeGap >= .250)
        return 7;
    if(relativeGap >= .200)
        return 6;
    if(relativeGap >= .150)
        return 5;
    if(relativeGap >= .100)
        return 4;
    if(relativeGap >= .050)
        return 3;
    if(relativeGap >= .001)
        return 2;
    return 1;
}

int utils::GetSimStep(float * simValue)
{

    if(*simValue == 0.F)
        return 11;
    if(*simValue == 1.F)
        return 0;
    if(*simValue >= .750)
        return 10;
    if(*simValue >= .500)
        return 9;
    if(*simValue >= .250)
        return 8;
    if(*simValue >= .100)
        return 7;
    if(*simValue >= .010)
        return 6;
    if(*simValue >= .001)
        return 5;
    if(*simValue >= 1e-4)
        return 4;
    if(*simValue >= 1e-5)
        return 3;
    if(*simValue >= 1e-6)
        return 2;
    return 1;
}

int utils::GetConsStep(float * consValue)
{
    // Special cases. Upper and Lower limits.
    if(*consValue == 1.F)
        return 11;
    if(*consValue == 0.F)
        return 0;
    
    if(*consValue >= .750)
        return 10;
    if(*consValue >= .500)
        return 9;
    if(*consValue >= .350)
        return 8;
    if(*consValue >= .250)
        return 7;
    if(*consValue >= .200)
        return 6;
    if(*consValue >= .150)
        return 5;
    if(*consValue >= .100)
        return 4;
    if(*consValue >= .050)
        return 3;
    if(*consValue >= .001)
        return 2;
    return 1;
}

void utils::printAccSVG(int * x, float * y, int num, float threshold, std::string title, std::string out)
{
//     threshold = 0.8F;
//     num = 24;
//     x = new int[24] {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23/*,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47*/};
//     y = new float[24] 
//       {0.6354759646 , 12.31123498 , 36.02644687 , 45.07530835 , 48.42959361 , 50.10772454 , 51.37472328 , 
//         52.09420462 , 52.88089026 , 53.57665244 , 54.25067204 , 55.00177894 , 55.61155914 , 56.14425206 , 
//         56.77478653 , 57.37073055 , 57.92022454 , 58.45687065 , 59.24454459 , 59.97292062 , 61.02348197 , 
//         63.64939121 , 100, 100};
/*        
    num = 144;
    x = new int[144] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143};
    y = new float[144] 
      {0.6354759646 , 12.31123498 , 36.02644687 , 45.07530835 , 48.42959361 , 50.10772454 , 51.37472328 , 
        52.09420462 , 52.88089026 , 53.57665244 , 54.25067204 , 55.00177894 , 55.61155914 , 56.14425206 , 
        56.77478653 , 57.37073055 , 57.92022454 , 58.45687065 , 59.24454459 , 59.97292062 , 61.02348197 , 
        63.64939121 , 100, 100, 0.6354759646 , 12.31123498 , 36.02644687 , 45.07530835 , 48.42959361 , 50.10772454 , 51.37472328 , 
        52.09420462 , 52.88089026 , 53.57665244 , 54.25067204 , 55.00177894 , 55.61155914 , 56.14425206 , 
        56.77478653 , 57.37073055 , 57.92022454 , 58.45687065 , 59.24454459 , 59.97292062 , 61.02348197 , 
        63.64939121 , 100, 100,
          0.6354759646 , 12.31123498 , 36.02644687 , 45.07530835 , 48.42959361 , 50.10772454 , 51.37472328 , 
        52.09420462 , 52.88089026 , 53.57665244 , 54.25067204 , 55.00177894 , 55.61155914 , 56.14425206 , 
        56.77478653 , 57.37073055 , 57.92022454 , 58.45687065 , 59.24454459 , 59.97292062 , 61.02348197 , 
        63.64939121 , 100, 100, 0.6354759646 , 12.31123498 , 36.02644687 , 45.07530835 , 48.42959361 , 50.10772454 , 51.37472328 , 
        52.09420462 , 52.88089026 , 53.57665244 , 54.25067204 , 55.00177894 , 55.61155914 , 56.14425206 , 
        56.77478653 , 57.37073055 , 57.92022454 , 58.45687065 , 59.24454459 , 59.97292062 , 61.02348197 , 
        63.64939121 , 100, 100,
        0.6354759646 , 12.31123498 , 36.02644687 , 45.07530835 , 48.42959361 , 50.10772454 , 51.37472328 , 
        52.09420462 , 52.88089026 , 53.57665244 , 54.25067204 , 55.00177894 , 55.61155914 , 56.14425206 , 
        56.77478653 , 57.37073055 , 57.92022454 , 58.45687065 , 59.24454459 , 59.97292062 , 61.02348197 , 
        63.64939121 , 100, 100, 0.6354759646 , 12.31123498 , 36.02644687 , 45.07530835 , 48.42959361 , 50.10772454 , 51.37472328 , 
        52.09420462 , 52.88089026 , 53.57665244 , 54.25067204 , 55.00177894 , 55.61155914 , 56.14425206 , 
        56.77478653 , 57.37073055 , 57.92022454 , 58.45687065 , 59.24454459 , 59.97292062 , 61.02348197 , 
        63.64939121 , 100, 100
    };*/
//     cout << "CUT POINT " << threshold << endl;
    float xMultiplier = 300.F / num;
    int spacer = num / 10;
    cout << "xMultiplier = " << xMultiplier << endl << "Spacer" << spacer << endl;;
    
    ofstream file;
    file.open(out);  

    // svg header  
    file    << "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" version=\"1.1\" height=\"175\" width=\"350px\">" << endl;
    // Gray box
    file    << "<rect x=\"0\" width=\"350\" y=\"0\" height=\"175\" style=\"fill:gray\"/>" << endl;
    // White box
    file    << "<rect x=\"25\" width=\"300\" y=\"49\" height=\"102\" style=\"fill:white\"/>" << endl;
    // Header text
    file    << "<text text-anchor=\"middle\" x=\"175\" y=\"31.25\">" << title << "</text>" << endl;
    // Y Axis
    for (int i = 0; i < 110; i+= 10)
    {
        // Y Axis lines
        file    << "<line x1=\"25\" y1=\"" << 50 + i << "\" x2=\"325\" y2=\"" << 50 + i << "\" style=\"stroke:black;stroke-width:0.5\" stroke-dasharray=\"1, 1\" opacity=\"0.25\"/>" << endl;
        // Y Axis labels
        file    << "<text x=\"22\" y= \"" << 50 + i + 1.25F << "\" text-anchor=\"end\" font-size=\"5\">" << 100 - i << "%" << "</text>" << endl;
    }
    // X Axis
    for (int i = 0; i < num; i++)
    {
        // X axis lines
        file    << "<line " <<
                    "x1=\""  << 25 + i * xMultiplier <<"\" y1=\"50\" " <<
                    "x2=\""  << 25 + i * xMultiplier <<"\" y2=\"150\" " <<
                    "style=\"stroke:black;stroke-width:0.5\" " <<
                    "stroke-dasharray=\"1, 1\" " <<
                    "opacity=\"0.125\"/>" << endl;
        
        // X axis labels
        if (i % spacer == 0)
        {
            // Indicator rectangle
            file    << "<rect " <<
                    "x=\""  << 25 + i * xMultiplier <<"\" " <<
                    "y=\"49\" " <<
                    "width=\""<< std::max(1.F, xMultiplier) <<"\" " <<
                    "height=\"102\" " <<
                    "style=\"fill:gray\" " <<
                    "opacity=\"0.125\"/>" << endl;
                    
            // Labels
            file    << "<text x=\""  << 25 + (i + 0.5F) * xMultiplier <<"\" " 
                    << "y=\"" << 160 
                    << "\" text-anchor=\"middle\" font-size=\"5\">"
                    << x[i]
                    << "</text>" << endl;
        }
    }
    
    // Points Line
    file    << "<polyline stroke-linecap=\"round\" style=\"fill:none;stroke:black;stroke-width:1\" opacity=\"0.8\" points=\"";
    for (int i = 0; i < num; i++)
    {
        file    << (25 + (x[i] + 0.5F) * xMultiplier) << "," << 50 + 100 - y[i] << " ";
    }
    file    << "\"/>" << endl;
    
    // Kept
    file    << "<rect " <<
                "x=\""  << 25 <<"\" " <<
                "y=\"49\" " <<
                "width=\""<< 300 * threshold <<"\" " <<
                "height=\""<< 102  <<"\" " <<
                "style=\"fill:green\" " <<
                "opacity=\"0.125\"/>" << endl;
                
    // Rejected
    file    << "<rect " <<
                "x=\""  << 25 + 300 * threshold <<"\" " <<
                "y=\"49\" " <<
                "width=\""<< 300 * (1.F - threshold) <<"\" " <<
                "height=\""<< 102  <<"\" " <<
                "style=\"fill:red\" " <<
                "opacity=\"0.125\"/>" << endl;
    
    file << "</svg>";
    file.close();
    
//     delete [] x;
//     delete [] y;
}

void utils::printColSVG(float * gapScore, int num, float threshold, std::string title, std::string out)
{

    float xMultiplier = 300.F / num;
    int spacer = num / 10;
//     cout << "xMultiplier = " << xMultiplier << endl << "Spacer" << spacer << endl;;
    cout << "N = " << num << endl;
    
    ofstream file;
    file.open(out);  

    // svg header  
    file    << "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" version=\"1.1\" height=\"175\" width=\"350px\">" << endl;
    // Gray box
    file    << "<rect x=\"0\" width=\"350\" y=\"0\" height=\"175\" style=\"fill:gray\"/>" << endl;
    // White box
    file    << "<rect x=\"25\" width=\"300\" y=\"49\" height=\"102\" style=\"fill:white\"/>" << endl;
    // Header text
    file    << "<text text-anchor=\"middle\" x=\"175\" y=\"31.25\">" << title << "</text>" << endl;
    // Y Axis
    for (int i = 0; i < 11; i+= 1)
    {
        // Y Axis lines
        file    << "<line x1=\"25\" y1=\"" << 50 + i * 10 << "\" x2=\"325\" y2=\"" << 50 + i * 10 << "\" style=\"stroke:black;stroke-width:0.5\" stroke-dasharray=\"1, 1\" opacity=\"0.25\"/>" << endl;
        // Y Axis labels
        file    << "<text x=\"22\" y= \"" << 50 + i * 10 + 1.25F << "\" text-anchor=\"end\" font-size=\"5\">" << 1 - i * 0.1F << "</text>" << endl;
    }
    // X Axis
    for (int i = 0; i < num; i++)
    {
        // X axis lines
        file    << "<line " <<
                    "x1=\""  << 25 + i * xMultiplier <<"\" y1=\"50\" " <<
                    "x2=\""  << 25 + i * xMultiplier <<"\" y2=\"150\" " <<
                    "style=\"stroke:black;stroke-width:0.5\" " <<
                    "stroke-dasharray=\"1, 1\" " <<
                    "opacity=\"0.125\"/>" << endl;
        
        // X axis labels
        if (i % spacer == 0)
        {
            // Indicator rectangle
            file    << "<rect " <<
                    "x=\""  << 25 + i * xMultiplier <<"\" " <<
                    "y=\"49\" " <<
                    "width=\""<< std::max(1.F, xMultiplier) <<"\" " <<
                    "height=\"102\" " <<
                    "style=\"fill:gray\" " <<
                    "opacity=\"0.125\"/>" << endl;
                    
            // Labels
            file    << "<text x=\""  << 25 + (i + 0.5F) * xMultiplier <<"\" " 
                    << "y=\"" << 160 
                    << "\" text-anchor=\"middle\" font-size=\"5\">"
                    << i
                    << "</text>" << endl;
        }
    }
    
    // Points Line
    file    << "<polyline stroke-linecap=\"round\" style=\"fill:none;stroke:black;stroke-width:1\" opacity=\"0.8\" points=\"";
    for (int i = 0; i < num; i++)
    {
        file    << 25 + (i + 0.5F) * xMultiplier << "," << 50 + 100 - gapScore[i] << " ";
    }
    file    << "\"/>" << endl;
    
    // Kept
    file    << "<rect " <<
                "x=\"25\" " <<
                "y=\""<< 49 <<"\" " <<
                "width=\"300\" " <<
                "height=\""<< 1 + 100 * threshold  <<"\" " <<
                "style=\"fill:green\" " <<
                "opacity=\"0.125\"/>" << endl;
                
    // Rejected
    file    << "<rect " <<
                "x=\"25\" " <<
                "y=\""<< 50 + 100 * threshold <<"\" " <<
                "width=\"300\" " <<
                "height=\""<< 101 - 100 * threshold <<"\" " <<
                "style=\"fill:red\" " <<
                "opacity=\"0.125\"/>" << endl;
    
    file << "</svg>";
    file.close();
//     
//     delete [] x;
//     delete [] y;
}
