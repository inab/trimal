/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 

    trimAl v1.3: a tool for automated alignment trimming in large-scale 
                 phylogenetics analyses.

    readAl v1.2: a tool for automated alignment conversion among different
                 formats.

    Copyright (C) 2009 Capella-Gutierrez S. and Gabaldon, T.
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

 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#include "utils.h"
#include "values.h"

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

/* ****************************************************************************************************************** */
/* ****************************************************************************************************************** */


bool utils::checkFile(ifstream &file) {

  long begin, end;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check file */
  if(!file) return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** Check the alignment content ***** */
  begin = file.tellg();
  file.seekg(0, ios::end);
  end = file.tellg();
	file.seekg(0, ios::beg);
  if(!(end - begin)) return false;
  return true;
  /* ***** ***** ***** ***** ***** ***** ***** */
}

char* utils::readLine(ifstream &file) {	

  char c;
  int state;
  string nline;
  static char *line;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
    nline.resize(nline.size() + 1, c);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 

  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  if(file.eof())
    return NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  state = nline.find("\r", 0);
  while(state != (int) string::npos) {
    nline.erase(state, 1);
    state = nline.find("\r", state);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 

  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  line = new char[nline.size() + 1];  
  strcpy(line, nline.c_str());
  return line;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
}


char* utils::readLine2(ifstream &file) {	

  char c;
  int state;
  string nline;
  static char *line;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  for(nline.clear(), file.read(&c, 1); (c != '\n') && ((!file.eof())); file.read(&c, 1))
    nline.resize(nline.size() + 1, c);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 

  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  cerr << endl << nline << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 

  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  if(file.eof())
    return NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  state = nline.find("\r", 0);
  while(state != (int) string::npos) {
    nline.erase(state, 1);
    state = nline.find("\r", state);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 

  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  line = new char[nline.size() + 1];  
  strcpy(line, nline.c_str());
  return line;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
}


char* utils::trimLine(string nline) {

  int state, pos, npos, next;
  static char *line;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  next = nline.size() + 1; 
  pos = -1;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  /* We need to look for comments along the sequence */
  while(true) {
    pos  = nline.find("\"", (pos + 1));
    npos = nline.find("\"", (pos + 1));
    next = nline.rfind("\"", next - 1);

    if(pos != (int) string::npos) {
      if(npos == next) {
        nline.erase(pos, (next - pos + 1));
        next = nline.size() + 1;
        pos = -1; 
      }
    } else break;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  do {
    pos = -1; 
    npos = pos;

    while((pos = nline.find("[", (pos + 1))) != (int) string::npos)
      npos = pos;

    if((int) nline.find("]", (npos + 1)) != (int) string::npos)
      nline.erase(npos, (nline.find("]", (npos + 1)) - npos + 1));

  } while(npos != -1);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 

  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  state = nline.find("\r", 0);
  while(state != (int) string::npos) {
    nline.erase(state, 1);
    state = nline.find("\r", state);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 

  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  line = new char[nline.size() + 1];  
  strcpy(line, nline.c_str());
  return line;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
}

char* utils::readLineMEGA(ifstream &file) {	
  char c;
  int state;
  string nline;
  static char *line;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

	nline += c;
  while((c != '\n') && (!file.eof())) {
		/* ***** ***** ***** ***** ***** ***** ***** ***** */
		for(file.read(&c, 1); (c != '\n') && (!file.eof()) && (c != '#'); file.read(&c, 1))
			nline += c;
		/* ***** ***** ***** ***** ***** ***** ***** ***** */

		/* ***** ***** ***** ***** ***** ***** ***** ***** */
   if(file.eof() || (c == '#')) break;
   nline += c;
   file.read(&c, 1);
		/* ***** ***** ***** ***** ***** ***** ***** ***** */  
  }

	/* ***** ***** ***** ***** ***** ***** ***** ***** */  
  state = nline.find("\r", 0);
  while(state != (int) string::npos) {
    nline.erase(state, 1);
    state = nline.find("\r", state);
  }
	/* ***** ***** ***** ***** ***** ***** ***** ***** */ 

	/* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  line = new char[nline.size() + 1];  
  strcpy(line, nline.c_str());
  return line;
	/* ***** ***** ***** ***** ***** ***** ***** ***** */ 
}

string utils::getReverse(string toReverse) {

  string line;
  int i;
  
  /* ***** ***** ***** ***** ***** ***** ***** ***** */ 
  for(i = toReverse.size(); i >= 0; i--)
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

  int i, j, k, l, hitDNA, hitRNA, gDNA, gRNA;
  char listRNA[6] = "AGCUN";
  char listDNA[6] = "AGCTN";

  /* For each sequences, this method locks at the 100 letters (excluding gaps). If 95% or more of those letters are 
     valid nucleotides, then the files is treated as nucleotides. The method also recognizes between ADN and ARN. */
  for(i = 0, gDNA = 0, gRNA = 0; i < seqNumber; i++) {

    /* Looks at the 100 letters (excluding gaps) while doesn's get the sequence's end */
    for(j = 0, k = 0, hitDNA = 0, hitRNA = 0; ((j < residNumber) && (k  < 100)); j++)
      if(sequences[i][j] != '-') {
        k++;

        /* Recognizes between DNA and RNA. */
        for(l = 0; l < (int) strlen(listDNA); l++)
          if(listDNA[l] == sequences[i][j])
            hitDNA++;

        for(l = 0; l < (int) strlen(listRNA); l++)
          if(listRNA[l] == sequences[i][j])
            hitRNA++;
      }

    /* If for an alignment's sequences the nucleotides don't achieve the threshold, then the method finish and fix
       the alignment's datatype as AminoAcids. */
    if((((float) hitDNA/k) < 0.95) && (((float) hitRNA/k) < 0.95)) {
      return AAType;
    }

    /* Computes the greater value between DNA's nucleotides and RNA's nucleotides */
    else if(hitRNA > hitDNA)
      gRNA++;

    else
      gDNA++;
 }
  /* Return the datatype with the greater value */
  if(gRNA > gDNA)
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
