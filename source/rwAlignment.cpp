/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 

    trimAl v1.3: a tool for automated alignment trimming in large-scale 
                 phylogenetics analyses.

    readAl v1.3: a tool for automated alignment conversion among different
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

#include "alignment.h"
#include "utils.h"

extern int errno;
#include <errno.h>
#include <ctype.h>
#include <string>

#define DELIMITERS     "   \n"
#define OTHDELIMITERS  "   \n,:*"
#define OTH2DELIMITERS "   \n,:;*"
#define PHYLIPDISTANCE 10

using namespace std;

bool alignment::fillMatrices(bool aligned) {
  int i, j;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  residuesNumber = new int[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(i = 0; i < sequenNumber; i++)
    residuesNumber[i] = sequences[i].size();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(i = 1; i < sequenNumber; i++)
    if(residuesNumber[i] != residuesNumber[i-1])
    break;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
   if (i != sequenNumber) isAligned = false;
   else                   isAligned = true;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** Fill some info ***** ***** ***** */  
  if(residNumber == 0)
    residNumber = residuesNumber[0];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** Check content ***** ***** ***** */
  if(aligned) {
    for(i = 0; i < sequenNumber; i++) {
      if(residuesNumber[i] != residNumber) {
        cerr << endl << "ERROR: The sequence \"" << seqsName[i] << "\" (" << residuesNumber[i] 
             << ") does not have the same number of residues fixed by the alignment (" << residNumber << ").";
        return false;
      }
    }
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if((aligned) || (isAligned)) {

    /* Asign its position to each column. */
    saveResidues = new int[residNumber];
    for(i = 0; i < residNumber; i++)
    saveResidues[i] = i;

    /* Asign its position to each sequence. */
    saveSequences = new int[sequenNumber];
    for(i = 0; i < sequenNumber; i++)
    saveSequences[i] = i;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(i = 0; i < sequenNumber; i++)
    for(j = 0; j < residuesNumber[i]; j++)
      if((!isalpha(sequences[i][j])) && (!ispunct(sequences[i][j]))) {
        cerr << endl << "ERROR: The sequence \"" << seqsName[i] << "\" has an unknow (" << sequences[i][j]
             <<") character.";
        return false;
      }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  return true;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

int alignment::formatInputAlignment(char *alignmentFile) {

  char c, *firstWord , *line = NULL;
  int format = 0, blocks = 0;
  ifstream file;
  string nline;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file)) return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Read lines in a safer way */
  line = utils::readLine(file);
  firstWord = strtok(line, OTHDELIMITERS);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Clustal Format */
  if((!strcmp(firstWord, "CLUSTAL")) || (!strcmp(firstWord, "clustal")))
    format = 1;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* NBRF/PIR Format */
  else if(firstWord[0] == '>' && firstWord[3] == ';')
    format = 3;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Fasta Format */
  else if(firstWord[0] == '>')
    format = 8;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Nexus Format */
  else if((!strcmp(firstWord, "#NEXUS")) || (!strcmp(firstWord, "#nexus")))
    format = 17;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Mega Format */
  else if((!strcmp(firstWord, "#MEGA")) || (!strcmp(firstWord, "#mega"))) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    blocks = 0;
    do { file.read(&c, 1); } while((c != '#') && (!file.eof()));
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    do {
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      while((c != '\n') && (!file.eof())) file.read(&c, 1);
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      file.read(&c, 1); 
      if(c == '#') blocks++;
    } while((c != '\n') && (!file.eof()));
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* Mega NonInterleaved */
    if(!blocks) format = 22;

    /* Mega Interleaved */
    else format = 21;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Phylip Format */
  else {
    sequenNumber = atoi(firstWord);
    firstWord = strtok(NULL, DELIMITERS);

    if(firstWord != NULL)
      residNumber = atoi(firstWord);

    if((sequenNumber == 1) && (residNumber != 0))
	  format = 12;

    else if((sequenNumber != 0) && (residNumber != 0)) {

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* Read lines in a safer way */
      delete [] line;
      line = utils::readLine(file);
      firstWord = strtok(line, DELIMITERS);
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      if(firstWord != NULL) blocks = 1;

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      while(true) {
        firstWord = strtok(NULL, DELIMITERS);
        if(firstWord != NULL) blocks++;
        else break;
      }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* Read lines in a safer way */
      delete [] line;
      line = utils::readLine(file);
      firstWord = strtok(line, DELIMITERS);
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      if(firstWord != NULL) blocks--;

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      while(true) {
        firstWord = strtok(NULL, DELIMITERS);
        if(firstWord != NULL) blocks--;
        else break;
      }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* Phylip Format */
      if(!blocks) format = 12;

      /* Phylip3.2 Format */
      else format = 11;
    }
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Close the input file and delete dinamic memory */
  file.close();
  delete [] line;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Return the input alignment format */
  return format;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}


bool alignment::loadPhylipAlignment(char *alignmentFile) {

  char *str, *line = NULL;
  ifstream file;
  int i;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file)) return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* We store the file name */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Read lines in a safer way */
  line = utils::readLine(file);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  str = strtok(line, DELIMITERS);
  if(str != NULL) sequenNumber = atoi(str);
  else return false;

  str = strtok(NULL, DELIMITERS);
  if(str != NULL) residNumber = atoi(str);
  else return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the parameters */
  if((sequenNumber == 0) || (residNumber == 0)) 
    return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate memory  */
  sequences  = new string[sequenNumber];
  seqsName   = new string[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Read the first block of lines that contains two tokens */
  for(i = 0; i < sequenNumber; i++) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Read lines in a safer way */
    delete [] line;
    line = utils::readLine(file);
    if(file.eof()) break;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* First token: Sequence name */
    str = strtok(line, DELIMITERS);
    seqsName[i].append(str, strlen(str));

    /* Next token: first N aminoacid block */
    str = strtok(NULL, DELIMITERS);
    sequences[i].append(str, strlen(str));
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* we trim the line of spaces and tabulations */
    /* and store it in the alignment matrix       */
    while(true) {
      str = strtok(NULL, DELIMITERS);
      if(str != NULL)
        sequences[i].append(str, strlen(str));
      else break;
    }
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Read the next aminoacids blocks */
  while(!file.eof()) {
    for(i = 0; i < sequenNumber; i++) {

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* Read lines in a safe way */
      delete [] line;
      line = utils::readLine(file);
      if(file.eof()) break;
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* we append this sequence fragment to the previous */
      /* stored sequence fragment */
      str = strtok(line, DELIMITERS);
      sequences[i].append(str, strlen(str));
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* we trim the line of spaces and tabulations */
      /* and append it to the sequence */
      while(true) {
        str = strtok(NULL, DELIMITERS);
        if(str != NULL) 
          sequences[i].append(str, strlen(str));
        else break;
      }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Close the input file and delete dinamic memory */
  file.close();
  delete [] line;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the matrix's content */
  return fillMatrices(true);
}

bool alignment::loadPhylip3_2Alignment(char *alignmentFile) {

  char c, *str, *line = NULL;
  int i, firstLine = true;
  ifstream file;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file)) return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* We store the file name */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Read lines in a safer way */
  line = utils::readLine(file);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  str = strtok(line, DELIMITERS);
  if(str != NULL) sequenNumber = atoi(str);
  else return false;

  str = strtok(NULL, DELIMITERS);
  if(str != NULL) residNumber = atoi(str);
  else return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the parameters */
  if((sequenNumber == 0) || (residNumber == 0)) 
    return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate memory  */
  sequences  = new string[sequenNumber];
  seqsName   = new string[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  i = 0;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  do {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Read lines in a safer way */
    delete [] line;
    line = utils::readLine(file);
    if(file.eof()) break;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* First block: Sequence Name + Sequence fragment */
    str = strtok(line, OTHDELIMITERS);
    if(firstLine) {
      firstLine = false;
      seqsName[i].append(str, strlen(str));
      str = strtok(NULL, OTHDELIMITERS);
    }

    /* Sequence fragment */
    while(true) {
      if(str == NULL) break;
      else sequences[i].append(str, strlen(str));
      str = strtok(NULL, OTHDELIMITERS);
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* to detect a new sequence, we read next token and */
    /* if it is a blank space then there is a new sequence */
    file.read(&c, 1); 
    if(c == '\n') {
      firstLine = true;
      i++; 
    }

    /* Move one position back to continue */
    /* with the reading flow */
    file.clear();
    file.seekg (-1, ios_base::cur);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

  } while(!file.eof());

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Close the input file and delete dinamic memory */
  file.close();
  delete [] line;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the matrix's content */
  return fillMatrices(true);
}

bool alignment::loadFastaAlignment(char *alignmentFile) {

  char *str, *line = NULL;
  ifstream file;
  int i;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file)) return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */   

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* We store the file name */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  while(true) {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Read lines in a safe way */
    line = utils::readLine(file);
    if(file.eof()) break;
    str = strtok(line, DELIMITERS);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if(str[0] == '>') sequenNumber++;
    delete [] line;
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Finish to preprocess the input file. */
  file.clear();
  file.seekg(0);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate memory for sequenNumber names vector */
  seqsName  = new string[sequenNumber];
  sequences = new string[sequenNumber];

  for(i = -1; i < sequenNumber; ) {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Read lines in a safe way */
    line = utils::readLine(file);
    if(file.eof()) break;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    str = strtok(line, OTHDELIMITERS);

    /* Sequence Name */
    if(str[0] == '>') {
      if(strlen(str) == 1) str = strtok(NULL, OTHDELIMITERS);
      else str = str + 1;
      seqsName[++i].append(str, strlen(str));
    }

    /* Sequence */
    else {
      while(true) {
        if(str == NULL) break;
        else sequences[i].append(str, strlen(str));
        str = strtok(NULL, DELIMITERS);
      }
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    delete [] line;
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Close the input file */
  file.close();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the matrix's content */
  return fillMatrices(false);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

bool alignment::loadClustalAlignment(char *alignmentFile) {

  int i, state, length = 0, store, firstBlock = true;
  char c, *str, *line = NULL;
  ifstream file;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file)) return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* We store the file name */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");
  /* ***** ***** ***** ***** ***** ***** ***** ***** */   

  /* ***** ***** ***** ***** Title ***** ***** ***** ***** */
  /* We ignore the title line. We are only interested in the */
  /* sequences number from the input alignment */
  do { file.read(&c, 1); } while((c != '\n') && (!file.eof()));
  do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  file.seekg(-1, ios::cur);
  /* ***** ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** ***** */
  /* We want to read the first block only to know the */
  /* sequences number */
  while(true) {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Read lines in a safe way */
    line = utils::readLine2(file);
    if(file.eof()) break;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    state = (int) strlen(line);
    for(i = 0, length = 0; i < state; i++)
      if((isalpha(line[i])) || (line[i] == '-')) length++;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if(!length) break;
    sequenNumber++;
    delete [] line;
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Finish to preprocess the input file. */
  file.clear();
  file.seekg(0);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** Allocate memory ***** ***** ***** */
  seqsName  = new string[sequenNumber];
  sequences = new string[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Read the title line and store it */
  line = utils::readLine(file);
  aligInfo.append(line, strlen(line));
  delete [] line;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  file.seekg(-1, ios::cur);
  /* ***** ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  firstBlock = true;
  i = 0;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  while(true) {
     if(i >= sequenNumber) i = 0;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Read line  and analyze it*/
    line = utils::readLine2(file);
    if(file.eof()) break;

    store = i;
    state = (int) strlen(line);
    for(i = 0, length = 0; i < state; i++)
      if((isalpha(line[i])) || (line[i] == '-')) length++;
    i = store;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* if there is a valid line, we need to split it in  */
    /* order to extract the sequence name, only from the */
    /* first block, as well the sequences, from each block */
    if(length) {
      str = strtok(line, OTHDELIMITERS);

      if(str != NULL) {

        /* Extract the sequence name */
        if(firstBlock)
          seqsName[i].append(str, strlen(str));
        str = strtok(NULL, OTHDELIMITERS);
        if(str != NULL)
          sequences[i].append(str, strlen(str));
        i++;
      }
    /* Start a new block in the input alignment */
    } else
        firstBlock = false;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    delete [] line;
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Close the input file */
  file.close();

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the matrix's content */
  return fillMatrices(true);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

bool alignment::loadNexusAlignment(char *alignmentFile) {

  char *frag, *str, *line = NULL;
  int i, state = false, firstBlock = true;
  ifstream file;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file)) return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* We store the file name */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");
  /* ***** ***** ***** ***** ***** ***** ***** ***** */   

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* Read lines in a safe way */
  utils::readLine(file);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  do {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */  
    /* Read lines in a safe way */
    line = utils::readLine(file);
    str = strtok(line, DELIMITERS);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */  

    /* ***** ***** ***** ***** ***** ***** ***** ***** */  
    if(str != NULL) {

      /* ***** ***** ***** ***** ***** ***** ***** ***** */  
      for(i = 0; i < (int) strlen(str); i++)
        str[i] = toupper(str[i]);

      /* ***** ***** ***** ***** ***** ***** ***** ***** */  
      if(!strcmp(str, "MATRIX")) break;

      /* ***** ***** ***** ***** ***** ***** ***** ***** */  
      if(!strcmp(str, "BEGIN")) state = true;

      /* ***** ***** ***** ***** ***** ***** ***** ***** */  
      if((!strcmp(str, "DIMENSIONS")) && state) {
        str = strtok(NULL, DELIMITERS);
        frag = strtok(NULL, DELIMITERS);

        str = strtok(str, "=;");
        sequenNumber = atoi(strtok(NULL, "=;"));

        frag = strtok(frag, "=;");
        residNumber = atoi(strtok(NULL, "=;"));
      }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }

  } while(!file.eof());

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* Check correct read parameters */
  if(strcmp(str, "MATRIX") || (sequenNumber == 0) || (residNumber == 0))
    return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** Allocate memory ***** ***** ***** */
  seqsName  = new string[sequenNumber];
  sequences = new string[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  i = 0;
  while(true) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */  
    /* Read lines in a safe way */
    line = utils::readLine(file);
    if(file.eof()) break;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */  

    /* ***** ***** ***** ***** ***** ***** ***** ***** */  
    if((!strncmp(line, "end;", 4)) || (!strncmp(line, "END;", 4))) 
      break;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */  

    /* ***** ***** ***** ***** ***** ***** ***** ***** */  
    str = strtok(line, OTH2DELIMITERS);
    if(str != NULL) {

      /* Store the sequence name, only from the first block */
      if(firstBlock)
        seqsName[i].append(str, strlen(str));

      /* Store the sequence */
      while(true) {
        str = strtok(NULL, OTH2DELIMITERS);
        if(str != NULL)
          sequences[i].append(str, strlen(str));
        else break;
      }
      i++;
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if(i >= sequenNumber) {
      firstBlock = false;
      i = 0;
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    delete [] line;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Close the input file */
  file.close();

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the matrix's content */
  return fillMatrices(true);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

bool alignment::loadNBRF_PirAlignment(char *alignmentFile) {

  int i, firstLine = true, seqLines = false;
  char *str, *line = NULL;
  ifstream file;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file)) return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* We store the file name */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  while(true) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */  
    line = utils::readLine(file);
    if(file.eof()) break;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */  

    /* ***** ***** ***** ***** ***** ***** ***** ***** */  
    if((line[0] == '>') && (line[3] == ';'))
      sequenNumber++;
    delete [] line;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }

  /* Finish to preprocess the input file. */
  file.clear();
  file.seekg(0);

  /* ***** ***** ***** Allocate memory ***** ***** ***** */
  sequences = new string[sequenNumber];
  seqsName  = new string[sequenNumber];
  seqsInfo  = new string[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  i = 0;
  while(true) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */  
    line = utils::readLine(file);
    if(file.eof()) break;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */  

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if((line[0] == '>') && (line[3] == ';') && (firstLine)) {
      firstLine = false;

      str = strtok(line, ">;"); 
      seqsInfo[i].append(str, strlen(str));

      str = strtok(NULL, ">;");
      seqsName[i].append(str, strlen(str));
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    else if((!firstLine) && (!seqLines)) 
      seqLines = true;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    else if(seqLines) {
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      if(line[strlen(line)-1] == '*') {
        seqLines = false; 
        firstLine = true;
      }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      str = strtok(line, OTHDELIMITERS);
      while(true) {
        if(str == NULL) break;
        sequences[i].append(str, strlen(str));
        str = strtok(NULL, OTHDELIMITERS);
      }
    }

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if((firstLine) && (!seqLines))
      i++;
    delete [] line;
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Close the input file */
  file.close();

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the matrix's content */
  return fillMatrices(false);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}


bool alignment::loadMegaInterleavedAlignment(char *alignmentFile) {

  char c, *str, *line = NULL;
  int i, firstBlock = true;
  string nline;
  ifstream file;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file)) return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* We store the file name */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  do { file.read(&c, 1); } while((c != '\n') && (!file.eof()));
  line = utils::readLine(file);
  nline.append(line, strlen(line));
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  str = strtok(line, "!: ");
  for(i = 0; i < (int) strlen(str); i++)
    str[i] = toupper(str[i]);

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* Read the input alignment title */
  if(!strcmp(str, "TITLE")) {
    filename.clear();

    if(nline[0] != '!') 
      filename = "!";
    filename += nline;
    cerr << endl << endl << filename << endl << endl;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  line = utils::readLineMEGA(file);
  aligInfo.append(line, strlen(line));
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  while((c != '#') && (!file.eof())) { file.read(&c, 1); }
  for(nline.clear(); (c != '\n') && (!file.eof()); file.read(&c, 1))
    nline += c;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Count the sequences number from the input alignment */
  while(!file.eof()) {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if(nline[0] != '!') {
      line = utils::trimLine(nline);
      if(strcmp(line, "")) {
        if(line[0] == '#')  sequenNumber++;
        else break;
      }
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    nline = utils::readLine(file);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* Finish to preprocess the input file. */
  file.clear();
  file.seekg(0);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** Allocate memory ***** ***** ***** */
  seqsName  = new string[sequenNumber];
  sequences = new string[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  do { file.read(&c, 1); } while((c != '#') && (!file.eof()));
  for(nline.clear(); (c != '\n') && (!file.eof()); file.read(&c, 1))
    nline += c;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  i = 0;
  while(!file.eof()) {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if(nline[0] != '!') {
      line = utils::trimLine(nline);
      str = strtok(line, " #\n");

      /* First Block */
      if(firstBlock) seqsName[i].append(str, strlen(str));

      /* Sequence */
      while(true) {
        str = strtok(NULL, " \n");
        if(str == NULL) break;
        sequences[i].append(str, strlen(str));
      }
      i++;
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    nline = utils::readLine(file);
    if(i >= sequenNumber) {
      firstBlock = false;
      i = 0; 
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Close the input file and delete dinamic memory */
  file.close();
  delete [] line;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the matrix's content */
  return fillMatrices(true);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

bool alignment::loadMegaNonInterleavedAlignment(char *alignmentFile) {

  int i, firstLine = true;
  char c, *str, *line = NULL;
  string nline;
  ifstream file;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file)) return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* We store the file name */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  do { file.read(&c, 1); } while((c != '\n') && (!file.eof()));
  line = utils::readLine(file);
  nline.append(line, strlen(line));
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  str = strtok(line, "!: ");
  for(i = 0; i < (int) strlen(str); i++)
    str[i] = toupper(str[i]);

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* Read the input alignment title */
  if(!strcmp(str, "TITLE")) {
    filename.clear();

    if(nline[0] != '!') 
      filename = "!";
    filename += nline;
    cerr << endl << endl << filename << endl << endl;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  line = utils::readLineMEGA(file);
  aligInfo.append(line, strlen(line));
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  while((c != '#') && (!file.eof())) { file.read(&c, 1); }
  for(nline.clear(); (c != '\n') && (!file.eof()); file.read(&c, 1))
    nline += c;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Count the sequences number from the input alignment */
  while(!file.eof()) {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if(nline[0] != '!') {
      line = utils::trimLine(nline);
      if(line[0] == '#')  sequenNumber++;
      }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    nline = utils::readLine(file);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  /* Finish to preprocess the input file. */
  file.clear();
  file.seekg(0);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** Allocate memory ***** ***** ***** */
  seqsName  = new string[sequenNumber];
  sequences = new string[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  do { file.read(&c, 1); } while((c != '#') && (!file.eof()));
  for(nline.clear(); (c != '\n') && (!file.eof()); file.read(&c, 1))
    nline += c;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  i = 0;
  while(!file.eof()) {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if(nline[0] != '!') {
      line = utils::trimLine(nline);
      str = strtok(line, " #\n");

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* Sequence Name */
      if(firstLine) {
        seqsName[i].append(str, strlen(str));
        str = strtok(NULL, " \n");
        firstLine = false;
      }

      while(str != NULL) {
        sequences[i].append(str, strlen(str));
        str = strtok(NULL, " \n");
      }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    nline = utils::readLine(file);
    if(nline[0] == '#') {
      firstLine = true;
      i++; 
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Close the input file and delete dinamic memory */
  file.close();
  delete [] line;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check the matrix's content */
  return fillMatrices(true);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
| void alignment::alignmentToFile(ostream &file)  |
|    Private method that put the alignment on the    |
|    parameter stream in PHYLIP2               |
+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::alignmentPhylipToFile(ostream &file) {
  
  string *tmpMatrix;
  int i, j;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not ";
    cerr << "have the sequences aligned." << endl << endl;
    return;
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate memory to output alignment. */   
  tmpMatrix = new string[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {
    if(!reverse) tmpMatrix[i] = sequences[i];
    else tmpMatrix[i] = utils::getReverse(sequences[i]);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Include in the first line the sequenNumber and the aminoacids of the alignment */
  file << " " << sequenNumber << " " << residNumber << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Put on the stream the sequenNumber names and the first 60 aminoacids block */
  for(i = 0; i < sequenNumber; i++) {
    file << setw(PHYLIPDISTANCE + 3) << left << seqsName[i].substr(0, 10);
    file << tmpMatrix[i].substr(0, 60)  << endl;
  }
  file << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Put on the stream the rest of the blocks */
  for(i = 60; i < residNumber; i += 60) {
    for(j = 0; j < sequenNumber; j++) {
      file << tmpMatrix[j].substr(i, 60) << endl;
    }
    file << endl;
  }
  file << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate local memory */
  delete [] tmpMatrix;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

void alignment::alignmentPhylip3_2ToFile(ostream &file) {

  string *tmpMatrix;
  int i, j, k;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not ";
    cerr << "have the sequences aligned." << endl << endl;
    return;
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate memory to output alignment. */   
  tmpMatrix = new string[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {
    if(!reverse) tmpMatrix[i] = sequences[i];
    else tmpMatrix[i] = utils::getReverse(sequences[i]);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Include in the first line the sequenNumber and the aminoacids of the alignment */
  file << " " << sequenNumber << " " << residNumber << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Put on the stream the sequenNumber names and the first 60 aminoacids block */
  for(i = 0; i < sequenNumber; i++) {
    file << setw(PHYLIPDISTANCE + 3) << left << seqsName[i].substr(0, 10);

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    for(j = 0; j < residNumber; j += 50) {
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      for(k = j; k < residNumber && k < j + 50; k += 10)
        file << sequences[i].substr(k, 10) << " ";
      file << endl;

      if(j + 50 < residNumber) 
        file << setw(PHYLIPDISTANCE + 3) << " ";
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    file << endl;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate local memory */
  delete [] tmpMatrix;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}


void alignment::alignmentClustalToFile(ostream &file) {

  int i, j, maxLongName = 0;
  string *tmpMatrix;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not ";
    cerr << "have the sequences aligned." << endl << endl;
    return;
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate memory to output alignment. */   
  tmpMatrix = new string[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {
    if(!reverse) tmpMatrix[i] = sequences[i];
    else tmpMatrix[i] = utils::getReverse(sequences[i]);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if((aligInfo.size() != 0)  && (iformat == oformat))
    file << aligInfo;
  else 
    file << "CLUSTAL W (1.8) multiple sequence alignment";
  file << endl << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Look for the maximum sequenNumber name size */
  for(i = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size()); 
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Print the sequence name as well as 60 residue blocks */
  for(j = 0; j < residNumber; j += 60) {
    file << endl;
    for(i = 0; i < sequenNumber; i++) {
      file << setw(maxLongName + 5) << left << seqsName[i];
      file << tmpMatrix[i].substr(j, 60) << endl;
    }
    file << endl;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate local memory */
  delete [] tmpMatrix;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}


void alignment::alignmentNBRF_PirToFile(ostream &file) {

  string *tmpMatrix;
  int i, j, k;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate memory to output alignment. */   
  tmpMatrix = new string[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {
    if(!reverse) tmpMatrix[i] = sequences[i];
    else tmpMatrix[i] = utils::getReverse(sequences[i]);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(i = 0; i < sequenNumber; i++) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Write the sequence name as well as the sequence type */
    file << ">";
/*
    if((seqsInfo != NULL) && (iformat == oformat))
      file << seqsInfo[i];

    else {
      getTypeAlignment();
      switch(dataType) {
        case DNAType: file << "DL"; break;
        case RNAType: file << "RL"; break;
        case AAType:  file << "P1"; break;
      }
    }
*/
    file << ";" << seqsName[i] << endl;
    file << seqsName[i] << " " << residuesNumber[i];
    file << " bases" << endl;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Write the sequence */
    for(j = 0; j < residuesNumber[i]; j += 50) {
      for(k = j; k < residuesNumber[i] && k < j + 50; k += 10)
        file << " " << tmpMatrix[i].substr(k, 10);

      if((j + 50) >= residNumber)
        file << "*"; 
      file << endl;
    }
    file << endl;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate local memory */
  delete [] tmpMatrix;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

void alignment::alignmentFastaToFile(ostream &file) {

  string *tmpMatrix;
  int i, j;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate memory to output alignment. */   
  tmpMatrix = new string[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {
    if(!reverse) tmpMatrix[i] = sequences[i];
    else tmpMatrix[i] = utils::getReverse(sequences[i]);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(i = 0; i < sequenNumber; i++) {
    /* Sequence Name */
    file << ">" << seqsName[i] << " " << dec << residuesNumber[i] << " bp" << endl;
    /* Sequence Residues */
    for(j = 0; j < residuesNumber[i]; j += 60)
      file << tmpMatrix[i].substr(j, 60) << endl;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate local memory */
  delete [] tmpMatrix;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}



void alignment::alignmentNexusToFile(ostream &file) {

  int i, j, k, maxLongName = 0;
  string *tmpMatrix;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not ";
    cerr << "have the sequences aligned." << endl << endl;
    return;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate memory to output alignment. */   
  tmpMatrix = new string[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {
    if(!reverse) tmpMatrix[i] = sequences[i];
    else tmpMatrix[i] = utils::getReverse(sequences[i]);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  file << "#NEXUS" << endl;
  file << "BEGIN DATA;" << endl;
  file << " DIMENSIONS NTAX=" << sequenNumber << " NCHAR=" << residNumber <<";" << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  getTypeAlignment();

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  file << " FORMAT";
  switch(dataType) {
    case DNAType: file << " DATATYPE=DNA "; break;
    case RNAType: file << " DATATYPE=RNA "; break;
    case AAType:  file << " DATATYPE=PROTEIN "; break;
  }
  file << "INTERLEAVE=yes GAP=-;" << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Print the sequence name and the residues number 
     for each sequence */
  for(i = 0; i < sequenNumber; i++)
    file << "[Name: " << seqsName[i] << "      Len: " << residNumber << " Check: 0]" << endl;
  file << endl << "MATRIX" << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Look for the maximum sequenNumber name size */
  for(i = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size()); 
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Put on the stream the sequenNumber names and the first 60 aminoacids block */
  for(j = 0; j < residNumber; j += 50) {
    for(i = 0; i < sequenNumber; i++) {
      file << setw(maxLongName + 4) << left << seqsName[i];
      for(k = j; k < (j + 50) && k < residNumber; k += 10)
        file << " " << sequences[i].substr(k, 10);
      file << endl;
    }
    file << endl;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  file << endl << ";" << endl << "END;" << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate local memory */
  delete [] tmpMatrix;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

void alignment::alignmentMegaToFile(ostream &file) {

  int i, j, k;
  string *tmpMatrix;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not ";
    cerr << "have the sequences aligned." << endl << endl;
    return;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate memory to output alignment. */   
  tmpMatrix = new string[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {
    if(!reverse) tmpMatrix[i] = sequences[i];
    else tmpMatrix[i] = utils::getReverse(sequences[i]);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Phylemon Webserver Version */
  /* file << "#MEGA" << endl << "Alignment file" << endl; */
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Standard Version */
  file << "#MEGA" << endl << filename << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  getTypeAlignment();

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  switch(dataType) {
    case DNAType:
      file << "!Format DataType=DNA ";
      break;
    case RNAType:
      file << "!Format DataType=RNA ";
      break;
    case AAType:
      file << "!Format DataType=protein ";
      break;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  file << "NSeqs=" << sequenNumber << " Nsites=";
  file << residNumber << " indel=- CodeTable=Standard;";
  file << endl << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(i = 0; i < sequenNumber; i++) {
    file << "#" << seqsName[i] << endl;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
    for(j = 0; j < residNumber; j += 50) {
      for(k = j; ((k < residNumber) && (k < j + 50)); k += 10)
        file << tmpMatrix[i].substr(k, 10) << " ";
      file << endl;
    }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
    file << endl;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate local memory */
  delete [] tmpMatrix;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

bool alignment::alignmentSummaryHTML(char *destFile, int residues, int seqs, int *selectedRes, int *selectedSeq) {

  int i, j, k, maxLongName = 0;
  bool *res, *seq;
  string tmpName;
  ofstream file;

  
  /* ***** ***** ***** ***** ***** ***** ***** ***** */  
  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not ";
    cerr << "have the sequences aligned." << endl << endl;
    return false;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Open file and check the operations */
  file.open(destFile);
  if(!file) return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Look for the maximum sequenNumber name size */
  for(i = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size()); 
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  res = new bool[residNumber];
  for(i = 0; i < residNumber; i++) res[i] = false;
  for(i = 0; i < residues; i++) res[selectedRes[i]] = true;

  seq = new bool[sequenNumber];
  for(i = 0; i < sequenNumber; i++) seq[i] = false;
  for(i = 0; i < seqs; i++) seq[selectedSeq[i]] = true;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* File Header */
  file << "<!DOCTYPE html>" << endl << "<html><head>" << endl;
  file << "    <meta http-equiv=\"Content-Type\" content=\"text/html;charset=ISO-8859-1\" />" << endl;
  file << "    <title>trimAl v1.3 Summary</title>" << endl;

  file << "    <style type=\"text/css\">" << endl;
  file << "    .sel  { background-color: #C9C9C9; }\n";
  file << "    </style>\n  </head>\n\n";
  file << "  <body>\n";

  file << "  <pre>" << endl;
  file << "  <span class=sel>Selected Residue / Sequence</span>" << endl;
  file << "  Deleted Residue / Sequence";
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(j = 0; j < residNumber; j += 120) {

    file << endl;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Print the sequence number identifier */
    for(i = j + 1; ((i <= residNumber) && (i <= (j + 10))); i++)
      if(!((i + 1) % 10)) file << setw(maxLongName + 19) << right << (i + 1) << " ";

    for(i = j + 11; ((i <= residNumber) && (i <= (j + 120))); i++)
      if(!((i + 1) % 10)) file << setw(10) << right << (i + 1) << " ";
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* A pretty format */
    file << endl << setw(maxLongName + 10);
    for(i = j + 1; ((i <= residNumber) && (i <= (j + 120))); i++)
      if(!(i % 10)) file << "+" << " ";
      else file << "=";
    file << endl;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    for(i = 0; i < sequenNumber; i++) {
      tmpName = seqsName[i] + "</span>";
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      if(seq[i]) file << "    <span class=sel>";
      else       file << "    <span>";
    file << setw(maxLongName + 12) << left << tmpName;
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      for(k = j; ((k < residNumber) && (k < (j + 120))); k++) {
        if((seq[i]) && (res[k])) 
          file << "<span class=sel>" << sequences[i][k] << "</span>";
        else
          file << sequences[i][k];
    if(!((k + 1) % 10)) file << " ";
      }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      file << endl;
      }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }
  
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  file << "    </pre>" << endl;
  file << "  </body>" << endl << "</html>" << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  file.close();
  delete [] seq;
  delete [] res;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  return true;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

void alignment::getSequences(ostream &file) {

  string *tmpMatrix;
  int i, j;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate memory to output alignment. */   
  tmpMatrix = new string[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {
    if(!reverse) tmpMatrix[i] = utils::removeCharacter('-', sequences[i]);
    else tmpMatrix[i] = utils::getReverse(utils::removeCharacter('-', sequences[i]));
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(i = 0; i < sequenNumber; i++) {
    file << ">" << seqsName[i] << endl;
    for(j = 0; j < (int) tmpMatrix[i].size(); j += 60)
      file << tmpMatrix[i].substr(j, 60) << endl;
    file << endl;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  file << endl;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate local memory */
  delete [] tmpMatrix;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}
