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
#include <exception>

#define DELIMITERS     "   \t\n"
#define OTHDELIMITERS  "   \t\n,:*"
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
   if (i != sequenNumber)
    isAligned = false;
   else
    isAligned = true;
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

  char c, *firstWord = NULL, *line = NULL;
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
  try {
  line = utils::readLine(file);
  firstWord = strtok(line, OTHDELIMITERS);
  } catch (exception& e) {
    cerr << "Standard exception: " << e.what() << endl;
  }

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

/* *****************************************************************************
 *
 * START - Refactored code
 *
 * ************************************************************************** */

bool alignment::loadPhylipAlignment(char *alignmentFile) {

  char *str, *line = NULL;
  ifstream file;
  int i;

  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file))
    return false;

  /* Store some data about filename for possible uses in other formats */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");

  /* Read the first line in a safer way */
  line = utils::readLine(file);
  if (line == NULL)
    return false;

  /* Read the input sequences and residues for each sequence numbers */
  str = strtok(line, DELIMITERS);
  sequenNumber = 0;
  if(str != NULL)
    sequenNumber = atoi(str);

  str = strtok(NULL, DELIMITERS);
  residNumber = 0;
  if(str != NULL)
    residNumber = atoi(str);

  /* If something is wrong about the sequences or/and residues number,
   * return an error to warn about that */
  if((sequenNumber == 0) || (residNumber == 0))
    return false;

  /* Allocate memory  for the input data */
  sequences  = new string[sequenNumber];
  seqsName   = new string[sequenNumber];

  /* Read the lines block containing the sequences name + first fragment */
  for(i = 0; i < sequenNumber; i++) {

    /* Read lines in a safer way. Destroy previous stored information */
    if (line != NULL)
      delete [] line;

    line = utils::readLine(file);
    /* It the input line/s are blank lines, skip the loop iteration  */
    if(line == NULL)
      continue;

    /* First token: Sequence name */
    str = strtok(line, DELIMITERS);
    seqsName[i].append(str, strlen(str));

    /* Trim the rest of the line from blank spaces, tabs, etc and store it */
    str = strtok(NULL, DELIMITERS);
    while(str != NULL) {
      sequences[i].append(str, strlen(str));
      str = strtok(NULL, DELIMITERS);
    }
  }

  /* Read the rest of the input file */
  while(!file.eof()) {
    /* Try to get for each sequences its corresponding residues */
    for(i = 0; i < sequenNumber; i++) {
      /* Read lines in a safer way. Destroy previous stored information */
      if (line != NULL)
        delete [] line;

      line = utils::readLine(file);
      /* It the input line/s are blank lines, skip the loop iteration  */
      if(line == NULL)
        continue;

      /* Remove from the current line non-printable characters and add fragments
       * to previous stored sequence */
      str = strtok(line, DELIMITERS);
      while(str != NULL) {
        sequences[i].append(str, strlen(str));
        str = strtok(NULL, DELIMITERS);
      }
    }
  }

  /* Close the input file and delete dinamic memory */
  file.close();
  if (line != NULL)
    delete [] line;

  /* Check the matrix's content */
  return fillMatrices(true);
}

bool alignment::loadPhylip3_2Alignment(char *alignmentFile) {

  int i, blocksFirstLine, firstLine = true;
  char *str, *line = NULL;
  ifstream file;

  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file))
    return false;

  /* Store the file name for futher format conversion*/
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");

  /* Read the first line in a safer way */
  line = utils::readLine(file);
  if (line == NULL)
    return false;

  /* Get the sequences and residues numbers. If there is any mistake,
   * return a FALSE value to warn about the possible error */
  str = strtok(line, DELIMITERS);
  sequenNumber = 0;
  if(str != NULL)
    sequenNumber = atoi(str);

  str = strtok(NULL, DELIMITERS);
  residNumber = 0;
  if(str != NULL)
    residNumber = atoi(str);

  if((sequenNumber == 0) || (residNumber == 0))
    return false;

  /* Reserve memory according to the input parameters */
  sequences  = new string[sequenNumber];
  seqsName   = new string[sequenNumber];

  /* Point to the first sequence in the alignment. Since the alignment could not
   * have blank lines to separate the different sequences. Store the blocks size
   * for the first line including a sequence identifier */
  i = 0;
  blocksFirstLine = 0;

  do {
    /* Read lines in a safer way. Destroy previous stored information */
    if (line != NULL)
      delete [] line;

    line = utils::readLine(file);
    /* If there is nothing in the input line, skip the loop instructions */
    if(line == NULL)
      continue;

    str = strtok(line, OTHDELIMITERS);
    /* First block: Sequence Name + Sequence fragment. Count how many blocks
     * the first sequence line is divided. It could help to identify the
     * different sequences from the input file */
    if(firstLine) {
      seqsName[i].append(str, strlen(str));
      str = strtok(NULL, OTHDELIMITERS);
      firstLine = 1;
    }

    /* Sequence fragment */
    while(str != NULL) {
      sequences[i].append(str, strlen(str));
      str = strtok(NULL, OTHDELIMITERS);
      /* Count the blocks number for the sequences first line */
      if (firstLine)
        firstLine += 1;
    }

    /* Store the blocks number for the first sequence including the name */
    if ((blocksFirstLine == 0) and firstLine)
      blocksFirstLine = firstLine;

    /* If a false positive new sequence was detected, add stored information for
     * the current sequence to the previous one and clear the data structure
     * for the current sequence. Finally, move the sequence pointer to the
     * previous one. */
    if ((firstLine != false) and (firstLine != blocksFirstLine)) {
      sequences[i-1].append(seqsName[i]);
      seqsName[i].clear();
      sequences[i-1].append(sequences[i]);
      sequences[i].clear();
      i --;
    }

    firstLine = false;
    /* There are many ways to detect a new sequence. */
    /* One of them -experimental- is just to detect if the residues number for
     * the current entry is equal to the residues number for the whole align */
    if ((int) sequences[i].size() == residNumber) {
      firstLine = true;
      i++;
    }
  } while(!file.eof());

  /* Close the input file and delete dinamic memory */
  file.close();
  if (line != NULL)
    delete [] line;

  /* Check the matrix's content */
  return fillMatrices(true);
}

bool alignment::loadClustalAlignment(char *alignmentFile) {

  int i, seqLength, pos, firstBlock;
  char *str, *line = NULL;
  ifstream file;

  /* Check input file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file))
    return false;

  /* Store some details about input file to be used in posterior format
   * conversions */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");

  /* The first line corresponding to CLUSTAL label is ignored */
  line = utils::readLine(file);
  if (line == NULL)
    return false;

  /* Ignore blank lines before first sequence block starts */
  while(!file.eof()) {

    /* Deallocate previously used dynamic memory */
    if (line != NULL)
      delete [] line;

    /* Read lines in safe way */
    line = utils::readLine(file);

    if (line != NULL)
      break;
  }

  /* The program in only interested in the first blocks of sequences since
   * it wants to know how many sequences are in the input file */
  sequenNumber = 0;
  while(!file.eof()) {

    /* If a new line without any valid character is detected
     * means the first block is over */
    if (line == NULL)
      break;

    /* Count how many times standard characters as well as
     * gap symbol "-" is detected in current line. */
    seqLength = (int) strlen(line);
    for(pos = 0; pos < seqLength; pos++)
      if((isalpha(line[pos])) || (line[pos] == '-'))
        break;

    /* If not standard characters are detected in current line means that
     * the program has found typical line in clustal alignment files to mark
     * some scores for some columns. In that case, the first block is over */
    if(pos == seqLength)
      break;
    sequenNumber++;

    /* Deallocate previously used dynamic memory */
    if (line != NULL)
      delete [] line;

    /* Read lines in safe way */
    line = utils::readLine(file);
  }

  /* Finish to preprocess the input file. */
  file.clear();
  file.seekg(0);

  /* Allocate memory for the input alignmet */
  seqsName  = new string[sequenNumber];
  sequences = new string[sequenNumber];

  /* Read the title line and store it */
  line = utils::readLine(file);
  if (line == NULL)
    return false;
  aligInfo.append(line, strlen(line));

  /* Ignore blank lines before first sequence block starts */
  while(!file.eof()) {

    /* Deallocate previously used dynamic memory */
    if (line != NULL)
      delete [] line;

    /* Read lines in safe way */
    line = utils::readLine(file);

    if (line != NULL)
      break;
  }


  /* Set-up sequences pointer to the first one and the flag to indicate
   * the first blocks. That flag implies that sequences names have to be
   * stored */
  i = 0;
  firstBlock = true;

  while(!file.eof()) {

    if (line == NULL) {
      /* Sometimes, clustalw files does not have any marker after first block
       * to indicate conservation between its columns residues. In that cases,
       * mark the end of first block */
      if (i == 0)
        firstBlock = false;
      /* Read current line and analyze it*/
      line = utils::readLine(file);
      continue;
    }

    /* Check whteher current line is a standard line or it is a line to mark
     * quality scores for that alignment column */
    seqLength = (int) strlen(line);
    for(pos = 0; pos < seqLength; pos++)
      if((isalpha(line[pos])) || (line[pos] == '-'))
        break;

    /* Start a new block in the input alignment */
    if (pos == seqLength) {
      firstBlock = false;

      /* Deallocate dinamic memory if it has been used before */
      if (line != NULL)
        delete [] line;

      /* Read current line and analyze it*/
      line = utils::readLine(file);

      continue;
    }

    /* If it is a standard line, split it into two parts. The first one contains
     * sequence name and the second one the residues. If the "firstBlock" flag
     * is active then store the sequence name */
    str = strtok(line, OTHDELIMITERS);
    if(str != NULL) {
      if(firstBlock)
        seqsName[i].append(str, strlen(str));
      str = strtok(NULL, OTHDELIMITERS);
      if(str != NULL)
        sequences[i].append(str, strlen(str));

      /* Move sequences pointer in a circular way */
      i = (i + 1) % sequenNumber;
    }

    /* Deallocate dinamic memory if it has been used before */
    if (line != NULL)
      delete [] line;

    /* Read current line and analyze it*/
    line = utils::readLine(file);
  }

  /* Close the input file */
  file.close();

  /* Deallocate dinamic memory */
  if (line != NULL)
    delete [] line;

  /* Check the matrix's content */
  return fillMatrices(true);
}

bool alignment::loadFastaAlignment(char *alignmentFile) {

  char *str, *line = NULL;
  ifstream file;
  int i;

  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file))
    return false;

  /* Store input file name for posterior uses in other formats */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");

  /* Compute how many sequences are in the input alignment */
  sequenNumber = 0;
  while(!file.eof()) {

    /* Deallocate previously used dinamic memory */
    if (line != NULL)
      delete [] line;

    /* Read lines in a safe way */
    line = utils::readLine(file);
    if (line == NULL)
      continue;

    /* It the line starts by ">" means that a new sequence has been found */
    str = strtok(line, DELIMITERS);
    if (str == NULL)
      continue;

    /* If a sequence name flag is detected, increase sequences counter */
    if(str[0] == '>')
      sequenNumber++;
  }

  /* Finish to preprocess the input file. */
  file.clear();
  file.seekg(0);

  /* Allocate memory for the input alignmet */
  seqsName  = new string[sequenNumber];
  sequences = new string[sequenNumber];

  for(i = -1; (i < sequenNumber) && (!file.eof()); ) {

    /* Deallocate previously used dinamic memory */
    if (line != NULL)
      delete [] line;

    /* Read lines in a safe way */
    line = utils::readLine(file);
    if (line == NULL)
      continue;

    /* Cut the current line and check whether there are valid characters */
    str = strtok(line, OTHDELIMITERS);
    if (str == NULL)
      continue;

    /* Check whether current line belongs to the current sequence
     * or it is a new one. In that case, store the sequence name */
    if(str[0] == '>') {
      /* Move sequence name pointer until a valid string name is obtained */
      do {
        str = str + 1;
      } while(strlen(str) == 0);
      seqsName[++i].append(str, strlen(str));
      continue;
    }

    /* Sequence */
    while(str != NULL) {
      sequences[i].append(str, strlen(str));
      str = strtok(NULL, DELIMITERS);
    }
  }

  /* Close the input file */
  file.close();

  /* Deallocate previously used dinamic memory */
  if (line != NULL)
    delete [] line;

  /* Check the matrix's content */
  return fillMatrices(false);
}

bool alignment::loadNexusAlignment(char *alignmentFile) {

  char *frag = NULL, *str = NULL, *line = NULL;
  int i, pos, state, firstBlock;
  ifstream file;

  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file))
    return false;

  /* Store input file name for posterior uses in other formats */
  /* We store the file name */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");

  state = false;
  do {

    /* Destroy previous assigned memory */
    if (line != NULL)
      delete [] line;

    /* Read line in a safer way */
    line = utils::readLine(file);
    if (line == NULL)
      continue;

    /* Discard line where there is not information */
    str = strtok(line, DELIMITERS);
    if(str == NULL)
      continue;

    /* If the line has any kind of information, try to catch it */
    /* Firstly, convert to capital letters the input line */
    for(i = 0; i < (int) strlen(str); i++)
      str[i] = toupper(str[i]);

    /* ... and then compare it again specific tags */
    if(!strcmp(str, "BEGIN"))
      state = true;

    else if(!strcmp(str, "MATRIX"))
      break;

    /* Store information about input format file */
    else if(!strcmp(str, "FORMAT")) {
      str = strtok(NULL, DELIMITERS);
      while(str != NULL) {
        aligInfo.append(str, strlen(str));
        aligInfo.append(" ", strlen(" "));
        str = strtok(NULL, DELIMITERS);
      }
    }

    /* In this case, try to get matrix dimensions */
    else if((!strcmp(str, "DIMENSIONS")) && state) {
      str = strtok(NULL, DELIMITERS);
      frag = strtok(NULL, DELIMITERS);
      str = strtok(str, "=;");
      sequenNumber = atoi(strtok(NULL, "=;"));
      frag = strtok(frag, "=;");
      residNumber = atoi(strtok(NULL, "=;"));
    }
  } while(!file.eof());

  /* Check all parameters */
  if(strcmp(str, "MATRIX") || (sequenNumber == 0) || (residNumber == 0))
    return false;

  /* Allocate memory for the input alignmet */
  seqsName  = new string[sequenNumber];
  sequences = new string[sequenNumber];

  pos = 0;
  state = false;
  firstBlock = true;

  while(!file.eof()) {
    /* Destroy previous assigned memory */
    if (line != NULL)
      delete [] line;

    /* Read line in a safer way */
    line = utils::readLine(file);
    if (line == NULL)
      continue;

    /* Discard any comments from input file */
    for(i = 0; i < (int) strlen(line); i++) {
      if (line[i] == '[')
        state = true;
      else if (line[i] == ']' && state) {
        state = false;
        break;
      }
    }

    /* If there is a multi-line comments, skip it as well */
    if ((state) || (not state && i != (int) strlen(line)))
     continue;

    /* Check for a specific tag indicating matrix end */
    if((!strncmp(line, "end;", 4)) || (!strncmp(line, "END;", 4)))
      break;

    /* Split input line and check it if it is valid */
    str = strtok(line, OTH2DELIMITERS);
    if (str == NULL)
      continue;

    /* Store the sequence name, only from the first block */
    if(firstBlock)
      seqsName[pos].append(str, strlen(str));

    /* Store rest of line as part of sequence */
    str = strtok(NULL, OTH2DELIMITERS);
    while(str != NULL) {
      sequences[pos].append(str, strlen(str));
      str = strtok(NULL, OTH2DELIMITERS);
    }

    /* move sequences pointer to next one. It if it is last one, move it to
     * the beginning and set the first block to false for avoiding to rewrite
     * sequences name */
    pos = (pos + 1) % sequenNumber;
    if (not pos)
      firstBlock = false;
  }

  /* Deallocate memory */
  if (line != NULL)
    delete [] line;

  /* Close the input file */
  file.close();

  /* Check the matrix's content */
  return fillMatrices(true);
}

bool alignment::loadMegaNonInterleavedAlignment(char *alignmentFile) {

  int i, firstLine = true;
  char *frag = NULL, *str = NULL, *line = NULL;
  ifstream file;

  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file))
    return false;

  /* Filename is stored as a title for MEGA input alignment.
   * If it is detected later a label "TITLE" in input file, this information
   * will be replaced for that one */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");

  /* Skip first line */
  line = utils::readLine(file);

  /* Try to get input alignment information */
  while(!file.eof()) {

    /* Destroy previously allocated memory */
    if (line != NULL)
      delete [] line;

    /* Read a new line in a safe way */
    line = utils::readLine(file);
    if (line == NULL)
      continue;

    /* If a sequence name flag is found, go out from getting information loop */
    if(!strncmp(line, "#", 1))
      break;

    /* Destroy previously allocated memory */
    if (frag != NULL)
      delete [] frag;

    /* Create a local copy from input line */
    frag = new char[strlen(line) + 1];
    strcpy(frag, line);

    /* Split input line copy into pieces and analize it
     * looking for specific labels */
    str = strtok(frag, "!: ");
    for(i = 0; i < (int) strlen(str); i++)
      str[i] = toupper(str[i]);

    /* If TITLE label is found, replace previously stored information with
     * this info */
    if(!strcmp(str, "TITLE")) {
      filename.clear();
      if(strncmp(line, "!", 1))
        filename += "!";
      filename += line;
    }

    /* If FORMAT label is found, try to get some details from input file */
    else if(!strcmp(str, "FORMAT"))
      aligInfo.append(line, strlen(line));
  }

  /* Deallocate local memory */
  if (frag != NULL)
    delete [] frag;

  /* Count how many sequences are in input alignment file */
  do {

    /* Check whether input line is valid or not */
    if (line == NULL) {
      line = utils::readLine(file);
      continue;
    }

    /* If current line starts by a # means that it is a sequence name */
    if (!strncmp(line, "#", 1))
      sequenNumber++;

    /* Destroy previously allocated memory */
    if (line != NULL)
      delete [] line;

    /* Read a new line in a safe way */
    line = utils::readLine(file);

  } while(!file.eof());

  /* Move file pointer to the beginner */
  file.clear();
  file.seekg(0);

  /* Allocate memory */
  seqsName  = new string[sequenNumber];
  sequences = new string[sequenNumber];

  /* Skip first line */
  line = utils::readLine(file);

  /* Skip lines until first sequence name is found */
  while(!file.eof()) {

    /* Destroy previously allocated memory */
    if (line != NULL)
      delete [] line;

    /* Read a new line in a safe way */
    line = utils::readLine(file);
    if (line == NULL)
      continue;

    /* If sequence name label is found, go out from loop */
    if (!strncmp(line, "#", 1))
      break;
  }

  /* First sequence is already detected so its name should be stored */
  i = -1;
  firstLine = true;

  /* This loop is a bit tricky because first sequence name has been already
   * detected, so it is necessary to process it before moving to next line.
   * That implies that lines are read at loop ends */
  while(!file.eof()) {

    /* Skip blank lines */
    if (line == NULL) {
      line = utils::readLine(file);
      continue;
    }

    /* Skip lines with comments */
    if (!strncmp(line, "!", 1)) {
      /* Deallocate memory and read a new line */
      delete [] line;
      line = utils::readLine(file);
      continue;
    }

    /* Remove comments inside of sequences and split input line */
    frag = utils::trimLine(line);

    /* Skip lines with only comments */
    if (frag == NULL) {

      /* Deallocate memory and read a new line */
      if (line != NULL)
        delete [] line;
      line = utils::readLine(file);
      continue;
    }

    /* Otherwise, split it into fragments */
    str = strtok(frag, " #\n");

    /* Sequence Name */
    if (!strncmp(line, "#", 1)) {
      i += 1;
      seqsName[i].append(str, strlen(str));
      str = strtok(NULL, " #\n");
    }

    /* Sequence itself */
    while(str != NULL) {
      sequences[i].append(str, strlen(str));
      str = strtok(NULL, " \n");
    }

    /* Deallocate dynamic memory */
    if (frag != NULL)
      delete [] frag;

    if (line != NULL)
      delete [] line;

    /* Read a new line in a safe way */
    line = utils::readLine(file);
  }

  /* Close input file */
  file.close();

  /* Deallocate dynamic memory */
  if (line != NULL)
    delete [] line;

  /* Check the matrix's content */
  return fillMatrices(true);
}

bool alignment::loadMegaInterleavedAlignment(char *alignmentFile) {

  char *frag = NULL, *str = NULL, *line = NULL;
  int i, firstBlock = true;
  ifstream file;

  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file))
    return false;

  /* Filename is stored as a title for MEGA input alignment.
   * If it is detected later a label "TITLE" in input file, this information
   * will be replaced for that one */
  filename.append("!Title ");
  filename.append(alignmentFile);
  filename.append(";");

  /* Skip first line */
  line = utils::readLine(file);

  /* Try to get input alignment information */
  while(!file.eof()) {

    /* Destroy previously allocated memory */
    if (line != NULL)
      delete [] line;

    /* Read a new line in a safe way */
    line = utils::readLine(file);
    if (line == NULL)
      continue;

    /* If a sequence name flag is found, go out from getting information loop */
    if(!strncmp(line, "#", 1))
      break;

    /* Create a local copy from input line */
    frag = new char[strlen(line) + 1];
    strcpy(frag, line);

    /* Split input line copy into pieces and analize it
     * looking for specific labels */
    str = strtok(frag, "!: ");
    for(i = 0; i < (int) strlen(str); i++)
      str[i] = toupper(str[i]);

    /* If TITLE label is found, replace previously stored information with
     * this info */
    if(!strcmp(str, "TITLE")) {
      filename.clear();
      if(strncmp(line, "!", 1))
        filename += "!";
      filename += line;
    }

    /* If FORMAT label is found, try to get some details from input file */
    else if(!strcmp(str, "FORMAT"))
      aligInfo.append(line, strlen(line));

    /* Destroy previously allocated memory */
    if (frag != NULL)
      delete [] frag;
  }

  /* Count how many sequences are in input file */
  while(!file.eof()) {

    /* If a sequence name flag has been detected, increase counter */
    if(!strncmp(line, "#", 1))
      sequenNumber++;

    /* Deallocate dynamic memory */
    if(line != NULL)
      delete [] line;

    /* Read lines in a safe way */
    line = utils::readLine(file);

    /* If a blank line is detected means first block of sequences is over */
    /* Then, break counting sequences loop */
    if (line == NULL)
      break;
  }

  /* Finish to preprocess the input file. */
  file.clear();
  file.seekg(0);

  /* Allocate memory */
  seqsName  = new string[sequenNumber];
  sequences = new string[sequenNumber];

  /* Skip first line */
  line = utils::readLine(file);

  /* Skip lines until first # flag is reached */
  while(!file.eof()) {

    /* Deallocate previously used dynamic memory */
    if (line != NULL)
      delete [] line;

    /* Read line in a safer way */
    line = utils::readLine(file);

    /* Determine whether a # flag has been found in current string */
    if (line != NULL)
      if(!strncmp(line, "#", 1))
        break;
  }

  /* Read sequences and get from first block, the sequences names */
  i = 0;
  firstBlock = true;

  while(!file.eof()) {

    if (line == NULL) {
    /* Read line in a safer way */
    line = utils::readLine(file);
    continue;
    }

    if (!strncmp(line, "!", 1)) {
      /* Deallocate memory and read a new line */
      delete [] line;
      line = utils::readLine(file);
      continue;
    }

    /* Trim lines from any kind of comments and split it */
    frag = utils::trimLine(line);
    str = strtok(frag, " #\n");

    /* Check whether a line fragment is valid or not */
    if (str == NULL)
      continue;

    /* Store sequences names if firstBlock flag is TRUE */
    if(firstBlock)
      seqsName[i].append(str, strlen(str));

    /* Store sequence */
    str = strtok(NULL, " \n");
    while(str != NULL) {
      sequences[i].append(str, strlen(str));
      str = strtok(NULL, " \n");
    }

    /* Deallocate previously used dynamic memory */
    if (frag != NULL)
      delete [] frag;

    if (line != NULL)
      delete [] line;

    /* Read line in a safer way */
    line = utils::readLine(file);

    i = (i + 1) % sequenNumber;
    if (i == 0)
      firstBlock = false;
  }

  /* Close input file */
  file.close();

  /* Deallocate local memory */
  if (line != NULL)
    delete [] line;

  /* Check the matrix's content */
  return fillMatrices(true);
}

/* *****************************************************************************
 *
 * END - Refactored code
 *
 * ************************************************************************** */

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

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
| void alignment::alignmentToFile(ostream &file)  |
|    Private method that put the alignment on the    |
|    parameter stream in PHYLIP2               |
+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::alignmentPhylipToFile(ostream &file) {

  int i, j, maxLongName;
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
  /* Include in the first line the sequenNumber and the aminoacids of the alignment */
  file << " " << sequenNumber << " " << residNumber << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Look for the maximum sequenNumber name size */
  maxLongName = PHYLIPDISTANCE;
  for(i = 0; (i < sequenNumber) && (!shortNames); i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Put on the stream the sequenNumber names and the first 60 aminoacids block */
  for(i = 0; i < sequenNumber; i++) {
	if(!shortNames)    file << setw(maxLongName + 3) << left << seqsName[i];
	else file << setw(maxLongName + 3) << left << seqsName[i].substr(0, 10);
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

  int i, j, k, maxLongName;
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
  /* Include in the first line the sequenNumber and the aminoacids of the alignment */
  file << " " << sequenNumber << " " << residNumber << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Look for the maximum sequenNumber name size */
  maxLongName = PHYLIPDISTANCE;
  for(i = 0; (i < sequenNumber) && (!shortNames); i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Put on the stream the sequenNumber names and the first 60 aminoacids block */
  for(i = 0; i < sequenNumber; i++) {
	if(!shortNames)    file << setw(maxLongName + 3) << left << seqsName[i];
	else file << setw(maxLongName + 3) << left << seqsName[i].substr(0, 10);

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    for(j = 0; j < residNumber; j += 50) {
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      for(k = j; k < residNumber && k < j + 50; k += 10)
        file << sequences[i].substr(k, 10) << " ";
      file << endl;

      if(j + 50 < residNumber)
        file << setw(maxLongName + 3) << " ";
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


void alignment::alignmentPhylip_PamlToFile(ostream &file) {

  int i, maxLongName;
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
  /* Include in the first line the sequenNumber and the aminoacids of the alignment */
  file << " " << sequenNumber << " " << residNumber << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Look for the maximum sequenNumber name size */
  maxLongName = PHYLIPDISTANCE;
  for(i = 0; (i < sequenNumber) && (!shortNames); i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Put on the stream the sequence name, limited to
     the 10 first characters, and the sequence */
  for(i = 0; i < sequenNumber; i++) {
	if(!shortNames)    file << setw(maxLongName + 3) << left << seqsName[i];
	else file << setw(maxLongName + 3) << left << seqsName[i].substr(0, 10);
    file << sequences[i] << endl;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  file << endl;

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
	if(!shortNames)    file << ">" << seqsName[i] << endl;
	else file << ">" << seqsName[i].substr(0, 10) << endl;
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

  /* Compute input file datatype */
  getTypeAlignment();
  switch(dataType) {
    case DNAType:
      file << "FORMAT DATATYPE=DNA INTERLEAVE=yes GAP=-";
      break;
    case RNAType:
      file << "FORMAT DATATYPE=RNA INTERLEAVE=yes GAP=-";
      break;
    case AAType:
      file << "FORMAT DATATYPE=PROTEIN INTERLEAVE=yes GAP=-";
      break;
  }

  /* Try to use information from input alignment */

  /* Endding characters like ";" are removed from input information line */
  while((int) aligInfo.find(";") != (int) string::npos)
    aligInfo.erase(aligInfo.find(";"), 1);

  i = 0;
  j = aligInfo.find(" ", i);
  /* Scan information line looking for specific tags. No all available tags
   * are taking into account */
  while(j != (int) string::npos) {

    if((aligInfo.substr(i, j - i)).compare(0, 7, "MISSING") == 0 ||
       (aligInfo.substr(i, j)).compare(0, 7, "missing") == 0)
      file << " " << (aligInfo.substr(i, j - i));

    if((aligInfo.substr(i, j)).compare(0, 9, "MATCHCHAR") == 0 ||
       (aligInfo.substr(i, j)).compare(0, 9, "matchchar") == 0)
      file << " " << (aligInfo.substr(i, j - i));

    i = j + 1;
    j = aligInfo.find(" ", i);
  }
  file << ";" << endl;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* Look for the maximum sequenNumber name size */
  for(i = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Print the sequence name and the residues number
     for each sequence */
  for(i = 0; i < sequenNumber; i++)
    file << "[Name: " << setw(maxLongName + 4) << left << seqsName[i] << "Len: " << residNumber << " Check: 0]" << endl;
  file << endl << "MATRIX" << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */

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
  file << ";" << endl << "END;" << endl;
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


bool alignment::alignmentColourHTML(ostream &file, float *GapsVector, float *SimilVector, float *ConsVector, float *IdentVector) {

  int i, j, kj, upper, k = 0, maxLongName = 0;
  string tmpColumn;

  tmpColumn.reserve(sequenNumber);

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not ";
    cerr << "have the sequences aligned." << endl << endl;
    return false;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Look for the maximum sequenNumber name size */
  for(i = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* File Header */
  file << "<!DOCTYPE html>" << endl << "<html><head>" << endl;
  file << "    <meta http-equiv=\"Content-Type\" content=\"text/html;charset=ISO-8859-1\" />" << endl;
  file << "    <title>trimAl v1.3 Summary</title>" << endl;

  file << "    <style type=\"text/css\">" << endl;
  file << "    .w  { background-color: #FFFFFF; }\n";
  file << "    .b  { background-color: #3366ff; }\n";
  file << "    .r  { background-color: #cc0000; }\n";
  file << "    .g  { background-color: #33cc00; }\n";
  file << "    .p  { background-color: #ff6666; }\n";
  file << "    .m  { background-color: #cc33cc; }\n";
  file << "    .o  { background-color: #ff9900; }\n";
  file << "    .c  { background-color: #46C7C7; }\n";
  file << "    .y  { background-color: #FFFF00; }\n";
  file << "    </style>\n  </head>\n\n";
  file << "  <body>\n";

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  file << "  <pre>" << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(j = 0, upper = 120; j < residNumber; j += 120, upper += 120) {

    file << endl;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Print the sequence number identifier */
    file << setw(maxLongName + 20) << right << (j + 10);
    for(i = j + 20; ((i <= residNumber) && (i <= upper)); i += 10)
      file << setw(10) << right << i;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* A pretty format */
    file << endl << setw(maxLongName + 10);
	for(i = j + 10; ((i < residNumber) && (i <= upper)); i += 10) {
	 for(k = i - 9; k < i; k++) file << "=";
     file << "+";
	}
	for( ; ((k < residNumber) && (k < upper)); k++) file << "=";
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    for(i = 0; i < sequenNumber; i++) {
      file << endl << setw(maxLongName + 9) << left << seqsName[i];
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      for(k = j; ((k < residNumber) && (k < upper)); k++) {
		for(kj = 0, tmpColumn.clear(); kj < sequenNumber; kj++)
		  tmpColumn += sequences[kj][k];
		file << "<span class=" << utils::determineColor(sequences[i][k], tmpColumn)
		     << ">" << sequences[i][k] << "</span>";
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      }
    }
	file << endl;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  file << "    </pre>" << endl;
  file << "  </body>" << endl << "</html>" << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  return true;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

bool alignment::alignmentColourHTML(ostream &file) {

  int i, j, kj, upper, k = 0, maxLongName = 0;
  string tmpColumn;
  char type;

  tmpColumn.reserve(sequenNumber);

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not ";
    cerr << "have the sequences aligned." << endl << endl;
    return false;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Look for the maximum sequenNumber name size */
  for(i = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* File Header */
  file << "<!DOCTYPE html>" << endl << "<html><head>" << endl;
  file << "    <meta http-equiv=\"Content-Type\" content=\"text/html;charset=ISO-8859-1\" />" << endl;
  file << "    <title>trimAl v1.3 Summary</title>" << endl;

  file << "    <style type=\"text/css\">" << endl;
  file << "    #b  { background-color: #3366ff; }\n";
  file << "    #r  { background-color: #cc0000; }\n";
  file << "    #g  { background-color: #33cc00; }\n";
  file << "    #p  { background-color: #ff6666; }\n";
  file << "    #m  { background-color: #cc33cc; }\n";
  file << "    #o  { background-color: #ff9900; }\n";
  file << "    #c  { background-color: #46C7C7; }\n";
  file << "    #y  { background-color: #FFFF00; }\n";
  file << "    </style>\n  </head>\n\n";
  file << "  <body>\n";

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  file << "  <pre>" << endl;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(j = 0, upper = 100; j < residNumber; j += 100, upper += 100) {

    file << endl;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Print the sequence number identifier */
    file << setw(maxLongName + 20) << right << (j + 10);
    for(i = j + 20; ((i <= residNumber) && (i <= upper)); i += 10)
      file << setw(10) << right << i;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* A pretty format */
    file << endl << setw(maxLongName + 10);
	for(i = j + 10; ((i < residNumber) && (i <= upper)); i += 10) {
	 for(k = i - 9; k < i; k++) file << "=";
     file << "+";
	}
	for( ; ((k < residNumber) && (k < upper)); k++) file << "=";
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    for(i = 0; i < sequenNumber; i++) {
      file << endl << setw(maxLongName + 9) << left << seqsName[i];
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      for(k = j; ((k < residNumber) && (k < upper)); k++) {
        for(kj = 0, tmpColumn.clear(); kj < sequenNumber; kj++)
          tmpColumn += sequences[kj][k];
        type = utils::determineColor(sequences[i][k], tmpColumn);
        if (type != 'w')
          file << "<span id=" << type << ">" << sequences[i][k] << "</span>";
        else
          file << sequences[i][k];
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      }
    }
	file << endl;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  file << "    </pre>" << endl;
  file << "  </body>" << endl << "</html>" << endl;
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
