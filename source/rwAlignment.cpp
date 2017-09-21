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

#include "alignment.h"
#include "defines.h"
#include "utils.h"

extern int errno;
#include <errno.h>
#include <ctype.h>
#include <string>

using namespace std;

bool alignment::fillMatrices(bool aligned) {
  /* Function to determine if a set of sequences, that can be aligned or not,
   * have been correctly load and are free of errors. */
  int i, j;

  /* Initialize some variables */
  residuesNumber = new int[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {
    residuesNumber[i] = sequences[i].size();
  }

  /* Check whether there are any unknow/no allowed character in the sequences */
  for(i = 0; i < sequenNumber; i++)
    for(j = 0; j < residuesNumber[i]; j++)
      if((!isalpha(sequences[i][j])) && (!ispunct(sequences[i][j]))) {
        cerr << endl << "ERROR: The sequence \"" << seqsName[i] << "\" has an "
          << "unknown (" << sequences[i][j] << ") character." << endl;
        return false;
      }

  /* Check whether all sequences have same size or not */
  for(i = 1; i < sequenNumber; i++)
    if(residuesNumber[i] != residuesNumber[i-1])
      break;
  /* Set an appropriate flag for indicating if sequences are aligned or not */
  isAligned = (i != sequenNumber) ? false : true;

  /* Warm about those cases where sequences should be aligned
   * and there are not */
  if (aligned and !isAligned) {
    cerr << endl << "ERROR: Sequences should be aligned (all with same length) "
      << "and there are not. Check your input alignment" << endl;
    return false;
  }

  /* Full-fill some information about input alignment */
  if(residNumber == 0)
    residNumber = residuesNumber[0];

  /* Check whether aligned sequences have the length fixed for the input alig */
  for(i = 0; (i < sequenNumber) and (aligned); i++) {
    if(residuesNumber[i] != residNumber) {
      cerr << endl << "ERROR: The sequence \"" << seqsName[i] << "\" ("
        << residuesNumber[i] << ") does not have the same number of residues "
        << "fixed by the alignment (" << residNumber << ")." << endl;
      return false;
    }
  }

  /* If the sequences are aligned, initialize some additional variables.
   * These variables will be useful for posterior analysis */
  if((aligned) || (isAligned)) {

    /* Asign its position to each column. That will be used to determine which
     * columns should be kept in output alignment after applying any method
     * and which columns should not */
    saveResidues = new int[residNumber];
    for(i = 0; i < residNumber; i++)
      saveResidues[i] = i;

    /* Asign its position to each sequence. Similar to the columns numbering
     * process, assign to each sequence its position is useful to know which
     * sequences will be in the output alignment */
    saveSequences = new int[sequenNumber];
    for(i = 0; i < sequenNumber; i++)
      saveSequences[i] = i;
  }

  /* Return an flag indicating that everything is fine */
  return true;
}

void alignment::getSequences(ostream &file) {
  /* Get only residues sequences from input alignment */

  int i, j;
  string *tmpMatrix;

  /* Allocate local memory for generating output alignment */
  tmpMatrix = new string[sequenNumber];

  /* Depending on alignment orientation: forward or reverse. Copy directly
   * sequence information or get firstly the reversed sequences and then copy
   * it into local memory. Once orientation has been determined, remove gaps
   * from resuting sequences */
  for(i = 0; i < sequenNumber; i++)
    tmpMatrix[i] = (!reverse) ? utils::removeCharacter('-', sequences[i]) : \
      utils::removeCharacter('-', utils::getReverse(sequences[i]));

  /* Print output set of sequences in FASTA format*/
  for(i = 0; i < sequenNumber; i++) {
    file << ">" << seqsName[i] << endl;
    for(j = 0; j < (int) tmpMatrix[i].size(); j += 60)
      file << tmpMatrix[i].substr(j, 60) << endl;
    file << endl;
  }

  /* Deallocate local memory */
  delete [] tmpMatrix;
}

int alignment::formatInputAlignment(char *alignmentFile) {
  /* Guess input alignment format */

  char c, *firstWord = NULL, *line = NULL;
  int format = 0, blocks = 0;
  ifstream file;
  string nline;

  /* Check the file and its content */
  file.open(alignmentFile, ifstream::in);
  if(!utils::checkFile(file))
    return false;

  /* Read first valid line in a safer way */
  do {
    line = utils::readLine(file);
  } while ((line == NULL) && (!file.eof()));

  /* If the file end is reached without a valid line, warn about it */
  if (file.eof())
    return false;

  /* Otherwise, split line */
  firstWord = strtok(line, OTHDELIMITERS);

  /* Clustal Format */
  if((!strcmp(firstWord, "CLUSTAL")) || (!strcmp(firstWord, "clustal")))
    format = 1;

  /* NBRF/PIR Format */
  else if(firstWord[0] == '>' && firstWord[3] == ';')
    format = 3;

  /* Fasta Format */
  else if(firstWord[0] == '>')
    format = 8;

  /* Nexus Format */
  else if((!strcmp(firstWord, "#NEXUS")) || (!strcmp(firstWord, "#nexus")))
    format = 17;

  /* Mega Format */
  else if((!strcmp(firstWord, "#MEGA")) || (!strcmp(firstWord, "#mega"))) {

    /* Determine specific mega format: sequential or interleaved.
     * Counting the number of blocks (set of lines starting by "#") in
     * the input file. */
    blocks = 0;
    do {
      file.read(&c, 1);
    } while((c != '#') && (!file.eof()));

    do {
      while((c != '\n') && (!file.eof()))
        file.read(&c, 1);
      file.read(&c, 1);
      if(c == '#')
        blocks++;
    } while((c != '\n') && (!file.eof()));

    /* MEGA Sequential (22) or Interleaved (21) */
    format = (!blocks) ? 22 : 21;
  }

  /* Phylip Format */
  else {
    /* Determine specific phylip format: sequential or interleaved. */

    /* Get number of sequences and residues */
    sequenNumber = atoi(firstWord);
    firstWord = strtok(NULL, DELIMITERS);
    if(firstWord != NULL)
      residNumber = atoi(firstWord);

    /* If there is only one sequence, use by default sequential format since
     * it is impossible to determine exactly which phylip format is */
    if((sequenNumber == 1) && (residNumber != 0))
      format = 12;

    /* If there are more than one sequence, analyze sequences distribution to
     * determine its format. */
    else if((sequenNumber != 0) && (residNumber != 0)) {
      blocks = 0;

      /* Read line in a safer way */
      do {
        if (line != NULL)
          delete [] line;
        line = utils::readLine(file);
      } while ((line == NULL) && (!file.eof()));

      /* If the file end is reached without a valid line, warn about it */
      if (file.eof())
        return false;

      firstWord = strtok(line, DELIMITERS);
      while(firstWord != NULL) {
        blocks++;
        firstWord = strtok(NULL, DELIMITERS);
      }

      /* Read line in a safer way */
      do {
        if (line != NULL)
          delete [] line;
        line = utils::readLine(file);
      } while ((line == NULL) && (!file.eof()));

      firstWord = strtok(line, DELIMITERS);
      while(firstWord != NULL) {
        blocks--;
        firstWord = strtok(NULL, DELIMITERS);
      }

      /* If the file end is reached without a valid line, warn about it */
      if (file.eof())
        return false;

      /* Phylip Interleaved (12) or Sequential (11) */
      format = (!blocks) ? 12 : 11;
    }
  }

  /* Close the input file and delete dinamic memory */
  file.close();
  if (line != NULL)
    delete [] line;

  /* Return the input alignment format */
  return format;
}

bool alignment::loadPhylipAlignment(char *alignmentFile) {
  /* PHYLIP/PHYLIP 4 (Sequential) file format parser */

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

  /* Read first valid line in a safer way */
  do {
    line = utils::readLine(file);
  } while ((line == NULL) && (!file.eof()));

  /* If the file end is reached without a valid line, warn about it */
  if (file.eof())
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
  i = 0;
  while((i < sequenNumber) && (!file.eof())){

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
    i++;
  }

  /* Read the rest of the input file */
  while(!file.eof()) {

    /* Try to get for each sequences its corresponding residues */
    i = 0;
    while((i < sequenNumber) && (!file.eof())) {
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
      i++;
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
  /* PHYLIP 3.2 (Interleaved) file format parser */

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

  /* Read first valid line in a safer way */
  do {
    line = utils::readLine(file);
  } while ((line == NULL) && (!file.eof()));

  /* If the file end is reached without a valid line, warn about it */
  if (file.eof())
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

    str = strtok(line, DELIMITERS);
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
  /* CLUSTAL file format parser */

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

  /* The first valid line corresponding to CLUSTAL label is ignored */
  do {
    line = utils::readLine(file);
  } while ((line == NULL) && (!file.eof()));

  /* If the file end is reached without a valid line, warn about it */
  if (file.eof())
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
  /* FASTA file format parser */

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
  seqsInfo  = new string[sequenNumber];

  for(i = -1; (i < sequenNumber) && (!file.eof()); ) {

    /* Deallocate previously used dinamic memory */
    if (line != NULL)
      delete [] line;

    /* Read lines in a safe way */
    line = utils::readLine(file);
    if (line == NULL)
      continue;

    /* Store original header fom input sequences including non-standard
     * characters */
    if (line[0] == '>')
      seqsInfo[i+1].append(&line[1], strlen(line) - 1);

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

  /* NEXUS file format parser */
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

    /* Move sequences pointer to next one. It if it is last one, move it to
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
  /* MEGA sequential file format parser */

  char *frag = NULL, *str = NULL, *line = NULL;
  ifstream file;
  int i;

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

  /* Skip first valid line */
  do {
    line = utils::readLine(file);
  } while ((line == NULL) && (!file.eof()));

  /* If the file end is reached without a valid line, warn about it */
  if (file.eof())
    return false;

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
  /* MEGA interleaved file format parser */

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

  /* Skip first valid line */
  do {
    line = utils::readLine(file);
  } while ((line == NULL) && (!file.eof()));

  /* If the file end is reached without a valid line, warn about it */
  if (file.eof())
    return false;

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

bool alignment::loadNBRF_PirAlignment(char *alignmentFile) {
  /* NBRF/PIR file format parser */

  bool seqIdLine, seqLines;
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
  sequences = new string[sequenNumber];
  seqsName  = new string[sequenNumber];
  seqsInfo  = new string[sequenNumber];

  /* Initialize some local variables */
  seqIdLine = true;
  seqLines = false;
  i = -1;

  /* Read the entire input file */
  while(!file.eof()) {

    /* Deallocate local memory */
    if (line != NULL)
      delete [] line;

    /* Read lines in a safe way */
    line = utils::readLine(file);
    if (line == NULL)
      continue;

    /* Sequence ID line.
     * Identification of these kind of lines is based on presence of ">" and ";"
     * symbols at positions 0 and 3 respectively */
    if((line[0] == '>') && (line[3] == ';') && (seqIdLine)) {
      seqIdLine = false;
      i += 1;

      /* Store information about sequence datatype */
      str = strtok(line, ">;");
      seqsInfo[i].append(str, strlen(str));

      /* and the sequence identifier itself */
      str = strtok(NULL, ">;");
      seqsName[i].append(str, strlen(str));
    }

    /* Line just after sequence Id line contains a textual description of
     * the sequence. */
    else if((!seqIdLine) && (!seqLines)) {
      seqLines = true;
      seqsInfo[i].append(line, strlen(line));
    }

    /* Sequence lines itself */
    else if (seqLines) {

      /* Check whether a sequence end symbol '*' exists in current line.
       * In that case, set appropriate flags to read a new sequence */
      if (line[strlen(line) - 1] == '*') {
        seqLines = false;
        seqIdLine = true;
      }

      /* Process line */
      str = strtok(line, OTHDELIMITERS);
      while (str != NULL) {
        sequences[i].append(str, strlen(str));
        str = strtok(NULL, OTHDELIMITERS);
      }

      /* In case the end symbol '*' has been detected, remove it */
      if(sequences[i][sequences[i].size() - 1] == '*')
        sequences[i].erase(sequences[i].size()-1);

    }
  }
  /* Close the input file */
  file.close();

  /* Deallocate dinamic memory */
  if (line != NULL)
    delete [] line;

  /* Check the matrix's content */
  return fillMatrices(true);
}

void alignment::alignmentPhylipToFile(ostream &file) {
  /* Generate output alignment in PHYLIP/PHYLIP 4 format (sequential) */

  int i, j, maxLongName;
  string *tmpMatrix;

  /* Check whether sequences in the alignment are aligned or not.
   * Warn about it if there are not aligned. */
  if (!isAligned) {
    cerr << endl << "ERROR: Sequences are not aligned. Format (PHYLIP) "
      << "not compatible with unaligned sequences." << endl << endl;
    return ;
  }

  /* Allocate local memory for generating output alignment */
  tmpMatrix = new string[sequenNumber];

  /* Depending on alignment orientation: forward or reverse. Copy directly
   * sequence information or get firstly the reversed sequences and then
   * copy it into local memory */
  for(i = 0; i < sequenNumber; i++)
    tmpMatrix[i] = (!reverse) ? sequences[i] : utils::getReverse(sequences[i]);

  /* Depending on if short name flag is activated (limits sequence name up to
   * 10 characters) or not, get maximum sequence name length */
  maxLongName = PHYLIPDISTANCE;
  for(i = 0; (i < sequenNumber) && (!shortNames); i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());

  /* Generating output alignment */
  /* First Line: Sequences Number & Residued Number */
  file << " " << sequenNumber << " " << residNumber << endl;

  /* First Block: Sequences Names & First 60 residues */
  for(i = 0; i < sequenNumber; i++)
    file << setw(maxLongName + 3) << left << seqsName[i].substr(0, maxLongName)
      << tmpMatrix[i].substr(0, 60) << endl;
  file << endl;

  /* Rest of blocks: Print 60 residues per each blocks of sequences */
  for(i = 60; i < residNumber; i += 60) {
    for(j = 0; j < sequenNumber; j++)
      file << tmpMatrix[j].substr(i, 60) << endl;
    file << endl;
  }
  file << endl;

  /* Deallocate local memory */
  delete [] tmpMatrix;
}

void alignment::alignmentPhylip3_2ToFile(ostream &file) {
  /* Generate output alignment in PHYLIP 3.2 format (interleaved) */

  int i, j, k, maxLongName;
  string *tmpMatrix;

  /* Check whether sequences in the alignment are aligned or not.
   * Warn about it if there are not aligned. */
  if (!isAligned) {
    cerr << endl << "ERROR: Sequences are not aligned. Format (PHYLIP) "
      << "not compatible with unaligned sequences." << endl << endl;
    return ;
  }

  /* Allocate local memory for generating output alignment */
  tmpMatrix = new string[sequenNumber];

  /* Depending on alignment orientation: forward or reverse. Copy directly
   * sequence information or get firstly the reversed sequences and then
   * copy it into local memory */
  for(i = 0; i < sequenNumber; i++)
    tmpMatrix[i] = (!reverse) ? sequences[i] : utils::getReverse(sequences[i]);

  /* Depending on if short name flag is activated (limits sequence name up to
   * 10 characters) or not, get maximum sequence name length */
  maxLongName = PHYLIPDISTANCE;
  for(i = 0; (i < sequenNumber) && (!shortNames); i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());

  /* Generating output alignment */
  /* First Line: Sequences Number & Residued Number */
  file << " " << sequenNumber << " " << residNumber << endl;

  /* Alignment */
  /* For each sequence, print its identifier and then the sequence itself in
   * blocks of 50 residues */
  for(i = 0; i < sequenNumber; i++) {
    /* Sequence Name */
    file << setw(maxLongName + 3) << left << seqsName[i].substr(0, maxLongName);
    /* Sequence. Each line contains a block of 5 times 10 residues. */
    for(j = 0; j < residNumber; j += 50) {
      for(k = j; (k < residNumber) && (k < (j + 50)); k += 10)
        file << sequences[i].substr(k, 10) << " ";
      file << endl;
      /* If the sequences end has not been reached, print black spaces
       * to follow format specifications */
      if((j + 50) < residNumber)
        file << setw(maxLongName + 3) << " ";
    }
    /* Print a blank line to mark sequences separation */
    file << endl;
  }

  /* Deallocate local memory */
  delete [] tmpMatrix;
}

void alignment::alignmentPhylip_PamlToFile(ostream &file) {
  /* Generate output alignment in PHYLIP format compatible with PAML program */

  int i, maxLongName;
  string *tmpMatrix;

  /* Check whether sequences in the alignment are aligned or not.
   * Warn about it if there are not aligned. */
  if (!isAligned) {
    cerr << endl << "ERROR: Sequences are not aligned. Format (PHYLIP) "
      << "not compatible with unaligned sequences." << endl << endl;
    return ;
  }

  /* Allocate local memory for generating output alignment */
  tmpMatrix = new string[sequenNumber];

  /* Depending on alignment orientation: forward or reverse. Copy directly
   * sequence information or get firstly the reversed sequences and then
   * copy it into local memory */
  for(i = 0; i < sequenNumber; i++)
    tmpMatrix[i] = (!reverse) ? sequences[i] : utils::getReverse(sequences[i]);

  /* Depending on if short name flag is activated (limits sequence name up to
   * 10 characters) or not, get maximum sequence name length */
  maxLongName = PHYLIPDISTANCE;
  for(i = 0; (i < sequenNumber) && (!shortNames); i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());

  /* Generating output alignment */
  /* First Line: Sequences Number & Residued Number */
  file << " " << sequenNumber << " " << residNumber << endl;

  /* Print alignment */
  /* Print sequences name follow by the sequence itself in the same line */
  for(i = 0; i < sequenNumber; i++)
    file << setw(maxLongName + 3) << left << seqsName[i].substr(0, maxLongName)
      << sequences[i] << endl;
  file << endl;

  /* Deallocate local memory */
  delete [] tmpMatrix;
}

void alignment::alignmentClustalToFile(ostream &file) {
  /* Generate output alignment in CLUSTAL format */

  int i, j, maxLongName = 0;
  string *tmpMatrix;

  /* Check whether sequences in the alignment are aligned or not.
   * Warn about it if there are not aligned. */
  if (!isAligned) {
    cerr << endl << "ERROR: Sequences are not aligned. Format (CLUSTAL) "
      << "not compatible with unaligned sequences." << endl << endl;
    return ;
  }

  /* Allocate local memory for generating output alignment */
  tmpMatrix = new string[sequenNumber];

  /* Depending on alignment orientation: forward or reverse. Copy directly
   * sequence information or get firstly the reversed sequences and then
   * copy it into local memory */
  for(i = 0; i < sequenNumber; i++)
    tmpMatrix[i] = (!reverse) ? sequences[i] : utils::getReverse(sequences[i]);

  /* Compute maximum sequences name length */
  for(i = 0; (i < sequenNumber) && (!shortNames); i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());

  /* Print alignment header */
  if((aligInfo.size() != 0)  && (iformat == oformat))
    file << aligInfo << endl << endl;
  else
    file << "CLUSTAL multiple sequence alignment" << endl << endl;

  /* Print alignment itself */
  /* Print as many blocks as it is needed of lines composed
   * by sequences name and 60 residues */
  for(j = 0; j < residNumber; j += 60) {
    for(i = 0; i < sequenNumber; i++)
      file << setw(maxLongName + 5) << left << seqsName[i]
        << tmpMatrix[i].substr(j, 60) << endl;
    file << endl << endl;
  }

  /* Deallocate local memory */
  delete [] tmpMatrix;
}

void alignment::alignmentFastaToFile(ostream &file) {
  /* Generate output alignment in FASTA format. Sequences can be unaligned. */

  int i, j, maxLongName;
  string *tmpMatrix;

  /* Allocate local memory for generating output alignment */
  tmpMatrix = new string[sequenNumber];

  /* Depending on alignment orientation: forward or reverse. Copy directly
   * sequence information or get firstly the reversed sequences and then
   * copy it into local memory */
  for(i = 0; i < sequenNumber; i++)
    tmpMatrix[i] = (!reverse) ? sequences[i] : utils::getReverse(sequences[i]);

  /* Depending on if short name flag is activated (limits sequence name up to
   * 10 characters) or not, get maximum sequence name length. Consider those
   * cases when the user has asked to keep original sequence header */
  maxLongName = 0;
  for(i = 0; i < sequenNumber; i++)
    if (!keepHeader)
      maxLongName = utils::max(maxLongName, seqsName[i].size());
    else if (seqsInfo != NULL)
      maxLongName = utils::max(maxLongName, seqsInfo[i].size());

   if (shortNames && maxLongName > PHYLIPDISTANCE) {
    maxLongName = PHYLIPDISTANCE;
    if (keepHeader)
      cerr << endl << "WARNING: Original sequence header will be cut by charac"
        << "ter 10" << endl;
    }
  /* Print alignment. First, sequences name id and then the sequences itself */
  for(i = 0; i < sequenNumber; i++) {
    if (!keepHeader)
      file << ">" << seqsName[i].substr(0, maxLongName) << endl;
    else if (seqsInfo != NULL)
      file << ">" << seqsInfo[i].substr(0, maxLongName) << endl;
    for(j = 0; j < residuesNumber[i]; j+= 60)
      file << tmpMatrix[i].substr(j, 60) << endl;
  }

  /* Deallocate local memory */
  delete [] tmpMatrix;
}

void alignment::alignmentNexusToFile(ostream &file) {
  /* Generate output alignment in NEXUS format setting only alignment block */

  int i, j, k, maxLongName = 0;
  string *tmpMatrix;

  /* Check whether sequences in the alignment are aligned or not.
   * Warn about it if there are not aligned. */
  if (!isAligned) {
    cerr << endl << "ERROR: Sequences are not aligned. Format (NEXUS) "
      << "not compatible with unaligned sequences." << endl << endl;
    return ;
  }

  /* Allocate local memory for generating output alignment */
  tmpMatrix = new string[sequenNumber];

  /* Depending on alignment orientation: forward or reverse. Copy directly
   * sequence information or get firstly the reversed sequences and then
   * copy it into local memory */
  for(i = 0; i < sequenNumber; i++)
    tmpMatrix[i] = (!reverse) ? sequences[i] : utils::getReverse(sequences[i]);

  /* Compute maximum sequences name length */
  for(i = 0; (i < sequenNumber) && (!shortNames); i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());

  /* Compute output file datatype */
  getTypeAlignment();

  /* Remove characters like ";" from input alignment information line */
  while((int) aligInfo.find(";") != (int) string::npos)
    aligInfo.erase(aligInfo.find(";"), 1);

  /* Print Alignment header */
  file << "#NEXUS" << endl << "BEGIN DATA;" << endl << " DIMENSIONS NTAX="
    << sequenNumber << " NCHAR=" << residNumber <<";" << endl;

  /* Print alignment datatype */
  if ((dataType == DNAType) || (dataType == DNADeg))
    file << "FORMAT DATATYPE=DNA INTERLEAVE=yes GAP=-";
  else if ((dataType == RNAType) || (dataType == RNADeg))
    file << "FORMAT DATATYPE=RNA INTERLEAVE=yes GAP=-";
  else if (dataType == AAType)
    file << "FORMAT DATATYPE=PROTEIN INTERLEAVE=yes GAP=-";

  i = 0;
  /* Using information from input alignment. Use only some tags. */
  while((j = aligInfo.find(" ", i)) != (int) string::npos) {

    if((aligInfo.substr(i, j - i)).compare(0, 7, "MISSING") == 0 ||
       (aligInfo.substr(i, j)).compare(0, 7, "missing") == 0)
      file << " " << (aligInfo.substr(i, j - i));

    else if((aligInfo.substr(i, j)).compare(0, 9, "MATCHCHAR") == 0 ||
       (aligInfo.substr(i, j)).compare(0, 9, "matchchar") == 0)
      file << " " << (aligInfo.substr(i, j - i));

    i = j + 1;
  }
  file << ";" << endl;

  /* Print sequence name and sequence length */
  for(i = 0; i < sequenNumber; i++)
    file << "[Name: " << setw(maxLongName + 4) << left << seqsName[i] << "Len: "
      << residNumber << "]" << endl;
  file << endl << "MATRIX" << endl;

  /* Print alignment itself. Sequence name and 50 residues blocks */
  for(j = 0; j < residNumber; j += 50) {
    for(i = 0; i < sequenNumber; i++) {
      file << setw(maxLongName + 4) << left << seqsName[i];
      for(k = j; k < (j + 50) && k < residNumber; k += 10)
        file << " " << sequences[i].substr(k, 10);
      file << endl;
    }
    file << endl;
  }
  file << ";" << endl << "END;" << endl;

  /* Deallocate local memory */
  delete [] tmpMatrix;
}

void alignment::alignmentMegaToFile(ostream &file) {
  /* Generate output alignment in MEGA format */

  int i, j, k;
  string *tmpMatrix;

  /* Check whether sequences in the alignment are aligned or not.
   * Warn about it if there are not aligned. */
  if (!isAligned) {
    cerr << endl << "ERROR: Sequences are not aligned. Format (MEGA) "
      << "not compatible with unaligned sequences." << endl << endl;
    return ;
  }

  /* Allocate local memory for generating output alignment */
  tmpMatrix = new string[sequenNumber];

  /* Depending on alignment orientation: forward or reverse. Copy directly
   * sequence information or get firstly the reversed sequences and then
   * copy it into local memory */
  for(i = 0; i < sequenNumber; i++)
    tmpMatrix[i] = (!reverse) ? sequences[i] : utils::getReverse(sequences[i]);

  /* Compute output file datatype */
  getTypeAlignment();

  /* Print output alignment header */
  file << "#MEGA" << endl << filename << endl;

  /* Print alignment datatype */
  if ((dataType == DNAType) || (dataType == DNADeg))
    file << "!Format DataType=DNA ";
  else if ((dataType == RNAType) || (dataType == RNADeg))
    file << "!Format DataType=RNA ";
  else if (dataType == AAType)
    file << "!Format DataType=protein ";

  /* Print number of sequences and alignment length */
  file << "NSeqs=" << sequenNumber << " Nsites=" << residNumber
    << " indel=- CodeTable=Standard;" << endl << endl;

  /* Print sequences name and sequences divided into blocks of 50 residues */
  for(i = 0; i < sequenNumber; i++) {
    file << "#" << seqsName[i] << endl;
    for(j = 0; j < residNumber; j += 50) {
      for(k = j; ((k < residNumber) && (k < j + 50)); k += 10)
        file << tmpMatrix[i].substr(k, 10) << " ";
      file << endl;
    }
    file << endl;
  }

  /* Deallocate local memory */
  delete [] tmpMatrix;
}

void alignment::alignmentNBRF_PirToFile(ostream &file) {
  /* Generate output alignment in NBRF/PIR format. Sequences can be unaligned */

  int i, j, k;
  string alg_datatype, *tmpMatrix;

  /* Allocate local memory for generating output alignment */
  tmpMatrix = new string[sequenNumber];

  /* Depending on alignment orientation: forward or reverse. Copy directly
   * sequence information or get firstly the reversed sequences and then
   * copy it into local memory */
  for(i = 0; i < sequenNumber; i++)
    tmpMatrix[i] = (!reverse) ? sequences[i] : utils::getReverse(sequences[i]);

  /* Compute output file datatype */
  getTypeAlignment();
  if ((dataType == DNAType) || (dataType == DNADeg))
    alg_datatype = "DL";
  else if ((dataType == RNAType) || (dataType == RNADeg))
    alg_datatype = "RL";
  else if (dataType == AAType)
    alg_datatype = "P1";

  /* Print alignment */
  for(i = 0; i < sequenNumber; i++) {

    /* Print sequence datatype and its name */
    if((seqsInfo != NULL) && (iformat == oformat))
      file << ">" << seqsInfo[i].substr(0, 2) << ";" << seqsName[i]
        << endl << seqsInfo[i].substr(2) << endl;
    else
      file << ">" << alg_datatype << ";" << seqsName[i] << endl
        << seqsName[i] << " " << residuesNumber[i] << " bases" << endl;

    /* Write the sequence */
    for(j = 0; j < residuesNumber[i]; j += 50) {
      for(k = j; (k < residuesNumber[i]) && (k < (j + 50)); k += 10)
        file << " " << tmpMatrix[i].substr(k, 10);

      if(k >= residuesNumber[i]) {
        if((residuesNumber[i] % 50) == 0)
          file << endl << " ";
        else if((residuesNumber[i] % 10) == 0)
          file << " ";
        file << "*";
      }
      file << endl;
    }
    file << endl;
  }

  /* Deallocate local memory */
  delete [] tmpMatrix;
}

bool alignment::alignmentSummaryHTML(char *destFile, int residues, int seqs, \
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
  tmpColumn.reserve(sequenNumber);

  /* Check whether sequences in the alignment are aligned or not.
   * Warn about it if there are not aligned. */
  if (!isAligned) {
    cerr << endl << "ERROR: Sequences are not aligned." << endl << endl;
    return false;
  }

  /* Open output file and check that file pointer is valid */
  file.open(destFile);
  if(!file)
    return false;

  /* Compute maximum sequences name length. */
  maxLongName = 0;
  for(i = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());

  /* Compute HTML blank spaces */
  minHTML = utils::max(25, maxLongName + 10);

  /* Initialize local variables to control which columns/sequences
   * will be kept in the output alignment */
  res = new bool[residNumber];
  for(i = 0; i < residNumber; i++)
    res[i] = false;

  seq = new bool[sequenNumber];
  for(i = 0; i < sequenNumber; i++)
    seq[i] = false;

  /* Record which columns/sequences from original alignment
   * have been kept in the final one */
  for(i = 0; i < residues; i++)
    res[selectedRes[i]] = true;
  for(i = 0; i < seqs; i++)
    seq[selectedSeq[i]] = true;

  /* Recover some stats about different scores from current alignment */
  gapsValues = NULL;
  if (sgaps != NULL)
    gapsValues = sgaps -> getGapsWindow();
  simValues = NULL;
  if (scons != NULL)
    simValues = scons -> getMdkwVector();

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
    << sequenNumber - seqs << " /Deleted Residues:  " << setw(7) << right
    << residNumber - residues << "</span>" << endl;

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
  for(j = 0, upper = HTMLBLOCKS; j < residNumber; j += HTMLBLOCKS, upper += \
    HTMLBLOCKS) {

    /* Print main columns number */
    file << endl << setw(minHTML + 10) << right << (j + 10);
    for(i = j + 20; ((i <= residNumber) && (i <= upper)); i += 10)
      file << setw(10) << right << (i);

    /* Print special characters to delimit sequences blocks */
    file << endl << setw(minHTML + 1) << right;
    for(i = j + 1; ((i <= residNumber) && (i <= upper)); i++)
      file << (!(i % 10) ? "+" : "=");
    file << endl;

    /* Print sequences name */
    for(i = 0; i < sequenNumber; i++) {
      file << "    <span class=" << ((seq[i]) ? "sel>" : "nsel>") << seqsName[i]
        << "</span>" << setw(minHTML - 4 - seqsName[i].size()) << right << "";

      /* Print residues corresponding to current sequences block */
      for(k = j; ((k < residNumber) && (k < upper)); k++) {
        for(kj = 0, tmpColumn.clear(); kj < sequenNumber; kj++)
          tmpColumn += sequences[kj][k];
        /* Determine residue color based on residues across the alig column */
        type = utils::determineColor(sequences[i][k], tmpColumn);
        if (type == 'w')
          file << sequences[i][k];
        else
          file << "<span id=" << type << ">" << sequences[i][k] << "</span>";
      }
      file << endl;
    }

    file << endl << setw(minHTML) << left << "    Selected Cols:      ";
    for(k = j; ((k < residNumber) && (k < (j + HTMLBLOCKS))); k++)
      file << "<span class=" << (res[k] ? "sel" : "nsel") << "> </span>";
    file << endl;

    /* If there is not any score to print, skip this part of the function */
    if ((gapsValues == NULL) and (simValues == NULL) and (consValues == NULL))
      continue;

    /* Print score colors according to certain predefined thresholds */
    if (gapsValues != NULL) {
      file << endl << setw(minHTML) << left << "    Gaps Scores:        ";
      for(k = j; ((k < residNumber) && (k < (j + HTMLBLOCKS))); k++)
        if(gapsValues[k] == 0)
          file << "<span class=c12> </span>";
        else if(gapsValues[k] == sequenNumber)
          file << "<span class=c1> </span>";
        else if(1 - (float(gapsValues[k])/sequenNumber) >= .750)
          file << "<span class=c11> </span>";
        else if(1 - (float(gapsValues[k])/sequenNumber) >= .500)
          file << "<span class=c10> </span>";
        else if(1 - (float(gapsValues[k])/sequenNumber) >= .350)
          file << "<span  class=c9> </span>";
        else if(1 - (float(gapsValues[k])/sequenNumber) >= .250)
          file << "<span  class=c8> </span>";
        else if(1 - (float(gapsValues[k])/sequenNumber) >= .200)
          file << "<span  class=c7> </span>";
        else if(1 - (float(gapsValues[k])/sequenNumber) >= .150)
          file << "<span  class=c6> </span>";
        else if(1 - (float(gapsValues[k])/sequenNumber) >= .100)
          file << "<span  class=c5> </span>";
        else if(1 - (float(gapsValues[k])/sequenNumber) >= .050)
          file << "<span  class=c4> </span>";
        else if(1 - (float(gapsValues[k])/sequenNumber) >= .001)
          file << "<span  class=c3> </span>";
        else
          file << "<span  class=c2> </span>";
    }
    if (simValues != NULL) {
      file << endl << setw(minHTML) << left << "    Similarity Scores:  ";
      for(k = j; ((k < residNumber) && (k < (j + HTMLBLOCKS))); k++)
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
      for(k = j; ((k < residNumber) && (k < (j + HTMLBLOCKS))); k++)
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

bool alignment::alignmentColourHTML(ostream &file) {

  int i, j, kj, upper, k = 0, maxLongName = 0;
  string tmpColumn;
  char type;

  /* Allocate some local memory */
  tmpColumn.reserve(sequenNumber);

  /* Check whether sequences in the alignment are aligned or not.
   * Warn about it if there are not aligned. */
  if (!isAligned) {
    cerr << endl << "ERROR: Sequences are not aligned." << endl << endl;
    return false;
  }

  /* Compute maximum sequences name length */
  maxLongName = 0;
  for(i = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, seqsName[i].size());


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
  for(j = 0, upper = HTMLBLOCKS; j < residNumber; j += HTMLBLOCKS, upper += \
    HTMLBLOCKS) {

    file << endl;
    /* Print main columns number */
    file << setw(maxLongName + 19) << right << (j + 10);
    for(i = j + 20; ((i <= residNumber) && (i <= upper)); i += 10)
      file << setw(10) << right << i;

    /* Print special characters to delimit sequences blocks */
    file << endl << setw(maxLongName + 10);
    for(i = j + 1; ((i <= residNumber) && (i <= upper)); i++)
      file << (!(i % 10) ? "+" : "=");

    /* Print sequences themselves */
    for(i = 0; i < sequenNumber; i++) {

      /* Print sequences name */
      file << endl << setw(maxLongName + 9) << left << seqsName[i];

      /* Print residues corresponding to current sequences block */
      for(k = j; ((k < residNumber) && (k < upper)); k++) {
        for(kj = 0, tmpColumn.clear(); kj < sequenNumber; kj++)
          tmpColumn += sequences[kj][k];
        /* Determine residue color based on residues across the alig column */
        type = utils::determineColor(sequences[i][k], tmpColumn);
        if (type == 'w')
          file << sequences[i][k];
        else
          file << "<span id=" << type << ">" << sequences[i][k] << "</span>";
      }
    }
    file << endl;
  }

  /* Print HTML footer into output file */
  file << "    </pre>" << endl << "  </body>" << endl << "</html>" << endl;

  return true;
}

void alignment::printAlignmentInfo(ostream &file) {
  /* Print information about sequences number, average sequence length, maximum
   * and minimum sequences length, etc */

  int i, j, valid_res, max, min, max_pos, min_pos, total_res;

  /* Storage which sequences are the longest and shortest ones */
  max = 0;
  max_pos = 0;
  min_pos = 0;
  min = residuesNumber[0];

  for(i = 0, total_res = 0; i < sequenNumber; i++) {

    /* Discard gaps from current sequence and then compute real length */
    for(j = 0, valid_res = 0; j < residuesNumber[i]; j++)
      valid_res += (sequences[i][j] != '-' ? 1 : 0);

    /* Compute the total residues in the alignment to calculate avg. sequence
     * length */
    total_res += valid_res;

    /* Get values for the longest sequence */
    max_pos = (max > valid_res) ? max_pos : i;
    max = (max > valid_res) ? max : valid_res;
    /* Similarily, get values for the shortest sequence */
    min_pos = (min < valid_res) ? min_pos : i;
    min = (min < valid_res) ? min : valid_res;
  }

  file << "## Total sequences\t" << sequenNumber << endl;
  if (isFileAligned())
    file << "## Alignment length\t" << residNumber << endl;
  file  << "## Avg. sequence length\t" << (float) total_res / sequenNumber << endl
    << "## Longest seq. name\t'"   << seqsName[max_pos] << "'" << endl
    << "## Longest seq. length\t"  << max << endl
    << "## Shortest seq. name\t'"  << seqsName[min_pos] << "'" << endl
    << "## Shortest seq. length\t" << min << endl;
}
