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

#ifndef UTILS_H
#define UTILS_H

#include <string.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>

#define DNAType 1
#define RNAType 2
#define AAType  3

#define DNADeg  4
#define RNADeg  5



using namespace std;

/** \brief Utils class.
 *
 * This class implements util methods.
 */
class utils {

 public:
  /** \brief Vector initialization.
   * \param vector The vector that will be initializated.
   * \param tam The size of the vector.
   * \param valor The initialization value of all positions of the vector.
   *
   * This method is used to initializate all positions of a vector with a given value.
   */
  static void initlVect(int *vector, int tam, int valor);
    /** \brief Vector initialization.
   * \param vector The vector that will be initializated.
   * \param tam The size of the vector.
   * \param valor The initialization value of all positions of the vector.
   *
   * This method is used to initializate all positions of a vector with a given value.
   */
  static void initlVect(float *vector, int tam, float valor);

  /** \brief Integer vector copying.
   * \param vect1 Vector that we want to copy.
   * \param vect2 Destination vector of the copy.
   * \param tam Vectors size.
   *
   * This method copies integer vector 1 to integer vector 2.
   */
  static void copyVect(int *vect1, int *vect2, int tam);

  /** \brief Float vector copying.
   * \param vect1 Vector that we want to copy.
   * \param vect2 Destination vector of the copy.
   * \param tam Vectors size.
   *
   * This method copies float vector 1 to float vector 2.
   */
  static void copyVect(float *vect1, float *vect2, int tam);

  /** \brief Round double to inferior integer method.
   * \param number The number that will be rounded.
   * \return the rounded number.
   *
   * This method rounds a double number to the inferior integer.
   */
  static int roundToInf(double number);

  /** \brief Round double to integer method.
   * \param number The number that will be rounded.
   * \return the rounded number.
   *
   * This method rounds a double number to a integer.
   */
  static int roundInt(double number);

  /** \brief Round double to greater integer method.
   * \param number The number that will be rounded.
   * \return the rounded number.
   *
   * This method rounds a double number to the greater integer.
   */
  static int roundToSup(double number);

  /** \brief Maximum of two numbers method.
   * \param x The first number.
   * \param y The second number.
   * \return The maximum between the two given numbers.
   *
   * This method returns the maximum between the two numbers given as parameters.
   */
  static int max(int x, int y);
  /** \brief Maximum of two numbers method.
   * \param x The first number.
   * \param y The second number.
   * \return The maximum between the two given numbers.
   *
   * This method returns the maximum between the two numbers given as parameters.
   */
  static float max(float x, float y);
  /** \brief Maximum of two numbers method.
   * \param x The first number.
   * \param y The second number.
   * \return The maximum between the two given numbers.
   *
   * This method returns the maximum between the two numbers given as parameters.
   */
  static double max(double x, double y);
  /** \brief Minimum of two numbers method.
   * \param x The first number.
   * \param y The second number.
   * \return The minumum between the two given numbers.
   *
   * This method returns the minimum between the two numbers given as parameters.
   */
  static int min(int x, int y);
  /** \brief Minimum of two numbers method.
   * \param x The first number.
   * \param y The second number.
   * \return The minumum between the two given numbers.
   *
   * This method returns the minimum between the two numbers given as parameters.
   */
  static float min(float x, float y);
  /** \brief Minimum of two numbers method.
   * \param x The first number.
   * \param y The second number.
   * \return The minumum between the two given numbers.
   *
   * This method returns the minimum between the two numbers given as parameters.
   */
  static double min(double x, double y);

  /** \brief String-is-number checking.
   * \param num The string we want to check.
   * \return \b true if the string is a number, \b false if not.
   *
   * This method checks if the given string is a float number.
   */
  static bool isNumber(char *num);

  /** \brief String comparing method.
   * \param a The first string that will be compared.
   * \param b The second string that will be compared.
   * \return \b true if the two strings are the same, \b false if not.
   *
   * This method compares the two strings given, and returns \b true if the two strings are equal.
   */
  static bool compare(char *a, char *b);

  /** \brief Removing spaces method.
   * \param in The string that we want to clean.
   * \param[out] out The destination of the clean string.
   *
   * This method removes spaces in the input string and put the result in the output string.
   */
  static void removeSpaces(char *in, char *out);

  /** \brief Quicksort sorting method.
   * \param list The vector that we want to sort.
   * \param ini The first element of the vector.
   * \param fin The last element of the vector.
   *
   * This method sorts the vector using the quicksort method.
   */
  static void quicksort(float *list, int ini, int fin);

  /** \brief Swapping elements method
   * \param a One element to swap.
   * \param b Other element to swap.
   *
   * This method swaps the values in a and b.
   */
  static void swap(float *a, float *b);

  /** \brief Quicksort sorting method.
   * \param list The vector that we want to sort.
   * \param ini The first element of the vector.
   * \param fin The last element of the vector.
   *
   * This method sorts the vector using the quicksort method.
   */
  static void quicksort(int *list, int ini, int fin);

  /** \brief Swapping elements method
   * \param a One element to swap.
   * \param b Other element to swap.
   *
   * This method swaps the values in a and b.
   */
  static void swap(int *a, int *b);
  /**
   \brief Check if a given file exists and its size is greater than 0.
   \param file ifstream to check
   */
  static bool checkFile(ifstream &file);
    /**
     \brief Read a new line from current input stream.\n
     This function is better than standard one since cares of operative system compability.\n
     It is useful as well because removes tabs and blank spaces at lines at beginning/ending.\n
     \param file ifstream to read line from.
     \return \n
        NULL if there is nothing to read.\n
        Line that has been read.   
     */
  static char* readLine(ifstream &file);
      /**
     \brief Read a new line from current input stream.\n
     This function is better than standard one since cares of operative system compability.\n
     It is useful as well because removes tabs and blank spaces at lines at beginning/ending.\n
     \param file ifstream to read line from.
     \return \n
        NULL if there is nothing to read.\n
        Line that has been read.   
     */
  static char* readLine(istream &file);
    /**
     \brief Remove comments inside a biological sequence.\n
            Remove all content surrounded by ("") or ([]).\n
            It warns as well when a mismatch for these flags is found. \n
    \param nline Line to be trimmed.
    \return NULL if there has been a mismatch\n
            Line trimmed of comments.
     */
  static char* trimLine(string nline);
    /**
     \todo Implement this function.
     */
  static char* readLineMEGA(ifstream &file);
/**
 \brief Reverses a string
 \param toReverse String to get a reversed copy.
 \return Reversed string of toReverse.
 */
  static string getReverse(string toReverse);
/**
 \brief Removes a determined char from the string
 \param c Character to remove from line
 \param line String to remove c from.
 \return New string without c character
 */
  static string removeCharacter(char c, string line);
/**
 \brief Checks an alignment type 
 \param seqNumber Number of sequences to check it's type.
 \param residNumber Number of residues of the alignment.
 \param sequences Sequences pointer
 \return Integer that represents the alignment type.*/
  static int checkAlignmentType(int seqNumber, int residNumber, string *sequences);
/**
 \brief Reads a line and converts it to an array of number
 \param line Line to convert to array of ints
 \return Pointer to an array of numbers that contains line*/
  static int* readNumbers(string line);
  static int* readNumbers_StartEnd(string line);

/** \brief Quicksort sorting method.
* \param vect The vector that we want to sort.
* \param ini The first element of the vector.
* \param fin The last element of the vector.
*
* This method sorts the vector using the quicksort method.
*/
static void quicksort(int ** vect, int ini, int fin);
/** \brief Swaps pointers values
* \param a Pointer A
* \param b Pointer B
*
* This method swaps the values in a and b.
*/
  static void swap(int ** a , int ** b);
/**
 \brief Checks the color that has to be used on the output report
 \param res Resiude to check its color
 \param column Column to which this residue belongs.
 \return Char that represents the color to be used.
 */
  static char determineColor(char res, string column);
/**
 \brief Looks for a pattern
 \todo Give a good description for this function.*/
  static bool lookForPattern(string, string, float);
/**
 * \brief Function that replaces a substring with another substring in a string.
            It does not make a copy of the original string, but modifies it.
    \param subject String to be modified
    \param search Substring to search and change
    \param replace Substring to put in place of search
 */
  static void ReplaceStringInPlace(std::string& subject, const std::string& search, const std::string& replace);
  /**
 * \brief Function that replaces a substring with another substring in a string.
            It makes a copy of the original string.
    \param subject String to be modified
    \param search Substring to search and change
    \param replace Substring to put in place of search
 */
  static std::string ReplaceString(std::string subject, const std::string& search, const std::string& replace);
};
#endif
