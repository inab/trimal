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

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <string>
#include <cmath>
#include <map>
#include <sys/stat.h>

/**
 \brief Utilities class.
 This class contains shared methods to be used in multiple parts of the code.
 */
namespace utils {

// public:
    /**
     \brief Vector initialization.
     \param [out] vector The vector that will be initialized.
     \param tam The size of the vector.
     \param valor The initialization value that will of all positions of the vector.

     This method is used to initialize
     all positions of a vector with a given value.
     */
     void initlVect(int *vector, int tam, int valor);

    /**
     \brief Vector initialization.
     \param [out] vector The vector that will be initialized.
     \param tam The size of the vector.
     \param valor The initialization value of all positions of the vector.

     This method is used to initialize
     all positions of a vector with a given value.
     */
     void initlVect(float *vector, int tam, float valor);

    /**
     \brief Integer vector copying.
     \param vect1 Vector that we want to copy.
     \param [out] vect2 Destination vector of the copy.
     \param tam Vectors size.

     This method copies integer vector 1 to integer vector 2.
      */
     void copyVect(int *vect1, int *vect2, int tam);

    /**
     \brief Float vector copying.
     \param vect1 Vector that we want to copy.
     \param [out] vect2 Destination vector of the copy.
     \param tam Vectors size.

     This method copies float vector 1 to float vector 2.
     */
     void copyVect(float *vect1, float *vect2, int tam);

    /**
     \brief Round double to inferior integer method.
     \param number The number that will be rounded.
     \return the rounded number.

     This method rounds a double number to the inferior integer.
     */
     int roundToInf(double number);

    /**
     \brief Round double to integer method.
     \param number The number that will be rounded.
     \return the rounded number.

     This method rounds a double number to a integer.
     */
     int roundInt(double number);

    /**
     \brief Round double to greater integer method.
     \param number The number that will be rounded.
     \return the rounded number.

     This method rounds a double number to the greater integer.
     */
     int roundToSup(double number);

    /**
     \brief Maximum of two numbers method.
     \param x The first number.
     \param y The second number.
     \return The maximum between the two given numbers.

     This method returns the maximum between the two numbers given as parameters.
     */
     int max(int x, int y);

    /**
     \brief Maximum of two numbers method.
     \param x The first number.
     \param y The second number.
     \return The maximum between the two given numbers.

     This method returns the maximum between the two numbers given as parameters.
     */
     float max(float x, float y);

    /**
     \brief Maximum of two numbers method.
     \param x The first number.
     \param y The second number.
     \return The maximum between the two given numbers.

     This method returns the maximum between the two numbers given as parameters.
     */
     double max(double x, double y);

    /**
     \brief Minimum of two numbers method.
     \param x The first number.
     \param y The second number.
     \return The minumum between the two given numbers.

     This method returns the minimum between the two numbers given as parameters.
     */
     int min(int x, int y);

    /**
     \brief Minimum of two numbers method.
     \param x The first number.
     \param y The second number.
     \return The minumum between the two given numbers.

     This method returns the minimum between the two numbers given as parameters.
     */
     float min(float x, float y);

    /**
     \brief Minimum of two numbers method.
     \param x The first number.
     \param y The second number.
     \return The minumum between the two given numbers.

     This method returns the minimum between the two numbers given as parameters.
     */
     double min(double x, double y);

    /**
     \brief String-is-number checking.
     \param num The string we want to check.
     \return \b true if the string is a number, \b false if not.

     This method checks if the given string is a number
     (taking in mind the possibility of floating numbers and scientific notation)
     */
     bool isNumber(char *num);

    /**
     \brief String comparing method.
     \param a The first string that will be compared.
     \param b The second string that will be compared.
     \return \b true if the two strings are the same, \b false if not.

     This method compares the two strings given,
     and returns \b true if the two strings are equal.
     */
     bool compare(char *a, char *b);

    /**
     \brief Removing spaces method.
     \param in The string that we want to clean.
     \param[out] out The destination of the clean string.

     This method removes spaces in the input string
     and put the result in the output string.
     */
     void removeSpaces(char *in, char *out);

    /**
     \brief Quicksort sorting method.
     \param list The vector that we want to sort.
     \param ini The first element of the vector.
     \param fin The last element of the vector.

     This method sorts the vector using the quicksort method.
     */
     void quicksort(float *list, int ini, int fin);

    /**
     \brief Swapping elements method
     \param a One element to swap.
     \param b Other element to swap.

     This method swaps the values in a and b.
     */
     void swap(float *a, float *b);

    /**
     \brief Quicksort sorting method.
     \param list The vector that we want to sort.
     \param ini The first element of the vector.
     \param fin The last element of the vector.

     This method sorts the vector using the quicksort method.
     */
     void quicksort(int *list, int ini, int fin);

    /**
     \brief Swapping elements method
     \param a One element to swap.
     \param b Other element to swap.

     This method swaps the values in a and b.
     */
     void swap(int *a, int *b);

    /**
     \brief Check if a given file exists and its size is greater than 0.
     \param file ifstream to check
     */
     bool checkFile(std::ifstream &file);

    /**
     \brief Read a new line from current input stream.\n
     This function is better than standard one
      since cares of operative system compatibility.\n
     It is useful as well because removes tabs and blank spaces
      at lines at beginning/ending.\n
     \param file ifstream to read line from.
     \return \n
        Line that has been read or
        nullptr if there is nothing to read.\n
     */
     char *readLine(std::ifstream &file);

    /**
     \brief Read a new line from current input stream.\n
     This function is better than standard one
     since cares of operative system compatibility.\n
     It is useful as well because removes tabs
     and blank spaces at lines at beginning/ending.\n
     \param file ifstream to read line from.
     \return \n
      nullptr if there is nothing to read.\n
      Line that has been read.
    */
     char *readLine(std::istream &file);

    /**
     \brief Remove all content surrounded by ("") or ([]).\n
            It warns as well when a mismatch for these flags is found. \n
     \param nline Line to be trimmed.
     \return
            Line trimmed of comments or
            nullptr if there has been a mismatch\n
     */
     char *trimLine(std::string nline);

    /**
     \brief Reverses a string
     \param toReverse String to get a reversed copy.
     \return Reversed string of toReverse.
     */
     std::string getReverse(const std::string &toReverse);

    /**
     \brief Removes a determined char from the string
     \param c Character to remove from line
     \param line String to remove c from.
     \return New string without c character
     */
     std::string removeCharacter(char c, std::string line);

    /**
     \brief Checks an alignment type
     \param seqNumber Number of sequences to check it's type.
     \param residNumber Number of residues of the alignment.
     \param sequences Sequences pointer
     \return Integer that represents the alignment type.
     */
     int checkAlignmentType(int seqNumber,
                                  int residNumber,
                                  std::string *sequences);

    /**
     \brief Reads a line and converts it to an array of number
     \param line Line to convert to array of ints
     \return Pointer to an array of numbers that contains line
     */
     int *readNumbers(const std::string &line);

    /**
     \brief Quicksort sorting method.
     \param vect The vector that we want to sort.
     \param ini The first element of the vector.
     \param fin The last element of the vector.
    */
     void quicksort(int **vect, int ini, int fin);

    /**
     \brief Swaps double pointers.
     \param a Double pointer A
     \param b Double pointer B
    */
     void swap(int **a, int **b);

    /**
     \brief Checks the color that has to be used on the output report
     \param res Resiude to check its color
     \param column Column to which this residue belongs.
     \return Char that represents the color to be used.
    */
     char determineColor(char res, const std::string& column);

    /**
     \brief Method to check for a pattern in a string.\n
            The method will check, character by character of the first string if
             there is some equality for each character in the pattern.\n
            When done, it will calculate the fraction of characters
             present in the pattern and compare to the threshold argument.
     \param data string that will be compared against a pattern
     \param pattern string that contains the pattern.
     \param threshold minimum ratio of hits to consider the pattern valid
     */
     bool lookForPattern(const std::string& data,
                               const std::string& pattern,
                               const float threshold);

    /**
     \brief Function that replaces a substring
              with another substring in a string.
            It does not make a copy of the original string, but modifies it.
     \param [in,out] subject String to be modified
     \param search Substring to search and change
     \param replace Substring to put in place of search
     */
     void ReplaceStringInPlace(std::string &subject,
                                     const std::string &search,
                                     const std::string &replace);

    /**
     \brief Function that replaces a substring
             with another substring in a string.
            It makes a copy of the original string.
     \param [in] subject String to be modified
     \param search Substring to search and change
     \param replace Substring to put in place of search
     */
     std::string ReplaceString(std::string subject,
                                     const std::string &search,
                                     const std::string &replace);

    /**
     \brief Function that gives the gap classification of a column of values.
     \param gapValue Number of gaps present in the column.
     \param sequenNumber Number of sequences.
     \return Int representing the classification of this gap value.
     */
     int GetGapStep(int *gapValue, int sequenNumber);

    /**
     \brief Function that gives the gap classification of a column of values.\n
     This function should work faster than it's sister
       utils::GetGapStep(int * gapValue, int sequenNumber),
       as it uses a precomputed (by the user)
       inverseSequenNumber (1F / Alignment::sequenNumber),
       instead of calculating it over again each time the function
       is called (which is equal to number of residues). \n
     This comes with a precision cost that shouldn't be a problem.
     \param gapValue Number of gaps present in the column.
     \param inverseSequenNumber Inverse of number of sequences. (1F / sequenNumber)
     \return Int representing the classification of this gap value.
     */
     int GetGapStep(int *gapValue, float inverseSequenNumber);

    /**
     \brief Function that gives the
      similarity classification of a column of values.
     \param simValue Similarity value.
     \return Int representing the classification of this gap value.
     */
     int GetSimStep(float *simValue);

    /**
     \brief Function that gives the
      consistency classification of a column of values.
     \param consValue Consistency value.
     \return Int representing the classification of this gap value.
     */
     int GetConsStep(float *consValue);

     /**
      * \brief Method to check the existance of a file
      * @param path Path to the file to check
      * @return Wheter the file exists or not.
      */
     bool fileExists(std::string & path);

    /**
     * \brief Method to check the existance of a file.
     * Works exactly as fileExists(std::string & path),
     * but accepts r-value reference.
     * @param path Path to the file to check
     * @return Wheter the file exists or not.
     */
    bool fileExists(std::string && path);

    /***
     *  Method to transform a char to its upper version
     *  Will return the same char if its already and uppercase char
     *  Works using bit shifts, to avoid using locale.
     *
     *  If char is not a alpha-character, it will return the same char.
     * @param c Original character
     * @return upperCase version of the character
     */
     char toUpper(char c);

     namespace TerminalColors {
         enum terminalColor
         { RESET, BLACK, RED,
           GREEN, YELLOW, BLUE,
           MAGENTA, CYAN, WHITE,
           BOLD, UNDERLINE };

         extern std::map<terminalColor, const std::string> colors;
     }
};



#endif
