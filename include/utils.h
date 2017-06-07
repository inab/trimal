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

/** \brief Utils class.
 *
 * This class implements util methods.
 */

using namespace std;

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

  static float max(float x, float y);

  static double max(double x, double y);

  static int min(int x, int y);

  static float min(float x, float y);

  static double min(double x, double y);

  // static bool getArg(int argc, char *argv[], int *var, char *argument, char *abrevArg);
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
   * \param out The destination of the clean string.
   *
   * This method removes spaces in the input string and put the result in the output string.
   */
  static void removeSpaces(char *in, char *out);

  /** \brief Quicksort sorting method.
   * \a param list The vector that we want to sort.
   * \a param ini The first element of the vector.
   * \a param fin The last element of the vector.
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
   * \a param list The vector that we want to sort.
   * \a param ini The first element of the vector.
   * \a param fin The last element of the vector.
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

  static bool checkFile(ifstream &file);

  static char* readLine(ifstream &file);

  static char* trimLine(string nline);

  static char* readLineMEGA(ifstream &file);

  static string getReverse(string toReverse);

  static string removeCharacter(char c, string line);

  static int checkTypeAlignment(int, int, string *);

  static int* readNumbers(string);

  static void quicksort(int **, int, int);

  static void swap(int **, int **);

  static char determineColor(char res, string column);

  static bool lookForPattern(string, string, float);

};
#endif
