#ifndef MMA_FORMATTER_H_
#define MMA_FORMATTER_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

#include "mmaformatter.h"

const char* BASE_TEN_FORMAT = "*10^";
const char* ARRAY_DELIM_CHARS = ", ";
const char* ARRAY_OUTER_CHARS[] = {"{", "}"};





/* 
 * internal functions
 */




 // public but used internally, so needs pre-declaration here
 char* convertDoubleArrToMMA(double* array, int length, int precision);

void writeNestedDoubleList(
	FILE* file, void* arr, int numDimensions, int* lengths, int lengthInd, int precision
) {
	// base-case: write one list
	if (numDimensions == 1) {
		char* arrStr = convertDoubleArrToMMA((double*) arr, lengths[lengthInd], precision);
		fprintf(file, "%s", arrStr);
		free(arrStr);
		return;
	}
	
	// recursive case: wrap in list (no trailing newline)
	fprintf(file, "{\n");
	for (int i=0; i < lengths[lengthInd]; i++) {
		
		// and write each elem (possibly a pointer) as a list
		writeNestedDoubleList(file, ((double **) arr)[i], numDimensions-1, lengths, lengthInd+1, precision);
		
		// seperate elements
		fprintf(file, ",\n");
	}
	
	// remove trailing comma and newline
	fseek(file, -2, SEEK_END);
	
	// close list
	fprintf(file, "\n}");
}


void writeNestedDoubleArr(
	FILE* file, double* arr, int arrInd, int numDimensions, int* lengths, int lengthInd, int innerTrimLength, int precision
) {
	// base-case: write inner array
	if (lengthInd == numDimensions - 1) {
		char* arrStr = convertDoubleArrToMMA(&(arr[arrInd]), innerTrimLength, precision);
		fprintf(file, "%s", arrStr);
		free(arrStr);
		return;
	}
	
	// work out arrInd jump between outer list elements
	int lengthProd = 1;
	for (int l=lengthInd+1; l < numDimensions; l++)
		lengthProd *= lengths[l];
	
	// recursive case: iterate outer ind, wrapping inner arrays
	fprintf(file, "{\n");
	for (int i=0; i < lengths[lengthInd]; i++) {
		
		// write inner array
		writeNestedDoubleArr(file, arr, arrInd+i*lengthProd, numDimensions, lengths, lengthInd+1, innerTrimLength, precision);
		
		// seperate elements
		fprintf(file, ",\n");
	}
	
	// remove trailing comma and newline
	fseek(file, -2, SEEK_END);
	
	// close array
	fprintf(file, "\n}");
}




/* 
 * public functions
 */




/**
 * returns 'number' as a MMA compatible scientific notation string 
 * with 'precision' digits after the decimal point. The returned string
 * must be freed.
 * 
 * @brief converts double to a scientific notation string
 * @param number		the double to convert
 * @param precision	number of digits after decimal point in sci-not
 * @return 				the sci-not string e.g."-4.321*10^-03", which
 * 						should be freed
 */
char* getScientificNotation(double number, int precision) {
   
   // build sci-not format string
   int precDigits = 1 + ((precision>0)? floor(log10(abs(precision))) : 0);
   int formatStrSize = sizeof(char)*(precDigits + 3 + 1); // % d . e NULL
   char *formatStr = malloc(formatStrSize);
   sprintf(formatStr, "%%.%de", precision); // e.g. %.5e
   
   // build sci-not string with e
   ssize_t rawStrSize = snprintf(NULL, 0, formatStr, number);
   char *rawStr = malloc(rawStrSize + 1); // NULL char
   sprintf(rawStr, formatStr, number); // e.g. 1.23e+03
   
   // find the number of chars in our strings/substrings
   int numBase10Chars = strlen(BASE_TEN_FORMAT);
   int rawStrNumChars = rawStrSize/sizeof(char);
   ssize_t newStrSize = rawStrSize + sizeof(BASE_TEN_FORMAT) - sizeof(char);
   
   // build sci-not string with *10^ by...
   char *string = malloc(newStrSize + 1); // NULL
   memset(string, 0x00, newStrSize+1);
   
   // copying old chars up until e
   int rawInd=0;
   while (rawStr[rawInd] != 'e') {
       string[rawInd] = rawStr[rawInd];
       rawInd += 1;
   }
   
   // copying *10^ over
   for (int copyInd=0; copyInd < numBase10Chars; copyInd++)
       string[rawInd+copyInd] = BASE_TEN_FORMAT[copyInd];
       
   // copying over the rest of the original string (the exponent)
   for (rawInd++; rawInd < rawStrNumChars; rawInd++)
       string[rawInd + numBase10Chars - 1] = rawStr[rawInd];
   
   // clear 
   free(formatStr);
   free(rawStr);
   
   return string;
}

/**
 * returns a string containing a MMA array of numbers in
 * scientific notation with `precision` digits after each decimal point.
 * The returned string must be freed.
 * 
 * @brief converts double array to a MMA array string of sci-not decimals
 * @param array			an array of doubles to convert to sci-not
 * @param length		number of doubles in array
 * @param precision	number of digits after decimal point in sci-not
 * @return 				a MMA array of sci-not e.g."{1.2*10^+3, ...}",
 * 						which should be freed
 */
char* convertDoubleArrToMMA(double* array, int length, int precision) {
   
   // calculate array string length
   int sizeOfLeftBrace = strlen(ARRAY_OUTER_CHARS[0]);
   int sizeOfRightBrace = strlen(ARRAY_OUTER_CHARS[1]);
   int sizeOfDelim = strlen(ARRAY_DELIM_CHARS);
   int sizeOfNum = 11 + precision; // -1.(precision)*10^+(2 or 3)
   int sizeOfTerm = sizeOfNum + sizeOfDelim;
   int sizeOfArr = (
       + sizeOfLeftBrace			// {
       + sizeOfTerm*(length - 1)	// 	a, b, c,
       + sizeOfNum					//  d
       + sizeOfRightBrace);		// }
   
   char* arrString = malloc(sizeOfArr * sizeof(char));
   memset(arrString, '\0', sizeOfArr);
   
   // writte open-array character
   int strInd = 0;
   for (int outerInd=0; outerInd < sizeOfLeftBrace; outerInd++)
       arrString[strInd + outerInd] = ARRAY_OUTER_CHARS[0][outerInd];
   strInd += sizeOfLeftBrace;
   
   // for each number in the array
   for (int arrInd=0; arrInd < length; arrInd++) {
       
       // convert it to scientific notation
       char* numStr = getScientificNotation(array[arrInd], precision);
       int lenNumStr = strlen(numStr);
       
       // write the sci-not number to the string
       for (int charInd=0; charInd < lenNumStr; charInd++)
           arrString[strInd + charInd] = numStr[charInd];
       strInd += lenNumStr;
       
       // write the array delimiter to the string (except for last)
       if (arrInd < length - 1) {
           for (int delimInd=0; delimInd < sizeOfDelim; delimInd++)
               arrString[strInd + delimInd]  = ARRAY_DELIM_CHARS[delimInd];
           strInd += sizeOfDelim;
       }
       
       // free sci-not malloc
       free(numStr);
   }
   
   // write the close-array chars to the string
   for (int outerInd=0; outerInd < sizeOfRightBrace; outerInd++)
       arrString[strInd + outerInd] = ARRAY_OUTER_CHARS[1][outerInd];
   strInd += sizeOfRightBrace;
   
   return arrString;
}

/**
 * returns a file handle to be passed to subsequent functions for
 * writing data to a MMA Association. File handle must be eventually
 * passed to closeAssocWrite.
 * 
 * @brief begins writing an association to file
 * @param filename		name of the file to save the MMA association to.
 * 						The assocation is read into MMA by Get[filename]
 * @return 				a file handle to pass to subsequent Association
 * 						Write functions, and which must be eventually
 * 						passed to closeAssocWrite
 */
FILE* openAssocWrite(char* filename) {
   
   FILE* file = fopen(filename, "w");
   fprintf(file, "<|\n");
   return file;
}

/**
 * returns a file handle to be passed to subsequent functions for
 * writing data to a MMA Association, which appends elements to the
 * existing assoc in filename. File handle must be eventually
 * passed to closeAssocWrite or closeAssocAppend
 * 
 * @brief begins appending to an existing association in a file
 * @param filename		name of the file to containing the existing assoc
 * 						The assocation is read into MMA by Get[filename]
 * @return 				a file handle to pass to subsequent Association
 * 						Write functions, and which must be eventually
 * 						passed to closeAssocWrite or closeAssocAppend
 */
FILE* openAssocAppend(char* filename) {
   
   FILE* file = fopen(filename, "r+");
   
   // delete closing |> and newline
   fseek(file, -3, SEEK_END);
   
   // write comma for previous field
   fprintf(file, ",\n");
   
   return file;
}

/**
 * @brief formats and finalises the association, so it is ready
 * 		   to read by MMA's Get[filename] function
 * @param file	file handle returned by openAssocWrite
 */
void closeAssocWrite(FILE* file) {
	
	// delete comma (and newline necessarily)
	fseek(file, -2, SEEK_END);
	
	// restore newline and close Associonary
	fprintf(file, "\n|>");
	fclose(file);
}

/**
 * @brief equivalent to closeAssocWrite
 */
void closeAssocAppend(FILE* file) {
    closeAssocWrite(file);
}
 
/**
 * @brief adds an int to the association
 * @param file		file handle returned by openAssocWrite
 * @param keyname	key to add to the association
 * @param num		integer to add under key in association
 */
void writeIntToAssoc(FILE* file, char* keyname, int num) {
 	
    fprintf(file, "\"%s\" -> %d,\n", keyname, num);
}
 
/**
 * @brief adds a sci-not number to the association
 * @param file			file handle returned by openAssocWrite
 * @param keyname		key to add to the association
 * @param num			number to add in sci-not to the association
 * @param precision	number of digits after decimal point in sci-not
 */
void writeDoubleToAssoc(FILE* file, char* keyname, double num, int precision) {
	
	char* sciNotStr = getScientificNotation(num, precision);
	fprintf(file, "\"%s\" -> %s,\n", keyname, sciNotStr);
	free(sciNotStr);
}

/**
 * @brief adds a string to the association
 * @param file		file handle returned by openAssocWrite
 * @param keyname	key to add to the association
 * @param string	string to add under key in the association
 */
void writeStringToAssoc(FILE* file, char* keyname, char* string) {
	
	fprintf(file, "\"%s\" -> \"%s\",\n", keyname, string);
}

/**
 * @brief adds a MMA array of integers to the association
 * @param file		file handle returned by openAssocWrite
 * @param keyname	key to add to the association
 * @param arr		array of ints to add under key
 * @param length	length of the array of ints
 */
void writeIntArrToAssoc(FILE* file, char* keyname, int* arr, int length) {
	
	fprintf(file, "\"%s\" -> {", keyname);
	for (int i=0; i < length-1; i++)
		fprintf(file, "%d, ", arr[i]);
	fprintf(file, "%d},\n", arr[length-1]);
}

/**
 * @brief adds a once nested array of integers, where the inner arrays have
 * an inconsistent length, to the association
 * @param file				file handle returned by openAssocWrite
 * @param keyname			key to add to the association
 * @param arr				2D array of ints to add under key
 * @param outerLength		number of inner arrays
 * @param innerLengths		length of each inner array
 * @param innerSpace		the allocated size of the inner arrays 
 */
void writeUnevenOnceNestedIntArrToAssoc(FILE* file, char* keyname, int (*arr)[], int outerLength, int* innerLengths, int innerSpace) {
	
	// open outer list
	fprintf(file, "\"%s\" -> {\n", keyname);
	
	int* flatArr = (int *) arr;
	for (int outer=0; outer < outerLength; outer++) {
		
		// open inner list
		fprintf(file, "{");
		
		// write inner list
		for (int inner=0; inner < innerLengths[outer]; inner++)
			fprintf(file, "%d, ", *((flatArr+outer*innerSpace) + inner));
		
		// remove trailing comma and space (note, not written for empty lists)
		if (innerLengths[outer] > 0)
			fseek(file, -2, SEEK_END);
		
		// close inner list
		fprintf(file, "},\n");
	}
	
	// remove trailing comma and newline
	fseek(file, -2, SEEK_END);
	
	// close outer list
	fprintf(file, "\n},\n");
}

/**
 * @brief adds an array of long unsigned ints to the association 
 * 		  as a MMA array of integers
 * @param file		file handle returned by openAssocWrite
 * @param keyname	key to add to the association
 * @param arr		array of long unsigned ints to add under key
 * @param length	length of the array
 */
void writeUnsignedLongArrToAssoc(FILE* file, char* keyname, unsigned long *arr, int length) {
	
	fprintf(file, "\"%s\" -> {", keyname);
	for (int i=0; i < length-1; i++)
		fprintf(file, "%lu, ", arr[i]);
	fprintf(file, "%lu},\n", arr[length-1]);
}

/**
 * @brief adds a MMA array of sci-not numbers to the association
 * @param file			file handle returned by openAssocWrite
 * @param keyname		key to add to the association
 * @param arr			array of doubles to convert to sci-not
 * @param length		length of the array
 * @param precision	number of digits after decimal in sci-not
 */
void writeDoubleArrToAssoc(FILE* file, char* keyname, double* arr, int length, int precision) {
	
	char* arrStr = convertDoubleArrToMMA(arr, length, precision);
	fprintf(file, "\"%s\" -> %s,\n", keyname, arrStr);
	free(arrStr);
}

/**
 * @brief adds a once nested array of doubles, where the inner arrays have
 * an inconsistent length, to the association in scientific notation
 * @param file				file handle returned by openAssocWrite
 * @param keyname			key to add to the association
 * @param arr				2D array of doubles
 * @param outerLength		number of inner arrays in arr
 * @param innerLengths		length of each inner array
 * @param innerSpace		allocated space of each inner array
 * @param precision		number of digits after decimal in sci-not
 */
void writeUnevenOnceNestedDoubleArrToAssoc(
	FILE* file, char* keyname, double (*arr)[], int outerLength, int* innerLengths, int innerSpace, int precision) 
{
	// open outer list
	fprintf(file, "\"%s\" -> {\n", keyname);
	
	// flatten list
	double* flatArr = (double *) arr;
	
	// write each inner list to file with trailing comma and newline
	for (int outer=0; outer < outerLength; outer++) {
		char* arrStr = convertDoubleArrToMMA((flatArr+outer*innerSpace), innerLengths[outer], precision);
		fprintf(file, "%s,\n", arrStr);
		free(arrStr);
	}
	
	// delete trailing newline, comma
	fseek(file, -2, SEEK_END);
	
	// closer outer list
	fprintf(file, "\n},\n");
}
	
/**
 * @brief adds a once-nested MMA array of sci-not numbers to the association
 * @param file			file handle returned by openAssocWrite
 * @param keyname		key to add to the association
 * @param arr			nested array of doubles to convert to sci-not
 * @param outerLength	length of the outer list
 * @param innerLength	lenght of every inner list
 * @param precision	number of digits after decimal in sci-not
 */
void writeOnceNestedDoubleListToAssoc(
	FILE* file, char* keyname, double** arr, int outerLength, int innerLength, int precision) 
{
	// open outer list
	fprintf(file, "\"%s\" -> {\n", keyname);
	
	// write each inner list to file with trailing comma and newline
	for (int i=0; i < outerLength; i++) {
		char* arrStr = convertDoubleArrToMMA(arr[i], innerLength, precision);
		fprintf(file, "%s,\n", arrStr);
		free(arrStr);
	}
	
	// delete trailing newline, comma
	fseek(file, -2, SEEK_END);
	
	// closer outer list
	fprintf(file, "\n},\n");
}

/**
 * @brief adds a numDimensions-nested MMA array of sci-not numbers to the association
 * @param file			file handle returned by openAssocWrite
 * @param keyname		key to add to the association
 * @param arr			nested array of doubles to convert to sci-not
 * @param lengths		array of length of each dimension
 * @param precision	number of digits after decimal in sci-not
 */
void writeNestedDoubleListToAssoc(
	FILE* file, char* keyname, void* arr, int numDimensions, int* lengths, int precision
) {
	// start assoc
	fprintf(file, "\"%s\" -> ", keyname);
	
	writeNestedDoubleList(file, arr, numDimensions, lengths, 0, precision);

	// add comma and newline
	fprintf(file, ",\n");
}

/**
 * @brief adds a numDimensions-nested MMA array of sci-not numbers to the association
 * @param file			file handle returned by openAssocWrite
 * @param keyname		key to add to the association
 * @param arr			nested array of doubles to convert to sci-not
 * @param lengths		array of length of each dimension
 * @param innerTrimLength	
 * 						how many of each inner-most array to keep (to trim out place-holder data)
 * @param precision	number of digits after decimal in sci-not
 */
void writeNestedDoubleArrToAssoc(
	FILE* file, char* keyname, void* arr, int numDimensions, int* lengths, int innerTrimLength, int precision
) {
	// start assoc
	fprintf(file, "\"%s\" -> ", keyname);
	
	writeNestedDoubleArr(file, (double *) arr, 0, numDimensions, lengths, 0, innerTrimLength, precision);

	// add comma and newline
	fprintf(file, ",\n");
}
	
#endif // MMA_FORMATTER_H_

