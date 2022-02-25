#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <sys/time.h>



/* we absoutely insist our inline demands are obeyed */
#define FORCE_INLINE static inline __attribute__((always_inline))



/* timing (repeats must be inside brace scope) */

#define START_TIMING() \
    struct timeval tval_before, tval_after, tval_result; \
    gettimeofday(&tval_before, NULL);

#define STOP_TIMING() \
    gettimeofday(&tval_after, NULL); \
    timersub(&tval_after, &tval_before, &tval_result); \
    printf("\t\t%ld.%06ld (s)\n", \
        (long int) tval_result.tv_sec, \
        (long int) tval_result.tv_usec);
        
#define RECORD_TIMING(VAR) \
    gettimeofday(&tval_after, NULL); \
    timersub(&tval_after, &tval_before, &tval_result); \
    VAR = tval_result.tv_sec + (1.0/1000000) * tval_result.tv_usec;
    
    

/* analysis */

void getAverageAndVariance(double* data, int lenData, double *av, double *var) {
    // assumes data is similar size, to safely use two-pass algorithm
    
    double s = 0;
    for (int i=0; i<lenData; i++)
        s += data[i];
    *av = s/lenData;
    
    s = 0;
    for (int i=0; i<lenData; i++)
        s += (data[i] - *av)*(data[i] - *av);
    *var = s/(lenData - 1);
}



/* number types */

typedef long long unsigned int INDEX;

typedef double complex amp;

FORCE_INLINE double getAbsSquared(amp val) {
    double r = creal(val);
    double i = cimag(val);
    return r*r + i*i;
}

FORCE_INLINE amp expI(double phase) {
    double r = cos(phase);
    double i = sin(phase);
    return r + i*I;
}



/* bit twiddling */

FORCE_INLINE INDEX pow2(int p) {
    return (1ULL << p);
}

FORCE_INLINE INDEX flipBit(INDEX num, int i) {
    return num ^ pow2(i);
}

FORCE_INLINE int getBit(INDEX num, int i) {
    return (num >> i) & 1;
}

FORCE_INLINE INDEX insertZeroBit(INDEX num, int i) {
    INDEX l = (num >> i) << (i+1);
    INDEX r = num & (pow2(i)-1);
    return l | r;
}

FORCE_INLINE INDEX getBitMask(int* bits, int numBits) {
    INDEX mask = 0;
    for (int b=0; b<numBits; b++)
        mask = flipBit(mask, bits[b]);
    return mask;
}

FORCE_INLINE INDEX truncateBits(INDEX num, int numLowerBits) {
    return num & (pow2(numLowerBits) - 1);
}

FORCE_INLINE int bitsAreAllOne(INDEX i, INDEX mask) {
    return (mask & i) == mask;
}

FORCE_INLINE INDEX getZeroBitFromAffix(INDEX prefix, INDEX suffix, int i) {
    return (prefix << (i+1)) | suffix;
}

FORCE_INLINE INDEX getZeroBitsFromAffixes(INDEX prefix, INDEX infix, INDEX suffix, int t2, int t1) {
    return (prefix << (t2+1)) | (infix << (t1+1)) | suffix;
} 



/* randomness */

double getRandomDecimal(double min, double max) {
    
    double r = rand() / (double) RAND_MAX;
    return min + r*(max-min);
}

amp getRandomComplex(amp min, amp max) {
    
    double re = getRandomDecimal(creal(min), creal(max));
    double im = getRandomDecimal(cimag(min), cimag(max));
    return re + I*im;
}

int getRandomInt(int min, int max) {
    
    return round(getRandomDecimal(min, max));
}

INDEX getRandomBitMask(int len, int numOnes) {
    
    INDEX mask = 0;
    for (int n=0; n<numOnes; n++) {
        
        int i = getRandomInt(0, len-1);
        while (getBit(mask, i))
            i = getRandomInt(0, len-1);
            
        mask = flipBit(mask, i);
    }
    return mask;
}

void getSortedRandomSubReg(int* ret, int SubRegSize, int RegSize) {
    
    INDEX mask = getRandomBitMask(RegSize, SubRegSize);
    
    int q = 0;
    for (int i=0; i<SubRegSize; i++) {
        while (getBit(mask, q) != 1)
            q++;
        ret[i] = q++;
    }
}



/* state-vector management */

amp* createStatevector(int numQubits) {
    
    INDEX numAmps = pow2(numQubits);
    amp* vec = malloc(numAmps * sizeof *vec);
    return vec;
}

void initRandomStatevector(amp* vec, int numQubits) {
    
    INDEX numAmps = pow2(numQubits);
    double mag = 0;
    for (INDEX i=0; i<numAmps; i++) {
        vec[i] = getRandomComplex(-1-I, 1+I);
        mag += getAbsSquared(vec[i]);
    }
    
    mag = sqrt(mag);
    for (INDEX i=0; i<numAmps; i++)
        vec[i] /= mag;
}

void initOnesStatevector(amp* vec, int numQubits) {
    
    INDEX numAmps = pow2(numQubits);
    for (INDEX i=0; i<numAmps; i++)
        vec[i] = 1;
}



/* printing */

void printIntArray(char* label, int* arr, int len) {
    
    printf("%s[%d] = {", label, len);
    for (int i=0; i<len-1; i++)
        printf("%d, ", arr[i]);
    printf("%d}\n", arr[len-1]);
}

void printStatevector(amp* vec, int numQubits) {
    
    INDEX numAmps = pow2(numQubits);
    for (INDEX i=0; i<numAmps; i++)
        printf("psi[%llu] = %g + i(%g)\n", i, creal(vec[i]), cimag(vec[i]));
    printf("\n");
}

void printStatevectorForMMA(amp* vec, int numQubits) {
    
    INDEX numAmps = pow2(numQubits);
    printf("{");
    for (INDEX i=0; i<numAmps; i++) {
        printf("%.10f + I(%.10f)", creal(vec[i]), cimag(vec[i]));
        if (i < numAmps-1)
            printf(", ");
    }
    printf("}\n\n");
}



#endif // UTILITIES_H