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



/* number types */

typedef long long unsigned int index;

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

FORCE_INLINE index pow2(int p) {
    return (1ULL << p);
}

FORCE_INLINE index flipBit(index num, int i) {
    return num ^ pow2(i);
}

FORCE_INLINE int getBit(index num, int i) {
    return (num >> i) & 1;
}

FORCE_INLINE index insertZeroBit(index num, int i) {
    index l = (num >> i) << i;
    index r = num - l;
    return (l << 1ULL) ^ r;
}

FORCE_INLINE index getBitMask(int* bits, int numBits) {
    index mask = 0;
    for (int b=0; b<numBits; b++)
        mask = flipBit(mask, bits[b]);
    return mask;
}

FORCE_INLINE index truncateBits(index num, int numLowerBits) {
    return num & (pow2(numLowerBits) - 1);
}

FORCE_INLINE int bitsAreAllOne(index i, index mask) {
    return (mask & i) == mask;
}

FORCE_INLINE index getZeroBitFromAffix(index prefix, index suffix, int i) {
    return (prefix << (i+1)) | suffix;
}

FORCE_INLINE index getZeroBitsFromAffixes(index prefix, index infix, index suffix, int t2, int t1) {
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



/* state-vector management */

amp* createStatevector(int numQubits) {
    
    index numAmps = pow2(numQubits);
    amp* vec = malloc(numAmps * sizeof *vec);
    return vec;
}

void initRandomStatevector(amp* vec, int numQubits) {
    
    index numAmps = pow2(numQubits);
    double mag = 0;
    for (index i=0; i<numAmps; i++) {
        vec[i] = getRandomComplex(-1-I, 1+I);
        mag += getAbsSquared(vec[i]);
    }
    
    mag = sqrt(mag);
    for (index i=0; i<numAmps; i++)
        vec[i] /= mag;
}

void initOnesStatevector(amp* vec, int numQubits) {
    
    index numAmps = pow2(numQubits);
    for (index i=0; i<numAmps; i++)
        vec[i] = 1;
}

void printStatevector(amp* vec, int numQubits) {
    
    index numAmps = pow2(numQubits);
    for (index i=0; i<numAmps; i++)
        printf("psi[%llu] = %g + i(%g)\n", i, creal(vec[i]), cimag(vec[i]));
    printf("\n");
}

void printStatevectorForMMA(amp* vec, int numQubits) {
    
    index numAmps = pow2(numQubits);
    printf("{");
    for (index i=0; i<numAmps; i++) {
        printf("%.10f + I(%.10f)", creal(vec[i]), cimag(vec[i]));
        if (i < numAmps-1)
            printf(", ");
    }
    printf("}\n\n");
}



#endif // UTILITIES_H