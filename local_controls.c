
/* Comparing the performance of several local techniques for 
 * incorporating single (s_) and multiple (m_) control qubits into 
 * simulation of a single-target unitary gate. Note the extension
 * to distributed simulation requires edge cases for all methods (except 
 * A & B), to handle the scenarios when none or all local amplitudes
 * pass the control condition. This file does not include multithreading (yet!)
 *
 * As a matter of convenience, this function will use real arrays in lieu of 
 * complex statevectors.
 *
 * It is important that inline functions are actually inlined, so pass optimisation
 * flags (especially with Clang)
 * run with:
 *      gcc local_controls.c -O1 -o test; ./test
 */

#include "utilities.h"



/* array management */

void printSubArray(double* amps, INDEX numAmps) {
    for (INDEX i=0; i<numAmps; i++) {
        printf("amp[%lld] = %g\n", i, amps[i]);
    }
    printf("\n");
}

void initArray(double* amps, INDEX numAmps) {
    for (INDEX i=0; i<numAmps; i++)
        amps[i] = 1;
}



/* a stand-in function for modifying amplitudes */

#define f(amp) (1.5 * pow(amp - .1,2))



/* single control methods */

void s_methodA(double* amps, INDEX numAmps, int c) {
    for (INDEX i=0; i<numAmps; i++)
        if (getBit(i, c))
            amps[i] = f(amps[i]);
}

void s_methodB(double* amps, INDEX numAmps, int c) {
    for (INDEX i=0; i<numAmps; i++) {
        int b = getBit(i, c);
        amps[i] = (1-b)*amps[i] + b*f(amps[i]);
    }
}

void s_methodC(double* amps, INDEX numAmps, int c) {
    const INDEX jNum = numAmps >> (c+1);
    const INDEX iNum = pow2(c);
    for (INDEX j=0; j<jNum; j++) {
        for (INDEX i=0; i<iNum; i++) {
            INDEX j0i = getZeroBitFromAffix(j, i, c);
            INDEX j1i = flipBit(j0i, c);
            amps[j1i] = f(amps[i]);
        }
    }
}

void s_methodD(double* amps, INDEX numAmps, int c) {
    const INDEX l1 = numAmps >> 1;
    for (INDEX m=0; m<l1; m++) {
        INDEX i = flipBit(insertZeroBit(m, c), c);
        amps[i] = f(amps[i]);
    }
}


void (*s_methods[4]) (double* amps, INDEX numAmps, int c) = {
    s_methodA, s_methodB, s_methodC, s_methodD
};
char* s_methodNames[4] = {"A", "B", "C", "D"};



/* many control methods */

void m_methodA(double* amps, INDEX numAmps, int* ctrls, int numCtrls) {
    INDEX cMask = getBitMask(ctrls, numCtrls);
    for (INDEX i=0; i<numAmps; i++)
        if (bitsAreAllOne(i, cMask))
            amps[i] = f(amps[i]);
}

void m_methodB(double* amps, INDEX numAmps, int* ctrls, int numCtrls) {
    INDEX cMask = getBitMask(ctrls, numCtrls);
    for (INDEX i=0; i<numAmps; i++) {
        int b = bitsAreAllOne(i, cMask);
        amps[i] = (1-b)*amps[i] + b*f(amps[i]);
    }
}

void m_methodD(double* amps, INDEX numAmps, int* ctrls, int numCtrls) {
    const INDEX l1 = numAmps >> numCtrls;
    for (INDEX l=0; l<l1; l++) {
        INDEX j=l;
        for (int c=0; c<numCtrls; c++)
            j = flipBit(insertZeroBit(j, ctrls[c]), ctrls[c]);
        amps[j] = f(amps[j]);
    }
}

void (*m_methods[3]) (double* amps, INDEX numAmps, int* ctrls, int numCtrls) = {
    m_methodA, m_methodB, m_methodD
};
char* m_methodNames[3] = {"A", "B", "D"};




/* launch */

int main() {
    
    int numQubits = 27;
    INDEX numAmps = (1LL << numQubits);
    double* amps = malloc(numAmps * sizeof *amps);
    printf("[%d qubits]\n\n", numQubits);
    
    
    printf("single control\n");

    int c = 2;
    
    for (int m=0; m<4; m++) {
        
        printf("%s\n", s_methodNames[m]);
        
        START_TIMING()
        
        initArray(amps, numAmps);
        s_methods[m](amps, numAmps, c);

        STOP_TIMING()
    }
    
    
    printf("multiple controls\n");
    
    int numCtrls = 10;
    int ctrls[] = {0,2,4,6,7,15,16,20,21,22}; // must be increasing
    
    for (int m=0; m<3; m++) {
        
        printf("%s\n", m_methodNames[m]);
        
        START_TIMING()
        
        initArray(amps, numAmps);
        m_methods[m](amps, numAmps, ctrls, numCtrls);
        
        STOP_TIMING()
    }
    
    
    free(amps);
    
    return 0;
}