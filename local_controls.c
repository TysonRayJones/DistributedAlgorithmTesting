
/* Comparing the performance of several local techniques for 
 * incorporating single (s_) and multiple (m_) control qubits into 
 * simulation of a single-target unitary gate. Note the extension
 * to distributed simulation requires edge cases for all methods (except 
 * A & B), to handle the scenarios when none or all local amplitudes
 * pass the control condition.
 *
 * As a matter of convenience, this function will use real arrays in lieu of 
 * complex statevectors.
 *
 * It is important that inline functions are actually inlined, so pass optimisation
 * flags (especially with Clang)
 * run with:
 *      gcc local_controls.c -O3 -lm -fopenmp -o test
 *      ./test
 */

#include "utilities.h"
#include "mmaformatter.h"

#ifdef _OPENMP
#include <omp.h>
#endif



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

#define f(amp) (1.5 * pow(amp - .1, 2))



/* single control methods */

void s_methodA(double* amps, INDEX numAmps, int c) {
    INDEX i;
    #pragma omp parallel for shared(amps,numAmps,c) private(i) schedule(static)
    for (i=0; i<numAmps; i++)
        if (getBit(i, c))
            amps[i] = f(amps[i]);
}

void s_methodB(double* amps, INDEX numAmps, int c) {
    INDEX i; int b;
    #pragma omp parallel for shared(amps,numAmps,c) private(i,b) schedule(static)
    for (i=0; i<numAmps; i++) {
        b = getBit(i, c);
        amps[i] = (1-b)*amps[i] + b*f(amps[i]);
    }
}

void s_methodC(double* amps, INDEX numAmps, int c) {
    INDEX jNum = numAmps >> (c+1);
    INDEX iNum = pow2(c);
    INDEX j,i,j0i,j1i;
    #pragma omp parallel for shared(amps,numAmps,c,jNum,iNum) private(j,i,j0i,j1i) schedule(static) collapse(2)
    for (j=0; j<jNum; j++) {
        for (i=0; i<iNum; i++) {
            j0i = getZeroBitFromAffix(j, i, c);
            j1i = flipBit(j0i, c);
            amps[j1i] = f(amps[i]);
        }
    }
}

void s_methodD(double* amps, INDEX numAmps, int c) {
    INDEX l1 = numAmps >> 1;
    INDEX m,i;
    #pragma omp parallel for shared(amps,numAmps,c,l1) private(m,i) schedule(static)
    for (m=0; m<l1; m++) {
        i = flipBit(insertZeroBit(m, c), c);
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
    INDEX i;
    #pragma omp parallel for shared(amps,numAmps,ctrls,numCtrls,cMask) private(i) schedule(static)
    for (i=0; i<numAmps; i++)
        if (bitsAreAllOne(i, cMask))
            amps[i] = f(amps[i]);
}

void m_methodB(double* amps, INDEX numAmps, int* ctrls, int numCtrls) {
    INDEX cMask = getBitMask(ctrls, numCtrls);
    INDEX i; int b;
    #pragma omp parallel for shared(amps,numAmps,ctrls,numCtrls,cMask) private(i,b) schedule(static)
    for (i=0; i<numAmps; i++) {
        b = bitsAreAllOne(i, cMask);
        amps[i] = (1-b)*amps[i] + b*f(amps[i]);
    }
}

void m_methodD(double* amps, INDEX numAmps, int* ctrls, int numCtrls) {
    INDEX lNum = numAmps >> numCtrls;
    INDEX l,j; int c;
    #pragma omp parallel for shared(amps,numAmps,ctrls,numCtrls,lNum) private(l,j,c) schedule(static)
    for (l=0; l<lNum; l++) {
        j=l;
        for (c=0; c<numCtrls; c++)
            j = flipBit(insertZeroBit(j, ctrls[c]), ctrls[c]);
        amps[j] = f(amps[j]);
    }
}

void (*m_methods[3]) (double* amps, INDEX numAmps, int* ctrls, int numCtrls) = {
    m_methodA, m_methodB, m_methodD
};
char* m_methodNames[3] = {"A", "B", "D"};




/* launch */

void simpleTest() {
    
    int numQubits = 27;
    INDEX numAmps = (1LL << numQubits);
    double* amps = malloc(numAmps * sizeof *amps);
    printf("[%d qubits]\n\n", numQubits);
    
    initArray(amps, numAmps);
    
    
    printf("single control\n");

    int c = 2;
    
    for (int m=0; m<4; m++) {
        
        printf("%s\n", s_methodNames[m]);
        
        START_TIMING()
        s_methods[m](amps, numAmps, c);
        STOP_TIMING()
    }
    
    
    printf("multiple controls\n");
    
    int numCtrls = 10;
    int ctrls[] = {0,2,4,6,7,15,16,20,21,22}; // must be increasing
    
    for (int m=0; m<3; m++) {
        
        printf("%s\n", m_methodNames[m]);
        
        initArray(amps, numAmps);
        
        START_TIMING()
        m_methods[m](amps, numAmps, ctrls, numCtrls);
        STOP_TIMING()
    }
    
    
    free(amps);
}


void s_benchmarking(int numQubits, int numReps, char* outFN) {
    
    int outPrec = 5;
    INDEX numAmps = (1LL << numQubits);
    double* amps = malloc(numAmps * sizeof *amps);
    printf("[%d qubits]\n\n", numQubits);
    
    // you MUST init array before benchmarking, because the very first write to 
    // heap memory has an overhead on some platforms! Funky!
    initArray(amps, numAmps);
    
    double durs[4][numQubits];
    double vars[4][numQubits];
    
    for (int m=0; m<4; m++) {
        
        for (int c=0; c<numQubits; c++) {
            
            double totalDur = 0;
            double totalDurSquared = 0;
            
            for (int r=0; r<numReps; r++) {
                initArray(amps, numAmps);
                
                START_TIMING()
                s_methods[m](amps, numAmps, c);
                RECORD_TIMING(double dur);
                
                totalDur += dur;
                totalDurSquared += dur*dur;
            }
            
            durs[m][c] = totalDur/numReps;
            vars[m][c] = (totalDurSquared/numReps) - durs[m][c]*durs[m][c];
        }
    }
    
    
    FILE* file = openAssocWrite(outFN);
    writeStringToAssoc(file, "note", "timings are already per-rep");
    writeIntToAssoc(file, "numQubits", numQubits);
    writeIntToAssoc(file, "numReps", numReps);
    writeIntToAssoc(file, "outPrec", outPrec);
    for (int m=0; m<4; m++) {
        char buff[50];
        strcpy(buff, "dur_"); strcat(buff, s_methodNames[m]);
        writeDoubleArrToAssoc(file, buff, durs[m], numQubits, outPrec);
        strcpy(buff, "var_"); strcat(buff, s_methodNames[m]);
        writeDoubleArrToAssoc(file, buff, vars[m], numQubits, outPrec);
    }
    closeAssocWrite(file);
    
    free(amps);
}




void m_benchmarking(int numQubits, int numReps, char* outFN) {
    
    
    int outPrec = 10;
    INDEX numAmps = (1LL << numQubits);
    double* amps = malloc(numAmps * sizeof *amps);
    printf("[%d qubits]\n\n", numQubits);
    
    // you MUST init array before benchmarking, because the very first write to 
    // heap memory has an overhead on some platforms! Funky!
    initArray(amps, numAmps);
    
    
    /* we obviously cannot try every {ctrls}, since there are O(2^(numQubits))
     * choices! So we'll try every number of controls from 1 to all, and for each,
     * try #numReps different random qubit indices. We'll make sure that each method 
     * is trying the identical controls assignments too, so we don't have to worry about 
     * sampling thresholds.
     */
     
    double durs[3][numQubits+1];
    double vars[3][numQubits+1];
    for (int m=0; m<3; m++) {
        for (int i=0; i<(numQubits+1); i++) {
            durs[m][i] = -1;
            vars[m][i] = -1;
        }
    }

    
    // try every possible #numCtrls
    for (int numCtrls=2; numCtrls<=numQubits; numCtrls++) {
        int ctrls[numCtrls];
        
        // obtain data
        double totalDurs[3];
        double totalDursSquared[3];
        
        for (int r=0; r<numReps; r++) {            
            getSortedRandomSubReg(ctrls, numCtrls, numQubits);
            
            for (int m=0; m<3; m++) {                
                initArray(amps, numAmps);    
                            
                START_TIMING();
                m_methods[m](amps, numAmps, ctrls, numCtrls);
                RECORD_TIMING(double dur);
                
                totalDurs[m] += dur;
                totalDursSquared[m] += dur*dur;
            }
        }
        
        // process data
        for (int m=0; m<3; m++) {
            durs[m][numCtrls] = totalDurs[m]/numReps;
            vars[m][numCtrls] = (totalDursSquared[m]/numReps) - durs[m][numCtrls]*durs[m][numCtrls];
        }
    }
    
    
    FILE* file = openAssocWrite(outFN);
    writeStringToAssoc(file, "note", "timings are already per-rep");
    writeIntToAssoc(file, "numQubits", numQubits);
    writeIntToAssoc(file, "numReps", numReps);
    writeIntToAssoc(file, "outPrec", outPrec);
    for (int m=0; m<3; m++) {
        char buff[50];
        strcpy(buff, "dur_"); strcat(buff, m_methodNames[m]);
        writeDoubleArrToAssoc(file, buff, durs[m], numQubits+1, outPrec);
        strcpy(buff, "var_"); strcat(buff, m_methodNames[m]);
        writeDoubleArrToAssoc(file, buff, vars[m], numQubits+1, outPrec);
    }
    closeAssocWrite(file);
    
    
    free(amps);
}







int main(int argc, char* argv[]) {
    
    srand(123456789);
    
    if (argc == 1)
        simpleTest();
        
    else if (argc == 5 && argv[1][0] == 's') {
        int numQubits = atoi(argv[2]);
        int numReps = atoi(argv[3]);
        char* outFN = argv[4];
        s_benchmarking(numQubits, numReps, outFN);
    
    } else if (argc == 5 && argv[1][0] == 'm') {
        int numQubits = atoi(argv[2]);
        int numReps = atoi(argv[3]);
        char* outFN = argv[4];
        m_benchmarking(numQubits, numReps, outFN);
        
    } else
        printf("call as either:\n\t./exec\n\t./exec [s/m] numQubits numReps outFN\n");
    
    return 0;
}