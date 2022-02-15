
/* Comparing the performance of local simulation of the QFT via 
 * a direct evaluation of the gates, vs one where contiguous 
 * controlled-phase gates have been merged into a single diagonal.
 * Note that all operators (hadamards, swaps, phases) have been implemented
 * optimally, that is without branching or superfluous memory access.
 * This file does not include multithreading (yet!)
 *
 * It is important that inline functions are actually inlined, so pass optimisation
 * flags (especially with Clang)
 * run with:
 *      gcc local_qft.c -O1 -o test; ./test
 */

#include "utilities.h"



/* gates */

void applyHadamard(amp* psi, int t, int N) {
    
    const double fac = 1/sqrt(2);
    
    const index jNum = pow2(N-(t+1));
    const index kNum = pow2(t);
    
    for (index j=0; j<jNum; j++) {
        for (index k=0; k<kNum; k++) {
    
            // |j>|0>|k> and |j>|1>|k>
            index j0k = getZeroBitFromAffix(j, k, t);
            index j1k = flipBit(j0k, t);
            
            amp a1 = psi[j0k];
            amp a2 = psi[j1k];
            
            psi[j0k] = fac*a1 + fac*a2;
            psi[j1k] = fac*a1 - fac*a2;
        }
    }
}

void applyControlledPhase(amp* psi, int c, int t, double theta, int N) {
    const int t1 = (t < c)? t : c;
    const int t2 = (c > t)? c : t;
    
    const amp fac = expI(theta);
    
    // phase shift |j>|1>|k>|1>|l>
    const index jNum = pow2(N-(t2+1));
    const index kNum = pow2(t2-(t1+1));
    const index lNum = pow2(t1);
    
    for (index j=0; j<jNum; j++) {
        for (index k=0; k<kNum; k++) {
            for (index l=0; l<lNum; l++) {
                
                index j0k0l = getZeroBitsFromAffixes(j, k, l, t2, t1);
                index j1k1l = flipBit(flipBit(j0k0l, t2), t1);
                
                psi[j1k1l] *= fac;
            }
        }
    } 
}

void applySwap(amp* psi, int t1, int t2, int N) {
    if (t1 > t2) {
        int t = t1;
        t1 = t2;
        t2 = t;
    }
    
    // |j>|0>|k>|1>|l> <-> |j>|1>|k>|0>|l>
    const index jNum = pow2(N-(t2+1));
    const index kNum = pow2(t2-(t1+1));
    const index lNum = pow2(t1);
    
    for (index j=0; j<jNum; j++) {
        for (index k=0; k<kNum; k++) {
            for (index l=0; l<lNum; l++) {
                
                index j0k0l = getZeroBitsFromAffixes(j, k, l, t2, t1);
                index j0k1l = flipBit(j0k0l, t1);
                index j1k0l = flipBit(j0k0l, t2);

                amp tmp = psi[j0k1l];
                psi[j0k1l] = psi[j1k0l];
                psi[j1k0l] = tmp;
            }
        }
    }
}



/* QFT by circuit */

void applyMultiplePhases(amp* psi, int tMax, int N) {
    
    int m = 2;
    for (int t=tMax-1; t>=0; t--) {
        double theta = 2 * M_PI / (double) pow2(m);
        applyControlledPhase(psi, tMax, t, theta, N);
        m++;
    }
}

void applyQFTCircuit(amp* psi, int N) {
    
    for (int t=N-1; t>0; t--) {
        applyHadamard(psi, t, N);
        applyMultiplePhases(psi, t, N);
    }
    applyHadamard(psi, 0, N);
    
    for (int t=0; t<N/2; t++)
        applySwap(psi, t, N-t-1, N);
}



/* QFT by algorithm */

void applyMergedPhases(amp* psi, int tMax, int N) {

    // |j>|1>|k>
    const index jNum = pow2(N-(tMax+1));
    const index kNum = pow2(tMax);
    const index kMask = kNum-1;
    
    const double fac = (M_PI / (double) kNum);
    
    for (index j=0; j<jNum; j++) {
        for (index k=0; k<kNum; k++) {
            
            index j0k = getZeroBitFromAffix(j, k, tMax);
            index j1k = flipBit(j0k, tMax);
            
            double theta = fac * (j1k & kMask);
            psi[j1k] *= expI(theta);
        }
    }
}

void applyQFTAlgorithm(amp* psi, int N) {
    
    for (int t=N-1; t>0; t--) {
        applyHadamard(psi, t, N);
        applyMergedPhases(psi, t, N);
    }
    applyHadamard(psi, 0, N);
    
    for (int t=0; t<N/2; t++)
        applySwap(psi, t, N-t-1, N);
}



/* launch */

int main() {
    
    int N = 24;
    amp* psi = createStatevector(N);
    initRandomStatevector(psi, N);
    printf("[%d qubits]\n\n", N);    
    
    printf("contiguous phases\n");
    {
        printf("\tas N gates\n");
        START_TIMING()
        applyMultiplePhases(psi, N-1, N);
        STOP_TIMING()
    }
    {
        printf("\tas 1 merged gate\n");
        START_TIMING()
        applyMergedPhases(psi, N-1, N);
        STOP_TIMING()
    }
    
    printf("QFT\n");
    {
        printf("\tusing full circuit\n");
        START_TIMING()
        applyQFTCircuit(psi, N);
        STOP_TIMING()
    }
    {   
        printf("\tusing merged phases\n");
        START_TIMING()
        applyQFTAlgorithm(psi, N);
        STOP_TIMING()
    }
    

    free(psi);
    return 0;
}