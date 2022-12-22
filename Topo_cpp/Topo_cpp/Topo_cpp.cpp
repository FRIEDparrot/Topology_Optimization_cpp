// Topo_cpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include <iostream>
#include "Fun_assertion.h"    // attention that should

int main()
{
    FEM(4, 3);

    // we had "using namespace std" in the head file
    return 0;
}

void FEM(int xelm, int yelm) {
    // initialize the dynamic_matrix
    DyMat<double> K(2 * (xelm + 1) * (yelm + 1), 2 * (xelm + 1) * (yelm + 1), 0.);  // this is the method to initialize 
    DyMat<double> F(2 * (xelm + 1), 1, 0.);
    DyMat<double> U(2 * (xelm + 1), 1, 0.);
    DyMat<double> TestMat(2, 2, 5.3);
    // =========initialize the boudary condition =============
    F.set(2, 1, -1.);
    // set the value of the F matrix;
    DyArray<int> alldofs(2 * (xelm + 1) * (yelm + 1), 1);  // firstly we set all dimension as free
    Apply_fixed_constraint(alldofs, xelm, yelm);

    int xline[] {0,1};
    int yline[] {2,3};
    K.setChunk(xline, yline,TestMat);


    // alldofs.Display_Content();  ---> then we get the 


    // --------------- matlab code ---------------
    // fixeddofs = [1: 2 : 2 * (nely + 1)];
    // alldofs = [1:2 * (nely + 1) * (nelx + 1)];
    // freedofs = setdiff(alldofs, fixeddofs);
    // fixeddofs = union([1:2 : 2 * (nely + 1)], [2 * (nelx + 1) * (nely + 1)]);
};

DyArray<int> Apply_fixed_constraint(DyArray<int> freedof, int xelm, int yelm) {
    // apply a fixed constraint on the left of the beam;  --> that part is the definition of the matrix
    for (int i = 0; i < yelm + 1; i ++) {
        freedof.set(2*i*( xelm + 1)+1, 0);
    }
    freedof.set(2*(xelm + 1) * (yelm + 1) - 1, 0); 
    // attention the last element should be length-1
    return freedof;
}

