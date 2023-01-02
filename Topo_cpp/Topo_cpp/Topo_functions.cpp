#include "Fun_assertion.h"

DyMat<double> Stiffness_Update(DyMat<double> K, DyMat<double> K_reshape, int proj_seq[]) {
    // This funciton is only used in the step to calculate Data.K() in the FEMupdate function
    DyMat<double> K_mdf(K);
    for (int i = 0; i < K_reshape.xelm; i++) {
        for (int j = 0; j < K_reshape.yelm; j++) {
            K_mdf.set(proj_seq[i], proj_seq[j], K.get(proj_seq[i], proj_seq[j]) + K_reshape.get(i, j));
        }
    }
    return  K_mdf;
}

void FEMupdate(FEM_parameters Data, DyMat<double> Ke, DyMat<double> x, int nelx, int nely, double penal)
{
    // penal is the penalization factor
    for (int elx = 0; elx < nelx; elx++) {
        for (int ely = 0; ely < nely; ely++) {
            // calculate the number of the element  --> only for y
            int n1 = (nely + 1) * elx + ely + 1;   // left-top point of the element
            int n2 = (nely + 1) * (elx + 1) + ely + 1;
            int seq[] = { 2 * n1 - 2, 2 * n1 - 1 , 2 * n2 - 2, 2 * n2 - 1,  2 * n2, 2 * n2 + 1 , 2 * n1, 2 * n1 + 1};
            // define the projection sequencce
            // int columns[] = { 2 * n1 - 2, 2 * n1 - 1 , 2 * n2 - 2, 2 * n2 - 1,  2 * n2, 2 * n2 + 1 , 2 * n1, 2 * n1 + 1 };
            // then we update the stiffness matrix K
            DyMat<double> dK(8, 8, pow(x.get(ely,elx),penal));
            // Begin Matrix Projection here
            DyMat<double> K_pre = Stiffness_Update(Data.K(), DyMat_dotMul(Ke, pow(x.get(ely, elx), penal)) , seq);
            // ***** attention if we want to change the K to a totally new matrix, don't use data.K() = K2; 
            // ***** Use data.set_K(K2) instead; 
            // set the corresponding Chunk of the stiffness matrix;
        };
    };

    // matlab code ---->  U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
    // count the number of the freedofs and build an array that 

    int freeArrLength = Data.Boundaries().sum();
    int index = 0;

    DyArray<int> freeArr(freeArrLength, (int)0);
    for (int i = 0; i < Data.Boundaries().length; i++) {
        if (Data.Boundaries().get(i) != 0) {  // for every freedofs
            freeArr.set(index, i);
            index = index + 1;
        }
    }

    //  the problem is that we can't  use get() to obtain the value of the DyArray<int> type;
    // ***** Data.U() = DyMatMul(DyMat_Inv(Data.K()), Data.F());
    DyArray<int> Sequential(Data.U().yelm, 0);  // use zero to initialize the array
    int _p = 0;
    for (int i = 0; i < Sequential.length; i++) {
        Sequential.set(i, _p);
        _p += 1;    // note that if the program triggered a breakpoint is because the array exceed the boundary
    } // all the columns of the DyMat

    Data.U().setChunk(freeArr, Sequential, DyMatMul(DyMat_Inv(Data.K().getChunk(freeArr,freeArr)), Data.F().getChunk(freeArr, Sequential)));
    // it would be 0 if only U has 1 line
    // also the U(fixdoofs) = 0 should be set, but it's not necessary because U is initialized by zero
}

DyMat<double> ChangeUpdate(DyMat<double> U,DyMat<double> Ke, DyMat<double> x, double* c, DyMat<double> dc, double penal, double rmin) {
    // in this function we update the changes and make a filter of the net
    int nelx = dc.yelm;
    int nely = dc.xelm;
    // attention the dc(yelm, xelm);

    // I just use the original code written in Matlab
    for (int ely = 0; ely < nely; ely++) {
        for (int elx = 0; elx < nelx; elx++) {
            int n1 = (nely + 1) * elx + ely + 1;  
            int n2 = (nely + 1) * (elx + 1) + ely + 1;   // n2 = 5 firstly
            int lines  [] = { 2 * n1 - 2, 2 * n1 - 1, 2 * n2 - 2, 2 * n2 - 1, 2 * n2, 2 * n2 + 1, 2 * n1, 2 * n1 + 1};
            int columns[] = { 0 };
            // attention this should be modified if the dimension of the U is changed
            DyMat<double> Ue = U.getChunk(lines,columns,8,1);  // this is for get a Chunk of size (8x1)
            double delta = DyMatMul(DyMatMul(DyMat_transpose(Ue),Ke),Ue).get(0,0);
            // this time we can only calculate the Ue^T Ke Ue once
            *c += pow(x.get(ely, elx), penal) * delta;
            dc.set(ely, elx, -penal * pow(x.get(ely, elx), penal -1) * delta);
        }
    }
    return dc;
}

// Mesh-independency filter
DyMat<double> Grid_filter(DyMat<double> x, DyMat<double> dc, double rmin) {
    int nelx = dc.xelm; int nely = dc.yelm;

    DyMat<double> dc_new(nelx, nely, (double)0 );   // dc_new(elmy,elmx)
    int r_filter;
    if (rmin < 1) {
        r_filter = 1;
    }
    else {
        r_filter = (int)floor(rmin);
    }

    for (int i = 0; i < nelx; i++) {
        for (int j = 0 ; j < nely; j++) {
            double _sum = 0;   // i,j is the current coordination of the point
            for (int k = max(i - r_filter , 0); k <= min(nelx-1, i + r_filter); k++) {   // yelm
                for (int p = max( j - r_filter , 0); p <= min(nely-1, j + r_filter); p++) {  // xelm
                    // fac means the weight --> this should be over zero
                    double fac = max(rmin - sqrt( pow((i-k),2.) + pow((j-p),2.)), (double)0);
                    _sum += fac;
;                   dc_new.set(i, j, dc_new.get(i, j) + fac * x.get(k, p) * dc.get(k, p));
                    // Note that x and c,dc are all the dymat of elmy*elmx
                }
            }
            dc_new.set(i, j, dc_new.get(i, j) /(_sum * x.get(i, j )));
        }
    }
    return dc_new;
}

FEM_parameters FEMInit(int xelm, int yelm) {
    // initialize the dynamic_matrix of K, F and U
    // K -global stiffness
    DyMat<double> K(2 * (xelm + 1) * (yelm + 1), 2 * (xelm + 1) * (yelm + 1), 0.);
    DyMat<double> F(2 * (xelm + 1) * (yelm + 1), 1, 0.);
    DyMat<double> U(2 * (xelm + 1) * (yelm + 1), 1, 0.);
    // =========initialize the boudary condition =============
    F.set(1, 0, -1.);
    // set the force boundary condition 
    DyArray<int> alldofs(2 * (xelm + 1) * (yelm + 1), 1);  // firstly we set all dimension as free
    // apply the constraint in the cantilever beam or simplily supported beam; 
    alldofs = Apply_fixed_constraint(alldofs, xelm, yelm);
    FEM_parameters para(K, F, U, alldofs);
    return para;
    // --------------- matlab code ---------------
    // fixeddofs = [1: 2 : 2 * (nely + 1)];
    // alldofs = [1:2 * (nely + 1) * (nelx + 1)];
    // freedofs = setdiff(alldofs, fixeddofs);
    // fixeddofs = union([1:2 : 2 * (nely + 1)], [2 * (nelx + 1) * (nely + 1)]);
};

DyArray<int> Apply_fixed_constraint(DyArray<int> freedof, int xelm, int yelm) {
    // apply a fixed constraint on the left of the beam;  --> that part is the definition of the matrix
    for (int i = 0; i < 2*(yelm + 1) ; i+=2 ) {
        freedof.set(i , 0);
    }
    freedof.set(2 * (xelm + 1) * (yelm + 1) - 1, 0);  // apply fixed constraint
    // attention the last element should be length-1
    return freedof;
}

// Generate the Ke matrix ----> finished
DyMat<double> Elm_Stiff_Generate(double E, double nu) {
    double coef = E / (1. - pow(nu, 2.0));

    double k[] = { coef * (1. / 2 - nu / 6), coef * (1. / 8 + nu / 8),
                       coef * (-1. / 4 - nu / 12), coef * (-1. / 8 + 3. * nu / 8) ,
                       coef * (-1. / 4 + nu / 12), coef * (-1. / 8 - nu / 8),
                       coef * (nu / 6) , coef * (1. / 8 - 3 * nu / 8) };
    /*don't forget the . or the computer would consider 1, 2.... as the integral
    */
    DyMat<double> k_t(1, 8, (double*)k);

    double Ke_pre[8][8] = {
    {k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]},
    {k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]},
    {k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]},
    {k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]},
    {k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]},
    {k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]},
    {k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]},
    {k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]}
    };

    DyMat<double> Ke(8, 8, (double*)Ke_pre);
    // if we just input **Ke it like this it will regard it as a value 
    // input the pointer of the 2-D array
    return Ke;
};

DyMat<double> Optimization_Criterion(DyMat<double> x, DyMat<double> dc, double volfrac) {
    double Lim1 = 0; double Lim2 = 100000;  double MaxMove = 0.2;
    int nelx = x.xelm; int nely = x.yelm;
    double mid = 0.5 * (Lim2 + Lim1);   // every time we have to recalculate the mid
    // multiplication is faster than the division
    // use the bi-sectional algorithm to find the lambda (Lagrange Operator)
    DyMat<double> x_new = DyMat_dotMul(x, DyMat_dotPow(DyMat_dotMul(dc, -1/mid), 0.5));

    while (Lim2 - Lim1 > 0.0001) {
        mid = 0.5 * (Lim2 + Lim1);   // every time we have to recalculate the mid
        x_new = DyMat_dotMul(x, DyMat_dotPow(DyMat_dotMul(dc, -1 / mid), 0.5));
        // since mid is different every time in the iteration
        // then we judge what the x_new(i,j) should be set as
        for (int i = 0; i < nelx; i++) {
            for (int j = 0; j < nely; j++) {
                // In the iteration process, the data should be less than 1 but greater than 0
                if (x_new.get(i, j) < max(0.001, x.get(i,j) - MaxMove)) {
                    x_new.set(i, j , max(0.001, x.get(i, j) - MaxMove));
                }
                else if (x_new.get(i, j) > min(1., x.get(i, j) + MaxMove)) {
                    x_new.set(i, j, min(1., x.get(i, j) + MaxMove));
                }
            }
        }
        if ((DyMat_sum(x_new)) > (double)nelx * (double)nely * volfrac) {
            Lim1 = mid;
        }
        else {
            Lim2 = mid;
        }
        DyMat<double> x(x_new);   // renew the x Matrix
    }
    return x_new;
}
