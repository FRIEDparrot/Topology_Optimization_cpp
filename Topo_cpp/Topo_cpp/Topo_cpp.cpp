// Topo_cpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
// we had "using namespace std" in the head file
// all the function types that include it can use the relevant functions and types
#include "Fun_assertion.h";    
#include "Plot_Functions.h";
#include <conio.h>;



int main()
{
#pragma region  Parameter Definition
    double max_iteration = 1;
    double cut_error = 0.01;
    int elmx = 10;
    int elmy = 3;
    double volfrac = 0.4;
    double E = 1.;
    double nu = 0.3;
    double p = 3.;
    double rmin = 1.2;
    /// <summary>
    /// Parameters :
    ///     E = 1.0; nu = 0.3; volfrac = 0.4; penalization_factor = 3.0;
    ///     xelm = 50; yelm = 30; rmin = 1.2;
    /// </summary>
#pragma endregion

    DyMat<double> Ke = Elm_Stiff_Generate(E, nu);     // we firstly initialize the stiffness matrix
    // initialize the denisity matrix
    DyMat<double> x(elmy, elmx, volfrac);
    // the optimization  object is to minimize c(  c = U^T K_e U )
    int iteration = 0;
    double Error = 1.;

    while ( Error > cut_error && iteration < max_iteration) {
        iteration += 1;
        FEM_parameters FEM_Data = FEMInit(elmx, elmy);   // Init the K, F and U
        FEMupdate(FEM_Data, Ke, x, elmx, elmy, p);
        cout << "run successfully" << endl;
        double c = 0; // I set the breakpoint here


        // DyMat<double> dc(elmy, elmx, 0.);
        // dc = ChangeUpdate(FEM_Data.U(), Ke, x, &c, dc, p, rmin);

        // DyMat<double> dc_new = Grid_filter(x, dc, rmin);

        // x = Optimization_Criterion(x, dc, volfrac);
        // logs : FEMupdate will crush when the lines reach 125; ---> where ely = 30; ---> fix( < nely)
        // DyMat<double> dc( elmx, elmy , (double)0);   // every time we recalculate the dymap dc
    }

    // DyMat<double> Test1(3, 3, (double*) new double[3][3]{ {2,2,3},{1,-1, 0},{-1,2,1 } });
    // DyMat<double> Test2(3, 1, (double*) new double[3][1]{ {7},{-1},{4} });
    // DyMatMul(DyMat_Inv(Test1), Test2).Display_Content();
    // char b = _getch();
    return 0;
}