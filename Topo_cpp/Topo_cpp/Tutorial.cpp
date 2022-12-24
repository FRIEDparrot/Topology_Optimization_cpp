/// This cpp file is not relevant to other source files
void Te1(double* A) {
    *A = *A + 2;
    // * is the sign to modify the parameter that the main function gives
    // Te1(&c);
    // cout << c;
}


void Test() {
    int t = 5;
    int* Array = new int[t];
    for (int i = 0; i < t; i++)
    {
        *(Array + i) = 5;
    }
    // cout << *(Array + 2);
    ///
    /*
    DyMat<double> Test1(3, 3, (double*) new double[3][3]{ {2,2,3},{1,-1, 0},{-1,2,1 } });
    DyMat<double> Test2(3, 1, (double*) new double[3][1]{ {7},{-1},{4} });
    DyArray<int> A(3, new int[] {0, 1, 2});
    DyArray<int> B(1, new int[] {1});
    Test1.Display_Content();
    Test1.setChunk(A, B, Test2);
    Test1.Display_Content();
    */
}