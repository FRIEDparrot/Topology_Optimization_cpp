#ifndef FEM
#define FEM
#include <iostream>
#include <string.h>

#define PI 3.1415

using namespace std;

class DyMat {  // global parameters definition
public:
	DyMat(int xelm, int yelm) {
		//  initialize an array of dynamic size
		this->xelm = xelm;
		this->yelm = yelm;
		this->Mat = new int* [xelm];    
		// attention the initialization of the Matrix with a userdefined boudary should be placed in initialization function
		// define a pointer which is for a 2-D array -> 
		// in this case we can initialize an 2-D array with predefined-size
		
		for (int j = 0; j < yelm; j++) {
			*(this -> Mat + j) = new int[yelm];
		}
		for (int i = 0; i < xelm; i++) {
			for (int j = 0; j < yelm; j++) {
				*(*(this -> Mat + i) + j) = 0;   // initialize the corresponding element of the pointer array
			}
		}
	};
	void Display() {
		cout << "Haha";
	};
private:
	int xelm;
	int yelm;
	int** Mat;
};

// class DyMat(int,int);
// the function in the other cpp source file can be asserted


#endif