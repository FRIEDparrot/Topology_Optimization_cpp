#ifndef FUN_ASSERTION
#define FUN_ASSERTION
#include <iostream>
#include <string.h>
#include <assert.h>    // this package is used for raise error

#define PI 3.1415

using namespace std;

// The definition of the class should be put here. 
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
	void set(int i,int j, value) {
		if (i >= xelm || j >= yelm || i<0 || j < 0) {
			throw std::runtime_error("the index out of boundary");
		}
		else {
			
		}
	}
	void Display_Content(){
		for (int i = 1; i < xelm; i++) {
			for (int j = 1; j < yelm; j++) {
				cout << *(*(this->Mat + i) + j) << " ";
			}
			cout << "\n";
		}
	};
private:
	int xelm;
	int yelm;
	int** Mat;
};

// class DyMat(int,int);
// the function in the other cpp source file can be asserted


#endif