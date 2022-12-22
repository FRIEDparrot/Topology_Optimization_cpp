#pragma once
#ifndef FUN_ASSERTION
#define FUN_ASSERTION
#include <iostream>
#include <string.h>
#include <assert.h>    // this package is used for raise error
#define PI 3.1415

using namespace std;
void FEM(int xelm, int yelm);

// The definition of the class should be put here. 

template <typename ItemType>
class DyArray {
public:
	int length;
	DyArray(int length, ItemType InitValue) {
		this->length = length;
		this->Array = new ItemType[length];
		for (int i = 0; i < length; i++) {
			*(this->Array + i) = InitValue;
		}
	}
	void set(int index, ItemType value) {
		try {
			*(this->Array + index) = value;
		}
		catch (runtime_error &e) {
			throw std::runtime_error("an error occurred!");
			cout << e.what();
		}
	}
	ItemType get(int index) {
		return *(this -> Array + index);
	}
	void Display_Content() {
		for (int i = 0; i < length; i++) {
			cout << *(this -> Array + i) << " ";
		}
		cout << endl;
	}
private:
	ItemType* Array;
};

DyArray<int> Apply_fixed_constraint(DyArray<int> freedof, int nelx, int nely);

// -----------------dynamic array----------------------
template <typename ItemType>
class DyMat {  // global parameters definition
	public:
		int xelm;
		int yelm;
		DyMat(int xelm, int yelm, ItemType InitValue){
			//  initialize an array of dynamic size
			this->xelm = xelm;
			this->yelm = yelm;
			this->Mat = new ItemType* [xelm];
			// attention the initialization of the Matrix with a userdefined boudary should be placed in initialization function
			// define a pointer which is for a 2-D array -> 
			// in this case we can initialize an 2-D array with predefined-size
			for (int i = 0; i < xelm; i++) {
				*(this -> Mat + i) = new ItemType[yelm];
			}
			for (int i = 0; i < xelm; i++) {
				for (int j = 0; j < yelm; j++) {
					*(*(this -> Mat + i) + j) = InitValue;   // initialize the corresponding element of the pointer array
				}
			}
		};
		void set(int i, int j, ItemType value) {
			try {
				*(*(this->Mat + i) + j) = value;
			}
			catch (runtime_error& e) {
				throw std::runtime_error("the index out of boundary");
				cout << e.what();
			}
		};
		// This funciton is used to set a peice of chunk;
		void setChunk(int lines[], int cols[], DyMat<ItemType> dmatrix ) {
			// there we initilize the lines[] as a pointer --> so we can't get it's size

			int lineNum = dmatrix.xelm;
			int colNum = dmatrix.yelm;
			for (int i = 0; i < dmatrix.xelm; i++) {
				cout << i;
			};
			cout << "\n";
			for (int i = 0; i < lineNum ; i++) {
				for (int j = 0; j < colNum ; j++){
					int x_i = lines[i];
					int y_i = cols[j];
					cout << x_i << " " << y_i << endl;
				}
			}
			
		};

		ItemType get(int i, int j) {
			return *(*(this->Mat + i) + j);
		}
		void Display_Content(){
			for (int i = 0; i < xelm; i++) {
				for (int j = 0; j < yelm; j++) {
					cout << *(*(this->Mat + i) + j) << " ";
				}
				cout << "\n";
			}
		};
	private:
		ItemType** Mat;
};


// class DyMat(int,int);
// the function in the other cpp source file can be asserted


#endif