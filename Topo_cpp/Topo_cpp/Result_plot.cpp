#include "Fun_assertion.h";
#include <conio.h>;
#include <graphics.h>;
#include <stdio.h>;

// note that EasyX should be installed firstly before running
void PlotResult(DyMat<double> x) {
	initgraph(700, 300);
	// note that x is a matrix of size nely * nelx
	for (int i = 0; i < x.yelm; i++) {
		for (int j = 0; j < x.xelm; j++) {
			// total length --> 750  total height --> 250
			setfillcolor(RGBtoGRAY(RGB(-x.get(j,i) * 256, -x.get(j,i) * 256 ,-x.get(j,i) * 256 )));
			fillrectangle(25 + i * 650 / x.yelm , 25 + j * 250 / x.xelm , 25+(i+1)*650/x.yelm, 25+ (j+1) * 250/ x.xelm);
			
		}
	}
}
