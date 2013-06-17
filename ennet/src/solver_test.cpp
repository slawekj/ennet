/*
 * solver_test.cpp
 *
 *  Created on: May 4, 2013
 *      Author: Janusz Slawek
 */

#include <iostream>
#include <cstdlib>
#include "regression_stump.h"

#define N 100
#define P 5

using namespace std;

int main() {

	// dummy test
	double x[N * P];
	for (int i = 0; i < N * P; i++) {
		x[i] = rand() % 1000;
	}
	double y[N];
	for (int i = 0; i < N; i++) {
		y[i] = rand() % 1000;
	}

	train_regression_stump(N, P, x, y, 0.5, 0.5, 1, 0.01);
	cout << "done" << endl;
}
