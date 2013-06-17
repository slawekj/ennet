/*
 * test_regression_stump_R.cpp
 *
 *  Created on: May 14, 2013
 *      Author: Janusz Slawek
 */

#include "regression_stump.h"

extern "C" {

/*
 * This is a simple R interface
 * input variables are passed as constant pointers
 * results are returned as non-constant pointers
 */
void test_regression_stump_R(const int *N_test, const int *P_test,
		const int *P_train, const double *x_test, const double *y_test,
		const int *M_test, const int *M_train, const double *nu,
		const double *f0, const int *featI, const double * featT,
		const double *gamma_l, const double *gamma_r, double *loss, double *p) {

	/*
	 * R use:
	 * int*           - as.logical()
	 * int*           - as.integer()
	 * double*        - as.double()
	 * Rcomplex*      - as.complex()
	 * char**         - as.character()
	 * unsigned char* - as.raw()
	 */

	/*
	 * First step is to recreate Model from raw arrays
	 */
	Model m(*M_train, *P_test);
	m.setNu(*nu);
	m.setF0(*f0);
	for (int i = 0; i < *M_train; i++) {
		m.setFeatSplitI(i, featI[i]);
		m.setFeatSplitT(i, featT[i]);
		m.setGammaL(i, gamma_l[i]);
		m.setGammaR(i, gamma_r[i]);
		//leave importance blank
	}

	/*
	 * Now calculate a prediction
	 */
	Prediction prediction = test_regression_stump(m, *N_test, x_test, y_test,
			*M_test);

	for (int i = 0; i < *M_test; i++) {
		loss[i] = prediction.getLoss(i);
	}
	for (int i = 0; i < *N_test; i++) {
		p[i] = prediction.getPrediction(i);
	}
}

}

