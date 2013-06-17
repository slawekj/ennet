/*
 * ennet_R.cpp
 *
 *  Created on: May 3, 2013
 *      Author: Janusz Slawek
 */

#include "regression_stump.h"

extern "C" {

/*
 * This is a simple R interface
 * input variables are passed as constant pointers
 * results are returned as non-constant pointers
 */
void train_regression_stump_R(const int *N_train, const int *P_test,
		const double *x_train, const double *y_train, const double *s_f,
		const double *s_s, const int *M_train, const double *nu, double *I,
		double *f0, int *featI, double * featT, double *gamma_l,
		double *gamma_r) {

	/*
	 * R use:
	 * int*           - as.logical()
	 * int*           - as.integer()
	 * double*        - as.double()
	 * Rcomplex*      - as.complex()
	 * char**         - as.character()
	 * unsigned char* - as.raw()
	 */

	Model m = train_regression_stump(*N_train, *P_test, x_train, y_train, *s_f,
			*s_s, *M_train, *nu);
	for (int i = 0; i < *P_test; i++) {
		I[i] = m.getImportance(i);
	}
	*f0 = m.getF0();
	for (int i = 0; i < *M_train; i++) {
		featI[i] = (int) m.getFeatSplitI(i);
		featT[i] = m.getFeatSplitT(i);
		gamma_l[i] = m.getGammaL(i);
		gamma_r[i] = m.getGammaR(i);
	}
}

}
