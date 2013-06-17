/*
 * solver.h
 *
 *  Created on: May 1, 2013
 *      Author: Janusz Slawek
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "Model.h"
#include "Prediction.h"

const Model train_regression_stump(const int N_train, const int P_train,
		const double *x, const double *y, const double col_sampling_rate,
		const double row_sampling_rate, const int M_train, const double nu);

const Prediction test_regression_stump(const Model &m, const int N_test,
		const double *x, const double *y, const int M_test);

#endif /* SOLVER_H_ */
