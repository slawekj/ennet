/*
 * regression_stump.h
 *
 *  Created on: May 1, 2013
 *      Author: Janusz Slawek
 */

#ifndef REGRESSION_STUMP_H_
#define REGRESSION_STUMP_H_

#include "Model.h"
#include "Prediction.h"

const Model train_regression_stump(const int N_train, const int P_train,
		const double *x, const double *y, const double col_sampling_rate,
		const double row_sampling_rate, const int M_train, const double nu);

const Prediction test_regression_stump(const Model &m, const int N_test,
		const double *x, const double *y, const int M_test);

#endif /* REGRESSION_STUMP_H_ */
