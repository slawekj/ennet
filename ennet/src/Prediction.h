/*
 * Result.h
 *
 *  Created on: May 14, 2013
 *      Author: Janusz Slawek
 */

#ifndef RESULT_H_
#define RESULT_H_

#include <vector>

using namespace std;

class Prediction {
public:
	Prediction(size_t M, size_t N) {
		this->M = M;
		this->N = N;
		loss.resize(M);
		prediction.resize(N);
	}

	virtual ~Prediction() {
	}

	size_t getM() const {
		return M;
	}

	const double getLoss(size_t tree) const {
		return loss[tree];
	}

	void setLoss(size_t tree, double l) {
		loss[tree] = l;
	}

	const double getPrediction(size_t row) const {
		return prediction[row];
	}

	void setPrediction(size_t row, double p) {
		prediction[row] = p;
	}

private:
	size_t M;
	size_t N;
	vector<double> prediction;
	vector<double> loss;
};

#endif /* RESULT_H_ */
