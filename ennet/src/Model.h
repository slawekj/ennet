/*
 * Model.h
 *
 *  Created on: May 13, 2013
 *      Author: Janusz Slawek
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <vector>

using namespace std;

class Model {
public:
	Model(size_t M, size_t P) {
		this->M = M;
		this->P = P;
		this->f0 = 0.0;
		this->nu = 0.0;
		featSplitI.resize(M);
		featSplitT.resize(M);
		gamma_l.resize(M);
		gamma_r.resize(M);
		importance.resize(P);
	}

	virtual ~Model() {
	}

	size_t getM() const {
		return M;
	}

	size_t getP() const {
		return P;
	}

	double getF0() const {
		return f0;
	}

	void setF0(double f0) {
		this->f0 = f0;
	}

	double getNu() const {
		return nu;
	}

	void setNu(double nu) {
		this->nu = nu;
	}

	const double getGammaL(size_t tree) const {
		return gamma_l[tree];
	}

	void setGammaL(size_t tree, double g_l) {
		gamma_l[tree] = g_l;
	}

	const double getGammaR(size_t tree) const {
		return gamma_r[tree];
	}

	void setGammaR(size_t tree, double g_r) {
		gamma_r[tree] = g_r;
	}

	const size_t getFeatSplitI(size_t tree) const {
		return featSplitI[tree];
	}

	void setFeatSplitI(size_t tree, size_t ft) {
		featSplitI[tree] = ft;
	}

	const double getFeatSplitT(size_t tree) const {
		return featSplitT[tree];
	}

	void setFeatSplitT(size_t tree, double thr) {
		featSplitT[tree] = thr;
	}

	const double getImportance(size_t ft) const {
		return importance[ft];
	}

	void setImportance(size_t ft, double imp) {
		importance[ft] = imp;
	}

private:
	size_t M;
	size_t P;
	double f0;
	double nu;
	vector<size_t> featSplitI;
	vector<double> featSplitT;
	vector<double> gamma_l;
	vector<double> gamma_r;
	vector<double> importance;
};

#endif /* MODEL_H_ */
