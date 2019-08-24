/*
 * ATCTiedEgoCentricNetworkModel.h
 *
 *  Created on: Oct 28, 2010
 *      Author: duyvu
 */

#ifndef ATCTIEDEGOCENTRICNETWORKMODEL_H_
#define ATCTIEDEGOCENTRICNETWORKMODEL_H_

#include "TiedEgoCentricNetworkModel.h"

namespace ndip {

class ATCTiedEgoCentricNetworkModel: public ndip::TiedEgoCentricNetworkModel {

protected:
	unsigned int ageTimeBetaIndex;
	void initializeMore() {
		ageTimeBetaIndex = data->numOfNodalNetworkStatistics;
	}

	// for computing the log-likelihood
	double computeLogLikelihoodWithLBF(const DoubleVector& beta);
	double computeLogLikelihoodWithoutLBF(const DoubleVector& beta);
	double computeSumNodalStatistics(unsigned int nodeID,
			const DoubleVector& beta, double currentCitingTime);
	double computeKappaOfNodalUpdateMap(int cacheIndex,
			const DoubleVector& beta, double currentCitingTime);
	double computeKappaOfComingVertices(int cacheIndex,
			const DoubleVector& beta, double currentCitingTime);

	// for computing the score vector and the Hessian matrix
	void computeNodalDerivatives(unsigned int nodeID, const DoubleVector& beta,
			double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa, double currentCitingTime);
	void computeKappaDerivativesOfNodalUpdateMap(int cacheIndex,
			const DoubleVector& beta, double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa, double currentCitingTime,
			bool resetKappas);
	void computeKappaDerivativesOfComingVertices(int cacheIndex,
			const DoubleVector& beta, double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa, double currentCitingTime,
			bool resetKappas);
	void
			updateUandI(DoubleVector& u, DoubleMatrix& I, const double& kappa,
					const DoubleVector& deltaKappa,
					const DoubleMatrix& squaredDeltaKappa,
					unsigned int numOfTiedEvents);

public:
	ATCTiedEgoCentricNetworkModel();
	virtual ~ATCTiedEgoCentricNetworkModel();

	double computeLogLikelihood(const DoubleVector& beta, bool isUsingLBF);
	void computeLogLikelihoodOverEvents(DoubleVector& cummulativeLogLLH,
			const DoubleVector& beta, bool isUsingLBF);
	void computeLogLikelihoodOverEdges(ofstream& byEdgesOutputFile,
				const DoubleVector& beta, bool isUsingLBF);

	void computeScoreVectorAndInformationMatrix(DoubleVector& u,
			DoubleMatrix& I, const DoubleVector& beta);

	double computeLogLikelihoodAndScoreVector(DoubleVector& u,
				const DoubleVector& beta);

	void computeOnlineGradientAscentMLE(DoubleVector& beta,
				double learningRate);

	void computeCummulativeBaselineHazard(VectorOfCummulativePoints& cumBH,
			const DoubleVector& beta);
	void computeSchoenfeldResiduals(DoubleMatrix& schoenfeldResiduals,
			const DoubleVector& beta, bool isPlusBeta);
	void getSumRanks(DoubleVector& sumRanks, const DoubleVector& beta,
			const std::string& tieMethod);
	void getRanks(VectorOfDoubleVector& ranks,
				const DoubleVector& beta, const std::string& tieMethod);
	void simulateNetwork(Graph& graph, DoubleMatrix& covariates,
			const DoubleVector& beta);
};

}

#endif /* ATCTIEDEGOCENTRICNETWORKMODEL_H_ */
