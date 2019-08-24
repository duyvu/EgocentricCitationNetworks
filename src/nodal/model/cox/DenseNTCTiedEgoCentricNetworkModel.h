/*
 * DenseNTCTiedEgoCentricNetworkModel.h
 *
 *  Created on: Jan 13, 2011
 *      Author: duyvu
 */

#ifndef DENSENTCTIEDEGOCENTRICNETWORKMODEL_H_
#define DENSENTCTIEDEGOCENTRICNETWORKMODEL_H_

#include "TiedEgoCentricNetworkModel.h"

namespace ndip {

class DenseNTCTiedEgoCentricNetworkModel: public ndip::TiedEgoCentricNetworkModel {
protected:
	// this class does need any specialized subclass initialization
	void initializeMore() {
		return;
	}
	// likelihood computation
	double computeLogLikelihoodWithoutLBF(const DoubleVector& beta);
	double computeSumNodalStatistics(unsigned int nodeID,
			const DoubleVector& beta);
	// for computing the score vector and the Hessian matrix
	void computeNodalDerivatives(unsigned int nodeID, const DoubleVector& beta,
			double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa);
	void
			updateUandI(DoubleVector& u, DoubleMatrix& I, const double& kappa,
					const DoubleVector& deltaKappa,
					const DoubleMatrix& squaredDeltaKappa,
					unsigned int numOfTiedEvents);

public:
	DenseNTCTiedEgoCentricNetworkModel();
	virtual ~DenseNTCTiedEgoCentricNetworkModel();

	double computeLogLikelihood(const DoubleVector& beta, bool isUsingLBF);
	void computeLogLikelihoodOverEvents(DoubleVector& cummulativeLogLLH,
			const DoubleVector& beta, bool isUsingLBF);
	void computeLogLikelihoodOverEdges(ofstream& byEdgesOutputFile,
			const DoubleVector& beta, bool isUsingLBF);

	void computeScoreVectorAndInformationMatrix(DoubleVector& u,
			DoubleMatrix& I, const DoubleVector& beta);
	double computeLogLikelihoodAndScoreVector(DoubleVector& u,
			const DoubleVector& beta);
	void
			computeOnlineGradientAscentMLE(DoubleVector& beta,
					double learningRate);

	void computeCummulativeBaselineHazard(VectorOfCummulativePoints& cumBH,
			const DoubleVector& beta);

	void computeSchoenfeldResiduals(DoubleMatrix& schoenfeldResiduals,
			const DoubleVector& beta, bool isPlusBeta);

	void getSumRanks(DoubleVector& sumRanks, const DoubleVector& beta,
			const std::string& tieMethod);
	void getRanks(VectorOfDoubleVector& ranks, const DoubleVector& beta,
			const std::string& tieMethod);
	void getRanks(ofstream& byEdgesOutputFile, ofstream& byEventsOutputFile,
				ofstream& accumlativeOutputFile, const DoubleVector& beta,
				const std::string& tieMethod);

	void simulateNetwork(Graph& graph, DoubleMatrix& covariates,
			const DoubleVector& beta);

};

}

#endif /* DENSENTCTIEDEGOCENTRICNETWORKMODEL_H_ */
