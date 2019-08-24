/*
 * NTCTiedEgoCentricNetworkModel.h
 *
 *  Created on: Jul 19, 2010
 *      Author: duyvu
 */

#ifndef NTCTIEDEGOCENTRICNETWORKMODEL_H_
#define NTCTIEDEGOCENTRICNETWORKMODEL_H_

#include "TiedEgoCentricNetworkModel.h"

namespace ndip {

class NTCTiedEgoCentricNetworkModel: public ndip::TiedEgoCentricNetworkModel {
protected:
	// this class does need any specialized subclass initialization
	void initializeMore() {
		return;
	}

	// for computing the log-likelihood
	double computeLogLikelihoodWithLBF(const DoubleVector& beta);
	void computeLogLikelihoodOverEventsWithLBF(DoubleVector& cummulativeLogLLH,
			const DoubleVector& beta);
	double computeLogLikelihoodWithoutLBF(const DoubleVector& beta);
	double computeSumNodalStatistics(unsigned int nodeID,
			const DoubleVector& beta);
	double computeKappaOfNodalUpdateMap(int cacheIndex,
			const DoubleVector& beta);
	double computeKappaOfComingVertices(int cacheIndex,
			const DoubleVector& beta);

	// for computing the gradient only
	void computeNodalGradients(unsigned int nodeID, const DoubleVector& beta,
			double& kappa, DoubleVector& deltaKappa);
	void computeKappaGradientsOfNodalUpdateMap(int cacheIndex,
			const DoubleVector& beta, double& kappa, DoubleVector& deltaKappa,
			bool resetKappas, unsigned int citingNode);
	void computeKappaGradientsOfComingVertices(int cacheIndex,
			const DoubleVector& beta, double& kappa, DoubleVector& deltaKappa,
			bool resetKappas);
	void updateU(DoubleVector& u, const double& kappa,
			const DoubleVector& deltaKappa, unsigned int numOfTiedEvents);

	// for computing the score vector and the Hessian matrix
	void computeNodalDerivatives(unsigned int nodeID, const DoubleVector& beta,
			double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa);
	void computeKappaDerivativesOfNodalUpdateMap(int cacheIndex,
			const DoubleVector& beta, double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa, bool resetKappas,
			unsigned int citingNode);
	void computeKappaDerivativesOfComingVertices(int cacheIndex,
			const DoubleVector& beta, double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa, bool resetKappas);
	void
			updateUandI(DoubleVector& u, DoubleMatrix& I, const double& kappa,
					const DoubleVector& deltaKappa,
					const DoubleMatrix& squaredDeltaKappa,
					unsigned int numOfTiedEvents);

	void
			computeSchoenfeldResiduals(DoubleVector& residuals, double kappa,
					const DoubleVector& deltaKappa,
					const DoubleMatrix& squaredDeltaKappa,
					unsigned int numOfTiedEvents);
public:
	NTCTiedEgoCentricNetworkModel();
	virtual ~NTCTiedEgoCentricNetworkModel();

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

#endif /* NTCTIEDEGOCENTRICNETWORKMODEL_H_ */
