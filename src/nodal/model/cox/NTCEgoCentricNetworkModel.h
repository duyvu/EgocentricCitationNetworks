/*
 * NTCEgoCentricNetworkModel.h
 *
 *  Created on: May 27, 2010
 *      Author: duyvu
 */

#ifndef NTCEGOCENTRICNETWORKMODEL_H_
#define NTCEGOCENTRICNETWORKMODEL_H_

#include "EgoCentricNetworkModel.h"

namespace ndip {

class NTCEgoCentricNetworkModel: public ndip::EgoCentricNetworkModel {
protected:
	void initializeMore() {
		return;
	}

	double computeLogLikelihoodWithLBF(const DoubleVector& beta);
	double computeLogLikelihoodWithoutLBF(const DoubleVector& beta);
	double computeKappaOfNodalUpdateMap(int cacheIndex,
			const DoubleVector& beta);
	double computeKappaOfComingVertices(int cacheIndex,
			const DoubleVector& beta);
	void computeKappaDerivativesOfNodalUpdateMap(int cacheIndex,
			const DoubleVector& beta, double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa, bool resetKappas);
	void computeKappaDerivativesOfComingVertices(int cacheIndex,
			const DoubleVector& beta, double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa, bool resetKappas);
public:
	NTCEgoCentricNetworkModel();
	virtual ~NTCEgoCentricNetworkModel();

	double computeLogLikelihood(const DoubleVector& beta, bool isUsingLBF);
	double computeLogLikelihoodAndScoreVector(DoubleVector& u,
			const DoubleVector& beta);
	void computeScoreVectorAndInformationMatrix(DoubleVector& u,
			DoubleMatrix& I, const DoubleVector& beta);

	void getNodalCovariateMatrix(DoubleMatrix& covariates, double time);
	void computeCummulativeBaselineHazard(
			VectorOfCummulativePoints& cumBH, const DoubleVector& beta);
	void computeCoxSnellResiduals(
			VectorOfCummulativePoints& csResiduals,
			const DoubleVector& beta);
	void computeMartingaleResiduals(
			VectorOfCummulativePoints& martingaleResiduals,
			VectorOfCummulativePoints& timeMartingaleResiduals,
			const DoubleVector& beta);
	void computeSchoenfeldResiduals(DoubleMatrix& schoenfeldResiduals,
			const DoubleVector& beta, bool isPlusBeta);
	void getTrainingRanks(DoubleVector& ranks, const DoubleVector& beta,
			const std::string& tieMethod);
	void simulateNetwork(Graph& graph, DoubleMatrix& covariates,
			const DoubleVector& beta);
	void simulateNetworkWithBothEnds(Graph& graph, DoubleMatrix& covariates,
			const DoubleVector& betaOut, const DoubleVector& betaIn);
};

}

#endif /* NTCEGOCENTRICNETWORKMODEL_H_ */
