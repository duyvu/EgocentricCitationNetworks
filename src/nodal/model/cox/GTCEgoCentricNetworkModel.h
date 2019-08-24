/*
 * GTCEgoCentricNetworkModel.h
 *
 *  Created on: Jun 1, 2010
 *      Author: duyvu
 */

#ifndef GTCEGOCENTRICNETWORKMODEL_H_
#define GTCEGOCENTRICNETWORKMODEL_H_

#include "EgoCentricNetworkModel.h"

namespace ndip {

class GTCEgoCentricNetworkModel: public ndip::EgoCentricNetworkModel {
protected:
	unsigned int gapTimeBetaIndex;
	void initializeMore() {
		gapTimeBetaIndex = data->numOfNodalNetworkStatistics;
	}

	double computeLogLikelihoodWithLBF(const DoubleVector& beta);
	double computeLogLikelihoodWithoutLBF(const DoubleVector& beta);

	double computeKappaOfNodalUpdateMap(int cacheIndex,
			const DoubleVector& beta, double currentEdgeTime);
	double computeKappaOfComingVertices(int cacheIndex,
			const DoubleVector& beta, double currentEdgeTime);

	void computeKappaDerivativesOfNodalUpdateMap(int cacheIndex,
			const DoubleVector& beta, double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa, double currentEdgeTime,
			bool resetKappas);
	void computeKappaDerivativesOfComingVertices(int cacheIndex,
			const DoubleVector& beta, double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa, double currentEdgeTime,
			bool resetKappas);
	void updateUandI(DoubleVector& u, DoubleMatrix& I, double kappa,
			const DoubleVector& deltaKappa,
			const DoubleMatrix& squaredDeltaKappa);

	void computeSchoenfeldResiduals(DoubleVector& residuals, double kappa,
			const DoubleVector& deltaKappa,
			const DoubleMatrix& squaredDeltaKappa);
public:
	GTCEgoCentricNetworkModel();
	virtual ~GTCEgoCentricNetworkModel();

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

#endif /* GTCEGOCENTRICNETWORKMODEL_H_ */
