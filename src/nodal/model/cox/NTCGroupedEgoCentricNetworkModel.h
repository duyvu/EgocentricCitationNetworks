/*
 * NTCGroupedEgoCentricNetworkModel.h
 *
 *  Created on: Jan 18, 2011
 *      Author: duyvu
 */

#ifndef NTCGROUPEDEGOCENTRICNETWORKMODEL_H_
#define NTCGROUPEDEGOCENTRICNETWORKMODEL_H_

#include "GroupedEgoCentricNetworkModel.h"

namespace ndip {

class NTCGroupedEgoCentricNetworkModel: public ndip::GroupedEgoCentricNetworkModel {
protected:
	// this class does need any specialized subclass initialization
	void initializeMore() {
		return;
	}

	// for computing the log-likelihood
	double computeLogLikelihoodWithLBF(const DoubleVector& beta);
	double computeSumNodalStatistics(unsigned int nodeID,
			const DoubleVector& beta);
	double computeKappaOfNodalUpdateMap(int cacheIndex,
			const DoubleVector& beta);
	double computeKappaOfComingVertices(int cacheIndex,
			const DoubleVector& beta);

	// for computing the score vector and the Hessian matrix
	void computeNodalDerivatives(unsigned int nodeID, const DoubleVector& beta,
			double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa);
	void computeKappaDerivativesOfNodalUpdateMap(int cacheIndex,
			const DoubleVector& beta, double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa, bool resetKappas,
			const NodeSet& citingNodes);
	void computeKappaDerivativesOfComingVertices(int cacheIndex,
			const DoubleVector& beta, double& kappa, DoubleVector& deltaKappa,
			DoubleMatrix& squaredDeltaKappa, bool resetKappas);
	void
			updateUandI(DoubleVector& u, DoubleMatrix& I, const double& kappa,
					const DoubleVector& deltaKappa,
					const DoubleMatrix& squaredDeltaKappa,
					unsigned int numOfTiedEvents);
public:
	NTCGroupedEgoCentricNetworkModel();
	virtual ~NTCGroupedEgoCentricNetworkModel();

	double computeLogLikelihood(const DoubleVector& beta, bool isUsingLBF);
	void computeLogLikelihoodOverEvents(DoubleVector& cummulativeLogLLH,
			const DoubleVector& beta, bool isUsingLBF);
	void computeLogLikelihoodOverEdges(ofstream& outputFile,
			const DoubleVector& beta, bool isUsingLBF);

	void computeScoreVectorAndInformationMatrix(DoubleVector& u,
			DoubleMatrix& I, const DoubleVector& beta);

	void getSumRanks(DoubleVector& sumRanks, const DoubleVector& beta,
			const std::string& tieMethod);
	void getRanks(VectorOfDoubleVector& ranks, const DoubleVector& beta,
			const std::string& tieMethod);
	void getRanks(ofstream& byEdgesOutputFile, ofstream& byEventsOutputFile,
			ofstream& accumlativeOutputFile, const DoubleVector& beta,
			const std::string& tieMethod);

};

}

#endif /* NTCGROUPEDEGOCENTRICNETWORKMODEL_H_ */
