/*
 * DataTypes.h
 *
 *  Created on: Dec 27, 2009
 *      Author: duyvu
 */

#define DEBUG_LEVEL_1
#define DEBUG_EGO_CENTRIC_NETWORK_DATA
#define DEBUG_AP_EGO_CENTRIC_NETWORK_DATA

#ifndef DATATYPES_H_
#define DATATYPES_H_

#include <cmath>

#include <map>
#include <vector>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>

#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>

using namespace std;
using namespace boost;
namespace ublas = boost::numeric::ublas;

namespace ndip {

#define MIN_TAU 1e-15

#define NODAL_POPULARITY_ACTIVITY_MODEL 1
#define NODAL_TRANSITIVITY_POPULARITY_ACTIVITY_MODEL 2
#define NODAL_CLUSTERING_POPULARITY_ACTIVITY_MODEL 3

#define DYADIC_TRANSITIVITY_RECIPROCITY_MODEL 4

#define PREFERENTIAL_ATTACHMENT_EGO_CENTRIC_MODEL	0
#define ACTIVITY_POPULARITY_EGO_CENTRIC_MODEL	1
#define TRANSITIVITY_ACTIVITY_POPULARITY_EGO_CENTRIC_MODEL	2
#define NODAL_TRANSITIVITY_ACTIVITY_POPULARITY_EGO_CENTRIC_MODEL	3
#define SIMULATING_TRANSITIVITY_ACTIVITY_POPULARITY_EGO_CENTRIC_MODEL	4
#define SIMULATING_NODAL_TRANSITIVITY_ACTIVITY_POPULARITY_EGO_CENTRIC_MODEL	5
#define P2P_EGO_CENTRIC_MODEL	6
#define P2PT_EGO_CENTRIC_MODEL	7
#define P2PTR_EGO_CENTRIC_MODEL	8
#define TextualP2PTR_EGO_CENTRIC_MODEL	9
#define LDA_P2PTR_EGO_CENTRIC_MODEL	10
#define LDA_EGO_CENTRIC_MODEL	11
#define PA_GAP_EGO_CENTRIC_MODEL	12
#define SIMULATING_P2PTR_EGO_CENTRIC_MODEL	13

#define TRANSITIVITY_RECIPROCITY_DEGREE_RELATIONAL_EVENT_NETWORK_MODEL	1
#define FAVOR_POST_TRANSITIVITY_RECIPROCITY_DEGREE_RELATIONAL_EVENT_NETWORK_MODEL	2

#define DEGREE_CONSTRAINED_RELATIONAL_EVENT_NETWORK_MODEL	1
#define SHARED_PARTNER_DEGREE_CONSTRAINED_RELATIONAL_EVENT_NETWORK_MODEL	2
#define EGO_POSTS_SHARED_PARTNER_DEGREE_CONSTRAINED_RELATIONAL_EVENT_NETWORK_MODEL	3
#define RECENCY_SHARED_PARTNER_DEGREE_CONSTRAINED_RELATIONAL_EVENT_NETWORK_MODEL	4
#define RECENCY_EGO_POSTS_SHARED_PARTNER_DEGREE_CONSTRAINED_RELATIONAL_EVENT_NETWORK_MODEL	5

#define UNIFORM	1
#define TRIANGULAR	2
#define EPANECHNIKOV	3

typedef vector<double> ExposureTimeVectorType;
typedef pair<int, int> Edge;
typedef pair<int, unsigned int> StatValue; // the stat index and the value
typedef pair<Edge, StatValue> UpdateElement;
typedef pair<Edge, double> EdgeTime;
typedef vector<EdgeTime> EdgeTimeVectorType;
typedef list<EdgeTime> EdgeTimeList;
typedef pair<int, double> NodalActivityTime;
typedef ublas::vector<NodalActivityTime> NodalActivityTimeVector;
typedef ublas::vector<NodalActivityTimeVector> VectorOfNodalActivityTimeVectors;
typedef pair<Edge, double> DyadicActivityTime;
typedef ublas::vector<DyadicActivityTime> DyadicActivityTimeVector;
typedef ublas::vector<DyadicActivityTimeVector>
		VectorOfDyadicActivityTimeVectors;
typedef ublas::vector<unsigned short> UnsignedShortVector;
typedef ublas::matrix<unsigned short> UnsignedShortMatrix;
typedef ublas::vector<unsigned int> UnsignedIntVector;
typedef ublas::matrix<unsigned int> UnsignedIntMatrix;
typedef ublas::vector<unsigned long> UnsignedLongVector;
typedef ublas::matrix<unsigned long> UnsignedLongMatrix;
typedef std::set<int> NodeSet;
typedef ublas::vector<NodeSet> VectorOfNodeSets;
typedef std::set<Edge> EdgeSet;
typedef pair<double, double> CummulativePoint;
typedef ublas::vector<CummulativePoint> VectorOfCummulativePoints;
typedef ublas::matrix<CummulativePoint> MatrixOfCummulativePoints;
typedef ublas::vector<int> IntVector;
typedef ublas::matrix<int> IntMatrix;
typedef ublas::vector<long> LongVector;
typedef ublas::matrix<long> LongMatrix;
typedef ublas::vector<float> FloatVector;
typedef ublas::matrix<float> FloatMatrix;
typedef ublas::vector<double> DoubleVector;
typedef ublas::matrix<double> DoubleMatrix;

typedef pair<double, DoubleVector> TimeParameters;
typedef ublas::vector<TimeParameters> VectorOfTimeParameters;

typedef ublas::matrix<UnsignedShortVector> MatrixOfUnsignedShortVector;
typedef ublas::matrix<DoubleVector> MatrixOfDoubleVector;
typedef ublas::vector<UnsignedShortMatrix> VectorOfUnsignedShortMatrix;
typedef ublas::vector<DoubleMatrix> VectorOfDoubleMatrix;
typedef ublas::vector<DoubleVector> VectorOfDoubleVector;

typedef adjacency_list<vecS, vecS, undirectedS> UndirectedGraph;
typedef graph_traits<UndirectedGraph> UndirectedGraphTraits;

typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
typedef graph_traits<Graph> GraphTraits;

typedef property<vertex_name_t, std::string> VertexUsernameProperty;
typedef property<edge_weight_t, double> EdgeTimeProperty;
typedef adjacency_list<vecS, vecS, bidirectionalS, VertexUsernameProperty,
		EdgeTimeProperty> VertexNameEdgeTimeGraph;
typedef adjacency_list<vecS, vecS, bidirectionalS, no_property,
		EdgeTimeProperty> EdgeTimeGraph;

typedef ublas::generalized_vector_of_vector<unsigned short, ublas::row_major,
		ublas::vector<ublas::compressed_vector<unsigned short> > >
		SparseUnsignedShortMatrix;
typedef ublas::vector<SparseUnsignedShortMatrix>
		VectorOfSparseUnsignedShortMatrix;

typedef ublas::generalized_vector_of_vector<unsigned int, ublas::row_major,
		ublas::vector<ublas::compressed_vector<unsigned int> > >
		SparseUnsignedIntMatrix;
typedef ublas::vector<SparseUnsignedIntMatrix> VectorOfSparseUnsignedIntMatrix;

typedef ublas::generalized_vector_of_vector<long, ublas::row_major,
		ublas::vector<ublas::compressed_vector<long> > > SparseLongMatrix;
typedef ublas::vector<SparseLongMatrix> VectorOfSparseLongMatrix;

typedef ublas::vector<ublas::compressed_vector<unsigned int> >
		SparseUnsignedIntVector;

typedef ublas::compressed_matrix<unsigned char, ublas::row_major>
		SparseBoolMatrix;
typedef ublas::generalized_vector_of_vector<unsigned char, ublas::row_major,
		//ublas::vector<ublas::coordinate_vector<unsigned char> > >
		ublas::vector<ublas::mapped_vector<unsigned char> > >
		SparseGVOCVBoolMatrix;

typedef list<UpdateElement> ListOfUpdateElements;
typedef ublas::vector<ListOfUpdateElements> VectorOfListsOfUpdateElements;

typedef struct {
	unsigned int source;
	unsigned int target;
	unsigned int index;
} EdgeStatistic;
typedef list<EdgeStatistic> ListOfEdgeStatistics;
typedef pair<Edge, unsigned int> EdgeStatisticPair;
typedef map<EdgeStatisticPair, int> MapOfEdgeStatistics;

typedef list<Edge> ListOfUpdatedEdges;
typedef ublas::vector<ListOfUpdatedEdges> VectorOfListOfUpdatedEdges;

typedef list<StatValue> NodalUpdates;
typedef map<int, NodalUpdates> NodalUpdateMap;
typedef ublas::vector<NodalUpdateMap> VectorOfNodalUpdateMaps;

typedef list<StatValue> EdgeUpdates;
typedef map<Edge, EdgeUpdates> EdgeUpdateMap;
typedef ublas::vector<EdgeUpdateMap> VectorOfEdgeUpdateMaps;

typedef map<Edge, unsigned int> Edge2UnsignedIntMap;
typedef list<double> DoubleList;
typedef map<Edge, DoubleList> Edge2PostTimesMap;

typedef list<int> ListOfVertices;
typedef ublas::vector<ListOfVertices> VectorOfListsOfVertices;

typedef list<Edge> ListOfEdges;
typedef ublas::vector<ListOfEdges> VectorOfListsOfEdges;

typedef vector<int> NodeVector;
typedef pair<int, NodeVector> EdgeEvents;
typedef vector<EdgeEvents> VectorOfEdgeEvents;
typedef list<EdgeEvents> ListOfEdgeEvents;

typedef pair<double, VectorOfEdgeEvents> DiscreteTimeVectorOfEdgeEvents;
typedef vector<DiscreteTimeVectorOfEdgeEvents>
		VectorOfDiscreteTimeVectorOfEdgeEvents;
typedef list<DiscreteTimeVectorOfEdgeEvents>
		ListOfDiscreteTimeVectorOfEdgeEvents;

typedef set<int> SetOfVertices;
typedef map<int, SetOfVertices> MatchMap;

/* Matrix inversion routine.
 Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool InvertMatrix(const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
	using namespace boost::numeric::ublas;

	// create a working copy of the input
	matrix<T> A(input);

	// create a permutation matrix for the LU-factorization
	permutation_matrix<std::size_t> pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(ublas::identity_matrix<T>(A.size1()));

	// back substitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}

inline double min(double a, double b) {
	return (a < b) ? a : b;
}

inline double max(double a, double b) {
	return (a > b) ? a : b;
}

}

#endif /* DATATYPES_H_ */
