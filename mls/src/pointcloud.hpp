#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/Geometry>
#include <nanoflann.hpp>
#include <cstdlib>
#include <memory>
#include <string>
#include <boost/math/constants/constants.hpp>

template <typename REAL, unsigned int DIM>
class PointCloud {
public:
    /// Define the type used to store point data
    typedef Eigen::Matrix<REAL, DIM, 1> VectorDr;
    typedef VectorDr PointType;

    /// Store our matrix in a big 'ol std::vector
    typedef std::vector<PointType> PointVector;

    /// A type to represent hyperplanes (handy for projections)
    typedef Eigen::Hyperplane<REAL,DIM> Hyperplane;

    /// A type to represent dynamic matrices
    typedef Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> MatrixXr;

private:
    /// A KDTreeIndex type
    typedef nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<REAL, PointCloud<REAL,DIM> > ,
        PointCloud<REAL,DIM>,
        DIM> KDTreeIndex;

    /// A type to store the results of a nearest neighbour query
    typedef nanoflann::RadiusResultSet<REAL,size_t> NeighbourSet;

    /// This contains the indices and distance from query point for a NN query
    typedef std::vector<std::pair<size_t,REAL> > ResultsVector;

public:
    /// Constructor
    PointCloud();

    /// Initialise with random data
    void initRandomSine(const size_t /*number of points*/);

    /// Initialise from an obj file
    void initFromObj(const std::string /*filename*/) {}

    /// Destructor
    ~PointCloud() {}

    /// Perform a moving least squares projection of the given point
    PointType mlsProject(const PointType& /*pt*/,
                         const REAL& /*radius*/,
                         const REAL& /*tol*/ = REAL(0.00001));

    /// Testy McTesterson!
    void test_neighbour_query(const PointType /*query point*/, const REAL radius);

    // *** Needed for nanoflann ***
    /// Returns the number of points in the point cloud
    inline size_t kdtree_get_point_count() const {
        return m_data.size();
    }

    // *** Needed for nanoflann ***
    /// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
    inline REAL kdtree_distance(const REAL *p1_data, const size_t idx_p2, size_t size) const {
        Eigen::Map<const PointType> p1_map(p1_data,DIM);
        return (m_data[idx_p2] - p1_map).norm();
    }

    // *** Needed for nanoflann ***
    /// Returns the dim'th component of the idx'th point in the class (random point access)
    inline REAL kdtree_get_pt(const size_t idx, int dim) const {
        return m_data[idx][dim];
    }

    // *** Needed for nanoflann ***
    /// Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const {return false;}

public:
    /// Dump the data to the screen
    void debug_dump_data();

private:
    /// A big old data store for our structure
    PointVector m_data;

    /// A KDTree structure
    KDTreeIndex m_index;

    /// Keep track on whether the class has been appropriately initialized
    bool m_init;

    /// Compute a weighted LS projection of a given point onto it's neighbours
    size_t weightedLS(const PointType & /*pt*/,
                      const REAL & /*radius*/,
                      Hyperplane &/*plane*/,
                      ResultsVector &/*index weights pair*/);

    /// A private initialisation function
    void _init();

    /// A private deletion function
    void _clear();
};

/**
 * Build an empty point cloud with an empty KDTree.
 */
template <typename REAL, unsigned int DIM>
PointCloud<REAL,DIM>::PointCloud()
    : m_index(DIM, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10)) {
}

/**
 * We assume that data has been loaded into m_data before this can be called.
 */
template <typename REAL, unsigned int DIM>
void PointCloud<REAL,DIM>::_init() {
    if (m_init) return; // sanity check
    m_index.buildIndex();
    m_init = true;
}

template <typename REAL, unsigned int DIM>
void PointCloud<REAL,DIM>::_clear() {
    m_data.clear();
    m_index.freeIndex();
    m_init = false;
}

/**
 * Initialise the point cloud from random data
 */
template <typename REAL, unsigned int DIM>
void PointCloud<REAL,DIM>::initRandomSine(size_t N) {
    srand(time(NULL));
    if (m_init) _clear();
    int i,j;
    m_data.resize(N);
    for (i = 0; i<N; ++i) {
        PointType vec;
        vec[0]= (rand() % 1000) / REAL(1000);
        for (j=1; j<DIM; ++j) {
            vec[j]= sin(j * boost::math::constants::pi<REAL>() * vec[0]);
        }
        m_data[i] = vec;
    }
    _init();
}

/**
 * Dump the contents of the data structure to the screen
 */
template <typename REAL, unsigned int DIM>
void PointCloud<REAL,DIM>::debug_dump_data() {
    if (!m_init) return;

    typename PointVector::iterator it;
    for (it = m_data.begin(); it != m_data.end(); ++it) {
        std::cout << (*it).transpose() << "\n";
    }
}

/**
 *
 */
template <typename REAL, unsigned int DIM>
void PointCloud<REAL,DIM>::test_neighbour_query(const PointType pt, const REAL radius) {
    if (!m_init) return;

    // Create a structure to store the indices and distances for neighbours
    ResultsVector idxDist;

    // The results structure will contain our results (not sure of it's relationship with idxDist?)
    NeighbourSet results(radius,idxDist);

    // Actually perform the query using the input data (which is of a raw data point type)
    std::cout << "Results of Query Test at "<<pt.transpose()<<" with radius "<<radius<<":\n";
    if (m_index.findNeighbors(results, pt.data(), nanoflann::SearchParams())) {
        // The result of the call is true if a neighbour could be found
        typename ResultsVector::iterator it;
        for (it = idxDist.begin(); it != idxDist.end(); ++it) {
            std::cout << " - ("<<m_data[(*it).first].transpose()<<"), distance "<<(*it).second<<"\n";
        }
    }
}

/**
 * Perform a moving least squares projection of pt onto the unknown MLS surface defined
 * by the points in the point cloud. This uses the weightedLS function to determine the
 * best fit hyperplane defined by the points within the given radius.
 */
template <typename REAL, unsigned int DIM>
typename PointCloud<REAL,DIM>::PointType PointCloud<REAL,DIM>::mlsProject(const PointType& pt,
                                           const REAL& radius,
                                           const REAL& tol) {
    if (!m_init) return PointType();

    // Initialise our error to be something big
    REAL err = std::numeric_limits<REAL>::max();

    // Keep our current point handy
    PointType cpt = pt;
    PointType tmp;

    // Keep our current hyperplane for projection purposes
    Hyperplane plane;

    // Our results are not used in mlsProject yet
    ResultsVector idxDist;

    // Continue this loop until we're close enough to the surface
    while (err > tol) {
        // Perform a least squares projection from the current point onto a best
        // fit hyperplane. The result will tell us if any neighbours were found.
        if (weightedLS(cpt, radius, plane, idxDist) == 0) {
            // What should we do here? throw an error?
            std::cout << "mlsProject(["<<pt.transpose()<<"],"<<radius<<") - no points within radius!\n";
            return cpt;
        } else {
            std::cout << "mlsProject(["<<pt.transpose()<<"],"<<radius<<") - "<<idxDist.size()<<" points within radius!\n";
        }

        // If it all went well, we'll have a plane onto which the point can be projected
        tmp = plane.projection(cpt);

        // Our error is just the distance between the input point and our current point
        err = (tmp-cpt).norm();

        // Update our current point to the projected point
        cpt = tmp;
    }

    // If it all went well, we can assume our current point is the best one available
    return cpt;
}

/**
 *
 */
template <typename REAL, unsigned int DIM>
size_t PointCloud<REAL,DIM>::weightedLS(const PointType &pt,
                                        const REAL &radius,
                                        Hyperplane &h,
                                        ResultsVector &idxDist) {
    if (!m_init) return 0;

    // Make sure the input results vector contains nothing
    idxDist.clear();

    // The results structure will contain our results (not sure of it's relationship with idxDist?)
    NeighbourSet results(radius,idxDist);

    if (m_index.findNeighbors(results, pt.data(), nanoflann::SearchParams())) {
        // Weights matrix
        MatrixXr W(idxDist.size(), idxDist.size()); W.setZero();
        // The matrix of the actual points
        MatrixXr K(idxDist.size(), DIM);

        // The result of the call is true if a neighbour could be found
        typename ResultsVector::iterator it;
        int i,j;
        for (it = idxDist.begin(), i=0; it != idxDist.end(); ++it,++i) {
            // First check to see if this evaluates to 0 exactly (BAD)
            if ((*it).second == REAL(0)) {
                // Set the weight to half the max REAL value (risky)
                W(i,i) = std::numeric_limits<REAL>::max() * REAL(0.5);
            } else {
                // Set the weight to the inverse distance
                W(i,i) = 1.0 / (*it).second;
            }
            // Build the matrix K out of the point data
            for (j=0;j<DIM;++j) K(i,j) = m_data[(*it).first][j];

        }
        // Normalise the weights to make sure they sum to 1
        W = W / W.sum();

        // Now solve the system using the following formula from MATLAB
        // N' = (K'*diag(W)*K) \ (K')*diag(W)*(-ones(size(I,1),1));
        // d  = 1/norm(N);
        // N  = N' * d;
        VectorDr ans;
        MatrixXr A = K.transpose() * W * K;
        VectorDr b = -K.transpose() * W.diagonal();
        ans =  A.llt().solve(b);
        REAL d = REAL(1) / ans.norm();
        Hyperplane _h(ans * d, d); // Should the scalar be negative??
        h = _h;

        // For the results structure to be useful we need it to store the weights from W, so
        // we copy these back into it
        for (it = idxDist.begin(), i=0; it != idxDist.end(); ++it,++i) {
            (*it).second = W(i,i);
        }
        // We return the number of effective matches
        return idxDist.size();
    } else {
        // No neighbours found! What do we do?
        return 0;
    }
}

#endif // POINTCLOUD_H
