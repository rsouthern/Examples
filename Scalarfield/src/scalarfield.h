#ifndef SCALARFIELD_H
#define SCALARFIELD_H

#include <stdlib.h>

// Needed to pass a function to this class when constructing the initial grid of scalars in initialise
#include <functional>

// Needed for std::array
#include <array>

#include <fstream> // file i/o
#include <string.h> // string stuff
#include <iostream> // cout
#include <algorithm> //std::copy
#include "catmulrom.h" // Catmulrom splines


// This could cause problems, as uint's are popular typedefs
typedef unsigned int uint;

/**
 * Construct a spacial scalar field in any dimension, with the goal of interpolating across the space using
 * an appropriate interpolation scheme. This particular class uses Catmull-Rom splines formed from
 * neighbouring points to the input space.
 *
 * This is general dimensional, but only interpolates a scalar field in that space. The dimension is
 * specified on construction - this could be made a template parameter I guess, but doesn't help much
 * at the moment. Obvious the return types and precision of the scalar field could be user specified.
 *
 * For portability reasons I've not depended on STL structures to hold the data - everything is
 * stored in plain-old data chunks. The n-D scalar values are stored in m_data, which is the flattened
 * array of values. The min and max values of the space is stored in m_minValue and m_maxValue, and
 * the resolution is stored in m_resolution. The partition size (useful for internal computations,
 * calculated from the previous values) is stored in m_partSize and m_invPartSize (the inverse).
 *
 * The general dimensional functionality is achieved using recursion. There are two functions: initRecurse
 * and evalRecurse, which use this approach. The main API for this class is the init() with
 * user specified function, max and min grid ranges in each dimension, and the resolution, which
 * initialises the grid using the scalar function (calling again will clear existing data). Also the
 * eval() function is the key interface, which recursively evaluates the interpolated scalar value
 * using catmull-rom splines.
 *
 * Note that at some stage in the future it would be great to extend this as a base class and offer
 * alternative interpolation techniques, such as RBF's, for comparison. Also, note that it is limited
 * in that the grid is regular (required by Catmull-Rom splines).
 */
template <uint DIM>
class ScalarField {
public:
    /// Unsigned int array
    typedef std::array<uint, DIM> uintArray;

    /// Double array
    typedef std::array<double, DIM> doubleArray;

public:
    /// Construct the empty scalar field with the specified dimension
    ScalarField() : m_data(nullptr), m_isInit(false) {;}

    /// Destroy the scalar field
    ~ScalarField() {cleanup();}

    /// Initialise the grid
    void init(std::function<double (double*)> _func,
              const doubleArray & _minValue,
              const doubleArray & _maxValue,
              const uintArray & _resolution);

    /// Retrieve the value at the specified point location using a preferred interpolation technique
    double eval(const doubleArray & pt) const;

    /// For debugging, dump the entire structure as an obj file (only works if field is 2D)
    void dumpObj() const;

    /// Resample the scalar field using the catmul rom goodies
    void resampleObj(const uintArray &_resolution) const;

private:
    /// Recursively build the n-Dimensional grid using the specified function
    void initRecurse(std::function<double (double*)> _func,
                     const uint &dim,
                     uintArray &idx,
                     doubleArray &pos);

    /// Recursively evaluate the n-Dimensional spline using catmull rom spline interpolation
    double evalRecurse(uint dim,
                       const uintArray & idx,
                       uintArray &currentIdx,
                       doubleArray & x) const ;

    /// Delete memory associated with this structure
    void cleanup();

    /// Store the data in a fat array
    double *m_data;

    /// Keep track of the min, max and resolution values
    doubleArray m_minValue, m_maxValue, m_partSize, m_invPartSize;
    uintArray m_resolution;

    /// Set to true when this class is ready to be used
    bool m_isInit;
};


/**
 * @brief ScalarField::dumpObj
 */
template<uint DIM>
void ScalarField<DIM>::dumpObj() const {
    std::cout << "ScalarField::dumpObj() not implemented for dimension "<<DIM<<"\n";
}

template <>
void ScalarField<2>::dumpObj() const {
    if (!m_isInit) return;
    uint i,j;
    for (i=0; i<m_resolution[0]; ++i) {
        for (j=0; j<m_resolution[1]; ++j) {
            std::cout  << "v "<<m_minValue[0] + i*m_partSize[0]
                    <<" "  <<m_minValue[1] + j*m_partSize[1]
                    <<" "  <<m_data[i * m_resolution[0] + j] << "\n";
        }
    }
}

template<uint DIM>
void ScalarField<DIM>::resampleObj(const uintArray &_resolution) const {
    std::cout << "ScalarField::resampleObj() not implement for dimension "<<DIM<<"\n";
}

template<>
void ScalarField<2>::resampleObj(const uintArray &_resolution) const {
    if (!m_isInit) return;
    uint i,j;
    doubleArray _partSize;
    doubleArray _pos;

    for (i=0; i<2; ++i) _partSize[i] = (m_maxValue[i]-m_minValue[i]) / double(_resolution[i]);

    for (i=0; i<_resolution[0]; ++i) {
        _pos[0] = m_minValue[0] + i*_partSize[0];
        for (j=0; j<_resolution[1]; ++j) {
            _pos[1] = m_minValue[1] + j*_partSize[1];
            std::cout << "v "<< _pos[0] << " " << _pos[1] << " " << eval(_pos) << "\n";
        }
    }
}


/**
 * @brief ScalarField::eval Evaluate the scalar value at the specified point.
 * @param pt The position at which to evaluate the scalar field
 * @return
 */
template<uint DIM>
double ScalarField<DIM>::eval(const doubleArray& pt) const {
    if (!m_isInit) return 0.0;

    // Determine the index of the bottom corner of the evaluation grid (if outside, extrapolation
    // may be used, so it's capped at the bounds of the grid)
    uintArray idx;
    uintArray currentIdx;
    doubleArray x;

    uint i;
    // Deduce the index of the points from the min, max grid values and total resolution
    for (i=0; i<DIM; ++i) {
        // Calculate the index accounting for bounds
        int _idx = int((pt[i] - (m_minValue[i]+m_partSize[i])) * m_invPartSize[i]);
        if (_idx >= (int(m_resolution[i])-3)) {
            idx[i] = int(m_resolution[i]) - 4;
        } else if (_idx < 0) {
            idx[i] = 0;
        } else {
            idx[i] = _idx;
        }
        // Now compute the interpolation value (this may be negative)

        x[i] = (pt[i] - (m_minValue[i] + (double(idx[i]+1) * m_partSize[i]))) * m_invPartSize[i];
    }

    // Evaluate the position recursively over all dimensions
    double ret_val = evalRecurse(0, idx, currentIdx, x);

    // Clear memory and return
    return ret_val;
}

/**
 * @brief ScalarField::evalRecurse This is a rather elegant recursive method for evaluating
 * n-D scalar catmull-rom splines. It recursively evaluates the spline in each dimension, combining
 * all the results in each pass.
 * @param dim The current dimension (for initial execution, use 0)
 * @param idx The initial index (corner of the subgrid from which the value is interpolated)
 * @param currentIdx The index of this particular evaluation
 * @param x The parameter value to use for interpolation
 * @return
 */
template<uint DIM>
double ScalarField<DIM>::evalRecurse(uint dim,
                                     const uintArray & idx,
                                     uintArray & currentIdx,
                                     doubleArray & x) const  {
    // The termination condition - just return the scalar value at the specified position
    if (dim >= DIM) {

        // First determine the index in the data
        uint dataIdx = currentIdx[DIM-1];
        uint dimSize = m_resolution[DIM-1];
        int i;
        for (i = DIM-2; i >= 0; --i) {
            dataIdx += currentIdx[i] * dimSize;
            dimSize *= m_resolution[i];
        }
        // Just return the data value at this position!
        return m_data[dataIdx];
    } else {
        // Evaluate the catmull rom spline recursively by combining the current indices
        // along all dimensions
        double v0, v1, v2, v3;
        currentIdx[dim] = idx[dim] + 0; v0 = evalRecurse(dim+1, idx, currentIdx, x);
        currentIdx[dim] = idx[dim] + 1; v1 = evalRecurse(dim+1, idx, currentIdx, x);
        currentIdx[dim] = idx[dim] + 2; v2 = evalRecurse(dim+1, idx, currentIdx, x);
        currentIdx[dim] = idx[dim] + 3; v3 = evalRecurse(dim+1, idx, currentIdx, x);
        return catmullRomSpline(x[dim], v0, v1, v2, v3);
    }
}

/**
 * @brief ScalarField::initialise Inits the scalar field
 * @param _func A boost-ified function declaration which returns the scalar value
 * @param _minValue The minimum values of our given field.
 * @param _maxValue The maximum values of our given field.
 * @param _resolution The resolution of the field in each dimension.
 */
template<uint DIM>
void ScalarField<DIM>::init(std::function<double (double*)> _func,
                            const doubleArray &_minValue,
                            const doubleArray &_maxValue,
                            const uintArray& _resolution) {
    // Check we've not been here before already
    if (m_isInit) cleanup();

    // Copy over the bounds (should be able to std move this?)
    std::copy(_minValue.cbegin(), _minValue.cend(), m_minValue.begin());
    std::copy(_maxValue.cbegin(), _maxValue.cend(), m_maxValue.begin());
    std::copy(_resolution.cbegin(), _resolution.cend(), m_resolution.begin());

    // Precompute the partition size in each case
    int i;
    for (i=0; i<DIM; ++i) {
        m_partSize[i] = (m_maxValue[i] - m_minValue[i]) / double(m_resolution[i] - 1);
        m_invPartSize[i] = 1.0 / m_partSize[i];
    }

    // Now we can allocate some memory for our data
    uint datasize = _resolution[0];
    for (i=1; i<DIM; ++i) datasize *= _resolution[i];
    m_data = (double*) malloc(sizeof(double) * datasize);

    // Now for each value in the field we can use the provided function to derive the scalar value
    // this is done recursively (could easily be done in parallel
    doubleArray pos;
    uintArray idx;
    initRecurse(_func, 0u, idx, pos);

    // Set this to initialised
    m_isInit = true;
}

/**
 * @brief ScalarField::initRecurse Recursively initialise the scalar field along each dimension.
 * @param _func The boostified function to execute to determine the scalar value
 * @param dim The current dimension (start this at 0)
 * @param idx The current index to evaluate (empty when you start)
 * @param pos The position to evaluate (empty when you start)
 */
template<uint DIM>
void ScalarField<DIM>::initRecurse(std::function<double (double*)> _func,
                                   const uint &dim,
                                   uintArray &idx,
                                   doubleArray &pos) {
    int i;
    // At the end of the recursion chain we can calculate the value using the function
    if (dim >= DIM) {
        // First determine the index in the data
        uint dataIdx = idx[DIM-1];
        uint dimSize = m_resolution[DIM-1];

        for (i = DIM-2; i >= 0; --i) {
            dataIdx += idx[i] * dimSize;
            dimSize *= m_resolution[i];
        }

        // Now set the value of the data in the right place to the result of the function
        m_data[dataIdx] = _func(&pos[0]);

        return;
    }

    // If we didn't get to the end of the chain, we run another iteration
    for (i = 0 ; i < m_resolution[dim]; ++i) {
        idx[dim] = i;
        //std::cout << "ScalarField::initRecurse(): idx["<<dim<<"]="<<idx[dim]<<"\n";
        pos[dim] = m_minValue[dim] + i * m_partSize[dim];
        initRecurse(_func, dim+1, idx, pos);
    }
}

/**
 * @brief ScalarField::cleanup
 */
template<uint DIM>
void ScalarField<DIM>::cleanup() {
    if (m_isInit) {
        delete [] m_data;
        m_isInit = false;
    }
}

#endif // SCALARFIELD_H
