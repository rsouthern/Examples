#include <iostream>

#include "pointcloud.hpp"

typedef double real;
typedef PointCloud<real,2> MyPointCloud;

/**
 * @brief buildQuerySet
 * @param pts
 */
void buildQuerySet(MyPointCloud::PointVector &pts, size_t N, real pert) {
    srand(time(NULL));

    MyPointCloud::PointType query;
    real x;
    real part = 2.0f * boost::math::constants::pi<real>() / real(N);
    size_t i;

    for (i = 0; i < N; ++i) {
        x = real(i) * part;
        query << x, sin(x)+ pert * ((rand() % 1000)/1000.0f);
        pts.push_back(query);
    }
}

int main(int, char **) {
    MyPointCloud pc;
    pc.initRandomSine(100);
    MyPointCloud::PointVector pts;
    buildQuerySet(pts, 1000, 0.1);

    //pc.debug_dump_data();

    MyPointCloud::PointVector::iterator it;
    for (it = pts.begin(); it != pts.end(); ++it) {
        MyPointCloud::PointType result = pc.mlsProject((*it), 0.2);
        std::cout << result.transpose()<<"\n";
    }

    return EXIT_SUCCESS;
}

