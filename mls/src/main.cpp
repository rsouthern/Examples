#include <iostream>

#include "pointcloud.hpp"

typedef PointCloud<float,2> MyPointCloud;

/**
 * @brief buildQuerySet
 * @param pts
 */
void buildQuerySet(MyPointCloud::PointVector &pts, size_t N, float pert) {
    srand(time(NULL));

    MyPointCloud::PointType query;
    float x;
    float part = 2.0f * boost::math::constants::pi<float>() / float(N);
    size_t i;

    for (i = 0; i < N; ++i) {
        x = float(i) * part;
        query << x, sin(x)+ pert * ((rand() % 1000)/1000.0f);
        pts.push_back(query);
    }
}

int main(int, char **) {
    MyPointCloud pc;
    pc.initRandomSine(100);
    MyPointCloud::PointVector pts;
    buildQuerySet(pts, 1000, 0.01);

    MyPointCloud::PointVector::iterator it;
    for (it = pts.begin(); it != pts.end(); ++it) {
        //MyPointCloud::PointType result = pc.mlsProject((*it), 0.4);
        std::cout << (*it).transpose()<<"\n";
    }

    return EXIT_SUCCESS;
}

