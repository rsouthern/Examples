#include <iostream>

#include "pointcloud.hpp"

int main(int, char **) {
    PointCloud<float,2> pc;
    pc.initRandomSine(100);

    PointCloud<float,2>::PointType query; query << 0.5, 0.5;
    PointCloud<float,2>::PointType result = pc.mlsProject(query, 0.4);

    std::cout << "Query point "<<query.transpose()<<" projected result "<<result.transpose()<<"\n";

    return EXIT_SUCCESS;
}

