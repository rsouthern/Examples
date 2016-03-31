#include "scalarfield.hpp"
#include <iostream>
#include <noise/src/noise.h>

// Silly test function
double func2(double *x) {
    return 0.05 * sin(20.0 * x[0]+ 20.0 * x[1]);
}

// Silly test function
double func3(double *x) {
    return 0.05 * sin(20.0 * x[0] + 20.0 * x[1] + 20.0 * x[2]);
}

noise::module::Perlin perlin;

double perlin_noise_2D(double *x) {
    return 0.2 * perlin.GetValue(x[0], x[1], 0.0);
}

typedef ScalarField<2u> ScalarField2;
typedef ScalarField<3u> ScalarField3;


int main() {

    ScalarField2 sfield;
    ScalarField2::doubleArray min, max;
    ScalarField2::uintArray res, res2;
    min[0] = min[1] = 0.0;
    max[0] = max[1] = 1.0;
    res[0] = res[1] = 50;
    res2[0]=res2[1]=500;

    sfield.init(perlin_noise_2D, min, max, res);
    sfield.resampleObj(res2);
//    std::cout << "Evaluating scalarfield "<< sfield.eval(max)<<"\n";

}
