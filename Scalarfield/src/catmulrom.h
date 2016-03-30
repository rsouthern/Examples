#ifndef CATMULROM_H
#define CATMULROM_H

/**
  * The interpolation values (where is the best place to leave these?)
  */
static double M[16] = {
    0.0, 1.0, 0.0, 0.0,
    -0.5, 0.0, 0.5, 0.0,
    1.0,-2.5, 2.0,-0.5,
    -0.5, 1.5,-1.5, 0.5};

/**
 * @brief catmullRomSpline
 * @param x Interpolation parameter, preferably between 0 and 1, but extrapolation is theoretically possible
 * @param v0 scalar value
 * @param v1 scalar value
 * @param v2 scalar value
 * @param v3 scalar value
 * @return The interpolated value
 */
inline double catmullRomSpline(const double x,
                        const double v0,
                        const double v1,
                        const double v2,
                        const double v3) {
    double c1,c2,c3,c4;
    c1 =               M[0*4+1]*v1;
    c2 = M[1*4+0]*v0               + M[1*4+2]*v2;
    c3 = M[2*4+0]*v0 + M[2*4+1]*v1 + M[2*4+2]*v2 + M[2*4+3]*v3;
    c4 = M[3*4+0]*v0 + M[3*4+1]*v1 + M[3*4+2]*v2 + M[3*4+3]*v3;
    return (((c4 * x + c3)*x + c2)*x + c1);
}

#endif // CATMULROM_H
