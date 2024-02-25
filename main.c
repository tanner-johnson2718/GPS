// * Use GSL for all data types and math operations
// * All vector data types will be a gsl matrix
// * All (non-trivial) functions return an error type
// * Based on book GPS Theory, Algorithms, and applications

#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <math.h>

//*****************************************************************************
// Util
//*****************************************************************************

#define SPATIAL_DIM 3
#define LAT_MIN -M_PI_2
#define LAT_MAX M_PI_2
#define LON_MIN -M_PI
#define LON_MAX M_PI
#define A_EARTH_METERS 6378137.0
#define B_EARTH_METERS 6356752.0
#define E_EARTH 0.08181979099
#define TAN_EPS .0001
#define ITERATIVE_EPS .00000001
#define ITERATIVE_MAX_ITERS 1000

typedef enum
{
    OK,
    INVALID_DIM,
    ITERATIVE_SOLN_FAILED
} err_t;

static err_t print_matrix(gsl_matrix* m)
{
    size_t i, j;
    for(i = 0; i < m->size1; ++i)
    {
        for(j = 0; j < m->size2; ++j)
        {
            printf("%0.4f ", gsl_matrix_get(m, i, j));
        }
        printf("\n");
    }

    return OK;
}

//*****************************************************************************
// ECEF and Ellipsoidial cordinates. ECEF or Earth Centered, Earth Fixed is a 
// frame a measurement as it sounds expressed in cartesian coordinates. The z
// axis extends through the poles with North being pos and origin at center of 
// earth. The X-plane is zero-ed at the Greenwhich median. Y is perdicualar
// east ward. We map x,y,z to our gsl_mat indicies as such:
//
//    gsl_mat index 0 -> X
//    gsl_mat index 1 -> Y
//    gsl_mat index 2 -> Z
//
// Spherical. Since Eart is not a perfect sphere, we generally do no map
// cartesian to sperical. As a prelude to ellipsoidal cords, spherical 
// cordinates can be expressed as)
//
//    x = r * cos(theta) * cos(lamdba) 
//    y = r * cos(theta) * sin(lambda) 
//    z = r * sin(theta)
//
// Theta is the angle between the x/y plane and line from your current pos to
// the origin. R is the length of that line and lambda is angle of your pos
// projected onto the x/y plane and the x-axis. To derive this yourself, take
// a sphere, keep Z facing up and take a circular cut at an arbitrary angle
// lamda. Now you can see cos(theta) projects you onto the the X/Y plane. From
// there you have a circular projection of radius r which can be parameterized 
// in lamda the usual way i.e. (rcos(lamdba), rsin(lambda)). We can also see from
// here that Z is just the verical projection of the line shown or r*sin(thata)
//
//                                 Z
//                                 |
//                          ooo OOO|OOO ooo  *
//                      oOO        |        /
//                  oOO            |       /    OOo
//               oOO               |      /        OOo
//             oOO                 |     /           OOo
//           oOO                   |    /              OOo
//          oOO                    |   /                OOo
//         oOO                     |  /                  OOo
//        oOO                      | /                    OOo
//        oOO                      |/  theta              OOo
//   ------------------------------|--------------------------------- Lambda
//        oOO                      |                      OOo
//        oOO                      |                      OOo
//         oOO                     |                     OOo
//          oOO                    |                    OOo
//           oOO                   |                   OOo
//             oOO                 |                 OOo
//               oO                |               OOo
//                  oOO            |            OOo
//                      oOO        |        OOo
//                          ooo OOO|OOO ooo
//                                 |
//       
//
// Ellipsoidal. Now we can extend this type of reasoning to an ellipsoid which
// is an ellipse rotated about the Z-axis showb below. We denote the vertical
// axis length of thw ellipse as b (this is the smaller i.e. minor axis). We
// denote the major or horizontal axis as a. Throughout the code and our system
// we will use a and the parameter e, the eccentricty, which is defined as shown 
// below and gives a positive such that the larger the number the more "oblong"
// the ellipsoid is with e=0 being the case of sphere.
//
//    e = sqrt(a*a - b*b) / a
//
//
//
//                                               *   P(h, lambda, theta)
//                                 Z            /
//                                 |         h /  
//                          ooo OOO|OOO ooo   /  
//                      oOO        |         *
//                  oOO            |        /   OOo
//               oOO               |       /       OOo
//             oOO                 |      /          OOo
//           oOO                   |   N /            OOo
//          oOO                    |    /               OOo
//         oOO                     |   / theta           OOo
//   ------------------------------|--/------------------------------ Lambda
//         oOO                     | /                   OOo
//          oOO                    |/                   OOo
//           oOO                   *                   OOo
//             oOO                 |                 OOo
//               oO                |               OOo
//                  oOO            |            OOo
//                      oOO        |        OOo
//                          ooo OOO|OOO ooo
//                                 |
//
// As you can see we need to generalize the radial line, which in ellipsoidail
// coordinates is given by N
//
//    N - is the prime normal vertical i.e. I am sitting on the surface of an 
//        ellisoid. I draw a tangent. The normal to the tangent plane extented
//        from my location to the Z-axis. The length of this line is N.
//
// We also replace the radial distance with simply a height above the surface,
// giving an above sea level altitude. The formal coordinate transformation is 
// given by the following:
//
//    x = (N+h)         * cos(theta) * cos(lamdba) 
//    y = (N+h)         * cos(theta) * sin(lambda) 
//    z = (N(1-e^2) +h) * sin(theta)
//    N = a / sqrt(1 - (e^2)(sin^2 (theta)))
//
// Deriving x and y and is a simple extension of the logic used in the
// spherical case. To compute both N and Z one needs to do the following.
//
//    1) Take an ellipse slice at some angle lamba
//    2) Parameterize ellipse: {x = acos(theta), y = b sin(theta)}
//    3) Compute Tangent Line slope dy/dx
//    4) Compute normal line to tangent as y = (-dx/dy)x + A
//         4a) A is the y intercept
//    5) N is distance between point {acos(theta), b sin(theta)} and y intercept
//    6) Z is N - the portion of the normal verical below the y axis. 
//
// We use the following mapping between ellipsoidal cords and gsl mat indicies
//
//    gsl_mat index 0 -> h (height above ellipse)
//    gsl_mat index 1 -> lamda or X/Y plane angle or longitutde
//    gsl_mat index 2 -> theta or Angle made with the Z plane
//
//*****************************************************************************

static inline double N(double a, double e, double theta)
{
    return  a / sqrt(1.0 - (e*e*sin(theta) * sin(theta)));
}

static inline double deg2rad(double d)
{
    return ((d * M_PI) / 180.0);
}


// in place conversion 
static err_t ellipsodial2cartesian(gsl_matrix* m, double a, double e)
{
    if(m->size1 != SPATIAL_DIM || m->size2 != 1)
    {
        return INVALID_DIM;
    }
    
    double h = gsl_matrix_get(m, 0,0);
    double lambda = gsl_matrix_get(m, 1,0);
    double theta = gsl_matrix_get(m, 2, 0);
    

    double _N = N(a, e, theta);
    // printf("N = %0.4f\n", _N);

    double x = (_N+h) * cos(theta) * cos(lambda);
    double y = (_N+h) * cos(theta) * sin(lambda);
    double z = ((_N*(1.0 - (e*e))) + h ) * sin(theta);

    gsl_matrix_set(m, 0,0, x);
    gsl_matrix_set(m, 1,0, y);
    gsl_matrix_set(m, 2,0, z); 
 
    return OK;
}

// theta iteratice func s.t. roots are soln
typedef struct
{
    double a;
    double e;
    double c0;
    double z;
} theta_params_t;
static double theta_roots(double theta, void* params)
{
    theta_params_t* p = (theta_params_t*) params; 
    double a = p->a;
    double e = p->e;
    double c0 = p->c0;
    double z = p->z;
    double _N = N(a, e, theta);

    return (z / (c0 - ((e*e)*_N*cos(theta)) )) - tan(theta);
}

// Requires solving iteratively as ellipsodial cordinates cannot be expressed
// in a closed form soln of x,y,z.
static err_t cartesian2ellipsodial(gsl_matrix* m, double a, double e)
{
    if(m->size1 != SPATIAL_DIM || m->size2 != 1)
    {
        return INVALID_DIM;
    }
    
    double x = gsl_matrix_get(m,0,0);
    double y = gsl_matrix_get(m,1,0);
    double z = gsl_matrix_get(m,2,0);

    double lambda;
    double theta;
    double h;

    if(fabs(x) < TAN_EPS )
    {
        lambda = x;
    }
    else
    {
        lambda = atan2(y, x);
    }

    double c0 = sqrt( (x*x) + (y*y) );
    theta_params_t p;
    p.a = a;
    p.c0 = c0;
    p.e = e;
    p.z = z;

    const gsl_root_fsolver_type * T = gsl_root_fsolver_bisection;
    gsl_root_fsolver * s = gsl_root_fsolver_alloc (T);
    gsl_function F;
    F.function = &theta_roots;
    F.params = &p;
    gsl_root_fsolver_set(s, &F, LAT_MIN, LAT_MAX);

    int i = 0;
    double prev;
    double curr;

    gsl_root_fsolver_iterate(s);
    curr = gsl_root_fsolver_root(s);
    ++i;
    
    for(i = 1; i < ITERATIVE_MAX_ITERS; ++i)
    {
        if(gsl_root_fsolver_iterate(s))
        {
            i = ITERATIVE_MAX_ITERS;
            break;
        }

        prev = curr;
        curr = gsl_root_fsolver_root(s);
        // printf("Iter = %d   Curr = %0.4f   Delta = %0.4f\n", i, curr, fabs(prev-curr));

        if(fabs(prev-curr) < ITERATIVE_EPS)
        {
            break;
        }

    }

    gsl_root_fsolver_free(s);

    if(i == ITERATIVE_MAX_ITERS)
    {
        printf("WARING!!\n");
        return ITERATIVE_SOLN_FAILED;
    }

    theta = curr;
    h = (c0 / cos(theta)) - N(a, e, theta);

    gsl_matrix_set(m, 0,0, h);
    gsl_matrix_set(m, 1, 0, lambda);
    gsl_matrix_set(m, 2, 0, theta);
    
    return OK;
}

//*****************************************************************************
// Local Coordinate Transformation
//*****************************************************************************

//*****************************************************************************
// Unit Tests
//*****************************************************************************

double l1_error(gsl_matrix* m1, gsl_matrix* m2)
{
    double e1 = fabs(gsl_matrix_get(m1, 0,0) - gsl_matrix_get(m2, 0,0));
    double e2 = fabs(gsl_matrix_get(m1, 1,0) - gsl_matrix_get(m2, 1,0));
    double e3 = fabs(gsl_matrix_get(m1, 2,0) - gsl_matrix_get(m2, 2,0));
    return (e1 + e2 + e3);
}

double test_l1_error(double h, double lon, double lat)
{
    gsl_matrix* m1 = gsl_matrix_calloc(3,1);
    gsl_matrix* m2 = gsl_matrix_calloc(3,1);
    gsl_matrix_set(m1, 0, 0, h);
    gsl_matrix_set(m1, 1, 0, deg2rad(lon));
    gsl_matrix_set(m1, 2, 0, deg2rad(lat));
    gsl_matrix_set(m2, 0, 0, h);
    gsl_matrix_set(m2, 1, 0, deg2rad(lon));
    gsl_matrix_set(m2, 2, 0, deg2rad(lat));

    ellipsodial2cartesian(m1, A_EARTH_METERS, E_EARTH);
    cartesian2ellipsodial(m1, A_EARTH_METERS, E_EARTH);
    
    double error = l1_error(m1, m2);
    gsl_matrix_free(m1);
    gsl_matrix_free(m2); 

    return error;
}

void ellipsoid_test(char* prompt, double h, double lon, double lat)
{
    double e = test_l1_error(h,lon,lat);
    printf("%-10s (%0.4f,%0.4f,%0.4f) |%0.4f|\n", prompt, h, lon, lat, e);
}

int main(int argc, char** argv)
{
    ellipsoid_test("Nominal", 2041.5504, deg2rad(-104.80121891303988), deg2rad(38.996328766277756));
    ellipsoid_test("Median 0", 1220, 0, M_PI/4.0);
    ellipsoid_test("Median 1", 1440, 0, 0);
    ellipsoid_test("Neg h", -100, M_PI/2.0, M_PI/4.0);
    ellipsoid_test("~Median", 1000, M_PI, 0.2);
    ellipsoid_test("~Median", 1000, -M_PI, 0.2);
    ellipsoid_test("N Pole", 0, 0, M_PI/2.0);
    ellipsoid_test("S Pole", 0, 0, -M_PI/2.0);
    ellipsoid_test("Large h", 20200000, 0,0);
    ellipsoid_test("Large h", 20200000, -1.3,0.7);
    ellipsoid_test("Large h", 20200000, M_PI,0);
    ellipsoid_test("Large h", 20200000, -M_PI,0);
    ellipsoid_test("Large h", 20200000, 0,M_PI/2.0);
    ellipsoid_test("Large h", 20200000, 0,-M_PI/2.0);
    ellipsoid_test("Large h", 20200000, M_PI,M_PI/2.0);
    return 0;
}