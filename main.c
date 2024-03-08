// * Math Based on book GPS Theory, Algorithms, and applications
// * leaflet mapping JS front end
// * wsServer C websocket library as communication
// * Use int32_t / uint16_t where appropriate
// * Use LOG instead of printf
// * Nontrivial functions shall return bool indicating success / failures.
//     * On failure nontrivial functions shall log error details
// * in ws.h change max clients to 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>
#include <time.h>
#include <pthread.h>
#include <stdint.h>
#include <unistd.h>
#include "ws.h"

#define ASYNC_TEST 0
#define GEOCORD_TEST 0
#define PATH_DB_TEST 0
#define WS_SERVER_IP "127.0.0.1"
#define WS_SERVER_PORT 6969
#define WS_RESPONSE_MSG_SIZE 128
#define NUM_SVS 30
#define MAX_SV_PATH_SIZE 32
#define AVG_SV_HEIGHT 20200000

//*****************************************************************************
// Util
//*****************************************************************************

#define LOG(tag, format, ...) _log(tag, __LINE__, __func__, format, ##__VA_ARGS__)

static void _log(char* tag, int line_no, const char* func, char* format, ...)
{
    time_t now;
    time(&now);
    va_list args;
    va_start(args, format);

    printf("[ %5s ][ %20s:%04d ] ", tag, func, line_no);
    vprintf(format, args);
    printf("\n");

    va_end(args);
}

//*****************************************************************************
// ECEF and Ellipsoidial cordinates. ECEF or Earth Centered, Earth Fixed is a 
// frame a measurement as it sounds expressed in cartesian coordinates. The z
// axis extends through the poles with North being pos and origin at center of 
// earth. The X-plane is zero-ed at the Greenwhich median. Y is perdicualar
// east ward. We map x,y,z to our gsl_mat indicies as such:
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
// As one can see going from ECEF -> WGS is much harder and requires numerics
// or approximation. Included in the repo is a paper on the "improved Zhus"
// algorithm for doing this fast and accurately. We will not go into depth
// about it as it really isn't super important.
//
//*****************************************************************************

static const double  a = +6.37813700000000000000e+0006; /* a */
static const double  invaa = +2.45817225764733181057e-0014; /* 1/(a^2) */
static const double  ee = +6.69437999014131705734e-0003; /* e^2 */
static const double  l = +3.34718999507065852867e-0003; /* (e^2)/2 */
static const double  ll4 = +4.48147234524044602618e-0005; /* e^4 */
static const double  ll = +1.12036808631011150655e-0005; /* (e^4)/4 */
static const double  p1mee = +9.93305620009858682943e-0001; /* 1-(e^2) */
static const double  p1meedaa = +2.44171631847341700642e-0014; /* (1-(e^2))/(a^2) */
static const double  Hmin = +2.25010182030430273673e-0014; /* (e^12)/4 */
static const double  invcbrt2 = +7.93700525984099737380e-0001; /* 1/(2^(1/3)) */
static const double  inv3 = +3.33333333333333333333e-0001; /* 1/3 */
static const double  inv6 = +1.66666666666666666667e-0001; /* 1/6 */
static const double  d2r = +1.74532925199432957691e-0002; /* pi/180 */
static const double  r2d = +5.72957795130823208766e+0001; /* 180/pi */

typedef struct
{
    double x;
    double y;
    double z;
} ECEF_t;

typedef struct
{
    double lat;
    double lon;
    double alt;
} WGS84_t;


static inline double N(double theta)
{
    return  a / sqrt(1.0 - (ee*sin(theta) * sin(theta)));
}

double rad(double deg)
{
    return deg * d2r;
}

double deg(double rad)
{
    return rad * r2d;
}

ECEF_t ECEF(WGS84_t in)
{
    ECEF_t p;
    double _N = N(in.lat);

    p.x = (_N+in.alt) * cos(in.lat) * cos(in.lon);
    p.y = (_N+in.alt) * cos(in.lat) * sin(in.lon);
    p.z = ((_N*(1.0 - (ee))) + in.alt ) * sin(in.lat);

    return p;
}

// Zhu's algorith (see docs in repo)
WGS84_t WGS84(ECEF_t in)
{
    double x, y, z;
    double lat, lon, alt;
    // The variables below correspond to symbols used in the paper
    // "Accurate Conversion of Earth-Centered, Earth-Fixed Coordinates
    // to Geodetic Coordinates"
    double beta;
    double C;
    double dFdt;
    double dt;
    double dw;
    double dz;
    double F;
    double G;
    double H;
    double i;
    double k;
    double m;
    double n;
    double p;
    double P;
    double t;
    double u;
    double v;
    double w;
    // Intermediate variables
    double j;
    double ww; // w^2
    double mpn; // m+n
    double g;
    double tt; // t^2
    double ttt; // t^3
    double tttt; // t^4
    double zu; // z * u
    double wv; // w * v
    double invuv; // 1 / (u * v)
    double da;
    double t1, t2, t3, t4, t5, t6, t7;
    x = in.x;
    y = in.y;
    z = in.z;
    ww = x * x + y * y;
    m = ww * invaa;
    n = z * z * p1meedaa;
    mpn = m + n;
    p = inv6 * (mpn - ll4);
    G = m * n * ll;
    H = 2 * p * p * p + G;
    if (H < Hmin)
    {
    return (WGS84_t) {0,0,0};
    }
    C = pow(H + G + 2 * sqrt(H * G), inv3) * invcbrt2;
    i = -ll - 0.5 * mpn;
    P = p * p;
    beta = inv3 * i - C - P / C;
    k = ll * (ll - mpn);
    // Compute left part of t
    t1 = beta * beta - k;
    t2 = sqrt(t1);
    t3 = t2 - 0.5 * (beta + i);
    t4 = sqrt(t3);
    // Compute right part of t
    t5 = 0.5 * (beta - i);
    // t5 may accidentally drop just below zero due to numeric turbulence
    // This only occurs at latitudes close to +- 45.3 degrees
    t5 = fabs(t5);
    t6 = sqrt(t5);
    t7 = (m < n) ? t6 : -t6;
    // Add left and right parts
    t = t4 + t7;
    // Use Newton-Raphson's method to compute t correction
    j = l * (m - n);
    g = 2 * j;
    tt = t * t;
    ttt = tt * t;
    tttt = tt * tt;
    F = tttt + 2 * i * tt + g * t + k;
    dFdt = 4 * ttt + 4 * i * t + g;
    dt = -F / dFdt;
    // compute latitude (range -PI/2..PI/2)
    u = t + dt + l;
    v = t + dt - l;
    w = sqrt(ww);
    zu = z * u;
    wv = w * v;
    lat = atan2(zu, wv);
    // compute altitude
    invuv = 1 / (u * v);
    dw = w - wv * invuv;
    dz = z - zu * p1mee * invuv;
    da = sqrt(dw * dw + dz * dz);
    alt = (u < 1) ? -da : da;
    // compute longitude (range -PI..PI)
    lon = atan2(y, x);
    return (WGS84_t) {lat, lon, alt};
}

#if GEOCORD_TEST

static void ellipsoid_test(char* prompt, double h, double lon, double lat)
{
    WGS84_t w1 = {lat, lon, h};
    WGS84_t w2 = WGS84(ECEF(w1));
    printf("%-10s W1 = (%012.4f,%012.4f,%012.4f) W2=(%012.4f,%012.4f,%012.4f)\n", prompt, w1.lat, w1.lon, w1.alt, w2.lat, w2.lon, w2.alt);
}

static void run_ellipsoid_test()
{
    ellipsoid_test("Nominal", 2041.5504, rad(-104.80121891303988), rad(38.996328766277756));
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
}

#endif

//*****************************************************************************
// Async Queue / Runners. Primitives:
//    * Q - This is the base data structure. Currently operates as a stack
//          with controled access. Posters must wait the poster cv if the q is
//          full. This cv is signaled when a comsumer or getter pops an
//          element off the queue. Elements of the queue are of the same length
//          and *should* be of the same type. Access to the queue both post
//          and get only happens atomically with the lock held.
//
//    * Post - Copies an element from the passed buffer on the queue. If full
//             wait for poster cv signal. On completion signal the getter cv.
//             Is reentrant i.e. can have multiple threads posting at the same
//             time
//
//    * Get - Copies an element form the head of the queue to the passed buffer
//            When empty will block until it recieves a getter cv signal. Upon
//            successful get, will signal the poster cv. Is reentrant.
//
//    * Runners - 
//*****************************************************************************

typedef struct
{
    pthread_cond_t poster_cv;
    pthread_cond_t getter_cv;
    pthread_mutex_t lock;
    uint32_t element_size;
    uint32_t n_max;
    uint32_t n;
    void* q;
} async_q_t;

typedef void (*async_consumer_t)(void*);

typedef struct
{
    pthread_t pthread;
    uint32_t thread_num;
    uint16_t n_runners;
    async_q_t* q;
    async_consumer_t consume;
    bool killed;
} async_runner_t;

static bool async_q_create(async_q_t* q, uint32_t element_size, uint32_t n_elements)
{
    q->poster_cv = (pthread_cond_t) PTHREAD_COND_INITIALIZER;
    q->getter_cv = (pthread_cond_t) PTHREAD_COND_INITIALIZER;
    q->lock = (pthread_mutex_t) PTHREAD_MUTEX_INITIALIZER;
    q->element_size = element_size;
    q->n_max = n_elements;  
    q->n = 0;
    q->q = calloc(element_size, n_elements);
    if(!(q->q))
    {
        LOG("ERROR", "In async_q_create calloc failure");
        return false;
    }

    return true;
}

// Unsafe. Make sure all threads waiting on the queue are joined and killed
static bool aysnc_q_destroy(async_q_t* q)
{
    bool ret = true;
    if(pthread_mutex_destroy(&(q->lock)))
    {
        LOG("ERROR", "In async_q_destroy failed to destroy lock");
        ret = false;
    }

    if(pthread_cond_destroy(&(q->poster_cv)))
    {
        LOG("ERROR", "In async_q_destroy failed to destroy cv");
        ret = false;
    }

    if(pthread_cond_destroy(&(q->getter_cv)))
    {
        LOG("ERROR", "In async_q_destroy failed to destroy cv");
        ret = false;
    }

    free(q->q);

    return ret;
}

static bool async_post(async_q_t* q, void* elem)
{
    bool ret = true;
    if(pthread_mutex_lock(&(q->lock)))
    {
        LOG("ERROR", "In async_post mutex lock failed");
        return false;
    }

    while(q->n == q->n_max)
    {
        if(pthread_cond_wait(&(q->poster_cv), &(q->lock)))
        {
            LOG("ERROR", "In async_post cond wait failed");
        }
    }

    memcpy(q->q + ((q->n)* (q->element_size)), elem, q->element_size);
    ++(q->n);
    if(pthread_cond_signal(&(q->getter_cv)))
    {
        LOG("ERROR", "In async_post cond signal failed");
        ret = false;
    }


    if(pthread_mutex_unlock(&(q->lock)))
    {
        LOG("ERROR", "In async_post mutex unlock failed");
        ret = false;
    }

    return ret;
    
}

// blocking
static bool async_get(async_q_t* q, void* elem)
{
    bool ret = true;
    if(pthread_mutex_lock(&(q->lock)))
    {
        LOG("ERROR", "In async_get mutex lock failed");
        return false;
    }

    while(q->n == 0)
    {
        if(pthread_cond_wait(&(q->getter_cv), &(q->lock)))
        {
            LOG("ERROR", "In async_get cond wait failed");
        }
    }

    memcpy(elem, q->q + (((q->n)-1)* (q->element_size)), q->element_size);
    --(q->n);

    if(pthread_cond_signal(&(q->poster_cv)))
    {
        LOG("ERROR", "In async_get cond signal failed");
        ret = false;
    }

    if(pthread_mutex_unlock(&(q->lock)))
    {
        LOG("ERROR", "In async_get mutex unlock failed");
        ret = false;
    }

    return ret;

}

static void* _async_thread_func(void* args)
{
    async_runner_t* me = (async_runner_t*) args;
    uint8_t elem[me->q->element_size];

    LOG("INFO", "Async Runner %d starting", me->thread_num);

    while(!(me->killed))
    {
        if(!async_get(me->q, (void*) &elem))
        {
            LOG("ERROR", "In async_runner %d async_get failed");
        }

        if(me->killed)
        {
            break;
        }

        me->consume((void*) &elem);
    }

    LOG("INFO", "Async Runner %d killed", me->thread_num);
    return NULL;
}


static bool async_launch_runners(async_runner_t* runners, 
                          uint16_t n_runners,
                          async_q_t* alloced_queue,
                          async_consumer_t consumer_func)
{

    if(n_runners >= alloced_queue->n_max)
    {
        LOG("ERROR", "In async_launch_runners we do not allow there to be more runners than the size of the queue");
        return false;
    }

    uint16_t i;
    for(i = 0; i < n_runners; ++i)
    {
        async_runner_t* runner = &(runners[i]) ;

        runner->thread_num = i;
        runner->consume = consumer_func;
        runner->q = alloced_queue;
        runner->killed = false;
        runner->n_runners = n_runners;

        if( pthread_create(&(runner->pthread), 
                       NULL, 
                       _async_thread_func,
                       (void*) runner)
        )
        {
            LOG("ERROR", "In async_launch_runners pthread create failed");
            return false;
        }
    }

    return true;

}

static bool async_kill_runners(async_runner_t* runners)
{
    uint16_t i;
    uint16_t n_runners = runners[0].n_runners;
    bool ret = true;

    for(i = 0; i < n_runners; ++i)
    {
        runners[i].killed = true;
    }

    for(i = 0; i < n_runners; ++i)
    {
        uint8_t elem[runners[i].q->element_size];
        if(!async_post((runners[i].q), elem))
        {
            LOG("ERROR", "In async_kill_runners Failed to post");
            ret = false;
        }
    }

    for(i = 0; i < n_runners; ++i)
    {
        if(pthread_join((runners[i].pthread), NULL))
        {
            LOG("ERROR", "In async_kill_runners Failed to join");
            ret = false;
        }
    }

    return ret;

}

#if ASYNC_TEST

static bool check_off[1000] = {0};

static void consumer(void* args)
{
    int arg = *((int*) args);

    if(check_off[arg])
    {
        LOG("ERROR", "TEST FAILED - duplicate elements: %d", arg);
    }
    else
    {
        check_off[arg] = true;
    }
}

static void async_q_test()
{
    uint16_t n_runners = 3;
    uint32_t n_elem = 10;
    async_q_t q;
    async_runner_t runners[n_runners];
    uint32_t i;

    async_q_create(&q, sizeof(int), n_elem);
    async_launch_runners(runners, n_runners, &q, consumer);

    for(i = 0; i < 1000; ++i)
    {
        async_post(&q, (void*) (&i));
    }

    sleep(1);


    for(i = 0; i < 1000; ++i)
    {
        if(check_off == false)
        {
            LOG("ERROR", "TEST FAILED - elem not consumed");
            break;
        }
    }

    LOG("INFO", "TEST PASSED");

    async_kill_runners(runners);
    aysnc_q_destroy(&q);

}

#endif

//*****************************************************************************
// SV Real time data and location prediction and regression. For each GPS sat
// i.e. sv we keep only its most up to date location in WGS84 cordinates. We 
// also keep its SV ID and a time for which that location was recorded in a 
// floating value in UTC time. Finally we keep its calculated instantaneous
// angluar velocity. This how we can predict its path.
//
// We provide two functions. 1) is the SV in some given bounds. This is so we
// can determine which SVs we need to send the client data on. We also provide
// generate path function which given an SV in some given bounds, we predict
// and regress its location to the boundary and fill in a path buffer giving
// its location at uniform time intervals while in the boundary. 
//*****************************************************************************

typedef struct
{   uint16_t sv_id;
    WGS84_t p;    // rad
    WGS84_t v;    // rad / s
    double  t;
} PATH_t;

static PATH_t sv_locs[NUM_SVS];

static bool is_sv_on_map(uint16_t id, 
                         double lat_min, 
                         double lat_max,
                         double lon_min,
                         double lon_max)
{
    if(id >= NUM_SVS)
    {
        LOG("ERROR", "invalid id");
        return false;
    }

    if(lat_min < -M_PI_2 || lat_max > M_PI_2 || lat_max <= lat_min)
    {
        LOG("ERROR", "lat range invalid: [%0.2lf,%0.2lf]", lat_min, lat_max);
        return false;
    }

    if(lon_min < -M_PI || lon_max > M_PI || lon_max <= lon_min)
    {
        LOG("ERROR", "lon range invalid: [%0.2lf,%0.2lf]", lon_min, lon_max);
        return false;
    }

    WGS84_t w = sv_locs[id].p;

    if(w.lat < lat_min || w.lat > lat_max || w.lon < lon_min || w.lon > lon_max)
    {
        return false;
    }

    return true;
}



bool generate_path(uint16_t id, 
                   double lat_min, 
                   double lat_max,
                   double lon_min,
                   double lon_max,
                   PATH_t* path,
                   uint16_t path_size)
{
    if(!is_sv_on_map(id, lat_min,lat_max,lon_min,lon_max) )
    {
        LOG("ERROR", "SV not on map");
        return false;
    }

    if(path == NULL)
    {
        LOG("ERROR", "NULL path passed");
        return false;
    }

    if(path_size < 2)
    {
        LOG("ERROR", "Path size must be greater than 1");
        return false;
    }

    //*************************************************************************
    // We compute the distance and time from the point to each boundary using 
    // the lat and lon rates of change. The two closest boundary points are
    // those whose time to reach them are closest to zero. This allows us to 
    // compute the two points of the path on the boundary
    //*************************************************************************

    WGS84_t p = sv_locs[id].p;
    WGS84_t v = sv_locs[id].v;
    double lat_d_north = lat_max - p.lat;
    double lat_d_south = lat_min - p.lat;
    double lon_d_east = lon_max - p.lon;
    double lon_d_west = lon_min - p.lon;

    double lat_t_north;
    double lat_t_south;
    double lon_t_east;
    double lon_t_west;

    double smallest, next_smallest;

    lat_t_north = ((v.lat == 0.0) ? INFINITY : (lat_d_north / v.lat));
    lat_t_south = ((v.lat == 0.0) ? INFINITY : (lat_d_south / v.lat));
    if(fabs(lat_t_north) < fabs(lat_t_south))
    {
        smallest = lat_t_north;
        next_smallest = lat_d_south;
    }
    else
    {
        smallest = lat_t_south;
        next_smallest = lat_t_north;
    }

    lon_t_east = ((v.lon == 0.0) ? INFINITY : (lon_d_east / v.lon));
    if(fabs(lon_t_east) < fabs(smallest))
    {
        next_smallest = smallest;
        smallest = lon_t_east;
    }
    else if(fabs(lon_t_east) < fabs(next_smallest))
    {
        next_smallest = lon_t_east;
    }

    lon_t_west = ((v.lon == 0.0) ? INFINITY : (lon_d_west / v.lon));
    if(fabs(lon_t_east) < fabs(smallest))
    {
        next_smallest = smallest;
        smallest = lon_t_east;
    }
    else if(fabs(lon_t_east) < fabs(next_smallest))
    {
        next_smallest = lon_t_east;
    }

    WGS84_t path_bd_first;
    WGS84_t path_bd_second;

    if(smallest == lat_t_north)
    {
        path_bd_first.lat = lat_max;
        path_bd_first.lon = v.lon*lat_t_north;
    }
    else if(smallest == lat_t_south)
    {
        path_bd_first.lat = lat_min;
        path_bd_first.lon = v.lon*lat_t_south;
    }
    else if(smallest == lon_t_east)
    {
        path_bd_first.lon = lon_max;
        path_bd_first.lat = v.lat*lon_t_east;
    }
    else if(smallest == lon_t_west)
    {
        path_bd_first.lon = lon_min;
        path_bd_first.lat = v.lat*lon_t_west;
    }

    if(next_smallest == lat_t_north)
    {
        path_bd_second.lat = lat_max;
        path_bd_second.lon = v.lon*lat_t_north;
    }
    else if(next_smallest == lat_t_south)
    {
        path_bd_second.lat = lat_min;
        path_bd_second.lon = v.lon*lat_t_south;
    }
    else if(next_smallest == lon_t_east)
    {
        path_bd_second.lon = lon_max;
        path_bd_second.lat = v.lat*lon_t_east;
    }
    else if(next_smallest == lon_t_west)
    {
        path_bd_second.lon = lon_min;
        path_bd_second.lat = v.lat*lon_t_west;
    }

    #if PATH_DB_TEST
        LOG("DEBUG", "");
        LOG("DEBUG", "Path Boundaries for id=%d", id);
        LOG("DEBUG", "lat_bd = [%0.2lf x %0.2lf]   lon_bd = [%0.2lf x %0.2lf]", deg(lat_min), deg(lat_max), deg(lon_min), deg(lon_max));
        LOG("DEBUG", "P = [%0.2lf , %0.2lf]   V = [%0.2lf , %0.2lf]", deg(p.lat), deg(p.lon), deg(v.lat), deg(v.lon));
        LOG("DEBUG", "BD_1 = [%0.2lf , %0.2lf]   BD_2 = [%0.2lf , %0.2lf]", deg(path_bd_first.lat), deg(path_bd_first.lon), deg(path_bd_second.lat), deg(path_bd_second.lon));
        LOG("DEBUG", "lat_t_north = %0.2lf", lat_t_north);
        LOG("DEBUG", "lat_t_south = %0.2lf", lat_t_south);
        LOG("DEBUG", "lon_t_east = %0.2lf", lon_t_east);
        LOG("DEBUG", "lon_t_west = %0.2lf", lon_t_west);
        LOG("DEBUG", "");
    #endif

    //*************************************************************************
    // We now have the boundary points where the projected path of the SV 
    // intersects the boundary. We need to compute the distance between these
    // points in each dir and divide by the number points in the passed path 
    // buffer. This gives us the step size in each dir of the path points.
    // Using the step size populate the buffer.
    //*************************************************************************

    double lat_step_size = fabs(path_bd_first.lat - path_bd_second.lat);
    lat_step_size = lat_step_size / (((double)path_size) -1.0);

    double lon_step_size = fabs(path_bd_first.lon - path_bd_second.lon);
    lon_step_size = lon_step_size / (((double)path_size) -1.0);

    double t_step_size = fabs(smallest - next_smallest);
    t_step_size = t_step_size / (((double)path_size) -1.0);

    if(smallest < 0.0)
    {
        path[0].p.lat = path_bd_first.lat;
        path[0].p.lon = path_bd_first.lon;
        path[0].p.alt = sv_locs[id].p.alt;
        path[0].t = smallest + sv_locs[id].t;
    }
    else
    {
        path[0].p.lat = path_bd_second.lat;
        path[0].p.lon = path_bd_second.lon;
        path[0].p.alt = sv_locs[id].p.alt;
        path[0].t = next_smallest +  + sv_locs[id].t;
    }

    path[0].sv_id = id;

    if(sv_locs[id].v.lat < 0.0)
    {
        lat_step_size *= -1.0;
    }

    if(sv_locs[id].v.lon < 0.0)
    {
        lon_step_size *= -1.0;
    }

    uint16_t i;
    for(i = 1; i < path_size; ++i)
    {
        path[i].p.lat = path[0].p.lat + (i*lat_step_size);
        path[i].p.lon = path[0].p.lon + (i*lon_step_size);
        path[i].t = path[0].t + (i*t_step_size);
        path[i].sv_id = id;
    }

    return true;
}

// for now just put some random bs in there
void pull_init_sv_locs()
{
    uint16_t i;
    for(i = 0; i < NUM_SVS; ++i)
    {
        sv_locs[i].p.lat = M_PI_2 * (((double) rand() / (double) RAND_MAX) - 1.0 );
        sv_locs[i].p.lon = M_PI * (((double) rand() / (double) RAND_MAX) - 1.0);
        sv_locs[i].p.alt = 0.0;

        sv_locs[i].v.lat = (((double) rand() / (double) 0x7fffffff) - 1.0) / 10.0;
        sv_locs[i].v.lon = (((double) rand() / (double) 0x7fffffff) - 1.0) / 10.0;

        LOG("DEBUG", "Sat %d  p=%0.2lf, %0.2lf  v=%0.2lf,%0.2lf", 
            deg(sv_locs[i].p.lat),
            deg(sv_locs[i].p.lon),
            deg(sv_locs[i].v.lat),
            deg(sv_locs[i].v.lon));

        sv_locs[i].t = 0.0;

        sv_locs[i].sv_id = i;
    }
}

#if PATH_DB_TEST
void path_gen_test()
{
    uint16_t i;
    PATH_t path[50];
    for(i = 0; i < NUM_SVS; ++i)
    {
        generate_path(i, 
                      rad(-50.0), 
                      rad(50.0),
                      rad(-120.0),
                      rad(120.0),
                      &path,
                      50);

    }
}
#endif

//*****************************************************************************
// Web Socket interface. Our front end UI is running leaflet js. It has an
// open websocket to this backend. The backend has 3 threads:
//
// -------------
// | Listening |
// -------------
//      |        
//      | Conn recv, Spawn     ---------------
//      |--------------------->| Recv Thread |
//      |                      ---------------
//      | (blocks incoming)         |
//      |                           | Spawn         ----------------
//      |                           |-------------->| Async Runner |
//      |                           |               ----------------
//      |                           |                       |
//      |                           | Process req, Q-post   |
//      |                           |---------------------->|
//      |                           |                       | Q-get -> Send
//      |                           |                       |----------------->
//      |                           |    conn closed, kill  |
//      |   die, client count -> 0  |---------------------->|
//      |<--------------------------|
//      |
//      V        
//
// The listening thread replaces the main thread. It opens a TCP listening
// socket. When an incoming connection is accepted it launches a receive thread
// to handle input from the client. In ws.h we set MAX CLIENTS to 1 so when 
// a client is connected, the listening thread will block all other incoming
// connections till that one client is closed. The listneing thread logic is
// wholey contained in the wsServer.
//
// The receive immediatly calls the registered call back "onopen" once
// spawned. It creates an async_q and a single async_runner. These async 
// objects are created to handle output from the backend to the front end 
// client. The async queue takes a tuple (sv ID, WGS84 point, utc time in
// seconds) as its q elements ,packages them, and sends them to the front end 
// client. The receivers main purpose is to handle requests from the client. 
// On input, the onmessage call back is called which processes the request and 
// posts to the async Q ID, point, UTC time tuples to be sent to the client. 
// When a client closes a session, the onclose call back is invoked which 
// cleans the async resource. The recv thread then gracefully exits. Most the 
// logic besides the onopen, onclose and onmessage are contained in wsServer.
//*****************************************************************************

static uint32_t path_q_size = 100;
static ws_cli_conn_t *ui_client = NULL;
static async_runner_t path_runner;
static async_q_t path_q;


void consome_path(void* arg)
{
    PATH_t path = *((PATH_t*) arg);
    WGS84_t w = path.p;
    char send_str[WS_RESPONSE_MSG_SIZE];
    snprintf(send_str, WS_RESPONSE_MSG_SIZE, "%d,%0.10lf,%0.10lf,%0.10lf", 
        path.sv_id, 
        deg(w.lat), 
        deg(w.lon),
        path.t
    );
    LOG("INFO", "Sending: %s", send_str);
    ws_sendframe_txt(ui_client, send_str);
}

void onopen(ws_cli_conn_t *client)
{
	char *cli, *port;
	cli  = ws_getaddress(client);
	port = ws_getport(client);
	LOG("INFO", "WS Connection opened, addr: %s, port: %s", cli, port);

    if(ui_client != NULL)
    {
        LOG("WARN", "Multiple clients connected .. closing duplicate client socket");
        ws_close_client(client);
    }
    ui_client = client;

    async_q_create(&path_q, sizeof(PATH_t), path_q_size);
    async_launch_runners(&path_runner, 1, &path_q, consome_path);
}

void onclose(ws_cli_conn_t *client)
{
	char *cli, *port;
	cli = ws_getaddress(client);
    port = ws_getport(client);
	LOG("INFO", "WS Connection closed, addr: %s, port: %s", cli, port);

    async_kill_runners(&path_runner);
    aysnc_q_destroy(&path_q);

    ui_client = NULL;
}

void onmessage(ws_cli_conn_t *client,
	const unsigned char *msg, uint64_t size, int type)
{
	char *cli, *port;
	cli = ws_getaddress(client);
    port = ws_getport(client);
	LOG("INFO", "I receive a message: %s (size: %" PRId64 ", type: %d), from: %s:%s",
		msg, size, type, cli, port);

    double lat_min, lat_max, lon_min, lon_max;
    int parsed = sscanf(msg, "%lf,%lf,%lf,%lf", 
        &lat_min,
        &lat_max,
        &lon_min,
        &lon_max
    );

    lat_min = rad(lat_min);
    lat_max = rad(lat_max);
    lon_min = rad(lon_min);
    lon_max = rad(lon_max);

    if(parsed != 4)
    {
        LOG("ERROR", "Incoming message parse error");
        return;
    }

    uint16_t i;
    PATH_t path[2];
    for(i = 0; i < NUM_SVS; ++i)
    {
        if(is_sv_on_map(i, lat_min, lat_max, lon_min, lon_max))
        {
            generate_path(i, 
                lat_min, 
                lat_max, 
                lon_min, 
                lon_max, 
                path, 
                2
            );

            async_post(&path_q, &path[0]);
            async_post(&path_q, &sv_locs[i]);
            async_post(&path_q, &path[1]);
        }
    }

}
 
//*****************************************************************************
// Entry
//*****************************************************************************

int main(void)
{

    pull_init_sv_locs();

    ws_socket(&(struct ws_server){
		.host = WS_SERVER_IP,
		.port = WS_SERVER_PORT,
		.thread_loop   = 0,
		.timeout_ms    = 1000,
		.evs.onopen    = &onopen,
		.evs.onclose   = &onclose,
		.evs.onmessage = &onmessage
	});

    // Does not return gets replace with the listening thread
    
    return 0;
}