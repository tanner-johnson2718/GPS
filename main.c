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

#define MAX_CLIENTS 1
#include "ws.h"

#define ASYNC_TEST 0
#define GEOCORD_TEST 0


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
// Path Generator. All points in ECEF
//*****************************************************************************

typedef struct
{
    ECEF_t p;      // Pos (m)
    ECEF_t v;      // Vel (m/s)
    ECEF_t a;      // Acc (m/s2)
    double  t;     // Time in simulation (sec).(ns)
} PATH_t;

// We assume path[0] and path[n-1] has valid data giving ini and fini cond.
// returns success or not.
bool path_gen(PATH_t* path, int n)
{
    if(path[0].t >= path[n-1].t )
    {
        LOG("ERROR", "In path_gen end point has timestamp before start point");
        return false;
    }

    if(n < 2)
    {
        LOG("ERROR", "In path_gen must have a start and end data point in passed points");
        return false;
    }

    PATH_t delta;
    delta.t = (path[n-1].t - path[0].t) / (double) n;
    delta.p = (ECEF_t) {
        (path[n-1].p.x - path[0].p.x) / (double) n,
        (path[n-1].p.y - path[0].p.y) / (double) n,
        (path[n-1].p.z - path[0].p.z) / (double) n
    };

    delta.v = (ECEF_t) {
        delta.p.x / delta.t,
        delta.p.y / delta.t,
        delta.p.z / delta.t
    };

    delta.a = (ECEF_t) {0.0,0.0,0.0};

    // Overwrite init cond v and a
    path[0].v = delta.v;
    path[0].a = delta.a;
    path[n-1].v = delta.v;
    path[n-1].a = delta.a;

    int i;
    for(i = 1; i < n-1; ++i)
    {
        path[i].v = delta.v;
        path[i].a = delta.a;

        path[i].p.x = path[0].p.x + ((delta.p.x) * i);
        path[i].p.y = path[0].p.y + ((delta.p.y) * i);
        path[i].p.z = path[0].p.z + ((delta.p.z) * i);
    }

    return 0;

}

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
// Web Socket interface. 
//*****************************************************************************

static uint32_t path_q_size = 100;
static ws_cli_conn_t *ui_client = NULL;
static async_runner_t path_runner;
static async_q_t path_q;

void consome_path(void* arg)
{
    WGS84_t point = *((WGS84_t*) arg);
    char send_str[48];
    snprintf(send_str, 100, "%15lf,%15lf", deg(point.lat), deg(point.lon));
    LOG("INFO", "Sending: %s", send_str);
    ws_sendframe_txt(ui_client, send_str);
}

void post_path()
{
    PATH_t path[10];
    path[0].p = ECEF((WGS84_t) {rad(38.0), rad(-105.0), 2000});
    path[0].v = (ECEF_t) {0,0,0};
    path[0].a = (ECEF_t) {0,0,0};
    path[0].t = 0.0;
    path[9].p = ECEF((WGS84_t) {rad(38.0), rad(-115.0), 2000});
    path[9].v = (ECEF_t) {0,0,0};
    path[9].a = (ECEF_t) {0,0,0};
    path[9].t = 10.0;
    path_gen(path, 10);

    int i;
    for(i = 0; i < 10; i++)
    {
        WGS84_t w = WGS84(path[i].p);
        async_post(&path_q, &w);
        sleep(1);
    }
    
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


    async_q_create(&path_q, sizeof(WGS84_t), path_q_size);
    async_launch_runners(&path_runner, 1, &path_q, consome_path);

    post_path();

    async_kill_runners(&path_runner);
    aysnc_q_destroy(&path_q);

    ws_close_client(client);
    ui_client = NULL;

}

void onclose(ws_cli_conn_t *client)
{
	char *cli, *port;
	cli = ws_getaddress(client);
    port = ws_getport(client);
	LOG("INFO", "WS Connection closed, addr: %s, port: %s", cli, port);
}

void onmessage(ws_cli_conn_t *client,
	const unsigned char *msg, uint64_t size, int type)
{
	char *cli, *port;
	cli = ws_getaddress(client);
    port = ws_getport(client);
	LOG("INFO", "I receive a message: %s (size: %" PRId64 ", type: %d), from: %s:%s",
		msg, size, type, cli, port);
}


//*****************************************************************************
// Entry
//*****************************************************************************

int main(void)
{

    ws_socket(&(struct ws_server){
		.host = "127.0.0.1",
		.port = 6969,
		.thread_loop   = 0,
		.timeout_ms    = 1000,
		.evs.onopen    = &onopen,
		.evs.onclose   = &onclose,
		.evs.onmessage = &onmessage
	});

    // Does not return gets replace with the listening thread
    
    return 0;
}