#pragma once

// GTSAM related includes.
#include <gtsam/gnssNavigation/GnssTools.h>
#include <gtsam/gnssNavigation/PhaseFactor.h>
#include <gtsam/gnssNavigation/PseudorangeFactor.h>
#include <gtsam/gnssNavigation/nonBiasStates.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/ISAM2.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/BatchFixedLagSmoother.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/slam/PriorFactor.h>

// STD
#include <unistd.h>
#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include "rtklib.h"

using namespace std;
using namespace gtsam;
namespace NM = gtsam::noiseModel;
typedef noiseModel::Diagonal DiagNoise;

// Intel Threading Building Block
#ifdef GTSAM_USE_TBB
#include <tbb/tbb.h>
#undef max  // TBB seems to include windows.h and we don't want these macros
#undef min
#endif

/************************************************************************/
#define MAXFILE 16 /* max number of input files */

// using symbol_shorthand::G;  // bias states ( Phase Biases )
// using symbol_shorthand::X;  // nonBiasStates ( dx, dy, dz, trop, cb )
#define C(x) (300000 + x)   // C: 3
#define N(x) (1400000 + x)  // N: 14
#define T(x) (2000000 + x)  // T: 20
#define X(x) (2500000 + x)  // X: 25
#define THRES_MW_JUMP 10.0
#define SQR(x) ((x) * (x))
#define ROUND(x) (int)floor((x) + 0.5)

#define VAR_POS SQR(60.0)  /* init variance receiver position (m^2) */
#define VAR_CLK SQR(60.0)  /* init variance receiver clock (m^2) */
#define VAR_ISB SQR(60.0)  /* init variance inter-system bias (m^2) */
#define VAR_ZTD SQR(0.6)   /* init variance ztd (m^2) */
#define VAR_GRA SQR(0.01)  /* init variance gradient (m^2) */
#define VAR_DCB SQR(30.0)  /* init variance dcb (m^2) */
#define VAR_BIAS SQR(60.0) /* init variance phase-bias (m^2) */
#define VAR_IONO SQR(60.0) /* init variance iono-delay */

#define ERR_SAAS 0.3      /* saastamoinen model error std (m) */
#define ERR_BRDCI 0.5     /* broadcast iono model error factor */
#define ERR_CBIAS 0.3     /* code bias error std (m) */
#define REL_HUMI 0.7      /* relative humidity for saastamoinen model */
#define GAP_RESION 120    /* default gap to reset ionos parameters (ep) */
#define EFACT_GPS_L5 10.0 /* error factor of GPS/QZS L5 */

/* number and index of states */
#define NF(opt) ((opt)->ionoopt == IONOOPT_IFLC ? 1 : (opt)->nf)
#define NP(opt) ((opt)->dynamics ? 9 : 3)
#define NC(opt) (NSYS)
#define NT(opt) \
  ((opt)->tropopt < TROPOPT_EST ? 0 : ((opt)->tropopt == TROPOPT_EST ? 1 : 3))
#define NI(opt) ((opt)->ionoopt == IONOOPT_IFLC ? 0 : MAXSAT)
#define ND(opt) ((opt)->nf >= 3 ? 1 : 0)
#define NR(opt) (NP(opt) + NC(opt) + NT(opt) + NI(opt) + ND(opt))
#define NB(opt) (NF(opt) * MAXSAT)
#define NX(opt) (NR(opt) + NB(opt))
#define IC(s, opt) (NP(opt) + (s))
#define IT(opt) (NP(opt) + NC(opt))
#define II(s, opt) (NP(opt) + NC(opt) + NT(opt) + (s)-1)
#define ID(opt) (NP(opt) + NC(opt) + NT(opt) + NI(opt))
#define IB(s, f, opt) (NR(opt) + MAXSAT * (f) + (s)-1)

static int nepoch = 0; /* number of observation epochs */
static int iobsu = 0;  /* current rover observation data arIndex_InFiles */
static int iobsr = 0;  /* current reference observation data arIndex_InFiles */
static char proc_rov[64] = "";  /* rover for current processing */
static char proc_base[64] = ""; /* base station for current processing */

/* CP data structure */
typedef int SatPRN;  // Satellite PRN Number
typedef int SatSYS;  // Satellite System

typedef struct { /* preprocessd GNSS data */
  double GNSS_time;
  gtime_t obsTime, solTime;
  int total_sv, prn, sys;    /* satellite total/cur prn/system */
  double SNR[NFREQ];         /* signal strength (1 dBHz) */
  unsigned char LLI[NFREQ];  /* loss of lock indicator */
  double P[NFREQ], L[NFREQ]; /* raw pseudorange(m)/carrier-phase(cycle) */
  double lamda[NFREQ];       /* carrier wave lengths (m) */
  double code[NFREQ];        /* code indicator (CODE_???) */
  int slip;                  // cycle-slip
  int phasebreak;            // cycle-slip or new
  double newBias;            // new phase bias
  double elevation, azimuth; /* elevation/azimuth(rad) */
  double err_tropo, err_iono, SatClk;
  Vector3 SatXYZ, SatVel; /* sat position (ecef, m) */
  double SatVar;
  double P_crr; /* pseudorange corrected measurement (m) */
  // double GNSSCovariance;
  int visable;  // 0-Not sure 1-visable 2-invisable;
  // string sat_system;
  int visable3DMA;  // ground truth visibility from 3DMA;
  double prE3dMA;   // ground truth pseudorange error from 3DMA;
} PreGNSS;

typedef struct {
  rtk_t rtk;
  list<obsd_t> lCurRawObs;
  map<SatPRN, PreGNSS> mCurPntData;
  void ClearCurdata() {
    lCurRawObs.clear();
    mCurPntData.clear();
  }
} CurData;

// // /************************* functions *************************/
int Showmsg(const char *format, ...);
int Checkbrk(const char *format, ...);
void Readpreceph(const char *Sp3File, const char *ClkFile, char *FcbFile,
                 const prcopt_t *prcopt, nav_t *nav);
int Readerp(const char *file, erp_t *erp);
void Readotl(prcopt_t *sProcOpt, const char *file, const sta_t *sta);
int Readobsnav(gtime_t sGTimeStart, gtime_t sGTimeEnd, double PrcsInterval,
               const char *ObsFile, const char *NavFile, const prcopt_t *prcopt,
               obs_t *obs, nav_t *nav, sta_t *sta);

int Nextobsf(const obs_t *obs, int *i, int rcv);
int Nextobsb(const obs_t *obs, int *i, int rcv);

void Corr_phase_bias_fcb(obsd_t *obs, int n, const nav_t *nav);
void Corr_phase_bias_ssr(obsd_t *obs, int n, const nav_t *nav);

double Gfmeas(const PreGNSS &CurPreGNSS, const nav_t *nav);
double Mwmeas(const PreGNSS &CurPreGNSS, const nav_t *nav);
void Detslp_ll(rtk_t *rtk, const map<SatPRN, PreGNSS> &mCurPntData);
void Detslp_gf(rtk_t *rtk, const map<SatPRN, PreGNSS> &mCurPntData,
               const nav_t *nav);
void Detslp_mw(rtk_t *rtk, const map<SatPRN, PreGNSS> &mCurPntData,
               const nav_t *nav);

void Testeclipse(const gtime_t &CurTime, map<SatPRN, PreGNSS> &mCurPntData,
                 const nav_t *nav);
void Corr_meas(const PreGNSS &CurPreGNSS, const nav_t *nav, const double *azel,
               const prcopt_t *opt, const double *dantr, const double *dants,
               double phw, double *L, double *P, double *Lc, double *Pc);
               
// double Sdobs(const PreGNSS &uCurSatPre, const PreGNSS &rCurSatPre, int f);
// void Detslp_ll(rtk_t *rtk, const PreGNSS &CurSatPre, int rcv);
// void Detslp_gf_L1L2(rtk_t *rtk, const PreGNSS &uCurSatPre,
//                     const PreGNSS &rCurSatPre, const nav_t *nav);
// void Detslp_gf_L1L5(rtk_t *rtk, const PreGNSS &uCurSatPre,
//                     const PreGNSS &rCurSatPre, const nav_t *nav);
// double Varerr(int sys, double el, double bl, double dt, int f,
//               const prcopt_t *opt);
// double SumNoise(SatSYS sys, double ele_i, double ele_j, double bl, double dt,
//                 int f, const prcopt_t *opt);