/*
 * @file pppBayesTree.cpp
 * @brief Iterative GPS Range/Phase Estimator with collected data
 * @author Ryan Watson & Jason Gross
 */

#include "PPPSamlib.h"

static obs_t gObss = {0};
static nav_t gNavs = {0};
static sta_t gStas[64];           /* station infomation */
static pcvs_t gRecAntParas = {0}; /* receiver antenna parameters */
static pcvs_t gSatAntParas = {0}; /* satellite antenna parameters */

vector<int> gAllGNKey_(NSATGPS + 2, 0);  // index-2-prn
vector<double> gvNavTrue, gvInitPosNoise;
double gInitClkNoise, gInitAmNoise, gBetwTrpNoise, gPhNoiseFactor;

bool gbFirstEpoch(true);
int gbCurBreak;
Values gNewInitValues;
Values gResultValues;
NonlinearFactorGraph *gNewGraph;

/* open procssing session ----------------------------------------------------*/
int Openses(const filopt_t &sFileOpt) {
  int i;

  trace(3, "Openses :\n");

  /* read satellite antenna parameters */
  if (sFileOpt.satantp && !(readpcv(sFileOpt.satantp, &gSatAntParas))) {
    Showmsg("error : no sat ant pcv in %s", sFileOpt.satantp);
    trace(1, "sat antenna pcv read error: %s\n", sFileOpt.satantp);
    return 0;
  }
  /* read receiver antenna parameters */
  if (sFileOpt.rcvantp && !(readpcv(sFileOpt.rcvantp, &gRecAntParas))) {
    Showmsg("error : no rec ant pcv in %s", sFileOpt.rcvantp);
    trace(1, "rec antenna pcv read error: %s\n", sFileOpt.rcvantp);
    return 0;
  }

  /* use satellite L2 offset if L5 offset does not exists */
  for (i = 0; i < gRecAntParas.n; i++) {
    if (norm(gRecAntParas.pcv[i].off[2], 3) > 0.0) continue;
    matcpy(gRecAntParas.pcv[i].off[2], gRecAntParas.pcv[i].off[1], 3, 1);
    matcpy(gRecAntParas.pcv[i].var[2], gRecAntParas.pcv[i].var[1], 19, 1);
  }
  for (i = 0; i < gSatAntParas.n; i++) {
    if (norm(gSatAntParas.pcv[i].off[2], 3) > 0.0) continue;
    matcpy(gSatAntParas.pcv[i].off[2], gSatAntParas.pcv[i].off[1], 3, 1);
    matcpy(gSatAntParas.pcv[i].var[2], gSatAntParas.pcv[i].var[1], 19, 1);
  }
  return 1;
}
/* set antenna parameters ----------------------------------------------------*/
void Setpcv(gtime_t time, prcopt_t *sProcOpt, nav_t *nav, const pcvs_t *pcvs,
            const pcvs_t *pcvr, const sta_t *sta) {
  pcv_t *pcv;
  double pos[3], del[3];
  int i, j,
      mode = PMODE_DGPS <= sProcOpt->mode && sProcOpt->mode <= PMODE_FIXED;
  char id[64];

  /* set satellite antenna parameters */
  for (i = 0; i < MAXSAT; i++) {
    if (!(satsys(i + 1, NULL) & sProcOpt->navsys)) continue;
    if (!(pcv = searchpcv(i + 1, "", time, pcvs))) {
      satno2id(i + 1, id);
      trace(3, "no satellite antenna pcv: %s\n", id);
      continue;
    }
    nav->pcvs[i] = *pcv;
  }
  for (i = 0; i < (mode ? 2 : 1); i++) {
    if (!strcmp(sProcOpt->anttype[i], "*")) { /* set by station parameters */
      strcpy(sProcOpt->anttype[i], sta[i].antdes);
      if (sta[i].deltype == 1) { /* xyz */
        if (norm(sta[i].pos, 3) > 0.0) {
          ecef2pos(sta[i].pos, pos);
          ecef2enu(pos, sta[i].del, del);
          for (j = 0; j < 3; j++) sProcOpt->antdel[i][j] = del[j];
        }
      } else { /* enu */
        for (j = 0; j < 3; j++) sProcOpt->antdel[i][j] = gStas[i].del[j];
      }
    }
    if (!(pcv = searchpcv(0, sProcOpt->anttype[i], time, pcvr))) {
      trace(2, "no receiver antenna pcv: %s\n", sProcOpt->anttype[i]);
      *sProcOpt->anttype[i] = '\0';
      continue;
    }
    strcpy(sProcOpt->anttype[i], pcv->type);
    sProcOpt->pcvr[i] = *pcv;
  }
}
/* input obs data, navigation messages and sbas correction -------------------*/
int Inputobs(obsd_t *obs, int solq, const prcopt_t *sProcOpt) {
  gtime_t time = {0};
  int i, nu, nr, n = 0;

  if (0 <= iobsu && iobsu < gObss.n) {
    time = gObss.data[iobsu].time;
    if (Checkbrk("processing : %s Q=%d", time_str(time, 0), solq)) {
      Showmsg("aborted");
      return -1;
    }
  }

  /* input forward data */
  if ((nu = Nextobsf(&gObss, &iobsu, 1)) <= 0) return -1;
  if (sProcOpt->intpref) {
    for (; (nr = Nextobsf(&gObss, &iobsr, 2)) > 0; iobsr += nr)
      if (timediff(gObss.data[iobsr].time, gObss.data[iobsu].time) > -DTTOL)
        break;
  } else {
    for (i = iobsr; (nr = Nextobsf(&gObss, &i, 2)) > 0; iobsr = i, i += nr)
      if (timediff(gObss.data[i].time, gObss.data[iobsu].time) > DTTOL) break;
  }
  nr = Nextobsf(&gObss, &iobsr, 2);
  for (i = 0; i < nu && n < MAXOBS * 2; i++) obs[n++] = gObss.data[iobsu + i];
  for (i = 0; i < nr && n < MAXOBS * 2; i++) obs[n++] = gObss.data[iobsr + i];
  iobsu += nu;

  return n;
}
/* temporal update of phase biases -------------------------------------------*/
void Udbias_ppp(gtime_t CurTime, rtk_t *rtk, const nav_t *nav,
                Values &ResultValues, map<SatPRN, PreGNSS> &mCurPntData) {
  const double *lam;
  double L[NFREQ], P[NFREQ], Lc, Pc, offset = 0.0, pos[3] = {0};
  double ion, vion, dantr[NFREQ] = {0}, dants[NFREQ] = {0};
  int i, j, k, l, f, sat, clk_jump = 0;

  /* handle day-boundary clock jump */
  if (rtk->opt.posopt[5]) {
    clk_jump = ROUND(time2gpst(CurTime, NULL) * 10) % 864000 == 0;
  }
  for (i = 0; i < MAXSAT; i++)
    for (j = 0; j < rtk->opt.nf; j++) {
      rtk->ssat[i].slip[j] = 0;
    }
  /* detect cycle slip by LLI */
  Detslp_ll(rtk, mCurPntData);

  /* detect cycle slip by geometry-free phase jump */
  Detslp_gf(rtk, mCurPntData, nav);

  /* detect slip by Melbourne-Wubbena linear combination jump */
  Detslp_mw(rtk, mCurPntData, nav);

  ecef2pos(rtk->sol.rr, pos);

  for (f = 0; f < NF(&rtk->opt); f++) {
    k = 0;
    for (auto &mpit : mCurPntData) {
      sat = mpit.second.prn;

      /* reset phase-bias if expire obs outage counter */
      if (++rtk->ssat[sat - 1].outc[f] > (unsigned int)rtk->opt.maxout ||
          rtk->opt.modear == ARMODE_INST || clk_jump) {
        mpit.second.phasebreak = 1;
      }

      Corr_meas(mpit.second, nav, rtk->ssat[sat - 1].azel, &rtk->opt, dantr,
                dants, 0.0, L, P, &Lc, &Pc);

      mpit.second.newBias = 0.0;

      mpit.second.slip = 0;
      if (rtk->opt.ionoopt == IONOOPT_IFLC) {
        mpit.second.newBias = Lc - Pc;
        mpit.second.slip =
            rtk->ssat[sat - 1].slip[0] || rtk->ssat[sat - 1].slip[1];
      } else if (L[f] != 0.0 && P[f] != 0.0) {
        mpit.second.slip = rtk->ssat[sat - 1].slip[f];
        l = satsys(sat, NULL) == SYS_GAL ? 2 : 1;
        lam = nav->lam[sat - 1];
        /*if (obs[i].P[0] == 0.0 || obs[i].P[l] == 0.0 ||
                lam[0] == 0.0 || lam[l] == 0.0 || lam[f] == 0.0) continue;
        ion = (obs[i].P[0] - obs[i].P[l]) / (1.0 - SQR(lam[l] / lam[0]));*/
        if (!ionocorr(CurTime, nav, sat, pos, rtk->ssat[sat - 1].azel,
                      rtk->opt.ionoopt, &ion, &vion))
          continue;
        mpit.second.newBias = L[f] - P[f] + 2.0 * ion * SQR(lam[f] / lam[0]);
      }
      if (!ResultValues.exists(N(gAllGNKey_[sat])) || mpit.second.slip ||
          mpit.second.newBias == 0.0)
        continue;

      offset +=
          mpit.second.newBias - ResultValues.at<double>(N(gAllGNKey_[sat]));
      k++;
    }
    /* correct phase-code jump to ensure phase-code coherency */
    if (k >= 2 && fabs(offset / k) > 0.0005 * CLIGHT) {
      for (auto &mpit : mCurPntData) {
        sat = mpit.second.prn;

        if (ResultValues.exists(N(gAllGNKey_[sat])))
          mpit.second.newBias += offset / k;
      }
    }
    /* calculate phase bias */
    for (auto &mpit : mCurPntData) {
      sat = mpit.second.prn;

      mpit.second.phasebreak = 0;
      if (mpit.second.newBias == 0.0 ||
          (ResultValues.exists(N(gAllGNKey_[sat])) && !mpit.second.slip))
        continue;

      /* reinitialize phase-bias if detecting cycle slip */
      mpit.second.phasebreak = 1;

      /* reset fix flags */
      for (k = 0; k < MAXSAT; k++) rtk->ambc[sat - 1].flags[k] = 0;
    }
  }
}
/* measurement error variance ------------------------------------------------*/
double Varerr(int sat, int sys, double el, int freq, int type,
              const prcopt_t *opt) {
  double fact = 1.0, sinel = sin(el);

  if (type == 1) fact *= opt->eratio[freq == 0 ? 0 : 1];
  // fact *= sys == SYS_GLO ? EFACT_GLO : (sys == SYS_SBS ? EFACT_SBS :
  // EFACT_GPS);
  switch (sys) {
    case SYS_GPS:
      fact *= EFACT_GPS;
      break;
    case SYS_GLO:
      fact *= EFACT_GLO;
      break;
    case SYS_GAL:
      fact *= EFACT_GAL;
      break;
    case SYS_QZS:
      fact *= EFACT_QZS;
      break;
    case SYS_SBS:
      fact *= EFACT_SBS;
      break;
    case SYS_CMP:
      fact *= EFACT_CMP;
      break;
    case SYS_IRN:
      fact *= EFACT_IRN;
      break;
  }

  if (sys == SYS_GPS || sys == SYS_QZS) {
    if (freq == 2) fact *= EFACT_GPS_L5; /* GPS/QZS L5 error factor */
  }
  if (opt->ionoopt == IONOOPT_IFLC) fact *= 3.0;
  return SQR(fact * opt->err[1]) + SQR(fact * opt->err[2] / sinel);
}
/* single-point positioning ----------------------------------------------------
 * compute receiver position, velocity, clock bias by single-point positioning
 * with pseudorange and doppler observables
 * args   : obsd_t *obs      I   observation data
 *          int    n         I   number of observation data
 *          nav_t  *nav      I   navigation data
 *          prcopt_t *opt    I   processing options
 *          sol_t  *sol      IO  solution
 *          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
 *          ssat_t *ssat     IO  satellite status              (NULL: no output)
 *          char   *msg      O   error message for error exit
 * return : status(1:ok,0:error)
 * notes  : assuming sbas-gps, galileo-gps, qzss-gps, compass-gps time offset
 *and receiver bias are negligible (only involving glonass-gps time offset and
 *receiver bias)
 *-----------------------------------------------------------------------------*/
int Pntpos(obsd_t *obs, int *ntemp, const nav_t *nav, const prcopt_t *opt,
           sol_t *sol, double *azel, ssat_t *ssat, char *msg,
           map<int, PreGNSS> &mCurPntData) {
  int n = *ntemp;
  prcopt_t opt_ = *opt;
  double *rs, *dts, *var, *azel_, *resp;
  int i, stat, vsat[MAXOBS] = {0}, svh[MAXOBS];

  sol->stat = SOLQ_NONE;

  if (n <= 0) {
    strcpy(msg, "no observation data");
    return 0;
  }

  trace(3, "Pntpos  : tobs=%s n=%d\n", time_str(obs[0].time, 3), n);

  sol->time = obs[0].time;
  msg[0] = '\0';

  rs = mat(6, n);
  dts = mat(2, n);
  var = mat(1, n);
  azel_ = zeros(2, n);
  resp = mat(1, n);

  if (opt_.mode != PMODE_SINGLE) { /* for precise positioning */
#if 0
		opt_.sateph =EPHOPT_BRDC;
#endif
    // opt_.ionoopt = IONOOPT_BRDC;
    if (opt_.nf > 1 && obs[0].L[0] && obs[0].L[1]) {
      opt_.ionoopt = IONOOPT_IFLC;
    } else {
      opt_.ionoopt = IONOOPT_BRDC;
    }

    opt_.tropopt = TROPOPT_SAAS;
  }
  /* satellite positons, velocities and clocks */
  satposs(sol->time, obs, n, nav, opt_.sateph, rs, dts, var, svh);

  /* estimate receiver position with pseudorange */
  stat = estpos(obs, n, rs, dts, var, svh, nav, &opt_, sol, azel_, vsat, resp,
                msg);

  // /* Jin: restore the preprocessed information via converged SPP */
  double CurTow = time2gpst(obs[0].time, 0), elemask = opt->elmin;
  for (int s_i = 0; s_i < n; s_i++) {
    PreGNSS GNSSPre{0};
    GNSSPre.GNSS_time = CurTow;
    GNSSPre.obsTime = obs[s_i].time;
    GNSSPre.solTime = sol->time;
    GNSSPre.total_sv = n;
    GNSSPre.prn = obs[s_i].sat;
    GNSSPre.sys = satsys(obs[s_i].sat, NULL);

    for (int k = 0; k < NFREQ; k++) {
      /* validate the obs */
      if (obs[s_i].L[k] == 0.0) continue;

      GNSSPre.SNR[k] = obs[s_i].SNR[k] * 0.25;
      GNSSPre.LLI[k] = obs[s_i].LLI[k];
      GNSSPre.P[k] = obs[s_i].P[k];
      GNSSPre.L[k] = obs[s_i].L[k];
      GNSSPre.lamda[k] = nav->lam[obs[s_i].sat - 1][k];
      GNSSPre.code[k] = obs[s_i].code[k];
    }

    /* get azi and ele */
    GNSSPre.azimuth = azel_[0 + s_i * 2];
    GNSSPre.elevation = azel_[1 + s_i * 2];

    double dion, dtrp, vmeas, vion, vtrp, rr[3], pos[3], e[3], P, lam_L1;
    /* psudorange with code bias correction */
    if ((P = prange(obs + s_i, nav, azel_ + s_i * 2, 2, &opt_, &vmeas)) == 0.0)
      continue;

    /* ionospheric corrections */
    for (i = 0; i < 3; i++) rr[i] = sol->rr[i];
    ecef2pos(rr, pos);
    if (!ionocorr(obs[s_i].time, nav, obs[s_i].sat, pos, azel_ + s_i * 2,
                  opt->ionoopt, &dion, &vion))
      continue;

    /* GPS-L1 -> L1/B1 */
    if ((lam_L1 = nav->lam[obs[s_i].sat - 1][0]) > 0.0) {
      dion *= pow(lam_L1 / lam_carr[0], 2);
    }

    /* tropospheric corrections */
    if (!tropcorr(obs[s_i].time, nav, pos, azel_ + s_i * 2, opt->tropopt, &dtrp,
                  &vtrp)) {
      continue;
    }
    GNSSPre.err_tropo = dtrp;
    GNSSPre.err_iono = dion;
    GNSSPre.SatClk = dts[0 + s_i * 2];
    GNSSPre.SatXYZ << rs[0 + s_i * 6], rs[1 + s_i * 6], rs[2 + s_i * 6];
    GNSSPre.SatVel << rs[3 + s_i * 6], rs[4 + s_i * 6], rs[5 + s_i * 6];
    GNSSPre.SatVar = var[i];

    /* remove the satellite clock bias, atmosphere error here */
    GNSSPre.P_crr = P + GNSSPre.SatClk * CLIGHT - dion - dtrp;

    if (GNSSPre.elevation > elemask) mCurPntData[GNSSPre.prn] = GNSSPre;
  }

  /* raim fde */
  if (!stat && n >= 6 && opt->posopt[4]) {
    stat = raim_fde(obs, &n, rs, dts, var, svh, nav, &opt_, sol, azel_, vsat,
                    resp, msg);
  }
  /* estimate receiver velocity with doppler */
  if (stat) estvel(obs, n, rs, dts, nav, &opt_, sol, azel_, vsat);

  if (azel) {
    for (i = 0; i < n * 2; i++) azel[i] = azel_[i];
  }
  if (ssat) {
    for (i = 0; i < MAXSAT; i++) {
      ssat[i].vs = 0;
      ssat[i].azel[0] = ssat[i].azel[1] = 0.0;
      ssat[i].resp[0] = ssat[i].resc[0] = 0.0;
      ssat[i].snr[0] = 0;
    }
    for (i = 0; i < n; i++) {
      ssat[obs[i].sat - 1].azel[0] = azel_[i * 2];
      ssat[obs[i].sat - 1].azel[1] = azel_[1 + i * 2];
      ssat[obs[i].sat - 1].snr[0] = obs[i].SNR[0];
      if (!vsat[i]) continue;
      ssat[obs[i].sat - 1].vs = 1;
      ssat[obs[i].sat - 1].resp[0] = resp[i];
    }
  }
  free(rs);
  free(dts);
  free(var);
  free(azel_);
  free(resp);
  return stat;
}
/* precise positioning preprocess
 *----------------------------------------------- input observation data and
 *navigation message, compute rover position by precise positioning return :
 *status (0:no solution,1:valid solution) notes  : before calling function, base
 *station position rtk->sol.rb[] should be properly set for relative mode except
 *for moving-baseline
 *-----------------------------------------------------------------------------*/
void Preproc_ppp(rtk_t *rtk, obsd_t *obs, int n, const nav_t *nav,
                 map<SatPRN, PreGNSS> &mCurPntData, Vector3 &dr) {
  prcopt_t *opt = &rtk->opt;
  gtime_t time = rtk->sol.time; /* previous epoch */
  int i;
  char msg[128] = "";

  /* rover position by single point positioning */
  if (!Pntpos(obs, &n, nav, &rtk->opt, &rtk->sol, NULL, rtk->ssat, msg,
              mCurPntData)) {
    printf("point pos error (%s)\n", msg);
  }
  if (time.time != 0) rtk->tt = timediff(rtk->sol.time, time);

  /* precise point positioning */
  if (opt->mode >= PMODE_PPP_KINEMA) {
    /* exclude measurements of eclipsing satellite (block IIA) */
    if (rtk->opt.posopt[3]) {
      Testeclipse(obs[0].time, mCurPntData, nav);
    }
    /* earth tides correction */
    if (opt->tidecorr) {
      double x[] = {0, 0, 0}, dr0[3];
      for (int i = 0; i < 3; i++) x[i] = rtk->sol.rr[i];
      tidedisp(gpst2utc(obs[0].time), x, opt->tidecorr == 1 ? 1 : 7, &nav->erp,
               opt->odisp[0], dr0);
      dr = Vector3(dr0[0], dr0[1], dr0[2]);
    }
  }
}
void UdStates(gtime_t CurTime, int CurCount, rtk_t *rtk,
              map<SatPRN, PreGNSS> &mCurPntData, Vector3 &NomXYZ) {
  if (gbFirstEpoch) {
    // Only one position state for static
    /* temporal update of pos */
    // Values
    NomXYZ = Vector3(rtk->sol.rr[0], rtk->sol.rr[1], rtk->sol.rr[2]);
    gNewInitValues.insert(X(CurCount), NomXYZ);
    // Prior facror
    DiagNoise::shared_ptr InitPosNoise = DiagNoise::Variances(
        Vector3(gvInitPosNoise[0], gvInitPosNoise[1], gvInitPosNoise[2]));
    gNewGraph->add(PriorFactor<Vector3>(X(CurCount), NomXYZ, InitPosNoise));

    /* temporal update of receiver clock (s) */
    // Values (no need for clock)

    /* temporal update of tropospheric deley (m) */
    double rr[] = {NomXYZ(0), NomXYZ(1), NomXYZ(2)};
    double pos[3], azel[] = {0.0, PI / 2.0}, ztd, var;
    ecef2pos(rr, pos);
    ztd = sbstropcorr(rtk->sol.time, pos, azel, &var);
    // Values
    gNewInitValues.insert(T(CurCount), ztd);
    // Prior facror
    DiagNoise::shared_ptr InitTrpNoise = DiagNoise::Variances(Vector1(var));
    gNewGraph->add(PriorFactor<double>(T(CurCount), ztd, InitTrpNoise));

    gbFirstEpoch = false;

  } else {
    // /* temporal update of pos */
    // // Prior facror
    // DiagNoise::shared_ptr InitPosNoise = DiagNoise::Variances(
    //     Vector3(gvInitPosNoise[0], gvInitPosNoise[1], gvInitPosNoise[2]));
    // gNewGraph->add(PriorFactor<Vector3>(X(0), gResultValues.at<Vector3>(X(0)),
    //                                     InitPosNoise));
    /* temporal update of tropospheric deley (m) */
    // Values
    gNewInitValues.insert(T(CurCount),
                          gResultValues.at<double>(T(CurCount - 1)));
    // Between facror
    gNewGraph->add(
        BetweenFactor<double>(T(CurCount), T(CurCount - 1), 0.0,
                              DiagNoise::Variances(Vector1(gBetwTrpNoise))));
  }

  /* temporal update of receiver clock (s) */
  gNewInitValues.insert(C(CurCount), rtk->sol.dtr[0] * CLIGHT);
  // Prior facror
  DiagNoise::shared_ptr InitClkNoise =
      DiagNoise::Variances(Vector1(gInitClkNoise));
  gNewGraph->add(
      PriorFactor<double>(C(CurCount), rtk->sol.dtr[0] * CLIGHT, InitClkNoise));

  /* temporal update of phase-bias */
  Udbias_ppp(CurTime, rtk, &gNavs, gResultValues, mCurPntData);
  gbCurBreak = 0;
  for (auto &mpit : mCurPntData) {
    int PRN = mpit.second.prn;

    // cycle-slip or new
    if (mpit.second.phasebreak) {
      gbCurBreak++;

      if (!gbFirstEpoch) {
        gAllGNKey_[PRN] = gAllGNKey_[PRN] + 1;
      }
      // Values
      gNewInitValues.insert(N(gAllGNKey_[PRN]), mpit.second.newBias);
      // Prior facror
      gNewGraph->add(
          PriorFactor<double>(N(gAllGNKey_[PRN]), mpit.second.newBias,
                              DiagNoise::Variances(Vector1(gInitAmNoise))));
    }
  }
}
/* Generate factors data */
void GenFactorData(const PreGNSS &CurPreGNSS, const nav_t *nav,
                   const prcopt_t *opt, CarrierPhaseObs &PhObs,
                   PseudorangeObs &PrObs, SatelliteData &SatData,
                   AntennaData &AntData, double &VarPh, double &VarPr) {
  int i, j, k;

  // PhObs
  PhObs.obsTime_.time = CurPreGNSS.obsTime.time;
  PhObs.obsTime_.sec = CurPreGNSS.obsTime.sec;
  PhObs.solTime_.time = CurPreGNSS.solTime.time;
  PhObs.solTime_.sec = CurPreGNSS.solTime.sec;
  PhObs.sys_ = CurPreGNSS.sys;
  PhObs.prn_ = CurPreGNSS.prn;
  for (i = 0; i < NFREQ; i++) {
    PhObs.Ph_[i] = CurPreGNSS.L[i];
    PhObs.lam_[i] = nav->lam[CurPreGNSS.prn - 1][i];
  }

  // PrObs
  PrObs.obsTime_.time = CurPreGNSS.obsTime.time;
  PrObs.obsTime_.sec = CurPreGNSS.obsTime.sec;
  PrObs.solTime_.time = CurPreGNSS.solTime.time;
  PrObs.solTime_.sec = CurPreGNSS.solTime.sec;
  PrObs.sys_ = CurPreGNSS.sys;
  PrObs.prn_ = CurPreGNSS.prn;
  for (i = 0; i < NFREQ; i++) {
    PrObs.Pr_[i] = CurPreGNSS.P[i];
    PrObs.code_[i] = CurPreGNSS.code[i];
    PrObs.cbias_[i] = nav->cbias[CurPreGNSS.prn - 1][i];
    PrObs.lam_[i] = nav->lam[CurPreGNSS.prn - 1][i];
  }

  // SatData
  SatData.SatClk_ = CurPreGNSS.SatClk;
  SatData.SatXYZ_ = CurPreGNSS.SatXYZ;
  SatData.SatVel_ = CurPreGNSS.SatVel;

  // AntData
  int sat = CurPreGNSS.prn;
  for (i = 0; i < NFREQ; i++) {
    AntData.Antdel_[i] = opt->antdel[0][i];
    for (j = 0; j < 19; j++) {
      AntData.PcvS_var_[i][j] = nav->pcvs[sat - 1].var[i][j];
      AntData.PcvR_var_[i][j] = opt->pcvr[0].var[i][j];
    }
    for (k = 0; k < 3; k++) AntData.PcvR_off_[i][k] = opt->pcvr[0].off[i][k];
  }

  VarPh = gPhNoiseFactor * Varerr(CurPreGNSS.prn, CurPreGNSS.sys,
                                  CurPreGNSS.elevation, 0, 0, opt) +
          SQR(0.01) + CurPreGNSS.SatVar;
  VarPr =
      Varerr(CurPreGNSS.prn, CurPreGNSS.sys, CurPreGNSS.elevation, 0, 1, opt) +
      SQR(0.01) + CurPreGNSS.SatVar;
}
/* process positioning -------------------------------------------------------*/
void Procpos(const prcopt_t *sProcOpt, const solopt_t *sSolutionOption,
             string OutFile) {
  sol_t sol = {{0}};
  rtk_t rtk;
  obsd_t obs[MAXOBS]; /* for rover */
  int i, nobs, n;

  vector<int> vPRNVec;
  int PRN, PreCount = 0, CurCount = 0;

  // 5维状态量：位置、钟差、对流层
  Vector3 NomXYZ, PropXYZ;
  double Clkm, dTrp;
  for (int i = 1; i < NSATGPS + 2; i++)
    gAllGNKey_[i] = gAllGNKey_[i - 1] + 10000;

  // 设置优化器，正片开始
  ISAM2DoglegParams doglegParams;
  ISAM2Params ISAM2Parameters;
  ISAM2Parameters.relinearizeThreshold = 0.1;
  ISAM2Parameters.relinearizeSkip = 100;
  ISAM2 Isam2(ISAM2Parameters);
  BatchFixedLagSmoother BFSmoother;
  // initialize factor gNewGraph
  gNewGraph = new NonlinearFactorGraph();

  ofstream ofXYZ;
  ofXYZ.open(OutFile);
  ofXYZ << std::fixed;

  rtkinit(&rtk, sProcOpt);

  while ((nobs = Inputobs(obs, rtk.sol.stat, sProcOpt)) >= 0) {
    /* exclude satellites */
    for (i = n = 0; i < nobs; i++) {
      if ((satsys(obs[i].sat, NULL) & sProcOpt->navsys) &&
          sProcOpt->exsats[obs[i].sat - 1] != 1)
        obs[n++] = obs[i];
    }
    if (n <= 0) continue;

    /* carrier-phase bias correction */
    if (gNavs.nf > 0) {
      Corr_phase_bias_fcb(obs, n, &gNavs);
    } else if (!strstr(sProcOpt->pppopt, "-DIS_FCB")) {
      Corr_phase_bias_ssr(obs, n, &gNavs);
    }

    /* ppp pre-process */
    map<SatPRN, PreGNSS> mCurPntData;
    Vector3 dr;
    Preproc_ppp(&rtk, obs, n, &gNavs, mCurPntData, dr);

    /* Graph construction */
    /* temporal update of states */
    UdStates(obs[0].time, CurCount, &rtk, mCurPntData, NomXYZ);

    /* temporal update of graph */
    // sat-by-sat
    for (auto &mpit : mCurPntData) {
      PRN = mpit.second.prn;

      CarrierPhaseObs PhObs;
      PseudorangeObs PrObs;
      SatelliteData SatData;
      AntennaData AntData;
      double VarPh, VarPr;
      GenFactorData(mpit.second, &gNavs, &rtk.opt, PhObs, PrObs, SatData,
                    AntData, VarPh, VarPr);

      // Generate Phase factor
      PhaseFactor PhFactor(X(0), C(CurCount), T(CurCount), N(gAllGNKey_[PRN]),
                           PhObs, SatData, dr, AntData,
                           DiagNoise::Variances(Vector1(VarPh)));

      // Generate pseudorange factor
      PseudorangeFactor PrFactor(X(0), C(CurCount), T(CurCount), PrObs, SatData,
                                 dr, AntData,
                                 DiagNoise::Variances(Vector1(VarPr)));

      gNewGraph->add(PhFactor);
      gNewGraph->add(PrFactor);
      // vPRNVec.push_back(PRN);
    }

    // gNewGraph->print("Graph:");
    // gNewInitValues.print("Values:");

    Isam2.update(*gNewGraph, gNewInitValues);
    gResultValues = Isam2.calculateEstimate();
    // BFSmoother.update(*gNewGraph, gNewInitValues);
    // gResultValues = BFSmoother.calculateEstimate();

    // 位置估值及显示
    PropXYZ = gResultValues.at<Vector3>(X(0));
    char strTimeYMDHMS[64];
    time2str(obs[0].time, strTimeYMDHMS, 3);
    ofXYZ << string(strTimeYMDHMS) << " " << PropXYZ.x() - gvNavTrue[0] << " "
          << PropXYZ.y() - gvNavTrue[1] << " " << PropXYZ.z() - gvNavTrue[2]
          << " " << gbCurBreak << endl;

    gNewGraph->resize(0);
    gNewInitValues.clear();
    vPRNVec.clear();

    PreCount = CurCount;
    CurCount++;
  }

  rtkfree(&rtk);
  ofXYZ.close();
}

int main(int argc, char *argv[]) {
  /*********** I. Configuration ***********/
  prcopt_t sProcOpt = prcopt_default;
  solopt_t sSolutionOption = solopt_default;
  filopt_t sFileOpt = {""};
  gtime_t sGTimeStart = {0}, sGTimeEnd = {0};
  double PrcsInterval = 0.0;
  int i;  //, j
  char *archInFiles[MAXFILE] = {NULL};

  sProcOpt.mode = PMODE_KINEMA;
  sProcOpt.navsys = 0;
  sProcOpt.refpos = 1;
  sProcOpt.glomodear = 1;
  sSolutionOption.timef = 0;

  /* load options from configuration file */
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-c") && i + 1 < argc) {
      resetsysopts();
      if (!loadopts(argv[++i], sysopts)) return -1;
      getsysopts(&sProcOpt, &sSolutionOption, &sFileOpt);
    }
  }
  if (!sProcOpt.navsys) {
    sProcOpt.navsys = SYS_GPS | SYS_GLO;
  }

  /* load options from yaml file */
  YAML::Node config = YAML::LoadFile(argv[3]);
  // Start and End Time
  vector<double> vecdata = config["TimeStart"].as<vector<double>>();
  double arTimeStart[] = {vecdata[0], vecdata[1], vecdata[2],   // YMD
                          vecdata[3], vecdata[4], vecdata[5]};  // HMS
  vecdata = config["TimeEnd"].as<vector<double>>();
  double arTimeEnd[] = {vecdata[0], vecdata[1], vecdata[2],   // YMD
                        vecdata[3], vecdata[4], vecdata[5]};  // HMS
  sGTimeStart = epoch2time(arTimeStart);
  sGTimeEnd = epoch2time(arTimeEnd);

  // RINEX files, muli-agent obs-files and nav-file
  string InFilepath = config["InFilepath"].as<string>();

  /* open processing session */
  // 天线文件在这里读取
  if (!Openses(sFileOpt)) return -1;

  // read sp3, clk
  string Sp3File = InFilepath + config["Infile-orbitfile"].as<string>();
  string ClkFile = InFilepath + config["Infile-clkfile"].as<string>();
  Readpreceph(Sp3File.c_str(), ClkFile.c_str(), NULL, &sProcOpt, &gNavs);
  /* read ionosphere data file */
  string IonFile = InFilepath + config["Infile-ionofile"].as<string>();
  if (IonFile[0] != '\0') readtec(IonFile.c_str(), &gNavs, 1);

  /* read erp data */
  string ErpFile = InFilepath + config["Infile-eopfile"].as<string>();
  if (ErpFile[0] != '\0') {
    if (!Readerp(ErpFile.c_str(), &gNavs.erp)) {
      printf("error : no erp data %s", ErpFile.c_str());
    }
  }
  /* read obs and nav data */
  string ObsFile = config["Infile-obsfile"].as<string>();
  string Obspath = InFilepath + ObsFile;
  string NavFile = InFilepath + config["Infile-navfile"].as<string>();
  if (!Readobsnav(sGTimeStart, sGTimeEnd, PrcsInterval, Obspath.c_str(),
                  NavFile.c_str(), &sProcOpt, &gObss, &gNavs, gStas))
    return 0;
  /* read dcb parameters */
  string DCBFile = InFilepath + config["Infile-dcbfile"].as<string>();
  if (DCBFile[0] != '\0') readdcb(DCBFile.c_str(), &gNavs, gStas);

  /* set antenna paramters */
  if (sProcOpt.mode != PMODE_SINGLE) {
    Setpcv(gObss.n > 0 ? gObss.data[0].time : timeget(), &sProcOpt, &gNavs,
           &gRecAntParas, &gSatAntParas, gStas);
  }
  /* read ocean tide loading parameters */
  string BLQFile = InFilepath + config["Infile-blqfile"].as<string>();
  if (sProcOpt.mode > PMODE_SINGLE && BLQFile[0] != '\0') {
    Readotl(&sProcOpt, BLQFile.c_str(), gStas);
  }

  // Output file
  string OutFilepath = config["OutFilepath"].as<string>();
  if (access(OutFilepath.c_str(), 0)) {
    string command = "mkdir -p " + OutFilepath;
    system(command.c_str());
  }
  string str1 = ObsFile.substr(0, ObsFile.find('.'));
  string OutFile =
      OutFilepath + str1 + config["Outfile-file_id"].as<string>() + ".pos";

  gvNavTrue = config["vNavTrue"].as<vector<double>>();
  gvInitPosNoise = config["vInitPosNoise"].as<vector<double>>();
  gInitClkNoise = config["InitClkNoise"].as<double>();
  gInitAmNoise = config["InitAmNoise"].as<double>();
  gBetwTrpNoise = config["BetwTrpNoise"].as<double>();
  gPhNoiseFactor = config["PhNoiseFactor"].as<double>();

  Procpos(&sProcOpt, &sSolutionOption, OutFile);

  return 0;
}
