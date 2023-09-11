#include "PPPSamlib.h"

/* show message and check break ----------------------------------------------*/
int Showmsg(const char *format, ...) {
  va_list arg;
  va_start(arg, format);
  vfprintf(stderr, format, arg);
  va_end(arg);
  fprintf(stderr, "\r");
  return 0;
}
/* show message and check break ----------------------------------------------*/
int Checkbrk(const char *format, ...) {
  va_list arg;
  char buff[1024], *p = buff;
  if (!*format) return Showmsg("");
  va_start(arg, format);
  p += vsprintf(p, format, arg);
  va_end(arg);
  if (*proc_rov && *proc_base)
    sprintf(p, " (%s-%s)", proc_rov, proc_base);
  else if (*proc_rov)
    sprintf(p, " (%s)", proc_rov);
  else if (*proc_base)
    sprintf(p, " (%s)", proc_base);
  return Showmsg(buff);
}
/* read prec ephemeris, sbas data, lex data, tec grid and open rtcm ----------*/
void Readpreceph(const char *Sp3File, const char *ClkFile, char *FcbFile,
                 const prcopt_t *prcopt, nav_t *nav) {
  char *ext;

  nav->ne = nav->nemax = 0;
  nav->nc = nav->ncmax = 0;
  nav->nf = nav->nfmax = 0;

  /* read precise ephemeris files */
  if (Sp3File) readsp3(Sp3File, nav, 0);

  /* read precise clock files */
  if (ClkFile) readrnxc(ClkFile, nav);

  /* read satellite fcb files */
  if (FcbFile) {
    if ((ext = strrchr(FcbFile, '.')) &&
        (!strcmp(ext, ".fcb") || !strcmp(ext, ".FCB"))) {
      readfcb(FcbFile, nav);
    }
  }
}
/* read earth rotation parameters ----------------------------------------------
 * read earth rotation parameters
 * args   : char   *file       I   IGS ERP file (IGS ERP ver.2)
 *          erp_t  *erp        O   earth rotation parameters
 * return : status (1:ok,0:file open error)
 *-----------------------------------------------------------------------------*/
int Readerp(const char *file, erp_t *erp) {
  FILE *fp;
  erpd_t *erp_data;
  double v[14] = {0};
  char buff[256];

  free(erp->data);
  erp->data = NULL;
  erp->n = erp->nmax = 0;

  if (!(fp = fopen(file, "r"))) {
    trace(2, "erp file open error: file=%s\n", file);
    return 0;
  }
  while (fgets(buff, sizeof(buff), fp)) {
    if (sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               v, v + 1, v + 2, v + 3, v + 4, v + 5, v + 6, v + 7, v + 8, v + 9,
               v + 10, v + 11, v + 12, v + 13) < 5) {
      continue;
    }
    if (erp->n >= erp->nmax) {
      erp->nmax = erp->nmax <= 0 ? 128 : erp->nmax * 2;
      erp_data = (erpd_t *)realloc(erp->data, sizeof(erpd_t) * erp->nmax);
      if (!erp_data) {
        free(erp->data);
        erp->data = NULL;
        erp->n = erp->nmax = 0;
        fclose(fp);
        return 0;
      }
      erp->data = erp_data;
    }
    erp->data[erp->n].mjd = v[0];
    erp->data[erp->n].xp = v[1] * 1E-6 * AS2R;
    erp->data[erp->n].yp = v[2] * 1E-6 * AS2R;
    erp->data[erp->n].ut1_utc = v[3] * 1E-7;
    erp->data[erp->n].lod = v[4] * 1E-7;
    erp->data[erp->n].xpr = v[12] * 1E-6 * AS2R;
    erp->data[erp->n++].ypr = v[13] * 1E-6 * AS2R;
  }
  fclose(fp);
  return 1;
}
/* read ocean tide loading parameters ----------------------------------------*/
void Readotl(prcopt_t *sProcOpt, const char *file, const sta_t *sta) {
  int i, mode = PMODE_DGPS <= sProcOpt->mode && sProcOpt->mode <= PMODE_FIXED;

  for (i = 0; i < (mode ? 2 : 1); i++) {
    readblq(file, sta[i].name, sProcOpt->odisp[i]);
  }
}
/* search next observation data
 * ----------------------------------------*/
int Nextobsf(const obs_t *obs, int *i, int rcv) {
  double tt;
  int n;

  for (; *i < obs->n; (*i)++)
    if (obs->data[*i].rcv == rcv) break;
  for (n = 0; *i + n < obs->n; n++) {
    tt = timediff(obs->data[*i + n].time, obs->data[*i].time);
    if (obs->data[*i + n].rcv != rcv || tt > DTTOL) break;
  }
  return n;
}
int Nextobsb(const obs_t *obs, int *i, int rcv) {
  double tt;
  int n;

  for (; *i >= 0; (*i)--)
    if (obs->data[*i].rcv == rcv) break;
  for (n = 0; *i - n >= 0; n++) {
    tt = timediff(obs->data[*i - n].time, obs->data[*i].time);
    if (obs->data[*i - n].rcv != rcv || tt < -DTTOL) break;
  }
  return n;
}
/* read obs and nav data -----------------------------------------------------*/
int Readobsnav(gtime_t sGTimeStart, gtime_t sGTimeEnd, double PrcsInterval,
               const char *ObsFile, const char *NavFile, const prcopt_t *prcopt,
               obs_t *obs, nav_t *nav, sta_t *sta) {
  int i, j, ind = 0, nobs = 0, rcv = 1;

  obs->data = NULL;
  obs->n = obs->nmax = 0;
  nav->eph = NULL;
  nav->n = nav->nmax = 0;
  nav->geph = NULL;
  nav->ng = nav->ngmax = 0;
  nav->seph = NULL;
  nav->ns = nav->nsmax = 0;
  nepoch = 0;

  /* read rinex obs file */
  if (readrnxt(ObsFile, rcv, sGTimeStart, sGTimeEnd, PrcsInterval,
               prcopt->rnxopt[rcv <= 1 ? 0 : 1], obs, NULL, sta) < 0) {
    Checkbrk("error : insufficient memory");
    return 0;
  }

  /* read rinex obs and nav file */
  if (readrnxt(NavFile, rcv, sGTimeStart, sGTimeEnd, PrcsInterval,
               prcopt->rnxopt[rcv <= 1 ? 0 : 1], NULL, nav, NULL) < 0) {
    Checkbrk("error : insufficient memory");
    return 0;
  }

  if (obs->n <= 0) {
    Checkbrk("error : no obs data");
    return 0;
  }
  if ((nav->n <= 0 && nav->ng <= 0 && nav->ns <= 0)) {
    Checkbrk("error : no nav data");
    return 0;
  }
  /* sort observation data */
  nepoch = sortobs(obs);

  /* delete duplicated ephemeris */
  uniqnav(nav);

  /* set time span for progress display */
  if (sGTimeStart.time == 0 || sGTimeEnd.time == 0) {
    for (i = 0; i < obs->n; i++)
      if (obs->data[i].rcv == 1) break;
    for (j = obs->n - 1; j >= 0; j--)
      if (obs->data[j].rcv == 1) break;
    if (i < j) {
      if (sGTimeStart.time == 0) sGTimeStart = obs->data[i].time;
      if (sGTimeEnd.time == 0) sGTimeEnd = obs->data[j].time;
    }
  }
  return 1;
}
/* carrier-phase bias correction by fcb --------------------------------------*/
void Corr_phase_bias_fcb(obsd_t *obs, int n, const nav_t *nav) {
  int i, j, k;

  for (i = 0; i < nav->nf; i++) {
    if (timediff(nav->fcb[i].sGTimeEnd, obs[0].time) < -1E-3) continue;
    if (timediff(nav->fcb[i].sGTimeStart, obs[0].time) > 1E-3) break;
    for (j = 0; j < n; j++) {
      for (k = 0; k < NFREQ; k++) {
        if (obs[j].L[k] == 0.0) continue;
        obs[j].L[k] -= nav->fcb[i].bias[obs[j].sat - 1][k];
      }
    }
    return;
  }
}
/* carrier-phase bias correction by ssr --------------------------------------*/
void Corr_phase_bias_ssr(obsd_t *obs, int n, const nav_t *nav) {
  double lam;
  int i, j, code;

  for (i = 0; i < n; i++)
    for (j = 0; j < NFREQ; j++) {
      if (!(code = obs[i].code[j])) continue;
      if ((lam = nav->lam[obs[i].sat - 1][j]) == 0.0) continue;

      /* correct phase bias (cyc) */
      obs[i].L[j] -= nav->ssr[obs[i].sat - 1].pbias[code - 1] / lam;
    }
}
/* geometry-free phase measurement -------------------------------------------*/
double Gfmeas(const PreGNSS &CurPreGNSS, const nav_t *nav) {
  const double *lam = nav->lam[CurPreGNSS.prn - 1];
  int i = (satsys(CurPreGNSS.prn, NULL) & (SYS_GAL | SYS_SBS)) ? 2 : 1;

  if (lam[0] == 0.0 || lam[i] == 0.0 || CurPreGNSS.L[0] == 0.0 ||
      CurPreGNSS.L[i] == 0.0)
    return 0.0;
  return lam[0] * CurPreGNSS.L[0] - lam[i] * CurPreGNSS.L[i];
}
/* Melbourne-Wubbena linear combination --------------------------------------*/
double Mwmeas(const PreGNSS &CurPreGNSS, const nav_t *nav) {
  const double *lam = nav->lam[CurPreGNSS.prn - 1];
  int i = (satsys(CurPreGNSS.prn, NULL) & (SYS_GAL | SYS_SBS)) ? 2 : 1;

  if (lam[0] == 0.0 || lam[i] == 0.0 || CurPreGNSS.L[0] == 0.0 ||
      CurPreGNSS.L[i] == 0.0 || CurPreGNSS.P[0] == 0.0 ||
      CurPreGNSS.P[i] == 0.0)
    return 0.0;
  return lam[0] * lam[i] * (CurPreGNSS.L[0] - CurPreGNSS.L[i]) /
             (lam[i] - lam[0]) -
         (lam[i] * CurPreGNSS.P[0] + lam[0] * CurPreGNSS.P[i]) /
             (lam[i] + lam[0]);
}
/* detect cycle slip by LLI --------------------------------------------------*/
void Detslp_ll(rtk_t *rtk, const map<SatPRN, PreGNSS> &mCurPntData) {
  int i, j;

  for (auto mpit : mCurPntData)
    for (j = 0; j < rtk->opt.nf; j++) {
      if (mpit.second.L[j] == 0.0 || !(mpit.second.LLI[j] & 3)) continue;

      rtk->ssat[mpit.second.prn - 1].slip[j] = 1;
    }
}
/* detect cycle slip by geometry free phase jump -----------------------------*/
void Detslp_gf(rtk_t *rtk, const map<SatPRN, PreGNSS> &mCurPntData,
               const nav_t *nav) {
  double g0, g1;
  int i, j;

  for (auto mpit : mCurPntData) {
    if ((g1 = Gfmeas(mpit.second, nav)) == 0.0) continue;

    g0 = rtk->ssat[mpit.second.prn - 1].gf;
    rtk->ssat[mpit.second.prn - 1].gf = g1;

    if (g0 != 0.0 && fabs(g1 - g0) > rtk->opt.thresslip) {
      for (j = 0; j < rtk->opt.nf; j++)
        rtk->ssat[mpit.second.prn - 1].slip[j] |= 1;
    }
  }
}
/* detect slip by Melbourne-Wubbena linear combination jump ------------------*/
void Detslp_mw(rtk_t *rtk, const map<SatPRN, PreGNSS> &mCurPntData,
               const nav_t *nav) {
  double w0, w1;
  int i, j;

  for (auto mpit : mCurPntData) {
    if ((w1 = Mwmeas(mpit.second, nav)) == 0.0) continue;

    w0 = rtk->ssat[mpit.second.prn - 1].mw;
    rtk->ssat[mpit.second.prn - 1].mw = w1;

    /*double b0 = w0*(FREQ1 - FREQ2) / CLIGHT;
    double b1 = w1*(FREQ1 - FREQ2) / CLIGHT;
    double deltb = b1 - b0;*/

    if (w0 != 0.0 && fabs(w1 - w0) > THRES_MW_JUMP) {
      for (j = 0; j < rtk->opt.nf; j++)
        rtk->ssat[mpit.second.prn - 1].slip[j] |= 1;
    }
  }
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
void Testeclipse(const gtime_t &CurTime, map<SatPRN, PreGNSS> &mCurPntData,
                 const nav_t *nav) {
  double rsun[3], esun[3], r, ang, erpv[5] = {0}, cosa;
  const char *type;

  /* unit vector of sun direction (ecef) */
  sunmoonpos(gpst2utc(CurTime), erpv, rsun, NULL, NULL);
  normv3(rsun, esun);

  for (auto &mpit : mCurPntData) {
    type = nav->pcvs[mpit.second.prn - 1].type;

    Vector3 rs = mpit.second.SatXYZ;
    if ((r = rs.norm()) <= 0.0) continue;

    /* only block IIA */
    if (*type && !strstr(type, "BLOCK IIA")) continue;

    /* sun-earth-satellite angle */
    cosa = rs.dot(Vector3(esun[0], esun[1], esun[2])) / r;
    cosa = cosa < -1.0 ? -1.0 : (cosa > 1.0 ? 1.0 : cosa);
    ang = acos(cosa);

    /* test eclipse */
    if (ang < PI / 2.0 || r * sin(ang) > RE_WGS84) continue;

    mpit.second.SatXYZ = Vector3::Zero();
  }
}
/* antenna corrected measurements --------------------------------------------*/
void Corr_meas(const PreGNSS &CurPreGNSS, const nav_t *nav, const double *azel,
               const prcopt_t *opt, const double *dantr, const double *dants,
               double phw, double *L, double *P, double *Lc, double *Pc) {
  const double *lam = nav->lam[CurPreGNSS.prn - 1];
  double C1, C2;
  int i, sys;

  for (i = 0; i < NFREQ; i++) {
    L[i] = P[i] = 0.0;
    if (lam[i] == 0.0 || CurPreGNSS.L[i] == 0.0 || CurPreGNSS.P[i] == 0.0)
      continue;
    if (testsnr(0, 0, azel[1], CurPreGNSS.SNR[i] * 0.25, &opt->snrmask))
      continue;

    /* antenna phase center and phase windup correction */
    L[i] = CurPreGNSS.L[i] * lam[i] - dants[i] - dantr[i] - phw * lam[i];
    P[i] = CurPreGNSS.P[i] - dants[i] - dantr[i];

    /* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
    if (CurPreGNSS.code[i] == CODE_L1C) {
      P[i] += nav->cbias[CurPreGNSS.prn - 1][1];
    } else if (CurPreGNSS.code[i] == CODE_L2C ||
               CurPreGNSS.code[i] == CODE_L2X ||
               CurPreGNSS.code[i] == CODE_L2L ||
               CurPreGNSS.code[i] == CODE_L2S) {
      P[i] += nav->cbias[CurPreGNSS.prn - 1][2];
#if 0
			L[i]-=0.25*lam[i]; /* 1/4 cycle-shift */
#endif
    }
  }
  /* iono-free LC */
  *Lc = *Pc = 0.0;
  sys = satsys(CurPreGNSS.prn, NULL);
  i = (sys & (SYS_GAL | SYS_SBS)) ? 2 : 1; /* L1/L2 or L1/L5 */
  if (lam[0] == 0.0 || lam[i] == 0.0) return;

  C1 = SQR(lam[i]) / (SQR(lam[i]) - SQR(lam[0]));  /*  f1^2/(f1^2-fi^2) */
  C2 = -SQR(lam[0]) / (SQR(lam[i]) - SQR(lam[0])); /* -fi^2/(f1^2-fi^2) */

#if 1
  /* P1-P2 dcb correction (P1->Pc,P2->Pc) */
  if (sys & (SYS_GPS | SYS_GLO | SYS_QZS)) {
    if (P[0] != 0.0) P[0] -= C2 * nav->cbias[CurPreGNSS.prn - 1][0];
    if (P[1] != 0.0) P[1] += C1 * nav->cbias[CurPreGNSS.prn - 1][0];
  }
#endif
  if (L[0] != 0.0 && L[i] != 0.0) *Lc = C1 * L[0] + C2 * L[i];
  if (P[0] != 0.0 && P[i] != 0.0) *Pc = C1 * P[0] + C2 * P[i];
}