/*------------------------------------------------------------------------------
 * postpos.c : post-processing positioning
 *
 *          Copyright (C) 2007-2014 by T.TAKASU, All rights reserved.
 *
 * version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
 * history : 2007/05/08  1.0  new
 *           2008/06/16  1.1  support binary inputs
 *           2009/01/02  1.2  support new rtk positioing api
 *           2009/09/03  1.3  fix bug on combined mode of moving-baseline
 *           2009/12/04  1.4  fix bug on obs data buffer overflow
 *           2010/07/26  1.5  support ppp-kinematic and ppp-static
 *                            support multiple sessions
 *                            support sbas positioning
 *                            changed api:
 *                                postpos()
 *                            deleted api:
 *                                postposopt()
 *           2010/08/16  1.6  fix bug sbas message synchronization (2.4.0_p4)
 *           2010/12/09  1.7  support qzss lex and ssr corrections
 *           2011/02/07  1.8  fix bug on sbas navigation data conflict
 *           2011/03/22  1.9  add function reading g_tec file
 *           2011/08/20  1.10 fix bug on freez if solstatic=single and combined
 *           2011/09/15  1.11 add function reading stec file
 *           2012/02/01  1.12 support keyword expansion of rtcm ssr corrections
 *           2013/03/11  1.13 add function reading otl and erp data
 *           2014/06/29  1.14 fix problem on overflow of # of satellites
 *           2015/03/23  1.15 fix bug on ant type replacement by rinex header
 *                            fix bug on combined filter for moving-base mode
 *           2015/04/29  1.16 fix bug on reading rtcm ssr corrections
 *                            add function to read satellite fcb
 *                            add function to read stec and troposphere file
 *                            add keyword replacement in dcb, erp and ionos file
 *           2015/11/13  1.17 add support of L5 antenna phase center paramters
 *                            add *.stec and *.trp file for ppp correction
 *           2015/11/26  1.18 support opt->freqopt(disable L2)
 *           2016/01/12  1.19 add carrier-phase bias correction by ssr
 *           2016/07/31  1.20 fix error message problem in rnx2rtkp
 *           2016/08/29  1.21 suppress warnings
 *-----------------------------------------------------------------------------*/
#include "rtklib.h"

static const char rcsid[] =
    "$Id: postpos.c,v 1.1 2008/07/17 21:48:06 ttaka Exp $";

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define SQRT(x) ((x) <= 0.0 ? 0.0 : sqrt(x))

#define MAXPRCDAYS 100 /* max days of continuous processing */
#define MAXINFILE 1000 /* max number of input files */

/* constants/global variables ------------------------------------------------*/

static pcvs_t RecAntParas = {0}; /* receiver antenna parameters */
static pcvs_t SatAntParas = {0}; /* satellite antenna parameters */
static obs_t obss = {0};         /* observation data */
static nav_t navs = {0};     /* navigation data */
static sbs_t sbss = {0};         /* sbas messages */
static lex_t lexs = {0};         /* lex messages */
static sta_t stas[MAXRCV];       /* station infomation */
static int nepoch = 0;           /* number of observation epochs */
static int iobsu = 0;  /* current rover observation data arIndex_InFiles */
static int iobsr = 0;  /* current reference observation data arIndex_InFiles */
static int isbs = 0;   /* current sbas message arIndex_InFiles */
static int ilex = 0;   /* current lex message arIndex_InFiles */
static int revs = 0;   /* analysis direction (0:forward,1:backward) */
static int aborts = 0; /* abort status */
static sol_t *solf;    /* forward solutions */
static sol_t *solb;    /* backward solutions */
static double *rbf;    /* forward base positions */
static double *rbb;    /* backward base positions */
static int isolf = 0;  /* current forward solutions arIndex_InFiles */
static int isolb = 0;  /* current backward solutions arIndex_InFiles */
static char proc_rov[64] = "";    /* rover for current processing */
static char proc_base[64] = "";   /* base station for current processing */
static char rtcm_file[1024] = ""; /* rtcm data file */
static char rtcm_path[1024] = ""; /* rtcm data path */
static rtcm_t rtcm;               /* rtcm control struct */
static FILE *fp_rtcm = NULL;      /* rtcm data file pointer */

static int showmsg(char *format, ...) {
  va_list arg;
  va_start(arg, format);
  vfprintf(stderr, format, arg);
  va_end(arg);
  fprintf(stderr, "\r");
  return 0;
}

/* show message and check break ----------------------------------------------*/
static int checkbrk(const char *format, ...) {
  va_list arg;
  char buff[1024], *p = buff;
  if (!*format) return printf("");
  va_start(arg, format);
  p += vsprintf(p, format, arg);
  va_end(arg);
  if (*proc_rov && *proc_base)
    sprintf(p, " (%s-%s)", proc_rov, proc_base);
  else if (*proc_rov)
    sprintf(p, " (%s)", proc_rov);
  else if (*proc_base)
    sprintf(p, " (%s)", proc_base);
  return showmsg(buff);
}
/* output reference position -------------------------------------------------*/
static void outrpos(FILE *fp, const double *r, const solopt_t *opt) {
  double pos[3], dms1[3], dms2[3];
  const char *sep = opt->sep;

  trace(3, "outrpos :\n");

  if (opt->posf == SOLF_LLH || opt->posf == SOLF_ENU) {
    ecef2pos(r, pos);
    if (opt->degf) {
      deg2dms(pos[0] * R2D, dms1, 5);
      deg2dms(pos[1] * R2D, dms2, 5);
      fprintf(fp, "%3.0f%s%02.0f%s%08.5f%s%4.0f%s%02.0f%s%08.5f%s%10.4f",
              dms1[0], sep, dms1[1], sep, dms1[2], sep, dms2[0], sep, dms2[1],
              sep, dms2[2], sep, pos[2]);
    } else {
      fprintf(fp, "%13.9f%s%14.9f%s%10.4f", pos[0] * R2D, sep, pos[1] * R2D,
              sep, pos[2]);
    }
  } else if (opt->posf == SOLF_XYZ) {
    fprintf(fp, "%14.4f%s%14.4f%s%14.4f", r[0], sep, r[1], sep, r[2]);
  }
}
/* output header -------------------------------------------------------------*/
static void outheader(FILE *fp, char **file, int NumInFiles,
                      const prcopt_t *sProcessOption,
                      const solopt_t *sSolutionOption) {
  const char *s1[] = {"GPST", "UTC", "JST"};
  gtime_t sGTimeStart, sGTimeEnd;
  double t1, t2;
  int i, j, w1, w2;
  char s2[32], s3[32];

  trace(3, "outheader: NumInFiles=%d\n", NumInFiles);

  if (sSolutionOption->posf == SOLF_NMEA ||
      sSolutionOption->posf == SOLF_STAT) {
    return;
  }
  if (sSolutionOption->outhead) {
    if (!*sSolutionOption->prog) {
      fprintf(fp, "%s program   : RTKLIB ver.%s\n", COMMENTH, VER_RTKLIB);
    } else {
      fprintf(fp, "%s program   : %s\n", COMMENTH, sSolutionOption->prog);
    }
    for (i = 0; i < NumInFiles; i++) {
      fprintf(fp, "%s inp file  : %s\n", COMMENTH, file[i]);
    }
    for (i = 0; i < obss.n; i++)
      if (obss.data[i].rcv == 1) break;
    for (j = obss.n - 1; j >= 0; j--)
      if (obss.data[j].rcv == 1) break;
    if (j < i) {
      fprintf(fp, "\n%s no rover obs data\n", COMMENTH);
      return;
    }
    sGTimeStart = obss.data[i].time;
    sGTimeEnd = obss.data[j].time;
    t1 = time2gpst(sGTimeStart, &w1);
    t2 = time2gpst(sGTimeEnd, &w2);
    if (sSolutionOption->times >= 1) sGTimeStart = gpst2utc(sGTimeStart);
    if (sSolutionOption->times >= 1) sGTimeEnd = gpst2utc(sGTimeEnd);
    if (sSolutionOption->times == 2)
      sGTimeStart = timeadd(sGTimeStart, 9 * 3600.0);
    if (sSolutionOption->times == 2) sGTimeEnd = timeadd(sGTimeEnd, 9 * 3600.0);
    time2str(sGTimeStart, s2, 1);
    time2str(sGTimeEnd, s3, 1);
    fprintf(fp, "%s obs start : %s %s (week%04d %8.1fs)\n", COMMENTH, s2,
            s1[sSolutionOption->times], w1, t1);
    fprintf(fp, "%s obs end   : %s %s (week%04d %8.1fs)\n", COMMENTH, s3,
            s1[sSolutionOption->times], w2, t2);
  }
  if (sSolutionOption->outopt) {
    outprcopt(fp, sProcessOption);
  }
  if (PMODE_DGPS <= sProcessOption->mode &&
      sProcessOption->mode <= PMODE_FIXED &&
      sProcessOption->mode != PMODE_MOVEB) {
    fprintf(fp, "%s ref pos   :", COMMENTH);
    outrpos(fp, sProcessOption->rb, sSolutionOption);
    fprintf(fp, "\n");
  }
  if (sSolutionOption->outhead || sSolutionOption->outopt)
    fprintf(fp, "%s\n", COMMENTH);

  outsolhead(fp, sSolutionOption);
}
/* search next observation data arIndex_InFiles
 * ----------------------------------------*/
static int nextobsf(const obs_t *obs, int *i, int rcv) {
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
static int nextobsb(const obs_t *obs, int *i, int rcv) {
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
/* update rtcm ssr correction ------------------------------------------------*/
static void update_rtcm_ssr(gtime_t time) {
  char path[1024];
  int i;

  /* open or swap rtcm file */
  reppath(rtcm_file, path, time, "", "");

  if (strcmp(path, rtcm_path)) {
    strcpy(rtcm_path, path);

    if (fp_rtcm) fclose(fp_rtcm);
    fp_rtcm = fopen(path, "rb");
    if (fp_rtcm) {
      rtcm.time = time;
      input_rtcm3f(&rtcm, fp_rtcm);
      trace(2, "rtcm file open: %s\n", path);
    }
  }
  if (!fp_rtcm) return;

  /* read rtcm file until current time */
  while (timediff(rtcm.time, time) < 1E-3) {
    if (input_rtcm3f(&rtcm, fp_rtcm) < -1) break;

    /* update ssr corrections */
    for (i = 0; i < MAXSAT; i++) {
      if (!rtcm.ssr[i].update || rtcm.ssr[i].iod[0] != rtcm.ssr[i].iod[1] ||
          timediff(time, rtcm.ssr[i].t0[0]) < -1E-3)
        continue;
      navs.ssr[i] = rtcm.ssr[i];
      rtcm.ssr[i].update = 0;
    }
  }
}
/* input obs data, navigation messages and sbas correction -------------------*/
static int inputobs(obsd_t *obs, int solq, const prcopt_t *sProcessOption) {
  gtime_t time = {0};
  int i, nu, nr, n = 0;

  trace(3, "infunc  : revs=%d iobsu=%d iobsr=%d isbs=%d\n", revs, iobsu, iobsr,
        isbs);

  if (0 <= iobsu && iobsu < obss.n) {
    time = obss.data[iobsu].time;
    if (checkbrk("processing : %s Q=%d", time_str(time, 0), solq)) {
      aborts = 1;
      printf("aborted");
      return -1;
    }
  }
  if (!revs) { /* input forward data */
    if ((nu = nextobsf(&obss, &iobsu, 1)) <= 0) return -1;
    if (sProcessOption->intpref) {
      for (; (nr = nextobsf(&obss, &iobsr, 2)) > 0; iobsr += nr)
        if (timediff(obss.data[iobsr].time, obss.data[iobsu].time) > -DTTOL)
          break;
    } else {
      for (i = iobsr; (nr = nextobsf(&obss, &i, 2)) > 0; iobsr = i, i += nr)
        if (timediff(obss.data[i].time, obss.data[iobsu].time) > DTTOL) break;
    }
    nr = nextobsf(&obss, &iobsr, 2);
    for (i = 0; i < nu && n < MAXOBS * 2; i++) obs[n++] = obss.data[iobsu + i];
    for (i = 0; i < nr && n < MAXOBS * 2; i++) obs[n++] = obss.data[iobsr + i];
    iobsu += nu;

    /* update sbas corrections */
    while (isbs < sbss.n) {
      time = gpst2time(sbss.msgs[isbs].week, sbss.msgs[isbs].tow);

      if (getbitu(sbss.msgs[isbs].msg, 8, 6) != 9) { /* except for geo nav */
        sbsupdatecorr(sbss.msgs + isbs, &navs);
      }
      if (timediff(time, obs[0].time) > -1.0 - DTTOL) break;
      isbs++;
    }
    /* update lex corrections */
    while (ilex < lexs.n) {
      if (lexupdatecorr(lexs.msgs + ilex, &navs, &time)) {
        if (timediff(time, obs[0].time) > -1.0 - DTTOL) break;
      }
      ilex++;
    }
    /* update rtcm ssr corrections */
    if (*rtcm_file) {
      update_rtcm_ssr(obs[0].time);
    }
  } else { /* input backward data */
    if ((nu = nextobsb(&obss, &iobsu, 1)) <= 0) return -1;
    if (sProcessOption->intpref) {
      for (; (nr = nextobsb(&obss, &iobsr, 2)) > 0; iobsr -= nr)
        if (timediff(obss.data[iobsr].time, obss.data[iobsu].time) < DTTOL)
          break;
    } else {
      for (i = iobsr; (nr = nextobsb(&obss, &i, 2)) > 0; iobsr = i, i -= nr)
        if (timediff(obss.data[i].time, obss.data[iobsu].time) < -DTTOL) break;
    }
    nr = nextobsb(&obss, &iobsr, 2);
    for (i = 0; i < nu && n < MAXOBS * 2; i++)
      obs[n++] = obss.data[iobsu - nu + 1 + i];
    for (i = 0; i < nr && n < MAXOBS * 2; i++)
      obs[n++] = obss.data[iobsr - nr + 1 + i];
    iobsu -= nu;

    /* update sbas corrections */
    while (isbs >= 0) {
      time = gpst2time(sbss.msgs[isbs].week, sbss.msgs[isbs].tow);

      if (getbitu(sbss.msgs[isbs].msg, 8, 6) != 9) { /* except for geo nav */
        sbsupdatecorr(sbss.msgs + isbs, &navs);
      }
      if (timediff(time, obs[0].time) < 1.0 + DTTOL) break;
      isbs--;
    }
    /* update lex corrections */
    while (ilex >= 0) {
      if (lexupdatecorr(lexs.msgs + ilex, &navs, &time)) {
        if (timediff(time, obs[0].time) < 1.0 + DTTOL) break;
      }
      ilex--;
    }
  }
  return n;
}
/* carrier-phase bias correction by fcb --------------------------------------*/
static void corr_phase_bias_fcb(obsd_t *obs, int n, const nav_t *nav) {
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
static void corr_phase_bias_ssr(obsd_t *obs, int n, const nav_t *nav) {
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
/* process positioning -------------------------------------------------------*/
static void procpos(FILE *fp, const prcopt_t *sProcessOption,
                    const solopt_t *sSolutionOption, int mode) {
  gtime_t time = {0};
  sol_t sol = {{0}};
  rtk_t rtk;
  obsd_t obs[MAXOBS * 2]; /* for rover and base */
  double rb[3] = {0};
  int i, nobs, n, solstatic, pri[] = {0, 1, 2, 3, 4, 5, 1, 6};

  trace(3, "procpos : mode=%d\n", mode);

  solstatic =
      sSolutionOption->solstatic && (sProcessOption->mode == PMODE_STATIC ||
                                     sProcessOption->mode == PMODE_PPP_STATIC);

  rtkinit(&rtk, sProcessOption);
  rtcm_path[0] = '\0';

  while ((nobs = inputobs(obs, rtk.sol.stat, sProcessOption)) >= 0) {
    char pause;
    if (!strcmp(
            time_str(obs[0].time, 2),
            "2019/04/20 00:38:16.00"))  //||			!rtk.sol.stat
    {
      printf("error next!");
      scanf("%c", &pause);
    }
    /* exclude satellites */
    for (i = n = 0; i < nobs; i++) {
      if ((satsys(obs[i].sat, NULL) & sProcessOption->navsys) &&
          sProcessOption->exsats[obs[i].sat - 1] != 1)
        obs[n++] = obs[i];
    }
    if (n <= 0) continue;

    /* carrier-phase bias correction */
    if (navs.nf > 0) {
      corr_phase_bias_fcb(obs, n, &navs);
    } else if (!strstr(sProcessOption->pppopt, "-DIS_FCB")) {
      corr_phase_bias_ssr(obs, n, &navs);
    }
    /* disable L2 */
#if 0
		if (sProcessOption->freqopt==1) {
			for (i=0;i<n;i++) obs[i].L[1]=obs[i].P[1]=0.0;
		}
#endif
    if (!rtkpos(&rtk, obs, n, &navs)) continue;

    if (mode == 0) { /* forward/backward */
      if (!solstatic) {
        outsol(fp, &rtk.sol, rtk.rb, sSolutionOption);
      } else if (time.time == 0 || pri[rtk.sol.stat] <= pri[sol.stat]) {
        sol = rtk.sol;
        for (i = 0; i < 3; i++) rb[i] = rtk.rb[i];
        if (time.time == 0 || timediff(rtk.sol.time, time) < 0.0) {
          time = rtk.sol.time;
        }
      }
    } else if (!revs) { /* combined-forward */
      if (isolf >= nepoch) return;
      solf[isolf] = rtk.sol;
      for (i = 0; i < 3; i++) rbf[i + isolf * 3] = rtk.rb[i];
      isolf++;
    } else { /* combined-backward */
      if (isolb >= nepoch) return;
      solb[isolb] = rtk.sol;
      for (i = 0; i < 3; i++) rbb[i + isolb * 3] = rtk.rb[i];
      isolb++;
    }
  }
  if (mode == 0 && solstatic && time.time != 0.0) {
    sol.time = time;
    outsol(fp, &sol, rb, sSolutionOption);
  }
  rtkfree(&rtk);
}
/* validation of combined solutions ------------------------------------------*/
static int valcomb(const sol_t *solf, const sol_t *solb) {
  double dr[3], var[3];
  int i;
  char tstr[32];

  trace(3, "valcomb :\n");

  /* compare forward and backward solution */
  for (i = 0; i < 3; i++) {
    dr[i] = solf->rr[i] - solb->rr[i];
    var[i] = solf->qr[i] + solb->qr[i];
  }
  for (i = 0; i < 3; i++) {
    if (dr[i] * dr[i] <= 16.0 * var[i]) continue; /* ok if in 4-sigma */

    time2str(solf->time, tstr, 2);
    trace(2, "degrade fix to float: %s dr=%.3f %.3f %.3f std=%.3f %.3f %.3f\n",
          tstr + 11, dr[0], dr[1], dr[2], SQRT(var[0]), SQRT(var[1]),
          SQRT(var[2]));
    return 0;
  }
  return 1;
}
/* combine forward/backward solutions and output results ---------------------*/
static void combres(FILE *fp, const prcopt_t *sProcessOption,
                    const solopt_t *sSolutionOption) {
  gtime_t time = {0};
  sol_t sols = {{0}}, sol = {{0}};
  double tt, Qf[9], Qb[9], Qs[9], rbs[3] = {0}, rb[3] = {0}, rr_f[3], rr_b[3],
                                  rr_s[3];
  int i, j, k, solstatic, pri[] = {0, 1, 2, 3, 4, 5, 1, 6};

  trace(3, "combres : isolf=%d isolb=%d\n", isolf, isolb);

  solstatic =
      sSolutionOption->solstatic && (sProcessOption->mode == PMODE_STATIC ||
                                     sProcessOption->mode == PMODE_PPP_STATIC);

  for (i = 0, j = isolb - 1; i < isolf && j >= 0; i++, j--) {
    if ((tt = timediff(solf[i].time, solb[j].time)) < -DTTOL) {
      sols = solf[i];
      for (k = 0; k < 3; k++) rbs[k] = rbf[k + i * 3];
      j++;
    } else if (tt > DTTOL) {
      sols = solb[j];
      for (k = 0; k < 3; k++) rbs[k] = rbb[k + j * 3];
      i--;
    } else if (solf[i].stat < solb[j].stat) {
      sols = solf[i];
      for (k = 0; k < 3; k++) rbs[k] = rbf[k + i * 3];
    } else if (solf[i].stat > solb[j].stat) {
      sols = solb[j];
      for (k = 0; k < 3; k++) rbs[k] = rbb[k + j * 3];
    } else {
      sols = solf[i];
      sols.time = timeadd(sols.time, -tt / 2.0);

      if ((sProcessOption->mode == PMODE_KINEMA ||
           sProcessOption->mode == PMODE_MOVEB) &&
          sols.stat == SOLQ_FIX) {
        /* degrade fix to float if validation failed */
        if (!valcomb(solf + i, solb + j)) sols.stat = SOLQ_FLOAT;
      }
      for (k = 0; k < 3; k++) {
        Qf[k + k * 3] = solf[i].qr[k];
        Qb[k + k * 3] = solb[j].qr[k];
      }
      Qf[1] = Qf[3] = solf[i].qr[3];
      Qf[5] = Qf[7] = solf[i].qr[4];
      Qf[2] = Qf[6] = solf[i].qr[5];
      Qb[1] = Qb[3] = solb[j].qr[3];
      Qb[5] = Qb[7] = solb[j].qr[4];
      Qb[2] = Qb[6] = solb[j].qr[5];

      if (sProcessOption->mode == PMODE_MOVEB) {
        for (k = 0; k < 3; k++) rr_f[k] = solf[i].rr[k] - rbf[k + i * 3];
        for (k = 0; k < 3; k++) rr_b[k] = solb[j].rr[k] - rbb[k + j * 3];
        if (smoother(rr_f, Qf, rr_b, Qb, 3, rr_s, Qs)) continue;
        for (k = 0; k < 3; k++) sols.rr[k] = rbs[k] + rr_s[k];
      } else {
        if (smoother(solf[i].rr, Qf, solb[j].rr, Qb, 3, sols.rr, Qs)) continue;
      }
      sols.qr[0] = (float)Qs[0];
      sols.qr[1] = (float)Qs[4];
      sols.qr[2] = (float)Qs[8];
      sols.qr[3] = (float)Qs[1];
      sols.qr[4] = (float)Qs[5];
      sols.qr[5] = (float)Qs[2];
    }
    if (!solstatic) {
      outsol(fp, &sols, rbs, sSolutionOption);
    } else if (time.time == 0 || pri[sols.stat] <= pri[sol.stat]) {
      sol = sols;
      for (k = 0; k < 3; k++) rb[k] = rbs[k];
      if (time.time == 0 || timediff(sols.time, time) < 0.0) {
        time = sols.time;
      }
    }
  }
  if (solstatic && time.time != 0.0) {
    sol.time = time;
    outsol(fp, &sol, rb, sSolutionOption);
  }
}
/* read prec ephemeris, sbas data, lex data, tec grid and open rtcm ----------*/
static void readpreceph(char **archInFiles, int n, const prcopt_t *prcopt,
                        nav_t *nav, sbs_t *sbs, lex_t *lex) {
  seph_t seph0 = {0};
  int i;
  char *ext;

  trace(2, "readpreceph: n=%d\n", n);

  nav->ne = nav->nemax = 0;
  nav->nc = nav->ncmax = 0;
  nav->nf = nav->nfmax = 0;
  sbs->n = sbs->nmax = 0;
  lex->n = lex->nmax = 0;

  /* read precise ephemeris files */
  for (i = 0; i < n; i++) {
    if (strstr(archInFiles[i], "%r") || strstr(archInFiles[i], "%b")) continue;
    readsp3(archInFiles[i], nav, 0);
  }
  /* read precise clock files */
  for (i = 0; i < n; i++) {
    if (strstr(archInFiles[i], "%r") || strstr(archInFiles[i], "%b")) continue;
    readrnxc(archInFiles[i], nav);
  }
  /* read satellite fcb files */
  for (i = 0; i < n; i++) {
    if (strstr(archInFiles[i], "%r") || strstr(archInFiles[i], "%b")) continue;
    if ((ext = strrchr(archInFiles[i], '.')) &&
        (!strcmp(ext, ".fcb") || !strcmp(ext, ".FCB"))) {
      readfcb(archInFiles[i], nav);
    }
  }
  /* read solution status files for ppp correction */
  for (i = 0; i < n; i++) {
    if (strstr(archInFiles[i], "%r") || strstr(archInFiles[i], "%b")) continue;
    if ((ext = strrchr(archInFiles[i], '.')) &&
        (!strcmp(ext, ".stat") || !strcmp(ext, ".STAT") ||
         !strcmp(ext, ".stec") || !strcmp(ext, ".STEC") ||
         !strcmp(ext, ".trp") || !strcmp(ext, ".TRP"))) {
      pppcorr_read(&nav->pppcorr, archInFiles[i]);
    }
  }
  /* read sbas message files */
  for (i = 0; i < n; i++) {
    if (strstr(archInFiles[i], "%r") || strstr(archInFiles[i], "%b")) continue;
    sbsreadmsg(archInFiles[i], prcopt->sbassatsel, sbs);
  }
  /* read lex message files */
  for (i = 0; i < n; i++) {
    if (strstr(archInFiles[i], "%r") || strstr(archInFiles[i], "%b")) continue;
    lexreadmsg(archInFiles[i], 0, lex);
  }
  /* allocate sbas ephemeris */
  nav->ns = nav->nsmax = NSATSBS * 2;
  if (!(nav->seph = (seph_t *)malloc(sizeof(seph_t) * nav->ns))) {
    printf("error : sbas ephem memory allocation");
    trace(1, "error : sbas ephem memory allocation");
    return;
  }
  for (i = 0; i < nav->ns; i++) nav->seph[i] = seph0;

  /* set rtcm file and initialize rtcm struct */
  rtcm_file[0] = rtcm_path[0] = '\0';
  fp_rtcm = NULL;

  for (i = 0; i < n; i++) {
    if ((ext = strrchr(archInFiles[i], '.')) &&
        (!strcmp(ext, ".rtcm3") || !strcmp(ext, ".RTCM3"))) {
      strcpy(rtcm_file, archInFiles[i]);
      init_rtcm(&rtcm);
      break;
    }
  }
}
/* free prec ephemeris and sbas data -----------------------------------------*/
static void freepreceph(nav_t *nav, sbs_t *sbs, lex_t *lex) {
  int i;

  trace(3, "freepreceph:\n");

  free(nav->peph);
  nav->peph = NULL;
  nav->ne = nav->nemax = 0;
  free(nav->pclk);
  nav->pclk = NULL;
  nav->nc = nav->ncmax = 0;
  free(nav->fcb);
  nav->fcb = NULL;
  nav->nf = nav->nfmax = 0;
  free(nav->seph);
  nav->seph = NULL;
  nav->ns = nav->nsmax = 0;
  free(sbs->msgs);
  sbs->msgs = NULL;
  sbs->n = sbs->nmax = 0;
  free(lex->msgs);
  lex->msgs = NULL;
  lex->n = lex->nmax = 0;
  for (i = 0; i < nav->nt; i++) {
    free(nav->tec[i].data);
    free(nav->tec[i].rms);
  }
  free(nav->tec);
  nav->tec = NULL;
  nav->nt = nav->ntmax = 0;

  if (fp_rtcm) fclose(fp_rtcm);
  free_rtcm(&rtcm);
}
/* read obs and nav data -----------------------------------------------------*/
static int readobsnav(gtime_t sGTimeStart, gtime_t sGTimeEnd,
                      double PrcsInterval, char **archInFiles,
                      const int *arIndex_InFiles, int n, const prcopt_t *prcopt,
                      obs_t *obs, nav_t *nav, sta_t *sta) {
  int i, j, ind = 0, nobs = 0, rcv = 1;

  trace(3, "readobsnav: sGTimeStart=%s n=%d\n", time_str(sGTimeStart, 0), n);

  obs->data = NULL;
  obs->n = obs->nmax = 0;
  nav->eph = NULL;
  nav->n = nav->nmax = 0;
  nav->geph = NULL;
  nav->ng = nav->ngmax = 0;
  nav->seph = NULL;
  nav->ns = nav->nsmax = 0;
  nepoch = 0;

  for (i = 0; i < n; i++) {
    if (checkbrk("")) return 0;
    if (strstr(archInFiles[i], ".sp3") || strstr(archInFiles[i], ".clk"))
      continue;

    if (arIndex_InFiles[i] != ind) {
      if (obs->n > nobs) rcv++;
      ind = arIndex_InFiles[i];
      nobs = obs->n;
    }
    /* read rinex obs and nav file */
    if (readrnxt(archInFiles[i], rcv, sGTimeStart, sGTimeEnd, PrcsInterval,
                 prcopt->rnxopt[rcv <= 1 ? 0 : 1], obs, nav,
                 rcv <= 2 ? sta + rcv - 1 : NULL) < 0) {
      checkbrk("error : insufficient memory");
      trace(1, "insufficient memory\n");
      return 0;
    }
  }
  if (obs->n <= 0) {
    checkbrk("error : no obs data");
    trace(1, "\n");
    return 0;
  }
  if (nav->n <= 0 && nav->ng <= 0 && nav->ns <= 0) {
    checkbrk("error : no nav data");
    trace(1, "\n");
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
/* free obs and nav data -----------------------------------------------------*/
static void freeobsnav(obs_t *obs, nav_t *nav) {
  trace(3, "freeobsnav:\n");

  free(obs->data);
  obs->data = NULL;
  obs->n = obs->nmax = 0;
  free(nav->eph);
  nav->eph = NULL;
  nav->n = nav->nmax = 0;
  free(nav->geph);
  nav->geph = NULL;
  nav->ng = nav->ngmax = 0;
  free(nav->seph);
  nav->seph = NULL;
  nav->ns = nav->nsmax = 0;
}
/* average of single position ------------------------------------------------*/
static int avepos(double *ra, int rcv, const obs_t *obs, const nav_t *nav,
                  const prcopt_t *opt) {
  obsd_t data[MAXOBS];
  gtime_t sGTimeStart = {0};
  sol_t sol = {{0}};
  int i, j, n = 0, m, iobs;
  char msg[128];

  trace(3, "avepos: rcv=%d obs.n=%d\n", rcv, obs->n);

  for (i = 0; i < 3; i++) ra[i] = 0.0;

  for (iobs = 0; (m = nextobsf(obs, &iobs, rcv)) > 0; iobs += m) {
    for (i = j = 0; i < m && i < MAXOBS; i++) {
      data[j] = obs->data[iobs + i];
      if ((satsys(data[j].sat, NULL) & opt->navsys) &&
          opt->exsats[data[j].sat - 1] != 1)
        j++;
    }
    if (j <= 0 || !screent(data[0].time, sGTimeStart, sGTimeStart, 1.0))
      continue; /* only 1 hz */

    if (!pntpos(data, j, nav, opt, &sol, NULL, NULL, msg)) continue;

    for (i = 0; i < 3; i++) ra[i] += sol.rr[i];
    n++;
  }
  if (n <= 0) {
    trace(1, "no average of base station position\n");
    return 0;
  }
  for (i = 0; i < 3; i++) ra[i] /= n;
  return 1;
}
/* station position from file ------------------------------------------------*/
static int getstapos(const char *file, char *name, double *r) {
  FILE *fp;
  char buff[256], sname[256], *p, *q;
  double pos[3];

  trace(3, "getstapos: file=%s name=%s\n", file, name);

  if (!(fp = fopen(file, "r"))) {
    trace(1, "station position file open error: %s\n", file);
    return 0;
  }
  while (fgets(buff, sizeof(buff), fp)) {
    if ((p = strchr(buff, '%'))) *p = '\0';

    if (sscanf(buff, "%lf %lf %lf %s", pos, pos + 1, pos + 2, sname) < 4)
      continue;

    for (p = sname, q = name; *p && *q; p++, q++) {
      if (toupper((int)*p) != toupper((int)*q)) break;
    }
    if (!*p) {
      pos[0] *= D2R;
      pos[1] *= D2R;
      pos2ecef(pos, r);
      fclose(fp);
      return 1;
    }
  }
  fclose(fp);
  trace(1, "no station position: %s %s\n", name, file);
  return 0;
}
/* antenna phase center position ---------------------------------------------*/
static int antpos(prcopt_t *opt, int rcvno, const obs_t *obs, const nav_t *nav,
                  const sta_t *sta, const char *posfile) {
  double *rr = rcvno == 1 ? opt->ru : opt->rb, del[3], pos[3], dr[3] = {0};
  int i, postype = rcvno == 1 ? opt->rovpos : opt->refpos;
  char *name;

  trace(3, "antpos  : rcvno=%d\n", rcvno);

  if (postype == POSOPT_SINGLE) { /* average of single position */
    if (!avepos(rr, rcvno, obs, nav, opt)) {
      printf("error : station pos computation");
      return 0;
    }
  } else if (postype == POSOPT_FILE) { /* read from position file */
    name = stas[rcvno == 1 ? 0 : 1].name;
    if (!getstapos(posfile, name, rr)) {
      printf("error : no position of %s in %s", name, posfile);
      return 0;
    }
  } else if (postype == POSOPT_RINEX) { /* get from rinex header */
    if (norm(stas[rcvno == 1 ? 0 : 1].pos, 3) <= 0.0) {
      printf("error : no position in rinex header");
      trace(1, "no position position in rinex header\n");
      return 0;
    }
    /* antenna delta */
    if (stas[rcvno == 1 ? 0 : 1].deltype == 0) { /* enu */
      for (i = 0; i < 3; i++) del[i] = stas[rcvno == 1 ? 0 : 1].del[i];
      del[2] += stas[rcvno == 1 ? 0 : 1].hgt;
      ecef2pos(stas[rcvno == 1 ? 0 : 1].pos, pos);
      enu2ecef(pos, del, dr);
    } else { /* xyz */
      for (i = 0; i < 3; i++) dr[i] = stas[rcvno == 1 ? 0 : 1].del[i];
    }
    for (i = 0; i < 3; i++) rr[i] = stas[rcvno == 1 ? 0 : 1].pos[i] + dr[i];
  }
  return 1;
}
/* open procssing session ----------------------------------------------------*/
static int openses(const prcopt_t *sProcessOption,
                   const solopt_t *sSolutionOption, const filopt_t *sFileOption,
                   nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr) {
  int i;

  trace(3, "openses :\n");

  /* read satellite antenna parameters */
  if (*sFileOption->satantp && !(readpcv(sFileOption->satantp, pcvs))) {
    printf("error : no sat ant pcv in %s", sFileOption->satantp);
    trace(1, "sat antenna pcv read error: %s\n", sFileOption->satantp);
    return 0;
  }
  /* read receiver antenna parameters */
  if (*sFileOption->rcvantp && !(readpcv(sFileOption->rcvantp, pcvr))) {
    printf("error : no rec ant pcv in %s", sFileOption->rcvantp);
    trace(1, "rec antenna pcv read error: %s\n", sFileOption->rcvantp);
    return 0;
  }
  /* open geoid data */
  if (sSolutionOption->geoid > 0 && *sFileOption->geoid) {
    if (!opengeoid(sSolutionOption->geoid, sFileOption->geoid)) {
      printf("error : no geoid data %s", sFileOption->geoid);
      trace(2, "no geoid data %s\n", sFileOption->geoid);
    }
  }
  /* use satellite L2 offset if L5 offset does not exists */
  for (i = 0; i < pcvs->n; i++) {
    if (norm(pcvs->pcv[i].off[2], 3) > 0.0) continue;
    matcpy(pcvs->pcv[i].off[2], pcvs->pcv[i].off[1], 3, 1);
    matcpy(pcvs->pcv[i].var[2], pcvs->pcv[i].var[1], 19, 1);
  }
  for (i = 0; i < pcvr->n; i++) {
    if (norm(pcvr->pcv[i].off[2], 3) > 0.0) continue;
    matcpy(pcvr->pcv[i].off[2], pcvr->pcv[i].off[1], 3, 1);
    matcpy(pcvr->pcv[i].var[2], pcvr->pcv[i].var[1], 19, 1);
  }
  return 1;
}
/* close procssing session ---------------------------------------------------*/
static void closeses(nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr) {
  trace(3, "closeses:\n");

  /* free antenna parameters */
  free(pcvs->pcv);
  pcvs->pcv = NULL;
  pcvs->n = pcvs->nmax = 0;
  free(pcvr->pcv);
  pcvr->pcv = NULL;
  pcvr->n = pcvr->nmax = 0;

  /* close geoid data */
  closegeoid();

  /* free erp data */
  free(nav->erp.data);
  nav->erp.data = NULL;
  nav->erp.n = nav->erp.nmax = 0;

  /* close solution statistics and debug trace */
  rtkclosestat();
  traceclose();
}
/* set antenna parameters ----------------------------------------------------*/
static void setpcv(gtime_t time, prcopt_t *sProcessOption, nav_t *nav,
                   const pcvs_t *pcvs, const pcvs_t *pcvr, const sta_t *sta) {
  pcv_t *pcv;
  double pos[3], del[3];
  int i, j,
      mode = PMODE_DGPS <= sProcessOption->mode &&
             sProcessOption->mode <= PMODE_FIXED;
  char id[64];

  /* set satellite antenna parameters */
  for (i = 0; i < MAXSAT; i++) {
    if (!(satsys(i + 1, NULL) & sProcessOption->navsys)) continue;
    if (!(pcv = searchpcv(i + 1, "", time, pcvs))) {
      satno2id(i + 1, id);
      trace(3, "no satellite antenna pcv: %s\n", id);
      continue;
    }
    nav->pcvs[i] = *pcv;
  }
  for (i = 0; i < (mode ? 2 : 1); i++) {
    if (!strcmp(sProcessOption->anttype[i],
                "*")) { /* set by station parameters */
      strcpy(sProcessOption->anttype[i], sta[i].antdes);
      if (sta[i].deltype == 1) { /* xyz */
        if (norm(sta[i].pos, 3) > 0.0) {
          ecef2pos(sta[i].pos, pos);
          ecef2enu(pos, sta[i].del, del);
          for (j = 0; j < 3; j++) sProcessOption->antdel[i][j] = del[j];
        }
      } else { /* enu */
        for (j = 0; j < 3; j++) sProcessOption->antdel[i][j] = stas[i].del[j];
      }
    }
    if (!(pcv = searchpcv(0, sProcessOption->anttype[i], time, pcvr))) {
      trace(2, "no receiver antenna pcv: %s\n", sProcessOption->anttype[i]);
      *sProcessOption->anttype[i] = '\0';
      continue;
    }
    strcpy(sProcessOption->anttype[i], pcv->type);
    sProcessOption->pcvr[i] = *pcv;
  }
}
/* read ocean tide loading parameters ----------------------------------------*/
static void readotl(prcopt_t *sProcessOption, const char *file,
                    const sta_t *sta) {
  int i, mode = PMODE_DGPS <= sProcessOption->mode &&
                sProcessOption->mode <= PMODE_FIXED;

  for (i = 0; i < (mode ? 2 : 1); i++) {
    readblq(file, sta[i].name, sProcessOption->odisp[i]);
  }
}
/* write header to output file -----------------------------------------------*/
static int outhead(const char *archOutFiles, char **archInFiles, int NumInFiles,
                   const prcopt_t *sProcessOption,
                   const solopt_t *sSolutionOption) {
  FILE *fp = stdout;

  trace(3, "outhead: archOutFiles=%s NumInFiles=%d\n", archOutFiles,
        NumInFiles);

  if (*archOutFiles) {
    createdir(archOutFiles);

    if (!(fp = fopen(archOutFiles, "w"))) {
      printf("error : open output file %s", archOutFiles);
      return 0;
    }
  }
  /* output header */
  outheader(fp, archInFiles, NumInFiles, sProcessOption, sSolutionOption);

  if (*archOutFiles) fclose(fp);

  return 1;
}
/* open output file for append -----------------------------------------------*/
static FILE *openfile(const char *archOutFiles) {
  trace(3, "openfile: archOutFiles=%s\n", archOutFiles);

  return !*archOutFiles ? stdout : fopen(archOutFiles, "a");
}
/* execute processing session ------------------------------------------------*/
static int execses(gtime_t sGTimeStart, gtime_t sGTimeEnd, double PrcsInterval,
                   const prcopt_t *sProcessOption,
                   const solopt_t *sSolutionOption, const filopt_t *sFileOption,
                   int OpenFileFlag, char **archInFiles,
                   const int *arIndex_InFiles, int NumInFiles,
                   char *archOutFiles) {
  FILE *fp;
  prcopt_t popt_ = *sProcessOption;
  char posfile[100], ambfile[100], resfile[100], azlfile[100], tracefile[100],
      statfile[100], outposmat[100], path[100], *ext;

  trace(3, "execses : NumInFiles=%d archOutFiles=%s\n", NumInFiles,
        archOutFiles);

  /* open debug trace */
  if (OpenFileFlag && sSolutionOption->trace > 0) {
    if (*archOutFiles) {
      strcpy(tracefile, archOutFiles);
      // strcat(tracefile, ".trace");
    } else {
      strcpy(tracefile, sFileOption->tracepath);
      strcat(tracefile, sFileOption->file_id);
      strcat(tracefile, sFileOption->prc_id4opname);
      strcat(tracefile, ".trace");
    }
    traceclose();
    traceopen(tracefile);
    tracelevel(sSolutionOption->trace);
  }
  /* read ionosphere data file */
  if (*sFileOption->iono && (ext = strrchr(sFileOption->iono, '.'))) {
    if (strlen(ext) == 4 && (ext[3] == 'i' || ext[3] == 'I')) {
      reppath(sFileOption->iono, path, sGTimeStart, "", "");
      readtec(path, &navs, 1);
    }
  }
  /* read erp data */
  if (*sFileOption->eop) {
    free(navs.erp.data);
    navs.erp.data = NULL;
    navs.erp.n = navs.erp.nmax = 0;
    reppath(sFileOption->eop, path, sGTimeStart, "", "");
    if (!readerp(path, &navs.erp)) {
      printf("error : no erp data %s", path);
      trace(2, "no erp data %s\n", path);
    }
  }
  /* read obs and nav data */
  if (!readobsnav(sGTimeStart, sGTimeEnd, PrcsInterval, archInFiles,
                  arIndex_InFiles, NumInFiles, &popt_, &obss, &navs, stas))
    return 0;

  /* read dcb parameters */
  if (*sFileOption->dcb) {
    reppath(sFileOption->dcb, path, sGTimeStart, "", "");
    readdcb(path, &navs, stas);
  }
  /* set antenna paramters */
  if (popt_.mode != PMODE_SINGLE) {
    setpcv(obss.n > 0 ? obss.data[0].time : timeget(), &popt_, &navs,
           &RecAntParas, &SatAntParas, stas);
  }
  /* read ocean tide loading parameters */
  if (popt_.mode > PMODE_SINGLE && *sFileOption->blq) {
    readotl(&popt_, sFileOption->blq, stas);
  }
  /* rover/reference fixed position */
  if (popt_.mode == PMODE_FIXED) {
    if (!antpos(&popt_, 1, &obss, &navs, stas, sFileOption->stapos)) {
      freeobsnav(&obss, &navs);
      return 0;
    }
  } else if (PMODE_DGPS <= popt_.mode && popt_.mode <= PMODE_STATIC) {
    if (!antpos(&popt_, 2, &obss, &navs, stas, sFileOption->stapos)) {
      freeobsnav(&obss, &navs);
      return 0;
    }
  }

  return 1;
  
  /* open solution status */
  if (OpenFileFlag && sSolutionOption->sstat > 0) {
    strcpy(statfile, sFileOption->solstatpath);
    strcat(statfile, sFileOption->file_id);
    strcat(statfile, sFileOption->prc_id4opname);
    strcat(statfile, ".pos.stat");
    rtkclosestat();
    rtkopenstat(statfile, sSolutionOption->sstat);
  }
  /* open solution ambiguity */
  if (OpenFileFlag && sFileOption->ambpath[0] != '\0') {
    strcpy(ambfile, sFileOption->ambpath);
    strcat(ambfile, sFileOption->file_id);
    strcat(ambfile, sFileOption->prc_id4opname);
    strcat(ambfile, ".amb");
    rtkcloseamb();
    rtkopenamb(ambfile, sSolutionOption->samb);
  }
  /* open solution residual */
  if (OpenFileFlag && sFileOption->respath[0] != '\0') {
    strcpy(resfile, sFileOption->respath);
    strcat(resfile, sFileOption->file_id);
    strcat(resfile, sFileOption->prc_id4opname);
    strcat(resfile, ".resd");
    rtkcloseres();
    rtkopenres(resfile, sSolutionOption->sres);
  }
  /* open solution azimuth/elevation angle */
  if (OpenFileFlag && sFileOption->azlpath[0] != '\0') {
    strcpy(azlfile, sFileOption->azlpath);
    strcat(azlfile, sFileOption->file_id);
    strcat(azlfile, sFileOption->prc_id4opname);
    strcat(azlfile, ".azl");
    rtkcloseazl();
    rtkopenazl(azlfile, sSolutionOption->sazl);
  }

  iobsu = iobsr = isbs = ilex = revs = aborts = 0;

  strcpy(posfile, sFileOption->pospath);
  strcat(posfile, sFileOption->file_id);
  strcat(posfile, sFileOption->prc_id4opname);
  strcpy(outposmat, posfile);

  archOutFiles = strcat(posfile, ".pos");
  /* write header to output file */
  if (OpenFileFlag && !outhead(archOutFiles, archInFiles, NumInFiles, &popt_,
                               sSolutionOption)) {
    freeobsnav(&obss, &navs);
    return 0;
  }

  if (popt_.mode == PMODE_SINGLE || popt_.soltype == 0) {
    if ((fp = openfile(archOutFiles))) {
      procpos(fp, &popt_, sSolutionOption, 0); /* forward */
      fclose(fp);
    }
  } else if (popt_.soltype == 1) {
    if ((fp = openfile(archOutFiles))) {
      revs = 1;
      iobsu = iobsr = obss.n - 1;
      isbs = sbss.n - 1;
      ilex = lexs.n - 1;
      procpos(fp, &popt_, sSolutionOption, 0); /* backward */
      fclose(fp);
    }
  } else { /* combined */
    solf = (sol_t *)malloc(sizeof(sol_t) * nepoch);
    solb = (sol_t *)malloc(sizeof(sol_t) * nepoch);
    rbf = (double *)malloc(sizeof(double) * nepoch * 3);
    rbb = (double *)malloc(sizeof(double) * nepoch * 3);

    if (solf && solb) {
      isolf = isolb = 0;
      procpos(NULL, &popt_, sSolutionOption, 1); /* forward */
      revs = 1;
      iobsu = iobsr = obss.n - 1;
      isbs = sbss.n - 1;
      ilex = lexs.n - 1;
      procpos(NULL, &popt_, sSolutionOption, 1); /* backward */

      /* combine forward/backward solutions */
      if (!aborts && (fp = openfile(archOutFiles))) {
        combres(fp, &popt_, sSolutionOption);
        fclose(fp);
      }
    } else
      printf("error : memory allocation");
    free(solf);
    free(solb);
    free(rbf);
    free(rbb);
  }
  /* free obs and nav data */
  freeobsnav(&obss, &navs);

  return aborts ? 1 : 0;
}

/* execute processing session for each rover ---------------------------------*/
static int execses_r(gtime_t sGTimeStart, gtime_t sGTimeEnd,
                     double PrcsInterval, const prcopt_t *sProcessOption,
                     const solopt_t *sSolutionOption,
                     const filopt_t *sFileOption, int OpenFileFlag,
                     char **archInFiles, const int *arIndex_InFiles,
                     int NumInFiles, char *archOutFiles, const char *pchRover) {
  gtime_t t0 = {0};
  int i, stat = 0;
  char *ifile[MAXINFILE], archReplacedOutFile[1024], *rov_, *p, *q, s[64] = "";

  trace(3, "execses_r: NumInFiles=%d archOutFiles=%s\n", NumInFiles,
        archOutFiles);

  for (i = 0; i < NumInFiles; i++)
    if (strstr(archInFiles[i], "%r")) break;

  if (i < NumInFiles) { /* include rover keywords */
    if (!(rov_ = (char *)malloc(strlen(pchRover) + 1))) return 0;
    strcpy(rov_, pchRover);

    for (i = 0; i < NumInFiles; i++) {
      if (!(ifile[i] = (char *)malloc(1024))) {
        free(rov_);
        for (; i >= 0; i--) free(ifile[i]);
        return 0;
      }
    }
    for (p = rov_;; p = q + 1) { /* for each rover */
      if ((q = strchr(p, ' '))) *q = '\0';

      if (*p) {
        strcpy(proc_rov, p);
        if (sGTimeStart.time)
          time2str(sGTimeStart, s, 0);
        else
          *s = '\0';
        if (checkbrk("reading    : %s", s)) {
          stat = 1;
          break;
        }
        for (i = 0; i < NumInFiles; i++)
          reppath(archInFiles[i], ifile[i], t0, p, "");
        reppath(archOutFiles, archReplacedOutFile, t0, p, "");

        /* execute processing session */
        stat = execses(sGTimeStart, sGTimeEnd, PrcsInterval, sProcessOption,
                       sSolutionOption, sFileOption, OpenFileFlag, ifile,
                       arIndex_InFiles, NumInFiles, archReplacedOutFile);
      }
      if (stat == 1 || !q) break;
    }
    free(rov_);
    for (i = 0; i < NumInFiles; i++) free(ifile[i]);
  } else {
    /* execute processing session */
    stat = execses(sGTimeStart, sGTimeEnd, PrcsInterval, sProcessOption,
                   sSolutionOption, sFileOption, OpenFileFlag, archInFiles,
                   arIndex_InFiles, NumInFiles, archOutFiles);
  }
  return stat;
}
/* execute processing session for each base station --------------------------*/
static int execses_b(gtime_t sGTimeStart, gtime_t sGTimeEnd,
                     double PrcsInterval, const prcopt_t *sProcessOption,
                     const solopt_t *sSolutionOption,
                     const filopt_t *sFileOption, int OpenFileFlag,
                     char **archInFiles, const int *arIndex_InFiles,
                     int NumInFiles, char *archOutFiles, const char *pchRover,
                     const char *pchBase) {
  gtime_t t0 = {0};
  int i, stat = 0;
  char *ifile[MAXINFILE], archReplacedOutFile[1024], *base_, *p, *q, s[64];

  trace(3, "execses_b: NumInFiles=%d archOutFiles=%s\n", NumInFiles,
        archOutFiles);

  /* read prec ephemeris and sbas data */
  readpreceph(archInFiles, NumInFiles, sProcessOption, &navs, &sbss, &lexs);

  for (i = 0; i < NumInFiles; i++)
    if (strstr(archInFiles[i], "%b")) break;

  if (i < NumInFiles) { /* include base station keywords */
    if (!(base_ = (char *)malloc(strlen(pchBase) + 1))) {
      freepreceph(&navs, &sbss, &lexs);
      return 0;
    }
    strcpy(base_, pchBase);

    for (i = 0; i < NumInFiles; i++) {
      if (!(ifile[i] = (char *)malloc(1024))) {
        free(base_);
        for (; i >= 0; i--) free(ifile[i]);
        freepreceph(&navs, &sbss, &lexs);
        return 0;
      }
    }
    for (p = base_;; p = q + 1) { /* for each base station */
      if ((q = strchr(p, ' '))) *q = '\0';

      if (*p) {
        strcpy(proc_base, p);
        if (sGTimeStart.time)
          time2str(sGTimeStart, s, 0);
        else
          *s = '\0';
        if (checkbrk("reading    : %s", s)) {
          stat = 1;
          break;
        }
        for (i = 0; i < NumInFiles; i++)
          reppath(archInFiles[i], ifile[i], t0, "", p);
        reppath(archOutFiles, archReplacedOutFile, t0, "", p);

        stat = execses_r(sGTimeStart, sGTimeEnd, PrcsInterval, sProcessOption,
                         sSolutionOption, sFileOption, OpenFileFlag, ifile,
                         arIndex_InFiles, NumInFiles, archReplacedOutFile,
                         pchRover);
      }
      if (stat == 1 || !q) break;
    }
    free(base_);
    for (i = 0; i < NumInFiles; i++) free(ifile[i]);
  } else {
    stat = execses_r(sGTimeStart, sGTimeEnd, PrcsInterval, sProcessOption,
                     sSolutionOption, sFileOption, OpenFileFlag, archInFiles,
                     arIndex_InFiles, NumInFiles, archOutFiles, pchRover);
  }
  /* free prec ephemeris and sbas data */
  freepreceph(&navs, &sbss, &lexs);

  return stat;
}
/* post-processing positioning -------------------------------------------------
 * post-processing positioning
 * args   : gtime_t sGTimeStart        I   processing start time
 *(sGTimeStart.time==0: no limit) : gtime_t sGTimeEnd          I   processing
 *end time   (sGTimeEnd.time==0: no limit) double PrcsInterval        I
 *processing interval  (s) (0:all) double PrcsUnitTime        I   processing
 *unit time (s) (0:all) prcopt_t *sProcessOption   I   processing options
 *          solopt_t *sSolutionOption   I   solution options
 *          filopt_t *sFileOption   I   file options
 *          char   **archInFiles  I   input files (see below)
 *          int    n         I   number of input files
 *          char   *archOutFiles  I   output file ("":stdout, see below)
 *          char   *pchRover      I   rover id list        (separated by " ")
 *          char   *pchBase     I   base station id list (separated by " ")
 * return : status (0:ok,0>:error,1:aborted)
 * notes  : input files should contain observation data, navigation data,
 *precise ephemeris/clock (optional), sbas log file (optional), ssr message log
 *file (optional) and tec grid file (optional). only the first observation data
 *file in the input files is recognized as the rover data.
 *
 *          the type of an input file is recognized by the file extention as ]
 *          follows:
 *              .sp3,.SP3,.eph*,.EPH*: precise ephemeris (sp3c)
 *              .sbs,.SBS,.ems,.EMS  : sbas message log files (rtklib or ems)
 *              .lex,.LEX            : qzss lex message log files
 *              .rtcm3,.RTCM3        : ssr message log files (rtcm3)
 *              .*i,.*I              : tec grid files (ionex)
 *              .fcb,.FCB            : satellite fcb
 *              others               : rinex obs, nav, gnav, hnav, qnav or clock
 *
 *          inputs files can include wild-cards (*). if an file includes
 *          wild-cards, the wild-card expanded multiple files are used.
 *
 *          inputs files can include keywords. if an file includes keywords,
 *          the keywords are replaced by date, time, rover id and base station
 *          id and multiple session analyses run. refer reppath() for the
 *          keywords.
 *
 *          the output file can also include keywords. if the output file does
 *          not include keywords. the results of all multiple session analyses
 *          are output to a single output file.
 *
 *          ssr corrections are valid only for forward estimation.
 *-----------------------------------------------------------------------------*/
extern int postpos(gtime_t sGTimeStart, gtime_t sGTimeEnd, double PrcsInterval,
                   double PrcsUnitTime, const prcopt_t *sProcessOption,
                   const solopt_t *sSolutionOption, const filopt_t *sFileOption,
                   char **archInFiles, int NumInFiles, char *archOutFiles,
                   const char *pchRover, const char *pchBase) {
  gtime_t tts, tte, ttte;
  double tunit, tss;
  int i, j, k, nf, stat = 0, week, OpenFileFlag = 1,
                   arIndex_InFiles[MAXINFILE] = {0};
  char *ifile[MAXINFILE], archReplacedOutFile[1024], *ext;

  trace(3,
        "postpos : PrcsInterval=%.0f PrcsUnitTime=%.0f NumInFiles=%d "
        "archOutFiles=%s\n",
        PrcsInterval, PrcsUnitTime, NumInFiles, archOutFiles);

  /* open processing session */
  if (!openses(sProcessOption, sSolutionOption, sFileOption, &navs,
               &RecAntParas, &SatAntParas))
    return -1;

  if (sGTimeStart.time != 0 && sGTimeEnd.time != 0 && PrcsUnitTime >= 0.0) {
    if (timediff(sGTimeEnd, sGTimeStart) < 0.0) {
      printf("error : no period");
      closeses(&navs, &RecAntParas, &SatAntParas);
      return 0;
    }
    for (i = 0; i < MAXINFILE; i++) {
      if (!(ifile[i] = (char *)malloc(1024))) {
        for (; i >= 0; i--) free(ifile[i]);
        closeses(&navs, &RecAntParas, &SatAntParas);
        return -1;
      }
    }
    if (PrcsUnitTime == 0.0 || PrcsUnitTime > 86400.0 * MAXPRCDAYS)
      PrcsUnitTime = 86400.0 * MAXPRCDAYS;
    tunit = PrcsUnitTime < 86400.0 ? PrcsUnitTime : 86400.0;
    tss = tunit * (int)floor(time2gpst(sGTimeStart, &week) / tunit);

    for (i = 0;; i++) { /* for each periods */
      tts = gpst2time(week, tss + i * PrcsUnitTime);
      tte = timeadd(tts, PrcsUnitTime - DTTOL);
      if (timediff(tts, sGTimeEnd) > 0.0) break;
      if (timediff(tts, sGTimeStart) < 0.0) tts = sGTimeStart;
      if (timediff(tte, sGTimeEnd) > 0.0) tte = sGTimeEnd;

      strcpy(proc_rov, "");
      strcpy(proc_base, "");
      if (checkbrk("reading    : %s", time_str(tts, 0))) {
        stat = 1;
        break;
      }
      for (j = k = nf = 0; j < NumInFiles; j++) {
        ext = strrchr(archInFiles[j], '.');

        if (ext && (!strcmp(ext, ".rtcm3") || !strcmp(ext, ".RTCM3"))) {
          strcpy(ifile[nf++], archInFiles[j]);
        } else {
          /* include next day precise ephemeris or rinex brdc nav */
          ttte = tte;
          if (ext && (!strcmp(ext, ".sp3") || !strcmp(ext, ".SP3") ||
                      !strcmp(ext, ".eph") || !strcmp(ext, ".EPH"))) {
            ttte = timeadd(ttte, 3600.0);
          } else if (strstr(archInFiles[j], "brdc")) {
            ttte = timeadd(ttte, 7200.0);
          }
          nf += reppaths(archInFiles[j], ifile + nf, MAXINFILE - nf, tts, ttte,
                         "", "");
        }
        while (k < nf) arIndex_InFiles[k++] = j;

        if (nf >= MAXINFILE) {
          trace(2, "too many input files. trancated\n");
          break;
        }
      }
      if (!reppath(archOutFiles, archReplacedOutFile, tts, "", "") && i > 0)
        OpenFileFlag = 0;

      /* execute processing session */
      stat = execses_b(tts, tte, PrcsInterval, sProcessOption, sSolutionOption,
                       sFileOption, OpenFileFlag, ifile, arIndex_InFiles, nf,
                       archReplacedOutFile, pchRover, pchBase);

      if (stat == 1) break;
    }
    for (i = 0; i < MAXINFILE; i++) free(ifile[i]);
  } else if (sGTimeStart.time != 0) {
    for (i = 0; i < NumInFiles && i < MAXINFILE; i++) {
      if (!(ifile[i] = (char *)malloc(1024))) {
        for (; i >= 0; i--) free(ifile[i]);
        return -1;
      }
      reppath(archInFiles[i], ifile[i], sGTimeStart, "", "");
      arIndex_InFiles[i] = i;
    }
    reppath(archOutFiles, archReplacedOutFile, sGTimeStart, "", "");

    /* execute processing session */
    stat = execses_b(sGTimeStart, sGTimeEnd, PrcsInterval, sProcessOption,
                     sSolutionOption, sFileOption, 1, ifile, arIndex_InFiles,
                     NumInFiles, archReplacedOutFile, pchRover, pchBase);

    for (i = 0; i < NumInFiles && i < MAXINFILE; i++) free(ifile[i]);
  } else {
    for (i = 0; i < NumInFiles; i++) arIndex_InFiles[i] = i;

    /* execute processing session */
    stat =
        execses_b(sGTimeStart, sGTimeEnd, PrcsInterval, sProcessOption,
                  sSolutionOption, sFileOption, 1, archInFiles, arIndex_InFiles,
                  NumInFiles, archOutFiles, pchRover, pchBase);
  }
  /* close processing session */
  closeses(&navs, &RecAntParas, &SatAntParas);

  return stat;
}
