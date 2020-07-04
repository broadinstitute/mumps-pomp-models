/* pomp C snippet file: snippets */
/* Time: 2020-06-30 21:52:26.476 +0900 */
/* Salt: B543EB1AD7720C61EEA31B29 */

#include <pomp.h>
#include <R_ext/Rdynload.h>

 


/* C snippet: 'rinit' */
#define beta		(__p[__parindex[0]])
#define rho		(__p[__parindex[1]])
#define psi		(__p[__parindex[2]])
#define w		(__p[__parindex[3]])
#define pop		(__p[__parindex[4]])
#define S_0		(__p[__parindex[5]])
#define E_0		(__p[__parindex[6]])
#define I_0		(__p[__parindex[7]])
#define R_0		(__p[__parindex[8]])
#define sigma		(__p[__parindex[9]])
#define gamma		(__p[__parindex[10]])
#define intervention		(__p[__parindex[11]])
#define S		(__x[__stateindex[0]])
#define E		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define R		(__x[__stateindex[3]])
#define C		(__x[__stateindex[4]])

void __pomp_rinit (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)
{
 
  double m = pop/(S_0 + E_0 + I_0 + R_0);

  S = nearbyint(m*S_0);
  E = nearbyint(m*E_0);
  I = nearbyint(m*I_0);
  R = nearbyint(m*R_0);

  C = 0;
 
}

#undef beta
#undef rho
#undef psi
#undef w
#undef pop
#undef S_0
#undef E_0
#undef I_0
#undef R_0
#undef sigma
#undef gamma
#undef intervention
#undef S
#undef E
#undef I
#undef R
#undef C

/* C snippet: 'step.fn' */
#define beta		(__p[__parindex[0]])
#define rho		(__p[__parindex[1]])
#define psi		(__p[__parindex[2]])
#define w		(__p[__parindex[3]])
#define pop		(__p[__parindex[4]])
#define S_0		(__p[__parindex[5]])
#define E_0		(__p[__parindex[6]])
#define I_0		(__p[__parindex[7]])
#define R_0		(__p[__parindex[8]])
#define sigma		(__p[__parindex[9]])
#define gamma		(__p[__parindex[10]])
#define intervention		(__p[__parindex[11]])
#define S		(__x[__stateindex[0]])
#define E		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define R		(__x[__stateindex[3]])
#define C		(__x[__stateindex[4]])

void __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t, double dt)
{
 
  double dSE, dEI, dIR;
  double betaf;

  if (intervention <= t) {
    betaf = w * beta;
  } else {
    betaf = beta;
  }
  
  double rate[3], trans[3];

  rate[0] = betaf * I/pop;
  rate[1] = sigma;
  rate[2] = gamma;

  reulermultinom(1, S, &rate[0], dt, &trans[0]);
  reulermultinom(1, E, &rate[1], dt, &trans[1]);
  reulermultinom(1, I, &rate[2], dt, &trans[2]);

  dSE = trans[0];
  dEI = trans[1];
  dIR = trans[2];

  S += -dSE;
  E += dSE - dEI;
  I += dEI - dIR;
  R += dIR;

  C += dIR;
 
}

#undef beta
#undef rho
#undef psi
#undef w
#undef pop
#undef S_0
#undef E_0
#undef I_0
#undef R_0
#undef sigma
#undef gamma
#undef intervention
#undef S
#undef E
#undef I
#undef R
#undef C

/* C snippet: 'rmeasure' */
#define beta		(__p[__parindex[0]])
#define rho		(__p[__parindex[1]])
#define psi		(__p[__parindex[2]])
#define w		(__p[__parindex[3]])
#define pop		(__p[__parindex[4]])
#define S_0		(__p[__parindex[5]])
#define E_0		(__p[__parindex[6]])
#define I_0		(__p[__parindex[7]])
#define R_0		(__p[__parindex[8]])
#define sigma		(__p[__parindex[9]])
#define gamma		(__p[__parindex[10]])
#define intervention		(__p[__parindex[11]])
#define S		(__x[__stateindex[0]])
#define E		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define R		(__x[__stateindex[3]])
#define C		(__x[__stateindex[4]])
#define cases		(__y[__obsindex[0]])

void __pomp_rmeasure (double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 
  double m = rho * C;
  double v = m * (1.0 - rho + psi * psi * m);
  double tol = 1.0e-18;
  cases = rnorm(m, sqrt(v) + tol);
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
 
}

#undef beta
#undef rho
#undef psi
#undef w
#undef pop
#undef S_0
#undef E_0
#undef I_0
#undef R_0
#undef sigma
#undef gamma
#undef intervention
#undef S
#undef E
#undef I
#undef R
#undef C
#undef cases

/* C snippet: 'dmeasure' */
#define beta		(__p[__parindex[0]])
#define rho		(__p[__parindex[1]])
#define psi		(__p[__parindex[2]])
#define w		(__p[__parindex[3]])
#define pop		(__p[__parindex[4]])
#define S_0		(__p[__parindex[5]])
#define E_0		(__p[__parindex[6]])
#define I_0		(__p[__parindex[7]])
#define R_0		(__p[__parindex[8]])
#define sigma		(__p[__parindex[9]])
#define gamma		(__p[__parindex[10]])
#define intervention		(__p[__parindex[11]])
#define S		(__x[__stateindex[0]])
#define E		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define R		(__x[__stateindex[3]])
#define C		(__x[__stateindex[4]])
#define cases		(__y[__obsindex[0]])
#define lik		(__lik[0])

void __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 
  double m = rho * C;
  double v = m * (1.0 - rho + psi * psi * m);
  double tol = 1.0e-18;
  if (cases > 0.0) {
    lik = pnorm(cases + 0.5, m, sqrt(v) + tol, 1, 0) - pnorm(cases - 0.5, m, sqrt(v)  + tol, 1, 0) + tol;
  } else {
    lik = pnorm(cases + 0.5, m, sqrt(v) + tol, 1, 0) + tol;
  }
  if (give_log) lik = log(lik);
 
}

#undef beta
#undef rho
#undef psi
#undef w
#undef pop
#undef S_0
#undef E_0
#undef I_0
#undef R_0
#undef sigma
#undef gamma
#undef intervention
#undef S
#undef E
#undef I
#undef R
#undef C
#undef cases
#undef lik

/* C snippet: 'toEst' */
#define beta		(__p[__parindex[0]])
#define rho		(__p[__parindex[1]])
#define psi		(__p[__parindex[2]])
#define w		(__p[__parindex[3]])
#define pop		(__p[__parindex[4]])
#define S_0		(__p[__parindex[5]])
#define E_0		(__p[__parindex[6]])
#define I_0		(__p[__parindex[7]])
#define R_0		(__p[__parindex[8]])
#define sigma		(__p[__parindex[9]])
#define gamma		(__p[__parindex[10]])
#define intervention		(__p[__parindex[11]])
#define T_beta		(__pt[__parindex[0]])
#define T_rho		(__pt[__parindex[1]])
#define T_psi		(__pt[__parindex[2]])
#define T_w		(__pt[__parindex[3]])
#define T_pop		(__pt[__parindex[4]])
#define T_S_0		(__pt[__parindex[5]])
#define T_E_0		(__pt[__parindex[6]])
#define T_I_0		(__pt[__parindex[7]])
#define T_R_0		(__pt[__parindex[8]])
#define T_sigma		(__pt[__parindex[9]])
#define T_gamma		(__pt[__parindex[10]])
#define T_intervention		(__pt[__parindex[11]])

void __pomp_to_trans (double *__pt, const double *__p, const int *__parindex)
{
 	T_beta = log(beta);
	T_psi = log(psi);
	T_rho = logit(rho);
	T_w = logit(w); 
}

#undef beta
#undef rho
#undef psi
#undef w
#undef pop
#undef S_0
#undef E_0
#undef I_0
#undef R_0
#undef sigma
#undef gamma
#undef intervention
#undef T_beta
#undef T_rho
#undef T_psi
#undef T_w
#undef T_pop
#undef T_S_0
#undef T_E_0
#undef T_I_0
#undef T_R_0
#undef T_sigma
#undef T_gamma
#undef T_intervention

/* C snippet: 'fromEst' */
#define beta		(__p[__parindex[0]])
#define rho		(__p[__parindex[1]])
#define psi		(__p[__parindex[2]])
#define w		(__p[__parindex[3]])
#define pop		(__p[__parindex[4]])
#define S_0		(__p[__parindex[5]])
#define E_0		(__p[__parindex[6]])
#define I_0		(__p[__parindex[7]])
#define R_0		(__p[__parindex[8]])
#define sigma		(__p[__parindex[9]])
#define gamma		(__p[__parindex[10]])
#define intervention		(__p[__parindex[11]])
#define T_beta		(__pt[__parindex[0]])
#define T_rho		(__pt[__parindex[1]])
#define T_psi		(__pt[__parindex[2]])
#define T_w		(__pt[__parindex[3]])
#define T_pop		(__pt[__parindex[4]])
#define T_S_0		(__pt[__parindex[5]])
#define T_E_0		(__pt[__parindex[6]])
#define T_I_0		(__pt[__parindex[7]])
#define T_R_0		(__pt[__parindex[8]])
#define T_sigma		(__pt[__parindex[9]])
#define T_gamma		(__pt[__parindex[10]])
#define T_intervention		(__pt[__parindex[11]])

void __pomp_from_trans (double *__p, const double *__pt, const int *__parindex)
{
 	beta = exp(T_beta);
	psi = exp(T_psi);
	rho = expit(T_rho);
	w = expit(T_w); 
}

#undef beta
#undef rho
#undef psi
#undef w
#undef pop
#undef S_0
#undef E_0
#undef I_0
#undef R_0
#undef sigma
#undef gamma
#undef intervention
#undef T_beta
#undef T_rho
#undef T_psi
#undef T_w
#undef T_pop
#undef T_S_0
#undef T_E_0
#undef T_I_0
#undef T_R_0
#undef T_sigma
#undef T_gamma
#undef T_intervention

static int __pomp_load_stack = 0;

void __pomp_load_stack_incr (void) {++__pomp_load_stack;}

void __pomp_load_stack_decr (int *val) {*val = --__pomp_load_stack;}

void R_init_snippets (DllInfo *info)
{
R_RegisterCCallable("snippets", "__pomp_load_stack_incr", (DL_FUNC) __pomp_load_stack_incr);
R_RegisterCCallable("snippets", "__pomp_load_stack_decr", (DL_FUNC) __pomp_load_stack_decr);
R_RegisterCCallable("snippets", "__pomp_rinit", (DL_FUNC) __pomp_rinit);
R_RegisterCCallable("snippets", "__pomp_stepfn", (DL_FUNC) __pomp_stepfn);
R_RegisterCCallable("snippets", "__pomp_rmeasure", (DL_FUNC) __pomp_rmeasure);
R_RegisterCCallable("snippets", "__pomp_dmeasure", (DL_FUNC) __pomp_dmeasure);
R_RegisterCCallable("snippets", "__pomp_to_trans", (DL_FUNC) __pomp_to_trans);
R_RegisterCCallable("snippets", "__pomp_from_trans", (DL_FUNC) __pomp_from_trans);
}
