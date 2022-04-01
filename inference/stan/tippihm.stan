data{
  int NC;                       /* number of countries */
  int NST;                      /* total number of sites */
  vector[NST] STi;              /* site time in intervention */
  vector[NST] STb;              /* site time in baseline */
  int DXi[NST];                 /* counts in intervention */
  int DXb[NST];                 /* counts in baseline */
  int st[NC];                   /* county index starts */
  int nd[NC];                   /* county index ends */
  /* prior parms */
  real lcmsig;                  /* effect mean prior variance */
  real lcssig;                  /* effect variance prior variance */
  real ratessig;                /* rate prior variance */
  real sigsig;                  /* variance prior variance */
}
parameters{
  real<lower=0> sig;            /* within country variance */
  real lcm;                     /* effect mean (log) */
  real<lower=0> lcs;            /* effect variance parm */
  vector[NC] lcrr;                /* country effects (log) */
  vector[NST] lsrr;                /* site effects (log) */
  vector<lower=0>[NST] rates;      /* site rates */
}
transformed parameters{
  vector<lower=0>[NST] lamb;
  vector<lower=0>[NST] lami;
  lamb = rates .* STb;
  lami = exp(lsrr) .* rates .* STi;
}
model{
  /* site rate prior */
  rates ~ cauchy(0,ratessig);
  /* top level prior */
  lcm ~ normal(0,lcmsig);       /* effect prior */
  lcs ~ cauchy(0,lcssig);       /* variance prior */
  sig ~ cauchy(0,sigsig);       /* variance prior */
  /* country rate ratios */
  lcrr ~ normal(lcm,lcs);
  /* site rate ratios */
  for(i in 1:NC){                       /* TODO consider country sig */
    lsrr[st[i]:nd[i]] ~ normal(lcrr[i],sig); /* site in country i */
  }
  /* data likelihood */
  DXb ~ poisson(lamb);          /* baseline */
  DXi ~ poisson(lami);          /* intervention */
}
