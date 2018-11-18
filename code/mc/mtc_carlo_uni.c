/*
  Compton scattering in a uniform, hot  plasma cloud.
  Based on Gorecki & Wilczewski and Pozdnyakov, Sobol & Sunyaev

   22 I 1992
   modifications: Jan/Feb 1998, Durham
                  - sigma_mean corrected - now it seems to be OK.
                  - light curve constructed
		  - externally irradiated spherical source geometry added
		  - Sunaev-Titarchuk spatial distribution added
		  Dec 1998
		  - cold reflection added
		  Jan 1999
		  - major change: distances no longer in tau units but in Rg !!
		  March 1999
                  - bulk motion added
		  March 1999
		  - active region geometry added (hemisphere, H, R)
                           (spat_distr = 5)
                  June 1999
                  - bulk motion corrected: <sigma> computed properly
                  Dec 1999 
                  - arbitrary opacites for reflection  
 
   Has to be compiled with quadrat.o


   2010: Lense-Thirring precesion project
		  

   Description of the parameters:
   spatial_distribution: 
   1 - sphere
   2 - slab
   3 - disk + sphere
   4 - spherical shell (note that it might be incorrect, because dempty is 
                                                    always added to lambda
   5 - active region
   6 - externally illuminated torus (Lense-Thirring)

*/

#include "mtc_incl_def.c"

double  tTh1Rg;
  
/* for irradiation emissivity */
double rhist[1000],r_hit_max = 2;
int    nrh = 200;


int debug,debug0,debug_cond;

double reflect[1000][10][5], incid[1000];
double sigma_KN();
double sigma_bf();

/* for cold opacities */ 
#define MAXA 1393 
 
/* for other opacities (from mansig) */ 
/* #define MAXA 459  */


double enop[MAXA],opac[MAXA],elop,coop;
int    jlo_op = 1;

double Rinner;

double  total_output_up = 0,
        total_output_down = 0,
        total_input = 0,
        total_refl = 0;


main(argc,argv)
  int    argc;
  char   *argv[];
{
  long int  il,istep;
  double    position[3],momentum[3],
            energy,weight,speed,r;
  int       runnum,i,j,k;

  double    ener0,ener_tot,ener_esc,ener0_tot,ener_esc_tot; 
  double    xf[200],yf[200],a,b,siga,sigb,chi2,e,tnew;
  int       nos;
  FILE *fp;

  for(i=0; i<1000; i++) {
    rhist[i] = 0;
  }

  for(i=0; i<1000; i++) {
    incid[i] = 0;
    for(j=0; j<10; j++) {
      for(k=0; k<5; k++)
	reflect[i][j][k] = 0;
    }
  }

  runnum = 0;
  if (argc == 2)   {
    runnum = atoi(argv[1]);
  }

  input_data(runnum);
  set_accuracy(1e-5);

  if (do_reflection) {
    read_opacity();
  }

  debug0 = 98097;

  define_variables();
  istep = mc_steps/5;
  if (istep==0) istep=1;
  debug = 0;
  debug_cond = 0;
  for (il=0; il<=mc_steps; il++) {

#ifdef DEBUG
    if (debug_cond)
      printf("%ld   =================================== \n",il);
    debug = il;
#endif

    if ((il % istep) == 0) { 
      /*        output_results(runnum,il); */
      writelog(runnum,il); 
    }

    do {
      generate_photon(position,momentum,&energy,&weight);
    } while (!inside(position));
    
#ifdef DEBUG
      if (debug_cond) {
	printf("generated: \n");
	writeposdir(position,momentum);
      }
#endif

    next_scattering(position,momentum,energy,&weight,&ener_esc);
#ifdef DEBUG
      if (debug_cond) {
	printf("next scat: \n");
	writeposdir(position,momentum);
      }
#endif

    do  {
      scattering(position,momentum,&energy,weight);

#ifdef DEBUG
      if (debug_cond) {
	printf("scat: \n");
	writeposdir(position,momentum);
      }
#endif

      next_scattering(position,momentum,energy,&weight,&ener_esc);

#ifdef DEBUG
      if (debug_cond) {
	printf("next scat: \n");
	writeposdir(position,momentum);
      }
#endif

    } while ((weight > wghtmin) && (energy > emin_o) && (energy < emax_o));
  }
  
  output_results(runnum,il);
  speed = (1.0*mc_steps)/(time(NULL)-start_time);
  printf(" Speed of computations: %.1f photons/sec \n",speed);
  
  /*  fp = fopen("rhist.dat","w");
  for(i=0; i<=nrh; i++) {
    r = pow(10.,(i+0.5)*r_hit_max/nrh);
    fprintf(fp,"%lf %lf \n",r,rhist[i]/mc_steps*nrh/r);
  }
  fclose(fp);
  */

  /* 2-20 spectral index */

  nos = num_bin_o*log10(emax_o/emin_o);

  j = 0;
  for(i=0; i<nos; i++) {
    e = 511*emin_o*pow(10.,(i+0.5)/num_bin_o );
    if ( (e>2) && (e<20) ) {
      j ++;
      xf[j] = log(e);
      yf[j] = log(out_s[i]);
    }
  }
  fit(xf,yf,j,&a,&b,&siga,&sigb,&chi2);
  printf("2-20 energy spectral index: %.3f, chi^2  %.3e\n",b,chi2);
  

  /*  printf("i: %.3e,   o_u: %.3e,   o_d: %.3e,   r: %.3e \n",total_input,
	 total_output_up,total_output_down,total_refl);
  tnew = rad_temp*511*pow((total_output_down-total_refl)/total_input,0.25);
  printf("old temp: %.3f,  new temp:  %.3f \n",rad_temp*511,tnew);
  */

  return 0;
}

/*end main*/


writelog(runnum,noph)
   long int noph;
        int runnum;
{
    FILE   *fl;
    double speed;
    char   log_file[20],num[10];

    strcpy(log_file,"mtc__log.000000");
    sprintf(num,"%06i",runnum);
    memcpy(log_file+9,num,6);

    if (noph == 0) {
      fl = fopen(log_file,"w");
      fprintf(fl,"Uniform density run.  %s \n",
	      asctime(localtime(&start_time)));
      fprintf(fl,
	      "Spatial distribution: %d (%.1f),  Rin: %.2f,  rndseed: %ld \n",
	      spat_distr,r_emission,Rin,idum0);
      fprintf(fl,"Tau: %.2f,  plasma temp (keV): %.1f,  T0 (keV): %.3f \n",
	      tau_max,511./med_temp,rad_temp*511.);
    } else {
      fl = fopen(log_file,"a");
      speed = (1.0*noph)/(time(NULL)-start_time);
      printf("%6li  %.1f \n",noph,speed);
      fprintf(fl,"%6li  %.1f \n",noph,speed);
    }
    fclose(fl);
    
    return 0;
}

/*end log*/



next_scattering(pos,mom,energy,weight,ener_esc)
   double pos[],mom[],energy,*weight,*ener_esc;
/* 
    determines position of the next scattering
    sets global variable:   energ2
*/
{
   double p_esc,p_n_esc,sigm_mean,wesc,
          l_i,x,lambda,ek,dempty,
          v1 = 1e-5,
          v2 = 1 - 1e-5;
   int    i,indx,dindex,iindx,nempty,findx;
   static int ins=1, inc=1;

   double d,r_r,mom_r[3],pos_r[3],te[5],ee,ww,e_sm,xe,xt,fi,post[3],momt[3];

   l_i = distance(pos,mom,&nempty,te); /* distance to the boundary, along 
                                        the direction of photon's motion */

   if (l_i < 1e-9) {
     printf("l_i = %.4lf \n",l_i);
   }

   energ2 = 2*energy;         /* global variable used by sigma_integr_   */

   /* very old method:
      sigm_mean = qsimp(sigma_integr_,v1,v2); */   /*    (8) in G&W    */

   if (fabs(b_bulk) > 1e-6) {
     e_sm = energy>sig_e[1]?
       (energy<sig_e[nop_se]?energy:sig_e[nop_se]*0.99999):
       sig_e[1]*1.000001;
     c_hunt(sig_e ,nop_se,e_sm  ,&ins);
     c_hunt(sig_ct,nop_ct,mom[2],&inc);
     xe = 1-log(e_sm/sig_e [ins])/log(sig_e [ins+1]/sig_e [ins]);
     xt = 1-(mom[2]-sig_ct[inc])/(sig_ct[inc+1]-sig_ct[inc]);
     
     if ( (xe<0) ||  (xe>1) ) printf("xe=%lf %lf %d\n",xe,e_sm*511,ins);
     if ( (xt<0) ||  (xt>1) ) printf("xt=%lf \n",xt);

     sigm_mean =   xe *   xt * sigma_mean[ins  ][inc  ] + 
                (1-xe)*   xt * sigma_mean[ins+1][inc  ] +
                (1-xe)*(1-xt)* sigma_mean[ins+1][inc+1] +
                   xe *(1-xt)* sigma_mean[ins  ][inc+1];
   } else {
     c_hunt(sig_e,nop_se,energy,&ins);
     sigm_mean = sigmm[ins]+(energy-sig_e[ins])*(sigmm[ins+1]-sigmm[ins])/
                                                (sig_e[ins+1]-sig_e[ins]);
   }
   

   /* l_i is now in Rg !! */
   p_esc = exp(-l_i*tTh1Rg*sigm_mean);        /*    (7) in G&W           */
#ifdef DEBUG
   if (debug_cond) 
     printf("d: %.3f, sig: %.3f  p_esc: %.3e\n",l_i,sigm_mean,p_esc);
#endif
   x = rnd()*(1 - p_esc);                     /* formula (13) in G&W     */

   wesc = p_esc*(*weight);
   *ener_esc = wesc*energy*511;

   if (nop_grid) {
     if ( (energy>ear[0]) && (energy<ear[nop_grid])) {
       c_hunt(ear,nop_grid,energy,&indx);
       photar[indx] += wesc;
     }
   }
   indx = num_bin_o*log10(energy/emin_o);     /* add escaping photons    */
   if (indx>=0)   {
     ek = 511*energy;
     dempty = te[2]-te[1] + te[4]-te[3];
     dindex = floor((total_dist+l_i+dempty)/total_dist_max*dist_points);
     if (dindex > MAXd) dindex = MAXd;

#ifdef DEBUG
     if (debug_cond) 
       printf("esc: dindx: %d, ek: %.2f num: %d\n",dindex,ek,num_dist_en);
#endif
     if ( (spat_distr != 2) || ((spat_distr==2) && (mom[2]>0))) {
       for(i=0; i<num_dist_en; i++) {
	 if ( (ek > dist_energies[0][i]) && (ek <= dist_energies[1][i])) {
	   dist_hist[i][dindex] += wesc;
	   scatt_hist[i][nscatt] += wesc;
	 }
       }
     }
#ifdef DEBUG
     if (debug_cond) printf("Escape: ek: %.2f  indx %d,  wesc  %.2e \n",ek,indx,wesc);
#endif

     switch (spat_distr) {
     case 0 : 
     case 1 :
     case 4 : {
       out_s[indx] += wesc;
       for(i=0; i<max_spsc; i++) {
	 if ((nscatt>=spsc[0][i]) && (nscatt<=spsc[1][i])) {
	   scatt_spec[indx][i] += wesc;
	 }
       }
       break ;
     }
     case 2 :   /* slab */ 
     case 5 :   /* active region */
       {
       if (mom[2]>0) {
	 out_s[indx] += wesc;    /* transmitted photons     */
	 /* 	angle of view dependence */
	 iindx = mom[2]*nangle;
	 if (nscatt>0) out_si[indx][iindx][0] += wesc;
	 for(i=0; i<max_spsc; i++) {
	   if ((nscatt>=spsc[0][i]) && (nscatt<=spsc[1][i])) {
	     scatt_spec[indx][i] += wesc;
	   }
	 }
       } else {
	 out_r[indx] += wesc;    /* reflected photons       */
       }
       break;
     }
     case 3 : {
       out_s[indx] += wesc; 
       /* 	angle of view dependence */
       iindx = fabs(mom[2])*nangle; /* both up- and downwards going photons */
       out_si[indx][iindx][0] += wesc;
       for(i=0; i<max_spsc; i++) {
	 if ((nscatt>=spsc[0][i]) && (nscatt<=spsc[1][i])) {
	   scatt_spec[indx][i] += wesc;
	 }
       }
       break ;
     }
     case 6 : {   /* torus */
       /* store photons escaping in MAXphi (4) different directions */

       if (precession_option<2) {
	 prectransfsimple(pos,mom,post,momt,prectheta,precphi,-1); 
       } else {
	 prectransf(pos,mom,post,momt,prectheta,precphi,-1); 
       }

       fi = atan2(momt[1],momt[0])*180/pi;
       if (fi < 0) { fi += 360; }
       if (fi >=360.) {
	 printf("FI > 360 !!! %6.2lf \n",fi);
       }
       if (momt[2]>0) { /* only upward-escaping photons are included */
	 if (!hits_outer(post,momt)) { /* shielding by the outer disk 
					  taken into account */

	   /* bins in \phi: 0-20, 90-110, 180-200, 270-290 degs */
	   /* precession axis at 10 degs ?? - is this true      */
	   findx = -1;
	   for(i=0; i<MAXphi; i++) {
	     if (fabs(fi-(i*90+10)) <= 10) {  
	       findx = fi/90.;
	     }
	   }
	   
	   out_s[indx] += wesc; 
	   /* 	angle of view dependence */
	   if (findx >=0) {  /* store only photons escaping within narrow 
				ranges of \phi angle   */
	     iindx = momt[2]*nangle; 
	     out_si[indx][iindx][findx] += wesc;
	     for(i=0; i<max_spsc; i++) {
	       if ((nscatt>=spsc[0][i]) && (nscatt<=spsc[1][i])) {
		 scatt_spec[indx][i] += wesc;
	       }
	     }
	   }
	 }
       }
       break ;
     }
     }

   }

   if (do_reflection) {

     /* cold reflection? */
     if ((spat_distr == 6) && (post[2]>0) && (momt[2] < 0)) {
       /* disk + TORUS geometry */
       d = -post[2]/momt[2];
       r_r = sqrt(pow(post[0]+d*momt[0],2.)+pow(post[1]+d*momt[1],2.));
       /*     printf("1: %6.3f  %6.3f  %6.3f \n",d,r_r,Rinner); */
       if (r_r > Rdisk) {
	 for(i=0; i<3; i++) {
	   pos_r[i] = post[i]+d*momt[i];
	   mom_r[i] = momt[i];
	 }
	 if (fabs(pos_r[2])>1e-5) {
	   printf("inaccurate pos_r:  %.3e \n",pos_r[2]);
	 }
	 /* local coordinate system in cold_reflection requires p_z<0 initially */
	 mom_r[2] = -fabs(mom_r[2]);
	 ee = energy;
	 if (r_r > Rout) {
	   ww = wesc;
	 } else {
	   ww = (*weight) * exp(-d*tTh1Rg*sigm_mean);
	 }
  	 indx = num_bin_o*log10(ee/emin_o);
	 incid[indx] += ww;
	 cold_reflection(pos_r,mom_r,ee,ww);
       }
     }
     if ((spat_distr == 3) && (pos[2]*mom[2] < 0)) {
       /* disk + sphere geometry */
       d = -pos[2]/mom[2];
       r_r = sqrt(pow(pos[0]+d*mom[0],2.)+pow(pos[1]+d*mom[1],2.));
       /*     printf("1: %6.3f  %6.3f  %6.3f \n",d,r_r,Rinner); */
       if (r_r > Rinner) {
	 for(i=0; i<3; i++) {
	   pos_r[i] = pos[i]+d*mom[i];
	   mom_r[i] = mom[i];
	 }
	 if (fabs(pos_r[2])>1e-5) {
	   printf("inaccurate pos_r:  %.3e \n",pos_r[2]);
	 }
	 /* local coordinate system in cold_reflection requires p_z<0 initially */
	 mom_r[2] = -fabs(mom_r[2]);
	 ee = energy;
	 if (r_r > Rout) {
	   ww = wesc;
	 } else {
	   ww = (*weight) * exp(-d*tTh1Rg*sigm_mean);
	 }
  	 indx = num_bin_o*log10(ee/emin_o);
	 incid[indx] += ww;
	 cold_reflection(pos_r,mom_r,ee,ww);
       }
     }
     
     /* cold reflection for slab/active region geometry - 
	no subsequent comptonization of the reflected photons */
     if ( ( (spat_distr == 2) && (mom[2] < 0) ) || 
	  ( (spat_distr == 5) && (mom[2] < 0) ) ) {
       d = -pos[2]/mom[2];
       for(i=0; i<3; i++) {
	 pos_r[i] = pos[i]+d*mom[i];
	 mom_r[i] = mom[i];
       }
       if (fabs(pos_r[2])>1e-5) {
	 printf("inaccurate pos_r:  %.3e \n",pos_r[2]);
       }
       
       ee = energy;
       ww = wesc;
       indx = num_bin_o*log10(ee/emin_o);
       incid[indx] += ww;
       cold_reflection(pos_r,mom_r,ee,ww);
     }

   }

   
   p_n_esc = 1 - p_esc;                       /* remain probability      */
   *weight *= p_n_esc;                        /* new weight              */
   lambda = -log(1-x)/(sigm_mean*tTh1Rg);     /* free path of photon     */

   if (lambda > te[1]) {                     /* add possible contribution */
     lambda += te[2]-te[1];                    /*  from the inner void      */
     if (lambda > te[3]) {
       lambda += te[4]-te[3];
     }
   }
					
   total_dist += lambda;

   for (i=0; i<3; i++)
     pos[i] += lambda*mom[i];                 /* position of the next    */
                                              /* collision               */
   return 0;
}
/*end next_scatter*/




hits_outer(pos,mom)
double pos[],mom[];
{
  int hit;
  double t,xh,yh,rh,zh;

  if (fabs(mom[2]) < 1e-18) {
    return 1;
  }

  if (pos[2]*mom[2]>0) {
    hit = 0;
  } else {
    t = fabs(pos[2]/mom[2]);
    xh = pos[0]+t*mom[0];
    yh = pos[1]+t*mom[1];
    /*    zh = pos[2]+t*mom[2]; */
    rh = sqrt(xh*xh+yh*yh);
    hit = rh>Rdisk;
  }

  return hit;
}
/*end hits_outer*/



/*define */
define_variables()
{
  double Mass, Rg;

  Mass = 10.;
  Rg = 1.5e5*Mass;


  if (spat_distr != 6) {
    if (rho0 > 1e6) { /* treat rho0 as n_0 */
      tTh1Rg = Rg*rho0*sigma_T;
      Rout = Rin + tau_max/tTh1Rg;
    } else {     /* rho0 means Rout in Rg */
      Rout = rho0;
      tTh1Rg = tau_max/(Rout-Rin);
    }
    Rdisk = Rout;
    printf("Rout: %.1lf,   tau_T/1Rg: %.4lf \n",Rout,tTh1Rg);
  }
  /*  Rout = tau_max + Rin; */
  Rq = Rout*Rout;
  Rin2 = Rin*Rin;
  if (r_emission < 0) {
    Rinner = fabs(r_emission)*Rout;
  } else {
    Rinner = 1e30;
  }
  
  g_bulk = 1/sqrt((1-b_bulk)*(1+b_bulk));

  return 0;
}
/*end define*/



#include "mtc_incl_code.c"

cold_reflection(pos,mom,ee,ww)
     double pos[],mom[],ee,ww;
{

  double mom_i[3],mom_f[3],pos_l[3],
         energy,weight,weight_min,mu_prime,sine2,term,bound,r_hit;
  int    i,rindx;
  
  energy = ee;
  weight = ww;

  for (i=0; i<3; i++) {
    mom_i[i] = mom[i];
    pos_l[i] = pos[i];
  }

  r_hit = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
  rindx = log10(r_hit/(fabs(r_emission)*Rout))/r_hit_max*nrh;
  if (rindx>=0) 
    rhist[rindx] += weight;

  weight_min = weight*1e-5;
  if (weight_min < wghtmin) weight_min = wghtmin;

  /*    printf("%.3e  %.3e \n",energy*511,weight/weight_min); */

  scatpos(pos_l,mom_i,energy,&weight);

  do {

    do {
      mu_prime = 1-(exp(log(1+2*energy)*rnd())-1)/energy;
      sine2 = 1-mu_prime*mu_prime;
      term  = 1/(1+energy*(1-mu_prime));
      bound = term*(term+1/term-sine2)/2;
    } while (rnd() > bound);
    
    transform(mom_i,mom_f,sqrt(sine2),mu_prime);
    for (i=0; i<3; i++)
      mom_i[i] = mom_f[i];
    
    energy *= term;
    
    scatpos(pos_l,mom_i,energy,&weight);
  }
  while ( (weight > weight_min) && (energy > emin_o) );

  return 0;
}

scatpos(pos,mom,energy,weight)
     double pos[],mom[],energy,*weight;
{
  double sigma_es,sigma_abs,sigma_tot,p_esc,wesc,lambda,d,fi;
  int    indx,iindx,i,findx;

  sigma_es = sigma_KN(energy);
  sigma_abs= sigma_bf(energy);
  sigma_tot= sigma_es+sigma_abs;
    
  if (mom[2] <= 0) { /* downwards - no escape */
    p_esc = 0;
  } else {
    d = fabs(pos[2]/mom[2]);
    p_esc = exp(-d*sigma_tot);
    wesc = p_esc*(*weight);
    
    indx = num_bin_o*log10(energy/emin_o);     /* add escaping photons    */
    if (indx>=0) {
      iindx = mom[2]*nangle + 1;

      fi = atan2(mom[1],mom[0])*180/pi;
      if (fi < 0) { fi += 360; }
      findx = -1;
      for(i=0; i<MAXphi; i++) {
	if (fabs(fi-(i*90+10)) <= 10) {  
	  findx = fi/90.;
	  findx ++;
	}
      }

      if (findx > 0) {
	reflect[indx][0][0] += wesc;   /* fi-indx=0 - all photons */
	reflect[indx][0][findx] += wesc;   /* fi-indx > 0 - fi-bins */
	/* 	angle of view dependence */
	reflect[indx][iindx][findx] += wesc;
      }
    }
  }
  lambda = -log(1-(1-p_esc)*rnd())/sigma_tot;
  
  for (i=0; i<3; i++)
    pos[i] += lambda*mom[i];                 /* position of the next    
						scattering               */  
  (*weight) *= (1-p_esc)*sigma_es/sigma_tot;
    
    return 0;
}

transform(omega_i,omega_f,st,ct)
    double omega_i[],omega_f[],
           st,ct;
{
    double x,y,r2,o3,dd,cph,sph,norm;
    int    i;

    do {
      x = rnd()-0.5;
      y = rnd()-0.5;
      r2 = x*x+y*y;
    }
    while (r2 > 0.25);
    o3 = 1 - omega_i[2]*omega_i[2];
    dd = sqrt(o3*r2);
    cph = x/dd;
    sph = y/dd;
    omega_f[0] = st*( omega_i[1]*sph-omega_i[0]*omega_i[2]*cph)+omega_i[0]*ct;
    omega_f[1] = st*(-omega_i[0]*sph-omega_i[1]*omega_i[2]*cph)+omega_i[1]*ct;
    omega_f[2] = st*cph*o3+omega_i[2]*ct;

    norm = 1/sqrt(omega_f[0]*omega_f[0]+omega_f[1]*omega_f[1]+omega_f[2]*omega_f
[2]);
    if (fabs(norm-1.0) > 0.001)
    {
      for (i=0; i<3; i++)
         omega_f[i] *= norm;
    }

    return 0;
}


double sigma_bf(e)
   double e;
{
   double  sigma;

   if (e >= elop) {
     sigma = coop/(e*e*e);
   } else {
     if (e > enop[0]) {
       c_hunt(enop,MAXA-1,e,&jlo_op);
       sigma = opac[jlo_op]+(e-enop[jlo_op])*(opac[jlo_op+1]-opac[jlo_op])/
	                                     (enop[jlo_op+1]-enop[jlo_op]);
     } else {
       sigma = opac[0];
     }
   }

   return(sigma);
}

double sigma_KN(e)
     double e;
{
  double  sigma,z;
/*
   Klein-Nishina cross section in units of sigma(Thomson),
   e is energy in units of the electron rest mass
   Per hydrogen atom! It means multiplying it by 1.21, which is the
   number of electrons per H atom for cosmic abundances

*/  
  z = 1+2*e;
  if (e < 0.02)
    sigma = (1+e*(2+e*(1.2-0.5*e)))/z/z;  /*  Hubbell 1980  */
  else
    sigma = 0.375/e*((1-2/e-2/e/e)*log(z)+0.5+4/e-0.5/z/z);

  sigma *= 1.21; 
  return(sigma);
}

read_opacity()
{
   double ee,dum,op;
   int    i;
   FILE  *fp,*fopen();

   if ((fp=fopen("absorp.dat","r")) != NULL)
   {
      for (i=0; i<MAXA; i++)
      {
        fscanf(fp,"%lf %lf %lf",&ee,&op,&dum);
  	enop[i] = ee/511;
        opac[i] = op;
      }
      fclose(fp);

      elop = enop[MAXA-1];
      coop = opac[MAXA-1]*pow(elop,3.0);
      /*      printf("opac: %.3e  %.4e \n",elop*511,coop); */
   }
   else
   {
      printf("Cannot find file with opacity data \n");
      exit(13);
   }

   return 0;
}

read_opacity_ion() 
{ 
   double ee,op; 
   int    i; 
   FILE  *fp,*fopen(); 
 
   if ((fp=fopen("opac_xi1e4.dat","r")) != NULL) 
   { 
      for (i=0; i<MAXA; i++) 
      { 
        fscanf(fp,"%lf %lf ",&ee,&op); 
        enop[i] = ee/511; 
        opac[i] = op; 
      } 
      fclose(fp); 
 
      elop = enop[MAXA-1]; 
      coop = opac[MAXA-1]*pow(elop,3.0); 
      /*      printf("opac: %.3e  %.4e \n",elop*511,coop); */ 
   } 
   else 
   { 
      printf("Cannot find file with opacity data \n"); 
      exit(13); 
   } 
 
   return 0; 
} 



fit(double x[], double y[], int ndata, double *a,
	 double *b, double *siga, double *sigb, double *chi2)
{
  int i;
  double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
  
  *b=0.0;
  for (i=1;i<=ndata;i++) {
    sx += x[i];
    sy += y[i];
  }
  ss=ndata;

  sxoss=sx/ss;
  for (i=1;i<=ndata;i++) {
    t=x[i]-sxoss;
    st2 += t*t;
    *b += t*y[i];
  }

  *b /= st2;
  *a=(sy-sx*(*b))/ss;
  *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
  *sigb=sqrt(1.0/st2);
  *chi2=0.0;

  for (i=1;i<=ndata;i++)
    *chi2 += (y[i]-(*a)-(*b)*x[i])*(y[i]-(*a)-(*b)*x[i]);
  sigdat=sqrt((*chi2)/(ndata-2));
  *siga *= sigdat;
  *sigb *= sigdat;
 
  return 0;
}

/* (C) Copr. 1986-92 Numerical Recipes Software . */

