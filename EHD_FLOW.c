/*THIS SERIES OF UDFs IMPLEMENTS TWO USER DEFINED SCALAR TRANSPORT EQUATIONS AND ONE MOMENTUM SOURCE TERM, THE TRANSPORT
  EQUATIONS ARE THE TRANSPORT OF SPACE CHARGE IN A CONTINUUM AND THE POISSON EQUATION FOR VOLTAGE, THE ADDED SOURCE ON
  THE MOMENTUM EQUATION IS THE LORENTZ FORCE TERM, DIELECTRIC AND ELECTROSTRICTIVE FORCES ARE NEGLECTED*/

/*A.KOURMATZIS*/

#include "udf.h"
#include "sg.h"


static real kappa = 0.00000000571075;			/*ionic mobility*/
static real epsilon = .00000000000855;  /*electric permittivity of free space*/
static real epsilon_r=3.54;			   /*relative permittivity*/



/*POISSON EQUATION*/
/*=========================================================================================================*/
DEFINE_SOURCE(electric_source,c,t,dS,eqn)
{
dS[eqn] = 0;
return (C_R(c,t)*C_UDSI(c,t,1)/epsilon);
}


DEFINE_DIFFUSIVITY(electric_diff, c, t, i)
{
return epsilon_r;
}

DEFINE_UDS_UNSTEADY(spcharge_unsteady,c,t,i,apu,su)
{
real physical_dt, vol, rho, phi_old;
physical_dt = RP_Get_Real("physical-time-step");
vol = C_VOLUME(c,t);
rho = C_R_M1(c,t);
*apu = -rho*vol / physical_dt;/*implicit part*/
phi_old = C_STORAGE_R(c,t,SV_UDSI_M1(i));
*su = rho*vol*phi_old/physical_dt;/*explicit part*/
}

/*This is a default scalar flux term used by fluent, see fluent 6.3 udf manual for details*/

DEFINE_UDS_FLUX(spcharge_flux,f,t,i)
{

cell_t c0, c1 = -1;
Thread *t0, *t1 = NULL;
real NV_VEC(psi_vec), NV_VEC(A), flux = 0.0;
c0 = F_C0(f,t);
t0 = F_C0_THREAD(f,t);
F_AREA(A, f, t);
/* If face lies at domain boundary, use face values; */
/* If face lies IN the domain, use average of adjacent cells. */
if (BOUNDARY_FACE_THREAD_P(t)) /*Most face values will be available*/
{
real dens;

/* Depending on its BC, density may not be set on face thread*/
if (NNULLP(THREAD_STORAGE(t,SV_DENSITY)))
dens = F_R(f,t); /* Set dens to face value if available */
else
dens = C_R(c0,t0); /* else, set dens to cell value */
NV_DS(psi_vec, =, F_U(f,t), F_V(f,t), F_W(f,t), *, dens);
flux = NV_DOT(psi_vec, A); /* flux through Face */
}
else
{
c1 = F_C1(f,t); /* Get cell on other side of face */
t1 = F_C1_THREAD(f,t);
NV_DS(psi_vec, =, C_U(c0,t0),C_V(c0,t0),C_W(c0,t0),*,C_R(c0,t0));
NV_DS(psi_vec, +=, C_U(c1,t1),C_V(c1,t1),C_W(c1,t1),*,C_R(c1,t1));
flux = NV_DOT(psi_vec, A)/2.0; /* Average flux through face */
}
/* Fluent will multiply the returned value by phi_f (the scalar's
value at the face) to get the ``complete'' advective term. */
return flux;
}


DEFINE_SOURCE(spcharge_lin,c,t,dS,eqn)
{

real source,source_a,source_b;

source_a=C_R(c,t)*kappa*(C_UDSI_G(c,t,0)[0])*(C_UDSI_G(c,t,1)[0]);

source_b=C_R(c,t)*kappa*(C_UDSI_G(c,t,0)[1])*(C_UDSI_G(c,t,1)[1]);

//source_c=C_R(c,t)*kappa*C_UDSI_G(c,t,0)[2]*C_UDSI_G(c,t,1)[2];

source=source_a+source_b;

dS[eqn] = 0;

return source;

}


DEFINE_SOURCE(spcharge_quad,c,t,dS,eqn)
{
real source;
source = -(pow(C_R(c,t),2)*pow(C_UDSI(c,t,1),2)*kappa)/(epsilon*epsilon_r);
dS[eqn] = -(2*pow(C_R(c,t),2)*C_UDSI(c,t,1)*kappa)/(epsilon*epsilon_r);
return source;
}


DEFINE_DIFFUSIVITY(spcharge_diff, c, t, i)
{
//return .00000001*kappa / 40;  
return C_R(c,t)*kappa / 40;
}


DEFINE_SOURCE(xmom_source,c,t,dS,eqn)
{
real source;
dS[eqn]=0;
source=-C_R(c,t)*C_UDSI(c,t,1)*C_UDSI_G(c,t,0)[0];
return source;

}

DEFINE_SOURCE(ymom_source,c,t,dS,eqn)
{

real source;
dS[eqn]=0;
source=-C_R(c,t)*C_UDSI(c,t,1)*C_UDSI_G(c,t,0)[1];
return source;
}

DEFINE_SOURCE(zmom_source,c,t,dS,eqn)
{

real source;
dS[eqn]=0;
source=-C_R(c,t)*C_UDSI(c,t,1)*C_UDSI_G(c,t,0)[2];
return source;
}








	
