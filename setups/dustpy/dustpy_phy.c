
#include "fargo3d.h"

//TODO: check with astropy
#define kB_cgs 1.380649e-16 // Boltzmman constant
//#define mp_cgs 1.67262192369e-24  // proton mass in g 
#define mu  2.3 // mean molecular weight in units of proton mass
#define mu_code 1.934e-57 // mean molecular weight in code unit (m_mu=2.3m_p)
//#define G  6.674299999999999e-08 // gravitational constants in cgs
//#define msun 1.988409870698051e+33 // solar mass
#define sigma_H2_code 3.31e-45 // the cross section for molecular hydrogen in code units (considering 1 length unit = distance to jupiter)

//define a structure for easier computation of v_visc and pressure_eta
typedef struct {
    int l_left, l_right;

    real r_left, r, r_right;
    real r_face_left, r_face_right;
    real omega_r_left, omega_r, omega_r_right;

} RadGeo;


RadGeo Read_Radial_Geometry(int i, int j, int k){
    
    RadGeo g;

    #ifdef GHOSTSX
        g.l_right = ((i)+(j+1)*(Nx+2*NGHX)+((k)*Stride));
        g.l_left  = ((i)+(j-1)*(Nx+2*NGHX)+((k)*Stride));
    #else
        g.l_right = ((i)+(j+1)*(Nx)+((k)*Stride));
        g.l_left  = ((i)+(j-1)*(Nx)+((k)*Stride));
    #endif

    g.r_left = Ymed(j-1);
    g.r      = Ymed(j);
    g.r_right= Ymed(j+1);

    g.r_face_left = Ymin(j);
    g.r_face_right= Ymax(j);

    g.omega_r_left  = sqrt(G*MSTAR/g.r_left/g.r_left/g.r_left);
    g.omega_r       = sqrt(G*MSTAR/g.r/g.r/g.r);
    g.omega_r_right = sqrt(G*MSTAR/g.r_right/g.r_right/g.r_right);

    return g;
}

real interp(real x_left, real x_right, real y_left, real y_right, real xi){
    // interpolate the value at xi (which is in thie middle of x_right and x_left)
    //  from the known values of y_right and y_left at x_right and x_left

    real yi, slope; 
    slope = (y_right - y_left) / (x_right-x_left);
    yi  = y_left + slope * (xi- x_left);

    return yi;

}

real Particle_Mass(real s_cgs) {

    real m_cgs;
    real m_code;

    m_cgs = 4.0/3.0 * M_PI * RHOSOLID * pow(s_cgs, 3.0);
    m_code = m_cgs * MSTAR/ MSTAR_CGS;

    return  m_code;
}

real Gas_Temperature(real soundspeed_code){
    
    real gas_temp_code;

    gas_temp_code  =  pow(soundspeed_code, 2) / R_MU;
    return gas_temp_code;
}

real Stokes_Number(real s_cgs, real sigma_g_code){

    real s_code;
    real St;

    s_code = R0/R0_CGS * s_cgs; 
    St = M_PI/ 2.0 * s_code * rhograin/ sigma_g_code;
    
    return St;
}

real Gas_Scaleheight(real cs_code, real omega_code){

    real h_g_code;

    h_g_code = cs_code/omega_code; 

    return h_g_code;

}

real Gas_Nu(real turb_alpha, real cs_code, real omega_code){

    real nu_code;

    nu_code = turb_alpha * cs_code * cs_code / omega_code;

    return nu_code;

}

real Dust_Scaleheight(real st, real hg_code, real turb_alpha, real cs_code, real omega_code){

    real hd_code;// st for stokes number

    hd_code = hg_code * sqrt(turb_alpha / (turb_alpha + st));

    return hd_code;
}

real Volume_Rho_Midplane(real sigma_g, real hg_code){

    real rho_midplane, hg_code;

    rho_midplane = sigma_g / (sqrt(2.0*M_PI) * hg_code);

    return rho_midplane;

}

real Pgas(real rho_g_mid, real cs_code){

    real pgas

    pgas = rho_g_mid * pow(cs_code,2);

    return pgas;
}

real Pressure_Eta(real *sigma, real *cs, real hg_code, RadGeo g){

    real pgas_left, pgas, pgas_right;
    real pgas_left_face, pgas_right_face;
    real dlogp_dlogr, eta; 

    pgas_left = Pgas(sigma[g.l_left], g.omega_r_left, cs[g.l_left]);
    pgas      = Pgas(sigma[l], g.omega_r, cs[l]);
    pgas_right= Pgas(sigma[g.l_right], g.omega_r_right, cs[g.l_right]);

    pgas_left_face = interp(g.r_left, g.r, pgas_left, pgas, g.r_face_left);
    pgas_right_face= interp(g.r, g.r_right, pgas, pgas_right, g.r_face_right);

    dlogp_dlogr = (log(pgas_right_face) - log(pgas_left_face)) / (log(g.r_face_right)-log(g.r_face_left));
    
    eta = -0.50 * pow(hg_code/g.r, 2) * dlogp_dlogr; 

    return eta;

}

real Vel_Settle(real st, real hd_code, real cs_code, real omega_code, real turb_alpha){

    real st, small; // st for stokes number
    real v_settle_code;

    small = fmin(0.5, st);
    v_settle_code = - hd_code * omega_code * small;
    
    return v_settle_code;
}

real Vel_Visc(real *sigma, real *cs, real turb_alpha, RadGeo g){

    real arg_left, arg, arg_right;
    real arg_left_face, arg_right_face;
    real grad_arg, v_visc;

    arg_left = sigma[g.l_left] * sqrt(g.r_left) * Gas_nu(turb_alpha, cs[g.l_left], g.omega_r_left);
    arg      = sigma[l] * sqrt(g.r) * Gas_nu(turb_alpha, cs[l], g.omega_r);
    arg_right= sigma[g.l_right] * sqrt(g.r_right) * Gas_nu(turb_alpha, cs[g.l_right], g.omega_r_right);

    // the value at the face cell on the left and right of the center cell l
    arg_left_face = interp(g.r_left, g.r, arg_left, arg, g.r_face_left);
    arg_right_face = interp(g.r, g.r_right, arg, arg_right, g.r_face_right);

    grad_arg = (arg_right_face - arg_left_face) / (g.r_face_right-g.r_face_right);

    v_visc = -3.0 / (sigma_g[l] * sqrt(g.r)) * grad_arg;

    return v_visc;

}

real Vel_Gas_Rad(real v_visc, real p_eta, real vk, real A, real B){


    // A and B are two parameters to account for the dust back-reaction;
    // TODO: how is the dust back-reaction accounted for in FARGO3D?
    // If no back reaction, then A = 1.0 and B = 0.0;
    real v_gas_r;

    v_gas_r= A * v_visc + 2.0 * B * p_eta * vk;

    return v_gas_r;

}

real v_drift_max(real v_visc, real p_eta, real A, real B, real vk){

    real v_drift_max;

    v_drift_max = 0.50 * B * v_visc - A * p_eta * vk;

    return v_drift_max;
}

real Vel_Dust_Rad(real v_gas_r, real v_drift_max, real vk, real A, real B, real turb_alpha, real st){

    real v_dust_r; 

    v_dust_r = (v_gas_r  + 2 * v_drift_max * st) / (1 + st*st);

    return v_dust_r;

}
real Rel_Vel_Brownian(real s1_cgs, real s2_cgs, real cs_code) {

    real m1_code, m2_code, T_gas_code;
    real kB_code;
    real dv_brown_code;


    m1_code = Particle_Mass(s1_cgs);
    m2_code = Particle_Mass(s2_cgs);

    T_gas_code = Gas_Temperature(cs_code);

    kB_code = kB_cgs * R_MU/R_MU_CGS * MSTAR/MSTAR_CGS;

    dv_brown_code = sqrt(8.0 * kB_code * T_gas_code/M_PI * (1.0/m1_code+1.0/m2_code));

    return dv_brown_code;
}

real Rel_Vel_Settle(real v_settle1, real v_settle2){

    real rel_vel_settle;

    rel_vel_settle = v_settle1 - v_settle2;

    return rel_vel_settle;

}

real Rel_Vel_Rad(real st1 , real st2, real v_gas_r, real v_drift_max){

    real rel_vel_rad;

    rel_vel_rad = v_gas_r * (1 / (1 + pow(st1, 2)) - 1 / (1 + pow(st2, 2)))
    + 2 * v_drift_max * (st1 / (1 + pow(st1, 2)) - st2 / (1 + pow(st2, 2)));

    return rel_vel_rad
   
}


// azimuthal velocity offset (after cancelling out the keplerian velocity and the frame velocity)
real Vel_Dust_Azi_Offset(real v_drift_max, real st){

    real st, v_dust_a;

    v_dust_a = v_drift_max / (1+ st*st);

    return v_dust_a;
}

real Rel_Vel_Azi(real st1, real st2, real v_drift_max){

    real rel_vel_azi;

    rel_vel_azi = v_drift_max * (1 / pow(st1, 2) - 1 / pow(st2, 2));

    return rel_vel_azi;
}

real Rel_Vel_Turb(real omega_code, real sigma_g_code, real cs_code, real turb_alpha, 
    real st1, real st2){

    const real c0  = 1.6015125;
    const real c1  = -0.63119577;
    const real c2  =  0.32938936;
    const real c3  = -0.29847604;
    const real ya  =  1.6;
    
    real yap1inv  = 1.0/ (1.0 + ya);
    real OmkInv   = 1.0 / omega_code;
    //TODO: check the definition of sigma_H2_code!
    real Re       = 0.5 * turb_alpha * sigma_g_code * sigma_H2_code / mu_code;
    real ReInvSqrt= sqrt(1.0/Re);

    real StL, StS, tauL, tauS, eps;
    real vn, vs, ts, vg2;
    real ys, h1, h2;
    real vrel_d_turb2, vrel_d_turb;
    
    vn            = sqrt(turb_alpha) * cs_code;
    vs            = pow(Re, -0.25) * vn;
    ts            = OmkInv * ReInvSqrt;
    vg2           = 1.5 * pow(vn,2);
    
    StL = fmax(st1, st2);
    StS = fmin(st1, st2);
    eps = StS/StL;

    tauL = StL * OmkInv;
    tauS = StS * OmkInv;

    ys = c0 + c1*StL + c2* pow(StL, 2.0) + c3* pow(StL, 3.0);
    

    if (tauL < 0.2 * ts){
        vrel_d_turb2 = 1.5 * pow(vs / ts * (tauL - tauS), 2.0);
    }
    else if (tauL* ya < ts) {
        vrel_d_turb2 = vg2 * (StL - StS) / (StL + StS) 
        * (pow(StL, 2.0) / (StL + ReInvSqrt) 
        - pow(StS, 2.0)  / (StS + ReInvSqrt));
    }
    else if (tauL < 5.0*ts){
        h1 = (StL - StS) / (StL + StS) * (StL * yap1inv - pow(StS, 2.0) / (StS + ya * StL));
        h2 = 2.0 * (ya * StL - ReInvSqrt) + StL * yap1inv - pow(StL, 2.0) / (StL + ReInvSqrt) + 
        pow(StS, 2.0) / (ya * StL + StS) - pow(StS, 2.0) / (StS + ReInvSqrt);
        
        vrel_d_turb2 = vg2 * (h1 + h2);
    }
    else if (tauL < 0.2 * OmkInv){
        vrel_d_turb2 = vg2 * StL * 
        (2.0*ya - 1.0 - eps + 2.0/ (1.0 + eps) 
        * (yap1inv + pow(eps, 3.0) / (ya+eps)));
    }
    else if (tauL < OmkInv){
        vrel_d_turb2 = vg2 * StL 
        * (2.0 * ys - 1.0 - eps 
        + 2.0/(1.0 + eps) * (1.0/(1.0 + ys)
        + pow(eps,3) / (ys + eps)));
    }

    else if (tauL >= OmkInv){
        vrel_d_turb2 = vg2 * (2.0 + StL + StS)/(1.0 + StL + StS + StL * StS);
    }
    
    vrel_d_turb = sqrt(vrel_d_turb2);

    return vrel_d_turb;
}

real Rel_Vel_Tot(real rel_vel_azi, real rel_vel_brownian, real rel_vel_rad, 
    real rel_vel_turb, real rel_vel_settle){

        real rel_vel_tot;

        rel_vel_tot = sqrt(pow(rel_vel_azi, 2.0) + pow(rel_vel_brownian, 2.0) 
        + pow(rel_vel_rad, 2.0) + pow(rel_vel_turb, 2.0) + pow(rel_vel_settle, 2.0));

        return rel_vel_tot;

}

