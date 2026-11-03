

typedef struct {

    int l_left, l_right;

    real r_left, r, r_right;
    real r_face_left, r_face_right;

    real omega_r_left, omega_r, omega_r_right;
} RadGeo;

RadGeo Read_Radial_Geometry(int i, int j, int k);

real interp(real x_left, real x_right,
            real y_left, real y_right, real xi);

real Particle_Mass(real s_cgs);

real Gas_Temperature(real soundspeed_code);

real Stokes_Number(real s_cgs, real sigma_g_code, real rhograin);

real Gas_Scaleheight(real cs_code, real omega_code);

real Gas_Nu(real turb_alpha, real cs_code, real omega_code);

real Dust_Scaleheight(real st, real hg_code, real turb_alpha, real cs_code, real omega_code);

real Volume_Rho_Midplane(real sigma_g_code, real hg_code);

real Pgas(real rho_g_mid, real cs_code);

real Pressure_Eta(real *sigma_g_code, real *cs, real hg_code, RadGeo g, int idx);

real Vel_Settle(real st, real hd_code, real cs_code, real omega_code, real turb_alpha);

real Vel_Visc(real *sigma_g_code, real *cs, real turb_alpha, RadGeo g, int idx);

real Vel_Gas_Rad(real v_visc, real p_eta, real vk, real A, real B);

real v_drift_max(real v_visc, real p_eta, real A, real B, real vk);

real Vel_Dust_Rad(real v_gas_r, real v_drift_max, real vk, real A, real B, real turb_alpha, real st);

real Vel_Dust_Azi_Offset(real v_drift_max, real st);

real Rel_Vel_Brownian(real s1_cgs, real s2_cgs, real cs_code);

real Rel_Vel_Settle(real v_settle1, real v_settle2);

real Rel_Vel_Rad(real st1, real st2, real v_gas_r, real v_drift_max);

real Rel_Vel_Azi(real st1, real st2, real v_drift_max);

real Rel_Vel_Turb(real omega_code, real sigma_g_code, real cs_code, real turb_alpha, real st1, real st2);

real Rel_Vel_Tot(real rel_vel_azi, real rel_vel_brownian, real rel_vel_rad, real rel_vel_turb, real rel_vel_settle);
