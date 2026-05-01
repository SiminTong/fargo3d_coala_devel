#include "coala_precision.h"
#include "math.h"
#include "coala_tension_module.h"
#include <stdbool.h>
#include "float.h"
#include <stdlib.h>
#include <stdio.h>

void TSPSI(unsigned int N, flt* COALA_X, flt* COALA_Y, flt* YP, flt* SIGMA, unsigned int IER, unsigned int SigErr){
	int IERR;
	IER = 0;

	if (N < 2){
		IER = -1;
	}else{
		YPC1(N,COALA_X,COALA_Y,YP,IERR);

		if (IERR != 0){
			IER = -4;
		}else{
			SIGS(N,COALA_X,COALA_Y,YP,SIGMA,IERR,SigErr);
		}
	}
}

void SIGS(unsigned int N,flt* COALA_X,flt* COALA_Y,flt* YP,flt* SIGMA, unsigned int IER,unsigned int SigErr){
	int I,ICNT,IP1,LUN,NIT,NM1;
	flt A,C1,C2,COSHM,COSHMM,D0,D1;
	flt D1D2, D1PD2, D2, DMAX, DSIG, DSM, DX;
    flt E, EMS, EMS2, F, F0, FMAX, FNEG, FP;
    flt FTOL, RTOL, S, S1, S2, SBIG, SCM;
    flt SGN, SIG, SIGIN, SINHM, SSINH, SSM;
    flt STOL, T, T0, T1, T2, TM, TP1;

    flt TOL, DSMAX;
    bool CONT,FLAG;

    FLAG = false;
    SBIG = ((flt)85);
    LUN = -1;
    NM1 = N - 2;
    if (NM1 < 0){
      DSMAX = ((flt)0);
      IER = -1;
    }else{
    	TOL = ((flt)0);
    	for (unsigned int I=0;I <= NM1; I++){
    		SIGMA[I] = ((flt)0);
    	}

		FTOL = fabs(TOL);
      	RTOL = ((flt)1);
      	while ( RTOL >= DBL_EPSILON  ){
      		RTOL = RTOL/((flt)2);
      	}
      	RTOL = RTOL*((flt)200);

		ICNT = 0;
      	DSM = ((flt)0);
      	for (unsigned int I=0;I <= NM1;I++){
      		IP1 = I + 1;
        	DX = COALA_X[IP1] - COALA_X[I];
        	if (DX <= ((flt)0)){
          		IER = -IP1;
          		exit(-1);
        	}
        	SIGIN = SIGMA[I];
        	if (SIGIN >= SBIG){
        		continue;
        	} 

        	S1 = YP[I];
	        S2 = YP[IP1];
	        S = (COALA_Y[IP1]-COALA_Y[I])/DX;
	        D1 = S - S1;
	        D2 = S2 - S;
	        D1D2 = D1*D2;

	        SIG = SBIG;
	        if ((D1D2 == ((flt)0)  &&  S1 != S2) || (S == ((flt)0) &&  S1*S2 > ((flt)0))){
	        	SIG = fmin(SIG,SBIG);
		        if (SIG > SIGIN){
		        	SIGMA[I] = SIG;
		            ICNT = ICNT + 1;
		            DSIG = SIG-SIGIN;
		            if (SIGIN > ((flt)0)){
		            	DSIG = DSIG/SIGIN;
		            } 
		            DSM = fmax(DSM,DSIG);
		        }
		        continue;
	        }

	        SIG = ((flt)0);
	        if (D1D2 >= ((flt)0)){
	        	if (D1D2 == ((flt)0)){
	        		SIG = fmin(SIG,SBIG);
            		if (SIG > SIGIN){
            			SIGMA[I] = SIG;
              			ICNT = ICNT + 1;
              			DSIG = SIG-SIGIN;
              			if (SIGIN > ((flt)0)){
              				DSIG = DSIG/SIGIN;
              			} 
              			DSM = fmax(DSM,DSIG);
            		}
              		continue;
	        	}

	        	T = fmax(D1/D2,D2/D1);
	        	if (T-((flt)2) <= DBL_EPSILON ){
	        		SIG = fmin(SIG,SBIG);
	        		if (SIG > SIGIN){
	        			SIGMA[I] = SIG;
		              	ICNT = ICNT + 1;
		              	DSIG = SIG-SIGIN;
		              	if (SIGIN > ((flt)0)){
		              		DSIG = DSIG/SIGIN;
		              	}
		              	DSM = fmax(DSM,DSIG);
	        		}
	        		continue;
	        	}
		            
		        TP1 = T + ((flt)1);
	            SIG = sqrt(((flt)10)*T-((flt)20));
	            NIT = 0;

	            while (NIT <= 10000){
	            	if (SIG <= ((flt)5e-1)){
	            		SNHCSH(SIG, &SINHM,&COSHM,&COSHMM);
	            		T1 = COSHM/SINHM;
              			FP = T1 + SIG*(SIG/SINHM - T1*T1 + ((flt)1));
	            	}else{
						EMS = exp(-SIG);
			            SSM = ((flt)1) - EMS*(EMS+SIG+SIG);
			            T1 = (((flt)1)-EMS)*(((flt)1)-EMS)/SSM;
			            FP = T1 + SIG*(((flt)2)*SIG*EMS/SSM - T1*T1 + ((flt)1));
	            	}

	            	F = SIG*T1 - TP1;
	            	NIT = NIT + 1;
	            	FLAG = false;
	            	if (FP <= ((flt)0)){
	            		FLAG = true;
	            		break;
	            	}
	            	DSIG = -F/FP;

	            	if((fabs(DSIG) <= RTOL*SIG) || (F >= ((flt)0) && F <= FTOL) || fabs(F) <= RTOL){
	            		FLAG = true;
	            		break;
	            	}
	            	SIG += DSIG;
	            }

	            if (NIT > 10000){
	            	SigErr=SigErr+1;
            		exit(-1);
	            }
	        }

	        if (FLAG){
	        	FLAG = false;
	        	SIG = fmin(SIG,SBIG);
	        	if (SIG > SIGIN){
	        		SIGMA[I] = SIG;
		            ICNT = ICNT + 1;
		            DSIG = SIG-SIGIN;
		            if (SIGIN > ((flt)0)){
		            	DSIG = DSIG/SIGIN;
		            }
		            DSM = fmax(DSM,DSIG);
	        	}
	        	continue;
	        }

			if (S1*S < ((flt)0) ||  S2*S < ((flt)0)){
				SIG = fmin(SIG,SBIG);
				if (SIG > SIGIN){
					SIGMA[I] = SIG;
		            ICNT = ICNT + 1;
		            DSIG = SIG-SIGIN;
		            if (SIGIN > ((flt)0)){
		            	DSIG = DSIG/SIGIN;
		            }
		            DSM = fmax(DSM,DSIG);
				}
				continue;
			}

			T0 = ((flt)3)*S - S1 - S2;
        	D0 = T0*T0 - S1*S2;

			if (D0 <= ((flt)0) || S*T0 >= ((flt)0)){
				SIG = fmin(SIG,SBIG);
				if (SIG > SIGIN){
					SIGMA[I] = SIG;
            		ICNT = ICNT + 1;
            		DSIG = SIG-SIGIN;
            		if (SIGIN > ((flt)0)){
            			DSIG = DSIG/SIGIN;
            		}
            		DSM = fmax(DSM,DSIG);
				}
				continue;
			}

			SGN = copysign(((flt)1),S);
        	SIG = SBIG;
        	FMAX = SGN*(SIG*S-S1-S2)/(SIG-((flt)2));
        	if (FMAX <= ((flt)0)){
        		SIG = fmin(SIG,SBIG);
        		if (SIG > SIGIN){
        			SIGMA[I] = SIG;
		            ICNT = ICNT + 1;
		            DSIG = SIG-SIGIN;
		            if (SIGIN > ((flt)0)){
		            	DSIG = DSIG/SIGIN;
		            }
		            DSM = fmax(DSM,DSIG);
        		}
        		continue;
        	}
        	STOL = RTOL*SIG;
	        F = FMAX;
	        F0 = SGN*D0/(((flt)3)*(D1-D2));
	        FNEG = F0;
	        DSIG = SIG;
	        DMAX = SIG;
	        D1PD2 = D1 + D2;
	        NIT = 0;

			while (!((fabs(DMAX) <= STOL) || (F >= ((flt)0) && F <= FTOL) || (fabs(F) <= RTOL))) {
				DSIG = -F*DSIG/(F-F0);
				if (fabs(DSIG) > fabs(DMAX) || DSIG*DMAX > ((flt)0)){
					DSIG = DMAX;
        			F0 = FNEG;
				}

				if (fabs(DSIG) < STOL/((flt)2)){
					DSIG = -copysign(STOL/((flt)2),DMAX);
				}

				SIG = SIG + DSIG;
          		F0 = F;
          		if (SIG <= ((flt)5e-1)){
					SNHCSH(SIG, &SINHM,&COSHM,&COSHMM);
					C1 = SIG*COSHM*D2 - SINHM*D1PD2;
		            C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2;
		            A = C2 - C1;
		            E = SIG*SINHM - COSHMM - COSHMM;
          		}else{
					EMS = exp(-SIG);
		            EMS2 = EMS + EMS;
		            TM = ((flt)1) - EMS;
		            SSINH = TM*(((flt)1)+EMS);
		            SSM = SSINH - SIG*EMS2;
		            SCM = TM*TM;
		            C1 = SIG*SCM*D2 - SSM*D1PD2;
		            C2 = SIG*SSINH*D2 - SCM*D1PD2;

		            F = FMAX;
            		CONT = true;
            		if (C1*(SIG*SCM*D1 - SSM*D1PD2) >= ((flt)0)){
            			CONT = false;
            		}
            		if (CONT){
            			A = EMS2*(SIG*TM*D2 + (TM-SIG)*D1PD2);
            		}
            		if(A*(C2+C1) < ((flt)0)){
            			CONT = false;
            		}
            		if (CONT){
            			E = SIG*SSINH - SCM - SCM;
            		}
          		}

          		if (CONT){
          			F = (SGN*(E*S2-C2) + sqrt(A*(C2+C1)))/E;
          		}

				NIT = NIT + 1;
				STOL = RTOL*SIG;
				if (fabs(DMAX) <= STOL || (F >= ((flt)0) && F <= FTOL) || fabs(F) <= RTOL ){
					exit(-1);
				}
				DMAX = DMAX + DSIG;
				if (F0*F > ((flt)0) && fabs(F) >= fabs(F0)){
					DSIG = DMAX;
            		F0 = FNEG;
				}

				if (F0*F <= ((flt)0)){
					T1 = DMAX;
		            T2 = FNEG;
		            DMAX = DSIG;
		            FNEG = F0;
		            if (fabs(DSIG) > fabs(T1)  && fabs(F) < fabs(T2)){
		            	DSIG = T1;
              			F0 = T2;
		            }
				}
			}

			SIG = fmin(SIG,SBIG);
			if (SIG > SIGIN){
				SIGMA[I] = SIG;
	            ICNT = ICNT + 1;
	            DSIG = SIG-SIGIN;
	            if (SIGIN > ((flt)0)){
	            	DSIG = DSIG/SIGIN;
	            }
	            DSM = fmax(DSM,DSIG);
			}
      	}
    }
}

void SNHCSH(flt COALA_X, flt* SINHM,flt* COSHM, flt* COSHMM){
	flt AX, EXPX, F, P, P1, P2, P3, P4, Q;
	flt Q1, Q2, Q3, Q4, XC, XS, XSD2, XSD4;

	P1 = -((flt)3.51754964808151394800e5);
	P2 = -((flt)1.15614435765005216044e4);
	P3 = -((flt)1.63725857525983828727e2);
	P4 = -((flt)7.89474443963537015605e-1);
	Q1 = -((flt)2.11052978884890840399e6);
	Q2 = ((flt)3.61578279834431989373e4);
	Q3 = -((flt)2.77711081420602794433e2);
	Q4 = ((flt)1);
	AX = fabs(COALA_X);
    XS = AX*AX;

    if (AX <= ((flt)5e-1)){
		XC = COALA_X*XS;
		P = ((P4*XS+P3)*XS+P2)*XS+P1;
		Q = ((Q4*XS+Q3)*XS+Q2)*XS+Q1;
		*SINHM = XC*(P/Q);
		XSD4 = ((flt)25e-2)*XS;
		XSD2 = XSD4 + XSD4;
		P = ((P4*XSD4+P3)*XSD4+P2)*XSD4+P1;
		Q = ((Q4*XSD4+Q3)*XSD4+Q2)*XSD4+Q1;
		F = XSD4*(P/Q);
		*COSHMM = XSD2*F*(F+((flt)2));
		*COSHM = *COSHMM + XSD2;
    }else{
    	EXPX = exp(AX);
      	*SINHM = -(((((flt)1)/EXPX+AX)+AX)-EXPX)/((flt)2);
      	if (COALA_X < ((flt)0)){
      		*SINHM = - *SINHM;
      	}
      	*COSHM = ((((flt)1)/EXPX-((flt)2))+EXPX)/((flt)2);
        *COSHMM = *COSHM - XS/((flt)2);
    }
}

void YPC1(unsigned int N, flt* COALA_X, flt* COALA_Y, flt* YP, unsigned int IER){
	int I, NM1;
	flt ASI, ASIM1, DX2, DXI, DXIM1, S2, SGN,SI, SIM1, T;

	NM1 = N - 2;
    I = 1;
    DXI = COALA_X[1] - COALA_X[0];
    if (DXI <= ((flt)0)){
      	IER = I + 1;
      	exit(-1);
    }

    SI = (COALA_Y[1]-COALA_Y[0])/DXI;
    if (NM1 == 0){
    	YP[0] = SI;
      	YP[1] = SI;
      	IER = 0;
      	exit(-1);
    }

	I = 2;
    DX2 = COALA_X[2] - COALA_X[1];
    if (DX2 <= ((flt)0)){
      	IER = I + 1;
      	exit(-1);
    }
    S2 = (COALA_Y[2]-COALA_Y[1])/DX2;
    T = SI + DXI*(SI-S2)/(DXI+DX2);
    if (SI >= ((flt)0)){
    	YP[0] = fmin(fmax(((flt)0),T), ((flt)3)*SI);
    }else{
    	YP[0] = fmax(fmin(((flt)0),T), ((flt)3)*SI);
    }

	for (unsigned int I=1;I <= NM1; I++){
		DXIM1 = DXI;
      	DXI = COALA_X[I+1] - COALA_X[I];
      	if (DXI <= ((flt)0)){
        	IER = I + 1;
        	break;
      	}
      	SIM1 = SI;
		SI = (COALA_Y[I+1]-COALA_Y[I])/DXI;
		T = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI);
		ASIM1 = fabs(SIM1);
		ASI = fabs(SI);
		SGN = copysign(((flt)1),SI);

		if (ASIM1 > ASI){
			SGN = copysign(((flt)1),SIM1);
		}

		if (SGN > ((flt)0)){
			YP[I] = fmin(fmax(((flt)0),T),((flt)3)*fmin(ASIM1,ASI));
		}else{
			YP[I] = fmax(fmin(((flt)0),T),-((flt)3)*fmin(ASIM1,ASI));
		}
	}

	T = SI + DXI*(SI-SIM1)/(DXIM1+DXI);
	if (SI >= ((flt)0)){
		YP[N-1] = fmin(fmax(((flt)0),T), ((flt)3)*SI);
	}else{
		YP[N-1] = fmax(fmin(((flt)0),T), ((flt)3)*SI);
	}
	IER = 0;
}

flt HVAL(flt T, unsigned int N, flt* COALA_X, flt* COALA_Y, flt* YP, flt* SIGMA, unsigned int IER){
	int I,IP1;
	flt B1, B2, CM, CM2, CMM, D1, D2, DUMMY;
	flt DX, E, E1, E2, EMS, S, S1, SB1, SB2;
	flt SBIG, SIG, SM, SM2, TM, TP, TS, U, Y1;
	flt resHVAL;

	SBIG = ((flt)85);
	if (N <= 2){
		resHVAL = ((flt)0);
		IER = -1;
	}else{
		if (T < COALA_X[0]){
			I = 1;
        	IER = 1;
		}else if (T > COALA_X[N-1]){
			I = N-2;
			IER = 1;
		}else{
			I = INTRVL(T,N,COALA_X);
			IER = 0;
		}

		IP1 = I + 1;
		DX = COALA_X[IP1] - COALA_X[I];
		if (DX <= ((flt)0)){
        	resHVAL = ((flt)0);
        	IER = -2;
		}else{
			U = T - COALA_X[I];
	        B2 = U/DX;
	        B1 = ((flt)1) - B2;
	        Y1 = COALA_Y[I];
	        S1 = YP[I];
	        S = (COALA_Y[IP1]-Y1)/DX;
	        D1 = S - S1;
	        D2 = YP[IP1] - S;
	        SIG = fabs(SIGMA[I]);
	        if (SIG < ((flt)1e-9)){
	        	resHVAL = Y1 + U*(S1 + B2*(D1 + B1*(D1-D2)));
	        }else if (SIG <= ((flt)0.5)){
				SB2 = SIG*B2;
				SNHCSH(SIG, &SM,&CM,&CMM);
				SNHCSH(SB2, &SM2,&CM2,&DUMMY);
				E = SIG*SM - CMM - CMM;
				resHVAL = Y1 + S1*U + DX*((CM*SM2-SM*CM2)*(D1+D2) + SIG*(CM*CM2-(SM+SIG)*SM2)*D1) / (SIG*E);
	        }else{
				SB1 = SIG*B1;
          		SB2 = SIG - SB1;
          		if (-SB1 > SBIG  ||  -SB2 > SBIG){
          			resHVAL = Y1 + S*U;
          		}else{
          			E1 = exp(-SB1);
		            E2 = exp(-SB2);
		            EMS = E1*E2;
		            TM = ((flt)1) - EMS;
		            TS = TM*TM;
		            TP = ((flt)1) + EMS;
		            E = TM*(SIG*TP - TM - TM);
		            resHVAL = Y1 + S*U + DX*(TM*(TP-E1-E2)*(D1+D2) +  
		                   SIG* ((E2+EMS*(E1-((flt)2))-B1*TS)*D1 +     
		                   (E1+EMS*(E2-((flt)2))-B2*TS)*D2)) / (SIG*E);
          		}
	        }
		}
	}
	return resHVAL;
}

flt HPVAL(flt T, unsigned int N, flt* COALA_X, flt* COALA_Y, flt* YP, flt* SIGMA, unsigned int IER){
	int I,IP1;
	flt B1, B2,   D1, D2;
	flt DX, E, E1, E2, EMS, S, S1, SB1, SB2;
	flt SBIG, SIG, SINH2, SM, SM2, TM;
	flt CM, CM2,CMM,DUMMY;
	flt resHPVAL;

	SBIG = ((flt)85);
	if (N <= 2){
		resHPVAL = ((flt)0);
		IER = -1;
	}else{
		if (T < COALA_X[0]){
			I = 1;
        	IER = 1;
		}else if (T > COALA_X[N-1]){
			I = N-2;
			IER = 1;
		}else{
			I = INTRVL (T,N,COALA_X);
			IER = 0;
		}

		IP1 = I + 1;
		DX = COALA_X[IP1] - COALA_X[I];
		if (DX <= ((flt)0)){
        	resHPVAL = ((flt)0);
        	IER = -2;
		}else{
			B1 = (COALA_X[IP1] - T)/DX;
	        B2 = ((flt)1) - B1;
	        S1 = YP[I];
	        S = (COALA_Y[IP1]-COALA_Y[I])/DX;
	        D1 = S - S1;
	        D2 = YP[IP1] - S;
	        SIG = fabs(SIGMA[I]);
	        if (SIG < ((flt)1e-9)){
	        	resHPVAL = S1 + B2*(D1 + D2 - ((flt)3)*B1*(D2-D1));
	        }else if (SIG <= ((flt)5e-1)){
				SB2 = SIG*B2;
				SNHCSH(SIG, &SM,&CM,&CMM);
				SNHCSH(SB2, &SM2,&CM2,&DUMMY);
				SINH2 = SM2 + SB2;
          		E = SIG*SM - CMM - CMM;
          		resHPVAL = S1 + ((CM*CM2-SM*SINH2)*(D1+D2) + SIG*(CM*SINH2-(SM+SIG)*CM2)*D1)/E;
	        }else{
				SB1 = SIG*B1;
          		SB2 = SIG - SB1;
          		if (-SB1 > SBIG  ||  -SB2 > SBIG){
          			resHPVAL = S;
          		}else{
          			E1 = exp(-SB1);
		            E2 = exp(-SB2);
		            EMS = E1*E2;
		            TM = ((flt)1) - EMS;
		            E = TM*(SIG*(((flt)1)+EMS) - TM - TM);
		            resHPVAL = S + (TM*((E2-E1)*(D1+D2) + TM*(D1-D2)) 
		                + SIG*((E1*EMS-E2)*D1 + (E1-E2*EMS)*D2))/E;
          		}
	        }
		}
	}
	return resHPVAL;
}

int INTRVL (flt T,unsigned int N,flt* COALA_X){
	int IH,IL,K;
	flt TT;
	bool FLAG;

	IL = 1;
	TT = T;
    FLAG = true;
    if (IL >= 1  &&  IL < N){
    	if (COALA_X[IL] <= TT && TT < COALA_X[IL+1]){
    		FLAG = false;
    	}
    }

    if (FLAG){
    	IL = 1;
        IH = N;
        while (IH > IL+1){
        	K = (IL+IH)/2;
        	if (TT < COALA_X[K]){
        		IH = K;
        	}else{
        		IL = K;
        	}
        }
    }
    return IL;
}