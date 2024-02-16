      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1 LAYER,KSPT)
      
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
      PARAMETER(ZERO=0.0E0)
      

      STATEV(1)=0.999999E0
      STATEV(2)=0.000001E0
     
      DO i = 3, NSTATV-1
       STATEV(i)=ZERO
      END DO
      RETURN
      END
      
      
      
C-------------------------------------------------------------------------
C-------Main UMAT SUBROUTINE for plain strain conditions----------------------------------------------
C234567
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)
      !implicit none
      INCLUDE 'ABA_PARAM.INC'
      
      CHARACTER*80 CMNAME
      
      DIMENSION STRESS(NTENS), STATEV(NSTATV), 
     1 DDSDDE(NTENS, NTENS),
     2 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     3 TIME(2), PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),
     4 DFGRD0(3, 3), DFGRD1(3, 3)
     
C LOCAL ARRAYS
C ----------------------------------------------------------------
C C_e1=Elastic stiffness matrix for phase 0 (original phase)
C C_e2=Elastic stiffness matrix for phase 1 (new/transformed phase)
C C_1_in=inverse of Elastic stiffness matrix for phase 0
C C_2_in=inverse of Elastic stiffness matrix for phase 1
C C_eff=Effective/total Elastic stiffness for entire body
C Cef_in=Inverse of Effective/total Elastic stiffness for entire body
C de_dev=deviatoric part of DSTRAN
C EG= Shear Modulus G/mu for both the phases
C STR=STRESS in a phase (sigma_ij)
C DEELAS - CHANGE IN ELASTIC STRAINS
C ESTRN0- SUM OF STRAIN AND PLASTIC STRAIN FOR PHASE0
C ESTRN1- SUM OF STRAIN AND PLASTIC STRAIN FOR PHASE1
C DPLAS-CHANGE IN PLASTIC STRAINS IN INDIVIDUAL PHASES
C DEPLAS-CHANGE IN TOTAL PLASTIC STRAINS
C EPLS0-PLASTIC STRAINS IN PHASE0
C EPLS1-PLASTIC STRAINS IN PHASE1
C EPLAS-TOTAL PLASTIC STRAINS
C S_dev=Deviatoric part of STR
C Sdev_t=Trial deviatoric stress 
C n=Flow direction in plasticity
C n_n=Dyadic product of two n vectors
C LOCAL ARRAYS
C ----------------------------------------------------------------
C EELAS - ELASTIC STRAINS
C 
C FLOW - DIRECTION OF PLASTIC FLOW
C ----------------------------------------------------------------
C ---------------Declaration of size of arrays and type of variables-------------------------------------------------
      DIMENSION EPLAS(NTENS),C_e2(NTENS,NTENS),dd1(NTENS,NTENS),
     1ESTRN0(NTENS),C_e1(NTENS,NTENS),dd2(NTENS,NTENS),DTHRM(NTENS),
     2 DEELAS(NTENS),DPLAS(NTENS),dd21(NTENS),dd22(NTENS),ETHRM1(NTENS),
     3EPLS0(NTENS),EPLS1(NTENS),ESTRN1(NTENS),DRDE(NTENS),ETHRM0(NTENS),
     4Cef_in(NTENS,NTENS),de_dev(NTENS),EG(2),STR(NTENS),DTHRM1(NTENS),
     5S_dev(NTENS),Sdev_t(NTENS),C_eff(NTENS,NTENS),D_n(NTENS),Cef_i(NTENS,NTENS),
     6DEPLAS(NTENS),Dn_n(NTENS,NTENS),Ce(NTENS,NTENS),STRN1(NTENS),
     7STRN0(NTENS),aLamd1(2),dlam(2),I_d(NTENS),dX_tr1(2),Cef_n(NTENS),
     8aL_tr(2),C_2_in(NTENS,NTENS),C_1_in(NTENS,NTENS),dX_tr(2),
     9dlamda(2),X_tr(2),Ce2(NTENS,NTENS),Ce1(NTENS,NTENS),X_tr1(2)
      DIMENSION dd3(NTENS,NTENS),DRDE2(NTENS),DMDE(NTENS),dsdm(NTENS),
     1ds_dl1(NTENS),de0dl(NTENS),de1dl(NTENS),II_d(NTENS,NTENS),
     2I3(NTENS,NTENS),DSTRN0(NTENS),STRN10(NTENS),DEV_I(NTENS,NTENS),
     3STRN00(NTENS),ESTN0(NTENS),ESTN1(NTENS),STRES1(NTENS),STR1(NTENS),
     4STRAN1(NTENS),dlm(2),STRAN2(NTENS),STRES2(NTENS),STRN12(NTENS),
     5STRN02(NTENS),ESTR0(NTENS),ESTR1(NTENS),STRS2(NTENS),STRS22(NTENS)
     6,STRS1(NTENS),STRS11(NTENS),eta1(NTENS),eta0(NTENS),eta_f(NTENS),
     7STR0(NTENS),Sdev0(NTENS),STRN_0(NTENS),Sdev1(NTENS),Sdev2(NTENS),
     8STN_0(NTENS),SN_A(NTENS),STN_1(NTENS),SN_B(NTENS),STN_E(NTENS),
     9OLDPL(6),OLDS(6),IkI(NTENS,NTENS),Sklij(NTENS,NTENS),DTHRM0(NTENS)
      
      DIMENSION ETHRM(NTENS),OLDTH(NTENS)
      DIMENSION EELAS(6),FLOW(6), HARD(3),DPSDDE(NTENS,NTENS)
      !double precision dX_tr1(2),dX_tr(2)
      LOGICAL,save:: firstcall=.true.
      PARAMETER(ZERO=0.0E0, ONE=1.0E0, TWO=2.0E0, 
     1THREE=3.0E0,TOLER=1.0E-7,TOLER1=1.0E-6,TOLER2=1.0E-3,
     2SIX=6.E0,ENUMAX=.4999E0, NEWTON=10, TOLER3=1.0E-6)
     
     !!for debug
      !if (firstcall==.true.) then
        ! firstcall=.false.
        ! read(*,*) dummyvar
         
      !end if
      !dummyvar=1
      !debug ends
C ----------------------------------------------------------------
C Material Properties
C Total number of material props=14
C PROPS(1) - E0, Elastic modulus
C PROPS(2) - E1, Elastic modulus
C PROPS(3) - NU0, poisson ratio
C PROPS(4) - NU1, poisson ratio
C PROPS(5) - f_c, energy barrier of phase transformation evolution 
C PROPS(6) - v_tr, a parametre of phase transformation evolution
C PROPS(7) - n_tr, a parametre of phase transformation evolution
C PROPS(8) - s_1, a parametre of phase transformation evolution
C PROPS(9) - s_2, a parametre of phase transformation evolution
C PROPS(10) - H_j, a parametre of plastic evolution
C PROPS(11) - n_pl, a parametre of plastic evolution
C PROPS(12) - v_pl,a parametre of plastic evolution
C PROPS(13)- Y_0,a parametre of plastic evolution yield stress of granulte
C PROPS(14)-Y_1, yield stress of eclogite
C PROPS(15)-inelastic heat fraction (n), for heating due to plastic dissipation 
C PROPS(16)-T_initial, initial temperature 973K
C PROPS(17)-alpha, coeff of thermal expansion olivine
C PROPS(18)-alpha2, coeff of thermal expansion spinel
C PROPS(19)-chem10, constant chemical force for olivine
C PROPS(20)-chem11, constant chemical force for spinel
C PROPS(21)-ch_t10, variable part of chemical force with temperature for olivine
C PROPS(22)-ch_t11, variable part of chemical force with temperature for spinel
C !PROPS(18)-Q_0, a parametre of plastic evolution (currently not used)
C !PROPS(19)-phi_ch, constant chemical energy (currently not used)
C !PROPS(20)-dphi_ch, changing part of chemical energy (currently not used)
C  PROPS(23)-PROPS(28):eta1
C--------------------------------------------------
C SOLUTION DEPENDENT VARIABLES
C 
C No. of SDVs=33
C STATEV(1)= aLmd_g, Volume fraction of granulite 10(phase 0)
C STATEV(2)=aLamda, Volume fraction of eclogite 11(phase 1)
C STATEV(3)=q_10, internal hardening variable for phase 0
C STATEV(4)=q_11, internal hardening variable for phase 0
C STATEV(5)=EPLS0(1), FIRST ELEMENT OF PLASTIC STRAIN IN PHASE 10
C STATEV(6)=EPLS1(1), FIRST ELEMENT OF PLASTIC STRAIN IN PHASE 11
C STATEV(7)=EPLAS(1), FIRST ELEMENT OF TOTAL PLASTIC STRAIN
C STATEV(8)=Mu_10, Viscoplasticity modulus for phase 10
C STATEV(9)=Mu_11, Viscoplasticity modulus for phase 11
C STATEV(10)=EPLS0(2), SECOND ELEMENT OF PLASTIC STRAIN IN PHASE 10
C STATEV(11)=EPLS0(3), THIRD ELEMENT OF PLASTIC STRAIN IN PHASE 10
C STATEV(12)=EPLS1(2), SECOND ELEMENT OF PLASTIC STRAIN IN PHASE 11
C STATEV(13)=EPLS1(3), THIRD ELEMENT OF PLASTIC STRAIN IN PHASE 11
C STATEV(14)=EPLAS(2), SECOND ELEMENT OF TOTAL PLASTIC STRAIN
C STATEV(15)=EPLAS(3), THIRD ELEMENT OF TOTAL PLASTIC STRAIN
C STATEV(16)=EPLAS(4), FOURTH ELEMENT OF TOTAL PLASTIC STRAIN
C STATEV(17)=EPLS0(4), FOURTH ELEMENT OF PLASTIC STRAIN IN PHASE 10
C STATEV(18)=EPLS1(4), FOURTH ELEMENT OF PLASTIC STRAIN IN PHASE 11
C STATEV(19)=gam3,lagrangian multipliers for lambda
C STATEV(20)=gam4,lagrangian multipliers for lambda
C STATEV(21)=EPLS0(5)
C STATEV(22)=EPLS0(6)     
C STATEV(23)=EPLS1(5)
C STATEV(24)=EPLS1(6)
C STATEV(25)=EPLAS(5)
C STATEV(26)=EPLAS(6)
C STATEV(27)=
C STATEV(34)=ETHRM(1)
C STATEV(35)=ETHRM(2)    
C STATEV(36)=ETHRM(3)
C STATEV(37)=ETHRM(4)
C STATEV(38)=ETHRM(5)
C STATEV(39)=ETHRM(6)

      
C ----------------------------------------------------------------
C ELASTIC PROPERTIES AND Ce MATRICES FOR PHASE 10 AND 11
C K1=1 stands for PHASE 10 (granulite) and K1=2 stands for PHASE11 (eclogite), C_e1 for granulite and C_e2 for eclogite
      
      DO K1=1,2
          EMOD=PROPS(K1)            ! Elastic modulus
          ENU=PROPS(K1+2)             ! Poisson ratio
          EBULK3=EMOD/(ONE-TWO*ENU)  ! 3 times K (bulk modulus)
          EG2=EMOD/(ONE+ENU)         ! 2 times G (shear modulus)
          EG(K1)=EG2/TWO                 ! G (shear modulus)
          ELAM=(EBULK3-EG2)/THREE    ! Lame parameter aLamda
          ELM=ELAM+TWO*EG(K1)
          if (K1 .EQ. 1) then
              CALL KC_e2(ELM,ELAM,EG(K1),NDI,NTENS,C_e1)  !CALLING SUBROUTINE KC_e to form stiffness matrix
          else 
              CALL KC_e1(ELM,ELAM,EG(K1),NDI,NTENS,C_e2)
          end if 
      END DO
      Ce1=C_e1
      Ce2=C_e2
     
      EQPLAS=STATEV(40)
      
      CALL inverse(Ce1,C_1_in,NTENS) ! Taking inverse of stiffness matrix
      CALL inverse(Ce2,C_2_in,NTENS)
      if (TIME(2).GT.5.0E0) then
      PROPS(5)=125.0E0
      end if

C-------------------------------------------------------
C ALGORITHM FOR PHASE TRANSFOMRATION
C 0 for granulite and 1 for eclogite
C Define trial state and elastic predictor
      aLamd1=(/STATEV(1),STATEV(2)/)!initial lambda values
      r_tr=100000.0E0 !residual initialization
      dX_tr(2)=.000E0!Unknown vector initialization
      dX_tr(1)=-dX_tr(2)
      gam1=ZERO !lagrangian multiplier lower bound
      gam2=ZERO !lagrangian multiplier upper bound
      s_tr=0.00000001E0!direction of search
      theta=TEMP
     
      !spinel to olivine
     
      chem10=PROPS(19)
      chem11=PROPS(20)
      ch_t10=PROPS(21)
      ch_t11=PROPS(22)
      
      eta0=(/0.0E0,0.0E0,0.0E0,0.0E0,0.0E0,0.0E0/)
      
       eta1=(/PROPS(23),PROPS(24),PROPS(25),PROPS(26),PROPS(27),PROPS(28)/)!U01
      
      
C define fixed scalars
      beta= .10E0; alpha0=1.0E0
C Line search   
C constraints on aLamda and dlamda
      AL1min=0.000001E0;AL1max=0.999999E0
      dL11mn=ZERO;dL11mx=0.999998E0
      q_10=STATEV(3) ! needs to come from abaqus as initial value of statev
      q_11=STATEV(4) ! needs to come from abaqus as initial value of statev
      CALL KPL_EN(STATEV(3),PROPS(10),PSI_10) 
      CALL KPL_EN(STATEV(4),PROPS(10),PSI_11)
      EPLS0(1)=STATEV(5)
      EPLS0(2)=STATEV(10)
      EPLS0(3)=STATEV(11)
      EPLS0(4)=STATEV(17)
      EPLS0(5)=STATEV(21)
      EPLS0(6)=STATEV(22)
      EPLS1(1)=STATEV(6)
      EPLS1(2)=STATEV(12)
      EPLS1(3)=STATEV(13)
      EPLS1(4)=STATEV(18)
      EPLS1(5)=STATEV(23)
      EPLS1(6)=STATEV(24)
      EPLAS(1)=STATEV(7)
      EPLAS(2)=STATEV(14)
      EPLAS(3)=STATEV(15)
      EPLAS(4)=STATEV(16)
      EPLAS(5)=STATEV(25)
      EPLAS(6)=STATEV(26)
      DEPLAS(1)=STATEV(28)
      DEPLAS(2)=STATEV(29)
      DEPLAS(3)=STATEV(30)
      DEPLAS(4)=STATEV(31)
      DEPLAS(5)=STATEV(32)
      DEPLAS(6)=STATEV(33)
      ETHRM(1)=STATEV(34)
      ETHRM(2)=STATEV(35)
      ETHRM(3)=STATEV(36)
      ETHRM(4)=STATEV(37)
      ETHRM(5)=STATEV(38)
      ETHRM(6)=STATEV(39)
     
      DO K1=1, NTENS
          OLDS(K1)=STRESS(K1)
          OLDPL(K1)=EPLAS(K1)
          OLDTH(K1)=ETHRM(K1)
      END DO
      ch_f0=chem10+STATEV(2)*ch_t10*theta
      ch_f1=chem11+STATEV(2)*ch_t11*theta
      aL_tr=aLamd1+dX_tr
      aL_tr(2)=projec(aL_tr(2),AL1min,AL1max)
      aL_tr(1)=1-aL_tr(2)
      
      Cef_in=aL_tr(2)*(C_2_in)+aL_tr(1)*(C_1_in)! inverse of effective C matrix with updated lamda
      call inverse(Cef_in,C_eff,NTENS)! effective C matrix
C       CALCULATE THERMAL EXPANSION
      DO K1=1,NDI
          ETHRM0(K1)=PROPS(17)*(TEMP-PROPS(16))
          ETHRM1(K1)=PROPS(18)*(TEMP-PROPS(16))
          DTHRM0(K1)=PROPS(17)*DTEMP
          DTHRM1(K1)=PROPS(18)*DTEMP
      END DO
      DO K1=NDI+1,NTENS
          ETHRM0(K1)=ZERO
          ETHRM1(K1)=ZERO
          DTHRM0(K1)=ZERO
          DTHRM1(K1)=ZERO
      END DO
      ETHRM=aL_tr(1)*ETHRM0+aL_tr(2)*ETHRM1
      DTHRM=ETHRM-OLDTH
  
C strain  and stress  at n+1 step
      STRAN2=STRAN+DSTRAN
      
       eta_f=aL_tr(2)*eta1
       STRES2=STRESS+MATMUL(C_eff,(DSTRAN-DTHRM)) 
       gam3=STATEV(19)
       gam4=STATEV(20)
       STR=aL_tr(2)*STRES2
       STRN12=MATMUL(C_2_in,STR)
       STRN12(4)=0.5*STRN12(4)
       STRN12(5)=0.5*STRN12(5)
       STRN12(6)=0.5*STRN12(6)

       STRN12=STRN12+EPLS1+eta1+ETHRM1
       STRAN2(4)=0.5*STRAN2(4)
       STRAN2(5)=0.5*STRAN2(5)
       STRAN2(6)=0.5*STRAN2(6)
       STRN02=(STRAN2-(aL_tr(2)*STRN12))/aL_tr(1)
       ESTR0=STRN02+EPLS0 ! sum of total strain array and plastic strain array for one phase
       ESTR1=STRN12+EPLS1+eta1

     
          aL_tr=(/aLamd1(1)+dX_tr(1),aLamd1(2)+dX_tr(2)/)
          aL_tr(2)=projec(aL_tr(2),AL1min,AL1max)
          aL_tr(1)=1-aL_tr(2)
C Calculate Thermodynamic force with obtained lambdsa trail 
          Cef_in=aL_tr(2)*(C_2_in)+aL_tr(1)*(C_1_in)! inverse of effective C matrix with updated lamda
          call inverse(Cef_in,C_eff,NTENS)! effective C matrix
      
          STRN1=MATMUL(C_2_in,STRES2)
          STRN1(4)=0.5*STRN1(4)
            STRN1(5)=0.5*STRN1(5)
            STRN1(6)=0.5*STRN1(6)
          STRN1=STRN1+EPLS1+eta1+ETHRM1
      
          STRN0=(STRAN2-(aL_tr(2)*STRN1))/aL_tr(1)
          ESTRN0=STRN0+EPLS0+ETHRM0 ! sum of total strain array and plastic strain array
          ESTRN1=STRN1+EPLS1+eta1+ETHRM1
          BLAM10=(-ONE*(0.50E0*(dot_product(STRES2,(MATMUL(C_1_in,STRES2)))
     1)+(dot_product(STRES2,ESTRN0))-(PSI_10)-(ch_f0)))
          BLAM11=(-ONE*(0.50E0*(dot_product(STRES2,(MATMUL(C_2_in,STRES2)))
     1)+(dot_product(STRES2,ESTRN1))-(PSI_11)-(ch_f1))) !! calculation of thermodynamic forces
          BL_fo=BLAM10-BLAM11-gam3+gam4
         
      aH_p=-(STRESS(1)+STRESS(2)+STRESS(3))/3.0E0
      aLine=aH_p-(3.571E0*theta)+3786.0E0
      
      if (aLine.GT.ZERO) then
          
          phi=(BL_fo-PROPS(5))
      else 
         phi=ZERO
          
      end if

      statev(47)=phi
     
      
C Check for yield
      if (phi .LE.ZERO) then 
          dlamda=(/ZERO,ZERO/)
          aLamd1=aLamd1+dlamda
          DEELAS=DSTRAN-DTHRM
          DDSDDE=C_eff
          STATEV(27)=ZERO
      else

          aLamd1=(/STATEV(1),STATEV(2)/)
          Rs1_tr=10000.0E0
          X_tr(1)=STATEV(1)
          X_tr(2)=STATEV(2)
          s_tr1=0.000001E0
          
C define fixed scalars
          beta1= 0.10E0; alph0=1.0E0
C Line search   
C constraints on Lambda
          
         
          r_0tr=100000.0E0 !residual
          ncont1=0
          r_nom1=1000000.0E0
          Rs0_n=sqrt(r_0tr*r_0tr)
          Rs0_n0=sqrt(r_0tr*r_0tr)
          R_old=Rs1_tr
          PJc_in=1.0E0
          CALL KMAC(phi,phi_m)
          CALL Ksatur(STATEV(2),STATEV(1),PROPS(8),
     1                    PROPS(9),H_str1) 
          dX_tr(2)=DTIME*((1/PROPS(6))*(H_str1**PROPS(7))*
     1                   ((phi_m)**PROPS(7)))
          dX_tr(1)=-(dX_tr(2))
          DO WHILE ((r_nom1) .GT. TOLER2 .AND. 
     1    ncont1.LT.40)
              if (ncont1.GT.38) then
                  print*,"lamda NR did not converge"
                  ncont1=ncont1+1
                  
              else
                  mp1=0 
                  
                  Rs0_n=sqrt(Rl_tr1*Rl_tr1)
                  Rs1_n=sqrt(Rs1_tr*Rs1_tr)
                  r_nom1=100000.0E0
                  DO WHILE ((r_nom1-Rs1_n) .GT. TOLER2 .AND. mp1.LT.40)
                      if (mp1.GT.38) then
                          print*,"lambda LS did not converge"
                          mp1=mp1+1
                      else
                          a_str1=(beta1**mp1)*alph0
                          
                          dX_tr1(2)=(dX_tr(2)+a_str1*s_tr1)
                      dX_tr1(2)=projec(dX_tr1(2),dL11mn,dL11mx)
                      dX_tr1(1)=-(dX_tr1(2))
                      X_tr1=(/X_tr(1)+dX_tr1(1),X_tr(2)+dX_tr1(2)/) 


Calculate Residual

                          CALL Ksatur(X_tr1(2),X_tr1(1),PROPS(8),
     1                    PROPS(9),H_str) 
                          q_10=STATEV(3) 
                          q_11=STATEV(4) 
                          CALL KPL_EN(STATEV(3),PROPS(10),PSI_10) 
                          CALL KPL_EN(STATEV(4),PROPS(10),PSI_11) 
                          
                          dlambd=dX_tr1(2)
                          !print*,X_tr,X_tr1
                          dlam=(/-dlambd,dlambd/)
                          
     
      
                          CALL KResid(BL_fo,PROPS(5),PROPS(6),
     1                    DTIME,PROPS(7),H_str,dlam(2),RsdL1)
                          Rl_tr1=RsdL1
                          
         
C Working set for lambda
                          I_al1=0;I_au1=0;I_p1=0;I_h1=0;I_b1=0
                          
                          if (X_tr1(2).EQ.0.000001E0 .and.dX_tr1(2).LE.
     1                         ZERO.and.Rl_tr1.
     2                    LT.ZERO) then
                                  
                              I_al1=1
                          else if (X_tr1(2).EQ.0.999999E0 .and.
     1                     dX_tr1(2).GE.ZERO.and.Rl_tr1.
     2                    GT.ZERO) then
                                  
                              I_au1=1
                          else 
                                  
                              I_p1=1
                          end if
                          
                  
                          if (I_p1==0) then
                              
                              I_h1=1
                          else
                              I_b1=1
                          end if
                
C Calculate projectetd residual 
                          
                          if (I_p1.EQ.0) then
                              P1_rtr=ZERO
                          else
                              P1_rtr=Rl_tr1
                          end if
                          P1rnom=sqrt(P1_rtr*P1_rtr)
                          r_nom1=P1rnom
                          mp1=mp1+1
                          
                          
                      end if
                  END DO
                 
                  mp1=mp1-1
                 
                  !!C Calculate  jacobian 
    
C Breydon's approximation
                  
                  
                  del_x=dX_tr1(2)-dX_tr(2)
                  d_f=RsdL1-R_old
                 
                  
     
                  df_n=sqrt(d_f*d_f)
                  AJ1=PJc_in + ((del_x-(PJc_in*d_f))/(df_n*df_n))*
     1            d_f
                 

  
                 
C Calculate projected Jacobian 
                      
                      
                  if (I_p1.EQ.0) then      
                      PJ1 =ONE            
                  else
                      PJ1=AJ1
                  end if
                 
                  
                  s_tr1=-PJ1*Rl_tr1
                  ncont1=ncont1+1
                  dX_tr=dX_tr1
                  Rs1_tr=Rl_tr1
                  R_old=RsdL1
                  
                  PJc_in=AJ1
                  
              end if
      end DO
         
      
          
                  
                  
          
          if (r_nom1.LE.TOLER2) then
              STATEV(1)=X_tr1(1)
              STATEV(2)=X_tr1(2)
              STATEV(27)=dlam(2)
          end if
          
          if (I_al1==0) then
              gam3=ZERO

          else 
              gam4=-Rl_tr1
           
          end if  
      
          if (I_au1==0) then
              gam4=ZERO
          else 
              gam3=Rl_tr1
             
          end if 
          STATEV(19)=gam3
          STATEV(20)=gam4
          aLamda=STATEV(2)
          aLmd_g=STATEV(1)
          Cef_in=aLamda*(C_2_in)+aLmd_g*(C_1_in)
          call inverse(Cef_in,C_eff,NTENS)
         
          
          
         
           
      
          dd1=C_eff
          
          DDSDDE=dd1
       
      end if 
      
      
      aLamda=STATEV(2)
      aLmd_g=STATEV(1)
      Cef_in=aLamda*(C_2_in)+aLmd_g*(C_1_in)
      Cef_i=Cef_in
      call inverse(Cef_in,C_eff,NTENS)
      ETHRM=aLmd_g*ETHRM0+aLamda*ETHRM1
          
          STATEV(34)=ETHRM(1)
          STATEV(35)=ETHRM(2)    
          STATEV(36)=ETHRM(3)
          STATEV(37)=ETHRM(4)
          STATEV(38)=ETHRM(5)
          STATEV(39)=ETHRM(6)
          DTHRM=ETHRM-OLDTH
           DEELAS=DSTRAN-(STATEV(27)*eta1)-DTHRM
      
      STRAN1=STRAN+DSTRAN
      
      eta_f=aLamda*eta1
C------------------------------------------------------------
C ALGORITHM FOR PLASTICITY
   !defining Id
      DO K1=1,NDI
          I_d(K1)=1
      END DO
      DO K1=NDI+1,NTENS
          I_d(K1)=0
      END DO
      !defining Sklij
      DO K1=1,NTENS
          DO K2=1,NTENS
              Sklij(K1,K2)=0.0
          END DO
      END DO
      DO K1=1,NDI
          Sklij(K1,K1)=1.0
      END DO
      
      DO K3=NDI+1,NTENS
          Sklij(K3,K3)=0.5E0
      END DO
      !defining IkI
      DO K1=1,NTENS
          DO K2=1,NTENS
              IkI(K1,K2)=0.0
          END DO
      END DO
      
      DO K1=1,NDI
          DO K2=1, NDI
              IkI(K1,K2)=1
          end DO
      END DO
  
    
      

      
C--------------------------------------------------------------------
C calculation for plasticity 
      DO K1=1,2
          
         
          STRES1=STRESS+MATMUL(DDSDDE,(DSTRAN-STATEV(27)*eta1-DTHRM))
          
          
          if (K1.EQ.1) then
            STN_E=MATMUL(Cef_i,STRES1)
            
            Ce=C_e1
          else
            STN_E=MATMUL(Cef_i,STRES1)
            
            Ce=C_e2
          end if
          
         
          
          CALL KTrace(STN_E,NTENS,NDI,tr_3)          
          CALL KTrace(STRES1,NTENS,NDI,sigm_p) ! sigm_p  is trace of stress 
       
         
          BULK=(ONE/THREE)*PROPS(K1)/(ONE-TWO*PROPS(K1+2))          
          Sdev1=STRES1-(ONE/THREE)*sigm_p*I_d
         
          Sdev_t=Sdev1
         
          CALL KNorm(Sdev_t,NTENS,NDI,sdev_n)! sdev_n is norm of Sdev_t
         
          S_vm_t=SQRT(THREE/TWO)*sdev_n   ! Trial von mises stress
         
          Q=PROPS(10)*STATEV(K1+2)
         
          Over_s=S_vm_t-(PROPS(K1+12)+Q)
          Syld_0=PROPS(K1+12)+Q
          
          
          
          CALL KMAC(Over_s,Ovr_sm)
          
      
          
C Check yielding
          if (Over_s. LE. ZERO) then
              del_mu=ZERO
              STATEV(K1+7)=STATEV(K1+7)+del_mu! mu_ij=mu_ij+del_mu
              
              DEELAS=(DSTRAN-STATEV(27)*eta1-DTHRM)
             
              
              
              DDSDDE=DDSDDE
              
              
             
              Syld_1=Syld_0
             
              if (K1.EQ.1) then
              STR=Sdev_t+BULK*tr_3*I_d
              else
              STR1=Sdev_t+BULK*tr_3*I_d
              end if
             
          else 
              del_mu=0.00001E0
              D_n=Sdev_t/sdev_n
              
              
              CALL KResd1(Sdev_t,del_mu,PROPS(12),PROPS(11),STATEV(K1),
     1        DTIME,EG(K1),C_eff,D_n,NTENS,NDI,PROPS(10),Q,PROPS(K1+12),Ce,resd1) !!calculate
              
              
              r_norm=SQRT(resd1*resd1)
              r_nor0=SQRT(resd1*resd1)
              
              
              
              iter8=0
              DO WHILE ((r_norm) .GT. TOLER2 .AND. iter8.LT.5000)
                  if (iter8.GT.4998) then
                      print*,"mu NR did not converge"
                      iter8=iter8+1
                      del_mu=0.00E0
                      
                       
                      
                  else
                  
                      
                          
                      
                      D_n=Sdev_t/sdev_n
                      CALL KDyad(D_n,NTENS,Dn_n)
                     
                      
                      CALL KDResd(Sdev_t,Ce,Dn_n,PROPS(10),del_mu,PROPS(12),
     1                PROPS(11),STATEV(K1),DTIME,NTENS,EG(1),C_eff,D_n,
     2                 NDI,Dresd1) 
                      
                      del_mu=del_mu-(resd1/Dresd1)
                      
                      CALL KResd1(Sdev_t,del_mu,PROPS(12),PROPS(11),STATEV(K1),
     1        DTIME,EG(K1),C_eff,D_n,NTENS,NDI,PROPS(10),Q,PROPS(K1+12),Ce,resd1) !!calculate
                      r_norm=SQRT(resd1*resd1)
                      iter8=iter8+1
                      
                      
                      
                      
                   
                      
                  end if
                  
              END DO
            
              
              D_n=Sdev_t/sdev_n
              Cef_n=MATMUL(C_eff,D_n)
              
             
              
              S_dev=Sdev_t-sqrt(3.0/2.0)*del_mu*Cef_n
              
              DPLAS=SQRT(THREE/TWO)*del_mu*D_n
              if (K1.EQ.1) then
              STR=S_dev+BULK*tr_3*I_d
              else
              STR1=S_dev+BULK*tr_3*I_d
              end if
              
              
              STATEV(K1+7)=STATEV(K1+7)+del_mu ! mu_ij=mu_ij+del_mu
              STATEV(K1+2)=STATEV(K1+2)+del_mu ! q_ij=q_ij+del_mu
              STATEV(K1+4)=STATEV(K1+4)+DPLAS(1) ! EPLS_ij=EPLS_ij+del_mu*n
              STATEV(2*(K1+4))=STATEV(2*(K1+4))+DPLAS(2)! EPLS_ij=EPLS_ij+del_mu*n
              STATEV(2*(K1+4)+1)=STATEV(2*(K1+4)+1)+DPLAS(3)
              STATEV(K1+16)=STATEV(K1+16)+DPLAS(4)
              STATEV(2*(K1+9)+1)=STATEV(2*(K1+9)+1)+DPLAS(5)
              STATEV(2*(K1+10))=STATEV(2*(K1+10))+DPLAS(6)
              
              
              
              
              

              Q=PROPS(10)*STATEV(K1+2)
         
              Syld_1=(PROPS(K1+12)+Q)
             
              
            
              
    
            
            end if
      END DO
                EPLAS(1)=aLmd_g*STATEV(5)+aLamda*STATEV(6)
              EPLAS(2)=aLmd_g*STATEV(10)+aLamda*STATEV(12)
              EPLAS(3)=aLmd_g*STATEV(11)+aLamda*STATEV(13)
              EPLAS(4)=aLmd_g*STATEV(17)+aLamda*STATEV(18)
              EPLAS(5)=aLmd_g*STATEV(21)+aLamda*STATEV(23)
              EPLAS(6)=aLmd_g*STATEV(22)+aLamda*STATEV(24)
              CALL KNorm(EPLAS,NTENS,NDI,EP_NM)
              STATEV(48)=EP_NM
              EPLAS(4)=2*EPLAS(4)
              EPLAS(5)=2*EPLAS(5)
              EPLAS(6)=2*EPLAS(6)
                STATEV(7)=EPLAS(1)
              STATEV(14)=EPLAS(2)
              STATEV(15)=EPLAS(3)
              STATEV(16)=EPLAS(4)
              STATEV(25)=EPLAS(5)
              STATEV(26)=EPLAS(6)
                 DEPLAS=EPLAS-OLDPL
              STATEV(28)=DEPLAS(1)
              STATEV(29)=DEPLAS(2)
              STATEV(30)=DEPLAS(3)
              STATEV(31)=DEPLAS(4)
              STATEV(32)=DEPLAS(5)
              STATEV(33)=DEPLAS(6)
                DEELAS=(DSTRAN-DEPLAS-STATEV(27)*eta1-DTHRM)
      
                DDSDDE=DDSDDE
     
                STRESS=STRESS+MATMUL(DDSDDE,DEELAS)

     
      
            
              EPLAS(4)=0.5*EPLAS(4)
              EPLAS(5)=0.5*EPLAS(5)
              EPLAS(6)=0.5*EPLAS(6)
              OLDPL(4)=0.5*OLDPL(4)
              OLDPL(5)=0.5*OLDPL(5)
              OLDPL(6)=0.5*OLDPL(6)
      
      
     
      

     
      

      RETURN
      END
 
      
C---------------------------------------------------------------


        SUBROUTINE UHARD(SYIELD,HARD,EQPLAS,EQPLASRT,TIME,DTIME,TEMP,
     1 DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,
     2 CMNAME,NSTATV,STATEV,NUMFIELDV,
     3 PREDEF,DPRED,NVALUE,TABLE)
        INCLUDE 'ABA_PARAM.INC'
        CHARACTER*80 CMNAME
        DIMENSION HARD(3),STATEV(NSTATV),TIME(*),
     1 PREDEF(NUMFIELDV),DPRED(*)
C
        DIMENSION TABLE(2, NVALUE)
C
        PARAMETER(ZERO=0.E0)
C
C SET YIELD STRESS TO LAST VALUE OF TABLE, HARDENING TO ZERO
C
        SYIELD=TABLE(1, NVALUE)
        HARD(1)=ZERO
C IF MORE THAN ONE ENTRY, SEARCH TABLE
C
        IF(NVALUE.GT.1) THEN
        DO K1=1, NVALUE-1
        EQPL1=TABLE(2,K1+1)
        IF(EQPLAS.LT.EQPL1) THEN
        EQPL0=TABLE(2, K1)
        IF(EQPL1.LE.EQPL0) THEN
        WRITE(7, 1)
  1     FORMAT(//, 30X, '***ERROR - PLASTIC STRAIN MUST BE ',
     1 'ENTERED IN ASCENDING ORDER')
        CALL XIT
        ENDIF
C
C CURRENT YIELD STRESS AND HARDENING
C
        DEQPL=EQPL1-EQPL0
        SYIEL0=TABLE(1, K1)
        SYIEL1=TABLE(1, K1+1)
        DSYIEL=SYIEL1-SYIEL0
        HARD(1)=DSYIEL/DEQPL
        SYIELD=SYIEL0+(EQPLAS-EQPL0)*HARD(1)
        GO TO 10
        ENDIF
        END DO
  10    CONTINUE
        ENDIF
        RETURN
        END
C---------------------------------------------------------------
C SUBROUTINES

C _________________________________________________________________
C 2. Subroutine for saturation function

        subroutine Ksatur(aLama,aLam_g,s1,s2,Hs)
C Performs a direct calculation of the saturation function for phase transformtion
          IMPLICIT REAL*8 (A-H,O-Z)
          
           Hs=((aLam_g)**s1 )*((aLama)**s2)
            RETURN
          END
C     _________________________________________________________________
C 3. Subroutine for platic energy function
C this is actually density times PSI_pl

        subroutine KPL_EN(qij,H_j,PSI_pl)
C Performs a direct calculation of the Rho times plastic energy
            IMPLICIT REAL*8 (A-H,O-Z)

            PSI_pl=(.5)*H_j*(qij**2.)
            RETURN
            END
C _________________________________________________________________
C 4. a Subroutine for Ce MATRIX CALCULATION (PLAIN STRAIN)


        subroutine KC_e(alm,alame,amue,numk,numb,Cek)
C Performs a direct calculation of the elastic stiffness matrix
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION Cek(numb,numb)    !! Elastic stiffness matrix for a phase Cek, output
           
            
  
            DO K1=1, numk
              DO K2=1, numk
                  if (K1==K2)then
                      Cek(K1, K2)=alm
                  else
                      Cek(K1, K2)=alame
                  end if
              END DO
            END DO
            DO K3=numk+1, numb
            DO K4=1,numb
            Cek(K3 ,K4)=0.
            Cek(K4 ,K3)=0.
            END DO
            Cek(K3 ,K3)=amue
            END DO
            RETURN
      END
C------------------------------------------------------------------
C 4. a Subroutine for Ce MATRIX CALCULATION (PLAIN STRAIN) spinel


        subroutine KC_e1(alm,alame,amue,numk,numb,Cek)
C Performs a direct calculation of the elastic stiffness matrix
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION Cek(numb,numb)    !! Elastic stiffness matrix for a phase Cek, output
           
            Cek(1,1)=300000.0E0 ;Cek(2,1)=118000.0E0;Cek(3,1)=118000.0E0
            Cek(1,2)=118000.0E0; Cek(2,2)=300000.0E0;Cek(3,2)=118000.0E0
            Cek(1,3)=118000.0E0;Cek(2,3)=118000.0E0;Cek(3,3)=300000.0E0
            Cek(1,4)=0.0E0;Cek(2,4)=0.0E0;Cek(3,4)=0.0E0
            Cek(1,5)=0.0E0;Cek(2,5)=0.0E0;Cek(3,5)=0.0E0
            Cek(1,6)=0.0E0;Cek(2,6)=0.0E0;Cek(3,6)=0.0E0
            Cek(4,1)=0.0E0 ;Cek(5,1)=0.0E0;Cek(6,1)=0.0E0
            Cek(4,2)=0.0E0; Cek(5,2)=0.0E0;Cek(6,2)=0.0E0
            Cek(4,3)=0.0E0;Cek(5,3)=0.0E0;Cek(6,3)=0.0E0
            Cek(4,4)=126000.0E0;Cek(5,4)=0.0E0;Cek(6,4)=0.0E0
            Cek(4,5)=0.0E0;Cek(5,5)=126000.0E0;Cek(6,5)=0.0E0
            Cek(4,6)=0.0E0;Cek(5,6)=0.0E0;Cek(6,6)=126000.0E0
        
            
  
            
            RETURN
      END
C----------------------------------------------
C4. a Subroutine for Ce MATRIX CALCULATION (PLAIN STRAIN)olivine


        subroutine KC_e2(alm,alame,amue,numk,numb,Cek)
C Performs a direct calculation of the elastic stiffness matrix
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION Cek(numb,numb)    !! Elastic stiffness matrix for a phase Cek, output
           
            Cek(1,1)=312000.0E0 ;Cek(2,1)=60000.0E0;Cek(3,1)=65000.0E0
            Cek(1,2)=60000.0E0; Cek(2,2)=187000.0E0;Cek(3,2)=66000.0E0
            Cek(1,3)=65000.0E0;Cek(2,3)=66000.0E0;Cek(3,3)=217000.0E0
            Cek(1,4)=0.0E0;Cek(2,4)=0.0E0;Cek(3,4)=0.0E0
            Cek(1,5)=0.0E0;Cek(2,5)=0.0E0;Cek(3,5)=0.0E0
            Cek(1,6)=0.0E0;Cek(2,6)=0.0E0;Cek(3,6)=0.0E0
            Cek(4,1)=0.0E0 ;Cek(5,1)=0.0E0;Cek(6,1)=0.0E0
            Cek(4,2)=0.0E0; Cek(5,2)=0.0E0;Cek(6,2)=0.0E0
            Cek(4,3)=0.0E0;Cek(5,3)=0.0E0;Cek(6,3)=0.0E0
            Cek(4,4)=71000.0E0;Cek(5,4)=0.0E0;Cek(6,4)=0.0E0
            Cek(4,5)=0.0E0;Cek(5,5)=66100.0E0;Cek(6,5)=0.0E0
            Cek(4,6)=0.0E0;Cek(5,6)=0.0E0;Cek(6,6)=57200.0E0
        
            
  
            
            RETURN
      END
C----------------------------------------------
C---------------------------------------------
C _________________________________________________________________
c SUBROUTINE FOR ROTATION MATRIX FOR EUER ANGLES
        SUBROUTINE ROT(Ang1,Ang2,Ang3,ROT_C)
         IMPLICIT REAL*8 (A-H,O-Z)   
        DIMENSION ROT_C(3,3)
        pi=acos(-1.0D0)
        Ang1=Ang1*(pi/180.0)
        Ang2=Ang2*(pi/180.0)
        Ang3=Ang3*(pi/180.0)
        ROT_C(1,1)=cos(Ang1)*cos(Ang3)-cos(Ang2)*sin(Ang1)*sin(Ang3)
        ROT_C(1,2)=cos(Ang1)*sin(Ang3)+cos(Ang3)*cos(Ang2)*sin(Ang1)
        ROT_C(1,3)=sin(Ang1)*sin(Ang2)
        ROT_C(2,1)=-cos(Ang3)*sin(Ang1)-cos(Ang1)*cos(Ang2)*sin(Ang3)
        ROT_C(2,2)=cos(Ang1)*cos(Ang3)*cos(Ang2)-sin(Ang1)*sin(Ang3)
        ROT_C(2,3)=cos(Ang1)*sin(Ang2)
        ROT_C(3,1)=sin(Ang3)*sin(Ang2)
        ROT_C(3,2)=-cos(Ang3)*sin(Ang2)
        ROT_C(3,3)=cos(Ang2)
        
        
        return
        
        end
c____________________________________________________________________
c DEFINING TRANSFROMATION m FOR 6X6 MATRIX
      SUBROUTINE trans(ROTA,AM_A)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ROTA(3,3)
      DIMENSION AM_A(6,6)
      DO J=1,3
          J1=1+J/3
          J2=2+J/2
          DO I=1,3
              I1=1+I/3
              I2=2+I/2
              AM_A(I,J)=ROTA(I,J)**2
              AM_A(I,J+3)=2.0E0*ROTA(I,J1)*ROTA(I,J2)
              AM_A(I+3,J)=ROTA(I1,J)*ROTA(I2,J)
              AM_A(I+3,J+3)=ROTA(I1,J1)*ROTA(I2,J2)+
     1        ROTA(I1,J2)*ROTA(I2,J1)
          END DO
      END DO
      RETURN
      END
      
c__________________________________________________________________
      SUBROUTINE TRANSP(XM_A,MSIZ,XM_AT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XM_AT(MSIZ,MSIZ)
      DIMENSION XM_A(MSIZ,MSIZ)
      DO I=1,MSIZ
          DO J=1,MSIZ
              XM_AT(I,J)=XM_A(J,I)
          END DO
      END DO
      RETURN
      END
      
c________________________________________________________________________
c---------------------------------------------------------------------
C 5. a Subroutine for calculting residual for phase transformation
        subroutine KResid(aL_frc,f_c,v_tr,DT,antr,Hs1,dlamd,
     1Resid1)
C Performs a direct calculation of the residual rij for given phase 
          IMPLICIT REAL*8 (A-H,O-Z)
         
            Resid1=(aL_frc-f_c)-((((v_tr/DT)**(1./antr))*((dlamd)
     1**(1./antr))*(1./Hs1)))

            RETURN
            END
           
C _________________________________________________________________
C 5. b Subroutine for calculting residual for phase transformation
        subroutine KRes(aL_frc,f_c,dlamd,Res1)
C Performs a direct calculation of the residual rij for given phase 
          IMPLICIT REAL*8 (A-H,O-Z)
         
            Res1=aL_frc-f_c-(((dlamd)))

            RETURN
            END
  
C _________________________________________________________________
C 6. a Subroutine for calculting Jacobian for phase transformation
        subroutine KJcob(aL_da,aL_g,an_tr,vtr,DTI,Hs2,dlad1,s11,s21,
     1dth01,dth11,aJac22)
          IMPLICIT REAL*8 (A-H,O-Z)
C Performs a direct calculation of the elastic stiffness matrix
            
    
  !! r11/lam11
      aJac22=(dth01-dth11)-((-(((vtr/DTI)**
     2(1.0E0/an_tr))*((dlad1)**(1.0E0/an_tr))*(1.0E0/Hs2/Hs2)*((aL_g)**
     3s11)*(s21*(aL_da**(s21-1.0E0)))))+(((vtr/DTI)**
     4(1.0E0/an_tr))*(1./Hs2)*(1.0E0/an_tr)*((dlad1)**((1.0E0-an_tr)/
     5an_tr))))
     
!     !!!  aJac22=(dth01-dth11)-((-(((vtr/DTI)**
     !!!2(1.0E0/an_tr))*((dlad1)**(1.0E0/an_tr))*(1.0E0/Hs2/Hs2)*((aL_g)**
     !!!3s11)*(s21*(aL_da**(s21-1.0E0))))))


     
            RETURN
            END
            
C _________________________________________________________________
C 6b. a Subroutine for calculting Jacobian for  delta lambda
        subroutine KdJcb(aL_da,aL_g,an_tr,vtr,DTI,Hs2,dlad1,s11,s21,
     1aJacdL)
     
C Performs a direct calculation of the elastic stiffness matrix
            IMPLICIT REAL*8 (A-H,O-Z)
   
    

      aJacdL=-(((1.0E0/an_tr)*((vtr/DTI)**(1.0E0/an_tr))*
     1(1.0E0/Hs2)*((dlad1)**((1.0E0-an_tr)/an_tr)))-(((vtr/DTI)**
     2(1.0E0/an_tr))*((dlad1)**(1.0E0/an_tr))*(1.0E0/Hs2/Hs2)*((aL_g)**
     3s11)*(s21*(aL_da**(s21-1.0E0)))))
             


            RETURN
            END
C-------------------------------------
C 6. c Subroutine for calculting Jacobian for phase transformation
        subroutine KJcb(dth01,dth11,aJac22)
          IMPLICIT REAL*8 (A-H,O-Z)
C Performs a direct calculation of the elastic stiffness matrix
            
    
  !! r11/lam11
      aJac22=(dth01-dth11)
     
!     !!!  aJac22=(dth01-dth11)-((-(((vtr/DTI)**
     !!!2(1.0E0/an_tr))*((dlad1)**(1.0E0/an_tr))*(1.0E0/Hs2/Hs2)*((aL_g)**
     !!!3s11)*(s21*(aL_da**(s21-1.0E0))))))


     
            RETURN
            END

C_________________________________________________________________
C 7. a Subroutine for calculting trace of stress o strain
        subroutine KTrace(AB,NTES,NI,trace1)
C Performs a  calculation of the trace of matrix for 2d case
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION AB(NTES)   !! Vector
           

            trace1=0
             DO K1=1,NI
                trace1=trace1+AB(K1)
             END DO
   
        RETURN
        END
  
C______________DO___________________________________________________
C 8. a Subroutine for calculting norm of vector
              Subroutine KNorm(xn,NTE,ND,anormx) 
            !! Performs a  calculation of the norm of a 3×1 matrix.
              IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION xn(NTE) !! actual value,i/p

          
            anomx1=0.
            DO K1=1,ND
            anomx1=anomx1+(xn(K1)**2.)
            END DO
            DO K1=(ND+1),NTE
            anomx1=anomx1+(2*(xn(K1)**2.))
            END DO
            
            anormx=SQRT(anomx1)
            RETURN
            END
C ______________________________________________________________
C 9. a Subroutine for calculting residual for  viscoplasticity 
              Subroutine KResd1(aSdv_t,ad_mu,v_pl,an_pl,aLa,
     1aDT,s_mu,Ceff,aDn,NTE,ND,H_1,aQ,aY,aCe,ared1)
              
            !! Performs a  calculation of the norm of a 3×1 matrix.
            IMPLICIT REAL*8 (A-H,O-Z)
              DIMENSION Ceff(NTE,NTE)
              DIMENSION aCe(NTE,NTE)
              DIMENSION aDn(NTE)
              DIMENSION Dmn(NTE)
              DIMENSION Cefn(NTE)
              DIMENSION Cefa(NTE)
              DIMENSION aSdv_t(NTE)
              DIMENSION ax(NTE)
              DIMENSION ax1(NTE)
              Dmn=sqrt(3./2.)*ad_mu*aDn
              Cefn=aLa*MATMUL(Ceff,Dmn)
              Cefa=MATMUL(aCe,Dmn)
              ax=aSdv_t-Cefn
              ax1=aSdv_t-Cefa
              call KNorm(Cefn,NTE,ND,cfnrm)
              call KNorm(ax,NTE,ND,xnrm)
              call KNorm(ax1,NTE,ND,xnrm1)
              
     
            ared1=(sqrt(3./2.)*xnrm1)-(H_1*ad_mu)-aQ-aY-
     1       ((ad_mu*v_pl/aDT/aLa)**(1./an_pl))

            !!!!!ared1=ov_s-(H_1*ad_mu)-
     !!!!!!!1       ((ad_mu*v_pl/aDT/aLa)**(1./an_pl))
           
            RETURN
            END
C _________________________________________________________________
C 10. a Subroutine for calculting derivative of residual for  viscoplasticity
    
               Subroutine KDResd(aSdv_t,C2,Ann,H_1,d_mu,vpl,an_pl,a_Lm,
     1aDTIM,kNTEN,s_mu,Ceff,aDn,ND,Dresd) 
               
             
           IMPLICIT REAL*8 (A-H,O-Z)
             
            
            DIMENSION C2(kNTEN,kNTEN) !! actual value,i/p
            DIMENSION Ann(kNTEN,kNTEN) !! actual value,i/p
            DIMENSION Ceff(kNTEN,kNTEN)
            DIMENSION Xnn(kNTEN,kNTEN)
              DIMENSION aDn(kNTEN)
              DIMENSION Cefn(kNTEN)
              DIMENSION Tax(kNTEN)
              DIMENSION Tax1(kNTEN)
              DIMENSION Tax2(kNTEN)
              DIMENSION Tadn(kNTEN)
              DIMENSION Tadmn(kNTEN)
              DIMENSION xvec(kNTEN)
              DIMENSION xvec1(kNTEN)
              DIMENSION xvec2(kNTEN)
              DIMENSION avec(kNTEN)
              DIMENSION avec1(kNTEN)
              DIMENSION aSdv_t(kNTEN)
              DIMENSION cfn(kNTEN)
              DIMENSION cfa(kNTEN)
              Tadmn=sqrt(3./2.)*d_mu*aDn
              !Tax=MATMUL(Ceff,Tadmn)
              !Tax1=(a_Lm)*MATMUL(Ceff,Tadmn)
              Tax2=MATMUL(C2,Tadmn)
              !xvec1=aSdv_t-Tax1
              xvec2=aSdv_t-Tax2
              !call KNorm(xvec1,kNTEN,ND,xnrm1)
              call KNorm(xvec2,kNTEN,ND,xnrm2)
              !Tadn=sqrt(3./2.)*aDn
              !Cefn=MATMUL(Ceff,Tadn)
              !call KNorm(Cefn,kNTEN,ND,cfnrm)
              !call KNorm(Tax,kNTEN,ND,xnrm)
              !xvec=Tax/xnrm
              !avec=xvec1/xnrm1
              avec1=xvec2/xnrm2
              !cfn=MATMUL(Ceff,aDn)
              cfa=MATMUL(C2,aDn)
              !call KDyad2(Tadn,xvec,kNTEN,Xnn)
            !CALL Kdcont(C2,Ann,kNTEN,dc1)
            !CALL Kdcont(Ceff,Xnn,kNTEN,Cnn)
            
     
        Dresd=(-(3./2.)*(dot_product(cfa,avec1)) - H_1) - (1./an_pl)*
     1((d_mu*vpl/aDTIM/a_Lm)**((1.-an_pl)/an_pl)) *(vpl/aDTIM/a_Lm)

            !!!!!Dresd= - H_1- (1./an_pl)*
     !!!!!!!1((d_mu*vpl/aDTIM/a_Lm)**((1.-an_pl)/an_pl)) *(vpl/aDTIM/a_Lm)
            
      
     
            RETURN
            END
C    _________________________________________________________________
C 11. a Subroutine for calculting dyadic product of vectors
              Subroutine KDyad(AK,KNS,BK) 
             IMPLICIT REAL*8 (A-H,O-Z)
     
            DIMENSION AK(KNS,1) !! actual value,i/p
             
            DIMENSION A_t(1,KNS)
             
            dimension  BK(KNS,KNS) !! output  
            Do K1=1,KNS
            A_t(1,K1)=AK(K1,1)
            End do
            BK=Matmul(AK,A_t)
            RETURN
            END
C----------------------------------------------------------------
C11.b
            Subroutine KDyad2(AK,CK,KNS,BK) 
             IMPLICIT REAL*8 (A-H,O-Z)
     
            DIMENSION AK(KNS,1) !! actual value,i/p
            DIMENSION CK(KNS,1) 
            DIMENSION C_t(1,KNS)
             
            dimension  BK(KNS,KNS) !! output  
            Do K1=1,KNS
            C_t(1,K1)=CK(K1,1)
            End do
            BK=Matmul(AK,C_t)
            RETURN
            END
    

C ______________________________________________________________________________________

C 12. a Subroutine for calculting derivative of thermodynamic force !! needs modification for plastic strain and total strainof variant 
          !!Subroutine Kd_thr(Cef,Cc_in,Cs_in,aSTRS,aSTRAN,aSTRN,aEPLAS
     !!1         ,aEPLSC,aEPLSs,KNTS,dthr) 
           Subroutine Kd_thr(Cef,Cc_in,Cs_in,aSTRS,aSTRAN,aSTRN,aEPLAS
     1         ,aEPLSC,aEPLSs,KNTS,etaj,etak,etaf,dthr)
     
                  IMPLICIT REAL*8 (A-H,O-Z)
                  DIMENSION Cef(KNTS,KNTS) !! actual value,i/p
                  DIMENSION Cc_in(KNTS,KNTS) !! actual value,i/p from def of Lam
                  DIMENSION Cs_in(KNTS,KNTS) !! actual value,i/p from diff wrt lamda
                  DIMENSION aSTRS(KNTS) !! actual value,i/p
                  DIMENSION aSTRAN(KNTS) !! actual value,i/p
                  DIMENSION aSTRN(KNTS) !! actual value,i/p
                  DIMENSION aEPLAS(KNTS) !! actual value,i/p
                  DIMENSION aEPLSC(KNTS) !! e plastc from def of lam
                  DIMENSION aEPLSs(KNTS) !! e plastic ffrom diff wrt small lamda
                  DIMENSION etaj(KNTS) !! e plastic ffrom diff wrt small lamda
                  DIMENSION etaf(KNTS) !! e plastic ffrom diff wrt small lamda
                  DIMENSION etak(KNTS) !! e plastic ffrom diff wrt small lamda
                  
                  
      
                  
                  CALL Kdcont(Cef,Cs_in,KNTS,dcont1) 
           
          
     !!             d_thr1=dot_product((MATMUL(aSTRS,Cc_in)),(MATMUL(Cef
     !!1            ,aEPLSs)+MATMUL((dcont1*Cef),(aSTRAN-aEPLAS))))
     !!             d_thr2= dot_product((MATMUL(Cef,aEPLSs)+MATMUL(
     !!1            (dcont1*Cef),(aSTRAN-aEPLAS))),(aSTRN+aEPLSC))
                  
                  
                  
                  d_thr1=dot_product((MATMUL(aSTRS,Cc_in)),(MATMUL(Cef
     1            ,aEPLSs)+MATMUL(Cef,etaj)+
     2            MATMUL((dcont1*Cef),(aSTRAN-aEPLAS-etaf))))
                  d_thr2= dot_product((MATMUL(Cef,aEPLSs)+
     1             MATMUL(Cef,etaj)+MATMUL((dcont1*Cef),
     2             (aSTRAN-aEPLAS-etaf))),(aSTRN+aEPLSC+etak))
     
                  dthr=d_thr1+d_thr2

    
     
    
            RETURN
            END
C _______________________________________________________________________________________
C 13. a Subroutine for calculting double contraction of matrices
              Subroutine kdcont(aNA,aNB,NS1,dcont) 
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION aNA(NS1,NS1) !! actual value,i/p
            DIMENSION aNB(NS1,NS1) !! actual value,i/p
            
            
            

            dcont=0
            Do K1=1,NS1
                DO K2=1,NS1
            dcont=dcont+aNA(K1,K2)*aNB(K1,K2)
                End DO
            END DO
    
            RETURN
            END
C_____________________________________________________________________________
C--------------------------------------------------------
         Subroutine Kdsdl(Cef,C_in,aSTRAN,aEPLAS,aEPLSs,KNTS,etaj,etaf,
     1     ds_dl) 
     
                  IMPLICIT REAL*8 (A-H,O-Z)
                  DIMENSION Cef(KNTS,KNTS) !! actual value,i/p
                  
                  DIMENSION C_in(KNTS,KNTS) !! actual value,i/p from diff wrt lamda
                  
                  DIMENSION aSTRAN(KNTS) !! actual value,i/p
                  
                  DIMENSION aEPLAS(KNTS) !! actual value,i/p
                  
                  DIMENSION aEPLSs(KNTS) !! e plastic ffrom diff wrt small lamda
                  
                  DIMENSION ds_dl(KNTS)
                  DIMENSION etaj(KNTS)
                  DIMENSION etaf(KNTS)
      
                  
                  CALL Kdcont(Cef,C_in,KNTS,dcont1) 
           
          

                  ds_dl= -(MATMUL(Cef,aEPLSs)+MATMUL(Cef,etaj)+MATMUL(
     1            (dcont1*Cef),(aSTRAN-aEPLAS-etaf)))


    
     
    
            RETURN
            END
C-------------------------------------------------
C subroutine for macauley brackets
            Subroutine KMAC(aP,bP) 
    
     
            IMPLICIT REAL*8 (A-H,O-Z)
            
            bP=(aP+abs(aP))/2.0E0
            
            RETURN
            END
C__________________________________________________________________________           
            
             Subroutine inverse(a,c,n1)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n1,n1) - array of coefficients for matrix A
! n1      - dimension
! output ...
! c(n1,n1) - inverse matrix of A
! comments ...
! the original matrix a(n1,n1) will be destroyed 
! during the calculation
!===========================================================
      !implicit none 
        
      IMPLICIT REAL*8 (A-H,O-Z)
      !integer n1 H here replaced L matrix in LU decompositin coz L is by default integer type in fortran
      DIMENSION a(n1,n1), c(n1,n1)
      DIMENSION H(n1,n1), U(n1,n1), b(n1), d(n1), x(n1)
      !Real(8) coeff
      !integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
      H=0.
      U=0.
      b=0.

! step 1: forward elimination
      do k=1, n1-1
       do i=k+1,n1
        coeff=a(i,k)/a(k,k)
        H(i,k) = coeff
      do j=k+1,n1
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
       end do
       end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
      do i=1,n1
      H(i,i) = 1.
      end do
! U matrix is the upper triangular part of A
      do j=1,n1
      do i=1,j
       U(i,j) = a(i,j)
      end do
      end do

! Step 3: compute columns of the inverse matrix C
      do k=1,n1
      b(k)=1.
      d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
       do i=2,n1
       d(i)=b(i)
       do j=1,i-1
       d(i) = d(i) - H(i,j)*d(j)
      end do
      end do
! Step 3b: Solve Ux=d using the back substitution
      x(n1)=d(n1)/U(n1,n1)
      do i = n1-1,1,-1
      x(i) = d(i)
      do j=n1,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
      end do
      x(i) = x(i)/u(i,i)
      end do
! Step 3c: fill the solutions x(n1) into column k of C
      do i=1,n1  
      c(i,k) = x(i)
      end do
      b(k)=0.
      end do
      end subroutine inverse
C--------------------------------------------------------------------------------------------

      
C------------------------------------------------------------------------
      function projec(v,vmin,vmax) result(vt)
    
      IMPLICIT REAL*8 (A-H,O-Z)
      

      vt=min(max(v,vmin),vmax)
      end function
            
   
C---------------------------------END----------------------