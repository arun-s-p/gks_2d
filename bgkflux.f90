!
! --- *******************************************************************
!     *                                                                 *
!     *   FILE: bgkflux.f90                                             *
!     *                                                                 *
! --- *******************************************************************
!
      SUBROUTINE BGKFLUX(WLP,WRP,W1P,W2P,DSL,DSR,DT,MU,FS,GAMMA,PRN)
!
! --- *******************************************************************
!     *                                                                 *
!     * -- DESCRIPTION:                                                 *
!     *    ------------                                                 *
!     *    COMPUTES THE FLUXES USING THE BGK METHOD                     *
!     *    INCLUDES NON EQUILIBRIUM TERMS                               *
!     *                                                                 *
! --- *******************************************************************
!
!!$      USE FLO_PARAM
!!$
!!$      USE DEFINES
!!$!
! --- *******************************************************************
!
      IMPLICIT NONE
!
! --- *******************************************************************
!
! --- IN
      REAL, INTENT(IN), DIMENSION(4) :: WLP,WRP,W1P,W2P ! LEFT/RIGHT AND CELL STATES
      REAL, INTENT(IN)               :: DSL,DSR     ! CC-FACE LENGTH LEFT/RIGHT
      REAL, INTENT(IN)               :: DT,MU       ! TIME STEP, KIN. VISC.
      REAL, INTENT(IN)               :: GAMMA,PRN
!
! --- OUT
      REAL, INTENT(OUT),DIMENSION(4) :: FS          ! GAS-KINETIC FLUX
!
! --- LOCAL
      INTEGER :: I

      REAL, DIMENSION(0:6)  :: UPL,UMR,UFL,UFR,VFL,VFR,WFL,WFR
      REAL, DIMENSION(0:6)  :: UP0,UM0,UF0,VF0,WF0
      REAL, DIMENSION(5)    :: W1,W2,WL,WR
      REAL, DIMENSION(5)    :: DWL,DWR,ACL,ACR,WW0
      REAL, DIMENSION(5)    :: DW0L,DW0R,AC0L,AC0R
      REAL, DIMENSION(5)    :: SLOPA,AA,ATL,ATR
      REAL, DIMENSION(5)    :: TERM1L,TERM2L,TERM3L,TERM5L,TERM6L,RHSL
      REAL, DIMENSION(5)    :: TERM1R,TERM2R,TERM3R,TERM5R,TERM6R,RHSR
      REAL, DIMENSION(5)    :: ATERML
      REAL, DIMENSION(5)    :: ATERMR
      REAL, DIMENSION(5)    :: TERM1,TERM2,TERM3,TERM4,TERM5,TERM6


      REAL    :: RHOL,RHOR,RHO1,RHO2
      REAL    :: UXL,UXR,UX1,UX2
      REAL    :: VXL,VXR,VX1,VX2
      REAL    :: WXL,WXR,WX1,WX2
      REAL    :: PL,PR,P1,P2
      REAL    :: EL,ER,E1,E2
      REAL    :: LAMBDAL,LAMBDAR,LAMINVL,LAMINVR

      REAL    :: XFL2,XFR2,XFL4,XFR4,XVWL2,XVWR2
      REAL    :: XF02,XF04,XVW02

      REAL    :: UX0,VX0,WX0,PP0
      REAL    :: LAMBDA0,LAMINV0

      REAL    :: TAU,VIS,eps,TAUI

      REAL    :: ALPHA1,ALPHA2,ALPHA3,ALPHA4,ALPHA5,ALPHA6,ALPHA7
      REAL    :: BETA1,BETA2
      REAL    :: EXTAU
      REAL    :: GAM0,GAM2,GAM5
      REAL    :: COEFF1,COEFF2,COEFF3

      REAL    :: DENSAV,UMOMAV,VMOMAV,WMOMAV,ENERAV
      REAL    :: UAV,VAV,WAV
      REAL    :: EKIN,Q

      REAL    :: GMG,GM1

      REAL    :: ERFC

      REAL, PARAMETER ::   CK  = 0.0
      REAL :: FLO_PI
!
! --- *******************************************************************
!

      GMG = (GAMMA -1.)/GAMMA
      GM1 = GAMMA -1.

      WL(1) = WLP(1) ;          WR(1) = WRP(1) ;                                         
      WL(2) = WLP(1)*WLP(2);    WR(2) = WRP(1)*WRP(2);                                 
      WL(3) = WLP(1)*WLP(3);    WR(3) = WRP(1)*WRP(3);                                 
      WL(4) = 0.;               WR(4) = 0.;                                            

      WL(5) = WLP(4)/GM1 + .5*WLP(1)*(WLP(2)**2 + WLP(3)**2);       
      WR(5) = WRP(4)/GM1 + .5*WRP(1)*(WRP(2)**2 + WRP(3)**2);


      W1(1) = W1P(1) ;          W2(1) = W2P(1) ;                                         
      W1(2) = W1P(1)*W1P(2);    W2(2) = W2P(1)*W2P(2);                                 
      W1(3) = W1P(1)*W1P(3);    W2(3) = W2P(1)*W2P(3);                                 
      W1(4) = 0.;               W2(4) = 0.;                                            

      W1(5) = W1P(4)/GM1 + .5*W1P(1)*(W1P(2)**2 + W1P(3)**2);       
      W2(5) = W2P(4)/GM1 + .5*W2P(1)*(W2P(2)**2 + W2P(3)**2);

      FLO_PI = 4.*ATAN(1.)
      

      RHOL  = WL(1);                 RHOR  = WR(1)        
      UXL   = WL(2)/WL(1);           UXR   = WR(2)/WR(1)  
      VXL   = WL(3)/WL(1);           VXR   = WR(3)/WR(1)  
      WXL   = WL(4)/WL(1);           WXR   = WR(4)/WR(1)  
      EL    = WL(5);                 ER    = WR(5)        
                                                          
      PL = GM1*(EL-.5*WL(1)*(UXL**2 + VXL**2 + WXL**2))
      PR = GM1*(ER-.5*WR(1)*(UXR**2 + VXR**2 + WXR**2))

!!$      RHO1  = W1(1);                 RHO2  = W2(1)        
!!$      UX1   = W1(2)/W1(1);           UX2   = W2(2)/W2(1)  
!!$      VX1   = W1(3)/W1(1);           VX2   = W2(3)/W2(1)  
!!$      WX1   = W1(4)/W1(1);           WX2   = W2(4)/W2(1)  
!!$      E1    = W1(5);                 E2    = W2(5)        
!!$                                                               
!!$      P1 = GM1*(E1-.5*W1(1)*(UX1**2 + VX1**2 + WX1**2))
!!$      P2 = GM1*(E2-.5*W2(1)*(UX2**2 + VX2**2 + WX2**2))
      
!!$      EL = HL - PL
!!$      ER = HR - PR
!!$      E1 = H1 - P1
!!$      E2 = H2 - P2

!!$      FR(1) = UXR*RHOR;               FL(1) = UXL*WL(1)          
!!$      FR(2) = UXR*WR(2)  + PR;        FL(2) = UXL*WL(2)  + PL    
!!$      FR(3) = UXR*WR(3);              FL(3) = UXL*WL(3)          
!!$      FR(4) = UXR*WR(4);              FL(4) = UXL*WL(4)          
!!$      FR(5) = UXR*HR;                 FL(5) = UXL*HL             
!!$
!!$      RETURN


! --------------------------------------------------------------
!     COMPUTE LAMBDA = m / (2 k T) = RHO / (2 P) 
!     THIS IS THE DAMPING FACTOR IN THE MAXWELLIAN DISTRIBUTION
! --------------------------------------------------------------
      LAMBDAL = .5*WL(1)/PL;        LAMINVL = 1./LAMBDAL  
      LAMBDAR = .5*WR(1)/PR;        LAMINVR = 1./LAMBDAR  
 
!----------------------------------------------------------------
!			MOMENTS
!	F- (-inf->inf)	M- (-inf->0)	P- (0->inf)
!----------------------------------------------------------------
      UPL(0) = .5*ERFC(-UXL*SQRT(LAMBDAL))
      UPL(1) = UXL*UPL(0) + 0.50*EXP(-LAMBDAL*UXL*UXL)/SQRT(LAMBDAL*FLO_PI)
      UMR(0) = .5*ERFC(UXR*SQRT(LAMBDAR))
      UMR(1) = UXR*UMR(0) - 0.50*EXP(-LAMBDAR*UXR*UXR)/SQRT(LAMBDAR*FLO_PI)

      UFL(0) = 1.;                  UFL(1) = UXL*UFL(0);     
      UFR(0) = 1.;                  UFR(1) = UXR*UFR(0);     

      VFL(0) = 1.;                  VFL(1) = VXL*VFL(0) 
      VFR(0) = 1.;                  VFR(1) = VXR*VFR(0) 

      WFL(0) = 1.;                  WFL(1) = WXL*WFL(0)
      WFR(0) = 1.;                  WFR(1) = WXR*WFR(0)

      DO I = 2,6
         UPL(I) = UXL*UPL(I-1) + (FLOAT(I)-1.)*0.5*UPL(I-2)*LAMINVL
         UFL(I) = UXL*UFL(I-1) + (FLOAT(I)-1.)*0.5*UFL(I-2)*LAMINVL
         VFL(I) = VXL*VFL(I-1) + (FLOAT(I)-1.)*0.5*VFL(I-2)*LAMINVL
         WFL(I) = WXL*WFL(I-1) + (FLOAT(I)-1.)*0.5*WFL(I-2)*LAMINVL

         UMR(I) = UXR*UMR(I-1) + (FLOAT(I)-1.)*0.5*UMR(I-2)*LAMINVR
         UFR(I) = UXR*UFR(I-1) + (FLOAT(I)-1.)*0.5*UFR(I-2)*LAMINVR
         VFR(I) = VXR*VFR(I-1) + (FLOAT(I)-1.)*0.5*VFR(I-2)*LAMINVR
         WFR(I) = WXR*WFR(I-1) + (FLOAT(I)-1.)*0.5*WFR(I-2)*LAMINVR
      END DO

! --- CONTRIBUTION FROM INTERNAL ENERGY (HERE JUST ROTATION)
      XFL2  =  0.5*CK*LAMINVL
      XFL4  =  0.25*CK*(CK+2.0)*LAMINVL**2
      XVWL2 =  XFL2 + VFL(2) + WFL(2)
     
      XFR2  =  0.5*CK*LAMINVR
      XFR4  =  0.25*CK*(CK+2.0)*LAMINVR**2
      XVWR2 =  XFR2 + VFR(2) + WFR(2)

!------------ Eqn 4.18, collapsing left and right at the cell interface      
      WW0(1)  = RHOL*UPL(0)        + RHOR*UMR(0)
      WW0(2)  = RHOL*UPL(1)        + RHOR*UMR(1)
      WW0(3)  = RHOL*UPL(0)*VFL(1) + RHOR*UMR(0)*VFR(1)
      WW0(4)  = RHOL*UPL(0)*WFL(1) + RHOR*UMR(0)*WFR(1)
      WW0(5)  = (RHOL*(UPL(2) + UPL(0)*XVWL2) + &
                 RHOR*(UMR(2) + UMR(0)*XVWR2)  ) * .5

      UX0  = WW0(2)/WW0(1)
      VX0  = WW0(3)/WW0(1)
      WX0  = WW0(4)/WW0(1)

      PP0  = GM1*(WW0(5) - .5*WW0(1)*(UX0**2 +VX0**2 +WX0**2))
      LAMBDA0 = .5*WW0(1)/PP0

!viscous      
!      VIS = 0.0005*DT!MU/PP0;    eps = 2.0
      VIS = MU/PP0;    eps = 2.0
      TAU = VIS+DT*ABS(PL-PR)/(PL+PR)    !  * eps


!inviscid
!	TAUI = 5.*(ABS((RHOL/LAMBDAL)-(RHOR/LAMBDAR)))/(ABS((RHOL/LAMBDAL)+(RHOR/LAMBDAR)))
!	TAU = .05*DT+DT*MIN(1.,TAUI)


      DO I=1,5
         DWL(I) = (WL(I) - W1(I))/DSL
         DWR(I) = (W2(I) - WR(I))/DSR
      ENDDO

      CALL DXE3(UXL,VXL,WXL,LAMBDAL,DWL,ACL,GAMMA)
      CALL DXE3(UXR,VXR,WXR,LAMBDAR,DWR,ACR,GAMMA)

! --- FIND NON-EQUILIBRIUM TERM COEFFS (SIMPLIFICATION POSSIBLE..LATER)

      RHSL(1) = ACL(1)*UFL(1) + ACL(2)*UFL(2) + ACL(3)*UFL(1)*VFL(1) +  &
                  ACL(4)*UFL(1)*WFL(1) +  ACL(5)*(UFL(3) + UFL(1)*XVWL2)
      RHSL(2) = ACL(1)*UFL(2) + ACL(2)*UFL(3) + ACL(3)*UFL(2)*VFL(1) +  &
                ACL(4)*UFL(2)*WFL(1) +  ACL(5)*(UFL(4) + UFL(2)*XVWL2)
      RHSL(3) = ACL(1)*UFL(1)*VFL(1) + ACL(2)*UFL(2)*VFL(1) + &
                ACL(3)*UFL(1)*VFL(2) + ACL(4)*UFL(1)*VFL(1)*WFL(1) +  &
                ACL(5)*(UFL(3)*VFL(1) + UFL(1)*VFL(3) + UFL(1)*VFL(1)*WFL(2) + &
                        UFL(1)*VFL(1)*XFL2)
      RHSL(4) = ACL(1)*UFL(1)*WFL(1) + ACL(2)*UFL(2)*WFL(1) + &
                ACL(3)*UFL(1)*VFL(1)*WFL(1) + ACL(4)*UFL(1)*WFL(2) +  &
                ACL(5)*(UFL(3)*WFL(1) + UFL(1)*WFL(3) + UFL(1)*WFL(1)*VFL(2) + &
                        UFL(1)*WFL(1)*XFL2)
      RHSL(5) = ACL(1)*(UFL(3) + UFL(1)*XVWL2) + &
                  ACL(2)*(UFL(4) + UFL(2)*XVWL2) + &
                  ACL(3)*(UFL(3)*VFL(1) + UFL(1)*VFL(3) + UFL(1)*VFL(1)*WFL(2) + &
                          UFL(1)*VFL(1)*XFL2) + &
                  ACL(4)*(UFL(3)*WFL(1) +UFL(1)*WFL(3) +UFL(1)*WFL(1)*VFL(2) + &
                          UFL(1)*WFL(1)*XFL2) + &
                  ACL(5)*(UFL(5) + UFL(1)*VFL(4) + UFL(1)*WFL(4) + UFL(1)*XFL4 + &
                          2.*UFL(3)*XVWL2 +  &
                          2.*UFL(1)*(VFL(2)*XFL2 +WFL(2)*XFL2 + VFL(2)*WFL(2)))
      RHSL(5) = .5*RHSL(5)

      RHSR(1) = ACR(1)*UFR(1) + ACR(2)*UFR(2) + ACR(3)*UFR(1)*VFR(1) +  &
                ACR(4)*UFR(1)*WFR(1) +  ACR(5)*(UFR(3) + UFR(1)*XVWR2)
      RHSR(2) = ACR(1)*UFR(2) + ACR(2)*UFR(3) + ACR(3)*UFR(2)*VFR(1) +  &
                ACR(4)*UFR(2)*WFR(1) +  ACR(5)*(UFR(4) + UFR(2)*XVWR2)
      RHSR(3) = ACR(1)*UFR(1)*VFR(1) + ACR(2)*UFR(2)*VFR(1) + &
                ACR(3)*UFR(1)*VFR(2) + ACR(4)*UFR(1)*VFR(1)*WFR(1) +  &
                ACR(5)*(UFR(3)*VFR(1) + UFR(1)*VFR(3) + UFR(1)*VFR(1)*WFR(2) + &
                        UFR(1)*VFR(1)*XFR2)
      RHSR(4) = ACR(1)*UFR(1)*WFR(1) + ACR(2)*UFR(2)*WFR(1) + &
                ACR(3)*UFR(1)*VFR(1)*WFR(1) + ACR(4)*UFR(1)*WFR(2) +  &
                ACR(5)*(UFR(3)*WFR(1) + UFR(1)*WFR(3) + UFR(1)*WFR(1)*VFR(2) + &
                        UFR(1)*WFR(1)*XFR2)
      RHSR(5) =ACR(1)*(UFR(3) + UFR(1)*XVWR2) + ACR(2)*(UFR(4) + UFR(2)*XVWR2) + &
                ACR(3)*(UFR(3)*VFR(1) + UFR(1)*VFR(3) + UFR(1)*VFR(1)*WFR(2) + &
                        UFR(1)*VFR(1)*XFR2) + &
                ACR(4)*(UFR(3)*WFR(1) + UFR(1)*WFR(3) + UFR(1)*WFR(1)*VFR(2) + &
                        UFR(1)*WFR(1)*XFR2) + &
                ACR(5)*(UFR(5) + UFR(1)*VFR(4) + UFR(1)*WFR(4) + UFR(1)*XFR4 + &
                        2.*UFR(3)*XVWR2 +  &
                        2.*UFR(1)*(VFR(2)*XFR2 +WFR(2)*XFR2 + VFR(2)*WFR(2)))
      RHSR(5) = .5*RHSR(5)

      DO I = 1,5
         RHSL(I) = -RHSL(I)!/RHOL
         RHSR(I) = -RHSR(I)!/RHOR
      END DO

      CALL DXE3(UXL,VXL,WXL,LAMBDAL,RHSL,ATL,GAMMA)
      CALL DXE3(UXR,VXR,WXR,LAMBDAR,RHSR,ATR,GAMMA)

      EXTAU  = EXP(-DT/TAU)
      
      ALPHA4  = TAU*(1.0-EXTAU)
      ALPHA5  = TAU*( DT*EXTAU - ALPHA4 )    
      ALPHA1  = DT - ALPHA4  
!!$      ALPHA2  = TAU*( 2.*ALPHA4 - DT*(1. + EXTAU) )
!!$      ALPHA3  = 0.5*DT*DT - TAU*( DT - ALPHA4 )

      ALPHA6  = TAU*( DT*EXTAU - 2.*ALPHA4 )    
!      alpha6 = alpha5
      ALPHA7  = TAU*ALPHA4


!!$      BETA1  = ALPHA2/(TAU*ALPHA1)  
!!$      BETA2  = ALPHA5/(TAU*ALPHA1)  

      TERM1L(1) = ACL(1)*UPL(1) + ACL(2)*UPL(2) + ACL(3)*UPL(1)*VFL(1) +  &
                  ACL(4)*UPL(1)*WFL(1) +  ACL(5)*(UPL(3) + UPL(1)*XVWL2)
      TERM1L(2) = ACL(1)*UPL(2) + ACL(2)*UPL(3) + ACL(3)*UPL(2)*VFL(1) +  &
                ACL(4)*UPL(2)*WFL(1) +  ACL(5)*(UPL(4) + UPL(2)*XVWL2)
      TERM1L(3) = ACL(1)*UPL(1)*VFL(1) + ACL(2)*UPL(2)*VFL(1) + &
                ACL(3)*UPL(1)*VFL(2) + ACL(4)*UPL(1)*VFL(1)*WFL(1) +  &
                ACL(5)*(UPL(3)*VFL(1) + UPL(1)*VFL(3) + UPL(1)*VFL(1)*WFL(2) + &
                        UPL(1)*VFL(1)*XFL2)
      TERM1L(4) = ACL(1)*UPL(1)*WFL(1) + ACL(2)*UPL(2)*WFL(1) + &
                ACL(3)*UPL(1)*VFL(1)*WFL(1) + ACL(4)*UPL(1)*WFL(2) +  &
                ACL(5)*(UPL(3)*WFL(1) + UPL(1)*WFL(3) + UPL(1)*WFL(1)*VFL(2) + &
                        UPL(1)*WFL(1)*XFL2)
      TERM1L(5) = ACL(1)*(UPL(3) + UPL(1)*XVWL2) + &
                  ACL(2)*(UPL(4) + UPL(2)*XVWL2) + &
                  ACL(3)*(UPL(3)*VFL(1) + UPL(1)*VFL(3) + UPL(1)*VFL(1)*WFL(2) + &
                          UPL(1)*VFL(1)*XFL2) + &
                  ACL(4)*(UPL(3)*WFL(1) +UPL(1)*WFL(3) +UPL(1)*WFL(1)*VFL(2) + &
                          UPL(1)*WFL(1)*XFL2) + &
                  ACL(5)*(UPL(5) + UPL(1)*VFL(4) + UPL(1)*WFL(4) + UPL(1)*XFL4 + &
                          2.*UPL(3)*XVWL2 +  &
                          2.*UPL(1)*(VFL(2)*XFL2 +WFL(2)*XFL2 + VFL(2)*WFL(2)))
      TERM1L(5) = .5*TERM1L(5)


      TERM1R(1) = ACR(1)*UMR(1) + ACR(2)*UMR(2) + ACR(3)*UMR(1)*VFR(1) +  &
                ACR(4)*UMR(1)*WFR(1) +  ACR(5)*(UMR(3) + UMR(1)*XVWR2)
      TERM1R(2) = ACR(1)*UMR(2) + ACR(2)*UMR(3) + ACR(3)*UMR(2)*VFR(1) +  &
                ACR(4)*UMR(2)*WFR(1) +  ACR(5)*(UMR(4) + UMR(2)*XVWR2)
      TERM1R(3) = ACR(1)*UMR(1)*VFR(1) + ACR(2)*UMR(2)*VFR(1) + &
                ACR(3)*UMR(1)*VFR(2) + ACR(4)*UMR(1)*VFR(1)*WFR(1) +  &
                ACR(5)*(UMR(3)*VFR(1) + UMR(1)*VFR(3) + UMR(1)*VFR(1)*WFR(2) + &
                        UMR(1)*VFR(1)*XFR2)
      TERM1R(4) = ACR(1)*UMR(1)*WFR(1) + ACR(2)*UMR(2)*WFR(1) + &
                ACR(3)*UMR(1)*VFR(1)*WFR(1) + ACR(4)*UMR(1)*WFR(2) +  &
                ACR(5)*(UMR(3)*WFR(1) + UMR(1)*WFR(3) + UMR(1)*WFR(1)*VFR(2) + &
                        UMR(1)*WFR(1)*XFR2)
      TERM1R(5) =ACR(1)*(UMR(3) + UMR(1)*XVWR2) + ACR(2)*(UMR(4) + UMR(2)*XVWR2) + &
                ACR(3)*(UMR(3)*VFR(1) + UMR(1)*VFR(3) + UMR(1)*VFR(1)*WFR(2) + &
                        UMR(1)*VFR(1)*XFR2) + &
                ACR(4)*(UMR(3)*WFR(1) + UMR(1)*WFR(3) + UMR(1)*WFR(1)*VFR(2) + &
                        UMR(1)*WFR(1)*XFR2) + &
                ACR(5)*(UMR(5) + UMR(1)*VFR(4) + UMR(1)*WFR(4) + UMR(1)*XFR4 + &
                        2.*UMR(3)*XVWR2 +  &
                        2.*UMR(1)*(VFR(2)*XFR2 +WFR(2)*XFR2 + VFR(2)*WFR(2)))
      TERM1R(5) = .5*TERM1R(5)

      DO I=1,5
!         TERM1(I) =  ALPHA5*(TERM1L(I) + TERM1R(I))
         TERM1(I) =  ALPHA6*(TERM1L(I) + TERM1R(I))
      END DO


! -----------------------------------------------------------------------------
!     TERMS FOR f_0 (g_L AND g_R EXPANSIONS)
! -----------------------------------------------------------------------------
      TERM2L(1)= RHOL*ALPHA4*UPL(1) + ALPHA6 * ( &
                 ACL(1)*UPL(2) + ACL(2)*UPL(3) + ACL(3)*UPL(2)*VFL(1) +  &
                 ACL(4)*UPL(2)*WFL(1) +  ACL(5)*(UPL(4) + UPL(2)*XVWL2) )
      TERM2L(2)= RHOL*ALPHA4*UPL(2) + ALPHA6 * ( &
                 ACL(1)*UPL(3) + ACL(2)*UPL(4) + ACL(3)*UPL(3)*VFL(1) +  &
                 ACL(4)*UPL(3)*WFL(1) +  ACL(5)*(UPL(5) + UPL(3)*XVWL2))
      TERM2L(3)= RHOL*ALPHA4*UPL(1)*VFL(1) + ALPHA6 * ( &
                 ACL(1)*UPL(2)*VFL(1) + ACL(2)*UPL(3)*VFL(1) + &
                 ACL(3)*UPL(2)*VFL(2) + ACL(4)*UPL(2)*VFL(1)*WFL(1) +  &
                 ACL(5)*(UPL(4)*VFL(1) +UPL(2)*(VFL(3)+VFL(1)*WFL(2)+VFL(1)*XFL2)))
      TERM2L(4)= RHOL*ALPHA4*UPL(1)*WFL(1) + ALPHA6 * ( &
                 ACL(1)*UPL(2)*WFL(1) + ACL(2)*UPL(3)*WFL(1) + &
                 ACL(3)*UPL(2)*VFL(1)*WFL(1) + ACL(4)*UPL(2)*WFL(2) +  &
                 ACL(5)*(UPL(4)*WFL(1) +UPL(2)*(WFL(3) +WFL(1)*VFL(2)+WFL(1)*XFL2)))
      TERM2L(5)= RHOL*ALPHA4*.5*(UPL(3) + UPL(1)*XVWL2) + .5*ALPHA6 * ( &
                 ACL(1)*(UPL(4) + UPL(2)*XVWL2) + ACL(2)*(UPL(5) + UPL(3)*XVWL2) + &
                 ACL(3)*(UPL(4)*VFL(1) + UPL(2)*VFL(3) + UPL(2)*VFL(1)*WFL(2) + &
                         UPL(2)*VFL(1)*XFL2) + &
                 ACL(4)*(UPL(4)*WFL(1) + UPL(2)*WFL(3) + UPL(2)*WFL(1)*VFL(2) + &
                 UPL(2)*WFL(1)*XFL2) + &
                 ACL(5)*(UPL(6) + UPL(2)*VFL(4) + UPL(2)*WFL(4) + UPL(2)*XFL4 + &
                         2.*UPL(4)*XVWL2 +  &
                         2.*UPL(2)*(VFL(2)*XFL2 +WFL(2)*XFL2 + VFL(2)*WFL(2))))


      TERM2R(1)= RHOR*ALPHA4*UMR(1) + ALPHA6 * ( &
                 ACR(1)*UMR(2) + ACR(2)*UMR(3) + ACR(3)*UMR(2)*VFR(1) +  &
                 ACR(4)*UMR(2)*WFR(1) +  ACR(5)*(UMR(4) + UMR(2)*XVWR2) )
      TERM2R(2)= RHOR*ALPHA4*UMR(2) + ALPHA6 * ( &
                 ACR(1)*UMR(3) + ACR(2)*UMR(4) + ACR(3)*UMR(3)*VFR(1) +  &
                 ACR(4)*UMR(3)*WFR(1) +  ACR(5)*(UMR(5) + UMR(3)*XVWR2))
      TERM2R(3)= RHOR*ALPHA4*UMR(1)*VFR(1) + ALPHA6 * ( &
                 ACR(1)*UMR(2)*VFR(1) + ACR(2)*UMR(3)*VFR(1) + &
                 ACR(3)*UMR(2)*VFR(2) + ACR(4)*UMR(2)*VFR(1)*WFR(1) +  &
                 ACR(5)*(UMR(4)*VFR(1) +UMR(2)*(VFR(3)+VFR(1)*WFR(2)+VFR(1)*XFR2)))
      TERM2R(4)= RHOR*ALPHA4*UMR(1)*WFR(1) + ALPHA6 * ( &
                 ACR(1)*UMR(2)*WFR(1) + ACR(2)*UMR(3)*WFR(1) + &
                 ACR(3)*UMR(2)*VFR(1)*WFR(1) + ACR(4)*UMR(2)*WFR(2) +  &
                 ACR(5)*(UMR(4)*WFR(1) +UMR(2)*(WFR(3)+WFR(1)*VFR(2)+WFR(1)*XFR2)))
      TERM2R(5)= RHOR*ALPHA4*.5*(UMR(3) + UMR(1)*XVWR2) + .5*ALPHA6 * ( &
                 ACR(1)*(UMR(4) + UMR(2)*XVWR2) + ACR(2)*(UMR(5) + UMR(3)*XVWR2) + &
                 ACR(3)*(UMR(4)*VFR(1) + UMR(2)*VFR(3) + UMR(2)*VFR(1)*WFR(2) + &
                         UMR(2)*VFR(1)*XFR2) + &
                 ACR(4)*(UMR(4)*WFR(1) + UMR(2)*WFR(3) + UMR(2)*WFR(1)*VFR(2) + &
                         UMR(2)*WFR(1)*XFR2) + &
                 ACR(5)*(UMR(6) + UMR(2)*VFR(4) + UMR(2)*WFR(4) + UMR(2)*XFR4 + &
                         2.*UMR(4)*XVWR2 +  &
                         2.*UMR(2)*(VFR(2)*XFR2 +WFR(2)*XFR2 + VFR(2)*WFR(2))))

      DO I=1,5
         TERM2(I) = TERM2L(I) + TERM2R(I)
      END DO

      DO I=1,5
         DW0L(I) = (WW0(I) - W1(I))/DSL
         DW0R(I) = (W2(I) - WW0(I))/DSR
      END DO
      
!!$      DO I=1,5
!!$         DW0L(I) = (W2(I) - W1(I))/(DSL+DSR)
!!$         DW0R(I) = (W2(I) - W1(I))/(DSL+DSR)
!!$      END DO
      

      CALL DXE3(UX0,VX0,WX0,LAMBDA0,DW0L,AC0L,GAMMA)
      CALL DXE3(UX0,VX0,WX0,LAMBDA0,DW0R,AC0R,GAMMA)


      LAMINV0 = 1./LAMBDA0

      UP0(0) = 0.5*ERFC(-UX0*SQRT(LAMBDA0))
      UP0(1) = UX0*UP0(0) + 0.50*EXP(-LAMBDA0*UX0*UX0)/SQRT(LAMBDA0*FLO_PI)

      UM0(0) = 0.5*ERFC( UX0*SQRT(LAMBDA0))
      UM0(1) = UX0*UM0(0) - 0.50*EXP(-LAMBDA0*UX0*UX0)/SQRT(LAMBDA0*FLO_PI)

      UF0(0) = 1.;      UF0(1) = UX0*UF0(0)
      VF0(0) = 1.;      VF0(1) = VX0*VF0(0)
      WF0(0) = 1.;      WF0(1) = WX0*WF0(0)
      
      DO I = 2,6
         UP0(I) = UX0*UP0(I-1) + (FLOAT(I)-1.)*0.5*UP0(I-2)*LAMINV0
         UM0(I) = UX0*UM0(I-1) + (FLOAT(I)-1.)*0.5*UM0(I-2)*LAMINV0
         UF0(I) = UX0*UF0(I-1) + (FLOAT(I)-1.)*0.5*UF0(I-2)*LAMINV0
         VF0(I) = VX0*VF0(I-1) + (FLOAT(I)-1.)*0.5*VF0(I-2)*LAMINV0
         WF0(I) = WX0*WF0(I-1) + (FLOAT(I)-1.)*0.5*WF0(I-2)*LAMINV0
      END DO

      XF02  =  0.5*CK*LAMINV0
      XF04  =  0.25*CK*(CK+2.0)*LAMINV0**2
      XVW02 =  XF02 + VF0(2) + WF0(2)
     

      TERM3L(1)= AC0L(1)*UP0(1) + AC0L(2)*UP0(2) + AC0L(3)*UP0(1)*VF0(1) +  &
                AC0L(4)*UP0(1)*WF0(1) +  AC0L(5)*(UP0(3) + UP0(1)*XVW02)
      TERM3L(2) = AC0L(1)*UP0(2) + AC0L(2)*UP0(3) + AC0L(3)*UP0(2)*VF0(1) +  &
                AC0L(4)*UP0(2)*WF0(1) +  AC0L(5)*(UP0(4) + UP0(2)*XVW02)
      TERM3L(3) = AC0L(1)*UP0(1)*VF0(1) + AC0L(2)*UP0(2)*VF0(1) + &
                AC0L(3)*UP0(1)*VF0(2) + AC0L(4)*UP0(1)*VF0(1)*WF0(1) +  &
                AC0L(5)*(UP0(3)*VF0(1) + UP0(1)*VF0(3) + UP0(1)*VF0(1)*WF0(2) + &
                        UP0(1)*VF0(1)*XF02)
      TERM3L(4) = AC0L(1)*UP0(1)*WF0(1) + AC0L(2)*UP0(2)*WF0(1) + &
                AC0L(3)*UP0(1)*VF0(1)*WF0(1) + AC0L(4)*UP0(1)*WF0(2) +  &
                AC0L(5)*(UP0(3)*WF0(1) + UP0(1)*WF0(3) + UP0(1)*WF0(1)*VF0(2) + &
                        UP0(1)*WF0(1)*XF02)
      TERM3L(5) =AC0L(1)*(UP0(3) +UP0(1)*XVW02) +AC0L(2)*(UP0(4) + UP0(2)*XVW02) + &
                AC0L(3)*(UP0(3)*VF0(1) + UP0(1)*VF0(3) + UP0(1)*VF0(1)*WF0(2) + &
                        UP0(1)*VF0(1)*XF02) + &
                AC0L(4)*(UP0(3)*WF0(1) + UP0(1)*WF0(3) + UP0(1)*WF0(1)*VF0(2) + &
                        UP0(1)*WF0(1)*XF02) + &
                AC0L(5)*(UP0(5) + UP0(1)*VF0(4) + UP0(1)*WF0(4) + UP0(1)*XF04 + &
                        2.*UP0(3)*XVW02 +  &
                        2.*UP0(1)*(VF0(2)*XF02 +WF0(2)*XF02 + VF0(2)*WF0(2)))
      TERM3L(5) = 0.5*TERM3L(5)


      TERM3R(1) = AC0R(1)*UM0(1) + AC0R(2)*UM0(2) + AC0R(3)*UM0(1)*VF0(1) +  &
                  AC0R(4)*UM0(1)*WF0(1) +  AC0R(5)*(UM0(3) + UM0(1)*XVW02)
      TERM3R(2) = AC0R(1)*UM0(2) + AC0R(2)*UM0(3) + AC0R(3)*UM0(2)*VF0(1) +  &
                  AC0R(4)*UM0(2)*WF0(1) +  AC0R(5)*(UM0(4) + UM0(2)*XVW02)
      TERM3R(3) = AC0R(1)*UM0(1)*VF0(1) + AC0R(2)*UM0(2)*VF0(1) + &
                  AC0R(3)*UM0(1)*VF0(2) + AC0R(4)*UM0(1)*VF0(1)*WF0(1) +  &
                  AC0R(5)*(UM0(3)*VF0(1) + UM0(1)*VF0(3) + UM0(1)*VF0(1)*WF0(2) + &
                           UM0(1)*VF0(1)*XF02)
      TERM3R(4) = AC0R(1)*UM0(1)*WF0(1) + AC0R(2)*UM0(2)*WF0(1) + &
                  AC0R(3)*UM0(1)*VF0(1)*WF0(1) + AC0R(4)*UM0(1)*WF0(2) +  &
                  AC0R(5)*(UM0(3)*WF0(1) + UM0(1)*WF0(3) + UM0(1)*WF0(1)*VF0(2) + &
                           UM0(1)*WF0(1)*XF02)
      TERM3R(5) = AC0R(1)*(UM0(3) +UM0(1)*XVW02) +AC0R(2)*(UM0(4) +UM0(2)*XVW02) + &
                  AC0R(3)*(UM0(3)*VF0(1) + UM0(1)*VF0(3) + UM0(1)*VF0(1)*WF0(2) + &
                           UM0(1)*VF0(1)*XF02) + &
                  AC0R(4)*(UM0(3)*WF0(1) + UM0(1)*WF0(3) + UM0(1)*WF0(1)*VF0(2) + &
                           UM0(1)*WF0(1)*XF02) + &
                  AC0R(5)*(UM0(5) + UM0(1)*VF0(4) + UM0(1)*WF0(4) + UM0(1)*XF04 + &
                           2.*UM0(3)*XVW02 +  &
                           2.*UM0(1)*(VF0(2)*XF02 +WF0(2)*XF02 + VF0(2)*WF0(2)))
      TERM3R(5) = 0.5*TERM3R(5)

      DO I=1,5
         TERM3(I) = TERM3L(I) + TERM3R(I)
      END DO
      

      ATERML(1) = ATL(1)*UPL(0) + ATL(2)*UPL(1) + ATL(3)*UPL(0)*VFL(1) +  &
                  ATL(4)*UPL(0)*WFL(1) +  ATL(5)*(UPL(2) + UPL(0)*XVWL2)
      ATERML(2) = ATL(1)*UPL(1) + ATL(2)*UPL(2) + ATL(3)*UPL(1)*VFL(1) +  &
                  ATL(4)*UPL(1)*WFL(1) +  ATL(5)*(UPL(3) + UPL(1)*XVWL2)
      ATERML(3) = ATL(1)*UPL(0)*VFL(1) + ATL(2)*UPL(1)*VFL(1) + &
                  ATL(3)*UPL(0)*VFL(2) + ATL(4)*UPL(0)*VFL(1)*WFL(1) +  &
                  ATL(5)*(UPL(2)*VFL(1) + UPL(0)*VFL(3) + UPL(0)*VFL(1)*WFL(2) + &
                        UPL(0)*VFL(1)*XFL2)
      ATERML(4) = ATL(1)*UPL(0)*WFL(1) + ATL(2)*UPL(1)*WFL(1) + &
                  ATL(3)*UPL(0)*VFL(1)*WFL(1) + ATL(4)*UPL(0)*WFL(2) +  &
                  ATL(5)*(UPL(2)*WFL(1) + UPL(0)*WFL(3) + UPL(0)*WFL(1)*VFL(2) + &
                        UPL(0)*WFL(1)*XFL2)
      ATERML(5) = ATL(1)*(UPL(2) + UPL(0)*XVWL2) + &
                  ATL(2)*(UPL(3) + UPL(1)*XVWL2) + &
                  ATL(3)*(UPL(2)*VFL(1) + UPL(0)*VFL(3) + UPL(0)*VFL(1)*WFL(2) + &
                          UPL(0)*VFL(1)*XFL2) + &
                  ATL(4)*(UPL(2)*WFL(1) +UPL(0)*WFL(3) +UPL(0)*WFL(1)*VFL(2) + &
                          UPL(0)*WFL(1)*XFL2) + &
                  ATL(5)*(UPL(4) + UPL(0)*VFL(4) + UPL(0)*WFL(4) + UPL(0)*XFL4 + &
                          2.*UPL(2)*XVWL2 +  &
                          2.*UPL(0)*(VFL(2)*XFL2 +WFL(2)*XFL2 + VFL(2)*WFL(2)))
      ATERML(5) = .5*ATERML(5)


      ATERMR(1) = ATR(1)*UMR(0) + ATR(2)*UMR(1) + ATR(3)*UMR(0)*VFR(1) +  &
                  ATR(4)*UMR(0)*WFR(1) +  ATR(5)*(UMR(2) + UMR(0)*XVWR2)
      ATERMR(2) = ATR(1)*UMR(1) + ATR(2)*UMR(2) + ATR(3)*UMR(1)*VFR(1) +  &
                  ATR(4)*UMR(1)*WFR(1) +  ATR(5)*(UMR(3) + UMR(1)*XVWR2)
      ATERMR(3) = ATR(1)*UMR(0)*VFR(1) + ATR(2)*UMR(1)*VFR(1) + &
                  ATR(3)*UMR(0)*VFR(2) + ATR(4)*UMR(0)*VFR(1)*WFR(1) +  &
                  ATR(5)*(UMR(2)*VFR(1) + UMR(0)*VFR(3) + UMR(0)*VFR(1)*WFR(2) + &
                        UMR(0)*VFR(1)*XFR2)
      ATERMR(4) = ATR(1)*UMR(0)*WFR(1) + ATR(2)*UMR(1)*WFR(1) + &
                  ATR(3)*UMR(0)*VFR(1)*WFR(1) + ATR(4)*UMR(0)*WFR(2) +  &
                  ATR(5)*(UMR(2)*WFR(1) + UMR(0)*WFR(3) + UMR(0)*WFR(1)*VFR(2) + &
                        UMR(0)*WFR(1)*XFR2)
      ATERMR(5) = ATR(1)*(UMR(2) + UMR(0)*XVWR2) + &
                  ATR(2)*(UMR(3) + UMR(1)*XVWR2) + &
                  ATR(3)*(UMR(2)*VFR(1) + UMR(0)*VFR(3) + UMR(0)*VFR(1)*WFR(2) + &
                          UMR(0)*VFR(1)*XFR2) + &
                  ATR(4)*(UMR(2)*WFR(1) +UMR(0)*WFR(3) +UMR(0)*WFR(1)*VFR(2) + &
                          UMR(0)*WFR(1)*XFR2) + &
                  ATR(5)*(UMR(4) + UMR(0)*VFR(4) + UMR(0)*WFR(4) + UMR(0)*XFR4 + &
                          2.*UMR(2)*XVWR2 +  &
                          2.*UMR(0)*(VFR(2)*XFR2 +WFR(2)*XFR2 + VFR(2)*WFR(2)))
      ATERMR(5) = .5*ATERMR(5)


      GAM0 = TAU*(DT - ALPHA4)
      GAM2 = -GAM0 - ALPHA5
      GAM5 = TAU*ALPHA4

      DO I=1,5
         SLOPA(I) = (TERM3(I)*GAM2 +TERM1(I) -GAM5*(ATERML(I) + ATERMR(I)))/GAM0
!         SLOPA(I) = (TERM3(I)*GAM2 +TERM1(I))/GAM0
      END DO
      
      CALL DXE3(UX0,VX0,WX0,LAMBDA0,SLOPA,AA,GAMMA)      

      TERM4(1)= AA(1)*UF0(1) + AA(2)*UF0(2) + AA(3)*UF0(1)*VF0(1) +  &
                AA(4)*UF0(1)*WF0(1) +  AA(5)*(UF0(3) + UF0(1)*XVW02)
      TERM4(2) = AA(1)*UF0(2) + AA(2)*UF0(3) + AA(3)*UF0(2)*VF0(1) +  &
                AA(4)*UF0(2)*WF0(1) +  AA(5)*(UF0(4) + UF0(2)*XVW02)
      TERM4(3) = AA(1)*UF0(1)*VF0(1) + AA(2)*UF0(2)*VF0(1) + &
                AA(3)*UF0(1)*VF0(2) + AA(4)*UF0(1)*VF0(1)*WF0(1) +  &
                AA(5)*(UF0(3)*VF0(1) + UF0(1)*VF0(3) + UF0(1)*VF0(1)*WF0(2) + &
                        UF0(1)*VF0(1)*XF02)
      TERM4(4) = AA(1)*UF0(1)*WF0(1) + AA(2)*UF0(2)*WF0(1) + &
                AA(3)*UF0(1)*VF0(1)*WF0(1) + AA(4)*UF0(1)*WF0(2) +  &
                AA(5)*(UF0(3)*WF0(1) + UF0(1)*WF0(3) + UF0(1)*WF0(1)*VF0(2) + &
                        UF0(1)*WF0(1)*XF02)
      TERM4(5) =AA(1)*(UF0(3) +UF0(1)*XVW02) +AA(2)*(UF0(4) + UF0(2)*XVW02) + &
                AA(3)*(UF0(3)*VF0(1) + UF0(1)*VF0(3) + UF0(1)*VF0(1)*WF0(2) + &
                        UF0(1)*VF0(1)*XF02) + &
                AA(4)*(UF0(3)*WF0(1) + UF0(1)*WF0(3) + UF0(1)*WF0(1)*VF0(2) + &
                        UF0(1)*WF0(1)*XF02) + &
                AA(5)*(UF0(5) + UF0(1)*VF0(4) + UF0(1)*WF0(4) + UF0(1)*XF04 + &
                        2.*UF0(3)*XVW02 +  &
                        2.*UF0(1)*(VF0(2)*XF02 +WF0(2)*XF02 + VF0(2)*WF0(2)))
      TERM4(5) = 0.5*TERM4(5)



      TERM5L(1) = AC0L(1)*UP0(2) + AC0L(2)*UP0(3) + AC0L(3)*UP0(2)*VF0(1) +  &
                  AC0L(4)*UP0(2)*WF0(1) +  AC0L(5)*(UP0(4) + UP0(2)*XVW02)
      TERM5L(2) = AC0L(1)*UP0(3) + AC0L(2)*UP0(4) + AC0L(3)*UP0(3)*VF0(1) +  &
                  AC0L(4)*UP0(3)*WF0(1) +  AC0L(5)*(UP0(5) + UP0(3)*XVW02)
      TERM5L(3) = AC0L(1)*UP0(2)*VF0(1) + AC0L(2)*UP0(3)*VF0(1) + &
                  AC0L(3)*UP0(2)*VF0(2) + AC0L(4)*UP0(2)*VF0(1)*WF0(1) +  &
                  AC0L(5)*(UP0(4)*VF0(1) + UP0(2)*VF0(3) + UP0(2)*VF0(1)*WF0(2) + &
                          UP0(2)*VF0(1)*XF02)
      TERM5L(4) = AC0L(1)*UP0(2)*WF0(1) + AC0L(2)*UP0(3)*WF0(1) + &
                  AC0L(3)*UP0(2)*VF0(1)*WF0(1) + AC0L(4)*UP0(2)*WF0(2) +  &
                  AC0L(5)*(UP0(4)*WF0(1) + UP0(2)*WF0(3) + UP0(2)*WF0(1)*VF0(2) + &
                          UP0(2)*WF0(1)*XF02)
      TERM5L(5) =AC0L(1)*(UP0(4) + UP0(2)*XVW02) +AC0L(2)*(UP0(5) +UP0(3)*XVW02) + &
                 AC0L(3)*(UP0(4)*VF0(1) + UP0(2)*VF0(3) + UP0(2)*VF0(1)*WF0(2) + &
                          UP0(2)*VF0(1)*XF02) + &
                 AC0L(4)*(UP0(4)*WF0(1) + UP0(2)*WF0(3) + UP0(2)*WF0(1)*VF0(2) + &
                         UP0(2)*WF0(1)*XF02) + &
                 AC0L(5)*(UP0(6) + UP0(2)*VF0(4) + UP0(2)*WF0(4) + UP0(2)*XF04 + &
                         2.*UP0(4)*XVW02 +  &
                         2.*UP0(2)*(VF0(2)*XF02 +WF0(2)*XF02 + VF0(2)*WF0(2)))
      TERM5L(5) = .5*TERM5L(5)



      TERM5R(1) = AC0R(1)*UM0(2) + AC0R(2)*UM0(3) + AC0R(3)*UM0(2)*VF0(1) +  &
                  AC0R(4)*UM0(2)*WF0(1) +  AC0R(5)*(UM0(4) + UM0(2)*XVW02)
      TERM5R(2) = AC0R(1)*UM0(3) + AC0R(2)*UM0(4) + AC0R(3)*UM0(3)*VF0(1) +  &
                  AC0R(4)*UM0(3)*WF0(1) +  AC0R(5)*(UM0(5) + UM0(3)*XVW02)
      TERM5R(3) = AC0R(1)*UM0(2)*VF0(1) + AC0R(2)*UM0(3)*VF0(1) + &
                  AC0R(3)*UM0(2)*VF0(2) + AC0R(4)*UM0(2)*VF0(1)*WF0(1) +  &
                  AC0R(5)*(UM0(4)*VF0(1) + UM0(2)*VF0(3) + UM0(2)*VF0(1)*WF0(2) + &
                          UM0(2)*VF0(1)*XF02)
      TERM5R(4) = AC0R(1)*UM0(2)*WF0(1) + AC0R(2)*UM0(3)*WF0(1) + &
                  AC0R(3)*UM0(2)*VF0(1)*WF0(1) + AC0R(4)*UM0(2)*WF0(2) +  &
                  AC0R(5)*(UM0(4)*WF0(1) + UM0(2)*WF0(3) + UM0(2)*WF0(1)*VF0(2) + &
                          UM0(2)*WF0(1)*XF02)
      TERM5R(5) =AC0R(1)*(UM0(4) + UM0(2)*XVW02) +AC0R(2)*(UM0(5) +UM0(3)*XVW02) + &
                 AC0R(3)*(UM0(4)*VF0(1) + UM0(2)*VF0(3) + UM0(2)*VF0(1)*WF0(2) + &
                          UM0(2)*VF0(1)*XF02) + &
                 AC0R(4)*(UM0(4)*WF0(1) + UM0(2)*WF0(3) + UM0(2)*WF0(1)*VF0(2) + &
                         UM0(2)*WF0(1)*XF02) + &
                 AC0R(5)*(UM0(6) + UM0(2)*VF0(4) + UM0(2)*WF0(4) + UM0(2)*XF04 + &
                         2.*UM0(4)*XVW02 +  &
                         2.*UM0(2)*(VF0(2)*XF02 +WF0(2)*XF02 + VF0(2)*WF0(2)))
      TERM5R(5) = .5*TERM5R(5)

      DO I=1,5
         TERM5(I) = TERM5L(I) + TERM5R(I)
      END DO


      TERM6L(1) = ATL(1)*UPL(1) + ATL(2)*UPL(2) + ATL(3)*UPL(1)*VFL(1) +  &
                  ATL(4)*UPL(1)*WFL(1) +  ATL(5)*(UPL(3) + UPL(1)*XVWL2)
      TERM6L(2) = ATL(1)*UPL(2) + ATL(2)*UPL(3) + ATL(3)*UPL(2)*VFL(1) +  &
                  ATL(4)*UPL(2)*WFL(1) +  ATL(5)*(UPL(4) + UPL(2)*XVWL2)
      TERM6L(3) = ATL(1)*UPL(1)*VFL(1) + ATL(2)*UPL(2)*VFL(1) + &
                  ATL(3)*UPL(1)*VFL(2) + ATL(4)*UPL(1)*VFL(1)*WFL(1) +  &
                  ATL(5)*(UPL(3)*VFL(1) + UPL(1)*VFL(3) + UPL(1)*VFL(1)*WFL(2) + &
                        UPL(1)*VFL(1)*XFL2)
      TERM6L(4) = ATL(1)*UPL(1)*WFL(1) + ATL(2)*UPL(2)*WFL(1) + &
                  ATL(3)*UPL(1)*VFL(1)*WFL(1) + ATL(4)*UPL(1)*WFL(2) +  &
                  ATL(5)*(UPL(3)*WFL(1) + UPL(1)*WFL(3) + UPL(1)*WFL(1)*VFL(2) + &
                        UPL(1)*WFL(1)*XFL2)
      TERM6L(5) = ATL(1)*(UPL(3) +UPL(1)*XVWL2) +ATL(2)*(UPL(4) + UPL(2)*XVWL2) + &
                  ATL(3)*(UPL(3)*VFL(1) + UPL(1)*VFL(3) + UPL(1)*VFL(1)*WFL(2) + &
                        UPL(1)*VFL(1)*XFL2) + &
                  ATL(4)*(UPL(3)*WFL(1) + UPL(1)*WFL(3) + UPL(1)*WFL(1)*VFL(2) + &
                        UPL(1)*WFL(1)*XFL2) + &
                  ATL(5)*(UPL(5) + UPL(1)*VFL(4) + UPL(1)*WFL(4) + UPL(1)*XFL4 + &
                        2.*UPL(3)*XVWL2 +  &
                        2.*UPL(1)*(VFL(2)*XFL2 +WFL(2)*XFL2 + VFL(2)*WFL(2)))
      TERM6L(5) = 0.5*TERM6L(5)


      TERM6R(1) = ATR(1)*UMR(1) + ATR(2)*UMR(2) + ATR(3)*UMR(1)*VFR(1) +  &
                  ATR(4)*UMR(1)*WFR(1) +  ATR(5)*(UMR(3) + UMR(1)*XVWR2)
      TERM6R(2) = ATR(1)*UMR(2) + ATR(2)*UMR(3) + ATR(3)*UMR(2)*VFR(1) +  &
                  ATR(4)*UMR(2)*WFR(1) +  ATR(5)*(UMR(4) + UMR(2)*XVWR2)
      TERM6R(3) = ATR(1)*UMR(1)*VFR(1) + ATR(2)*UMR(2)*VFR(1) + &
                  ATR(3)*UMR(1)*VFR(2) + ATR(4)*UMR(1)*VFR(1)*WFR(1) +  &
                  ATR(5)*(UMR(3)*VFR(1) + UMR(1)*VFR(3) + UMR(1)*VFR(1)*WFR(2) + &
                          UMR(1)*VFR(1)*XFR2)
      TERM6R(4) = ATR(1)*UMR(1)*WFR(1) + ATR(2)*UMR(2)*WFR(1) + &
                  ATR(3)*UMR(1)*VFR(1)*WFR(1) + ATR(4)*UMR(1)*WFR(2) +  &
                  ATR(5)*(UMR(3)*WFR(1) + UMR(1)*WFR(3) + UMR(1)*WFR(1)*VFR(2) + &
                          UMR(1)*WFR(1)*XFR2)
      TERM6R(5) = ATR(1)*(UMR(3) +UMR(1)*XVWR2) +ATR(2)*(UMR(4) + UMR(2)*XVWR2) + &
                  ATR(3)*(UMR(3)*VFR(1) + UMR(1)*VFR(3) + UMR(1)*VFR(1)*WFR(2) + &
                          UMR(1)*VFR(1)*XFR2) + &
                  ATR(4)*(UMR(3)*WFR(1) + UMR(1)*WFR(3) + UMR(1)*WFR(1)*VFR(2) + &
                          UMR(1)*WFR(1)*XFR2) + &
                  ATR(5)*(UMR(5) + UMR(1)*VFR(4) + UMR(1)*WFR(4) + UMR(1)*XFR4 + &
                          2.*UMR(3)*XVWR2 +  &
                          2.*UMR(1)*(VFR(2)*XFR2 +WFR(2)*XFR2 + VFR(2)*WFR(2)))
      TERM6R(5) = 0.5*TERM6R(5)


      DO I=1,5
         TERM6(I) = TERM6L(I) + TERM6R(I)
      END DO

      ALPHA3 = .5*DT*DT - TAU*DT + TAU*TAU*(1.-EXTAU)      ! alpha 3 
      ALPHA2 = TAU*(-DT + ALPHA4) - ALPHA5                 ! alpha 2
 
      FS(1) = ALPHA1*WW0(1)*UF0(1)        + ALPHA2*TERM5(1) + ALPHA3*TERM4(1) + &
              TERM2(1)    -ALPHA7*TERM6(1)
      FS(2) = ALPHA1*WW0(1)*UF0(2)        + ALPHA2*TERM5(2) + ALPHA3*TERM4(2) + &
              TERM2(2)    -ALPHA7*TERM6(2)
      FS(3) = ALPHA1*WW0(1)*UF0(1)*VF0(1) + ALPHA2*TERM5(3) + ALPHA3*TERM4(3) + &
              TERM2(3)    -ALPHA7*TERM6(3)
!!$      FS(4) = ALPHA1*WW0(1)*UF0(1)*WF0(1) + ALPHA2*TERM5(4) + ALPHA3*TERM4(4) + &
!!$              TERM2(4)    -ALPHA7*TERM6(4)
      FS(4) = .5*ALPHA1*WW0(1)*(UF0(3)+UF0(1)*XVW02) + &
                ALPHA2*TERM5(5) + ALPHA3*TERM4(5) + TERM2(5) -ALPHA7*TERM6(5)

! ---------------------------------------------
!     PRANDTL NUMBER FIX
! ---------------------------------------------
!!$      DENSAV = AA(1) +AA(2)*UF0(1) +AA(3)*VF0(1) +  &
!!$                AA(4)*WF0(1) +  AA(5)*(UF0(2) + XVW02)
!!$      UMOMAV = FS(1)
!!$!      test = TERM4(1)
!!$      VMOMAV = AA(1)*VF0(1) + AA(2)*UF0(1)*VF0(1) + AA(3)*VF0(2) +  &
!!$               AA(4)*VF0(1)*WF0(1) +  &
!!$               AA(5)*(VF0(1)*(XF02 +UF0(2) +WF0(2)) + VF0(3) )
!!$      WMOMAV = AA(1)*WF0(1) + AA(2)*UF0(1)*WF0(1) + AA(3)*VF0(1)*WF0(1) +  &
!!$               AA(4)*WF0(2) +  &
!!$               AA(5)*(WF0(1)*(XF02 +UF0(2) +VF0(2)) + WF0(3) )
!!$      ENERAV = AA(1)*(UF0(2) + XVW02) +AA(2)*(UF0(3) + UF0(1)*XVW02) + &
!!$               AA(3)*(UF0(2)*VF0(1) + VF0(3) + VF0(1)*WF0(2) + &
!!$                        VF0(1)*XF02) + &
!!$               AA(4)*(UF0(2)*WF0(1) + WF0(3) + WF0(1)*VF0(2) + &
!!$                        WF0(1)*XF02) + &
!!$               AA(5)*(UF0(4) + VF0(4) + WF0(4) + XF04 + &
!!$                        2.*UF0(2)*XVW02 +  &
!!$                        2.*(VF0(2)*XF02 +WF0(2)*XF02 + VF0(2)*WF0(2)))
!!$      ENERAV = 0.5*ENERAV
!!$      
!!$      DENSAV = DT*(WW0(1) + .5*DT*DENSAV)
!!$!      test = DT*(WW0(2) + .5*DT*test)
!!$!      UMOMAV = TEST
!!$      VMOMAV = DT*(WW0(3) + .5*DT*VMOMAV)
!!$      WMOMAV = DT*(WW0(4) + .5*DT*WMOMAV)
!!$      ENERAV = DT*(WW0(5) + .5*DT*ENERAV)
!!$
!!$      UAV  = UMOMAV/DENSAV
!!$      VAV  = VMOMAV/DENSAV
!!$      WAV  = WMOMAV/DENSAV
!!$      EKIN = .5*(UAV*UAV + VAV*VAV + WAV*WAV)
!!$      Q = FS(5) -UAV*FS(2) -VAV*FS(3) -WAV*FS(4) +EKIN*FS(1)  &
!!$           -UAV*(ENERAV -UAV*UMOMAV -VAV*VMOMAV -WAV*WMOMAV +EKIN*DENSAV)
!!$      FS(5) = FS(5) + (1./PRN-1.)*Q


      FS(1) = FS(1)/DT
      FS(2) = FS(2)/DT
      FS(3) = FS(3)/DT
      FS(4) = FS(4)/DT
!!$      FS(5) = FS(5)/DT

      END SUBROUTINE BGKFLUX


! subroutine which returns the slopes of al and ar
      SUBROUTINE DXE3(U,V,W,LAMBDA,DW,A,GAM)
        IMPLICIT NONE
        REAL, INTENT(IN) :: U,V,W,LAMBDA
        REAL, INTENT(IN) :: GAM
        REAL, DIMENSION(5), INTENT(IN)  :: DW
        REAL, DIMENSION(5), INTENT(OUT) :: A
        REAL :: DD,CC,BB,AA

        DD=2.0*DW(5)-(U*U+V*V+W*W+1./(GAM-1.)/LAMBDA)*DW(1)
        CC=DW(4)-W*DW(1)
        BB=DW(3)-V*DW(1)
        AA=DW(2)-U*DW(1)
        A(5)=(GAM-1.)*LAMBDA*LAMBDA*(DD-2.0*U*AA-2.0*V*BB-2.0*W*CC)
        A(2)=2.0*LAMBDA*(AA-U*A(5)/LAMBDA)
        A(3)=2.0*LAMBDA*(BB-V*A(5)/LAMBDA)
        A(4)=2.0*LAMBDA*(CC-W*A(5)/LAMBDA)
        A(1)=DW(1)-A(2)*U-A(3)*V-A(4)*W-A(5)*(U*U+V*V+W*W+1./(GAM-1.)/LAMBDA)
        RETURN
      END SUBROUTINE DXE3
