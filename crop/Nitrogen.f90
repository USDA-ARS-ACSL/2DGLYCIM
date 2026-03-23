!* * * * * * * * * * * * * * * * *  N I T R O  * * * * * * * * * * * * *
! Calculates the supply and demand for nitrogen in the whole plant and
! determines a nitrogen stress factor which is used to partition carbon
! to nodules.
!***********************************************************************
  SUBROUTINE NITRO
      use common_block
      INCLUDE 'common.h'
	  
	 
      IF ( INIT ) THEN
      INITR2 = 0
      PLNO3 = 0.0
	  PLANTN = 0.0065
	  END IF

! CALCULATE UPTAKE OF NO3 IN THE TRANSPIRATION STREAM.
! wsun Transfer VH2OC, water Nitrogen content, AWUP FROM 2DSOIL  
! VH2OC(NumNPD) volumetric water content cm3 cm-3 
! NO3_old(NumNPD) SOIL NITROGEN CONTENT
! CALCULATE MEAN ROOT GROWTH RATE AND MEAN WATER UPTAKE RATE  FOR THE PERIOD.    
! AWUP(NumELD)                                            *** Note 10.35
!     UPTH2O(L,K) = AWUP(L,K)*PERIOD
!	  SUPNO3 = 0.0
!      DO L = 1,NLR
!        DO K = 1,NKR
!            D11 = UPTH2O(L,K)*VNO3C(L,K)/VH2OC(L,K)
!            SUPNO3 = SUPNO3 + (D11/1000.0/POPSLB)
!         END DO
!      END DO
! --------------------------------------------------------------------
! THE RATIO OF THE ACTUAL NITROGEN CONTENT OF THE PLANT TO ITS maximum
! POSSIBLE NITROGEN CONTENT IS USED TO CALCULATE A NITROGEN STRESS FACTOR.
! ------------------------------------------------------------------
! CARBON ALLOCATED TO NODULE GROWTH AND ACTIVITY IS USED TO
! REDUCE ANY NITRATE NITROGEN PRESENT IN THE PLANT
! SIncrSink g slab-1 h-1
! Total nitrate-nitrogen in plant (g N).
	  PLNO3 = PLNO3 + NitrogenUptake/POPSLB ! THIS FROM 2DSOIL g plant-1
!      write(*,*) 'Time: ', time, 'PLNO3: ',PLNO3, 'NitrogenUptake: ',NitrogenUptake
      PLNH4 = 0.0
      IF (PLNO3.GT.0.0) THEN
! Note 11.01
         CNO3 = PLNO3*2.0/PERIOD     ! Rate of carbon use in nitrate reduction (g C plant-1 hr-1).
         IF (CNO3.GT.NODC) THEN
            USEDCS = USEDCS + (NODC*PERIOD)
!YAP        PLNH4 = NODC/2.0*PERIOD
            PLNH4 = PLNO3/CNO3*NODC
!YAPEND
            PLNO3 = PLNO3 - PLNH4
            NODC = 0.0
         ELSE
            USEDCS = USEDCS + (CNO3*PERIOD)
            PLNH4 = PLNO3
            PLNO3 = 0.0
            NODC = NODC - CNO3
         END IF
      END IF

! CALCULATE TOTAL NITROGEN IN REDUCED FORM IN THE PLANT.

       PLANTN = PLANTN + NODN + PLNH4 - ABSLPN - DEADRN - ABSPON
!	   PLANTN Total nitrogen in reduced form in plant (g N).!This is from 2dsoil
! RESET ACCUMULATORS

      NODN = 0.0
      PLNH4 = 0.0
      ABSLPN = 0.0
      DEADRN = 0.0
      ABSPON = 0.0

! CALCULATE NITROGEN REQUIRED BY PODS AND SEEDS

      PODN = PODWT*0.04
      SEEDN = SEEDWT*0.065
      D12 = PODN + SEEDN

! CALCULATE TOTAL NITROGEN REQUIRED BY PLANT AT MAXIMUM
! NITROGEN CONTENT                                                *** Note 11.02
      REQN = ((((FLWRWT + LEAFWT)*5.0) + ((STEMWT + PETWT)*2.0)        &
           + (ROOTWT*2.5))/100.) + PODN + SEEDN
! Total nitrogen required by plant at maximum nitrogen content (g).
	  Tgrowth = Tgrowth + Period                !To calculate how many hours for soybean growth
      HourlyNitrogenDemand =  REQN/Tgrowth

! CALCULATE NITROGEN SUPPLY/DEMAND FOR VEGETATIVE PARTS           *** Note 11.03
    NRATIO = (PLANTN - D12)/(REQN - D12)

      IF (NRATIO.LE.0.1) NRATIO = 0.1

! CALCULATE POTENTIAL CARBON SUPPLY TO NODULES FOR NEXT ITERATION *** Note 11.04
      PDNODC = ((REQN - PLANTN)*2.0/PERIOD)*0.25
      IF (PDNODC.LE.0.0) PDNODC = 0.0

      RETURN
      END