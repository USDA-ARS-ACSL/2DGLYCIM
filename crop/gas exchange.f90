!**********************************( PNET )****************************
!**                                                                  **
!**   This module uses single-leaf photosynethic characteristics to  **
!**   calculate crop canopy characterisitics and canopy gross        **
!**   photosynthetic rate. Photorespiration rate and maintenance rate**
!**   are calculated and subtracted to give net photosynthetic rate, **
!**   which is corrected for stomatal conductance.                   **
!**   From this are calculated the net carbon fixation rate          **
!**   and the rate of carbon translocation out of the leaves. Net    **
!**   photosynthesis is also calculated for the lowest layer of      **
!**   leaves in the canopy to see if they are "parasitic" on the     **
!**   plant.                                                         **
!**   Calculates the leaf area per unit ground area covered by the   **
!**   crop canopy (LCAI). If the current time period is the          **
!**   first of the day, canopy light utilization efficiency (CLUE)   **
!**   and canopy conductance to CO2 transfer (CCCT) are calculated.  **
!**   These are photosynthetic parameters.                           **
!**   During daylight periods,gross photosynthesis as limited by     **
!**   temperature, leaf nitrogen content and leaf age is calculated. **
!**   Photorespiration (LYTRES), and maintenance respiration (BMAIN) **
!**   are subtracted to leave net carbon fixation (FIXC).  From this **
!**   and the size of the shoot carbon pool (SCPOOL) potential carbon**
!**   translocation (TRANSC) is calculated.  The amount of carbon    **
!**   used by the plant in photorespiration and maintenance          **
!**   respiration is added to the sum of carbon used by the plant    **
!**   (USEDCS). If the vegetative stage of the plant is such that    **
!**   upper leaves may shade lower leaves to a significant degree,   **
!**   net carbon exchange rate is calculated for the leaves at the   **
!**   bottom of the canopy to see if they are "self-supporting" or   **
!**   "parasitic".                                                   **
!**********************************************************************
      SUBROUTINE PNET
	  
	  use common_block
      INCLUDE 'Common.h'
	 
      REAL    LLAREA, LLFWT, LYTRES
	  !WSUN change endday to dayend
      LOGICAL LIGHT, DAYEND       
      COMMON /PNET1/ D11, LLAREA, LLFWT
      
      CBAL1 = SCPOOL + USEDCS + XTRACS

! LEAF AREA PER UNIT GROUND AREA COVERED BY CROP CANOPY

      LCAI = LAREAT*POPROW/100.0/MIN(HEIGHT,ROWSP)

! WSUN calculate WATWK here (original weather section)	  
! AVERAGE DAYTIME LIGHT FLUX DENSITY ON CROP FOR PAST WEEK
      IF ( INIT ) THEN
         WATWK = RI*1.25E-4/DAYLNG
      ELSE
         WATWK = WATWK*6.0/7.0 + RI*1.7857E-5/DAYLNG
      END IF
! 1.25E-4 = 0.45/3600;    1.7857E-5 = 0.45/3600/7	  

! Note 7.01
! CANOPY PHOTOSYNTHETIC PARAMETERS ARE CALCULATED ONCE A DAY
! FROM LEAF PARAMETERS WEIGHTED AND INTEGRATED OVER THE CANOPY.
! ------------------------------------------------------------------
! CALCULATE CANOPY LIGHT UTILIZATION EFFICIENCY

      IF (ITIME.EQ.IDAWN) THEN
         D11 = EXP(-LCAI*CEC)
         CLUE = LLUE*(1.0 - D11)/(1.0 - LTC)

! CALCULATE CANOPY CONDUCTANCE TO CO2 TRANSFER
! Note 7.02-7.06
!        TAUM = 8.5E-5
         TAUM = 17.0E-05                         !YA
         TAUNIN = MAX(0.045 - (RSTAGE*0.01), 1.5E-2)
!wsun Change SINK TO SINKC
!Indicator of source/sink balance for carbon 
         TAUN = MAX(TAUN + (SINKC*2.5E-3), TAUNIN)
         TAUN = MIN(TAUN, 5.5E-2)
         D12 = TAUN*LLUE*CEC*WATWK
         D13 = LLUE*(1.0 - LTC)
         CCCT = ALOG((D12 + D13)/(D12*D11 + D13))*TAUM/TAUN/CEC

! OPEN THE STOMATA FULLY AT DAWN.

         SCF = 1.0
	  END IF 
      ! pass weather and plant common blocks to GasExchanger
	  Call gasexchanger(CDayofYear,NRATIO)
	  
      LIGHT = WATTSM(ITIME).GT.0.0
      IF (LIGHT) THEN
! ----------------------------------------------------------------------
! CANOPY GROSS PHOTOSYNTHETIC RATE AND ITS REDUCTION BY LOW TEMPERATURES,
! LOW LEAF NITROGEN CONTENT AND LEAF SENESCENCE ARE CALCULATED.


		IF (NRATIO.LE.0.4) THEN
            REDN = 0.
         ELSE
! Note 7.12
            REDN = ((NRATIO*8.25) - 3.3)/((NRATIO*6.6) - 1.64)
         ENDIF

! CALCULATE A PHOTOSYNTHESIS MULTIPLICATION FACTOR BASED ON
! REPRODUCTIVE STAGE.

         IF (RSTAGE.GT.6.0) THEN
! Note 7.13
!VRR        SEN = MAX (0.0, 1.0 - (RSTAGE - 6.0)*0.5)
            SEN = (9.0-RSTAGE)/3.0
         ELSE
            SEN = 1.0
		 ENDIF

! PHOTOSYNTHESIS IS EITHER LIMITED BY LEAF NITROGEN CONTENT
! OR BY SENESCENCE.                                               *** Note 7.14
         REDSEN = MIN(SEN, REDN)
! ------------------------------------------------------------------
! LIGHT RESPIRATION RATES AND MAINTENANCE RESPIRATION RATES ARE
! SUBTRACTED FROM GROSS PHOTOSYNTHETIC RATE. THE CARBON FIXED IS
! ALLOCATED TO STORAGE IN THE SHOOT OR IS TRANSLOCATED OUT OF
! THE LEAVES TO GROWING ORGANS.
! ------------------------------------------------------------------
! CALCULATE MAINTENANCE RESPIRATION RATE                          *** Note 7.16
         BMAIN = RESPM(LEAFWT,NRATIO,POPAREA,TAIR(ITIME))
! 
! CALCULATE NET PHOTOSYNTHETIC RATE 

		 PINDEX=PGR*REDSEN- (PGR-PN)-BMAIN

! CALCULATE NET PHOTOSYNTHETIC RATE ALLOWING FOR STOMATAL CLOSURE
! Note 7.17
!		 PN = MAX(PINDEX*SCF, 0.0)
         PN = MAX(PINDEX, 0.0) !mg CO2 GROUND m-2s-
! CALCULATE NET CARBON FIXATION RATE
! Note 7.18
         FIXC = PN*0.9818/POPAREA !g C plant h-
! 0.9818 = 0.27273*3600/1000
         FIXCS = FIXCS + FIXC*PERIOD
! FIXC   Carbon fixation rate allowing for stomatal closure (g plant-1 hr-1).
! FIXCS  Accumulated net carbon fixed (g plant-1).
! CALCULATE POTENTIAL RATE OF CARBON TRANSLOCATION FROM THE
! SHOOT TO GROWING ORGANS AND AMOUNT STORED IN SHOOT:
! The shoot's carbon pool may not exceed a fixed percentage of
! of the shoot's weight.  Keep track of how much carbon is "discarded."
! Note 7.19 - 7.20
         TRANSC = (FIXC*0.35) + (SCPOOL*(1.0 - (0.875**PERIOD))/PERIOD)
         SCPOOL = SCPOOL + ((FIXC-TRANSC)*PERIOD)
         SCPMAX = 0.4*SHTWT
         IF (SCPOOL.GT.SCPMAX) THEN
            XTRACS = XTRACS + SCPOOL - SCPMAX
            SCPOOL = SCPMAX
		 END IF
! XTRAC  Carbon lost because it cannot be used or stored (g plant-1 h-1).
! SCPOOL Carbon stored in the shoot (g plant-1).	 
! SCPMAX Maximum amount of carbon that can be stored in the shoot (g plant-1).	 
! TRANSC Potential rate of carbon translocation from the shoot to growing organs (g plant-1 hr-1).
		 
! AT NIGHT NO CARBON IS FIXED BUT MAINTENANCE RESPIRATION
! CONTINUES AND CARBON IS TRANSLOCATED OUT OF STORAGE.
      ELSE
!		 PGR = 0.0
!		 PGROSS = 0.0
         FIXC = 0.0
         BMAIN = RESPM(LEAFWT,NRATIO,POPAREA,TAIR(ITIME))
         TRANSC = (FIXC*0.35) + (SCPOOL*(1.0 - (0.875**PERIOD))/PERIOD)
         BMAINN = BMAIN*0.9818/POPAREA
! 0.9818 = 0.27273*3600/1000
         USEDCS = USEDCS + (MIN(TRANSC, BMAINN) * PERIOD)
         SCPOOL = SCPOOL - (TRANSC*PERIOD)
         TRANSC = MAX(TRANSC - BMAINN, 0.0)
      END IF
      
      
! CALCULATE FOR LEAVES AT THE BOTTOM OF THE CANOPY:
!   1) GROSS PHOTOSYNTHETIC RATE (PGLF)
!   2) PHOTORESPIRATION RATE (LYTRES)
!   3) WEIGHT AND AREA OF THE LOWEST LAYER OF LEAVES (LLFWT, LLAREA)
!   4) NET PHOTOSYNTHETIC RATE ACCUMULATED OVER THE DAY (PNLLFS)
! AND SET A FLAG (DROP) IF THE LOWEST LEAVES ARE "PARASITIC".
! THE LOGICAL ARGUEMENT LIGHT IS SET TRUE IF THERE IS LIGHT IN THE
! CURRENT TIME PERIOD.
! ------------------------------------------------------------------
! NET PHOTOSYNTHETIC RATE FOR LEAVES AT THE BOTTOM OF THE
! CANOPY IS INTEGRATED OVER 24 HOURS TO SEE IF THEY ARE STILL
! SELF-SUFFICIENT.
! ------------------------------------------------------------------
      IF (VSTAGE.GE.4.0) THEN
         DAYEND = ITIME.EQ.IPERD

! CALCULATE GROSS PHOTOSYNTHETIC RATE FOR LEAVES AT THE BOTTOM OF THE CANOPY

         IF (LIGHT) THEN                                        ! *** Note 7.21
            WATLF = PAR(ITIME)*CEC/(1.0 - LTC)*EXP(-CEC*LCAI)
            WWKLF = WATLF*WATWK/PAR(ITIME)
            CCTLF = TAUM*WWKLF/(1.0 + TAUN*WWKLF)
            PGLF = LLUE*WATLF*CCTLF*CO2*539.66/(273.15 + TAIR(ITIME))/ &
                  (LLUE*WATLF + CCTLF*CO2*539.66/(273.15 + TAIR(ITIME)))

! CALCULATE RESPIRATION RATES OF LEAVES AT THE BOTTOM OF THE CANOPY

            LYTRES = PGLF*REDSEN*151890.0/CO2*0.00012*                 &
			         EXP(0.0295*TAIR(ITIME))
         END IF

! CALCULATE THE WEIGHT AND AREA OF THE LOWEST LAYER OF LEAVES.

         IF (ILOWP .EQ. ILOW) THEN
            CONTINUE
         ELSE IF (ILOW.EQ.0) THEN
            LLFWT = ULEAFW
            LLAREA = ULAREA
         ELSE
            LLFWT = MLEAFW(ILOW)
            LLAREA = MLAREA(ILOW)
            IF (ILOW .NE. 1) THEN
               J = 1
			   DO WHILE (J .LE. NBRNCH )
                  IF (INITBR(J).GE.(ILOW)) EXIT
                  I = ILOWB(J)
                  LLFWT = LLFWT + BLEAFW(I,J)
                  LLAREA = LLAREA + BLAREA(I,J)
                  J = J + 1
			   END DO
            END IF
         END IF
         ILOWP = ILOW
         IF (.NOT.((LLAREA.GT.0.0).OR.(LLAREA.LT.0.0))) RETURN

         D15 = LLFWT/LLAREA*10000.
         D16 = 1.0
         BMAIN = RESPM(D15,NRATIO,D16,TAIR(ITIME))

! CALCULATE ACCUMULATED NET PHOTOSYNTHETIC RATE FOR LEAVES AT
! THE BOTTOM OF THE CANOPY
        
         IF (LIGHT) THEN                                        ! *** Note 7.22
            PNLLFS = PNLLFS + ((PGLF*REDSEN) - LYTRES -BMAIN)*SCF
         ELSE
            PNLLFS = PNLLFS - BMAIN
         END IF

! TEST TO SEE IF LOWEST LEAVES IN CANOPY ARE SELF SUPPORTING.
! IF NOT DROP THEM

         IF (DAYEND) THEN                                       ! *** Note 7.23
            IF (PNLLFS.LT.0.1) THEN
               DRPDAY = DRPDAY + 1
               IF (DRPDAY.GE.2) DROP = 1
            ELSE
               DRPDAY = 0
               DROP = 0
            END IF
            PNLLFS = 0.0
         END IF
      END IF

       
      END
!***********************************************************************
! ### Function RESPM calculates maintenance respiration rate ###

      REAL FUNCTION RESPM(LEAFWT,NRATIO,POPAREA,TAIRL)
      REAL LEAFWT, NRATIO, POPAREA, TAIRL
      RESPM = POPAREA*0.25*EXP(0.0693*TAIRL)*                           &
	          (11.733*LEAFWT + 366.67*LEAFWT*NRATIO*0.05) / 86400.
      RETURN
      END