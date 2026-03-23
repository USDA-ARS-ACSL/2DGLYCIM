!*************************(COMMON:Common Block)*************************
!**                                                                   **
!**                   GLYCIM with 2D soil Nov 09, 2020                **
!**                                                                   **
!***********************************************************************
      LOGICAL  INIT
      INTEGER Modnum
      REAL Period, NitrogenUptake
      COMMON /ABC/ Period, NitrogenUptake, Modnum, INIT
     
      REAL SEEDLB, FILL, IDET, PARM
      COMMON /VARIETY/ SEEDLB, FILL, IDET, PARM(25)

      REAL  CLUE, LLUE, LTC, SINKC,                                    &
      TAUM, TAUNIN, TAUN, WATWK, CCCT, PGLIMT, PGCOR,                  &
      WATLF, WWKLF, CCTLF, PGLIM, SEN, REDSEN, BMAIN,                  &
      PINDEX, PNLLFS, SUMESO, SUMRAIN    
      INTEGER ILOWP, ILOW, ILOWB, DRPDAY, DROP
      
      COMMON /PHN/ CLUE, LLUE, LTC,                                    &
      SINKC, TAUM, TAUNIN, TAUN, WATWK, CCCT, PGLIMT, PGCOR,           &
      WATLF, WWKLF, CCTLF, PGLIM, SEN, REDSEN, BMAIN,                  &
      PINDEX, ILOWP, ILOW, ILOWB(20), PNLLFS, DRPDAY, DROP,            &
      SUMESO, SUMRAIN
      
      REAL CBAL1, SCPOOL, USEDCS, XTRACS, FIXC, FIXCS,                 &
      TRANSC, SCPMAX
      COMMON /CARBON/ CBAL1, SCPOOL, USEDCS, XTRACS, FIXC, FIXCS,      &
      TRANSC, SCPMAX

      REAL SHTWT, ULEAFW, ULAREA, MLEAFW, MLAREA,                      &
      INITBR, BLAREA, BLEAFW
      INTEGER NBRNCH
      COMMON /GROWTH/ SHTWT, ULEAFW, ULAREA, MLEAFW(35), MLAREA(35),   &
      NBRNCH, INITBR(20), BLAREA(20, 20), BLEAFW(20, 20)
   
      INTEGER NTOPLF, MATURE, IPERD
      REAL RSTAGE, VSTAGE, TB, PDR
      COMMON /PHENOLOGY/ RSTAGE, VSTAGE, TB(20), NTOPLF, PDR, MATURE,  &
      IPERD

      INTEGER MM2, DD, YY
      CHARACTER*10  CROP_DATE1, CROP_DATE2
      CHARACTER*10  EDATE, R1DATE, R4DATE, R8DATE, R5DATE, R3DATE,     &
                    R7DATE

      COMMON / D_VARS / MM2, DD, YY, CROP_DATE1, CROP_DATE2,           &
           EDATE, R1DATE, R4DATE, R8DATE, R5DATE, R3DATE, R7DATE

      INTEGER INITR0, LIMIT, FILDAY, NFLWRS, NSEEDS
      REAL TBMAH, TNL, MPETL, BPETL, SDWMAX, RSTPC, BRNCHH
      REAL CONVSE, CONVPO
      COMMON /PGROWTH/ INITR0, LIMIT, TBMAH(20), TNL, FILDAY,          &
       PDV, PDPLM(3), PDLAM(3), PDPWM(3), MPETL(35), PDLAB(20, 20),    &
       PDPLB(20, 20), PDPWB(20, 20), BPETL(20, 20), PDMH,              &
       PDBRW(20), PCST, QCST, NFLWRS, PDFW, PDPODW, NPODS, PDPODC,     &
       NSEEDS, PDSW, PDSC, SDWMAX, RSTPC, BRNCHH(20),CONVSE, CONVPO
     
      REAL SCRTS1, NODC, FRUITC, SCRATO, EQSC, PRENOD
      COMMON /PCarbon/ SCRTS1, NODC, FRUITC, PRENOD, SCRATO,           &
             EQSC, SPCSP, PODCD, SVEGC, SPCSS, SEEDCD

      INTEGER INITR2
      REAL SUPNO3, CNO3
      REAL REDN, PLANTN, NODN, PLNH4, ABSLPN, DEADRN, ABSPON,          &
           PDNODC

      COMMON / NITROGEN / INITR2, PLNO3, SUPNO3, CNO3, REDN,           &
                  PLANTN, NODN, PLNH4, ABSLPN, DEADRN, ABSPON, PDNODC
    
      REAL PODWT, SEEDWT, FLWRWT, STEMWT, PETWT, ROOTWT, LEAFWT,YLD
      REAL SCRTS2, RSTPCS, RCPOOL, VSTMAH, SDWT, PARTRT, SGTR
      REAL UPETL, UPETW, MPETW, LAREAM, BPETW, LAREAB, MSTEMH, MSTEMW
      INTEGER INITR1, INITR3, INITR4, MBRNCH, NODES,FPODS, FSEEDS
      REAL BRNCHW, MG
      CHARACTER  LIMITF*1
      COMMON /AGROWTH/ PODWT, SEEDWT, FLWRWT, STEMWT, PETWT, ROOTWT,   &
            LEAFWT, SCRTS2, RSTPCS, RCPOOL, VSTMAH,                    &
            INITR1, INITR3, INITR4, MBRNCH, NODES,                     &
            PDWT(1000), PPDWT(1000), SDWT(2000), PSDWT(2000),          &
            YLD, PARTRT, SGTR, UPETL, UPETW, MPETW(35), LAREAM,        &
            BPETW(20, 20), LAREAB, MG, MSTEMH, MSTEMW, BRNCHW(20),     &
            FPODS, FSEEDS, LIMITF(24)

      REAL NCPOOL, POPSLB, NODNUM, NFIXNU, NFIXWT, NODIE, NNEWWT,      &
            NODNEW, DAYBR, RESPS, NFCF, PDNW
      COMMON /NODEGROW/ NCPOOL, POPSLB, NODNUM(4000), NFIXNU(4000),    &
             NFIXWT(4000), NODIE(4000), NNEWWT(10, 4000),              &
             NODNEW(10, 4000), DAYBR(20), RESPS(4000), NFCF(4000),     &
             PDNW(4000), CXT(4000)
      
      REAL NODWT, CONVN,PDMW, RNODC, CNFIXS, CMAXSG, CASG
      Integer NCR, NODWAR
      COMMON /sort1/ NCR, NODWAR(4000), NODWT, CONVN, PDMW, RNODC,     &
       CNFIXS, CMAXSG, CASG

      INTEGER  ABSCIS
      COMMON /LOST/ ABSCIS
      
      REAL Tgrowth, HourlyNitrogenDemand
      COMMON / MAINF / Tgrowth, HourlyNitrogenDemand
     
      INTEGER LONDAY
      REAL NRATS, ALLWT, TOPWT, LFWT, PW, ABSDW, PODWT6, RSTAGEold
      COMMON / tislost / NPA, NRATS, LONDAY, ALLWT, TOPWT, LFWT, PW,   &
           ABSDW, PODWT6, RSTAGEold

      REAL PFD, POTET, ACTET, WaterUptake,                             &
          Parave, Sradave, Tdayave, Tcanave,                           &
           gsmax, psilave, daypotet, dayactet,                         &
          wstressave, nstressave, NSTRESS, WSTRESS,                    &
          Cstress, PSNET, PGROSS, netsunlitleaf, AAgs, Avtemperature,  &
          howtime, whattime, DailyBio, Tleafmax, Biomass1
      Common / EVAPOR / PFD, POTET, ACTET, WaterUptake,                &
          Parave, Sradave, Tdayave, Tcanave,                           &
           gsmax, psilave, daypotet, dayactet,                         &
          wstressave, nstressave, NSTRESS, WSTRESS,                    &
          Cstress, PSNET, PGROSS, netsunlitleaf, AAgs, Avtemperature,  &
          howtime, whattime, DailyBio, Tleafmax, Biomass1

      REAL DDT, EDD, ADDT, AEDD, TRAIN, PPET,                          &
          wstressday, TDDT, TEDD, JJARAIN, TPET, TAET, JJAPET, JJAAET, &
          JJAVPD, TVPD
      INTEGER R1DAY, R4DAY, R5DAY, R8DAY
      Common / REGION / DDT, EDD, ADDT, AEDD, TRAIN, PPET,             &
          wstressday, TDDT, TEDD, JJARAIN, TPET, TAET, JJAPET, JJAAET, &
          JJAVPD, TVPD, R1DAY, R4DAY, R5DAY, R8DAY
      REAL  LeafNitrogen, StemNitrogen, PodNitrogen, SeedNitrogen,     &
            RootNitrogen, DeadNitrogen, PlantNitrogen, NCRatio

      Common /cropnitrogen/ LeafNitrogen, StemNitrogen, PodNitrogen,   &
          SeedNitrogen, RootNitrogen, DeadNitrogen, PlantNitrogen,     &
          NCRatio
      
      INTEGER CDayOfYear, CITIME, CIPERD
      REAL CWATTSM, CPAR, CTAIR, CCO2, CVPD, CWIND, CPSIL_,            &
           CLATUDE, CLAREAT, CLAI

      COMMON / Weather / CDayOfYear, CITIME, CIPERD, CWATTSM(24),      &
      CPAR(24), CTAIR(24), CCO2, CVPD(24), CWIND, CPSIL_,              &
      CLATUDE, CLAREAT, CLAI

      REAL NRATIO, PGR, PN, transpiration, temperature,  TLAI,         &
      SunlitLAI, ShadedLAI, LightIC, transpiration_sunlitleaf,         &
      transpiration_shadedleaf, temp1, Ags, ARH,                       &
      photosynthesis_netsunlitleaf, photosynthesis_netshadedleaf

      COMMON / Plant / NRATIO, PGR, PN, transpiration, temperature,    &
      TLAI, SunlitLAI, ShadedLAI, LightIC,                             &
      transpiration_sunlitleaf, transpiration_shadedleaf,              & 
      temp1, Ags, ARH,                                                 &
      photosynthesis_netsunlitleaf, photosynthesis_netshadedleaf