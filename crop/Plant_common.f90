!*************************(COMMON:Common Block)*************************
!**                                                                   **
!**                  2D SOIL GLYCIM october 20, 2020                  **
!**                                                                   **
!***********************************************************************
module common_block 

!DEC$ ATTRIBUTES DLLEXPORT:: / ShootR /
!DEC$ ATTRIBUTES DLLEXPORT:: / shtR_public /
!DEC$ ATTRIBUTES DLLEXPORT:: / Weath /
!DEC$ ATTRIBUTES DLLEXPORT:: / grid_public /
!DEC$ ATTRIBUTES DLLEXPORT:: / nodal_public /
!DEC$ ATTRIBUTES DLLEXPORT:: / elem_public /
!DEC$ ATTRIBUTES DLLEXPORT:: / bound_public /
!DEC$ ATTRIBUTES DLLEXPORT:: / time_public /
!DEC$ ATTRIBUTES DLLEXPORT:: / module_public /
!DEC$ ATTRIBUTES DLLEXPORT:: / DataFilenames /
!DEC$ ATTRIBUTES DLLEXPORT:: / Materials /

Parameter(NumNPD = 4000, NumElD = 3500, NumBPD = 600, NSeepD = 2,      &
    NumSPD = 30, NumSD = 10, NDrainD = 2, NumDR = 30,                  &
    NumGD = 3, NumPlD = 100,                                           &
    NMatD = 15, MNorth = 4,                                            &
    NumModD = 20, MBandD = 15, NumSurfDatD = 3 + NumGD + NumSD) 


    REAL PCRL, PCRQ, PCRS, ET_demand, LCAI, COVER,  CONVR,             &
         MaxRootDepth,SHADE, HEIGHT, LAI, AWUPS,                       &
         NitroDemand,xBSTEM, yBSTEM, SGT, PSIM,                        &
		  LAREAT, POPROW,ROWSP,ROWANG,PopArea, CEC,EORSCS,             &
          AWUPSS, SOLRAD, Total_Eor, Total_Pcrs, sincrsink, PSILD,     &
          OSMFAC, EOMULT, PSIL_,  NDemandError,                         &
               CumulativeNDemandError,                                 &
               InitialRootCarbo, PCRTS,                                &
               ConstI(2),constK(2), Cmin0(2) 
    
          
    Real(8) HourlyCarboUsed, TotalRootWeight
    integer isGerminated, isEmerged
    REAL cContentRootY,nContentRootY,cContentRootM,nContentRootM

    Common / ShootR / PCRL, PCRQ, PCRS, HourlyCarboUsed, ET_demand,     &
               LCAI, COVER, CONVR,                                      &
               MaxRootDepth, SHADE, HEIGHT, LAI, AWUPS,                 &
               NitroDemand,   xBSTEM, yBSTEM, SGT, PSIM,                &
               LAREAT, POPROW, ROWSP, ROWANG, PopArea,                  &
               CEC, EORSCS, AWUPSS, SOLRAD,                             &
               Total_Eor, Total_Pcrs, sincrsink, PSILD,                 &
               OSMFAC, EOMULT, PSIL_, NDemandError,                     &
               CumulativeNDemandError,TotalRootWeight,                  &
               InitialRootCarbo, PCRTS,                                 &
               ConstI,constK, Cmin0,                                    &
               isGerminated, isEmerged
               
    REAL RVR, AWUP, PSISM, PSILT, EOR, YRL, RGCF,                       &
               PSIS, TS, COND, PDWR, DRL, WUPM, WUPN, WUPT, TPL, PPDRL, &
        ADWR, AWR, ADRL, RTWL, ALPY,                                    &
      WidthE, HeightE, PILD, VMAX, vegsrc, potential_T, potential_T_EO, &
            AWUPS_old, WUPMS, EORSCF, WUPDS, WUP2S, WUP0S, DPSI02, SCF, &
            PSIST, PROPAR, VUP, FUP

    Common / shtR_public / RVR(NumNPD), AWUP(NumELD), PSISM, PSILT, EOR,&
               YRL(NumNPD), RGCF(NumNPD), iFavRoot(NumNPD),             &
               PSIS(NumNPD),                                            &
               TS(NumNPD), COND(NumNPD), PDWR(NumNPD),                  &
               DRL(NumNPD), WUPM(NumNPD),                               &
               WUPN(NumNPD), WUPT(NumNPD), TPL, PPDRL(NumNPD),          &
               ADWR(NumNPD), AWR(NumNPD), ADRL(NumNPD), RTWL,           &
               ALPY,                                                    &
               WidthE(NumEld), HeightE(NumEld),                         &
               PILD, VMAX, vegsrc, potential_T,                         &
               potential_T_EO,                                          &
               AWUPS_old, WUPMS, EORSCF, WUPDS, WUP2S,                  &
               WUP0S, DPSI02, SCF, PSIST, PROPAR,                       &
               VUP(NumNPD, 2), FUP(NumNPD, 2),                          &
               cContentRootM,cContentRootY,                             &
               nContentRootM, nContentRootY

    Integer MSW1, MSW2, MSW3, MSW4, MSW5, MSW6, MSW7
    REAL BSOLAR, ETCORR, BTEMP, ATEMP, ERAIN, BWIND, BIR, WINDA, IRAV
    Integer JDAY, NCD, JDLAST
    REAL CLDFAC, DEL, RINT, RNS, RNC, RAIN, IR
    REAL WIND, CO2, TDUSK, TDUSKY, CPREC, TAIR, VPD
    REAL ROUGH, RADINT, WATTSM, DIFINT
    REAL ROWINC, CLOUD, SHADOW, DIFWAT
    REAL DIRINT, WATACT, WATRAT, WATPOT, RNLU
    Integer NumF, NumFP
    REAL hFur, QF
    Integer IFUR
    REAL GAIR, PG, LATUDE, Longitude, Altitude, RI, par, parint, DAYLNG
    REAL DLNGMAX, AutoIrrigAmt, ESO, TMIN
    Integer AutoIrrigateF

    Common / Weath / MSW1, MSW2, MSW3, MSW4, MSW5, MSW6,               &
    MSW7, BSOLAR, ETCORR,                                              &
    BTEMP, ATEMP, ERAIN, BWIND, BIR, WINDA, IRAV, JDAY,                &
    NCD, JDLAST, CLDFAC, DEL(24), RINT(24), RNS,                       &
    RNC, RAIN, IR, WIND, CO2, TDUSK, TDUSKY,                           &
    CPREC(NumSD), TAIR(24), VPD(24), ROUGH,                            &
    RADINT(24), WATTSM(24), DIFINT(24),                                &
    ROWINC(24), CLOUD, SHADOW(24), DIFWAT(24),                         &
    DIRINT(24), WATACT, WATRAT, WATPOT, RNLU,                          &
    NumF(40), NumFP, hFur(40), QF, IFUR, GAIR(NumGD), PG,              &
    LATUDE, Longitude, Altitude, RI, PAR(24),                          &
    PARINT(24), daylng, DLNGMAX, AutoIrrigAmt,                         &
    AutoIrrigateF, ESO, TMIN

    Integer NumNP, NumEl, IJ, KAT, MBand, Nmat, KX
    Real x, y, Area, nodeArea

    Common / grid_public / NumNP, NumEl, IJ, KAT, MBand, Nmat,         &
    KX(NumElD, 4), x(NumNPD), y(NumNPD), Area(NumElD),                 &
    nodeArea(NumNPD)

    Integer NumSol, NumG, ListN, ListNE, MatNumN
    REAL hNew, ThNew, Vx, Vz, Q, Conc, g, Tmpr, Con, TcsXX, RO,        &
    hNew_org, QAct, ThetaAvail, ThetaFullRZ, ThAvail, ThFull,            &
    QGas, ThetaAir
    
    Logical*1 lOrt

    Common /nodal_public/ NumSol, NumG, ListN(NumNPD), ListNE(NumNPD), &
    MatNumN(NumNPD), hNew(NumNPD), ThNew(NumNPD), Vx(NumNPD),          &
    Vz(NumNPD), Q(NumNPD), Conc(NumNPD, NumSD),                        &
    g(NumNPD, NumGD), Tmpr(NumNPD), Con(NumNPD), TcsXX(NumNPD),        &
    RO(NumNPD), hNew_org(NumNPD), QAct(NumNPD), ThetaAvail,            &
    ThetaFullRZ, ThAvail(NumNPD), ThFull(NMatD), lOrt,                 &
    QGas(NumNPD, NumGD), ThetaAir(NumNPD)

    Integer MatNumE
	REAL  cSink, gSink               
    REAL  tSink, RTWT, RMassM, RDenM,                                  &
    RMassY, RDenY, gSink_OM, gSink_root, gsink_N2O,                    &
    gSink_rootY, gSink_rootM, FracClay, cSink_OM

    Common / elem_public / MatNumE(NumElD), Sink(NumNPD),              &  
      cSink(NumNPD, NumSD), gSink(NumNPD, NumGD), tSink(NumNPD),       &  
      RTWT(NumNPD), RMassM(NumNPD), RDenM(NumNPD),     & 
      RMassY(NumNPD), RDenY(NumNPD),                                   &
      gSink_OM(NumNPD, NumGD), cSink_OM(NumNPD,NumSD),                 &
      gSink_root(NumNPD, NumGD), gSink_rootY(NumNPD, NumGD),           &
      gSink_rootM(NumNPD, NumGD), gSink_N2O(NumNPD, NumGD)

    Integer NumBP, NSurf, NVarBW, NVarBS, NVarBT, NVarBG,              &
    NumSurfDat, NSeep, NSP, NP,                                        &
    NDrain, NDR, ND, KXB
    Integer CodeW, CodeS, CodeT, CodeG, PCodeW
    REAL Width, VarBW, VarBS, VarBT, VarBG, EO, Tpot

    Common /bound_public/ NumBP, NSurf, NVarBW, NVarBS, NVarBT,        &
    NVarBG, NumSurfDat, NSeep, NSP(NSeepD), NP(NSeepD, NumSPD),        &
    NDrain, NDR(NDrainD), ND(NDrainD, NumDR),                          &
    KXB(NumBPD),                                                       &
    CodeW(NumNPD), CodeS(NumNPD), CodeT(NumNPD), CodeG(NumNPD),        &
    PCodeW(NumNPD), Width(NumBPD),                                     &
    VarBW(NumBPD, 3),                                                  &
    VarBS(NumBPD, NumSD), VarBT(NumBPD, 4),                            &
	VarBG(NumBPD, NumGD, 3), EO, Tpot,                                 &
    pondingByHead, pondingByFlux

    Double precision Time, Step, tNext,                                &
        dtOpt, dtMin, dMul1, dMul2, dtMx, tTDB, tFin, tatm
    Double precision Starter, TimeStep
    REAL Tinit
    Integer lInput, Iter
    Integer DailyOutput, HourlyOutput, RunFlag, DailyWeather,          &
    HourlyWeather
    Integer beginDay, sowingDay, emergeDay, endDay, OutputSoilNo,      &
    OutPutSoilYes, Year
    Integer iTime, iDawn, iDusk
    Integer tmpday, DayOfYear

    Common /time_public/ tNext(NumModD), dtMx(4), Time, Step,        &
     dtOpt,                                                            &
     dtMin, dMul1, dMul2, tTDB(4), Tfin, tAtm, Tinit,                  &
     lInput, Iter, DailyOutput, HourlyOutput, RunFlag,                 &
     DailyWeather, HourlyWeather,                                      &
     beginDay, sowingDay, emergeDay, endDay,                           &
     OutputSoilNo, OutPutSoilYes, Year,                                &
     iTime, iDawn, iDusk, TimeStep, tmpday, DayOfYear

     Integer NumMod, Movers, NShoot
     Common /module_public/ NumMod, Movers(4), NShoot

         
     character WeatherFile*256, TimeFile*256, BiologyFile*256,        &
        ClimateFile*256, NitrogenFile*256, SoluteFile*256,            &
        ParamGasFile*256, SoilFile*256,                               &
        ManagementFile*256, IrrigationFile*256, DripFile*256,         &
        WaterFile*256, WaterBoundaryFile*256,                         &
        PlantGraphics*256, InitialsFile*256, VarietyFile*256,         &
        NodeGraphics*256, ElemGraphics*256,                           &
        NodeGeomFile*256,                                             &
        GeometryFile*256, SurfaceGraphics*256,                        &
        FluxGraphics*256, MassBalanceFile*256,                        &
        MassBalanceFileOut*256, LeafGraphics*256,                     &
        OrganicMatterGraphics*256,                                    &                               
        RunFile*256, MassBalanceRunoffFileOut*256,                    &
        MulchFile * 256, MassBalanceMulchFileOut * 256

     Common /DataFilenames/ Starter, WeatherFile, TimeFile,         &
        BiologyFile, ClimateFile, NitrogenFile, SoluteFile,           &
        ParamGasFile, SoilFile,                                       &
        ManagementFile, IrrigationFile, DripFile,                     &
        WaterFile, WaterBoundaryFile,                                 &
        PlantGraphics, InitialsFile, VarietyFile,                     &
        NodeGraphics, ElemGraphics, NodeGeomFile,                     &
        GeometryFile, SurfaceGraphics,                                &
        FluxGraphics, MassBalanceFile,                                &
        MassBalanceFileOut, LeafGraphics,                             &
        OrganicMatterGraphics,                                        &
        RunFile, MassBalanceRunoffFileOut,                            &
        MulchFile, MassBalanceMulchFileOut

      Common /Materials/ BlkDn(NMatD), FracSind(NMatD),             &
        FracClay(NMatD),                                              &
        FracOM(NMatD), TUpperLimit(NMatD), TLowerLimit(NMatD),        &
        soilair(NumNPD), thSat(NMatD), ThAMin(NMatD),                 &
        ThATr(NMatD)

     end module common_block
