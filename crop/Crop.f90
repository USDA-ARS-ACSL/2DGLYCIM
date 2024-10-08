! we will call this from the main program
! DEC$ statement tells compiler to expose the subroutine name and 
! common  block
!  Crop - subroutine 

    subroutine crop()
	use common_block
    Include 'common.h'
	
    !DEC$ATTRIBUTES DLLEXPORT :: crop
    Real CurrentNUptakeError, CumulativeNuptakeError

    common /crop_main/CumulativeNUptakeError
   
	!Write(*,*) Time 
    !***initializing crop module----------	
    If (lInput .eq. 1) then
		Write(*,*) 'Initializing crop module'
	    Period = (timeStep/60) !This is 1h period which only used in crop. 
		IPERD = 24
		INITR0 = 0
        EQSC   = 0.0
		POPSLB  = PopRow/100.0*EOMult!note: this is really-> poprow * (100/row-spacing) * (row-spacing / 10000) * EOMult
		!Convr = 1.0  !wsun all carbon currency should be on C basis both 2dsoil and crop
		AWUPS = 0.0   !initialize AWUPS, AWUPS_old and psil_ in 2DSOIL 
		psil_ = -0.5  !note, this is in bars
		ET_demand = 0.0 
		NitrogenUptake = 0 !should be zero at emergence since all plant N has come from seedpiece
		NDemandError = 0.0
		CumulativeNDemandError = 0.0
        CumulativeNUptakeError=0.0
	    RunFlag = 1
		Tgrowth = 0
        NumMod=NumMod+1
        ModNum=NumMod
        Mature=0
		tNext(ModNum) = emergeDay
        !INIT = .TRUE.
		
    !***A new plant model object is created and initialized (calls initialize function) here
        CALL PRI_INIT	! Iintial the common crop block
		CALL Initialize ! Read variety file and generate OUTPUT with plant related file (PlantGraphics,LeafGraphics) 
		MSTEMW = STEMWT
        SHTWT = LEAFWT + PETWT + STEMWT 
        LAREAM = LAREAT
		SDWMAX = 544.32/SEEDLB           !   544.32 = 453.6*1.2
    End if

    !***End Initialization
     If (Mature .eq. 1) return
    if (Abs(Time-sowingDay+1).lt.0.001) isGerminated = 1  
    !Simulate nitrogen and water uptake for each interval prior to 60min
   
    If (NShoot .GT. 0) then
		WaterUptake = WaterUptake + AWUPS*Step			   !g water per slab taken up in an hour
		NitrogenUptake = NitrogenUptake + SIncrSink/1.0e6  !Cumulative N (mass, g (plant slab)-1) in this time step
		HourlyCarboUsed = HourlyCarboUsed + PCRS*Step 
		
	End if 	 
   

 !If the emergence date has come then start the crop model
    ! code inside this if statement is executed every hour once the crop has emerged
  if (abs(time-tNext(modnum)).lt. 0.001*step) then
       if (NShoot .eq. 0) then 
	    NShoot=1
	    INIT = .TRUE.
        isEmerged = 1
      endif
  	   
 
    !check to make sure when crop is matured, dll is exited (avoids error message at end of simulation
    !If (Simulationdone .eq. 1) GOTO 1111 

        !Below is the values for gas exchanger.dll
		  CDayOfYear = DayOfYear
		  CITIME = ITIME
		  CIPERD = IPERD
		  CWATTSM(CITIME) = WATTSM(ITIME)
          CPAR(CITIME) = PAR(ITIME)
	      CTAIR(CITIME) = TAIR(ITIME)
	      CCO2 = CO2
	      CVPD(CITIME) =  VPD(ITIME)
	      CWIND = WIND
	      CPSIL_ = PSIL_
          CLATUDE = LATUDE
	      CLAREAT = LAREAT
	      CLAI = LAI
		  
!this should be comparison just for current hourly time-step - note NitrogenUptake reset to zero each hour, so above line was incorrect comparison          
		CurrentNUptakeError = NitrogenUptake/POPSLB  - HourlyNitrogenDemand
        CumulativeNUptakeError = CumulativeNUptakeError + CurrentNUptakeError

!Start here will call the main PLANT 
	        CALL PNET
            CALL PHEN
            CALL SHTGRO     ! WasPOTGRO
            CALL PARTIT  
			CALL ROOT
			CALL NITRO      ! WAS NUTRTS
			CALL NODULE
            CALL ACTGRO
            CALL TISLOS		
			INIT = .FALSE.
			
			!LCAI = LAREAT*PopArea/(100*100)
            !hieght is minimum of 2 cm when germinated
            Shade = Cover*RowSp
			Height=max(min(Shade,RowSp),2.0)
			LCAI = LAREAT*POPROW/100.0/MIN(HEIGHT,ROWSP)
			Cover = 1.0 - exp (-0.79*LCAI)
			Shade = Cover*RowSp
			Height=max(min(Shade,RowSp),2.0)
			ET_demand = transpiration*0.018*3600*24/PopArea !pass ET demand from shoot to root, g plant-1 d-1 //instanteous transpiration, mmol H2O m-2 ground s-1
			LAI = LAREAT*PopArea/(100*100)
			nitroDemand = HourlyNitrogenDemand*POPSLB*1e6*24 !pass the N demand into 2dsoil, units are ug N slab-1 d-1
			NDemandError = CurrentNUptakeError !but these units are g N plant-1 ?
			CumulativeNDemandError = CumulativeNUptakeError !these units are g N plant-1 - why is nitroDemand in ug N slab-1?
		
		    CALL Cropoutput
            
    			
	    if(Mature .eq. 1 ) then
		    !Simulationdone = 1
		    NShoot = 0             !tell 2dsoil that crops harvested
            tnext(ModNum)=1.e+32
        else
            tnext(Modnum) = Time + Period/24
			WaterUptake = 0
			NitrogenUptake = 0
			TotalPotentialRootWaterUptake = 0

			
	    end if
    end if  ! end of if statement for hourly time step calling crop model
    
!***End of Daily crop routine***
!1111 return
	end

