	!***********************************************************************
! This subroutine outputs to a graphics file.
! dt added leafwt and Stem weight 12/6/2006
!-----------------------------------------------------------------------
 SUBROUTINE Cropoutput
	 use common_block
     INCLUDE 'common.h'
	 
	 Character InString1*120
	 character * 10 date1
	 integer iday, iday1, hour
     REAL CWAD, GWAD, PWAD, LWAD, SWAD, AWAD, HIAD, SLAD, GWGD, GAD, difference
     REAL CumulativeNUptake, CumulativeNDemand,Pgday,Pnday

     common /output_G/ CumulativeNUptake, CumulativeNDemand,Pgday,    &
            Pnday, date1
  
	    If (lInput .eq. 1) then
          CumulativeNUptake=0.0
          CumulativeNDemand=0.0
        endif
        
	   PW=PodWt+SeedWt*1.13
	   TOPWT = STEMWT + PETWT + (SCPOOL+RSTPCS)*2.5*(STEMWT+PETWT)/SHTWT
	   LFWT = LEAFWT + (SCPOOL + RSTPCS)*2.5*LEAFWT/SHTWT
	   ALLWT  = LFWT + TOPWT + PW + ROOTWT + NODWT + ABSDW
	  ! ADD TO COMPARE AGMIP PROJECT
	   AWAD = (LFWT + TOPWT + PodWT + SeedWT*1.13 + ROOTWT) * PopArea * 10/1000  !Total ground biomass, dry (with grain)(kg/ha)
	   CWAD = (LFWT + TOPWT + PodWT + SeedWT*1.13) * PopArea * 10/1000  !Above ground biomass, dry (with grain)(kg/ha)
	   GWAD = SEEDWT * PopArea * 10/1000								!Grain dry weight (Mg/ha)
	   PWAD = (PodWT + SeedWT) * PopArea * 10/1000						!Pod biomass (seed + pod shell) (Mg/ha)
	   LWAD = LEAFWT * PopArea * 10									    !Leaf dry weight (kg/ha)
	   SWAD = (PETWT + STEMWT) * PopArea * 10                           !Stem dry weight (+ leaf petioles)(kg/ha)
	   IF (NSEEDS .GT. 0) GWGD = SEEDWT/NSEEDS*1000					    !Unit grain weight (mg)
	   HIAD = SEEDWT/(LEAFWT + PETWT + STEMWT + PodWT + SeedWT)			!Harvest index at R8
	   SLAD = LAREAT/LEAFWT											    !Specific leaf area(cm2/g)
	   GAD  = NSEEDS*PopArea											!Grain number(number/m2)
	   
	   
	   PGROSS = PGR*0.9818/POPAREA*(30/12)						  !conversion mg CO2 m-2 ground s-1	to g CH2O plant h-
	   PSNET = PN*0.9818/POPAREA*(30/12)						  !conversion mg CO2 m-2 ground s-1	to g CH2O plant h-
	   PFD = PAR(ITIME) * 4.55                    !conversion from PAR in W m-2 to umol s-1 m-2 
       POTET = transpiration*0.018*3600/PopArea   !for hourly output, Potential transpiration rate, in g h2o plant-1 h-1
       ACTET = WaterUptake/POPSLB                 !for hourly output, Actural transpiration rate, in g h2o plant-1 h-1
	  ! Calculate the plant nitrogen content
	   LeafNitrogen=LFWT*0.05*1000			  !leaf N content in plant, mg N plant-1
	   StemNitrogen=TOPWT*0.02*1000			  !Stem N content in plant, mg N plant-1
	   PodNitrogen=PodWt*0.04*1000			  !Pod N content in plant,  mg N plant-1
	   SeedNitrogen=SEEDWT*1.13*0.065*1000	  !Seed N content in plant, mg N plant-1
	   RootNitrogen=ROOTWT*0.025*1000		  !root N content in plant, mg N plant-1
	   DeadNitrogen=ABSDW*0.05*1000			  !dead N content in plant, mg N plant-1 
	   PlantNitrogen=LeafNitrogen+StemNitrogen+PodNitrogen  &
	   +SeedNitrogen+RootNitrogen+DeadNitrogen !total N content in plant, mg N plant-1
	   
	   NCRatio= PlantNitrogen/1000/USEDCS	  !N:C ratio of leaves, g N g-1
	   CumulativeNUptake=CumulativeNUptake+NitrogenUptake/POPSLB*1000	!total N uptake, cumulative, mg N plant-1
	   CumulativeNDemand=CumulativeNDemand+HourlyNitrogenDemand*1000	!total N demand, cumulative, mg N plant-1
	   
	   If (ACTET.eq.0) then
		   WSTRESS = 1.0
	   else
		   WSTRESS = ACTET/POTET		!Hourly water stress
	   end if 
	   If (WSTRESS.ge.1.0) WSTRESS = 1.0
	   If (Cstress.ge.1.0) Cstress = 1.0
	  !Nitrogen Stress 
	   If (NitrogenUptake.eq.0) then
		   NSTRESS = 1.0
	   else
		   NSTRESS = HourlyNitrogenDemand/(NitrogenUptake/POPSLB)	!Hourly Nitrogen stress
	   end if 
	   If (NSTRESS.ge.1.0) NSTRESS = 1.0
	  !Calucaltion of Daily Output
	   difference=abs(time-int(time))
	   If (difference.ge.0.999.or.difference.le.0.001) then         !first hour of day
		   Parave = 0.0    !mol m-2 timeinc-1 Daily average PFD in 24h
		   Sradave = 0.0   !MJ srad m-2 timeinc-1 Daily average Solrad in 24h
	       Tdayave = 0.0   !Daily average AIR temperature oC
		   Tcanave = 0.0   !Daily average canopy temperature oC
		   Pgday = 0.0     !daily PGROSS g CH2O plant d-1
		   Pnday = 0.0     !daily Photosynthesis g CH2O plant d-1
		   gsmax = 0.0     !maximum gs mmol H2O m-2 s-1
		   psilave = 0.0   !average psil leaf water potential bars
		   daypotet = 0.0  !potet g plant-1 h-1 daypotet g plant-1 day-1
		   dayactet = 0.0  !actet g plant-1 h-1 dayactet g plant-1 day-1
		   wstressave = 0.0!average waterstress from 0 to 1
		   nstressave = 0.0!average nitrogenstress from 0 to 1
		   Tleafmax = 0.0   !calculate the maximum leaf temperature
	   end if 
      
	   If (difference.ge.0.0) then         !then increment
		   Parave = Parave + PFD*timeStep*60/1000000    !mol m-2 timeinc-1
		   Sradave = Sradave + PFD*timeStep*60/4.57*2/1000000     !MJ srad m-2 timeinc-1
	       Tdayave = Tdayave + TAIR(ITIME)
		   Tcanave = Tcanave + temperature !Daily average canopy temperature oC
		   Pgday = Pgday + PGROSS          !daily PGROSS g CH2O plant d-1
		   Pnday = Pnday + PSNET           !daily Photosynthesis g CH2O plant d-1
		   if (Ags > gsmax) gsmax = Ags
		   if (TAIR(ITIME) > Tleafmax) Tleafmax = TAIR(ITIME)
		   psilave = psilave + psil_
		   daypotet = daypotet + POTET
		   dayactet = dayactet + ACTET
		   wstressave = wstressave + WSTRESS
		   nstressave = nstressave + NSTRESS
	   end if 
	   
	   If((difference.le.0.96).and.(difference.ge.0.95)) then  !last timestep prior to 24:00 or 0:00 next day
	       Tdayave = Tdayave/(24*60/timeStep)
		   Tcanave = Tcanave/(24*60/timeStep)!Daily average canopy temperature oC  
		   psilave = psilave/(24*60/timeStep)
		   daypotet = daypotet + POTET
		   dayactet = dayactet + ACTET
		   wstressave = wstressave/(24*60/timeStep)
		   nstressave = nstressave/(24*60/timeStep)
		end if    
		   
		   
	  !g01 Daily OUTPUT 
	   
       If((DailyOutput.eq.1).and.((difference.le.0.96).and.(difference.ge.0.95))) then
		iday=int(time)
		iday1=int(time)
        call caldat(iday,mm,id,iyyy) 
        write (date1,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
		Write(85,7) date1, iday1, 23, RSTAGE, VSTAGE, Parave, Sradave,    &
		              Tdayave, Tcanave, Pgday, Pnday, gsmax, psilave, LAI,		&
					  LAREAT, ALLWT, ROOTWT, STEMWT, LFWT, SeedWT, PodWt, ABSDW,&
			          daypotet, dayactet, wstressave, nstressave,               &
                        CumulativeNUptake/1000.0,                      &
                        CumulativeNDemand/1000.0,                      &
                        LIMITF(12) 
					  
	   end if 
	  
	  ! g01 Hourly output
       If(HourlyOutput.eq.1) then 
          If (difference.ge.0.999.or.difference.le.0.001) then	  
		   iday=int(time)
             call caldat(iday,mm,id,iyyy) 
             write (date1,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
           end if 
		  hour=int(difference*24 + 0.1)
		  iday1=int(time)
		  Write(85,7) date1, iday1, hour, RSTAGE, VSTAGE, PFD, WATTSM(ITIME), &
		              TAIR(ITIME), temperature, PGROSS, PSNET, Ags, psil_, LAI,	  &
					  LAREAT, ALLWT, ROOTWT, STEMWT, LFWT, SeedWT, PodWt, ABSDW,  &
					  POTET, ACTET, WSTRESS, NSTRESS,                &
                        CumulativeNUptake/1000.0,                      &
                        CumulativeNDemand/1000.0,                      &
                        LIMITF(ITIME)
	   end if 
	
	  ! daily g02 file
	  
      ! calculated average canopy temperature  
	  If (difference.ge.0.999.or.difference.le.0.001) THEN 
		  Avtemperature = 0.0 
		  howtime = 0.0
	  END IF  
	  If ((PFD.GE.50).AND.(LAI.GE.3)) THEN
		  Avtemperature = Avtemperature + temperature
		  howtime = howtime + 1.0
	  end if 
	 
	  If ((difference.le.0.96).and.(difference.ge.0.95)) THEN 
		  Avtemperature = Avtemperature/howtime
	  end if 
	  
		 
	  ! net sunlit photosynthesis   photosynthesis_netsunlitleaf
	  If (difference.ge.0.999.or.difference.le.0.001) THEN 
		  netsunlitleaf = 0.0 
	  END IF 
	  
	  If (difference.ge.0.0) then  
		  if (photosynthesis_netsunlitleaf > netsunlitleaf) netsunlitleaf = photosynthesis_netsunlitleaf
	  END IF
	  
	  If ((difference.le.0.39).and.(difference.ge.0.36)) THEN 
		  AAgs = Ags
	  END IF
	  
	  !Calculate the EDD(LEAF TEMPERATURE) from June 1 to August 30
	  If ((DayOfYear.GE.152).AND.(DayOfYear.LE.243)) Then
		  if (temperature.GE.30) then
			  DDT = (temperature-30)/24
			  EDD = EDD + DDT
		  end if
	  end if 
	  !Calculate the EDD(AIR TEMPERATURE) from June 1 to August 30
	  If ((DayOfYear.GE.152).AND.(DayOfYear.LE.243)) Then
		  if (TAIR(ITIME).GE.30) then
			  ADDT = (TAIR(ITIME)-30)/24
			  AEDD = AEDD + ADDT
		  end if
	  end if 
	  !Calculate the EDD from GROWING SEASON
		  if (TAIR(ITIME).GE.30) then
			  TDDT = (TAIR(ITIME)-30)/24
			  TEDD = AEDD + ADDT
		  end if
	  
	  !Calculated the total rain  from June 1 to August 30
	  If ((DayOfYear.GE.152).AND.(DayOfYear.LE.243)) Then
		  If (difference.ge.0.999.or.difference.le.0.001) THEN
			  JJARAIN =  JJARAIN + RAIN*24
		  end if
	  end if 
	  !Calculated the total rain  growing season
	  
	       If (difference.ge.0.999.or.difference.le.0.001) THEN
			  TRAIN =  TRAIN + RAIN*24
	       end if
	   
	   
	  !Calculate the TOTAL PET AND ET growing season
	  If ((difference.le.0.96).and.(difference.ge.0.95)) then
	              TPET = TPET + daypotet*poparea/1000 !mm TOATL potential ET
				  TAET = TAET + dayactet*poparea/1000 !mm TOATL actual  ET
	  End if
	  
	  !Calculate the JJAPET AND ET from June 1 to August 30
	  If ((DayOfYear.GE.152).AND.(DayOfYear.LE.243)) Then
	      If ((difference.le.0.96).and.(difference.ge.0.95)) then
	              JJAPET = JJAPET + daypotet*poparea/1000 !mm TOATL potential ET
				  JJAAET = JJAAET + dayactet*poparea/1000 !mm TOATL actual  ET
		  End if
	  End if
	  
	  !Calculate the average P-PET from June 1 to August 30
	  If ((DayOfYear.GE.152).AND.(DayOfYear.LE.243)) Then
		If ((difference.le.0.96).and.(difference.ge.0.95)) then
	              PPET = PPET + RAIN - daypotet*poparea/1000 !mm TOTAL RAIN minus potential ET
		End if
	  End if 
	  !Calculate the JJAVPD from June 1 to August 30
	  If ((DayOfYear.GE.152).AND.(DayOfYear.LE.243)) Then	
		    JJAVPD = JJAVPD + VPD(ITIME)/1000                     !JJA VPD(pa)	
	  End if 
	  !Calculate the TVPD growing season
			TVPD = TVPD + VPD(ITIME)/1000                        !TOTAL VPD(pa)	
	  
	  !Calculate the daily biomass (g/plant/day)
	   If (difference.ge.0.999.or.difference.le.0.001) then  !first hour of day
		   Biomass1 = ALLWT
	   END IF 
	   If((difference.le.0.96).and.(difference.ge.0.95)) then  
		   DailyBio	= ALLWT - Biomass1		!daily biomass grwoth rate.
		   if (wstressave.le.0.7) then 
		   wstressday = wstressday + 1
		   end if 
	   end if 
	   
	   
	   
	   !If (((DayOfYear.le.243).and.(DayOfYear.ge.182))) then  
	     !If((DailyOutput.eq.1).and.((difference.le.0.96).and.(difference.ge.0.95))) then
		 !iday=int(time)
		 !iday1=int(time)
           !call caldat(iday,mm,id,iyyy) 
           !write (date1,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy   
	       !Write(86,8) DayOfYear, Tleafmax, wstressday, daypotet, DailyBio, wstressave, dayactet
		!end if   
		 !date1,CWAD, GWAD, netsunlitleaf,    &
		 !            AAgs, Avtemperature, EDD, AEDD, PPET
	   !end if
	   
	   
	   
	  !For output of regional simulation
	  IF (MATURE .EQ. 1) WRITE (86,8) R1DAY, R4DAY, R5DAY, R8DAY,                                 &
									    GAD,GWAD,CWAD,PWAD,GWGD,HIAD,EDD,AEDD,TEDD,               &
	                                    JJARAIN,TRAIN,TPET,TAET,JJAPET,JJAAET,PPET,JJAVPD,TVPD,   &
										wstressday	
	  
	   

	 
      
	  !Nitrogen.crp file
	  !Nitrogen Daily OUTPUT 
      !changed to output N uptake and demand to be g /plant rather than mg/plant                                   
	   
       If((DailyOutput.eq.1).and.((difference.le.0.96).and.(difference.ge.0.95))) then
		iday=int(time)
		iday1=int(time)
        call caldat(iday,mm,id,iyyy) 
        write (date1,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
		Write(87,9) date1, iday1, 23, PlantNitrogen, LeafNitrogen, StemNitrogen, &
		      PodNitrogen, SeedNitrogen, RootNitrogen, DeadNitrogen, NCRatio,	&
              CumulativeNUptake/1000.0, CumulativeNDemand/1000.0, nstressave
					  
	   end if 
	  
	  ! nitrogen Hourly output
       If(HourlyOutput.eq.1) then 
	      If (difference.ge.0.999.or.difference.le.0.001) then
		  iday=int(time)
          call caldat(iday,mm,id,iyyy) 
          write (date1,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
		  end if 
		  hour=int(difference*24 + 0.1)
		  iday1=int(time)
		  Write(87,9) date1, iday1, hour, PlantNitrogen, LeafNitrogen, StemNitrogen, &
		      PodNitrogen, SeedNitrogen, RootNitrogen, DeadNitrogen, NCRatio,		&
              CumulativeNUptake/1000.0, CumulativeNDemand/1000.0, NSTRESS
		  
	   end if 
	   
	  !plantstress.crp file
	  !plantstress Daily OUTPUT 
	   
       If((DailyOutput.eq.1).and.((difference.le.0.96).and.(difference.ge.0.95))) then
		iday=int(time)
		iday1=int(time)
        call caldat(iday,mm,id,iyyy) 
        write (date1,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
		Write(88,10) date1, iday1, 23, wstressave, nstressave, Cstress, &
					NRATIO, SGT  
	   end if 
	 
	  ! plantstress Hourly output
       If(HourlyOutput.eq.1) then 
	      If (difference.ge.0.999.or.difference.le.0.001) then
		  iday=int(time)
          call caldat(iday,mm,id,iyyy) 
          write (date1,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
		  end if 
		  hour=int(difference*24 + 0.1)
		  iday1=int(time)
		  Write(88,10) date1, iday1, hour, WSTRESS, NSTRESS, Cstress,  &
			  NRATIO, SGT
		  
	   end if 
	  
7	   FORMAT (A12, (",",I8), (",",I4), 25(",",F15.3), (",",A5))
8	   FORMAT (I4, 3(",",I4),(",",F12.1), 18(",",F6.1))
9	   FORMAT (A12, (",",I8), (",",I4), 8(",",F15.3),2(",",F10.5),",",F15.3)	  
10     FORMAT (A12, (",",I8), (",",I4), 5(",",F15.3))
	   RETURN
 End   