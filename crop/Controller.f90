
!***********************************************************************
! This subroutine outputs to a graphics file.
! dt added leafwt and Stem weight 12/6/2006
!-----------------------------------------------------------------------
 SUBROUTINE Initialize
	 
	 use common_block
     INCLUDE 'common.h'
	 
     Character InString1*120,                                   &
	 nitrogencrp*256, plantstresscrp*256
     Character*256 extract_path, path
	 
  
   Write(*,*)'************************ GLYCIM ***********************'  
   Write(*,*)'*                    Version 9.1.1                    *'  
   Write(*,*)'*                                                     *'  
   Write(*,*)'*  A DYNAMIC SIMULATOR FOR SOYBEAN CROPS ORIGINALLY   *'  
   Write(*,*)'*  CREATED BY BASIL ACOCK AND SUBSEQUENTLY IMPROVED   *'  
   Write(*,*)'*  BY YAKOV PACHEPSKY, V.R. REDDY, DENNIS TIMLIN, AND *'  
   Write(*,*)'*  ANTHONY TRENT.                                     *'
   Write(*,*)'*  INCORPORATED INTO 2DSOIL AND LINKED WITH           *'
   Write(*,*)'*    FARQUHAR PHOTOSYNTHESIS MODEL BY                 *'
   write(*,*)'*    BY WENGUANG SUN                                  *'
   Write(*,*)'*  USDA-ARS,ADAPTIVE CROPPING SYSTEMS LABORATORY      *'  
   Write(*,*)'*  BELTSVILLE, MD  20705.   TEL:(301)504-5872         *'
   write(*,*)'*                                                     *'   
   Write(*,*)'*******************************************************'
   
   !Read variety file 

   Open(40,file = VarietyFile, status='old', ERR=13)

      im=220
      il=0
      Read(40,*, ERR=13)
      im=im+1
      il=il+1
      Read(40,*, ERR=13)
      im=im+1
      il=il+1
      Read(40,*, ERR=13)
      im=im+1
      il=il+1
      Read(40,*, ERR=13) MG, SEEDLB, FILL,						  &
	  PARM(2), PARM(3), PARM(4), PARM(5),                         &
	  PARM(6), PARM(7), PARM(8), PARM(9), PARM(10),               &
	  PARM(11), PARM(12), PARM(13), PARM(14), PARM(15),           &
	  PARM(16), PARM(17), PARM(18), PARM(19), PARM(20),           &
	  PARM(21), PARM(22), PARM(23), PARM(24), PARM(25)
      close(40)
	  
	  Open(85,file=PlantGraphics)  
	  Open(86,file=LeafGraphics)
	  
	  Path=extract_path(PlantGraphics)
      
	  nitrogencrp=trim(Path)//'nitrogen.crp'
	  plantstresscrp=trim(Path)//'plantstress.crp'
	  
      	   
	  Open(87,file=nitrogencrp)
	  Open(88,file=plantstresscrp)
	  
	  Write(85,5) "date", "jday", "time", "RSTAGE", "VSTAGE",         &
	  "PFD", "SolRad", "Tair", "Tcan", "Pgross", "Pnet", "gs", "PSIL",&
	  "LAI", "LAREAT", "totalDM", "rootDM", "stemDM", "leafDM",	      &
	  "seedDM", "podDM", "DeadDM", "Tr_pot", "Tr_act", "wstress", 	  &
	  "Nstress","N_Dem","NUPt", "Limit"
	 
	  !Write(86,6) "date", "DayOfYear", "Time", "CWAD", "GWAD",		  &
	  !"PN", "Ags", "Canopy T"
	  
	  Write(87,7) "date", "jday","time","total_N","leaf_N","stem_N",  &
		          "pod_N", "seed_N", "root_N", "dead_N", "plant_N_C", &
	              "N_uptake", "N_demand", "Nstress"
      Write(88,8) "date", "jday","time", "wstress", "Nstress",        &
	              "Cstress", "NEffect_veg", "wstress2"
	 
5	  FORMAT (A10, 28(",",A10))
6	  FORMAT (6(A10))
7	  FORMAT (A10, 13(",",A10))
8	  FORMAT (A10, 7(",",A10))
	  
	  RETURN
13	  Write(*,*) 'Error in VARIETY FILE'
	  GOTO 11
11	  CONTINUE	         

   END
   
   
       function extract_path(filename)
       character *256 filename, path, extract_path
       integer :: i, len

       len = len_trim(filename)
       path = ""
    ! Find the last occurrence of the directory separator '\'
       
       do i = len, 1, -1
          if ((filename(i:i) == '\').OR.(filename(i:i) == '/')) then
              path = filename(1:i)
              if (filename(i:i) == '/') then   ! if windows
                  path=path // '/'
               else 
                 path=path // '\'               ! if linux
              end if
             exit
          end if
       end do

       extract_path = path
       end function extract_path
       