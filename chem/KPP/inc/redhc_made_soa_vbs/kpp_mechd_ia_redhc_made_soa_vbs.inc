! RAR: modified to handle the new SOA mechanism based on the VBS approach and multi-generational
! VOC oxidation mechanism
! Correspondence between RACM and SAPRC-99 SOA precursors
! 1) OLT -> OLE1
! 2) OLI -> OLE2
! 3) TOL -> ARO1
! 4) XYL -> ARO2
! 5) CSL -> ARO2
! 6) HC5 -> ALK4
! 7) HC8 -> ALK5
! 8) ISO -> ISO
! 9) API -> TERP
!10) LIM -> TERP
!11) SESQ-> SESQ
!
if(p_nu0.gt.1)then

     ! OLE only, RACM_SOA is OLI and OLT
     roleho  = RCONST(00)
     roleo3  = RCONST(00)
     roleno3 = RCONST(00)

     ! TOL
     rtolho =   RCONST(00)

     ! NO XYL

     ! NO CSL

     ! PAR, RACM_SOA is HC5 and HC8 
     rparho  =  RCONST(10)

     ! ISO
     risoho  =  RCONST(10)
     risoo3  =  RCONST(10)
     risono3 =  RCONST(10)

     ! C10H16, RACM is seperate for API and LIM
     rc10h16ho  = RCONST(10)
     rc10h16o3  = RCONST(10)
     rc10h16no3 = RCONST(10)    
 
     ! SESQ
     rsesqho  =  RCONST(10)
     rsesqo3  =  RCONST(10)
     rsesqno3 =  RCONST(10)

     ! NO Isoprene radical, RACM uses ISOP + HO, NO, MO2, ACO3 to deterimine branching ratio, what whould we use?
! ---------------------    
  
     ! production from anthropogenic VOCs
     PRDROG(PPAR)= rparho*var(ind_par)*var(ind_ho)

     PRDROG(POLE)= roleho*var(ind_ole)*var(ind_ho) + roleo3*var(ind_ole)*var(ind_o3) + roleno3*var(ind_ole)*var(ind_no3)

     PRDROG(PARO)= rtolho*var(ind_tol)*var(ind_ho)

     ! Biogenic
     PRDROG(PISOP)= risoho*var(ind_iso)*var(ind_ho) + risoo3*var(ind_iso)*var(ind_o3) + risono3*var(ind_iso)*var(ind_no3)

     PRDROG(PTERP)= rc10h16ho*var(ind_c10h16)*var(ind_ho) + rc10h16o3*var(ind_c10h16)*var(ind_o3) + rc10h16no3*var(ind_c10h16)*var(ind_no3)

     PRDROG(PSESQ)= rsesqho*var(ind_sesq)*var(ind_ho) + rsesqo3*var(ind_sesq)*var(ind_o3) + rsesqno3*var(ind_sesq)*var(ind_no3)

! RAR: to calculate the branching ratios to determine high NOx versus low NOx

     PRDROG(PBRCH)= risopno*var(ind_no)

     ! VDROG carrying the branching ratios
!     if (PRDROG(PBRCH)>1.E-12) then
!        ro2loss= PRDROG(PBRCH) + risopho2*var(ind_ho2) + risopmo2*var(ind_mo2) + risopaco3*var(ind_aco3) +  &
!                               risopisop*var(ind_isop)
!        VDROG3_VBS( i,k,j,LDROG_VBS )= MIN( 1.D0,(PRDROG(PBRCH)/ro2loss) )
!     else
!         VDROG3_VBS( i,k,j,LDROG_VBS )= 0.
!     end if
!
     DO n = 1, LDROG_VBS-1
        VDROG3_VBS( i,k,j,n ) =  oconv* PRDROG( n ) * DTSTEPC
        VDROG3_VBS( i,k,j,n ) =  MAX( 0., VDROG3_VBS( i,k,j,n ) )
     ENDDO
endif

! RAR: debugging
!if (i==8 .AND. j==18) then
!   if (k==1) then
!        write(*,*)'rhch5ho',rhc5ho,'rhc8ho',rhc8ho,'rhc8ho',roltho,'roliho',roliho, &
!                  'rtolho',rtolho,'rxylho',rxylho,'rsesqno3',rsesqno3
!        write(*,*)'ind_tol',ind_tol,'var(ind_tol)',var(ind_tol)
!        write(*,*)'ind_ho',ind_ho,'var(ind_ho)',var(ind_ho)
!        write(*,*)'ind_iso',ind_iso,'risoho',risoho
!        write(*,*)'PRDROG(PBRCH)', PRDROG(PBRCH),'ro2loss=',ro2loss        
!        write(*,*)'VDROG3(8,1,18,:)', VDROG3(i,k,j,:)
!   end if
!end if
!
!if (j==18 .AND. k==1) then
!        write(*,*)'VDROG3(:,18,1,:)', VDROG3(i,k,j,:)
!end if
