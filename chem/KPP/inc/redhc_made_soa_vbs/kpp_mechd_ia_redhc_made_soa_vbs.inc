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
     roleho  = RCONST(56)
     roleo3  = RCONST(57)
     roleno3 = RCONST(58)

     ! TOL
     rtolho =   RCONST(75)

     ! NO XYL

     ! NO CSL

     ! PAR, RACM_SOA is HC5 and HC8
     rsoaalkho  =  RCONST(96)

     ! ISO
     risoho  =  RCONST(62)
     risoo3  =  RCONST(63)
     risono3 =  RCONST(64)

     ! C10H16, RACM is seperate for API and LIM
     rc10h16ho  = RCONST(77)
     rc10h16o3  = RCONST(76)
     rc10h16no3 = RCONST(78)

     ! SESQ
     rsesqho  =  RCONST(80)
     rsesqo3  =  RCONST(79)
     rsesqno3 =  RCONST(81)

     ! IVOC
     rivocho = RCONST(95)

     ! SESQ
     rsesqho  =  RCONST(10)
     rsesqo3  =  RCONST(10)
     rsesqno3 =  RCONST(10)

     ! IVOC 
     rivocho = RCONST(10)
     
     ! NO Isoprene radical, RACM uses ISOP + HO, NO, MO2, ACO3 to deterimine branching ratio, what whould we use?
! ---------------------    
  
     ! production from anthropogenic VOCs
     PRDROG(PSOAALK)= rsoaalkho*var(ind_soaalk)*var(ind_oh)

     PRDROG(POLE)= roleho*var(ind_ole)*var(ind_oh) + roleo3*var(ind_ole)*var(ind_o3) + roleno3*var(ind_ole)*var(ind_no3)

     PRDROG(PARO)= rtolho*var(ind_toluene)*var(ind_oh)

     PRDROG(PIVOC)= rivocho*var(ind_ivoc)*var(ind_oh)


     ! Biogenic
     PRDROG(PISOP)= risoho*var(ind_isop)*var(ind_oh) + risoo3*var(ind_isop)*var(ind_o3) + risono3*var(ind_isop)*var(ind_no3)

     PRDROG(PTERP)= rc10h16ho*var(ind_c10h16)*var(ind_oh) + rc10h16o3*var(ind_c10h16)*var(ind_o3) + rc10h16no3*var(ind_c10h16)*var(ind_no3)

     PRDROG(PSESQ)= rsesqho*var(ind_sesq)*var(ind_oh) + rsesqo3*var(ind_sesq)*var(ind_o3) + rsesqno3*var(ind_sesq)*var(ind_no3)

! RAR: to calculate the branching ratios to determine high NOx versus low NOx

!!!!!!     PRDROG(PBRCH)= risopno*var(ind_no)

     ! VDROG carrying the branching ratios
!     if (PRDROG(PBRCH)>1.E-12) then
!        ro2loss= PRDROG(PBRCH) + risopho2*var(ind_ho2) + risopmo2*var(ind_mo2) + risopaco3*var(ind_aco3) +  &
!                               risopisop*var(ind_isopp)
!        VDROG3_SOA2( i,k,j,LDROG_VBS )= MIN( 1.D0,(PRDROG(PBRCH)/ro2loss) )
!     else
!         VDROG3_SOA2( i,k,j,LDROG_VBS )= 0.
!     end if
!
     DO n = 1, ldrog_soa2-1
        VDROG3_SOA2( i,k,j,n ) =  oconv* PRDROG( n ) * DTSTEPC
        VDROG3_SOA2( i,k,j,n ) =  MAX( 0., VDROG3_SOA2( i,k,j,n ) )
     ENDDO
! All the produced OCV at bin-0 (CVASOA0) goes to the aerosol phase
! Below is the VDROG3_SOA2 component carrying CVASOA0 production
    VDROG3_SOA2(i,k,j,POCV0)= MAX( 0.0_dp, oconv* var(ind_cvasoa0) )
    var(ind_cvasoa0)= 0.0_dp

endif

! RAR: debugging
!if (i==8 .AND. j==18) then
!   if (k==1) then
!        write(*,*)'rhch5ho',rhc5ho,'rhc8ho',rhc8ho,'rhc8ho',roltho,'roliho',roliho, &
!                  'rtolho',rtolho,'rxylho',rxylho,'rsesqno3',rsesqno3
!        write(*,*)'ind_toluene',ind_toluene,'var(ind_toluene)',var(ind_toluene)
!        write(*,*)'ind_oh',ind_oh,'var(ind_oh)',var(ind_oh)
!        write(*,*)'ind_isop',ind_isop,'risoho',risoho
!        write(*,*)'PRDROG(PBRCH)', PRDROG(PBRCH),'ro2loss=',ro2loss        
!        write(*,*)'VDROG3(8,1,18,:)', VDROG3(i,k,j,:)
!   end if
!end if
!
!if (j==18 .AND. k==1) then
!        write(*,*)'VDROG3(:,18,1,:)', VDROG3(i,k,j,:)
!end if
