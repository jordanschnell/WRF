!
!
! CHECK, LUMPING
!
! Things to update:
!
! 1. rhc8oh --> set to some fraction of PAR
! 2. rOLI/OLT --> set to some fraction of OLE
!    --currently only set for 100% OLI for testing purposes
! 3. Partition API/LIM products
!

        if(p_nu0.gt.1)then


            rxylho  =   0._dp;

            rtolho  =   5.5e-11_dp;

            rcslho  =   0._dp;

            rcslno3 =   0._dp;

            rhc8ho  =   0._dp;

            roliho  =   ARR2( 5.2e-12_dp, -504._dp, TEMP);
            rolino3 =   7.7e-15_dp;
            rolio3  =   ARR2( 4.33e-15_dp, 1800._dp, TEMP);

            roltho  =   0._dp;
            roltno3 =   0._dp;
            rolto3  =   0._dp

! For API, set to folllow C10H16 Rxns but sorgam driver will not currently use since p_lim ==1
            rapiho = 6.e-11_dp;
            rapino3 = 4.e-12_dp;
            rapio3 = 1.e-16_dp;
! Maybe do same for LIM and set some partitioning ratio 
            rlimho = 0._dp;
            rlimno3 =  0._dp;
            rlimo3 = 0._dp;
!!!
!!!!!!!!!!!
            PRDROG(PXYL)  = 0._dp

            PRDROG(PTOL)  = rtolho * var(ind_toluene)*var(ind_oh)

            PRDROG(PCSL1) = 0._dp
            PRDROG(PCSL2) = 0._dp

            PRDROG(PHC8)  = rhc8ho * var(ind_par)*var(ind_oh)

            PRDROG(POLI1) = roliho * var(ind_ole)*var(ind_oh)
            PRDROG(POLI2) = rolino3 * var(ind_ole)*var(ind_no3)
            PRDROG(POLI3) = rolio3 * var(ind_ole)*var(ind_o3)
           
            PRDROG(POLT1) = roltho * var(ind_ole)*var(ind_oh)
            PRDROG(POLT2) = roltno3 * var(ind_ole)*var(ind_no3)
            PRDROG(POLT3) = rolto3 * var(ind_ole)*var(ind_o3)

            PRDROG(PAPI1) = 0._dp
            PRDROG(PAPI2) = 0._dp
            PRDROG(PAPI3) = 0._dp

            PRDROG(PLIM1) = 0._dp
            PRDROG(PLIM2) = 0._dp
            PRDROG(PLIM3) = 0._dp

            DO n = 1, LDROG
               VDROG3( i,k,j, n ) =  oconv * PRDROG( n ) * DTSTEPC
               VDROG3( i,k,j,n  ) = MAX( 0., VDROG3( i,k,j, n ) )
            ENDDO
        endif
