REAL(KIND=dp) :: roleho,roleo3,roleno3,rtolho,    &
                 rparho,rivocho,risoho,risoo3,risono3,           &
                 rc10h16ho,rc10h16o3,rc10h16no3,                  &
                 rsesqho,rsesqo3,rsesqno3

REAL(KIND=dp) , DIMENSION(ldrog_soa2)  :: PRDROG

REAL(kind=dp) :: es, qvs, rh

REAL(KIND=dp) :: RATE_CVPOA
REAL(KIND=dp) :: RATE_CVASOA
REAL(KIND=dp) :: RATE_CVASOA0
REAL(KIND=dp) :: RATE_CVBSOA

INTEGER, SAVE :: icall

RATE_CVPOA = REAL( config_flags%rate_cvpoa, KIND=dp )
RATE_CVASOA= REAL( config_flags%rate_cvasoa, KIND=dp )
RATE_CVASOA0= REAL( config_flags%rate_cvasoa0, KIND=dp )
RATE_CVBSOA= REAL( config_flags%rate_cvbsoa, KIND=dp )
