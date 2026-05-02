      program smulti(IFOLD)
      implicit real(a-h,o-z)
      PARAMETER (LM=06)
      PARAMETER (LMX=3*LM)
      PARAMETER (MAX_ANG=((LMX+1)*(LMX+2))/2)
      PARAMETER (MPX=MAX_ANG)
      PARAMETER (MXPR=MXPOISS)
      DIMENSION POLY(MPX,(LM+1)**2),SMULTI((LM+1)**2)


        PRINT *,'CALLING ANGMSH max_ang',MAX_ANG
        LMAX = LMX
        CALL ANGMSH(MAX_ANG,LMAX,NANG,ANG,DOMEGANG)
        CALL flush(6)
        PRINT *, 'NANG = ',NANG

        rad_sph = RSPHERE+20.0D0

        DO ipts = 1,NANG
          DO ix=1,3
            RMSH(ix,JMSH+1)=rad_sph*ANG(ix,ipts)
          ENDDO
          WMSH(JMSH+1)=0.0D0
          JMSH=JMSH+1
        ENDDO

        NMSH=JMSH

        print *,'rad_sph for SMULTI:',rad_sph
        print *,'NMSH, MMSH1, and NDCELL=',NMSH,MMSH1,NDCELL      

        CALL HARMONICS(MPX,NANG,LM,ANG,POLY,NPOLY)
!
! Checking non-orthonormality ERROR
!
        TOLER = 1.0D-4
        ERROR = 1.0d-30

        DO I=1,NPOLY
         DO J=I,NPOLY
           SSSS=0.0D0
           DO IANG=1,NANG
             SSSS=SSSS+POLY(IANG,I)*POLY(IANG,J)*DOMEGANG(IANG)
           ENDDO
           IF(I.NE.J)THEN
            ERROR=MAX(ERROR,ABS(SSSS))
            IF((ABS(SSSS).GT.TOLER))PRINT *,I,J,SSSS
            ELSE
            !PRINT *,I,J,SSSS
           ENDIF
         ENDDO
        ENDDO
        SOLANG=0.0D0
        DO IANG=1,NANG
          SOLANG=SOLANG+DOMEGANG(IANG)
        ENDDO
        !PRINT *,'SOLANG= ',SOLANG

        OPEN(199,FILE='SMULTI',form='formatted',status='unknown')
        REWIND(199)
     
        icount2=0
        DO I=1,((LM+1)**2)
          SMULTI(I)=0.0D0
          kk=1
          DO IPTS=MMSH1+NATMS+1,NMSH
            ZPOTTEMP=0.0D0
            DO ICNT=MMSH1+1,MMSH1+NATMS
              DIST=0.0D0
              DO j=1,3
                DIST=DIST+(RMSH(j,ICNT)-RMSH(j,IPTS))**2
              ENDDO
              DIST=1.0D0/DSQRT(DIST)
              ZPOTTEMP=ZPOTTEMP-ZATMS(ICNT-MMSH1)*DIST
            ENDDO
            !print *,'ZPOTTEMP',ZPOTTEMP

! COULOMB here is evaluated using the Wannier functions of central cell.            
            SMULTI(I)=SMULTI(I)+POLY(kk,I)*DOMEGANG(kk)
     &             *(COULOMB(IPTS)+ZPOTTEMP)
             kk=kk+1
          ENDDO
          ll=ceiling(dsqrt(dfloat(I))-1.0D0)
          SMULTI(I)=SMULTI(I)*(rad_sph**dfloat(ll+1))
          smulti_toll=1.0D-18 !Define
          IF(abs(SMULTI(I)).GT.smulti_toll)THEN
          icount2=icount2+1
          ENDIF
        END DO
        
        print *,'smulti_toll:',smulti_toll

        WRITE(199,*)icount2
        DO I=1,((LM+1)**2)
          IF(abs(SMULTI(I)).GT.smulti_toll)THEN
          WRITE(199,*)I,SMULTI(I)
          ENDIF
        ENDDO
        CLOSE(199)

      END



