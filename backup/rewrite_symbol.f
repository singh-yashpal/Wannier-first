        subroutine rewrite_symbol(ncalc,mode_run,imesh,ierr)
        implicit real*8 (a-h,o-z)
        !include 'commons.inc'
        character*100 line,mine,kine,jine(100)
        dimension prim(3,100),grp(3,3,121)
        dimension all_atoms(4,1000),v(3)
        logical exist
        dimension r(3,1000),dist2(4,1000),temp2(4)
        dimension inx(3),nxx(1000),nyy(1000),nzz(1000)
        dimension dist(1000),distt(1000),rcut2(1000)
        dimension tran(3),rtm(3,3)
        dimension ncon1(3),ncon2(3),mcon(3)
        dimension alpha(100),coef(100)
        dimension radmax1(10),radmin1(10),natmrad(10)
        dimension radmax2(10),radmin2(10),natmrad2(10)
        common/PTRANS/BBOUND(3,2),TVEC(3,3),TWIST(3,3,3),
     &       ROT(3,3),RSPHBAS,RSPHWD,NXYZ

!         COMMON/BASET/ZELC(MAX_FUSET),ZNUC(MAX_FUSET)
!     &  ,BFCON(MAX_BARE,MAX_CON,LDIM,MAX_FUSET)
!     &  ,BFALP(MAX_BARE,MAX_FUSET) 
!     &  ,N_BARE(MAX_FUSET),N_CON(LDIM,MAX_FUSET)
!     &  ,LSYMMAX(MAX_FUSET),N_POS(MAX_FUSET),NFNCT

        call system('rm SMULTI*')
        tvec = 0.0d0
        open(90,file='PTRANSLATIONS')
        do ii=1,3
        read(90,*)BBOUND(ii,1),BBOUND(ii,2)
        enddo
        read(90,*)NXYZ
        do iper=1,NXYZ
        read(90,*)(tvec(j,iper),j=1,3)
        !read(90,*)((twist(j,i,iper),j=1,3),i=1,3)   
        end do
        read(90,*)ucut
        read(90,*)twistang,twistax
        close(90)

        INQUIRE(FILE='ISYMTOT1',EXIST=EXIST)
        IF (EXIST) GOTO 999

        call choose_radius(ucut,radmax)
        print *,'after ucut'
        rsphbas = radmax!+0.2d0*radmax
        rsphwd  = radmax !+0.4d0*radmax
        print *,'rsphbas rsphwd:',rsphbas,rsphwd

        call system('rm TVECGEN TRANSLATIONS')
        call system('rm INRED')
        call system('rm CLUSTER-TUBE')
        call system('rm UTRANS ALL_CELLS')

        if(NXYZ.EQ.0)return
        call system('rm SYMPERTOT')
        call system('grep "ALL-" SYMBOL > CELL')
        call system('cat ISYMGEN > ISYMTOT')
        call system('XMOL_PER.DAT')
        open(60,file='SYMPERHLP')
        open(70,file='ALL_CELLS')
        open(71,file='UTRANS')
        open(80,file='XMOL_PER.DAT')
        open(90,file='CELL')

        ncel=0
        natm=0
        nnn=0

 100    continue     
        read(90,'(A)',end=200)line
        natm=natm+1
            read(line(15:100),*)(r(j,natm),j=1,3)
            print*,natm,line(5:7),(r(j,natm),j=1,3)
        goto 100
 200    continue
        close(90)

       do iatm=1,natm
          do j=1,3
          prim(j,iatm)=r(j,iatm)
          enddo
       end do

       open(90,file='GRPMAT')
       read(90,*)ngp
       read(90,*)(((grp(j,i,igp),j=1,3),i=1,3),igp=1,ngp)

       jpr=natm
       do itry=1,3
         kpr=jpr
         do ipr=1,kpr !loop over natm
         do igp=1,ngp

           do j=1,3
             prim(j,jpr+1)=0.0d0
           do k=1,3
             prim(j,jpr+1)=prim(j,jpr+1)+grp(j,k,igp)*prim(k,ipr)
           end do
           end do

           errmin=1.0d30
           do lpr=1,jpr
             err=abs(prim(1,jpr+1)-prim(1,lpr))
     &          +abs(prim(2,jpr+1)-prim(2,lpr))
     &          +abs(prim(3,jpr+1)-prim(3,lpr))
             errmin=min(err,errmin)
           end do
           if(errmin.gt.0.01)jpr=jpr+1
         end do ! ngp 
         end do !loop over natm
       end do

       do ipr=1,kpr
         print '(3F12.6)',(prim(j,ipr),j=1,3)
       end do
       close(90)

       inx(:)=0
       all_atoms(:,:)=0.0d0
       err=1.0d30

       if(NXYZ.ge.1)inx(1)=30
       if(NXYZ.ge.2)inx(2)=30
       if(NXYZ.eq.3)inx(3)=30
       maxcel = (2*inx(1)+1)*(2*inx(2)+1)*(2*inx(3)+1)
       ins=0
       do nx=-inx(1),inx(1)
       do ny=-inx(2),inx(2)
       do nz=-inx(3),inx(3)
         if(nx.eq.0.and.ny.eq.0.and.nz.eq.0)goto 499
         nin=0
         do ipr=1,kpr
           dst=0.0d0
           do j=1,3
           dst=dst+(prim(j,ipr)+
     &     tvec(j,1)*nx+tvec(j,2)*ny+tvec(j,3)*nz)**2
           end do
           dst=dsqrt(dst)
           if(dst.le.rsphbas)nin=nin+1
         if(dst.gt.rsphbas.and.dst.lt.rsphwd)then
         write(72,'(3I6)')nx,ny,nz
         endif
         enddo
         if(nin.eq.0)goto 499
         write(71,'(3I6)')nx,ny,nz
         do ipr=1,kpr
           ins=ins+1
           dst=0.0d0
           do j=1,3
             all_atoms(j,ins)=prim(j,ipr)+
     &       tvec(j,1)*nx+tvec(j,2)*ny+tvec(j,3)*nz
             dst=dst+all_atoms(j,ins)**2
           end do
           dst=dsqrt(dst)
           nxx(ins)=nx
           nyy(ins)=ny
           nzz(ins)=nz
           dist(ins)=dst
!           do jns=1,ins-1
!              err=abs(all_atoms(1,jns)-all_atoms(1,ins))
!     &           +abs(all_atoms(2,jns)-all_atoms(2,ins))
!     &           +abs(all_atoms(3,jns)-all_atoms(3,ins))
!           end do
         end do
  499  continue
       end do
       end do
       end do

        do jns=1,ins
           print '(3F12.6)',(all_atoms(j,jns),j=1,3)
           all_atoms(4,jns)=1.0d0
        end do

c purge sites that are equivalent by symmetry:
        print*,'inequivalent sites:'
        rmaxatm=0.0d0
        mns=0
        do jns=1,ins
           if(all_atoms(4,jns).gt.0.0d0)then
           do kns=jns+1,ins
           do igp=1,ngp
              do j=1,3
                v(j)=0.0d0
              do k=1,3
                v(j)=v(j)+grp(j,k,igp)*all_atoms(k,jns)
              enddo
              enddo
              err=abs(v(1)-all_atoms(1,kns))
     &           +abs(v(2)-all_atoms(2,kns))
     &           +abs(v(3)-all_atoms(3,kns))
              if(err.lt.0.01d0)all_atoms(4,kns)=-1.0d0
            enddo
            enddo
            mns=mns+1
            do j=1,4
              all_atoms(j,mns)=all_atoms(j,jns)
            enddo
            dist(mns)=dist(jns)
            nxx(mns)=nxx(jns)
            nyy(mns)=nyy(jns)
            nzz(mns)=nzz(jns)

            rmaxatm=max(rsphbas,dist(mns))
            !rmaxatm=max(rmaxatm,dist(mns))
            print '(1F12.4,3i4,3F12.6)',dist(mns),
     &                 nxx(mns),nyy(mns),nzz(mns),
     &                   (all_atoms(j,mns),j=1,3)

            print '(i5,3F12.6)',mns,(all_atoms(j,mns),j=1,3)
            write(70,'(F12.4,3I6,3F12.6,I3)')dist(mns),nxx(mns)
     &           ,nyy(mns),nzz(mns),(all_atoms(j,mns),j=1,3)
            write(60,'(F10.3," ",A,3F20.12)')dist(mns),
     &            line(1:13),(all_atoms(j,mns),j=1,3)
            write(80,*)' 6 ',(all_atoms(j,mns)*0.529177,j=1,3) 
           end if
        end do

        ncel=(mns/natm)+1
        print *,'ncel',ncel
       close(80)

       open(80,file='XMOL_PER_TMP')
       write(80,*)ncel
       write(80,*)' ',maxcel
       close(80)

       close(70)
       close(71)

       call system('cat XMOL_PER.DAT >>XMOL_PER_TMP')
       call system('mv XMOL_PER_TMP XMOL_PER.DAT')

       if(ncel.eq.maxcel)then
       print*,'cutoff too small in ...'
       stop        
       end if

       rewind(60)

       open(66,file='SYMPERHLP2')
       natmtot=(ncel-1)*natm
       print *,'ALL_ATMS:',natmtot

       dist2=0.0d0
       rcut2=0.0d0
       do i = 1,natmtot
        read(60,'(F10.3,A)')rcut,line
        dist2(4,i)=rcut
        read(line(15:100),*)(dist2(jj,i),jj=1,3)
       enddo
      ! Sorting the data based on the first column using bubble sort
      do i = 1, natmtot-1
        do j = i+1, natmtot
            if (dist2(4,i) .gt. dist2(4,j)) then
                temp2(:) = dist2(:,i)
                dist2(:,i) = dist2(:,j)
                dist2(:,j) = temp2(:)
            endif
        end do
      end do

      do i = 1,natmtot
            write(66,'(F10.3,"",A,3F20.12)')dist2(4,i),
     &            line(1:13),(dist2(j,i),j=1,3)

        !write(66,*)dist2(4,i),lab2,equ,(dist2(jj,i),jj=1,3)
        !print *,dist2(4,i),(dist2(jj,i),jj=1,3)
      end do
      close(66)
      close(60)

      call system('cp SYMPERHLP2 SYMPERHLP')
!================================================================================              
      open(60,file='SYMPERHLP')
      rewind(60)

      onept5= 1.6d0
      radmax=rmaxatm+0.01D0
      !radmin=radmax/onept5
      radmin=radmax-2.62D0
      print*,'natm,radmin,radmax',natm,radmin,radmax

      kk=0
      do irad=0,9
        rewind(60)
        nn=0
        do i=1,natmtot
           read(60,'(F10.3)')rcut
           if(rcut.ge.radmin.and.rcut.lt.radmax)then
           nn=nn+1
           endif
        enddo
        if(nn.ne.0)then
        kk=kk+1
        radmax2(kk)=radmax
        radmin2(kk)=radmin
        natmrad2(kk)=nn
        endif
        radmax=radmin
        !radmin=radmin/onept5
        if(irad.eq.0)then
        radmin=radmax-2.62D0
        endif
        if(irad.eq.1)then
        radmin=radmax-2.62D0
        endif
        if(irad.gt.1)then
        radmin =0.5d0
        endif
      enddo

      do jj=1,kk
        radmax1(jj)=radmax2(kk+1-jj)
        radmin1(jj)=radmin2(kk+1-jj)
        natmrad(jj)=natmrad2(kk+1-jj)
      enddo

      print *,'radmax(i)',(radmax1(ii),ii=1,kk)
      print *,'radmin(i)',(radmin1(ii),ii=1,kk)
      print *,'natmrad(i)',(natmrad(ii),ii=1,kk)

!================================================================================              
      rewind(60)

      open(67,file='XMOLWD.xyz',status='unknown')
      open(677,file='XMOLWD.DAT',status='unknown')
      write(67,*)natm+natmtot
      write(67,*)'xyz format'
      do ii=1,natm
        write(67,'(A,3F12.6)')'C',(prim(jj,ii)*0.529177,jj=1,3)
        write(677,'(A,3F12.6)')'C',(prim(jj,ii),jj=1,3)
      enddo
      temp2=0.0d0
      do ii=1,kk
        do jj=1,natmrad(ii)
          read(60,'(F10.3,A)')rcut,line
          read(line(15:100),*)(temp2(j),j=1,3)
          if(ii.lt.kk)then
          write(67,'(A,3F12.6)')'N',(temp2(j)*0.529177,j=1,3)
          write(677,'(A,3F12.6)')'C',(temp2(j),j=1,3)
          else
          write(67,'(A,3F12.6)')'O',(temp2(j)*0.529177,j=1,3)
          write(677,'(A,3F12.6)')'C',(temp2(j),j=1,3)
          endif
        enddo
      enddo
      close(67)
!================================================================================              
      open(70,file='SYMPERDEL')
      rewind(60)

      DO irad=1,kk
          n=natmrad(irad)
          print *,'LOOP over irad and N:',irad,n
          rewind(70)
          do i=1,n
            read(60,'(F10.3,A)')rcut,line
            write(70,'(A,I1,A,I3.3," = ",A,"$",F10.3)')'ALL',irad,
     &              line(6:8),i,trim(line(15:100)),rcut
          end do
          close(70)

          call system('cat SYMPERDEL >>SYMPERTOT')

          open(70,file='SYMPERDEL')
          rewind(70)
          rsum=0
          do m=1,n   
             read(70,'(A)')mine
             do l=1,100
               if(mine(l:l).eq.'$')then
               read(mine(l+1:100),*)rnow
               mine(l:l)=' '
               end if
             enddo
             rsum=rsum+rnow
             print*,trim(mine)                      
          end do

          ravg=rsum/n
! Decision barcut: exponent of the longest range gaussian
! Decision barsht: exponent of the shortest range gaussian
          
          !barcut=ucut/(rmaxatm-ravg)**2 !Decision: radmax vs rsphbas
          !barcut=min(barcut,0.4d0) !Decision
          !barcut=ucut/(rsphbas1-ravg)**2 !Decision: radmax vs rsphbas
          barcut=0.14d0 !Decision
          barsht=barcut+0.230D+05 ! Decision 
          open(40,file='ISYMGEN')
            rewind(40)
          print*,'opening ISYMADD',N
          open(45,file='ISYMADD')
**************************
          REWIND(40)
                       NEXEC=0
  10      CONTINUE
          READ(40,'(A)',END=15)KINE
             DO J=1,50
C            PRINT*,KINE(J:J+6)
             IF(KINE(J:J+6).EQ.'ATOMS O')READ(KINE,*)NINEQ     
             IF(KINE(J:J+6).EQ.' CHARGE')READ(KINE,*)NUCE,NUCN
                          NUCE=0
                          NUCN=0
             END DO
          IF(KINE(1:3).EQ.LINE(2:4).AND.
     &       KINE(5:7).EQ.LINE(6:8))THEN
          PRINT*,'KINE:',KINE(1:7),' ',LINE(2:8)
                JINE(1)=KINE
                       NEXEC=1
                DO IINEQU=2,NINEQ
                READ(40,'(A)')JINE(IINEQU)
                END DO
                DO IINEQU=1,NINEQ
                PRINT*,JINE(IINEQU)
                END DO
          READ(40,*)
          READ(40,*)NGAUSS
          READ(40,*)NCON1          
          READ(40,*)NCON2         
          READ(40,*)(ALPHA(IGAUSS),IGAUSS=1,NGAUSS)
              LGMIN=1      
              LGMAX=NGAUSS
              DO IGAUSS=NGAUSS,1,-1 
              IF(ALPHA(IGAUSS).LT.BARCUT)LGMAX=IGAUSS-1
              END DO
              DO IGAUSS=1,NGAUSS
              IF(ALPHA(IGAUSS).GT.BARSHT)LGMIN=IGAUSS+1
              END DO
              IF(LGMIN.GT.LGMAX)LGMAX=MAX(1,LGMIN-1)
              IF(LGMIN.GT.LGMAX)LGMAX=LGMIN+2          
              LGAUSS=LGMAX-LGMIN+1
              PRINT*,LGMIN,LGMAX,LGAUSS
          MCON=0
          WRITE(45,*)NUCE,NUCN,' ELECTRONIC AND NUCLEAR CHARGE'
          WRITE(45,'(A)')'ALL'
          !write(45,*)'RADMAX AND RAVG',radmax,ravg
          WRITE(45,*)N,' NUMBER OF ATOMS OF TYPE:',MINE(1:7)
          NNN=NNN+N
                    DO I=1,N
C                   PRINT*,MINE(1:7),I
                    WRITE(45,'(A7,I3.3)')MINE(1:7),I
                    END DO
          WRITE(45,'(A)')'EXTRABASIS'
          
          igaus=0
          if(irad.eq.kk-3)igaus=0
          if(irad.eq.kk-2)igaus=0
          if(irad.eq.kk-1)igaus=0
          if(irad.eq.kk)igaus=0
!              igaus=5
!              !if(irad.eq.kk)igaus=3
!          WRITE(45,'(I5,A)')LGAUSS-igaus, ' NUMBER OF BARE GAUSSIANS'
!          WRITE(45,'(A)')'PLACEHOLDER'
!          WRITE(45,'(5I5," !NS, NP, ND, NF, NG - EXTRA")')(0,LL=1,5)
!          WRITE(45,'(3G20.10)')(ALPHA(IG),IG=LGMIN,LGMAX-igaus)
!          else
          WRITE(45,'(I5,A)')LGAUSS-igaus, ' NUMBER OF BARE GAUSSIANS'
          WRITE(45,'(A)')'PLACEHOLDER'
          WRITE(45,'(5I5," !NS, NP, ND, NF, NG - EXTRA")')(0,LL=1,5)
          WRITE(45,'(3G20.10)')(ALPHA(IG),IG=LGMIN,LGMAX-igaus)
          WRITE(45,*)' '
!          endif

          DO L=0,2
          DO ICON=1,NCON1(L+1)+NCON2(L+1)
          READ(40,*)(COEF(IG),IG=1,NGAUSS)
          SUM1=-1.0D0
          SUM2= 1.0D0-NGAUSS
              DO IGAUSS=1,NGAUSS
              SUM1=SUM1+ABS(COEF(IGAUSS))
              SUM2=SUM2+ABS(COEF(IGAUSS)-1)
              END DO
          IF(ABS(SUM1).LT.0.001D0.AND.ABS(SUM2).LT.0.001D0)THEN
           ELSE
           if((MCON(L+1)+1).LE.LGAUSS)THEN
           MCON(L+1)=MCON(L+1)+1
           NCOEF=(LGMAX+1-LGMIN)
           !print *,'NCOEF AND LGAUSS',NCOEF,LGAUSS
           !if(NCOEF.NE.LGAUSS)then
           !        PRINT *,'NCOEF DOESNT MATCH LGAUSS',NCOEF,LGAUSS
           !  CALL STOPIT
           !endif

          igaus=0
          if(irad.eq.kk-3)igaus=0
          if(irad.eq.kk-2)igaus=0
          if(irad.eq.kk-1)igaus=0
          if(irad.eq.kk)igaus=0
!           if(irad.gt.kk-1)then
!              igaus=5
!              !if(irad.eq.kk)igaus=3
!           WRITE(45,'(3G20.10)') (COEF(IG),IG=LGMIN,LGMAX-igaus)
!           else
           WRITE(45,'(3G20.10)') (COEF(IG),IG=LGMIN,LGMAX-igaus)
!           endif
           !print *,'Hi Hi',(COEF(IG),IG=LGMIN,LGMAX)
           WRITE(45,*)' '
          END IF
          END IF

          END DO
          END DO
                GO TO 15
          ELSE
          GO TO 10
          END IF
 15       CONTINUE

          IF(NEXEC.EQ.0)STOP'NEXEC'
          CLOSE(40)
C         PRINT*,MCON,' NS, NP, ND'
C         PRINT*,' 0 0 0'
                   CLOSE(45)
                   OPEN(45,FILE='ISYMADD')
          REWIND(45)
          OPEN(46,FILE='DELISYMADD')
                 DO MM=1,100000
                 READ(45,'(A)',END=20)KINE
                 IF(KINE(1:11).NE.'PLACEHOLDER')THEN
                 WRITE(46,'(A)')KINE
                 ELSE
              WRITE(46,'(3I5,A)')MCON,' 0 0  NS, NP, ND, NF, NG'
                 END IF
                 END DO
 20       CONTINUE
          CLOSE(45)  
          CLOSE(46)  
          call system('mv DELISYMADD ISYMADD')
          call system('cat ISYMADD >> ISYMTOT')
          call system('rm    ISYMADD')
        ENDDO !irad

        close(60)
        close(70)
        call system(' rm SYMPERHLP')
        call system('grep "NUMBER OF ATOMS" ISYMTOT > SYMPERHLP')
          open(20,file='SYMPERHLP')
          open(25,file='ISYMTOT')
          READ(25,'(A)')LINE
          DO 300 I=1,10000
          READ(20,'(A)',END=210)LINE
          M=I
          PRINT*,LINE
  300     CONTINUE
  210     CONTINUE
          REWIND(20)
          WRITE(20,*)M,' !NUMBER OF FUNCTION SETS'
  220     CONTINUE
          READ(25,'(A)',END=225)LINE
          IF(LINE(1:5).NE.'WFOUT'.AND.LINE(1:9).NE.'ELECTRONS')THEN
          WRITE(20,'(A)')LINE
          END IF
          GO TO 220
  225     CONTINUE
          WRITE(20,'(A)')'ELECTRONS'
          WRITE(20,'(A)')'WFOUT'
          CLOSE(25)
          CLOSE(20)
          call system('mv SYMPERHLP ISYMTOT')
          open(20,file='SYMBOL')
          open(25,file='SYMHLP')
          DO I=1,3
          READ (20,'(A)')LINE
          WRITE(25,'(A)')TRIM(LINE)
          END DO
          WRITE(25,'(A)')' 2 NUMBER OF SYMBOLIC FILES'
          WRITE(25,'(A)')"ISYMTOT = INPUT"
          WRITE(25,'(A)')"TVECGEN = TRANSLATIONS"
               DO I=1,3
               READ(20,*)
               END DO
          READ(20,*)NSYMBOLS
          READ(20,*)NUCLEI
               MUCLEI=NUCLEI
               NUCLEI=NUCLEI+NNN
               NSYMBOLS=NSYMBOLS+NNN   
          WRITE(25,*)NSYMBOLS,' NSYMBOLS'
          WRITE(25,*)NUCLEI,' NUCLEI'
               DO I=1,NUCLEI
               WRITE(25,'(A)')'1.0 1 1 1'                     
               END DO
          WRITE(25,*)' 1 THIS IS ALWAYS THE FIRST CALCULATION'
          DO I=1,MUCLEI
          READ(20,*)
          END DO
          READ(20,*)
          DO I=1,MUCLEI
          READ (20,'(A)')LINE
          WRITE(25,'(A)')TRIM(LINE)
          END DO
          CLOSE(20)

!          open(66,file='SYMPERHLP2')
!          rewind(66)
!          do irad==1,kk
!            do I=1,natmrad(irad)
!             
!            enddo
!          enddo
!

          OPEN(20,FILE='SYMPERTOT')
          DO I=1,NNN          
          READ (20,'(A)')LINE
               K=100
               DO J=1,100
               IF(LINE(J:J).EQ.'$')K=J-1
               END DO
          WRITE(25,'(A)')TRIM(LINE(1:K))
          END DO
          CLOSE(20)
          CLOSE(25)
          call system('grep ELECTRONS SYMBOL |tail -1 >>SYMHLP')
          call system('grep EXTRA     SYMBOL |tail -1 >>SYMHLP')
          call system('mv SYMHLP SYMTOT')
        ierr=-1
        call symbol(mcalc,mode_run,imesh,ierr)
  999   CONTINUE
        print*,'END OF REWRITE_SYMBOL'
        call system('cp ISYMTOT1 ISYMTOT')
        call system('cp INPUT1 INPUT')
        call system('cp SYMTOT1 SYMTOT')
        call stopit
        end

        subroutine rottot(nx,twist,rot)
        dimension twist(3,3),rot(3,3),prd(3,3),rpm(3,3)
        rot=0.0d0
        prd=0.0d0
        do i=1,3
        rot(i,i)=1.0d0
        prd(i,i)=1.0d0
        end do
        if(nx.eq.0)return
        if(nx.gt.0)then
         rpm=twist
        else
          do i=1,3
          do j=1,3
          rpm(i,j)=twist(j,i)
          end do
          end do
        end if
        mx=abs(nx)
        rot=0.0d0
        do i=1,3
        rot(i,i)=1.0d0
        end do
        do ix=1,mx
         do i=1,3
         do j=1,3
         prd(i,j)=0.0d0
           do k=1,3
           prd(i,j)=prd(i,j)+rot(i,k)*rpm(k,j)
           end do
         end do
c        print '(3F12.4)',(prd(i,j),j=1,3)
         end do
         rot=prd
        end do
        return
        end
