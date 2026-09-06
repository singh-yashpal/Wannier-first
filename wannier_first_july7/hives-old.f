       include 'PARAMS' 
       include 'commons.inc'
       nproc=18-1
       irank=0
       call manage_hives(0,msplit)
c print*,'what is msplit?'
       read*,msplit
       call manage_hives(1,msplit)
       !call manage_hives(2,msplit)
       end 
       subroutine manage_hives(mode,nsplit)
       include 'PARAMS' 
       include 'commons.inc'
       common/hives/nqn,myqn(0:mx_proc),nbees(0:mx_proc),
     &         my_bees(0:mx_proc,0:mx_proc)
            msplit=min(nbees(irank),nsplit)
c mode=0 initialization
       if(mode.eq.0)then
           nbees=0
           nqn=0
           myqn=0
           nbees(0)=nproc
           do iproc=0,nbees(irank)
            my_bees(iproc,irank)=iproc
           end do
        print'(40I3)',(my_bees(iproc,irank),iproc=1,nproc)
        end if
        if(mode.eq.1)then
c mode=1: create mpslit new queens, notify new queens,notify bees of their new queen
        print*,' I am princess:',irank
        print*,' My queen is:',myqn(irank) 
        print "(' My bees:',40i3)",(my_bees(j,irank),j=0,nbees(irank))
        print*,' I am promoting',msplit,' new princesses'
              new=nbees(irank)/msplit
              nrm=nbees(irank)-(nbees(irank)/msplit)*msplit
              mew=new+1
              iqn=nqn+1
                 ibee=0
c                 print*,my_bees(jbee,irank)
c             nbees(nqn+msplit+1)=msplit
              isplit=0
              do jbee=1,nbees(irank)           
                  kbee=my_bees(jbee,irank)
                  ibee=ibee+1
                  if(ibee.le.mew)then
                     nbees(jrank)=mew
                     myqn(jrank)=irank
                  else
                  ibee=1
                  iqn=iqn+1
                   print*,' '
                  end if
                     if(ibee.eq.1)then
                        jrank=my_bees(kbee,irank)
                        isplit=isplit+1
                     end if
                  print*,'bee:',kbee,' is owned by:',jrank
                  my_bees(ibee ,jrank)=kbee
                             if(ibee.ne.1)myqn(kbee)=jrank
                  my_bees(isplit,40)=jrank
                  if(iqn-nqn.gt.nrm)mew=new
              end do
                 nbees(irank)=msplit
                 do isplit=1,msplit
                 my_bees(isplit,irank)=my_bees(isplit,40)
                 end do
                 nqn=iqn
              do jrank=0,nproc
              print"('My Boss:',I4,3I3,'   ',20I3)",myqn(jrank),
     & nrm,jrank,nbees(jrank),(my_bees(ibee,jrank),ibee=1,nbees(jrank))
              end do
        end if
c mode=2: return bees to former queen
         return
         end
