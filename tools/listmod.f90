!> \file listmod.f90

! AF 2003-12-17
! Write an extra 0 at the end of the first line to make the output with the V.2005 code.


!   Copyright 2002-2010 AstroFloyd - astrofloyd.org
!   
!   
!   This file is part of the eggleton-tools package.
!   
!   This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!   
!   This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License along with this code.  If not, see <http://www.gnu.org/licenses/>.


!> Reads an input or output structure model file for Eggeltons TWIN code and lists the properties of each model it contains.
!! One can then select a model to display its contents more precisely and optionally copy the model to a different file to serve as input. 
!<
program listmod
  use kinds
  use constants
  implicit none
  integer :: narg,iargc,blk,ans,nblk
  character :: fname*99,findfile*99

  
  call setconstants()
  
  
  !***   READ COMMAND LINE VARIABLES
  narg = iargc()
  if(narg.gt.0) then
     call getarg(1,fname)
  else
     write(6,'(A)')'  listmod: lists the contents of a mod-file to screen'
     write(6,'(A)')'           syntax:  listmod <filename>'
     write(6,*)''
     write(6,'(A)')"  I'll look in the current directory for a *.mod file..."
     fname = findfile('*.mod')
  end if
  write(6,'(/,A,/)')'  Reading file '//trim(fname)
  
  
  
  
  !***   LIST ALL STRUCTURE MODELS IN THE FILE AND THEIR MAIN PROPERTIES
3 continue
  call list_mod_file(fname,nblk)
  
  
  
  
  
  !***   CHOOSE STRUCTURE MODEL AND PRINT DETAILS

49 continue  
  blk = 0
  do while(blk.lt.1.or.blk.gt.nblk)
     write(6,'(A,I5,A4,$)')'  For which model do you want to print details (1 -',nblk,'):  '
     read*,blk
     if(blk.eq.0) then
        write(6,*)''
        stop
     end if
  end do
  
  call print_mod_details(fname,blk)
  
  
  !***   FINISH
  
  ans = -1
  do while(ans.gt.3.or.ans.lt.0)
     if(nblk.eq.1) then
        write(6,*)''
        stop
     end if
     
     write(6,*)''
     write(6,'(A)')'  You can:'
     write(6,'(A)')'    0) Quit'
     write(6,'(A)')'    1) Choose another model'
     write(6,'(A)')'    2) List all models again and choose another model'
     write(6,'(A)')'    3) Write this model to another file'
     write(6,*)''
     write(6,'(A28,$)')'  What do you want to do ?  '
     read*,ans
  end do
  
  if(ans.eq.0)  then
     close(10)
     write(6,*)''
     stop
  end if
  if(ans.eq.1) goto 49
  if(ans.eq.2) then
     close(10)
     goto 3
  end if
  
  
  
  
  
  !***   COPY MODEL TO DIFFERENT FILE
  call copy_mod(fname,blk)
  
  
  
  write(6,*)''
end program listmod
!***********************************************************************************************************************************




!***********************************************************************************************************************************
subroutine error_reading_header(b)
  implicit none
  integer, intent(in) :: b
  if(b.eq.0) then
     write(6,'(A,/)')'  Error reading header line, aborting.'
  else
     write(6,'(A,I5,A,/)')'  Error reading header line in block',b,', aborting.'
  end if
  close(10)
  stop
end subroutine error_reading_header
!***********************************************************************************************************************************


!***********************************************************************************************************************************
subroutine error_reading_block(l)
  implicit none
  integer, intent(in) :: l
  if(l.eq.0) then
     write(6,'(A,/)')'  Error reading (the first?) line of the data block, aborting.'
  else
     write(6,'(A,I5,A,/)')'  Error reading the data block, line',l,', aborting.'
  end if
  close(10)
  stop
end subroutine error_reading_block
!***********************************************************************************************************************************



!***********************************************************************************************************************************
subroutine list_mod_file(fname, nblk)
  use kinds
  use constants
  
  implicit none
  character,intent(in) :: fname*(*)
  integer, intent(out) :: nblk
  
  real(double) :: m1,dt,t,p,bms,ecc,p1,enc,horb
  real(double) :: lnf,lnt,x16,lnm,x1,dqdk,lnr,l,x4,x12,x20
  real(double) :: mhe,mco,mi,pr,phi,phis,e,f
  real(double) :: r1,l1,ts,hs,hes,zs,tc,hc,hec,zc
  real(double) :: dat(99)
  integer :: kh,kp,jmod,jb,jin,io
  integer :: bl,li
  
  open (unit=10,form='formatted',status='old',file=trim(fname))
!3 continue
  
  write(6,'(A)')'  Nr  Model Nmsh          Age       dT        M1    Mhe    Mco     Menv         R        L     Teff      Xs'// &
       '     Ys     Zs         Tc      Xc     Yc     Zc      Mtot     Porb     Prot'
  
  bl = 1  !Block/model number
  do
     read(10,*,iostat=io) m1,dt,t,p,bms,ecc,p1,enc,kh,kp,jmod,jb,jin
     if(io.lt.0) exit
     if(io.gt.0) then  !then you may be saving DH as well as H, if this is block 2
        write(0,'(A,I5,A)')'  Error reading the header line of block',bl,'.  Skipping the rest of the file.'
        exit
     end if
     
     mhe = 0.d0
     mco = 0.d0
     do li=1,kh  !Line/mesh point in the current block/model
        read(10,*,iostat=io)dat(1:jin)
        if(io.lt.0) exit
        if(io.gt.0) then
           write(0,'(A,2(I5,A))')'  Error reading line',li,' of block',bl,'.  Skipping the rest of the file.'
           exit
        end if
        
        lnf = dat(1)
        lnt = dat(2)
        x16 = dat(3)
        lnm = dat(4)
        x1 = dat(5)
        dqdk = dat(6)
        lnr = dat(7)
        l = dat(8)
        x4 = dat(9)
        x12 = dat(10)
        x20 = dat(11)
        mi = dat(12)
        pr = dat(13)
        phi = dat(14)
        phis = dat(15)
        horb = dat(17)
        e = dat(18)
        f = dat(19)
        
        if(li.eq.1) then
           r1  = exp(lnr)*1.e11/r0
           l1  = l*1.d33/l0
           ts  = exp(lnt)
           m1 = lnm*1.d33/m0
           hs  = x1
           hes = x4
           zs  = 1.d0 - x1 - x4
        end if
        if(li.eq.kh) then
           tc  = exp(lnt)
           hc  = x1
           hec = x4
           zc  = 1.d0 - x1 - x4
        end if
        if(mhe.eq.0.0.and.x1.lt.0.1) mhe = lnm*1.d33/m0
        if(mco.eq.0.0.and.x4.lt.0.1) mco = lnm*1.d33/m0
     end do !li
     
     if(mod(bl,50).eq.0) write(6,'(/,A)')'  Nr  Model Nmsh          Age       dT        M1    Mhe    Mco     Menv         R'// &
          '        L     Teff      Xs     Ys     Zs         Tc      Xc     Yc     Zc      Mtot     Porb     Prot'
     write(6,'(I4,I7,I5, ES13.5,ES9.2, F10.4,2F7.3,ES9.2, 1x,3ES9.2,1x,3F7.4, 2x,ES9.2,1x,3f7.4,1x,3ES9.2)') &
          bl,jmod,kh,t,dt,m1,mhe,mco,m1-mhe,r1,l1,ts,hs,hes,zs,tc,hc,hec,zc,bms,p,p1
     
     bl = bl+1
     
  end do 
  
  close(10)
  
  write(6,'(A)')'  Nr  Model Nmsh          Age       dT        M1    Mhe    Mco     Menv         R        L     Teff      Xs'// &
       '     Ys     Zs         Tc      Xc     Yc     Zc      Mtot     Porb     Prot'
  
  nblk = bl-1
  write(6,'(I5,A,/)')nblk,' blocks read.'
  
  if(nblk.eq.0) stop
  
end subroutine list_mod_file
!***********************************************************************************************************************************



!***********************************************************************************************************************************
subroutine print_mod_details(fname,blk)
  use kinds
  use constants
  
  implicit none
  character, intent(in) :: fname*(*)
  integer, intent(in) :: blk
  
  real(double) :: m1,dt,t,p,bms,ecc,p1,enc,horb
  real(double) :: lnf,lnt,x16,lnm,x1,dqdk,lnr,l,x4,x12,x20
  real(double) :: mi,pr,phi,phis,e,f
  real(double) :: m2,q1,q2,a,a1,a2,rl1,rl2,x
  real(double) :: r1,l1,ts,hs,hes,zs,cs,os,nes,tc,hc,hec,cc,oc,nec,zc
  real(double) :: mhe,mco,mhenv
  integer :: kh,kp,jmod,jb,jin,io
  integer :: bl,li
  
  
  !Read file, upto chosen model (blk-1)
  open (unit=10,form='formatted',status='old',file=trim(fname))
  do bl=1,blk-1  !Block/model
     read(10,*,iostat=io)m1,dt,t,p,bms,ecc,p1,enc,kh,kp,jmod,jb,jin
     if(io.ne.0) call error_reading_header(bl)
     do li=1,kh  !Line/mesh point in model
        read(10,*,iostat=io)lnf,lnt,x16,lnm,x1,dqdk,lnr,l,x4,x12,x20,mi,pr,phi,phis,x,horb,e,f,x,x,x,x,x
        if(io.ne.0) call error_reading_block(li)
     end do !li
  end do !bl
  
  
  
  
  !***   READ CHOSEN STRUCTURE MODEL AND GET VARIABLES TO PRINT
  
  read(10,*,iostat=io)m1,dt,t,p,bms,ecc,p1,enc,kh,kp,jmod,jb,jin   !jin = # columns
  if(io.ne.0) call error_reading_header(0)
  read(10,*,iostat=io)lnf,lnt,x16,lnm,x1,dqdk,lnr,l,x4,x12,x20,mi,pr,phi,phis,x,horb,e,f,x,x,x,x,x
  if(io.ne.0) call error_reading_block(0)
  
  m1  = lnm*1.d33/m0
  r1  = exp(lnr)*1.e11/r0
  l1  = l*1.d33/l0
  !l1  = l/3.844d0  !Peter's CLSN
  ts  = exp(lnt)
  hs  = x1
  hes = x4
  cs  = x12
  os  = x16
  nes = x20
  zs  = 1.d0 - hs - hes
  
  mhe = 0.d0
  mco = 0.d0
  do li=1,kh-1 !Number of Mesh points
     read(10,*,iostat=io)lnf,lnt,x16,lnm,x1,dqdk,lnr,l,x4,x12,x20,mi,pr,phi,phis,x,horb,e,f,x,x,x,x,x
     if(io.ne.0) call error_reading_block(li)
     if(mhe.eq.0.0.and.x1.lt.0.1) mhe = lnm*1.d33/m0
     if(mco.eq.0.0.and.x4.lt.0.1) mco = lnm*1.d33/m0
  end do
  
  close(10)
  
  
  mhenv = m1 - mhe
  
  tc  = exp(lnt)
  hc  = x1
  hec = x4
  cc  = x12
  oc  = x16
  nec = x20
  zc  = 1.d0 - hc - hec
  
  m2  = bms - m1
  q1  = m1/m2
  q2  = m2/m1
  
  a   = (p*day/((g*bms*m0)/(4.d0*pi**2))**(-.5d0))**(2.d0/3.d0)/r0
  
  a1  = a *m2/bms
  a2  = a *m1/bms
  
  rl1 = a*(0.49*q1**(2.d0/3.)/(0.6d0*q1**(2.d0/3.)+ dlog(1.d0+q1**(1.d0/3.))))
  rl2 = a*(0.49*q2**(2.d0/3.)/(0.6d0*q2**(2.d0/3.)+ dlog(1.d0+q2**(1.d0/3.))))
  
  
  
  

  
  
  
  
  !***   PRINT MODEL DETAILS
  
  write(6,'(A)')'  Properties of this model:'
  write(6,*)''
  write(6,81)jmod,m1,t,dt,zs
  write(6,82)kh,kp,jin,jb
  write(6,*)''
  write(6,83)m1,r1,l1,tc,ts
  write(6,84)mhe,mco,mhenv
  write(6,*)''
  write(6,85)m1,m2,bms,q1,q2
  write(6,*)''
  write(6,86)p,a,a1,a2,rl1,rl2
  write(6,87)horb*1.d50,e,pr
  write(6,*)''
  write(6,88)hs,hes,cs,os,nes,zs
  write(6,89)hc,hec,cc,oc,nec,zc
  
81 format ('  Model:        Model nr:',i5,',    Mass:',f7.2,' Mo,    Age: ',es10.4,' yr,  Time step:   ',es10.4,' yr,    Z ='f7.4)
82 format ('                Mesh pts: ',i4,',      Kp: ',i6,',       Jin: ',i3',            Binary component:',i2)
83 format ('  Primary:      M   =',f9.5,' Mo,  R   =',f9.5,' Ro,   L   =  ',es10.4' Lo,  Tc =  ',es10.4,' K,   Teff =',f8.0,' K')
84 format ('                Mhe =',f9.5,' Mo,  Mco =',f9.5,' Mo,   Menv=',f9.5,' Mo')
85 format ('  Binary:       M1  =',f9.5,' Mo,  M2  =',f9.5,' Mo,   Mb  =',f9.4' Mo,     q1 =',f9.5,',        q2   =',f9.5)
86 format ('  Orbit:        P  =',ES12.5,' d,   a  =',ES12.5,' Ro,   a1 =',ES12.5,' Ro,    a2 =',ES12.5,' Ro,   Rrl1 =',ES12.5, &
        ' Ro,   Rrl2 = ',ES12.5,' Ro')
87 format ('                J  =  ',es10.4,' erg s,                e  =',f9.5,',     Prot =',ES12.5,' days')
88 format ('  Composition:  Surface:  H: ',f6.4,',   He: ',f6.4,',   C: ',f6.4,',   O: ',f6.4,',   Ne: ',f6.4,',    Z: ',f6.4)
89 format ('  Composition:     Core:  H: ',f6.4,',   He: ',f6.4,',   C: ',f6.4,',   O: ',f6.4,',   Ne: ',f6.4,',    Z: ',f6.4)
  
  
  
  
  
end subroutine print_mod_details
!***********************************************************************************************************************************



!***********************************************************************************************************************************
subroutine copy_mod(fname,blk)
  use kinds
  implicit none
  character, intent(in) :: fname*(*)
  integer, intent(in) :: blk
  
  real(double) :: m1,dt,t,p,bms,ecc,p1,enc,horb
  real(double) :: lnf,lnt,x16,lnm,x1,dqdk,lnr,l,x4,x12,x20
  real(double) :: mi,pr,phi,phis,e,f,x
  real(double) :: dat1(8),dat2(24)
  integer :: kh,kp,jmod,jb,jin,io
  integer :: bl,li
  character :: outname*99
  logical :: ex
  
  !Read blocks before the desired one:
  open (unit=10,form='formatted',status='old',file=trim(fname))
  do bl=1,blk-1  !Block/model number
     read(10,*,iostat=io)m1,dt,t,p,bms,ecc,p1,enc,kh,kp,jmod,jb,jin
     if(io.ne.0) call error_reading_header(bl)
     do li=1,kh  !Line/mesh point in model
        read(10,*,iostat=io)lnf,lnt,x16,lnm,x1,dqdk,lnr,l,x4,x12,x20,mi,pr,phi,phis,x,horb,e,f,x,x,x,x,x
        if(io.ne.0) call error_reading_block(li)
     end do !li
  end do !bl
  
  
  !Read desired block:
  
  !Read header line:
  read(10,*,iostat=io)dat1,kh,kp,jmod,jb,jin
  if(io.ne.0) call error_reading_header(0)
  
  write(outname,'(I5.5,A4)')jmod,'.mod'
  inquire(file=trim(outname), exist=ex)
  if(ex) then
     write(6,'(A)')'  '//trim(outname)//' exists.'
     write(6,'(A,$)')'  Please enter a different name for the output file: '
     read*,outname
  end if
     
  open (unit=20,form='formatted',status='new',file=trim(outname),iostat=io)
  if(io.ne.0) call quit_program('Error opening output file '//trim(outname)//'.')

  
  !Write header line
  write(20,'(1X, 8ES23.15, 6I6)')dat1,kh,kp,jmod,jb,jin,0
  
  !Copy model block:
  do li=1,kh  !Line/mesh point
     read(10,*,iostat=io)dat2
     if(io.ne.0) call error_reading_block(li)
     write(20,'(1X, 24ES23.15)')dat2
  end do !li
  close(10)
  close(20)
  write(6,'(A)')'  Model written to '//trim(outname)//'.'
  
  
end subroutine copy_mod
!***********************************************************************************************************************************
