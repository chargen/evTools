!> \file plotmdln.f90   Plots the data contained in star.mdl[12], for some or all (n) structure models in the file, 
!!                      and for ONE set of variables


! AF, 19-05-2005


! Copyright 2002-2011 AstroFloyd - astrofloyd.org
! 
! 
! This file is part of the evTools package.
! 
! This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
! by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with this code.  If not, see 
! <http://www.gnu.org/licenses/>.


!***********************************************************************************************************************************
!> \brief Plot the data contained in star.mdl[12], for some or all (n) structure models in the file, and for ONE set of variables

program plotmdln  
  use constants
  use mdl_data !, mdl2 => mdl
  
  implicit none
  integer, parameter :: nbb=30  !,nn=1000,nq=30
  integer :: mdl(nbb),mdl1, nblk,nsel  !,nv, nm
  real :: dat(nbb,nq,nn), dat1(nq,nn), age(nbb),age1, x,xx(nbb,nn),yy(nbb,nn) !,dov
  real :: xmin,xmax,ymin,ymax,xmin0,xmax0,ymin0,ymax0
  real :: xsel(4),ysel(4)
  !real :: m1,r1,l1,ts,tc,mhe,mco
  !real :: hc,hec,cc,oc,zc
  !real :: hs,hes,cs,os,zs
  !real :: rhoc
  
  integer i,ii,j,nr,blk(nbb),vx,vy,plotagain,bl,nb,plt,nbmax !,system,status
  character findfile*(99),fname*(99),rng,log,mdlnr*(5)
  character :: lx*(50),ly*(50),title*(100)  !,labels(nq)*(60)
  
  
  !Set constants:
  call setconstants()
  write(6,*)
  call print_code_version(6)  !To screen
  
  call evTools_settings()
  
  
  plotagain = 0
  
  !Axis labels
  !labels(1)='\u\(2263) centre\d    Mesh point    \usurface \(2261)'
  !labels(2) = 'M (M\d\(2281)\u)'
  !labels(3) = 'R (R\d\(2281)\u)'
  !labels(4) = 'P (dyn cm\u-2\d)'
  !labels(5) = '\gr (g cm\u-3\d)'
  !labels(6) = 'T (K)'
  !labels(7) = 'k (cm\u2\d g\u-1\d)'
  !labels(8) = '\(2266)\dad\u'
  !labels(9) = '\(2266)\drad\u'
  !labels(10) = 'H abundance'
  !labels(11) = 'He abundance'
  !labels(12) = 'C abundance'
  !labels(13) = 'N abundance'
  !labels(14) = 'O abundance'
  !labels(15) = 'Ne abundance'
  !labels(16) = 'Mg abundance'
  !labels(17) = 'L (L\d\(2281)\u)'
  !labels(18) = '\ge\dth\u'
  !labels(19) = '\ge\dnucl\u'
  !labels(20) = '\ge\d\gn\u'
  !labels(21) = 'S (cgs)'
  !labels(22) = 'U\dint\u (erg g\u-1\d)'
  !labels(23) = '\(2266)\drad\u - \(2266)\dad\u'
  !labels(24) = 'm/M\d*\u'
  !labels(25) = 'r/R\d*\u'
  !labels(26) = 'C/O'
  !labels(27) = 'Ne/O'
  !labels(28) = '(Ne/O)/(Ne/O)\d0\u'
  
  ! Define variable labels:
  call set_mdl_labels()
  
  
3 continue
  
  ! Read currend path and use it as plot title:
  !status = system('pwd > tmppwd.txt')
  !open (unit=10,form='formatted',status='old',file='tmppwd.txt')
  !rewind 10
  !read(10,'(a100)')title
  !close(10)
  !status = system('rm tmppwd.txt')
  
  ! Read current path and use it as plot title:
  write(title,'(A)')trim(workdir)
  
  ! Get filename:
  if(command_argument_count().eq.1) then
     call get_command_argument(1,fname)
  else
     fname = findfile('*.mdl*')  ! Search for input file in current dir
  end if
  
  
  
  !************************************************************************      
  !***   READ ALL STRUCTURE MODELS IN THE FILE AND DISPLAY MAIN PROPERTIES
  !************************************************************************      
  
  write(6,*)''
4 continue
  call list_mdl_models(trim(fname),nblk)
  
!  write(6,'(A)')'Reading file '//trim(fname)
!  open (unit=10,form='formatted',status='old',file=trim(fname))
!  rewind 10
!  
!  read(10,5,err=11,end=11) nm,nv,dov
!5 format (2x,I4,4x,I2,F7.3)
!  write(6,*)''
!  write(6,'(A)')' Nr  Model Nmsh          Age        M1   Mhe   Mco     Menv         R        L     Teff       Tc     Rhoc'//  &
!       '      Xc     Yc     Cc     Oc     Xs    Ys    Zs'
!  
!  do b=1,nbb
!     read(10,6,err=12,end=15) mdl(b),age(b)
!6    format (I6,1x,E16.9)
!     mhe = 0.
!     mco = 0.
!     do j=1,nm
!        read(10,7,err=13,end=15) (dat(b,i+1,j),i=1,21)
!        if(j.eq.1) then
!           tc   = dat(b,6,j)
!           hc   = dat(b,10,j)
!           hec  = dat(b,11,j)
!           cc   = dat(b,12,j)
!           oc   = dat(b,14,j)
!           zc   = 1. - hc - hec
!           rhoc = dat(b,5,j)
!        end if
!        if(j.eq.nm) then
!           m1  = dat(b,2,j)
!           r1  = dat(b,3,j)
!           l1  = dat(b,17,j)
!           ts  = dat(b,6,j)
!           hs  = dat(b,10,j)
!           hes = dat(b,11,j)
!           cs  = dat(b,12,j)
!           os  = dat(b,14,j)
!           zs  = 1. - hs - hes
!        end if
!        if(mhe.eq.0.0.and.dat(b,10,j).gt.0.1) mhe = dat(b,2,j)
!        if(mco.eq.0.0.and.dat(b,11,j).gt.0.1) mco = dat(b,2,j)
!     end do !j
!7    format (1P,E13.6,4E11.4,16E11.3)
!     
!     if(mod(b,25).eq.0) then
!        write(6,*)''
!        write(6,'(A)')' Nr  Model Nmsh          Age        M1   Mhe   Mco     Menv         R        L     Teff       Tc'//  &
!             '     Rhoc      Xc     Yc     Cc     Oc     Xs    Ys    Zs'
!     end if
!     write(6,9)b,mdl(b),nm,age(b),m1,mhe,mco,m1-mhe,r1,l1,ts,tc,rhoc,hc,hec,cc,oc,hs,hes,zs!,bms,p,p1
!  end do !b
!  
!9 format(I4,I7,I5,ES13.5,f10.4,2f6.3,ES9.2,1x,4ES9.2,ES9.2,1x,4f7.4,1x,3f6.3)
!  
!  write(6,'(A)')'  Arrays are too small !'
!  goto 9999
!  
!11 write(6,'(A)')'  Error reading first line of file, aborting...'
!  close(10)
!  goto 9999
!12 write(6,'(A)')'  Error reading first line of block, aborting...'
!  close(10)
!  goto 15
!  goto 9999
!13 write(6,'(A,I4,A,I5,A)')'  Error reading block',i-1,' line',j-1,', aborting...'
!  close(10)
!  goto 9999
!15 close(10)
!  
!  write(6,'(A)')' Nr  Model Nmsh          Age        M1   Mhe   Mco     Menv         R        L     Teff       Tc     Rhoc'//  &
!       '      Xc     Yc     Cc     Oc     Xs    Ys    Zs'
!  write(6,*)''
!  nblk = b-1
!  write(6,'(A,I4,A)')'  EOF reached,',nblk,' blocks read.'
!  write(6,*)''
! 
  
  if(nblk.eq.0) goto 9999
  
  if(nblk.eq.1) then
     blk(1) = 1
     nb = 1
     goto 24
  end if
  
  
  !************************************************************************      
  !***   CHOOSE STRUCTURE MODELS
  !************************************************************************      
  
  write(6,'(A)')'Which structure models do you want to plot:' 
  write(6,'(A,I2,A3)')'  (press ENTER after each number, 0 for all models and -1 to end the list) (1-',nblk,'): '
  
  
  bl = 1
  nbmax = 30  !nblk   !Read max 30 models, so that all labels can be printed next to plot
  do j=1,nbmax
     read*,blk(bl)
     if(blk(bl).eq.-1) then
        nb = bl-1
        goto 24
     end if
     if(blk(bl).eq.0) then
        do i=1,nbmax
           blk(i) = i
           nb = min(nbmax,nblk)
        end do
        goto 24
     end if
     if(blk(bl).lt.0.or.blk(bl).gt.nblk) bl = bl-1
     bl = bl+1
  end do
  nb = bl-1
  
  
  
  ! Read the models again:
  call read_first_mdls(fname,0)  ! Read header
  do bl = 1,nb
     call read_chosen_mdl(bl, mdl1,age1,dat1)
     mdl(bl) = mdl1
     age(bl) = age1
     dat(bl,:,:) = dat1(:,:)
  end do
  
  
  
  !************************************************************************      
  !***   PRINT A LIST OF CHOSEN STRUCTURE MODELS
  !************************************************************************      
  
24 continue
  goto 30
  
  write(6,*)''
  write(6,'(I5,A)')nb,' blocks selected:'
  write(6,'(A)')'  Nr  Model       Age    Mass  Radius   Luminos      Teff        H     He    C+O'
  do bl=1,nb
     write(6,8)blk(bl),mdl(blk(bl)),age(blk(bl)),dat(blk(bl),2,nm),    dat(blk(bl),3,nm),dat(blk(bl),17,nm),dat(blk(bl),6,nm),  &
          dat(blk(bl),10,1),dat(blk(bl),11,1),    dat(blk(bl),12,1)+dat(blk(bl),14,1)
  end do
  write(6,*)''
8 format(2(I5,2x),1PE8.2,2x,2(0PF6.2,2x),2(1PE8.2,2x),3(0PF5.2,2x))
  
  
  
  
  
  
  ! Add file name to plot title:
  write(title,'(A)')trim(title)//'/'//trim(fname)
  
  
  !************************************************************************      
  !***   CREATE EXTRA PLOT VARIABLES
  !************************************************************************      
  
30 continue
  
  if(plotagain.eq.0) then
     !do i=1,nm
     !   dat(:,1,i) = real(i)
     !end do
     !
     !dat(:,23,1:nm) = dat(:,9,1:nm) - dat(:,8,1:nm)
     !dat(:,23,1:nm) = dat(:,23,1:nm)/abs(dat(:,23,1:nm))
     !do bl=1,nb
     !   dat(blk(bl),24,1:nm) = dat(blk(bl),2,1:nm)/dat(blk(bl),2,nm)
     !   dat(blk(bl),25,1:nm) = dat(blk(bl),3,1:nm)/dat(blk(bl),3,nm)
     !   dat(blk(bl),26,1:nm) = dat(blk(bl),12,1:nm)/dat(blk(bl),14,1:nm)
     !   dat(blk(bl),27,1:nm) = dat(blk(bl),15,1:nm)/dat(blk(bl),14,1:nm)
     !   dat(blk(bl),28,1:nm) = dat(blk(bl),27,1:nm)/dat(blk(bl),27,nm)
     !end do
     
     !***   COMPUTE ADDITIONAL PLOT VARIABLES
     do bl=1,nb
        dat1 = dat(bl,:,:)
        if(plotagain.eq.0) call compute_mdl_variables(dat1)
        dat(bl,:,:) = dat1
     end do
  end if !if(plotagain.eq.0) then
  
  
  
  !************************************************************************      
  !***   CHOOSE PLOT VARIABLES
  !************************************************************************      
  
32 continue
  
  !write(6,*)''
  !write(6,'(A)')'Variables:                       0: Quit'
  !write(6,*)''
  !write(6,'(A)')' 1: Mesh pt'
  !write(6,'(A)')' 2: M                                                    '
  !write(6,'(A)')' 3: R       10: H       17: L           24: M/M*       '
  !write(6,'(A)')' 4: P       11: He      18: Eth         25: R/R*       '
  !write(6,'(A)')' 5: Rho     12: C       19: Enuc        26: C/O       '
  !write(6,'(A)')' 6: T       13: N       20: Enu         27: Ne/O              '
  !write(6,'(A)')' 7: k       14: O       21: S           28: Ne/O change  '
  !write(6,'(A)')' 8: Nad     15: Ne      22: Uint                         '
  !write(6,'(A)')' 9: Nrad    16: Mg      23: Nrad-Nad                     '
  !write(6,*)''
  
  
  nr = 4 ! Number of variable columns
  ii = ceiling(real(nc)/real(nr)) !Number of rows
  write(6,'(A)')' Variables:                         0: Quit                   ' 
  do i=1,ii
     do j=0,nr-1
        if(pxnr(i+j*ii).eq.0) then
           write(6,'(A19)',advance='no')''
        else
           write(6,'(I9,A10,5x)',advance='no')i+j*ii,': '//pxns(pxnr(i+j*ii))
        end if
     end do
     write(6,*)
  end do
  
  
  ! Print derived variables, from number 201 on:
  write(6,'(A)')'                                                              '
  write(6,'(A)')'  Derived variables:                                          '
  
  nr = 4 !Number of variable columns
  ii = ceiling(real(nv_der)/real(nr)) !Number of rows
  do i=1,ii
     do j=0,nr-1
        if(pxnr(200+i+j*ii).eq.0) then
           write(6,'(A19)',advance='no')''
        else
           write(6,'(I9,A10,5x)',advance='no')200+i+j*ii,': '//pxns(pxnr(200+i+j*ii))
        end if
     end do
     write(6,*)
  end do
  
  
  ! Print special variables, from number 301 on:
  write(6,'(A)')'                                                              '
  write(6,'(A)')'  Special plots:                                          '
  
  nr = 2 !Number of variable columns
  ii = ceiling(real(nv_sp)/real(nr)) !Number of rows
  do i=1,ii
     do j=0,nr-1
        if(pxnr(300+i+j*ii).eq.0) then
           write(6,'(A39)',advance='no')''
        else
           write(6,'(I9,A30,5x)',advance='no')300+i+j*ii,': '//pxns(pxnr(300+i+j*ii))
        end if
     end do
     write(6,*)
  end do
  
  write(6,*)''
  write(6,*)''
  
  
  
!35 continue
  write(6,'(A36)',advance='no')' Choose the X-axis variable (1-28): '
  !read*,vx
  !if(vx.eq.0) goto 9999
  !if(vx.lt.1.or.vx.gt.28) goto 35
  vx = nq
  do while (pxnr(vx).le.0)
     write(6,'(A)',advance='no')' Choose the X-axis variable: '
     read*,vx
     vx = min(max(vx,0),nq)
     if(vx.eq.0) then
        write(6,'(A,/)')' Program finished'
        stop
     end if
  end do
  
!36 continue
  write(6,'(A36)',advance='no')' Choose the Y-axis variable (1-28): '
  !read*,vy
  !if(vy.eq.0) goto 9999
  !if(vy.lt.1.or.vy.gt.28) goto 36
  vy = nq
  do while (pxnr(vy).le.0)
     write(6,'(A)',advance='no')' Choose the Y-axis variable: '
     read*,vy
     vy = min(max(vy,0),nq)
     if(vy.eq.0) then
        write(6,'(A,/)')' Program finished'
        stop
     end if
  end do
  
  
41 continue
  lx = trim(labels(vx))
  ly = trim(labels(vy))
  
  xx(1:nb,1:nm) = dat(blk(1:nb),vx,1:nm)
  yy(1:nb,1:nm) = dat(blk(1:nb),vy,1:nm)
  
  
  
  
  
  !************************************************************************      
  !***   LIN/LOG AXES
  !************************************************************************      
  
  write(6,'(A68)',advance='no')' Do you want a logarithmic scale: (N)o, (X)-axis, (Y)-axis, (bl)oth: '
  read*,log
  if(log.eq.'X') log='x'
  if(log.eq.'Y') log='y'
  if(log.eq.'B') log='b'
  if(log.eq.'N') log='n'
  
  if(log.eq.'x'.or.log.eq.'b') then
     do bl=1,nb
        if(xx(bl,1).eq.0.) xx(bl,1) = xx(bl,2)
     end do
     xx(1:nb,1:nm) = log10(abs(xx(1:nb,1:nm))+1.e-30)
     lx = trim('log '//lx)
  end if
  if(log.eq.'y'.or.log.eq.'b') then
     do bl=1,nb
        if(yy(bl,1).eq.0.) yy(bl,1) = yy(bl,2)
     end do
     yy(1:nb,1:nm) = log10(abs(yy(1:nb,1:nm))+1.e-30)
     ly = trim('log '//ly)
  end if
  
  xmin = minval(xx(1:nb,1:nm))
  xmax = maxval(xx(1:nb,1:nm))
  ymin = minval(yy(1:nb,1:nm))
  ymax = maxval(yy(1:nb,1:nm))
  
  
  
  
  
  
  !************************************************************************      
  !***   PLOT RANGE
  !************************************************************************      
  
  xmin0 = xmin
  xmax0 = xmax
  ymin0 = ymin
  ymax0 = ymax
  
70 write(6,*)''
  write(6,'(A,ES12.3,A1,ES12.3)')'  X-range:',xmin,'-',xmax
  write(6,'(A,ES12.3,A1,ES12.3)')'  Y-range:',ymin,'-',ymax
  write(6,'(A)',advance='no')'  Do you want to change a plot range ?  (N)o, (X)-axis, (Y)-axis, (bl)oth:  '
  read*,rng
  if(rng.eq.'N') rng='n'
  if(rng.eq.'X') rng='x'
  if(rng.eq.'Y') rng='y'
  if(rng.eq.'B') rng='b'
  
  if(rng.eq.'n'.or.rng.eq.' ') goto 100
  
  
  if(rng.eq.'x'.or.rng.eq.'b') then
     write(6,'(A)')'Give the new range for the X-axis (Xmin, Xmax):'
     read*,xmin,xmax
     if(xmin.gt.xmax) then
        x = xmin
        xmin = xmax
        xmax = x
        write(6,'(A)')'  Swapped Xmin and Xmax'
     end if !if(xmin.gt.xmax)
     if(xmin.lt.xmin0) xmin = xmin0
     if(xmax.gt.xmax0) xmax = xmax0
  end if
  
  if(rng.eq.'y'.or.rng.eq.'b') then
     write(6,'(A)')'Give the new range for the Y-axis (Ymin, Ymax):'
     read*,ymin,ymax
     if(ymin.gt.ymax) then
        x = ymin
        ymin = ymax
        ymax = x
        write(6,'(A)')'  Swapped Ymin and Ymax'
     end if !if(ymin.gt.ymax)
     if(ymin.lt.ymin0) ymin = ymin0
     if(ymax.gt.ymax0) ymax = ymax0
  end if
  
  write(6,*)''
  write(6,'(A,ES12.3,A1,ES12.3)')'X-range:',xmin,'-',xmax
  write(6,'(A,ES12.3,A1,ES12.3)')'Y-range:',ymin,'-',ymax
  
  
100 continue
  x = 0.02*abs(xmax-xmin)
  if(x.eq.0.) x = 0.05*xmax
  xmin = xmin - x
  xmax = xmax + x
  x = 0.02*abs(ymax-ymin)
  if(x.eq.0.) x = 0.05*ymax
  ymin = ymin - x
  ymax = ymax + x
  
  
  
  
  
  
  
  
  
  
  
  
  !************************************************************************      
  !***   PLOT TO SCREEN, FILE
  !************************************************************************      
  
501 continue
  do plt=1,2    !Plot to PS, then screen
     
     if(plt.eq.1) call pgbegin(1,'plot_mdln.eps/cps',1,1)
     if(plt.eq.2) then
        call pgbegin(1,'3/xserve',1,1)
        call pgpap(scrsz,scrrat) !MacBook, Gentoo
        if(white_bg) then     !Create a white background; swap black (ci=0) and white (ci=1)
           call pgscr(0,1.,1.,1.)  !For some reason, this needs to be repeated for AquaTerm, see below
           call pgscr(1,0.,0.,0.)
           call pgsci(1)
           call pgsci(0)
           call pgsvp(0.,1.,0.,1.)
           call pgswin(-1.,1.,-1.,1.)
           call pgrect(-2.,2.,-2.,2.)
           call pgsci(1)
        end if
     end if
     call pgscf(1)
     call pgscr(7,0.7,0.7,0.7)  !Replace yellow by light grey
     
     call pgsvp(0.06,0.95,0.07,0.96)
     call pgswin(xmin,xmax,ymin,ymax)
     
     if(log.eq.'x') then
        call pgbox('BCLNTS',0.0,0,'BCNTS',0.0,0)  !log x, lin y
     else if(log.eq.'y') then
        call pgbox('BCNTS',0.0,0,'BCLNTS',0.0,0)  !lin x, log y
     else if(log.eq.'b') then
        call pgbox('BCLNTS',0.0,0,'BCLNTS',0.0,0) !log x,y
     else
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)   !lin x,y
     end if
     
     call pgmtxt('T',0.7,0.5,0.5,title(12:))  !13, to remove /home/user/
     call pgmtxt('B',2.4,0.5,0.5,lx)
     call pgmtxt('L',2.0,0.5,0.5,ly)
     
     if(plt.eq.1) call pgslw(2)
     
     call pgsch(0.7)
     if(vx.ne.1.and.vy.ne.1) then
        do bl=1,nb
           call pgsls(mod((bl-1)/12,5)+1)
           call pgsci(mod(bl,12)+1)
           call pgline(nm,xx(bl,1:nm),yy(bl,1:nm))
           write(mdlnr,'(I5)') mdl(blk(bl))
           call pgmtext('RV',0.5,1.-real(bl-1)/30.,0.,mdlnr)
        end do
     end if
     
     if(vx.eq.1.or.vy.eq.1) then
        do bl=1,nb
           call pgsci(mod(bl-1,12)+2)
           call pgpoint(nm,xx(bl,1:nm),yy(bl,1:nm),1)
        end do
     end if
     
     call pgsci(1)
     call pgsch(1.)
     call pgsls(2)
     if(vy.eq.22) call pgline(2,(/xmin,xmax/),(/0.,0./))
     
     if(plt.eq.1) call pgend
     
  end do  !do plt=1,2    !Plot to screen, then file
  
  
  
  
  
  !************************************************************************      
  !***   FINISH
  !************************************************************************      
  
900 write(6,*)''
  write(6,'(A)')' You can:'
  write(6,'(A)')'  0) quit'
  write(6,'(A)')'  1) change variables'
  write(6,'(A)')'  2) change lin/log axes'
  write(6,'(A)')'  3) change axis ranges'
  write(6,'(A)')'  4) select zoom region'
  write(6,'(A)')'  5) zoom out'
  write(6,'(A)')'  6) change structure model'
  write(6,'(A)')'  7) change input file'
  write(6,*)''
  write(6,'(A27)',advance='no')' What do you want to do ?  '
  read*,plotagain
  if(plotagain.lt.0.or.plotagain.gt.7) goto 900
  
  if(plotagain.ne.4) call pgend
  if(plotagain.eq.1) goto 32
  if(plotagain.eq.2) goto 41
  if(plotagain.eq.3) goto 70
  if(plotagain.eq.6) goto 4
  if(plotagain.eq.7) goto 3
  
  if(plotagain.eq.4) then  !Select region
941  call pgsci(1)
     xsel = 0.
     ysel = 0.
     write(6,'(A)')'Select 2-4 corner points with your left mouse button and press "x" to finish'
     nsel=0
     call pgolin(4,nsel,xsel,ysel,2)
     if(nsel.lt.2) then
        write(6,'(A)')' I need at least 2 corner points...'
        goto 941
     end if
     xmin = minval(xsel(1:nsel))  !The new window is drawn for the extreme values of these points
     xmax = maxval(xsel(1:nsel))
     ymin = minval(ysel(1:nsel))
     ymax = maxval(ysel(1:nsel))
     write(6,*)''
     write(6,'(A,ES12.3,A1,ES12.3)')'X-range:',xmin,'-',xmax
     write(6,'(A,ES12.3,A1,ES12.3)')'Y-range:',ymin,'-',ymax
     write(6,*)''
     call pgend
     goto 501
  end if
  
  if(plotagain.eq.5) then  !Zoom out
     xmin = (xmin+xmax)/2. - 2*abs((xmin+xmax)/2.-xmin) !Central value - 2x the 'radius', 'radius' = central value - minimum
     xmax = (xmin+xmax)/2. + 2*abs((xmin+xmax)/2.-xmin)
     ymin = (ymin+ymax)/2. - 2*abs((ymin+ymax)/2.-ymin)
     ymax = (ymin+ymax)/2. + 2*abs((ymin+ymax)/2.-ymin)
     write(6,*)''
     write(6,'(A,ES12.3,A1,ES12.3)')'X-range:',xmin,'-',xmax
     write(6,'(A,ES12.3,A1,ES12.3)')'Y-range:',ymin,'-',ymax
     write(6,*)''
     goto 501
  end if
  
9999 continue
  write(6,'(A,/)')'Program finished'
  
end program plotmdln
!***********************************************************************************************************************************




