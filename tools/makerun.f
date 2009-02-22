! Reads, optinally changes and (over!)writes an init.run (fort.23) input file
! AF January 21, 2004

program makerun  
  implicit none
  real*8 :: ct1(7),ct2(7),ct3(7)
  real*8 :: ml1,dml,ql1,dql,xl1,dxl
  real*8 :: rot,ex
  real*8 :: sm,dty,age,per,bms,ecc,p1,enc
  real*8 :: m2
  integer :: isb,ktw,ip1,im1,ip2,im2,kpt,kp
  integer :: kml,kql,kxl,kr,jmx
  integer :: io,narg,iargc
  character :: file*8,arg*10,bla*500
  
  write(6,*)''
  file = 'init.run'
  
  
  open(unit=10,form='formatted',status='old',file=file,iostat=io)
  if(io.ne.0) goto 90
  rewind 10
  read(10,*,err=91) isb,ktw,ip1,im1,ip2,im2,kpt,kp
  read(10,*,err=92) ml1,dml,kml
  read(10,*,err=93) ql1,dql,kql
  read(10,*,err=94) xl1,dxl,kxl
  read(10,*,err=95) rot,kr,ex
  read(10,*,err=96) sm,dty,age,per,bms,ecc,p1,enc,jmx
  read(10,*,err=97) ct1
  read(10,*,err=98) ct2
  read(10,*,err=99) ct3
  !close(10)
  
  
  kml = 1  !Only one iteration in mass
  m2 = 0.5d0*sm
  if(bms.gt.0) m2 = bms-sm
  
  narg = iargc()
  if(narg.eq.1) then
     call getarg(1,arg)
     read(arg,*)sm
  else if(narg.eq.2) then
     call getarg(1,arg)
     read(arg,*)sm
     call getarg(2,arg)
     read(arg,*)per
  else if(narg.eq.3) then
     call getarg(1,arg)
     read(arg,*)sm
     call getarg(2,arg)
     read(arg,*)m2
     call getarg(3,arg)
     read(arg,*)per
     write(6,'(A)')'  Synchronising binary...'
     p1 = per
  else if(narg.eq.4) then
     call getarg(1,arg)
     read(arg,*)sm
     call getarg(2,arg)
     read(arg,*)m2
     call getarg(3,arg)
     read(arg,*)per
     call getarg(4,arg)
     read(arg,*)p1
  else
     write(6,'(A)')'  Syntax: '
     write(6,'(A)')'    makerun <M1>'
     write(6,'(A)')'    makerun <M1> <Porb>'
     write(6,'(A)')'    makerun <M1> <M2> <Porb> (synchonise: Prot=Porb)'
     write(6,'(A,/)')'    makerun <M1> <M2> <Porb> <Prot>'
     stop
  end if
  
  if(bms.gt.0.d0.or.narg.eq.3) bms = sm + m2
  !dty = 1.d5
  !p1 = per  !Syncronise the initial binary
  !print*,sm,m2,bms
  
  ml1 = log10(sm)
  if(narg.eq.3) ql1 = log10(sm/m2)
  if(per.gt.0.d0) xl1 = log10(per)
  
  open(unit=20,form='formatted',file=file)
  write(20,50) isb,ktw,ip1,im1,ip2,im2,kpt,kp,  &
       ml1,dml,kml,ql1,dql,kql,xl1,dxl,kxl,  &
       rot,kr,ex,  &
       sm,dty,age,per,bms,ecc,p1,enc,jmx,  &
       ct1,ct2,ct3
  read(10,*)bla
  write(20,'(/,A)')'last five lines:'
  
  io = 0
  do while(io.eq.0)
     read(10,'(A500)',iostat=io)bla
     if(io.ne.0.or.len_trim(bla).eq.0.or.len_trim(bla).eq.500) exit
     write(20,'(A)')trim(bla)
  end do
  
  close(10)
  close(20)
  
50 format (6I6,1x,2I7,/,  3(2ES11.3,I5,/),  ES11.3,I3,ES10.2,/,   ES11.3,ES12.4,6ES10.2,I6,/,      3(7ES10.2,/))
  
  
  write(6,'(4(A,ES10.3))')'  M1:',sm,',  M2:',m2,',  Porb:',per,',  Prot,1:',p1
  write(6,'(A,/)')'  Program done'
  stop
  
90 write(6,'(A,/)')'  Error opening file: '//trim(file)
  stop
91 write(6,'(A,/)')'  Error reading file: '//trim(file)//', line 1'
  stop
92 write(6,'(A,/)')'  Error reading file: '//trim(file)//', line 2'
  stop
93 write(6,'(A,/)')'  Error reading file: '//trim(file)//', line 3'
  stop
94 write(6,'(A,/)')'  Error reading file: '//trim(file)//', line 4'
  stop
95 write(6,'(A,/)')'  Error reading file: '//trim(file)//', line 5'
  stop
96 write(6,'(A,/)')'  Error reading file: '//trim(file)//', line 6'
  stop
97 write(6,'(A,/)')'  Error reading file: '//trim(file)//', line 7'
  stop
98 write(6,'(A,/)')'  Error reading file: '//trim(file)//', line 8'
  stop
99 write(6,'(A,/)')'  Error reading file: '//trim(file)//', line 9'
  stop
  
end program makerun
