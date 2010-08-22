! getgrid.f, Find the parameters you need to get a grid with n values between x1 and x2

! Copyright 2002-2010 AstroFloyd - astrofloyd.org
! 
! 
! This file is part of the evTools package.
! 
! This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with this code.
! If not, see <http://www.gnu.org/licenses/>.

program getgrid      
  use kinds
  
  implicit none
  real(double) :: x1,x2,dlgx
  integer :: n,command_argument_count
  character :: bla*99
  
  if(command_argument_count().ne.3) then
     write(*,'(/,A)')'  This program returns the parameters you need for a certain grid of N models with values between'// &
          ' X1 and X2 in ev'
     write(6,'(A,/)')'  Syntax:  getgrid <X1> <X2> <N>'
     stop
  end if
  
  call get_command_argument(1,bla)
  read(bla,*)x1
  call get_command_argument(2,bla)
  read(bla,*)x2
  call get_command_argument(3,bla)
  read(bla,*)n
  
  dlgx = dlog10(x2/x1)/real(n-1)
  
  write(6,'(/,A,2ES11.3,I5)')'  Start at first value:  ',dlog10(x1),dlgx,n
  write(6,'(A,2ES11.3,I5,/)')'  Start at second value: ',dlog10(x1)+dlgx,dlgx,n-1
  
  write(*,*)
  call printgrid(dlog10(x1),dlgx,n)
  write(*,*)
  
  write(6,'(/,A,2ES11.3,I5)')'  Intermediate grid:     ',dlog10(x1)+dlgx/2.,dlgx,n-1
  write(*,*)
  call printgrid(dlog10(x1)+dlgx/2.,dlgx,n-1)
  write(*,*)
  
end program getgrid



!Shared with grid.f
subroutine printgrid(xi1,dx,n)
  use kinds
  
  implicit none
  real(double) :: xi1,xi,dx
  integer :: i,n
  
  xi = xi1
  write(6,'(A5,2A10)')'i','x','log x'
  do i=1,n
     write(6,'(I5,F10.5,ES12.3)')i,10.d0**xi,xi
     xi = xi + dx
  end do
  
end subroutine printgrid
