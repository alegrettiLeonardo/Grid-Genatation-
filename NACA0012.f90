
program naca0012
    implicit none
    integer,parameter :: cle = 360 ! x - direction
    integer,parameter :: row = 360 ! y - direction
    real :: x(0:cle,0:row),xx(0:cle,0:row)
    real :: y(0:cle,0:row),yy(0:cle,0:row)
    real :: mx(0:cle,0:row),my(0:cle,0:row)
    real :: f(0:cle,0:row),ff(0:cle,0:row),mf(0:cle,0:row),jj(0:cle,0:row),u(0:cle,0:row),v(0:cle,0:row),cp(0:cle,0:row)
    real :: aa(0:cle,0:row),bb(0:cle,0:row),rr(0:cle,0:row)
    real :: g = 0.0001, e = 0.001
    real :: pi = 3.1415926
    real :: n,A,B,R,H,H1,H2,mx1,my1,p,mf1,s
    integer :: c = 1, radius = 50
    integer :: m,i,j,k
    integer :: z = 0,w = 0
!--------------------------------------------------------------------------------------------------------------------------
! NACA 0012 Equations
    do i=0,cle/2
        n = i*2*pi/cle
        x(i,0) = 0.5*c*(1 + cos(n))
        y(i,0) = 0.6*(-0.1015*x(i,0)**4 + 0.2843*x(i,0)**3 - 0.3576*x(i,0)**2-0.1221*x(i,0) + 0.2969*x(i,0)**0.5)
        x(cle-i,0) = x(i,0)
        y(cle-i,0) = -y(i,0)

        x(i,row) = radius*c*cos(n)
        y(i,row) = radius*c*sin(n)
        x(cle-i,row) = x(i,row)
        y(cle-i,row) = -y(i,row)
    end do
!       x(cle,0) = x(0,0)
!	y(cle,0) = y(0,0)
!	x(cle,row) = x(0,row)
!	y(cle,row) = y(0,row)

  do i=0,cle
   do j=1,row-1
	  x(i,j) = x(i,0) + (x(i,row) - x(i,0))/row*j
	  y(i,j) = y(i,0) + (y(i,row) - y(i,0))/row*j
	end do
  end do

!-------------------------------------------------------------------------------------------------------------------------
  open(10,file='initial condition.plt')
  WRITE(10,*) 'VARIABLES= "Y", "Z"'
  WRITE(10,*) 'ZONE I=',cle+1,'J=',row+1
  do j=0,row
    do i=0,cle
	write(10,*) x(i,j),y(i,j)
    end do
  end do

  open(10,file='Chord_1.dat')
  !WRITE(10,*) cle+1
  do j=0,row
    do i=0,cle
	write(10,*) x(i,0)
    end do
  end do

  close(10)
  print *,"Finish Initial Condition!"
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
! SOR
    mx1 = 0.001
    my1 = 0.001
    do while(mx1>g.and.my1>g)
       do i=0,cle
         do j=1,row-1
            if(i==0)then
                A = (x(i,j+1) - x(i,j-1))**2 + (y(i,j+1) - y(i,j-1))**2
                B = (x(i+1,j) - x(cle-1,j))*(x(i,j+1) - x(i,j-1)) + (y(i+1,j)-y(cle-1,j))*(y(i,j+1)-y(i,j-1))
                R = (x(i+1,j) - x(cle-1,j))**2 + (y(i+1,j) - y(cle-1,j))**2
                xx(i,j) = 0.5*(A*(x(i+1,j) + x(cle-1,j)) + R*(x(i,j+1) + x(i,j-1)) - 0.5*B*(x(i+1,j+1) + x(cle-1,j-1) - x(i+1,j-1) - x(cle-1,j+1)))/(A + R)
                yy(i,j) = 0.5*(A*(y(i+1,j) + y(cle-1,j)) + R*(y(i,j+1) + y(i,j-1)) - 0.5*B*(y(i+1,j+1) + y(cle-1,j-1) - y(i+1,j-1) - y(cle-1,j+1)))/(A + R)
                mx(i,j) = abs(xx(i,j)-x(i,j))
                my(i,j) = abs(yy(i,j)-y(i,j))

            else if(i>0.and.i<cle)then
                A = (x(i,j+1) - x(i,j-1))**2 + (y(i,j+1)-y(i,j-1))**2
                B = (x(i+1,j) - x(i-1,j))*(x(i,j+1) - x(i,j-1)) + (y(i+1,j) - y(i-1,j))*(y(i,j+1) - y(i,j-1))
                R = (x(i+1,j) - x(i-1,j))**2 + (y(i+1,j)-y(i-1,j))**2
                xx(i,j) = 0.5*(A*(x(i+1,j) + x(i-1,j)) + R*(x(i,j+1) + x(i,j-1)) - 0.5*B*(x(i+1,j+1) + x(i-1,j-1) - x(i+1,j-1) - x(i-1,j+1)))/(A + R)
                yy(i,j) = 0.5*(A*(y(i+1,j) + y(i-1,j)) + R*(y(i,j+1) + y(i,j-1)) - 0.5*B*(y(i+1,j+1) + y(i-1,j-1) - y(i+1,j-1) - y(i-1,j+1)))/(A + R)
                mx(i,j) = abs(xx(i,j) - x(i,j))
                my(i,j) = abs(yy(i,j) - y(i,j))

            else
                xx(i,j) = xx(0,j)
                yy(i,j) = yy(0,j)
            end if
         end do
       end do

      mx1 = mx(1,1)
      my1 = my(1,1)
        do i=0,cle
          do j=1,row-1
          x(i,j) = xx(i,j)
          y(i,j) = yy(i,j)
          if(mx(i,j)>mx1) then
            mx1 = mx(i,j)
          end if
          if(my(i,j)>my1) then
            my1 = my(i,j)
            end if
          end do
        end do
        z = z + 1
        !print *,z,mx1,my1,g
    end do
    !-----------------------------------------------------------------
    open(10,file='grid.plt')
    WRITE(10,*) 'VARIABLES= "Y", "Z"'
    WRITE(10,*) 'ZONE I=',cle+1,'J=',row+1
    do j=0,row
        do i=0,cle
            write(10,*) x(i,j),y(i,j)
        end do
    end do
    close(10)

    open(10,file='NACA0012.tec')
    write(10,*) ' title="naca0012_grid"'
    write(10,*) ' variables= x, y'
    write(10,*) ' zone t="grid", i=',cle+1,'j=',row+1,', datapacking=point'
    do j=0,row
        do i=0,cle
            write(10,*) x(i,j), y(i,j)
        enddo
    enddo

    open(10,file='NACA0012.dat')
    !write(1,*) ' title="naca0012_grid"'
    !write(1,*) ' variables= x, y'
    !write(1,*) ' zone t="grid", i=',n,'j=',m,', datapacking=point'
    write(10,*) cle+1,row+1
    do j=0,row
       do i=0,cle
          write(10,*) x(i,j), y(i,j)
       enddo
    enddo
    close(10)
    print *,"Finish Grid"
end program
