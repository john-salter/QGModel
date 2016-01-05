program QGModel
  
  use NETCDF
  
  implicit none
  
  real :: height, g, f, d, omega, deltat, cf, beta, tes, AH
  integer :: i, j, k, n, m, ncattrID, ncfileID, nxDimID, nyDimID, DimID, gID, hId, scaleID
  integer :: fID, bID, dID, tID, sID, nID, zID, oID, uID, vID, nhID, cfID, tauID
  integer, dimension(2) :: DimsOut 
  integer :: status
  real :: coef, rho, arr, Ekno, Kv
  integer, parameter :: l = 1000, mid = 500, quart = 250 
  real, parameter :: err = 0.005, pi = 3.141592653589, lscale = 1.0d-4
  logical :: stp
  real, dimension(l) :: x
  real, dimension(l,l) :: psi, psinext, q, qnext, h, j1, forc, U, V, tau, qlast, tauz, diffu
  integer :: numsteps, numsaved, saveskip
  real :: deltah
  character(37) :: attrfile
  character(38) :: ncfile
  character(7) :: x1
  character(len=8) :: fmt ! format descriptor
  
  fmt = '(I7.7)' ! an integer of width 7 with zeros at the left 

  ! load the attributes for the run
  attrfile = '/Volumes/jsfronts/QG/Attributes.nc'
  !*********************************************************************************************
  !This block will change run-to-run
  !*********************************************************************************************
  !height of the water column
  height = 500.0
  !gravitational constant
  g = 9.8
  !coriolis parameter
  f = 0.000073
  ! d is related to the width of the domain (in meters)
  d = 10000.0
  ! one of the diffusivity coefficients
  Kv = f*height
  !dx = 15000 
  !dy = 30000
  !omega = 1/2*f
  !beta is the rate of change of f with latitude
  beta = 1.98026d-11
  ! cf is phase speed of rossby waves, one of the parameters which limits the time step
  cf = sqrt(Kv*f/2.0)/height!1.0d-7
  ! AH is an optional parameter which tunes the strength of some frictional forcing
  AH = 0.0
  rho = 1024.0
  !(-omega**2)/2/g
  ! time step was determined by standard stability criteria
  deltat = 0.3*(f**2.0)*d/(g*height*beta)
  !numsteps is how long the model is to run
  numsteps = 60*60*24*30*12/deltat
  print*, 'dt is ', deltat, ' numsteps is ', numsteps
  ! alternative stability criteria
  !deltat = 0.0025*1/sqrt(f**2/2+8*sqrt(2.0)*g*height/d**2)
  !deltah = -0.000504*2
  !coef is a coefficient used in the method
  coef = cos(pi/(l-1))
  coef = 2*(1 - sqrt(1-coef))/coef
  do i=1, l
     x(i) = 0.0 + (i-1)*d
  end do
  !tauz and tau are related to the wind forcing of the basin
  do i = 1,l
     do j = 1,l
        tauz(j,i) = -0.0001*pi/(d*l*height)*sin(pi*x(j)/(d*l))
     end do
  end do
  tau = tauz
  !*********************************************************************************************
  !end of block
  !*********************************************************************************************
  !  initialize the stream functions and model variables
  psi = 0.0
  psinext = 0.0
  q = 0.0
  qnext = 0.0
  h = 0.0
  numsaved = 0
  ! saveskip allows the model to cut down on storage space while keeping a small time step
  saveskip = int(12.0*3600.0/deltat)
  print *, 'Images saved every ', deltat*saveskip/60/60, ' hours'

  status = nf90_create(attrfile, 0, ncattrID)

  !     Create the dimensions for the output fields.

  status = nf90_def_dim(ncattrID, 'n', 1, DimID)
  status = nf90_def_dim(ncattrID, 'nx', l, nxDimID)
  status = nf90_def_dim(ncattrID, 'ny', l, nyDimID)
  
  DimsOut(1) = nxDimID
  DimsOut(2) = nyDimID
  ! create all the variables for the output fields
  status = nf90_def_var(ncattrID, 'g', nf90_real,DimID, gID)
  status = nf90_def_var(ncattrID, 'beta', nf90_real,DimID, bID)
  status = nf90_def_var(ncattrID, 'H', nf90_real,DimID, hID)
  status = nf90_def_var(ncattrID, 'f', nf90_real,DimID, fID)
  status = nf90_def_var(ncattrID, 'd', nf90_real,DimID, dID)
  status = nf90_def_var(ncattrID, 'dt', nf90_real,DimID, tID)
  status = nf90_def_var(ncattrID, 'dh', nf90_real,DimID, zID)
  status = nf90_def_var(ncattrID, 'skip', nf90_int,DimID, sID)
  status = nf90_def_var(ncattrID, 'domain', nf90_real,DimID, oID)
  status = nf90_def_var(ncattrID, 'numsaved', nf90_int,DimID, nID)
  status = nf90_def_var(ncattrID, 'cf', nf90_real,DimID, cfID)
  status = nf90_def_var(ncattrID, 'Psi', nf90_real,DimsOut, vID)
  status = nf90_def_var(ncattrID, 'Tau', nf90_real,DimsOut, tauID)
  status = nf90_def_var(ncattrID, 'Q0', nf90_real,DimsOut, uID)
  status = nf90_def_var(ncattrID, 'linear-scale', nf90_real, DimID, scaleID)
  !define attributes for the new variables
  status = nf90_put_att(ncattrID, gID, 'long name', 'gravity coefficient')

  status = nf90_put_att(ncattrID, cfID, 'long name', 'coefficient of bottom friction')

  status = nf90_put_att(ncattrID, bID, 'long name', 'beta included in run')

  status = nf90_put_att(ncattrID, hID, 'long name', 'Height of water column')

  status = nf90_put_att(ncattrID, fID, 'long name', 'coriolis parameter')

  status = nf90_put_att(ncattrID, dID, 'long name', 'resolution')
  status = nf90_put_att(ncattrID, dID, 'units', 'm/pixel')

  status = nf90_put_att(ncattrID, tID, 'long name', 'time step')
  status = nf90_put_att(ncattrID, tID, 'units', 's/step')

  status = nf90_put_att(ncattrID, zID, 'long name', 'pump flow rate parameter')
  status = nf90_put_att(ncattrID, zID, 'units', 'm/step')

  status = nf90_put_att(ncattrID, sID, 'long name', 'number of files between saved data')
  status = nf90_put_att(ncattrID, sID, 'units', 'number of files')

  status = nf90_put_att(ncattrID, oID, 'long name', 'min/max of domain')
  status = nf90_put_att(ncattrID, oID, 'units', 'm')

  status = nf90_put_att(ncattrID, nID, 'long name', 'number of files in directory')
  status = nf90_enddef(ncattrID)
  psi = 0.0

  !*********************************************************************************************
  !This block will change run to run
  !*********************************************************************************************
  ! in the following loop, the stream function is set to its initial condition, The comments here are various test runs I did
  do i = 2,l-1!75,l-75!mid-62, mid+62
     do j = 2,l-1!75,l-75!mid-62, mid+62
        !psi(j,i) = 1.3439d5*cos(pi*(x(i)-500.0*d)/1.213d6)*cos(pi*(x(j)-500.0*d)/1.213d6)
        !if(psi(j,i).le.0.0) psi(j,i) = 0.0
        !psi(j,i) = lscale*exp(-(x(i)-500.0*d)**2.0/2.0/10000000000.0)*exp(-(x(j)-500.0*d)**2.0/2.0/10000000000.0)
     end do
  end do
  ! the following loops initialize the variable Q, which is related to the potential vorticity of the flow
  do i = 2, l-1
     do j = 2, l-1     
        q(j,i) = (psi(j+1,i) - 2.0*psi(j,i) + psi(j-1,i))/d**2.0 + (psi(j,i+1) - 2.0*psi(j,i) + psi(j,i-1))/d**2.0 &
             - f**2.0/(g*height)*psi(j,i) + beta*x(j)! + f/height*h(j,i)
     end do
  end do
  !*********************************************************************************************
  !end block
  !*********************************************************************************************

  !first, save off the initialized variables
  status = nf90_put_var(ncattrID, vID, psi)
  status = nf90_put_var(ncattrID, scaleID, lscale)
  status = nf90_put_var(ncattrID, uID, q)
  status = nf90_put_var(ncattrID, tauID, tau)
  status = nf90_put_var(ncattrID, bID, beta)
  status = nf90_put_var(ncattrID, hID, height)
  status = nf90_put_var(ncattrID, fID, f)
  status = nf90_put_var(ncattrID, dID, d)
  status = nf90_put_var(ncattrID, tID, deltat)
  status = nf90_put_var(ncattrID, zID, deltah)
  status = nf90_put_var(ncattrID, sID, saveskip)
  status = nf90_put_var(ncattrID, oID, x(l))
  status = nf90_put_var(ncattrID, gID, g)
  status = nf90_put_var(ncattrID, cfID, cf)
  status = nf90_put_var(ncattrID, nID, numsaved)
  status = nf90_close(ncattrID)
 
  qlast = q

  !*********************************************************************************************
  !MAIN LOOP
  !*********************************************************************************************
  do n = 1, numsteps  
     if(n < 9010) then
        !this section is for a variable wind forcing field. In one run, it was set to increase after 9010 time steps
        !if(mod(n,100)==1) tau = tau + tauz
     end if     
     j1 = 0.0
     !this subroutine computes the jacobian of psi and q
     call Jac(psi, q, j1, l, d)
     !this subroutine computes the frictional terms for the stream function psi
     call fric(psi, l, diffu)
     if(mod(n,50).ne.0) then
        !I was playing around with running every 50th time step with more friction and/or a different method
        !so I could run the method with more stability for a given time step
        do i = 2,l-1
           do j = 2,l-1
              ! this is the main loop of the numerical method
              qnext(j,i) = qlast(j,i) - 2.0*deltat*j1(j,i) & 
                    - 2.0*cf*deltat*((psi(j+1,i) - 2.0*psi(j,i) + psi(j-1,i))/d**2.0 &
                   + (psi(j,i+1) - 2.0*psi(j,i) + psi(j,i-1))/d**2.0) + 2.0*deltat*tau(j,i) - 2.0*deltat*AH*diffu(j,i)/d**4
           end do
        end do
     else
        do i = 2,l-1
           do j = 2,l-1
              !this is the loop the program takes once every 50 iterations
              qnext(j,i) = q(j,i) - deltat*j1(j,i) & 
                   - cf*deltat*((psi(j+1,i) - 2.0*psi(j,i) + psi(j-1,i))/d**2.0 &
                   + (psi(j,i+1) - 2.0*psi(j,i) + psi(j,i-1))/d**2.0) + deltat*tau(j,i) - deltat*AH*diffu(j,i)/d**4
           end do
        end do
     end if
     qlast = q
     !q = 0.9*qnext + .1*q
     q = qnext
     qnext = 0.0
     ! initialize the while loop
     stp = .false.
     m = 1
     !initialize the forcing matrix
     do i = 1,l
        do j = 1,l
           forc(j,i) = q(j,i) - beta*x(j)! - f/height*h(j,i)
        end do
     end do
     ! this while loops solves the stream function iteratively by using the Q determined above
     do while(.not.stp)
        psinext = 0.0 
        do i = 1, l
           psinext(1,i) = 0.0
           psinext(l,i) = 0.0
           psinext(i,1) = 0.0
           psinext(i,l) = 0.0
        end do
        do j = 2,l-1
           do i = 2,l-1
              ! I was playing around with various methods of iteration, and the pros/cons of each
              !              psinext(j,i) = ((psi(j+1,i) + psi(j-1,i) + psi(j,i+1) + psi(j,i-1))/d**2.0 - (forc(j,i)))/&
              !(4.0/d**2.0+f**2.0/(g*height))!&
              psinext(j,i) = (psi(j+1,i) + psinext(j-1,i) + psi(j,i+1) + psinext(j,i-1))/4.0 - d**2.0/4.0*(forc(j,i)) &
              -d**2*psi(j,i)*f**2.0/(4.0*g*height)
              !print*, (1+d**2.0*f**2.0/(4.0*g*height))
!                   q(j,i) - beta*x(j) + f**2.0/(g*height)*psi(j,i))
              !psinext(j,i) = psi(j,i) + coef/4.0*(psi(j,i+1) + psinext(j,i-1) + (psi(j+1,i) - psinext(j-1,i)) &
               !    - 4.0*psi(j,i) - d**2.0*forc(j,i))
           end do
        end do
        stp = .true.
        ! the poorly named mee subroutine compares the iterated psi with psinext and determines the difference between the two, called arr
        call mee(psi,psinext,arr,l)
        !then, if arr is still higher than the tolerance, set the loop control variable to keep going
        if(arr.ge.err) stp = .false.

        ! every hundred iterations, I printed out the error to make sure it was converging
        if(mod(m,100)==0) print*, m, ' Error is ', arr
        if(m > 10000) then
           !if it hadn't converged by 10k iterations, something was wrong and I needed to know to debug it
           !usually it was a missed sign, or an unstable method (because time step was too large)
           stp = .true.
           print*, 'WARNING! Stream failed to converge in 10000 iterations! Error was ', tes
        end if                
        psi = psinext
        m = m + 1
     end do
     ! this is where I handled the boundary conditions
     do i = 1, l
        psi(1,i) = 0.0
        psi(l,i) = 0.0
        psi(i,1) = 0.0
        psi(i,l) = 0.0
        psinext(1,i) = 0.0
        psinext(l,i) = 0.0
        psinext(i,1) = 0.0
        psinext(i,l) = 0.0
     end do
     if(mod(n,saveskip) == 1) then
        !initialize the velocities
        U = 0.0
        V = 0.0
        do i = 2,l-1
           do j = 2, l-1
              !calculate the velocities that are about to be saved out
              U(j,i) = -(psi(j+1,i)-psi(j-1,i))/2.0/d
              V(j,i) = (psi(j,i+1)-psi(j,i-1))/2.0/d
           end do
        end do
        !here, I prepare to save off the nth iteration with the iteration number as the file name
        write (x1,fmt) n ! converting integer to string using an 'internal file'
        
        ncfile='/Volumes/jsfronts/QG/Data'//x1//'.nc'
        !this section is identical to the file created above for initial conditions
        status = nf90_create(ncfile, 0, ncfileID)
        
        !     Create the dimensions for the output fields.
        
        status = nf90_def_dim(ncfileID, 'nx', l, nxDimID)
        status = nf90_def_dim(ncfileID, 'ny', l, nyDimID)
        
        DimsOut(1) = nxDimID
        DimsOut(2) = nyDimID
        
        status = nf90_def_var(ncfileID, 'Psi', nf90_real,DimsOut, hID)
        status = nf90_def_var(ncfileID, 'U', nf90_real,DimsOut, uID)
        status = nf90_def_var(ncfileID, 'V', nf90_real,DimsOut, vID)
        
        status = nf90_put_att(ncfileID, uID, 'long name', 'Stream function')        
        
        status = nf90_enddef(ncfileID)
        
        status = nf90_put_var(ncfileID, hID, psi)
        status = nf90_put_var(ncfileID, uID, U)
        status = nf90_put_var(ncfileID, vID, V)
        
        status = nf90_close(ncfileID)
        
        numsaved = numsaved + 1
        status = nf90_open(attrfile, NF90_WRITE, ncattrID)
        status = nf90_put_var(ncattrID, nID, numsaved)
        status = nf90_close(ncattrID)
        print *, ' ...', numsaved, ' files saved so far'
        print *, n, ' of ', numsteps, ' steps so far', ' ... this time converged in ', m, ' steps with error ', arr  
     endif     
  end do
contains

  subroutine Jac(A,B,jacobian,c,d1)
    !this subroutine returns the jacobian of A and B
  implicit none
  integer, intent(in) :: c
  integer :: k, l
  real, intent(in) :: d1
  real, dimension(c,c), intent(in) :: A, B
  real, dimension(c,c), intent(out) :: jacobian
  real, dimension(c,c) :: J1, J2, J3
  jacobian = 0.0
  J1 = 0.0
  J2 = 0.0
  J3 = 0.0

  do l = 2, c-1
     do k = 2, c-1
        J1(k,l) = (A(k,l+1)-A(k,l-1))*(B(k+1,l)-B(k-1,l))/4.0/d1**2.0 - (B(k,l+1)-B(k,l-1))*(A(k+1,l)-A(k-1,l))/4.0/d1**2.0
        J2(k,l) = -(B(k,l+1)*(A(k+1,l+1)-A(k-1,l+1))/2.0/d1 - B(k,l-1)*(A(k+1,l-1)-A(k-1,l-1))/2.0/d1)/2.0/d1 &
             + (B(k+1,l)*(A(k+1,l+1)-A(k+1,l-1))/2.0/d1 - B(k-1,l)*(A(k-1,l+1)-A(k-1,l-1))/2.0/d1)/2.0/d1
        J3(k,l) = (A(k,l+1)*(B(k+1,l+1)-B(k-1,l+1))/2.0/d1 - A(k,l-1)*(B(k+1,l-1)-B(k-1,l-1))/2.0/d1)/2.0/d1 &
             - (A(k+1,l)*(B(k+1,l+1)-B(k+1,l-1))/2.0/d1 - A(k-1,l)*(B(k-1,l+1)-B(k-1,l-1))/2.0/d1)/2.0/d1
     end do
  end do

  jacobian = J1/3.0 + J2/3.0 + J3/3.0
end subroutine Jac

subroutine mee(A,B,me,c)
  !this subroutine gives the difference between A and B for error tolerance above
  implicit none
  integer :: i, j, m
  integer, intent(in) :: c
  real, dimension(c,c), intent(in) :: A,B
  real, dimension(c,c) :: D 
  real, intent(out) :: me
  real :: su
  me = 0.0
  m = 0
  su = 0.0
  D = 0.0
  do j = 1,c
     do i = 1,c
        if(A(j,i).ne.0.0) then
           D(j,i) = abs((A(j,i)-B(j,i))/A(j,i))
        else
           D(j,i) = -1.0
        end if
     end do
  end do
  do i = 1,c
     do j = 1,c
        if(D(j,i) .ge. 0.0)then
           su = su + D(j,i)
           m = m+1
        end if
     end do
  end do
  if(m.ne.0)me = su/real(m)
  
end subroutine mee
subroutine fric(A,c,fr)
  !This subroutine returns the friction terms for the differential equation. Nothing special, it was just kind of complicated so I
  !separated it out for easy debugging
  implicit none
  integer, intent(in) :: c
  integer :: jj, ii
  real, dimension(c,c), intent(in) :: A
  real, dimension(c,c), intent(out) :: fr
  fr = 0.0
  do jj = 3,c-2
     do ii = 3,c-2
        fr(jj,ii) = (A(jj,ii+2)-4.0*A(jj,ii+1)-4.0*A(jj,ii-1)+A(jj,ii-2)) + (A(jj+2,ii)-4.0*A(jj+1,ii)-4.0*A(jj-1,ii)+A(jj-2,ii)) &
             + A(jj-1,ii+1) - 2.0*A(jj-1,ii)+A(jj-1,ii-1) -2.0*(A(jj+1,ii) - 2.0*A(jj,ii)+A(jj-1,ii))+A(jj+1,ii-1)-2.0*A(jj,ii-1)&
             + A(jj-1,ii-1)
     end do
  end do
end subroutine fric
end program QGModel
