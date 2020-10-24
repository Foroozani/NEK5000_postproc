       program hpts_setup_rack
      implicit none
!------------------------------------------------------------------
!                      
!      To set hvirtual sensors along each beam line 
!      at 4 different heights, similar to the experinmental setup 
!      ifort -O2 -i4 -r8 hpts_setup_rack.f90 -o hpts_setup_rack.exe
!      ./hpts_setup_rack.exe my_file.his
!-----------------------------------------------------------------
      ! number of history points per beam
      integer*4, parameter :: npts_beam = 40
      ! number of points along the height
      integer*4, parameter :: npts_zz = 6
      real*8, parameter :: pi = 3.1415926535897932384626433832795d0

      ! geometry parameters
      real*8, save :: Radius = 0.25
      real*8, save :: Height = 1.00
      real*8, save :: Angle(1:4) = (/0.0, 45.0, 90.0, 135.0/)
      real*8, save :: Theta = 0.0 ! shift in degrees

      ! coordinates of history points
      real*8 :: xpos,ypos,zpos,dz

      integer*4 :: npts,ierr
      integer*4 :: i,j,k

      character*80 :: nek_his_file,fname

! --- get base name from command line ---
      call getarg(1,nek_his_file)

      open(unit=80, file=TRIM(nek_his_file), status='replace')

! --- start loop along the height ---
      dz = Height/float(npts_zz-1)
      zpos = 0.0

      do j = 1,npts_zz
         call hpts_plane_xy(npts_beam)
         zpos = zpos + dz
      enddo

      !write(*,*) 'zpos=',zpos

      close(80)

! --- for now, stop program here ---
      stop 'Finished !'

      contains

! --- create points in 2D plane ---
      subroutine hpts_plane_xy(np)
      ! dummy arguments
      integer*4 :: np ! npts per beam

      ! local variables
      integer*4 :: i,j,k
      real*8 :: d,r,dr,phi,cos_phi,sin_phi

      !write(80,'(A,1X,es13.6)') '# zpos=',zpos

      d = 2.0*Radius ! diameter
      dr = d/float(np-1)

! --- loop over phi-angles ---
      do k = 1,SIZE(Angle)
         phi = (Angle(k) + Theta) * pi/180.0

         cos_phi = cos(phi)
         if (abs(cos_phi).le.1.0d-9) cos_phi = 0.0d0

         sin_phi = sin(phi)
         if (abs(sin_phi).le.1.0d-9) sin_phi = 0.0d0

         r = Radius
         do i = 1,np
            xpos = r * cos_phi
            ypos = r * sin_phi

            write(80,'(3(1X,es13.6))') xpos,ypos,zpos

            r = r - dr
         enddo
      enddo

      !write(80,*) ! skip line

      return
      end subroutine hpts_plane_xy


   end program hpts_setup_rack
