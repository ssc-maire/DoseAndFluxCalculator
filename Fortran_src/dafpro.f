!=========================================================================
! program: dafpro
!> @file
!>      This is a program to apply response matricies for several dose
!>      types and for neutrons to incoming proton spectra hitting Earth's atmosphere. 
!>      The response matricies were calculated using the Geant4 based magnetocosmics 
!>      simulation software.
!>
!> This documentation is up-to-date as of 03/03/2021
!
!=========================================================================

!-----------------------------------------------------------
!> @param[in] spectrum.p a file containing the incoming spectrum of
!> protons
!> @param[in] flightfile the file containing the coordinates where dose
!> rates will be calculated
!> @param[in] temp.p the path to the output file which will be generated
!> and contain dose rates for each of the specified locations in the
!> flightfile.
!>
!> in addition to using the above input arguments, this program also accesses
!> four files in /data/mcnpx/proton/ relative to the \b $MAIREHOME environment
!> variable, which must be set externally prior to running this program.
!> The files are adose.rpf, edose.rpf, dosee.rpf and neutron.rpf, which
!> contain the response matricies for different dose types at many
!> different altitudes (altitudes in this code are also referred to as
!> layers).
!>
!> example usage as written in gcrdoseflux.sh which the program is
!> currently called from: $3$2$1.p $4 temp.p
!-----------------------------------------------------------

      program dafpro
*     
      character*1000  filnam
      character*1000  path
      integer i, j, H, layer, nlayer
*     
      real flux(50),fl(50),f1,f2
      real energy(51)

      real redose, radose, rdosee
      real tn1, tn2, tn3
*
*     
c     get inputs from command line
c
      INTEGER      ARGC
      CHARACTER    ARGV*1000
c
c     response matrix
c
      real adose(137,50), edose(137,50), dosee(137,50), rpf(137,150)
!the x value of rpf here must have the same length as the neutron.rpf
!file
c
      ARGC = IARGC()
      if (argc.ne.3) then 
         print*,
     +  ' to use: filghpro flux_fname flightpath_fname output_fname'
         call exit
      endif
c
      CALL GETARG( 1, ARGV)
c     open the analysis parameter data file
c     
      open(11,file=argv,form='formatted',status='old')

c     flux file is in units of particles/cm2/sr/mev/s  
      do i = 1,50
         read(11,*) aa, flux(i)
      enddo
      close (11)
c
      do i = 1,51
         energy(i) = 10**(0.1*(i-1)+1)
      enddo
c     convert to per cm2.s
      do i = 1,50
         flux(i) = flux(i)*(energy(i+1)-energy(i))*3.14159
      enddo
*     
c     open and read in the response matrixs
      call getenv('MAIREHOME',path)
      !print *, "MAIREHOME is ", path
      idx = index(path,' ')-1
      filnam = path(1:idx)//'/data/mcnpx/proton/adose.rpf'
      !print *, "path is ", path(1:idx)
      !print *, "filnam is ", filnam
      open(21,file=filnam,form='formatted',status='old')
      filnam = path(1:idx)//'/data/mcnpx/proton/edose.rpf'
      open(22,file=filnam,form='formatted',status='old')
      filnam = path(1:idx)//'/data/mcnpx/proton/dosee.rpf'
      open(23,file=filnam,form='formatted',status='old')
      filnam = path(1:idx)//'/data/mcnpx/proton/neutron.rpf'
      open(24,file=filnam,form='formatted',status='old')
c
      do i = 1,137
         read(21,*)(adose(i,j),j=1,50)
         read(22,*)(edose(i,j),j=1,50)
         read(23,*)(dosee(i,j),j=1,50)
         read(24,*)(rpf(i,j),j=1,150) 
         !rpf has 137 rows and 150 columns, 50 columns for each neutron flux
         !component
      enddo
      close (21)
      close (22)
      close (23)
      close (24)
c
      CALL GETARG( 2, ARGV)
c     open the flight data file
c     
      open(12,file=argv,form='formatted',status='old')
c
      CALL GETARG( 3, ARGV)
c     open the output file
c     
      open(13,file=argv,form='formatted',status='unknown')
      write(13,"(2(A))") "time,latitude,longitude",
     + ",altitude,adose,edose,dosee,tn1,tn2,tn3,SEU,SEL"
c     write(13,'(A10,9(",",A10))') 'time','latitude','longitude',
c     + 'altitude','adose', 'edose', 'dosee', 'tn1', 'tn2', 'tn3'
c
c     loop over flight path file
      do while (.true.)   
         read(12,*,end=100) rtime,rlatitude,rlongitude,altitude,rc,fr
         call findlayer (altitude, layer, f1)
c         print *,altitude, layer, f1
         nlayer = layer + 1
         if (nlayer .gt. 137) nlayer = 137
         f2 = 1. - f1
c     
         do i = 1,50
            fl(i) = flux(i)*fr
         enddo
         ene = 1000.*(sqrt(rc**2+0.938**2)-0.938)
         ie = 1
         do while (ene .gt. energy(ie))
            ie = ie + 1
         enddo
         if(ie.ne.1) then
            ie = ie-1
            et = energy(ie+1) - energy(ie)
            ed = energy(ie+1) - ene
            fe = ed/et
            fl(ie)  = fl(ie)*fe
         endif
c     
c     calculate the doses and fluxes
         radose = 0.
         redose = 0.
         rdosee = 0.
         tn1 = 0.
         tn2 = 0.
         tn3 = 0.
c
         do i = ie,50
            radose = radose + 
     +       adose(layer,i)*f1*fl(i) + adose(nlayer,i)*f2*fl(i)
            redose = redose + 
     +       edose(layer,i)*f1*fl(i) + edose(nlayer,i)*f2*fl(i)
            rdosee = rdosee +
     +       dosee(layer,i)*f1*fl(i) + dosee(nlayer,i)*f2*fl(i)
            tn1 = tn1 +
     +       rpf(layer,i)*f1*fl(i) + rpf(nlayer,i)*f2*fl(i)
            tn2 = tn2 +
     +       rpf(layer,50+i)*f1*fl(i) + rpf(nlayer,50+i)*f2*fl(i)
            tn3 = tn3 +
     +       rpf(layer,100+i)*f1*fl(i) + rpf(nlayer,100+i)*f2*fl(i)
         enddo
         SEU = tn2 * 1e-13
         SEL = tn2 * 1e-8
c
         write(13,'(23(g0))') 
     +          rtime,',',
     +          rlatitude,',',
     +          rlongitude,',',
     +          altitude,',',
     +          radose,',',
     +          redose,',',
     +          rdosee,',',
     +          tn1,',',
     +          tn2,',',
     +          tn3,',',
     +          SEU,',',
     +          SEL
c         write(13,'(f20.0,9(",",f20.6))') rtime,rlatitude,rlongitude,
c     +      altitude, radose, redose, rdosee, tn1, tn2, tn3
c     
      enddo
 100  close (12)
      close (13)
*     
      end


!-----------------------------------------------------------
! subroutine findlayer
!> @param[in] height the altitude in km at which the relevant response matrix
!> is to be determined
!> @param[out] layer the index for the layer in response matrix files that is directly
!> below the altitude given to this function
!> @param[out] f1 a number defining how the actual response matrix for
!> the altitude should be interpolated from the response matrix at the layer directly
!> below the altitude given to this function.
!>
!> This subroutine is a simple subroutine to return the index for the pre-calculated
!> response matrix for the layer directly below a given altitude, and a
!> number f1 describing how to interpolate from the pre-calculated
!> response matrix to the actual response matrix for the actual altitude  
!-----------------------------------------------------------

      subroutine findlayer (height,layer, f1)
      integer H
      if (height .gt. 100.)  then 
         print *, ' Height has to be less than 100 km! '
         call exit()
      endif
      H = int(height*1000)
      if (height .lt. 0.025) then
         layer = 1
         f1 = 1
      else if (height .lt. 1.025) then
         layer = int((H-25)/50) + 1
         hr = mod(H-25,50)
         f1 = 1. - hr/50.
      else if (height .lt. 1.15) then
         layer = 21
         f1 = 1. - (H-1025)/1025.
      else if (height .lt. 5.05) then
         layer = int((H-1150)/100) + 22
         hr = mod(H-1150,100)
         f1 = 1. - hr/100.
      else if (height .lt. 5.3) then
         layer = 61
         f1 = 1. - (height-5.05)/0.25
      else if (height .lt. 15.1) then
         layer = int((H-5300)/200) + 62
         hr = mod((H-5300),200)
         f1 = 1. - hr/200.
      else if (height .lt. 16.5) then
         layer = 111
         f1 = 1. - (height-15.1)/1.4
      else if (height .lt. 38.5) then
         layer = int((H-16500)/1000) + 112
         hr = mod((H-16500),1000)
         f1 = 1. - hr/1000.
      else if (height .lt. 40.5) then
         layer = 134
         f1 = 1. - (height-38.5)/2.
      else if (height .lt. 62.5) then
         layer = 135
         f1 = 1. - (height-40.5)/22.
      else if (height .lt. 97.5) then
         layer = 136
         f1 = 1. - (height-40.5)/57.5
      else 
         layer = 137
         f1 = 1
      end if 

      end
