!=====================================================================
!=====================================================================
!
!
!
!  halo_defs.f90
!
!
!
!=====================================================================
!=====================================================================
module fhalo_defs

  public

   real(kind=8) :: memory_used ! in GB
  !======================================================================
  ! parameters relative to the simulation analysis
  !======================================================================
!   character(len=8),parameter :: gravsoft  = 'cubsplin'         ! type of gravitational softening 
  integer(kind=4)            :: nbPes                          ! obsolete vars for reading treecode format
!   integer(kind=4)            :: nsteps                         ! number of timesteps
!   real(kind=8)               :: alpha,tnow                     ! 
  integer(kind=4)            :: nusedpart                        ! number of particles in the simulation 
  integer(kind=4)            :: nMembers                       ! minimal number of particles in a fof halo
!   real(kind=8)               :: b_init                         ! linking length parameter of the fof at z=0
!   character(len=4),parameter :: profile   = 'TSIS'             ! type of halo profile (only isothermal_sphere yet)
!   integer(kind=4),parameter  :: ninterp   = 30000              ! nb of bins for interp. of smoothed grav. field
!   integer(kind=4)            :: FlagPeriod                     ! flag for periodicity of boundary conditions 
!   real(kind=8)               :: fPeriod(3)
  !----------- For gadget format: ----------------------------------------
!   integer (kind=4)           :: nhr                            ! to read only the selected HR particles
!   real (kind=4)              :: minlrmrat                      ! to recognize contaminated halos
  !======================================================================


  !======================================================================
  ! Definitions specific to input/output
  !======================================================================
!   character(len=80)         :: output_dir
!   character(len=5)          :: file_num
!   integer(kind=4)           :: numstep
  integer(kind=4),parameter :: errunit = 0
!   logical(kind=4) :: write_resim_masses                        ! for writing resim_masses.dat file
!   integer(kind=4)           :: dchmod = int(o'755')
!   character(len=5)          :: dchmod_str = '755'
!   integer(kind=4)           :: fchmod = int(o'644')
!   character(len=5)          :: fchmod_str = '644'
!   integer(kind=4)           :: uid = -1
!   integer(kind=4)           :: gid = -1
  !======================================================================


  !======================================================================
  ! Constants
  !======================================================================
!   real(kind=8),parameter    :: gravconst = 430.1               ! G in units of (km/s)^2 Mpc/(10^11 Msol)
  real(kind=8),parameter    :: pi        = 3.141592654
  !======================================================================


  !======================================================================
  ! Global variables 
  !======================================================================
  real(kind=8),allocatable         :: pos(:,:)!,vel(:,:)
  real(kind=8)                     :: massp
!   real(kind=8),allocatable         :: epsvect(:)
  real(kind=8),allocatable         :: mass(:)
  logical, allocatable             :: refmask(:)
  real(kind=8)                     :: omega_t,omega_lambda_t,omega_f,omega_lambda_f,omega_c_f
  real(kind=8)                     :: rho_crit,aexp,Lboxp,mboxp,af,ai,Lf,H_f,H_i
  real(kind=8)                     :: age_univ,Lbox_pt,Lbox_pt2,Hub_pt,omega_0,hubble,omega_lambda_0
!   real(kind=8)                     :: vir_overdens,rho_mean
!   integer(kind=4),allocatable      :: linked_list(:)!, liste_parts(:)
!   integer(kind=4),allocatable      :: first_part(:),nb_of_parts(:)
!   type (halo),allocatable          :: liste_halos(:)
!   real(kind=8)                     :: phsmooth(0:1+ninterp)
!   integer(kind=4)                  :: nb_of_halos, nb_of_subhalos
!   integer(kind=4)                  :: numero_step  
  character(len=3)                 :: type
  !======================================================================


  !======================================================================
  ! defs for Adaptahop
  !======================================================================
   ! integer(kind=4), parameter :: nparmax=512*512*512
   ! integer(kind=4),parameter  :: nparbuffer=128**3
   integer(kind=4),parameter  :: ncellbuffermin=128**3
   ! integer(kind=4),parameter  :: nlevelmax=30,npartpercell=100
   integer(kind=4),parameter  :: npartpercell=100
   ! integer(kind=4)            :: lin=10,lsin=11,lin2=12
   ! integer(kind=4)            :: loudis=12,lounei=13,lounode=14,loupartnode=15,lounodedyn=16
   real(kind=8)            :: bignum=1.d30
   ! Physical constants (units : m s kg) ->
   real(kind=8), parameter :: gravitational_constant=6.6726d-11
   real(kind=8), parameter :: critical_density= 1.8788d-26
   real(kind=8), parameter :: mega_parsec=3.0857d22
   real(kind=8), parameter :: solar_mass=1.989d30
   real(kind=8), parameter :: convert_in_mps2=1.d6
   ! real(kind=8), allocatable    :: vxout(:),vyout(:),vzout(:),vdisout(:)
   integer(kind=4), allocatable :: mass_cell(:)
   ! real(kind=8), allocatable    :: tmass_cell(:)
   ! real(kind=8), allocatable    :: vcell(:,:)
   real(kind=8), allocatable    :: size_cell(:)
   real(kind=8), allocatable    :: pos_cell(:,:)
   integer(kind=4), allocatable :: sister(:)
   integer(kind=4), allocatable :: firstchild(:)
   integer(kind=4), allocatable :: idpart(:),idpart_tmp(:),idpart_adapt(:)
   integer(kind=4)              :: nactive_parts=0
   real(kind=8), allocatable    :: distance(:)
   ! real(kind=8), allocatable    :: density(:)
   integer(kind=4), allocatable :: firstpart(:)
   integer(kind=4), allocatable :: idgroup(:),idgroup_tmp(:)
   integer(kind=4), allocatable :: igroupid(:)
   integer(kind=4), allocatable :: color(:)
   ! integer(kind=4), allocatable :: partnode(:)
   real(kind=8), allocatable    :: densityg(:)
   real(kind=8), allocatable    :: zoombox(:)
   real(kind=8)    :: sizeroot
   real(kind=8)    :: xlong, ylong, zlong, boxsize, boxsize2, xyzlong
   real(kind=8)    :: xlongs2, ylongs2, zlongs2
   real(kind=8)    :: omega0,omegaL,mass_in_kg,GMphys
   real(kind=8)    :: aexp_hubble_const,aexp_mega_parsec_sqrt3
   real(kind=8)    :: aexp_mega_parsec,aexp_max
   ! integer(kind=4) :: npart
   integer(kind=8) :: npart
   integer(kind=4) :: nvoisins,nhop,ntype,nlevelmax
   ! integer(kind=4) :: ncellmx
   integer(kind=8) :: ncellmx
   integer(kind=4) :: ngroups,nmembthresh,nnodes,nnodesmax
   ! integer(kind=4) :: ncpu,nmpi,niterations
   ! integer(kind=4) :: ncellbuffer
   integer(kind=8) :: ncellbuffer
   real(kind=8)    :: rho_threshold
   logical         :: verbose,megaverbose,periodic,zoomin
   ! real(kind=8)    :: fgas
   real(kind=8)    :: fudge,alphap,epsilon,fudgepsilon
   real(kind=8)    :: pos_shift(3),pos_renorm,velrenorm

   integer(kind=4), allocatable :: dev_group_nhnei(:)
   integer(kind=4), allocatable :: dev_group_offset(:)
   integer(kind=4), allocatable :: dev_isad_all(:)
   real(kind=8),    allocatable :: dev_rho_saddle_all(:)

   type supernode
      sequence
      integer(kind=4) :: mass
      integer(kind=4) :: level
      integer(kind=4) :: mother
      integer(kind=4) :: firstchild
      integer(kind=4) :: nsisters
      integer(kind=4) :: sister
      real(kind=8)    :: rho_saddle
      real(kind=8)    :: density
      real(kind=8)    :: densmax
      real(kind=8)    :: radius
      real(kind=8)    :: truemass
      real(kind=8)    :: position(3)
   end type supernode
   type (supernode), allocatable :: node(:)

   
   

!======================================================================
! Flags for halo finder selection
!======================================================================
   character(len=3)    :: method       ! flag to notify which and how the halofinder is to be used
   logical             :: fsub         ! flag to notify whether subhaloes are included               
   logical             :: cdm          ! flag to select particle closest to the cdm instead of the one with the highest density
   logical             :: DPMMC=.false.! flag to select the densest particle in the most massive cell of the halo (not with FOF)
   logical             :: SC=.false.   ! flag to select the com within concentric spheres (not with FOF)
   real(kind=8)        :: dcell_min
   real(kind=8)        :: eps_SC
   logical             :: dump_dms=.false.

!======================================================================
! array to build the structure tree
!======================================================================
integer(kind=4), allocatable :: first_daughter(:), mother(:), first_sister(:), level(:)
integer(kind=4) :: nstruct

! used for the merger history method
! integer(kind=4), allocatable :: npfather(:),ex_liste_parts(:),removesub(:)
! integer(kind=4), allocatable :: ex_level(:),ex_nb_of_parts(:),ex_linked_list(:)
! integer(kind=4) ::  ex_nb_of_structs


! #ifdef ANG_MOM_OF_R
! integer(kind=4),parameter :: agor_unit = 19
! character(200)            :: agor_file
! integer(kind=4),parameter :: nshells = 100
! #endif

! logical :: isprinting=.false.
! integer :: tqdm_icalced=0

contains


!=====================================================================
! num_rec.f90
!=====================================================================
subroutine indexx(n,arr,indx)
   implicit none

   integer(kind=4)           :: n,indx(n)
   real(kind=8)              :: arr(n)
   integer(kind=4),parameter :: m=7,nstack=50
   integer(kind=4)           :: i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
   real(kind=8)              :: a

   do j = 1,n
      indx(j) = j
   enddo

   jstack = 0
   l      = 1
   ir     = n
   1 if (ir-l .lt. m) then
      do j = l+1,ir
         indxt = indx(j)
         a     = arr(indxt)
         do i = j-1,1,-1
            if (arr(indx(i)) .le. a) goto 2
            indx(i+1) = indx(i)
         enddo
         i         = 0
   2       indx(i+1) = indxt
      enddo
      if (jstack .eq. 0) return
      ir     = istack(jstack)
      l      = istack(jstack-1)
      jstack = jstack-2
   else
      k         = (l+ir)/2
      itemp     = indx(k)
      indx(k)   = indx(l+1)
      indx(l+1) = itemp
      if (arr(indx(l+1)) .gt. arr(indx(ir))) then
         itemp     = indx(l+1)
         indx(l+1) = indx(ir)
         indx(ir)  = itemp
      endif
      if (arr(indx(l)) .gt. arr(indx(ir))) then
         itemp    = indx(l)
         indx(l)  = indx(ir)
         indx(ir) = itemp
      endif
      if (arr(indx(l+1)) .gt. arr(indx(l))) then
         itemp     = indx(l+1)
         indx(l+1) = indx(l)
         indx(l)   = itemp
      endif
      i     = l+1
      j     = ir
      indxt = indx(l)
      a     = arr(indxt)
   3    continue
      i     = i+1
      if (arr(indx(i)) .lt. a) goto 3
   4    continue
      j     = j-1
      if (arr(indx(j)) .gt. a) goto 4
      if (j .lt. i) goto 5
      itemp   = indx(i)
      indx(i) = indx(j)
      indx(j) = itemp
      goto 3
   5    continue
      indx(l) = indx(j)
      indx(j) = indxt
      jstack  = jstack+2
      if (jstack .gt. nstack) stop 'nstack too small in indexx'
      if (ir-i+1 .ge. j-l) then
         istack(jstack)   = ir
         istack(jstack-1) = i
         ir               = j-1
      else
         istack(jstack)   = j-1
         istack(jstack-1) = l
         l                = i
      endif
   endif
   goto 1

end subroutine indexx


subroutine print_memory()
!=======================================================================
   implicit none
   if (verbose) write(errunit,*) '--- Memory used (in GB) :',memory_used
   if (memory_used .gt. 200.d0) then
      write(errunit,*) 'WARNING: Memory usage is very high. Consider increasing the size of integer variables to avoid overflow.'
      write(errunit,*) 'WARNING: Memory usage is very high. Consider increasing the size of integer variables to avoid overflow.'
      write(errunit,*) 'WARNING: Memory usage is very high. Consider increasing the size of integer variables to avoid overflow.'
      write(errunit,*) 'WARNING: Memory usage is very high. Consider increasing the size of integer variables to avoid overflow.'
      write(errunit,*) 'WARNING: Memory usage is very high. Consider increasing the size of integer variables to avoid overflow.'
      STOP
   endif
end subroutine print_memory

end module fhalo_defs























































module neiKDtree
   use fhalo_defs
   use omp_lib
   private

   public :: compute_adaptahop, sync_others, sync_from_change_pos, sync_from_init_adaptahop, &
             close, liste_parts, density, nnodes, real_table, integer_table

   ! KDtree variables
   integer(kind=4), allocatable :: icidpart_tmp(:)
   real(kind=8), parameter :: pos_ref_kdtree(3,0:7) = reshape( &
                           & [ -1.d0,-1.d0,-1.d0, &
                           &    1.d0,-1.d0,-1.d0, &
                           &   -1.d0, 1.d0,-1.d0, &
                           &    1.d0, 1.d0,-1.d0, &
                           &   -1.d0,-1.d0, 1.d0, &
                           &    1.d0,-1.d0, 1.d0, &
                           &   -1.d0, 1.d0, 1.d0, &
                           &    1.d0, 1.d0, 1.d0 ], [3,8] )
   integer(kind=4), allocatable :: queue_color(:)

   ! Output to python
   integer(kind=4),allocatable      :: liste_parts(:)
   real(kind=8),   allocatable      :: density(:)

   ! Convert node(:) to python-understandable array
   real(kind=8),    dimension(:,:),   allocatable :: real_table
   integer(kind=4), dimension(:,:),   allocatable :: integer_table

   contains

!=======================================================================
subroutine compute_adaptahop(pos_in, mass_in)
!=======================================================================
   implicit none
   real(kind=8), intent(in)         :: pos_in(:,:)
   real(kind=8), intent(in)         :: mass_in(:)

   memory_used = 0d0
   write(errunit,*) "[Memory] Initial"
   call print_memory()

   allocate(liste_parts(1:nusedpart))
   memory_used = memory_used + real(size(liste_parts), kind=8)*4/1.d9
   ! allocate(pos(1:npart,1:3))
   if (size(pos_in,1) .ne. 3) then
      write(errunit,*) 'ERROR: pos_in should have shape (3,npart)'
      stop
   endif
   allocate(pos(1:3,1:npart))
   memory_used = memory_used + real(size(pos), kind=8)*8/1.d9
   pos(:,:) = pos_in(:,:)
   allocate(mass(npart))
   memory_used = memory_used + real(size(mass), kind=8)*8/1.d9
   mass(:) = mass_in(:)
   write(errunit,*) "[Memory] Start create_tree_structure"
   call print_memory()

   call create_tree_structure
   write(errunit,*) "[Memory] Start compute_mean_density_and_np"
   call print_memory()
   call compute_mean_density_and_np
   write(errunit,*) "[Memory] Start find_local_maxima"
   call print_memory()
   call find_local_maxima
   write(errunit,*) "[Memory] Start create_group_tree"
   call print_memory()
   call create_group_tree
   write(errunit,*) "[Memory] Done create_group_tree"
   call print_memory()

   memory_used = memory_used - real(size(pos), kind=8)*8/1.d9
   deallocate(pos)
   memory_used = memory_used - real(size(mass), kind=8)*8/1.d9
   deallocate(mass)
   write(errunit,*) "[Memory] Done compute_adaptahop"
   call print_memory()
   
   return
end subroutine compute_adaptahop



!=====================================================================
! For f2py
!=====================================================================
subroutine node2table
   integer(kind=4) :: inodee
   allocate(real_table(8,nnodes)) !  rho_saddle, density, densmax, radius, truemass, position(3)
   memory_used = memory_used + real(size(real_table), kind=8)*8/1.d9
   allocate(integer_table(6,nnodes)) ! mass, level, mother, firstchild, nsisters, sister
   memory_used = memory_used + real(size(integer_table), kind=8)*4/1.d9
   do inodee=1,nnodes
      real_table(1,inodee) = node(inodee)%rho_saddle
      real_table(2,inodee) = node(inodee)%density
      real_table(3,inodee) = node(inodee)%densmax
      real_table(4,inodee) = node(inodee)%radius
      real_table(5,inodee) = node(inodee)%truemass
      real_table(6:8,inodee) = node(inodee)%position(:)
      integer_table(1,inodee) = node(inodee)%mass
      integer_table(2,inodee) = node(inodee)%level
      integer_table(3,inodee) = node(inodee)%mother
      integer_table(4,inodee) = node(inodee)%firstchild
      integer_table(5,inodee) = node(inodee)%nsisters
      integer_table(6,inodee) = node(inodee)%sister
   end do
   memory_used = memory_used - real(size(node)*((4*6)+(8*8)), kind=8)/1.d9
   deallocate(node)
      
end subroutine node2table


subroutine sync_from_change_pos(npart_in, nusedpart_in, epsilon_in, fudgepsilon_in, xlong_in, boxsize2_in)
   implicit none

   integer(kind=4), intent(in) :: nusedpart_in
   integer(kind=8), intent(in) :: npart_in
   real(kind=8), intent(in) :: epsilon_in, fudgepsilon_in, xlong_in, boxsize2_in

   npart = int(npart_in, kind=8)
   nusedpart = nusedpart_in
   epsilon = epsilon_in
   fudgepsilon = fudgepsilon_in
   xlong = xlong_in
   boxsize2 = boxsize2_in
end subroutine sync_from_change_pos

subroutine sync_from_init_adaptahop( &
            npart_in,nmembthresh_in,nMembers_in, zoomin_in, &
            omegaL_in, omega_lambda_f_in,  &
            omega0_in, omega_f_in,  &
            aexp_max_in, af_in,  &
            hubble_in, H_f_in,  &
            boxsize2_in, Lf_in,  &
            xlong_in, ylong_in, zlong_in,  &
            xlongs2_in, ylongs2_in, zlongs2_in, &
            Hub_pt_in,aexp_in, &
            pos_shift_in)
   implicit none

   integer(kind=4), intent(in) :: nmembthresh_in,nMembers_in, zoomin_in
   integer(kind=8), intent(in) :: npart_in
   real(kind=8), intent(in) :: omegaL_in, omega_lambda_f_in,  &
            omega0_in, omega_f_in,  &
            aexp_max_in, af_in,  &
            hubble_in, H_f_in,  &
            boxsize2_in, Lf_in,  &
            xlong_in, ylong_in, zlong_in,  &
            xlongs2_in, ylongs2_in, zlongs2_in, &
            Hub_pt_in,aexp_in
   real(kind=8), intent(in) :: pos_shift_in(3)

   npart = int(npart_in, kind=8)
   nmembthresh = nmembthresh_in
   nMembers = nMembers_in
   zoomin = .false.
   if (zoomin_in.gt.0) zoomin = .true.
   omegaL = omegaL_in
   omega_lambda_f = omega_lambda_f_in
   omega0 = omega0_in
   omega_f = omega_f_in
   aexp_max = aexp_max_in
   af = af_in
   hubble = hubble_in
   H_f = H_f_in
   boxsize2 = boxsize2_in
   Lf = Lf_in
   xlong = xlong_in
   ylong = ylong_in
   zlong = zlong_in
   xyzlong = xlong*ylong*zlong
   xlongs2 = xlongs2_in
   ylongs2 = ylongs2_in
   zlongs2 = zlongs2_in
   Hub_pt = Hub_pt_in
   aexp = aexp_in
   pos_shift(1:3) = pos_shift_in(1:3)
   fsub=.true.

end subroutine sync_from_init_adaptahop

subroutine sync_others( &
   verbose_in, npart_in, nbPes_in,&
   rho_threshold_in, massp_in, boxsize_in, &
   nhop_in, nvoisins_in, &
   fudge_in, alphap_in, &
   method_in, nlevelmax_in)

   implicit none

   logical, intent(in) :: verbose_in
   integer(kind=4), intent(in) :: npart_in, nbPes_in
   integer(kind=4), intent(in) :: nhop_in, nvoisins_in
   real(kind=8), intent(in) :: rho_threshold_in, massp_in, boxsize_in
   real(kind=8), intent(in) :: fudge_in, alphap_in
   character(len=*), intent(in) :: method_in
   integer(kind=4), intent(in) :: nlevelmax_in

   verbose         = verbose_in
   npart           = npart_in
   nbPes           = nbPes_in
   rho_threshold   = rho_threshold_in
   massp           = massp_in
   boxsize         = boxsize_in

   nhop            = nhop_in
   nvoisins        = nvoisins_in
   fudge           = fudge_in
   alphap          = alphap_in
   
   method          = method_in(1:min(len(method), len(method_in)))
   nlevelmax       = nlevelmax_in

end subroutine sync_others


!=======================================================================
subroutine compute_mean_density_and_np
!=======================================================================
   implicit none

   integer(kind=4)                     :: ipar,nactive
   real(kind=8), dimension(0:nvoisins) :: dist2
   integer, dimension(nvoisins)        :: iparnei
   real(kind=8)                        :: densav
   integer(kind=4) :: tttt0, tttt1, ttttrate
   real(kind=8)    :: dtdtdtdt
   integer(kind=4) :: vt0, vt1, vtrate
   real(kind=8)    :: dvt

   call system_clock(count=tttt0, count_rate=ttttrate)
   if (verbose) write(errunit,*) '    Compute mean density for each particle...'

   allocate(density(npart))
   memory_used = memory_used + real(size(density), kind=8)*8/1.d9
   nactive=0
   call print_memory()


   call omp_set_num_threads(nbPes)
   if(verbose) write(errunit,*) "    [OMP] compute density with ncore=",nbPes

   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipar,dist2,iparnei)
   !$OMP DO
   do ipar=1,npart
      call omp_find_nearest_parts(ipar,dist2,iparnei)
      call omp_compute_density(ipar,dist2,iparnei)
      if (density(ipar).gt.rho_threshold) then
         !$OMP ATOMIC
         nactive=nactive+1
      endif
   enddo
   !$OMP END DO
   !$OMP END PARALLEL

   ! Check for average density
   if (verbose) then
      call system_clock(count=vt0, count_rate=vtrate)
      write(errunit,*) '        Calc average density...'
      densav=0.d0

      do ipar=1,npart
         densav=(densav*dble(ipar-1)+dble(density(ipar)))/dble(ipar)
      enddo

      write(errunit,*) '    --> Average density :',densav
      call system_clock(count=vt1, count_rate=vtrate)
      dvt=real(vt1-vt0,8)/real(vtrate,8)
      if (verbose) write(errunit,'(A,F10.2,A)') "     --> ",dvt," seconds to compute average density"
   endif

   nactive_parts=nactive
   if (verbose) write(errunit,*) '    [dev3] active particles counted without neighbor table:', &
 &                               nactive_parts
   call print_memory()
   call system_clock(count=tttt1, count_rate=ttttrate)
   dtdtdtdt=real(tttt1-tttt0,8)/real(ttttrate,8)
   if (verbose) write(errunit,'(A,F10.2,A)') "     --> ",dtdtdtdt," seconds to compute mean density"

end subroutine compute_mean_density_and_np


!=======================================================================
subroutine find_local_maxima
!=======================================================================
   implicit none

   integer(kind=4)             :: ipar,idist,iparid,iparsel,igroup,nmembmax,nmembtot,nequal_density
   integer(kind=4),allocatable :: nmemb(:)
   real(kind=8), dimension(0:nvoisins) :: dist2
   integer, dimension(nvoisins)        :: iparnei
   real(kind=8)                :: densest
   integer(kind=4) :: groupid
   integer(kind=4) :: tttt0, tttt1, ttttrate
   real(kind=8)    :: dtdtdtdt

   call system_clock(count=tttt0, count_rate=ttttrate)
   if (verbose) write(errunit,*) '    Now Find local maxima...'

   liste_parts(1:nusedpart)=0
   allocate(idpart_adapt(npart))
   memory_used = memory_used + real(size(idpart_adapt), kind=8)*4/1.d9
   call print_memory()

   idpart_adapt=0
   ngroups=0
   nequal_density=0
   call omp_set_num_threads(nbPes)
   if(verbose) write(errunit,*) "    [OMP] find local maxima with ncore=",nbPes
   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipar,dist2,iparnei,densest,iparsel,idist,iparid) &
   !$OMP REDUCTION(+:nequal_density)
   !$OMP DO SCHEDULE(GUIDED,256)
   do ipar=1,npart
      if (density(ipar).le.rho_threshold) cycle
      call omp_find_nearest_parts(ipar,dist2,iparnei)
      densest=density(ipar)
      iparsel=ipar
      do idist=1,nhop
         iparid=iparnei(idist)
         if (density(iparid).gt.densest) then
            iparsel=iparid
            densest=density(iparid)
         elseif (density(iparid).eq.densest) then
            iparsel=min(iparsel,iparid)
            nequal_density=nequal_density+1
         endif
      enddo
      idpart_adapt(ipar)=iparsel
   enddo
   !$OMP END DO
   !$OMP END PARALLEL

   do ipar=1,npart
      if (density(ipar).le.rho_threshold) cycle
      if (idpart_adapt(ipar).eq.ipar) then
         ngroups=ngroups+1
         idpart_adapt(ipar)=-ngroups
      endif
   enddo
   if (verbose .and. nequal_density.gt.0) then
      write(errunit,*) '    WARNING : equal densities in find_local_maxima:', nequal_density
   endif
   
   if (verbose) write(errunit,*) '    --> Number of local maxima found :',ngroups

   ! Now Link the particles associated to the same maximum
   if (verbose) write(errunit,*) '    Create link list...'

   allocate(densityg(ngroups))
   memory_used = memory_used + real(size(densityg), kind=8)*8/1.d9
   allocate(firstpart(ngroups))
   memory_used = memory_used + real(size(firstpart), kind=8)*4/1.d9
   call print_memory()
   
   do ipar=1,npart
      if (density(ipar).le.rho_threshold) cycle
      iparid=idpart_adapt(ipar)

      if (iparid.lt.0) then
         groupid=-iparid
         densityg(groupid)=density(ipar)
         liste_parts(ipar)=groupid
      else
         groupid=0

         do while (iparid.gt.0)
            if (liste_parts(iparid).gt.0) then
               groupid=liste_parts(iparid)
               exit
            endif

            iparid=idpart_adapt(iparid)
         enddo

         if (iparid.lt.0) then
            groupid=-iparid
         endif

         liste_parts(ipar)=groupid
      endif
   enddo

   if (verbose) then
      allocate(nmemb(ngroups))
      memory_used = memory_used + real(size(nmemb), kind=8)*4/1.d9
      nmemb(1:ngroups)=0
      firstpart(1:ngroups)=0
      do ipar=1,npart
         if (density(ipar).le.rho_threshold) cycle
         igroup=liste_parts(ipar)
         if (igroup.gt.0) then
            idpart_adapt(ipar)=firstpart(igroup)   
            firstpart(igroup)=ipar
            nmemb(igroup)=nmemb(igroup)+1
         endif
      enddo

   
      nmembmax=0
      nmembtot=0
      do igroup=1,ngroups
         nmembmax=max(nmembmax,nmemb(igroup))
         nmembtot=nmembtot+nmemb(igroup)
      enddo
      write(errunit,*) '    --> Number of particles of the largest group :',nmembmax
      write(errunit,*) '    --> Total number of particles in groups ',nmembtot
      deallocate(nmemb)
   endif

   nactive_parts=0
   
   call system_clock(count=tttt1, count_rate=ttttrate)
   dtdtdtdt=real(tttt1-tttt0,8)/real(ttttrate,8)
   if (verbose) write(errunit,'(A,F10.2,A)') "     --> ",dtdtdtdt," seconds to find local maxima"
  
end subroutine find_local_maxima


!=======================================================================
subroutine create_group_tree
!=======================================================================
   implicit none

   integer(kind=4) :: mass_loc,igroupref
   integer(kind=4) :: igroup,igr1,igr2,inode
   integer(kind=4) :: ttt0, ttt1, tttrate
   real(kind=8)    :: dtdtdt
   real(kind=8)    :: rhot,truemass

   if (verbose) write(errunit,*) '    Create the tree of structures of structures'

   ! End of the branches of the tree
   call system_clock(count=ttt0, count_rate=tttrate)
   call dev_compute_saddle_list
   call system_clock(count=ttt1, count_rate=tttrate)
   dtdtdt=real(ttt1-ttt0,8)/real(tttrate,8)
   if (verbose) write(errunit,'(A,F10.2,A)') "     --> ",dtdtdt," seconds to saddle_list"

   if (verbose) write(errunit,*) '    Build the hierarchical tree'

   nnodesmax=2*ngroups
   ! Allocations
   allocate(node(0:nnodesmax))
   memory_used = memory_used + real(size(node)*((4*6)+(8*8)), kind=8)/1.d9
   allocate(idgroup(ngroups))
   memory_used = memory_used + real(size(idgroup), kind=8)*4/1.d9
   allocate(color(ngroups))
   memory_used = memory_used + real(size(color), kind=8)*4/1.d9
   allocate(igroupid(ngroups))
   memory_used = memory_used + real(size(igroupid), kind=8)*4/1.d9
   allocate(idgroup_tmp(ngroups))
   memory_used = memory_used + real(size(idgroup_tmp), kind=8)*4/1.d9

   ! Initializations
   liste_parts = 0

   ! Iterative loop to build the rest of the tree
   inode              = 0
   nnodes             = 0
   rhot               = rho_threshold
   node(inode)%mother = 0
   mass_loc           = 0
   truemass           = 0.d0
   igroupref          = 0
   do igroup=1,ngroups
      call treat_particles_linkonly(igroup,rhot)
   enddo

   node(inode)%mass          = mass_loc
   node(inode)%truemass      = truemass
   node(inode)%radius        = 0.d0
   node(inode)%density       = 0.d0
   node(inode)%position(1:3) = 0.d0
   node(inode)%densmax       = maxval(densityg(1:ngroups))
   node(inode)%rho_saddle    = 0.
   node(inode)%level         = 0
   node(inode)%nsisters      = 0
   node(inode)%sister        = 0
   igr1 = 1
   igr2 = ngroups
   do igroup=1,ngroups
      idgroup(igroup)=igroup
      igroupid(igroup)=igroup
   enddo
   allocate(queue_color(ngroups))
   memory_used = memory_used + real(size(queue_color), kind=8)*4/1.d9
   if (verbose) write(errunit,*) '    [OMP] Create nodes ...'
   call system_clock(count=ttt0, count_rate=tttrate)
   call dev_create_nodes_omp_root(rhot,inode,igr1,igr2)
   call system_clock(count=ttt1, count_rate=tttrate)
   dtdtdt=real(ttt1-ttt0,8)/real(tttrate,8)
   if (verbose) write(errunit,'(A,F10.2,A)') "     --> ",dtdtdt," seconds to create nodes"
   memory_used = memory_used - real(size(queue_color), kind=8)*4/1.d9
   deallocate(queue_color)

   memory_used = memory_used - real(size(idgroup), kind=8)*4/1.d9
   memory_used = memory_used - real(size(color), kind=8)*4/1.d9
   memory_used = memory_used - real(size(igroupid), kind=8)*4/1.d9
   memory_used = memory_used - real(size(idgroup_tmp), kind=8)*4/1.d9
   memory_used = memory_used - real(size(idpart_adapt), kind=8)*4/1.d9
   memory_used = memory_used - real(size(densityg), kind=8)*8/1.d9
   memory_used = memory_used - real(size(firstpart), kind=8)*4/1.d9
   deallocate(idgroup)
   deallocate(color)
   deallocate(igroupid)
   deallocate(idgroup_tmp)
   deallocate(idpart_adapt)
   deallocate(densityg)
   deallocate(firstpart)

   call reorder_nodes_deterministic
   call node2table
 
end subroutine create_group_tree


!=======================================================================
subroutine reorder_nodes_deterministic
!=======================================================================
   implicit none

   integer(kind=4) :: inode, new_id, old_id
   integer(kind=4) :: mapped
   integer(kind=4), allocatable :: new_to_old(:)
   integer(kind=4), allocatable :: old_to_new(:)
   type(supernode), allocatable :: node_sorted(:)

   if (nnodes.le.1) return

   allocate(new_to_old(1:nnodes))
   allocate(old_to_new(0:nnodes))
   allocate(node_sorted(0:nnodes))

   old_to_new(0)=0
   do inode=1,nnodes
      new_to_old(inode)=inode
   enddo

   call sort_node_ids(new_to_old,1,nnodes)

   do new_id=1,nnodes
      old_id=new_to_old(new_id)
      old_to_new(old_id)=new_id
   enddo

   node_sorted(0)=node(0)
   do new_id=1,nnodes
      old_id=new_to_old(new_id)
      node_sorted(new_id)=node(old_id)

      mapped=node_sorted(new_id)%mother
      if (mapped.gt.0) node_sorted(new_id)%mother=old_to_new(mapped)

      mapped=node_sorted(new_id)%firstchild
      if (mapped.gt.0) node_sorted(new_id)%firstchild=old_to_new(mapped)

      mapped=node_sorted(new_id)%sister
      if (mapped.gt.0) node_sorted(new_id)%sister=old_to_new(mapped)
   enddo

   do inode=1,nusedpart
      old_id=liste_parts(inode)
      if (old_id.gt.0) liste_parts(inode)=old_to_new(old_id)
   enddo

   node(0:nnodes)=node_sorted(0:nnodes)

   deallocate(new_to_old)
   deallocate(old_to_new)
   deallocate(node_sorted)

end subroutine reorder_nodes_deterministic

!=======================================================================
recursive subroutine sort_node_ids(ids,left,right)
!=======================================================================
   implicit none

   integer(kind=4), intent(inout) :: ids(:)
   integer(kind=4), intent(in)    :: left, right
   integer(kind=4)                :: i, j, pivot, tmp

   if (left.ge.right) return

   i=left
   j=right
   pivot=ids((left+right)/2)

   do
      do while (node_precedes(ids(i),pivot))
         i=i+1
      enddo
      do while (node_precedes(pivot,ids(j)))
         j=j-1
      enddo

      if (i.le.j) then
         tmp=ids(i)
         ids(i)=ids(j)
         ids(j)=tmp
         i=i+1
         j=j-1
      endif

      if (i.gt.j) exit
   enddo

   if (left.lt.j) call sort_node_ids(ids,left,j)
   if (i.lt.right) call sort_node_ids(ids,i,right)

end subroutine sort_node_ids

!=======================================================================
logical function node_precedes(ia,ib)
!=======================================================================
   implicit none

   integer(kind=4), intent(in) :: ia, ib
   integer(kind=4)             :: k

   node_precedes=.false.

   if (node(ia)%level.ne.node(ib)%level) then
      node_precedes = node(ia)%level.lt.node(ib)%level
      return
   endif

   if (node(ia)%truemass.ne.node(ib)%truemass) then
      node_precedes = node(ia)%truemass.gt.node(ib)%truemass
      return
   endif

   if (node(ia)%mass.ne.node(ib)%mass) then
      node_precedes = node(ia)%mass.gt.node(ib)%mass
      return
   endif

   if (node(ia)%densmax.ne.node(ib)%densmax) then
      node_precedes = node(ia)%densmax.gt.node(ib)%densmax
      return
   endif

   if (node(ia)%density.ne.node(ib)%density) then
      node_precedes = node(ia)%density.gt.node(ib)%density
      return
   endif

   if (node(ia)%rho_saddle.ne.node(ib)%rho_saddle) then
      node_precedes = node(ia)%rho_saddle.gt.node(ib)%rho_saddle
      return
   endif

   do k=1,3
      if (node(ia)%position(k).ne.node(ib)%position(k)) then
         node_precedes = node(ia)%position(k).lt.node(ib)%position(k)
         return
      endif
   enddo

   node_precedes = ia.lt.ib

end function node_precedes


!=======================================================================
subroutine dev_create_nodes_omp_root(rhot,inode,igr1,igr2)
!=======================================================================
   implicit none
   real(kind=8)    :: rhot
   integer(kind=4) :: inode, igr1, igr2
   !$OMP PARALLEL NUM_THREADS(nbPes)
   !$OMP SINGLE
   call omp_dev_create_nodes_task(rhot,inode,igr1,igr2,.false.)
   !$OMP END SINGLE
   !$OMP END PARALLEL

end subroutine dev_create_nodes_omp_root

!=======================================================================
subroutine reserve_nodes_block(nnew,inode_first)
!=======================================================================
   implicit none
   integer(kind=4), intent(in)  :: nnew
   integer(kind=4), intent(out) :: inode_first

   !$OMP CRITICAL(node_reserve)
   inode_first = nnodes + 1
   nnodes = nnodes + nnew
   if (nnodes.gt.nnodesmax) then
      write(errunit,*) 'ERROR in reserve_nodes_block :'
      write(errunit,*) 'nnodes > nnodes max'
      STOP
   endif
   !$OMP END CRITICAL(node_reserve)

end subroutine reserve_nodes_block

!=======================================================================
subroutine dev_do_colorize_local(igroup, igr, rhot, queue, color_id)
!=======================================================================
   implicit none

   integer(kind=4) :: igroup, igr, igroup_now, color_id
   integer(kind=4) :: ineig, igroup2, igr2, neig
   integer(kind=4) :: front, back
   integer(kind=4) :: base, idx
   integer(kind=4) :: queue(:)
   real(kind=8)    :: rhot

   queue(1) = igroup
   front = 1
   back = 2
   color(igr) = color_id

   do while (front < back)
      igroup_now = queue(front)
      front = front + 1
      neig = dev_group_nhnei(igroup_now)
      base = dev_group_offset(igroup_now)

      do ineig = 1, neig
         idx = base + ineig - 1

         if (dev_rho_saddle_all(idx) > rhot) then
            igroup2 = dev_isad_all(idx)
            igr2 = igroupid(igroup2)

            if (color(igr2) == 0) then
               color(igr2) = color_id

               if (back.gt.size(queue)) then
                  write(errunit,*) 'ERROR in dev_do_colorize_local: queue overflow'
                  write(errunit,*) 'size(queue), back =',size(queue),back
                  STOP
               endif

               queue(back) = igroup2
               back = back + 1

            elseif (color(igr2) /= color_id) then
               write(errunit,*) 'ERROR in dev_do_colorize_local : color(igr2) <> color_id'
               write(errunit,*) 'The connections are not symmetric or the task range is not isolated.'
               STOP
            endif

         else
            dev_group_nhnei(igroup_now) = ineig - 1
            exit
         endif
      enddo
   enddo

end subroutine dev_do_colorize_local

!=======================================================================
recursive subroutine omp_dev_create_nodes_task(rhot0,inode0,igr1_0,igr2_0,do_paint)
!=======================================================================
   implicit none

   real(kind=8),    intent(in) :: rhot0
   integer(kind=4), intent(in) :: inode0, igr1_0, igr2_0
   logical,         intent(in) :: do_paint

   integer(kind=4)              :: inode
   real(kind=8)                 :: rhot
   integer(kind=4)              :: icolor,igr_eff,ncolors
   integer(kind=4)              :: igroup, igr,igr1,igr2
   integer(kind=4)              :: inc_color_tot,inode1,mass_loc,masstmp,igroupref
   integer(kind=4)              :: inodeout,igr1out,igr2out,isisters,nsisters
   integer(kind=4)              :: mass_comp,icolor_ref
   real(kind=8)                 :: posg(3),posgtmp(3),posref(3),rhotout,rsquaretmp,rsquareg
   real(kind=8)                 :: densmoyg,densmoytmp
   real(kind=8)                 :: densmoy_comp_max,truemass,truemasstmp
   real(kind=8)                 :: posfin(3)
   real(kind=8)                 :: densmaxgroup

   integer(kind=4), allocatable :: igrpos(:),igrinc(:)
   integer(kind=4), allocatable :: igrposnew(:)
   integer(kind=4), allocatable :: massg(:)
   real(kind=8),    allocatable :: truemassg(:)
   real(kind=8),    allocatable :: densmaxg(:)
   integer(kind=4), allocatable :: mass_compg(:)
   real(kind=8),    allocatable :: posgg(:,:)
   real(kind=8),    allocatable :: rsquare(:)
   real(kind=8),    allocatable :: densmoy(:)
   logical,         allocatable :: ifok_tmp(:)
   integer(kind=4), allocatable :: ok_color(:)
   integer(kind=4), allocatable :: queue_local(:)

   logical                      :: ifok_this
   integer(kind=4), parameter   :: task_group_min = 32

   integer(kind=4) :: t_inodeout, t_igr1out, t_igr2out
   real(kind=8)    :: t_rhotout
   logical         :: make_task

   integer(kind=4) :: ipar

!  Initialize current job
   rhot  = rhot0
   inode = inode0
   igr1  = igr1_0
   igr2  = igr2_0

!  Paint current node if this is a child task.
!  This replaces parent-side paint_particles2.
   if (do_paint) then
      do igr=igr1,igr2
         igroup=idgroup(igr)
         ipar=firstpart(igroup)
         do while (ipar.gt.0)
            liste_parts(ipar)=inode
            ipar=idpart_adapt(ipar)
         enddo
      enddo
   endif
!  nsisters == 1 chain is handled inside this loop without recursive call.
   do

      allocate(queue_local(max(igr2-igr1+1,1)))

!     ------------------------------------------------------------
!     Percolate the groups
!     ------------------------------------------------------------
      color(igr1:igr2)=0

      ncolors=0
      do igr=igr1,igr2
         igroup=idgroup(igr)
         if (color(igr).eq.0) then
            ncolors=ncolors+1
            call dev_do_colorize_local(igroup,igr,rhot,queue_local,ncolors)
         endif
      enddo

      deallocate(queue_local)

!     ------------------------------------------------------------
!     Select groups above rhot and compact by color
!     ------------------------------------------------------------
      allocate(igrpos(0:ncolors))
      allocate(igrinc(1:ncolors))

      igrpos(0)=igr1-1
      igrpos(1:ncolors)=0

      do igr=igr1,igr2
         icolor=color(igr)
         igroup=idgroup(igr)
         if (densityg(igroup).gt.rhot) &
&           igrpos(icolor)=igrpos(icolor)+1
      enddo

      do icolor=1,ncolors
         igrpos(icolor)=igrpos(icolor-1)+igrpos(icolor)
      enddo

      if (igrpos(ncolors)-igr1+1.eq.0) then
         write(errunit,*) 'ERROR in omp_dev_create_nodes_task :'
         write(errunit,*) 'All subgroups are below the threshold.'
      endif

      igrinc(1:ncolors)=0

      do igr=igr1,igr2
         icolor=color(igr)
         igroup=idgroup(igr)
         if (densityg(igroup).gt.rhot) then
            igrinc(icolor)=igrinc(icolor)+1
            igr_eff=igrinc(icolor)+igrpos(icolor-1)
            idgroup_tmp(igr_eff)=igroup
            igroupid(igroup)=igr_eff
         endif
      enddo

      igr2=igrpos(ncolors)
      idgroup(igr1:igr2)=idgroup_tmp(igr1:igr2)

      inc_color_tot=0
      do icolor=1,ncolors
         if (igrinc(icolor).gt.0) inc_color_tot=inc_color_tot+1
      enddo

      allocate(igrposnew(0:inc_color_tot))
      igrposnew(0)=igrpos(0)

      inc_color_tot=0
      do icolor=1,ncolors
         if (igrinc(icolor).gt.0) then
            inc_color_tot=inc_color_tot+1
            igrposnew(inc_color_tot)=igrpos(icolor)
         endif
      enddo

      deallocate(igrpos)
      deallocate(igrinc)

!     ------------------------------------------------------------
!     Compute properties of each color component
!     ------------------------------------------------------------
      isisters=0

      allocate(posgg(1:3,1:inc_color_tot))
      allocate(massg(1:inc_color_tot))
      allocate(truemassg(1:inc_color_tot))
      allocate(densmaxg(1:inc_color_tot))
      allocate(mass_compg(1:inc_color_tot))
      allocate(rsquare(1:inc_color_tot))
      allocate(densmoy(1:inc_color_tot))
      allocate(ok_color(1:inc_color_tot))
      allocate(ifok_tmp(1:inc_color_tot))

      ifok_tmp(1:inc_color_tot)=.false.

!     Keep this loop serial inside each task.
!     Task-level parallelism is handled by child subtrees.
      do icolor=1,inc_color_tot
         posg(1:3)=0.d0
         mass_loc=0
         truemass=0.d0
         rsquareg=0.d0
         densmoyg=0.d0

         igr1out=igrposnew(icolor-1)+1
         igr2out=igrposnew(icolor)

         densmaxgroup=-1.d0
         mass_comp=0
         densmoy_comp_max=-1.d0
         igroupref=0

         do igr=igr1out,igr2out
            igroup=idgroup(igr)
            densmaxgroup=max(densmaxgroup,densityg(igroup))
            call treat_particles(igroup,rhot,posgtmp,masstmp, &
&                                 igroupref,posref,rsquaretmp, &
&                                 densmoytmp,truemasstmp)

            posg(1)=posg(1)+posgtmp(1)
            posg(2)=posg(2)+posgtmp(2)
            posg(3)=posg(3)+posgtmp(3)

            rsquareg=rsquareg+rsquaretmp
            mass_loc=mass_loc+masstmp
            truemass=truemass+truemasstmp
            densmoyg=densmoyg+densmoytmp
            mass_comp=max(mass_comp,masstmp)

            if (masstmp > 0) then
               densmoytmp=densmoytmp/dble(masstmp)
               densmoy_comp_max=max(densmoy_comp_max, &
&                 densmoytmp/(1.d0+fudge/sqrt(dble(masstmp))))
            endif
         enddo

         massg(icolor)=mass_loc
         truemassg(icolor)=truemass
         posgg(1:3,icolor)=posg(1:3)
         densmaxg(icolor)=densmaxgroup
         mass_compg(icolor)=mass_comp

         rsquare(icolor)=sqrt(abs( &
&           (truemass*rsquareg- &
&           (posg(1)**2+posg(2)**2+posg(3)**2) )/ &
&            truemass**2 ))

         densmoy(icolor)=densmoyg/dble(mass_loc)

         ifok_this=mass_loc.ge.nmembthresh.and. &
&           (densmoy(icolor).gt.rhot*(1.d0+fudge/sqrt(dble(mass_loc))).or. &
&            densmoy_comp_max.gt.rhot).and. &
&           densmaxg(icolor).ge.alphap*densmoy(icolor).and. &
&           rsquare(icolor).ge.epsilon

         ifok_tmp(icolor)=ifok_this
      enddo

      isisters=0
      icolor_ref=0

      do icolor=1,inc_color_tot
         if (ifok_tmp(icolor)) then
            isisters=isisters+1
            ok_color(isisters)=icolor
            icolor_ref=icolor
         endif
      enddo

      nsisters=isisters

!     ------------------------------------------------------------
!     Case 1: branching node, create child nodes and spawn tasks
!     ------------------------------------------------------------
      if (nsisters.gt.1) then

         call reserve_nodes_block(nsisters,inode1)

!        Create child node metadata in deterministic sibling order
         do isisters=1,nsisters
            icolor=ok_color(isisters)
            inodeout=inode1+isisters-1

            node(inodeout)%mother=inode
            node(inodeout)%densmax=densmaxg(icolor)

            if (isisters.gt.1) then
               node(inodeout)%sister=inodeout-1
            else
               node(inodeout)%sister=0
            endif

            node(inodeout)%nsisters=nsisters
            node(inodeout)%mass=massg(icolor)
            node(inodeout)%truemass=truemassg(icolor)

            if (massg(icolor).eq.0) then
               write(errunit,*) 'ERROR in omp_dev_create_nodes_task :'
               write(errunit,*) 'NULL mass for inode=',inodeout
               STOP
            endif

            posfin(1:3)=real(posgg(1:3,icolor)/truemassg(icolor),8)
            node(inodeout)%radius=real(rsquare(icolor),8)
            node(inodeout)%density=real(densmoy(icolor),8)

            if (posfin(1).ge.xlongs2) then
               posfin(1)=posfin(1)-xlong
            elseif (posfin(1).lt.-xlongs2) then
               posfin(1)=posfin(1)+xlong
            endif

            if (posfin(2).ge.ylongs2) then
               posfin(2)=posfin(2)-ylong
            elseif (posfin(2).lt.-ylongs2) then
               posfin(2)=posfin(2)+ylong
            endif

            if (posfin(3).ge.zlongs2) then
               posfin(3)=posfin(3)-zlong
            elseif (posfin(3).lt.-zlongs2) then
               posfin(3)=posfin(3)+zlong
            endif

            node(inodeout)%position(1:3)=posfin(1:3)
            node(inodeout)%rho_saddle=rhot
            node(inodeout)%level=node(inode)%level+1
         enddo

!        Original structure uses last child as firstchild head
         node(inode)%firstchild=inode1+nsisters-1

!        Spawn child subtree tasks
         do isisters=1,nsisters
            icolor=ok_color(isisters)

            t_inodeout=inode1+isisters-1
            t_igr1out=igrposnew(icolor-1)+1
            t_igr2out=igrposnew(icolor)
            t_rhotout=rhot*(1.d0+fudge/sqrt(dble(mass_compg(icolor))))

            make_task = (t_igr2out-t_igr1out+1 .ge. task_group_min)

            !$OMP TASK FIRSTPRIVATE(t_rhotout,t_inodeout,t_igr1out,t_igr2out) IF(make_task)
            call omp_dev_create_nodes_task(t_rhotout,t_inodeout,t_igr1out,t_igr2out,.true.)
            !$OMP END TASK
         enddo

         !$OMP TASKWAIT

         deallocate(igrposnew)
         deallocate(posgg)
         deallocate(massg)
         deallocate(truemassg)
         deallocate(densmaxg)
         deallocate(densmoy)
         deallocate(mass_compg)
         deallocate(rsquare)
         deallocate(ok_color)
         deallocate(ifok_tmp)

         exit

!     ------------------------------------------------------------
!     Case 2: only one surviving component.
!     Continue in the same task and same inode.
!     ------------------------------------------------------------
      elseif (nsisters.eq.1) then

         icolor=icolor_ref
         rhotout=rhot*(1.d0+fudge/sqrt(dble(mass_compg(icolor))))
         igr1out=igrposnew(0)+1
         igr2out=igrposnew(inc_color_tot)

         deallocate(igrposnew)
         deallocate(posgg)
         deallocate(massg)
         deallocate(truemassg)
         deallocate(densmaxg)
         deallocate(densmoy)
         deallocate(mass_compg)
         deallocate(rsquare)
         deallocate(ok_color)
         deallocate(ifok_tmp)

         if (igr2out.ne.igr1out) then
            rhot=rhotout
            igr1=igr1out
            igr2=igr2out
            cycle
         else
            node(inode)%firstchild=0
            exit
         endif

!     ------------------------------------------------------------
!     Case 3: no surviving component
!     ------------------------------------------------------------
      else

         node(inode)%firstchild=0

         deallocate(igrposnew)
         deallocate(posgg)
         deallocate(massg)
         deallocate(truemassg)
         deallocate(densmaxg)
         deallocate(densmoy)
         deallocate(mass_compg)
         deallocate(rsquare)
         deallocate(ok_color)
         deallocate(ifok_tmp)

         exit
      endif

   enddo

end subroutine omp_dev_create_nodes_task


!=======================================================================
subroutine treat_particles(igroup,rhot,posg,imass,igroupref,posref, &
&                          rsquare,densmoy,truemass)
!=======================================================================
   implicit none

   real(kind=8)    :: rhot, density_ipar
   real(kind=8)    :: posg(3),posref(3)
   real(kind=8)    :: posdiffx,posdiffy,posdiffz,rsquare,densmoy,truemass,xmasspart
   real(kind=8)    :: densmax,densmin
   integer(kind=4) :: imass,ipar,iparold,igroupref
   integer(kind=4) :: igroup
   logical         :: first_good

   imass=0
   truemass=0.d0
   rsquare=0.d0
   densmoy=0.d0
   posg(1:3)=0.d0
   ipar=firstpart(igroup)
   first_good=.false.
   do while (ipar.gt.0)
      density_ipar = density(ipar)
      if (density_ipar.gt.rhot) then
         if (.not.first_good) then
            if (igroupref.eq.0) then
               posref(1:3)=pos(1:3,ipar)
               igroupref=igroup
            endif
            first_good=.true.
            firstpart(igroup)=ipar
            densmin=density_ipar
            densmax=densmin
         else
            idpart_adapt(iparold)=ipar
         endif
         iparold=ipar
         imass=imass+1
         xmasspart=mass(ipar)
         truemass=truemass+xmasspart
         posdiffx=pos(1,ipar)-posref(1)
         posdiffy=pos(2,ipar)-posref(2)
         posdiffz=pos(3,ipar)-posref(3)
         if (posdiffx.ge.xlongs2) then
            posdiffx=posdiffx-xlong
         elseif (posdiffx.lt.-xlongs2) then
            posdiffx=posdiffx+xlong
         endif
         if (posdiffy.ge.ylongs2) then
            posdiffy=posdiffy-ylong
         elseif (posdiffy.lt.-ylongs2) then
            posdiffy=posdiffy+ylong
         endif
         if (posdiffz.ge.zlongs2) then
            posdiffz=posdiffz-zlong
         elseif (posdiffz.lt.-zlongs2) then
            posdiffz=posdiffz+zlong
         endif
         posdiffx=posdiffx+posref(1)
         posdiffy=posdiffy+posref(2)
         posdiffz=posdiffz+posref(3)
         posg(1)=posg(1)+posdiffx*xmasspart
         posg(2)=posg(2)+posdiffy*xmasspart
         posg(3)=posg(3)+posdiffz*xmasspart
         rsquare=rsquare+xmasspart*(posdiffx*posdiffx+posdiffy*posdiffy+posdiffz*posdiffz)
         densmoy=densmoy+dble(density_ipar)
         densmax=max(densmax,density_ipar)
         densmin=min(densmin,density_ipar)
      endif
      ipar=idpart_adapt(ipar)
   enddo
   ! if (.not.first_good) firstpart(igroup)=0
   if (first_good) then
      idpart_adapt(iparold)=0
   else
      firstpart(igroup)=0
   endif

   if (densmin.le.rhot.or.densmax.ne.densityg(igroup)) then
      write(errunit,*) 'ERROR in treat_particles'
      write(errunit,*) 'igroup, densmax, rhot=',igroup,densityg(igroup),rhot
      write(errunit,*) 'denslow, denshigh    =',densmin,densmax
      STOP
   endif

end subroutine treat_particles

!=======================================================================
subroutine treat_particles_linkonly(igroup,rhot)
!=======================================================================
   implicit none

   real(kind=8)    :: rhot
   integer(kind=4) :: ipar,iparold
   integer(kind=4) :: igroup
   logical         :: first_good

   ipar=firstpart(igroup)
   first_good=.false.
   do while (ipar.gt.0)
      if (density(ipar).gt.rhot) then
         if (.not.first_good) then
            first_good=.true.
            firstpart(igroup)=ipar
         else
            idpart_adapt(iparold)=ipar
         endif
         iparold=ipar
      endif
      ipar=idpart_adapt(ipar)
   enddo
   if (.not.first_good) firstpart(igroup)=0

end subroutine treat_particles_linkonly


!=======================================================================
subroutine dev_compute_saddle_list
!=======================================================================
! CSR version of compute_saddle_list.
! Instead of group(igroup)%isad_gr(:), use:
!   dev_group_nhnei(igroup)
!   dev_group_offset(igroup)
!   dev_isad_all(dev_group_offset(igroup):...)
!   dev_rho_saddle_all(dev_group_offset(igroup):...)
!=======================================================================

   implicit none

   integer(kind=4) :: ipar1, ipar2, igroup1, igroup2
   integer(kind=4) :: idist, ineig, ineig2, neig
   integer(kind=4) :: in1, in2, i, idx1, idx2
   integer(kind=4) :: idestroy, icon_count
   integer(kind=4) :: max_neig
   integer(kind=4) :: total_conn, new_total_conn
   integer(kind=4) :: base, base2, outbase
   integer(kind=4) :: tid
   real(kind=8), dimension(0:nvoisins) :: dist2
   integer, dimension(nvoisins)        :: iparnei
   real(kind=8)    :: density1, density2, rho_sad12
   logical         :: exist

   logical, allocatable         :: touch_thr(:,:)
   integer(kind=4), allocatable :: listg_thr(:,:)
   real(kind=8), allocatable    :: rho_sad(:)
   real(kind=8), allocatable    :: rho_sad_thr(:,:)
   integer(kind=4), allocatable :: isad(:)
   integer(kind=4), allocatable :: indx(:)

   integer(kind=4), allocatable :: new_group_nhnei(:)
   integer(kind=4), allocatable :: new_group_offset(:)
   integer(kind=4), allocatable :: new_isad_all(:)
   real(kind=8), allocatable    :: new_rho_saddle_all(:)

   if (verbose) write(errunit,*) '    Fill the end of the branches of the group tree'

   if (verbose) write(errunit,*) '    First count the number of neighbourgs'// &
&                                ' of each elementary group...'

   allocate(dev_group_nhnei(ngroups))
   memory_used = memory_used + real(size(dev_group_nhnei), kind=8)*4/1.d9
   allocate(dev_group_offset(ngroups+1))
   memory_used = memory_used + real(size(dev_group_offset), kind=8)*4/1.d9
   allocate(touch_thr(ngroups,nbPes))
   memory_used = memory_used + real(size(touch_thr), kind=8)*4/1.d9
   allocate(listg_thr(ngroups,nbPes))
   memory_used = memory_used + real(size(listg_thr), kind=8)*4/1.d9

   touch_thr(1:ngroups,1:nbPes)=.false.
   dev_group_nhnei(1:ngroups)=0

!=======================================================================
! First pass: count unique neighbouring groups for each group
!=======================================================================

   max_neig = 0

   call omp_set_num_threads(nbPes)
   if (verbose) write(errunit,*) '    [OMP] Count saddle neighbours with ncore=',nbPes
   !$OMP PARALLEL DEFAULT(SHARED) &
   !$OMP PRIVATE(tid,igroup1,ineig,ipar1,dist2,iparnei,idist,ipar2,igroup2,in1) &
   !$OMP REDUCTION(MAX:max_neig)
   tid=omp_get_thread_num()+1
   !$OMP DO SCHEDULE(DYNAMIC,8)
   do igroup1=1,ngroups
      ineig=0
      ipar1=firstpart(igroup1)

      do while (ipar1.gt.0)
         call omp_find_nearest_parts(ipar1,dist2,iparnei)
         do idist=1,nhop
            ipar2=iparnei(idist)
            igroup2=liste_parts(ipar2)

            if (igroup2.gt.0 .and. igroup2.ne.igroup1) then
               if (.not.touch_thr(igroup2,tid)) then
                  ineig=ineig+1
                  touch_thr(igroup2,tid)=.true.
                  listg_thr(ineig,tid)=igroup2
               endif
            endif
         enddo
         ipar1=idpart_adapt(ipar1)
      enddo

!     Reinitialize touch
      do in1=1,ineig
         igroup2=listg_thr(in1,tid)
         touch_thr(igroup2,tid)=.false.
      enddo

      dev_group_nhnei(igroup1)=ineig
      max_neig=max(max_neig,ineig)
   enddo
   !$OMP END DO
   !$OMP END PARALLEL

!=======================================================================
! Build CSR offset
!=======================================================================

   dev_group_offset(1)=1
   do igroup1=1,ngroups
      dev_group_offset(igroup1+1)=dev_group_offset(igroup1)+dev_group_nhnei(igroup1)
   enddo

   total_conn=dev_group_offset(ngroups+1)-1

   allocate(dev_isad_all(max(total_conn,1)))
   memory_used = memory_used + real(size(dev_isad_all), kind=8)*4/1.d9
   allocate(dev_rho_saddle_all(max(total_conn,1)))
   memory_used = memory_used + real(size(dev_rho_saddle_all), kind=8)*8/1.d9

   if (max_neig.gt.0) then
      allocate(rho_sad(max_neig))
      allocate(isad(max_neig))
      allocate(indx(max_neig))
      allocate(rho_sad_thr(max_neig,nbPes))
      memory_used = memory_used + real(size(rho_sad_thr), kind=8)*8/1.d9
   else
      allocate(rho_sad(1))
      allocate(isad(1))
      allocate(indx(1))
      allocate(rho_sad_thr(1,nbPes))
      memory_used = memory_used + real(size(rho_sad_thr), kind=8)*8/1.d9
   endif
   memory_used = memory_used + real(size(rho_sad), kind=8)*8/1.d9
   memory_used = memory_used + real(size(isad), kind=8)*4/1.d9
   memory_used = memory_used + real(size(indx), kind=8)*4/1.d9

   if (verbose) write(errunit,*) '    Compute lists of neighbourgs and saddle points...'

!=======================================================================
! Second pass: fill neighbour list and saddle densities
!=======================================================================

   if (verbose) write(errunit,*) '    [OMP] Fill saddle neighbours with ncore=',nbPes
   !$OMP PARALLEL DEFAULT(SHARED) &
   !$OMP PRIVATE(tid,igroup1,neig,ineig,ipar1,base,dist2,iparnei,density1, &
   !$OMP         idist,ipar2,igroup2,density2,rho_sad12,ineig2,in1)
   tid=omp_get_thread_num()+1
   !$OMP DO SCHEDULE(DYNAMIC,8)
   do igroup1=1,ngroups
      neig=dev_group_nhnei(igroup1)

      if (neig.gt.0) then
         ineig=0
         ipar1=firstpart(igroup1)
         base=dev_group_offset(igroup1)

         do while (ipar1.gt.0)
            call omp_find_nearest_parts(ipar1,dist2,iparnei)
            density1=density(ipar1)

            do idist=1,nhop
               ipar2=iparnei(idist)
               igroup2=liste_parts(ipar2)

               if (igroup2.gt.0 .and. igroup2.ne.igroup1) then
                  density2=density(ipar2)
                  rho_sad12=min(density1,density2)

                  if (.not.touch_thr(igroup2,tid)) then
                     ineig=ineig+1
                     touch_thr(igroup2,tid)=.true.
                     listg_thr(igroup2,tid)=ineig

                     rho_sad_thr(ineig,tid)=rho_sad12
                     dev_isad_all(base+ineig-1)=igroup2
                  else
                     ineig2=listg_thr(igroup2,tid)
                     rho_sad_thr(ineig2,tid)=max(rho_sad_thr(ineig2,tid),rho_sad12)
                  endif
               endif
            enddo

            ipar1=idpart_adapt(ipar1)
         enddo

         if (ineig.ne.neig) then
            write(errunit,*) 'ERROR in dev_compute_saddle_list :'
            write(errunit,*) 'The number of neighbourgs does not match.'
            write(errunit,*) 'ineig, neig =',ineig,neig
            STOP
         endif

         dev_rho_saddle_all(base:base+ineig-1)=rho_sad_thr(1:ineig,tid)

!        Reinitialize touch
         do in1=1,ineig
            igroup2=dev_isad_all(base+in1-1)
            touch_thr(igroup2,tid)=.false.
         enddo
      endif
   enddo
   !$OMP END DO
   !$OMP END PARALLEL

   memory_used = memory_used - real(size(touch_thr), kind=8)*4/1.d9
   deallocate(touch_thr)
   memory_used = memory_used - real(size(listg_thr), kind=8)*4/1.d9
   deallocate(listg_thr)
   memory_used = memory_used - real(size(rho_sad_thr), kind=8)*8/1.d9
   deallocate(rho_sad_thr)

!=======================================================================
! Establish symmetry in connections
!=======================================================================

   if (verbose) write(errunit,*) '    Establish symmetry in connections...'

   icon_count=0
   idestroy=0

   do igroup1=1,ngroups
      neig=dev_group_nhnei(igroup1)

      if (neig.gt.0) then
         base=dev_group_offset(igroup1)

         do in1=1,neig
            idx1=base+in1-1
            exist=.false.
            igroup2=dev_isad_all(idx1)

            if (igroup2.gt.0) then
               base2=dev_group_offset(igroup2)

               do in2=1,dev_group_nhnei(igroup2)
                  idx2=base2+in2-1

                  if (dev_isad_all(idx2).eq.igroup1) then
                     exist=.true.

                     rho_sad12=min(dev_rho_saddle_all(idx2), &
&                                  dev_rho_saddle_all(idx1))

                     dev_rho_saddle_all(idx2)=rho_sad12
                     dev_rho_saddle_all(idx1)=rho_sad12
                     exit
                  endif
               enddo
            endif

            if (.not.exist) then
               dev_isad_all(idx1)=0
               idestroy=idestroy+1
            else
               icon_count=icon_count+1
            endif
         enddo
      endif
   enddo

   if (verbose) write(errunit,*) '    --> Number of connections removed :',idestroy
   if (verbose) write(errunit,*) '    --> Total number of connections remaining :',icon_count

!=======================================================================
! Rebuild compact CSR list and sort saddle points by decreasing order
!=======================================================================

   if (verbose) then
      write(errunit,*) '    Rebuild groups with undesired connections removed...'
   endif

!  Count remaining connections per group
   allocate(new_group_nhnei(ngroups))
   memory_used = memory_used + real(size(new_group_nhnei), kind=8)*4/1.d9
   allocate(new_group_offset(ngroups+1))
   memory_used = memory_used + real(size(new_group_offset), kind=8)*4/1.d9

   new_group_nhnei(1:ngroups)=0

   do igroup1=1,ngroups
      neig=dev_group_nhnei(igroup1)
      if (neig.gt.0) then
         base=dev_group_offset(igroup1)
         do in1=1,neig
            idx1=base+in1-1
            if (dev_isad_all(idx1).gt.0) then
               new_group_nhnei(igroup1)=new_group_nhnei(igroup1)+1
            endif
         enddo
      endif
   enddo

!  New offsets
   new_group_offset(1)=1
   do igroup1=1,ngroups
      new_group_offset(igroup1+1)=new_group_offset(igroup1)+new_group_nhnei(igroup1)
   enddo

   new_total_conn=new_group_offset(ngroups+1)-1

   allocate(new_isad_all(max(new_total_conn,1)))
   memory_used = memory_used + real(size(new_isad_all), kind=8)*4/1.d9
   allocate(new_rho_saddle_all(max(new_total_conn,1)))
   memory_used = memory_used + real(size(new_rho_saddle_all), kind=8)*8/1.d9

!  Compact and sort each group
   do igroup1=1,ngroups
      neig=dev_group_nhnei(igroup1)

      if (neig.gt.0) then
         base=dev_group_offset(igroup1)
         ineig=0

         do in1=1,neig
            idx1=base+in1-1
            igroup2=dev_isad_all(idx1)

            if (igroup2.gt.0) then
               ineig=ineig+1
               rho_sad(ineig)=dev_rho_saddle_all(idx1)
               isad(ineig)=igroup2
            endif
         enddo

         if (ineig.ne.new_group_nhnei(igroup1)) then
            write(errunit,*) 'ERROR in dev_compute_saddle_list :'
            write(errunit,*) 'Rebuild count mismatch.'
            write(errunit,*) 'igroup, ineig, new_nhnei =', &
&                             igroup1, ineig, new_group_nhnei(igroup1)
            STOP
         endif

         if (ineig.gt.0) then
            call indexx(ineig,rho_sad,indx)

            outbase=new_group_offset(igroup1)

            do i=1,ineig
               ineig2=indx(i)

!              Store in decreasing order
               new_isad_all(outbase+ineig-i)=isad(ineig2)
               new_rho_saddle_all(outbase+ineig-i)=rho_sad(ineig2)
            enddo
         endif
      endif
   enddo

!  Replace old CSR arrays with compacted CSR arrays
   memory_used = memory_used - real(size(dev_group_nhnei), kind=8)*4/1.d9
   memory_used = memory_used - real(size(dev_group_offset), kind=8)*4/1.d9
   memory_used = memory_used - real(size(dev_isad_all), kind=8)*4/1.d9
   memory_used = memory_used - real(size(dev_rho_saddle_all), kind=8)*8/1.d9
   deallocate(dev_group_nhnei)
   deallocate(dev_group_offset)
   deallocate(dev_isad_all)
   deallocate(dev_rho_saddle_all)

   call move_alloc(new_group_nhnei,dev_group_nhnei)
   call move_alloc(new_group_offset,dev_group_offset)
   call move_alloc(new_isad_all,dev_isad_all)
   call move_alloc(new_rho_saddle_all,dev_rho_saddle_all)

   memory_used = memory_used - real(size(rho_sad), kind=8)*8/1.d9
   memory_used = memory_used - real(size(isad), kind=8)*4/1.d9
   memory_used = memory_used - real(size(indx), kind=8)*4/1.d9
   deallocate(rho_sad)
   deallocate(isad)
   deallocate(indx)
   memory_used = memory_used - real(size(mass_cell), kind=8)*4/1.d9
   memory_used = memory_used - real(size(size_cell), kind=8)*8/1.d9
   memory_used = memory_used - real(size(pos_cell), kind=8)*8/1.d9
   memory_used = memory_used - real(size(sister), kind=8)*4/1.d9
   memory_used = memory_used - real(size(firstchild), kind=8)*4/1.d9
   memory_used = memory_used - real(size(idpart), kind=8)*4/1.d9
   deallocate(mass_cell)
   deallocate(size_cell)
   deallocate(pos_cell)
   deallocate(sister)
   deallocate(firstchild)
   deallocate(idpart)
   nactive_parts=0

end subroutine dev_compute_saddle_list

!=======================================================================
subroutine omp_compute_density(ipar,dist2,iparnei)
!=======================================================================
   implicit none
   real(kind=8)          :: dist2(0:nvoisins)
   integer(kind=4)       :: iparnei(nvoisins)
   real(kind=8)          :: r,unsr,contrib,r3
   ! real(kind=8),external :: spline
   integer(kind=4)       :: idist,ipar

   r=sqrt(dist2(nvoisins))*0.5d0
   r3 = r*r*r
   unsr=1.d0/r
   contrib=0.d0
   do idist=1,nvoisins-1
      contrib=contrib+mass(iparnei(idist))*spline(sqrt(dist2(idist))*unsr)
   enddo
   ! Add the contribution of the particle itself and normalize properly
   ! to get a density with average unity (if computed on an uniform grid)
   ! note that this assumes that the total mass in the box is normalized to 1.
   density(ipar)=xyzlong*(contrib+mass(ipar)) &
         &              /(pi*r3)

end subroutine omp_compute_density

!=======================================================================
subroutine omp_find_nearest_parts(ipar,dist2,iparnei)
!=======================================================================
! For ID(ipar) particle, find `nvoisins` nearest particles
! and store their distances in dist2 and their IDs in iparnei
   implicit none

   integer(kind=4) :: ipar,idist,icell_identity,inccellpart
   real(kind=8)    :: dist2(0:nvoisins)
   integer(kind=4) :: iparnei(nvoisins)
   real(kind=8)    :: poshere(1:3)

   ! poshere(1:3)=pos(ipar,1:3)
   poshere(1:3)=pos(1:3,ipar)
   dist2(0)=0.
   do idist=1,nvoisins
      dist2(idist)=bignum
   enddo
   icell_identity =1
   inccellpart    =0
   call walk_tree2(icell_identity,poshere,dist2,ipar,inccellpart,iparnei)

end subroutine omp_find_nearest_parts



!=======================================================================
recursive subroutine walk_tree2(icellidin,poshere,dist2,    &
 &                   iparid,inccellpart,iparnei)
!=======================================================================
   implicit none


   integer(kind=4) :: icellidin,icell_identity,iparid,inccellpart,ic,iparcell
   real(kind=8)    :: poshere(3),dist2(0:nvoisins)
   real(kind=8)    :: dx,dy,dz,distance2,sizec,r2max
   integer(kind=4) :: idist,inc
   integer(kind=4) :: icellid_out
   real(kind=8)    :: discell2(0:8)
   integer(kind=4) :: iparnei(nvoisins)
   integer(kind=4) :: icid(8)

   integer(kind=4) :: i,npart_pos_this_node
   real(kind=8)    :: distance2p

   icell_identity=firstchild(icellidin)
   inc=1
   r2max = dist2(nvoisins)
   discell2(0)=0
   discell2(1:8)=1.d30
   ! Until icell_identity==0: (Final leaf of tree)
   ! Calc distance (poshere <-> cells)

   do while (icell_identity.ne.0)
      sizec=size_cell(icell_identity)

      dx=abs(pos_cell(1,icell_identity)-poshere(1))
      dy=abs(pos_cell(2,icell_identity)-poshere(2))
      dz=abs(pos_cell(3,icell_identity)-poshere(3))
      dx=max(0.d0,min(dx,xlong-dx)-sizec)
      dy=max(0.d0,min(dy,ylong-dy)-sizec)
      dz=max(0.d0,min(dz,zlong-dz)-sizec)
      distance2=dx*dx+dy*dy+dz*dz

      if (distance2.lt.r2max) then
         idist=inc-1
         do while (discell2(idist).gt.distance2)
            discell2(idist+1)=discell2(idist)
            icid(idist+1)=icid(idist)
            idist=idist-1
         enddo
         discell2(idist+1)=distance2
         icid(idist+1)=icell_identity
         inc=inc+1
      endif

      icell_identity=sister(icell_identity)
   enddo

   inccellpart=inccellpart+inc-1
   ! Loop for counted cells,
   ! Update the closest particle ID(in iparnei) and the distance to that part(in dist2)
   do ic=1,inc-1
      icellid_out=icid(ic)
      if (firstchild(icellid_out) < 0) then ! i.e., if leaf cell
         if (discell2(ic).lt.r2max) then ! Check if cell is closer than the farthest known particle
            npart_pos_this_node=-firstchild(icellid_out)-1

            do i=npart_pos_this_node+1,npart_pos_this_node+mass_cell(icellid_out)
               iparcell=idpart(i)
               if (iparcell .eq. iparid) cycle
               dx=abs(pos(1,iparcell)-poshere(1))
               dx=max(0.d0,min(dx,xlong-dx))
               distance2p=dx*dx
               if (distance2p .ge. r2max) cycle
               dy=abs(pos(2,iparcell)-poshere(2))
               dy=max(0.d0,min(dy,ylong-dy))
               distance2p=distance2p+dy*dy
               if (distance2p .ge. r2max) cycle
               dz=abs(pos(3,iparcell)-poshere(3))
               dz=max(0.d0,min(dz,zlong-dz))
               distance2p=distance2p+dz*dz
               ! distance2p = dx*dx + dy*dy + dz*dz
               if (distance2p .lt. r2max) then 
                  idist=nvoisins-1
                  do while (dist2(idist).gt.distance2p)
                     dist2(idist+1)=dist2(idist)
                     iparnei(idist+1)=iparnei(idist)
                     idist=idist-1
                  enddo
                  dist2(idist+1)=distance2p
                  iparnei(idist+1)=iparcell
                  r2max = dist2(nvoisins)
               endif
            enddo

         endif              
      elseif (discell2(ic).lt.r2max) then
         call walk_tree2(icellid_out,poshere,dist2,iparid,inccellpart,iparnei)
         r2max = dist2(nvoisins)
      endif
   enddo
end subroutine walk_tree2



!=======================================================================
subroutine create_tree_structure
!=======================================================================
   implicit none
   integer(kind=4) :: nlevel,inccell,idmother,ipar
   integer(kind=4) :: npart_this_node,npart_pos_this_node
   integer(kind=4) :: ncell
   real(kind=8)    :: pos_this_node(3)
   integer(kind=4) :: ttt0, ttt1, tttrate
   real(kind=8)    :: dtdtdt

   if (verbose) write(errunit,*) '    Create tree structure...'

   ! we modified to put 2*npart-1 instead of 2*npart so that AdaptaHOP can work on a 1024^3, 2*(1024^3)-1 is still an integer(kind=4), 2*(1024^3) is not 
   ncellmx=1*npart -1
   ncellbuffer=max(nint(0.1*npart),ncellbuffermin)
   allocate(idpart(npart))
   memory_used = memory_used + real(size(idpart, kind=8), kind=8)*4/1.d9
   allocate(idpart_tmp(npart))
   memory_used = memory_used + real(size(idpart_tmp, kind=8), kind=8)*4/1.d9
   allocate(mass_cell(ncellmx))
   memory_used = memory_used + real(size(mass_cell, kind=8), kind=8)*8/1.d9
   allocate(size_cell(ncellmx))
   memory_used = memory_used + real(size(size_cell, kind=8), kind=8)*8/1.d9
   allocate(pos_cell(3,ncellmx))
   memory_used = memory_used + real(size(pos_cell, kind=8), kind=8)*8/1.d9
   allocate(sister(1:ncellmx))
   memory_used = memory_used + real(size(sister, kind=8), kind=8)*4/1.d9
   allocate(firstchild(1:ncellmx))
   memory_used = memory_used + real(size(firstchild, kind=8), kind=8)*4/1.d9
   do ipar=1,npart
      idpart(ipar)=ipar
   enddo
   nlevel=0
   inccell=0
   idmother=0
   pos_this_node(1:3)=0.d0
   npart_this_node=npart
   npart_pos_this_node=0
   sizeroot=real(max(xlong,ylong,zlong),8)

   allocate(icidpart_tmp(npart))
   memory_used = memory_used + real(size(icidpart_tmp, kind=8), kind=8)*4/1.d9
   call print_memory()
   call system_clock(count=ttt0, count_rate=tttrate)

   call create_KDtree(nlevel,pos_this_node, &
   &                   npart_this_node,npart_pos_this_node,inccell, &
   &                   idmother,sizeroot*0.5d0)
   ncell=inccell

   call system_clock(count=ttt1, count_rate=tttrate)
   dtdtdt=real(ttt1-ttt0,8)/real(tttrate,8)
   
   if (verbose) write(errunit,*) '    --> total number of cells =',ncell
   if (verbose) write(errunit,'(A,F10.2,A)') "     --> ",dtdtdt," seconds to create the tree structure"

   memory_used = memory_used - real(size(idpart_tmp, kind=8), kind=8)*4/1.d9
   deallocate(idpart_tmp)
   memory_used = memory_used - real(size(icidpart_tmp, kind=8), kind=8)*4/1.d9
   call print_memory()
   deallocate(icidpart_tmp)

end subroutine create_tree_structure


!=======================================================================
recursive subroutine create_KDtree(nlevel,pos_this_node,npart_this_node, &
 &                npart_pos_this_node,inccell,idmother,size_my)
!=======================================================================
!  nlevel : level of the node in the octree. Level zero corresponds to 
!           the full box
!  pos_this_node : position of the center of this node
!  npart  : total number of particles 
!  idpart : array of dimension npart containing the id of each
!           particle. It is sorted such that neighboring particles in
!           this array belong to the same cell node.
!  idpart_tmp : temporary array of same size used as a buffer to sort
!           idpart.
!  npart_this_node : number of particles in the considered node
!  npart_pos_this_node : first position in idpart of the particles 
!           belonging to this node
!  pos :    array of dimension 3.npart giving the positions of 
!           each particle belonging to the halo
!  inccell : cell id number for the newly created structured grid site
!  pos_cell : array of dimension 3.npart (at most) corresponding
!           to the positions of each node of the structured grid
!  mass_cell : array of dimension npart (at most) giving the 
!           number of particles in each node of the structured grid
!  size_cell : array of dimension npart (at most) giving half the
!           size of the cube forming each node of the structured grid
!  sizeroot : size of the root cell (nlevel=0)
!  idmother : id of the mother cell
!  sister   : sister of a cell (at the same level, with the same mother)
!  firstchild : first child of a cell (then the next ones are found 
!             with the array sister). If it is a cell containing only
!           one particle, it gives the id of the particle.
!  ncellmx : maximum number of cells
!  megaverbose : detailed verbose mode
!=======================================================================
   implicit none

   integer(kind=4)           :: nlevel,npart_pos_this_node,npart_this_node,ipid
   real(kind=8)              :: size_my, size_child
   real(kind=8)              :: pos_this_node(3)
   integer(kind=4)           :: ipar,icid,j,inccell,nlevel_out
   integer(kind=4)           :: npart_pos_this_node_out,npart_this_node_out
   integer(kind=4)           :: incsubcell(0:7),nsubcell(0:7)
   real(kind=8)              :: pos_this_node_out(3)
   integer(kind=4)           :: idmother,idmother_out

   ! integer(kind=4)           :: ncellmx_old
   integer(kind=8)           :: ncellmx_old
   integer, allocatable      :: mass_cell_tmp(:),sister_tmp(:),firstchild_tmp(:)
   real(kind=8), allocatable :: size_cell_tmp(:),pos_cell_tmp(:,:)

   if (npart_this_node.gt.0) then
      inccell=inccell+1
      if (inccell.gt.2000000000) then
         write(errunit,*) 'WARNING: inccell is close to the maximum integer value.'
         write(errunit,*) 'Consider increasing the size of integer variables to avoid overflow.'
      endif         
      if (mod(inccell,1000000).eq.0.and.megaverbose) write(errunit,*) '    inccell=',inccell,nlevelmax
      if (inccell.gt.ncellmx) then
         ncellmx_old=ncellmx
         ncellmx=ncellmx+ncellbuffer
         if (megaverbose) write(errunit,*) &
 &          'ncellmx is too small. Increase it and reallocate arrays accordingly'
         allocate(mass_cell_tmp(ncellmx_old))
         memory_used = memory_used + real(size(mass_cell_tmp, kind=8), kind=8)*8/1.d9
         mass_cell_tmp(1:ncellmx_old)=mass_cell(1:ncellmx_old)
         memory_used = memory_used - real(size(mass_cell, kind=8), kind=8)*8/1.d9
         deallocate(mass_cell)
         allocate(mass_cell(ncellmx))
         memory_used = memory_used + real(size(mass_cell, kind=8), kind=8)*8/1.d9
         mass_cell(1:ncellmx_old)=mass_cell_tmp(1:ncellmx_old)
         memory_used = memory_used - real(size(mass_cell_tmp, kind=8), kind=8)*8/1.d9
         deallocate(mass_cell_tmp)
         allocate(sister_tmp(ncellmx_old))
         memory_used = memory_used + real(size(sister_tmp, kind=8), kind=8)*4/1.d9
         sister_tmp(1:ncellmx_old)=sister(1:ncellmx_old)
         memory_used = memory_used - real(size(sister, kind=8), kind=8)*4/1.d9
         deallocate(sister)
         allocate(sister(ncellmx))
         memory_used = memory_used + real(size(sister, kind=8), kind=8)*4/1.d9
         sister(1:ncellmx_old)=sister_tmp(1:ncellmx_old)
         memory_used = memory_used - real(size(sister_tmp, kind=8), kind=8)*4/1.d9
         deallocate(sister_tmp)
         allocate(firstchild_tmp(ncellmx_old))
         memory_used = memory_used + real(size(firstchild_tmp, kind=8), kind=8)*4/1.d9
         firstchild_tmp(1:ncellmx_old)=firstchild(1:ncellmx_old)
         memory_used = memory_used - real(size(firstchild, kind=8), kind=8)*4/1.d9
         deallocate(firstchild)
         allocate(firstchild(ncellmx))
         memory_used = memory_used + real(size(firstchild, kind=8), kind=8)*4/1.d9
         firstchild(1:ncellmx_old)=firstchild_tmp(1:ncellmx_old)
         firstchild(ncellmx_old+1:ncellmx)=0
         memory_used = memory_used - real(size(firstchild_tmp, kind=8), kind=8)*4/1.d9
         deallocate(firstchild_tmp)
         allocate(size_cell_tmp(ncellmx_old))
         memory_used = memory_used + real(size(size_cell_tmp, kind=8), kind=8)*8/1.d9
         size_cell_tmp(1:ncellmx_old)=size_cell(1:ncellmx_old)
         memory_used = memory_used - real(size(size_cell, kind=8), kind=8)*8/1.d9
         deallocate(size_cell)
         allocate(size_cell(ncellmx))
         memory_used = memory_used + real(size(size_cell, kind=8), kind=8)*8/1.d9
         size_cell(1:ncellmx_old)=size_cell_tmp(1:ncellmx_old)
         memory_used = memory_used - real(size(size_cell_tmp, kind=8), kind=8)*8/1.d9
         deallocate(size_cell_tmp)
         allocate(pos_cell_tmp(3,ncellmx_old))
         memory_used = memory_used + real(size(pos_cell_tmp, kind=8), kind=8)*8/1.d9
         pos_cell_tmp(1:3,1:ncellmx_old)=pos_cell(1:3,1:ncellmx_old)
         memory_used = memory_used - real(size(pos_cell, kind=8), kind=8)*8/1.d9
         deallocate(pos_cell)
         allocate(pos_cell(3,ncellmx))
         memory_used = memory_used + real(size(pos_cell, kind=8), kind=8)*8/1.d9
         pos_cell(1:3,1:ncellmx_old)=pos_cell_tmp(1:3,1:ncellmx_old)
         memory_used = memory_used - real(size(pos_cell_tmp, kind=8), kind=8)*8/1.d9
         deallocate(pos_cell_tmp)
      endif
      pos_cell(1:3,inccell)=pos_this_node(1:3)
      mass_cell(inccell)=npart_this_node
      ! size_cell(inccell)=2.**(-nlevel)*sizeroot*0.5d0
      size_cell(inccell)=size_my
      
      ! new initialization
      sister(inccell)=0
      firstchild(inccell)=0
      if (idmother.gt.0) then
         sister(inccell)=firstchild(idmother)
         firstchild(idmother)=inccell
      endif
!     If there is only one particle in the node or we have reach
!     maximum level of refinement, we are done
      if ((npart_this_node <= npartpercell).or.(nlevel.eq.nlevelmax)) then
         firstchild(inccell)=-(npart_pos_this_node+1)
         return
      endif
   else
      return
   endif


!  Count the number of particles in each subcell of this node
   incsubcell(0:7)=0
   do ipar=npart_pos_this_node+1,npart_pos_this_node+npart_this_node
      ipid = idpart(ipar)
      icid = 0
      if (pos(1,ipid) >= pos_this_node(1)) icid = icid + 1
      if (pos(2,ipid) >= pos_this_node(2)) icid = icid + 2
      if (pos(3,ipid) >= pos_this_node(3)) icid = icid + 4

      icidpart_tmp(ipar)=icid
      incsubcell(icid)=incsubcell(icid)+1
   enddo

!  Create the array of positions of the first particle of the lists
!  of particles belonging to each subnode
   nsubcell(0)=0
   do j=1,7
      nsubcell(j)=nsubcell(j-1)+incsubcell(j-1)
   enddo

!  Sort the array of ids (idpart) to gather the particles belonging
!  to the same subnode. Put the result in idpart_tmp.
   incsubcell(0:7)=0
   do ipar=npart_pos_this_node+1,npart_pos_this_node+npart_this_node
      icid=icidpart_tmp(ipar)
      incsubcell(icid)=incsubcell(icid)+1
      idpart_tmp(incsubcell(icid)+nsubcell(icid)+npart_pos_this_node)=idpart(ipar)
   enddo
   
!  Put back the sorted ids in idpart
   do ipar=npart_pos_this_node+1,npart_pos_this_node+npart_this_node
      idpart(ipar)=idpart_tmp(ipar)
   enddo

!  Call again the routine for the 8 subnodes:
!  Compute positions of subnodes, new level of refinement, 
!  positions in the array idpart corresponding to the subnodes,
!  and call for the treatment recursively.
   nlevel_out=nlevel+1
   idmother_out=inccell
   size_child = size_my*0.5d0
   do j=0,7
      if (incsubcell(j) <= 0) cycle

      pos_this_node_out(1)=pos_this_node(1) + size_child * pos_ref_kdtree(1,j)
      pos_this_node_out(2)=pos_this_node(2) + size_child * pos_ref_kdtree(2,j)
      pos_this_node_out(3)=pos_this_node(3) + size_child * pos_ref_kdtree(3,j)
      npart_pos_this_node_out=npart_pos_this_node+nsubcell(j)
      npart_this_node_out=incsubcell(j)
      call create_KDtree(nlevel_out,pos_this_node_out,npart_this_node_out, &
 &                npart_pos_this_node_out,inccell,idmother_out,size_child)
   enddo
end subroutine create_KDtree


!=======================================================================

!#####################################################################
subroutine close()
!#####################################################################
   ! Clean old data table is used more then once
   implicit none
   if(allocated(real_table)) deallocate(real_table)
   if(allocated(integer_table)) deallocate(integer_table)
   if(allocated(density)) deallocate(density)
   if(allocated(liste_parts)) deallocate(liste_parts)
end subroutine close




!=======================================================================
function spline(x)
!=======================================================================
   implicit none
   real(kind=8) :: spline,x,x2,x3

   if (x.le.1.d0) then
      x2 = x*x
      x3 = x2*x
      spline=1.d0 - 1.5d0*x2 + 0.75d0*x3
   elseif (x.le.2.d0) then
      x2 = (2.d0 - x)
      x3 = x2*x2*x2
      spline=0.25d0*x3
   else
      spline=0.d0
   endif

end function spline

end module neiKDtree

!=======================================================================
function icellid(xtest)
!=======================================================================
!  Compute cell id corresponding to the signs of coordinates of xtest
!  as follows :
!  (-,-,-) : 0
!  (+,-,-) : 1
!  (-,+,-) : 2
!  (+,+,-) : 3
!  (-,-,+) : 4
!  (+,-,+) : 5
!  (-,+,+) : 6
!  (+,+,+) : 7
!  For self-consistency, the array pos_ref should be defined exactly 
!  with the same conventions
!=======================================================================
   implicit none
   integer(kind=4) :: icellid,j,icellid3d(3)
   real(kind=8) :: xtest(3)
   do j=1,3
      if (xtest(j).ge.0) then
         icellid3d(j)=1
      else
         icellid3d(j)=0
      endif
   enddo
   icellid=icellid3d(1)+2*icellid3d(2)+4*icellid3d(3)
end function icellid
