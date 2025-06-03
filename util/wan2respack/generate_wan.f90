module qe_wfn_mod
  ! Read wfn in RESPACK format
  ! All the subroutines are originally from RESPACK.
  ! Modified for full k mesh with no symmetry.
  implicit none
  ! Global variables for wavefunction data
  integer:: NTK   ! Number of k-points in wannier90
  integer:: NTG   ! Maximum number of G-vectors across all k-points
  integer:: NTB   ! Number of bands
  real(8),allocatable::SK(:,:) ! SK(3,NTK) - k-point coordinates in reciprocal space
  integer,allocatable::NGI(:)  ! NGI(NTK) - Number of G-vectors for each k-point
  integer,allocatable::KGI(:,:,:)  ! KGI(3,NTG,NTK) - G-vectors for each k-point
  integer,allocatable::KG0(:,:,:)  ! KG0(3,NTG,NTK) - G-vectors for irreducible k-points
  !real(8),allocatable::E_EIGI(:,:) ! E_EIGI(NTB,NTK) - Eigenvalues (commented out)
  complex(8),allocatable::C0QE(:,:,:)! C0QE(NTG,NTB,NTK) - Wavefunction coefficients
  integer,allocatable::packing(:,:,:,:)! packing(-L1:L1,-L2:L2,-L3:L3,NTK) - Mapping array for G-vectors

  real(8)::a(3,3) ! Lattice vectors in real space

  contains

  ! Read lattice vectors from dat.lattice file
  ! Format: 3x3 matrix of lattice vectors (a1, a2, a3)
  subroutine rd_dat_lattice(file_dat_lattice)
    character(len=*), intent(in):: file_dat_lattice
    OPEN(105,FILE=file_dat_lattice)
    REWIND(105)
    READ(105,*) a(1,1),a(1,2),a(1,3)! a1 vector
    READ(105,*) a(2,1),a(2,2),a(2,3)! a2 vector
    READ(105,*) a(3,1),a(3,2),a(3,3)! a3 vector
    CLOSE(105)
  end subroutine

  ! Read k-point sampling information from dat.sample_k
  ! Format: First line is number of k-points, followed by k-point coordinates
  subroutine rd_dat_sample_k_full(file_dat_sample_k)
    character(len=*), intent(in):: file_dat_sample_k
    integer:: i, ik
    OPEN(101,FILE=file_dat_sample_k)
    rewind(101) 
    read(101,*) NTK
    allocate(SK(3,NTK))
    SK(:,:)=0.0D0 
    do ik=1,NTK
      read(101,*)(SK(i,ik),i=1,3) 
    enddo 
    close(101)
  end subroutine

  ! Read number of G-vectors for each k-point from dat.nkm
  ! Format: One number per line, representing number of G-vectors for each k-point
  subroutine rd_dat_nkm(file_dat_nkm)
    character(len=*), intent(in):: file_dat_nkm
    integer::ik 
    OPEN(132,FILE=file_dat_nkm)
    allocate(NGI(NTK))
    NGI(:)=0
    rewind(132)
    do ik=1,NTK
      read(132,*)NGI(ik) 
    enddo 
    close(132) 
    NTG=maxval(abs(NGI(:))) ! Find maximum number of G-vectors
  end subroutine

  ! Read G-vectors for full k-points from dat.kg
  ! Also creates packing array for efficient G-vector lookup
  subroutine rd_dat_kg_full(file_dat_kg)
    character(len=*), intent(in):: file_dat_kg
    integer:: i, ig, ik, i1, j1, k1
    integer:: L1, L2, L3
    integer:: NG_for_psi
    OPEN(104,FILE=file_dat_kg)
    rewind(104)
    allocate(KGI(3,NTG,NTK)); KGI(:,:,:)=0 
    do ik=1,NTK
      read(104,*) NG_for_psi 
      do ig=1,NGI(ik)
        read(104,*)(KGI(i,ig,ik),i=1,3) 
      enddo 
    enddo  
    close(104)

    ! Calculate maximum G-vector components for packing array
    L1=maxval(abs(KGI(1,:,:)))+1;write(6,*)'L1=',L1 
    L2=maxval(abs(KGI(2,:,:)))+1;write(6,*)'L2=',L2 
    L3=maxval(abs(KGI(3,:,:)))+1;write(6,*)'L3=',L3 
    allocate(packing(-L1:L1,-L2:L2,-L3:L3,NTK)); packing(:,:,:,:)=0 
    ! Create mapping from G-vector components to index
    do ik=1,NTK 
      do ig=1,NGI(ik) 
        i1=KGI(1,ig,ik);j1=KGI(2,ig,ik);k1=KGI(3,ig,ik) 
        packing(i1,j1,k1,ik)=ig 
      enddo 
    enddo 
  end subroutine

  ! Read G-vectors for irreducible k-points from dat.kg_respack
  ! Format: Similar to dat.kg but for irreducible k-points
  subroutine rd_dat_kg_irr(file_dat_kg)
    character(len=*), intent(in):: file_dat_kg
    integer:: NG_for_psi
    integer:: i, ik, ig
    OPEN(104,FILE=file_dat_kg)
    rewind(104)
    allocate(KG0(3,NTG,NTK)); KG0(:,:,:)=0 
    do ik=1,NTK
      read(104,*) NG_for_psi 
      do ig=1,NGI(ik)
        read(104,*)(KG0(i,ig,ik),i=1,3) 
      enddo 
      ! Verify number of G-vectors matches
      if(NG_for_psi/=NGI(ik)) then 
        write(6,*)'ERROR; STOP; NG_for_psi should be NGI(ik)'   
        write(6,*)'NG_for_psi=',NG_for_psi,'NG0(ik)=',NGI(ik)
        write(6,*)'ik',ik;STOP
      end if
    enddo  
    close(104)
  end subroutine

  ! Read number of bands from dat.eigenvalue
  ! Format: First line contains number of bands
  subroutine rd_dat_eigenvalue(file_dat_eigenvalue)
    character(len=*), intent(in):: file_dat_eigenvalue
    !integer::ik,ib 
    OPEN(111,FILE=file_dat_eigenvalue)
    rewind(111)
    read(111,*) NTB 
    !allocate(E_EIGI(NTB,NTK))
    !E_EIGI=0.0d0
    !do ik=1,NTK
    !  do ib=1,NTB
    !    read(111,*)E_EIGI(ib,ik)
    !  enddo!ib
    !enddo!ik
    close(111)
  end subroutine

  ! Read wavefunction coefficients from dat.wfn
  ! Format: Unformatted binary file containing complex coefficients
  subroutine rd_dat_wavefunction(file_dat_wfn)
    character(len=*), intent(in):: file_dat_wfn
    integer::ik,ib,ig,ncomp
    OPEN(102,FILE=file_dat_wfn,FORM='unformatted') 
    rewind(102)
    read(102)ncomp 
    if(ncomp/=1)then 
      write(6,*)'This program not suport ncomp/=1; then stop'
      stop
    endif 
    allocate(C0QE(NTG,NTB,NTK))
    C0QE=0.0d0
    do ik=1,NTK
      do ib=1,NTB
        read(102) (C0QE(ig,ib,ik), ig=1,NGI(ik))
      enddo!ib 
    enddo!ik          
    close(102) 
  end subroutine

  ! Convert wavefunction coefficients from QE format to RESPACK format
  ! Uses packing array to map G-vectors correctly
  subroutine make_dat_C0(C0QE, C0)
    complex(8)::C0QE(NTG,NTB,NTK),C0(NTG,NTB,NTK)
    integer:: ik, ig, jg, i, j, k
    do ik=1,NTK 
      do ig=1, NGI(ik)
        i = KG0(1,ig,ik); j = KG0(2,ig,ik); k = KG0(3,ig,ik)
        jg = packing(i,j,k,ik) 
        C0(ig,:,ik) = C0QE(jg,:,ik)
      end do
    end do
  end subroutine

end module

module wan90_chk_mod
  ! Module for reading and processing Wannier90 checkpoint files
  implicit none
  ! Constants
  integer, parameter :: dp = kind(1.0d0) ! double precision
  integer, parameter :: maxlen = 256 ! max length of line when reading file

  ! Variables from command line
  character(len=50) :: seedname ! Base name for Wannier90 files

  ! Variables from seedname.chk file
  character(len=33) :: header ! File header
  integer :: num_bands ! Total number of bands
  integer :: num_exclude_bands ! Number of bands to exclude
  integer :: num_kpts ! Number of k-points
  integer :: num_wann ! Number of Wannier functions
  integer :: nntot ! Number of nearest neighbors
  logical :: have_disentangled ! Whether disentanglement was performed
  integer :: mp_grid(3) ! Monkhorst-Pack grid dimensions
  real(kind=dp) :: real_lattice(3, 3) ! Real space lattice vectors
  real(kind=dp) :: recip_lattice(3, 3) ! Reciprocal space lattice vectors
  real(kind=dp) :: omega_invariant ! Invariant part of spread
  integer, allocatable :: exclude_bands(:) ! List of bands to exclude
  complex(kind=dp), allocatable :: u_matrix_opt(:, :, :) ! Optimal rotation matrix
  complex(kind=dp), allocatable :: u_matrix(:, :, :) ! Final rotation matrix
  real(kind=dp), allocatable :: wannier_centres(:, :) ! Wannier function centers
  real(kind=dp), allocatable :: kpt_latt(:, :) ! k-point coordinates
  character(len=20) :: checkpoint ! Checkpoint information
  integer, allocatable :: ndimwin(:) ! Number of bands in outer window
  logical, allocatable :: lwindow(:, :) ! Band window flags

contains
  ! Get seedname from command line arguments
  subroutine get_seedname()
    implicit none
    integer :: num_arg

    num_arg = COMMAND_ARGUMENT_COUNT()
    if (num_arg == 0) then
      seedname = "wannier" ! Default name if no argument provided
    elseif (num_arg == 1) then
      call GET_COMMAND_ARGUMENT(1, seedname)
    else
      print *, "command line argument error"
    end if
  end subroutine

  ! Read parameters from Wannier90 checkpoint file
  subroutine get_param_from_chk()
    implicit none
    integer :: chk_unit = 100
    integer :: i, j, k, nkp

    OPEN(unit=chk_unit, file=TRIM(seedname)//'.chk', form='unformatted')
    ! Read header and basic parameters
    READ (chk_unit) header
    READ (chk_unit) num_bands
    READ (chk_unit) num_exclude_bands
    allocate (exclude_bands(num_exclude_bands))
    READ (chk_unit) (exclude_bands(i), i=1, num_exclude_bands)
    READ (chk_unit) ((real_lattice(i, j), i=1, 3), j=1, 3)
    READ (chk_unit) ((recip_lattice(i, j), i=1, 3), j=1, 3)
    READ (chk_unit) num_kpts
    READ (chk_unit) (mp_grid(i), i=1, 3)
    allocate (kpt_latt(3, num_kpts))
    READ (chk_unit) ((kpt_latt(i, nkp), i=1, 3), nkp=1, num_kpts)
    READ (chk_unit) nntot
    READ (chk_unit) num_wann
    READ (chk_unit) checkpoint
    READ (chk_unit) have_disentangled

    ! Read disentanglement information if available
    if (have_disentangled) then
      READ (chk_unit) omega_invariant
      allocate (lwindow(num_bands, num_kpts))
      READ (chk_unit) ((lwindow(i, nkp), i=1, num_bands), nkp=1, num_kpts)
      allocate (ndimwin(num_kpts))
      READ (chk_unit) (ndimwin(nkp), nkp=1, num_kpts)
      ALLOCATE(u_matrix_opt(num_bands, num_wann, num_kpts))
      READ (chk_unit) (((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann), nkp=1, num_kpts)
    else
      ! Set default values if no disentanglement
      allocate (ndimwin(num_kpts)); ndimwin(:) = num_wann
      allocate (lwindow(num_bands, num_kpts)); lwindow(:, :) = .true.
      allocate (u_matrix_opt(num_wann, num_wann, num_kpts)); u_matrix_opt(:, :, :) = 0.0
      do i=1, num_wann
        u_matrix_opt(i, i, :) = 1.0
      enddo
    endif

    ! Read final rotation matrix and Wannier centers
    ALLOCATE(u_matrix(num_wann, num_wann, num_kpts))
    READ (chk_unit) (((u_matrix(i, j, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)
    READ (chk_unit) ! m_matrix
    ALLOCATE(wannier_centres(3, num_wann))
    READ (chk_unit) ((wannier_centres(i, j), i=1, 3), j=1, num_wann)

    CLOSE(chk_unit)
  end subroutine
end module

PROGRAM generate_wan
  ! Main program to generate Wannier functions from QE and Wannier90 data
  use qe_wfn_mod, only : rd_dat_lattice, rd_dat_sample_k_full, rd_dat_nkm, &
                         rd_dat_kg_full, rd_dat_kg_irr, rd_dat_eigenvalue, &
                         rd_dat_wavefunction
  use wan90_chk_mod, only : get_seedname, seedname, get_param_from_chk, dp
  implicit none

  ! Constants and variables
  integer, parameter:: iunit_log=148, iunit_file=100
  real(kind=dp), parameter :: bohr=0.529177d0 ! Bohr radius in Angstroms
  integer, allocatable :: ndim_exclude_low(:) ! Number of excluded bands at lower energy
  complex(kind=dp), allocatable :: UNT(:,:,:) ! Combined rotation matrix for RESPACK

  ! Initialize logging
  open(iunit_log, file='LOG.genwan', status='replace')

  ! Get seedname and read Wannier90 checkpoint file
  call get_seedname()
  write(iunit_log, '(" Reading ", a, ".chk")') trim(seedname)
  call get_param_from_chk()

  ! Read QE wavefunction data
  write(iunit_log, *) 'Reading ./dir-wfn/*'
  call rd_dat_lattice("./dir-wfn/dat.lattice")
  call rd_dat_sample_k_full("./dir-wfn/dat.sample-k")
  call rd_dat_nkm("./dir-wfn/dat.nkm")
  call rd_dat_kg_full("./dir-wfn/dat.kg")
  call rd_dat_eigenvalue("./dir-wfn/dat.eigenvalue")
  call rd_dat_wavefunction("./dir-wfn/dat.wfn")

  ! Read RESPACK G-vector data
  write(iunit_log, *) 'Reading dat.kg_respack'
  call rd_dat_kg_irr("dat.kg_respack")

  ! Verify compatibility between QE and Wannier90 data
  call check_qe_w90()

  ! Generate output files
  call write_wan_center() ! Write Wannier centers
  call write_nsnb() ! Write number of excluded bands and bands in window
  call write_umat() ! Write rotation matrices
  call write_dat_wan() ! Write Wannier functions

  close(iunit_log)

contains
  ! Check compatibility between QE and Wannier90 data
  subroutine check_qe_w90()
    use qe_wfn_mod, only : NTK, a
    use wan90_chk_mod, only : num_kpts, real_lattice
    implicit none
    ! Check number of k-points
    if (num_kpts /= NTK) then
      write(iunit_log, *) "Error: num of k points in dir-wfn and w90 are different." 
      write(iunit_log, *) NTK, num_kpts
      stop "Error: num of k points in dir-wfn and w90 are different." 
    end if
    ! Check lattice vectors
    if (any(abs(a(:,:)*bohr - real_lattice(:,:)) > 1e-4)) then
      write(iunit_log, *) "Error: lattice mismatch"
      write(iunit_log, *) a(:,:)
      write(iunit_log, *) real_lattice(:,:)
      stop "Error: lattice mismatch"
    end if
  end subroutine

  ! Write Wannier centers to file
  subroutine write_wan_center()
    use wan90_chk_mod, only : num_wann, wannier_centres
    implicit none
    real(kind=dp) :: lenconfac = 1/bohr ! Conversion factor for length units
    integer :: iw, j

    WRITE (iunit_log, *) "Writing ./dir-wan/dat.wan-center"
    OPEN(unit=iunit_file, file='./dir-wan/dat.wan-center', form='formatted', status='replace')
    WRITE (iunit_file, '(a)') '#Wannier center'
    WRITE (iunit_file, '(a)') '#1:x, 2:y, 3:z (in xyz coord)'
    do  iw = 1, num_wann
      WRITE (iunit_file, '(3F20.10)') (lenconfac*wannier_centres(j, iw), j=1, 3)
    end do
    CLOSE(iunit_file)
  end subroutine

  ! Write number of excluded bands and bands in window
  subroutine write_nsnb()
    use wan90_chk_mod, only : exclude_bands, num_exclude_bands, num_bands, num_kpts, lwindow, ndimwin
    implicit none
    integer :: n, ik
    integer:: nex_lower, nex_upper

    ! Count number of excluded bands at lower and upper energy
    nex_lower = 0
    do n=1, num_exclude_bands
      if(exclude_bands(n) /= n) exit
      nex_lower = n
    end do
    nex_upper = num_exclude_bands - nex_lower

    ! Verify band exclusion is valid
    if(nex_upper > 0) then
      if(exclude_bands(nex_lower+1) /= num_bands + nex_lower + 1) then
        write(iunit_log,*) exclude_bands(nex_lower+1)
        write(iunit_log,*) num_bands, nex_lower, nex_upper
        write(iunit_log,*) (exclude_bands(n), n=1, num_exclude_bands)
        write(iunit_log,*) "Error: Excluding intermediate bands does not work."
        stop "Error: Excluding intermediate bands does not work."
      end if
    end if

    ! Write number of excluded bands and bands in window for each k-point
    allocate(ndim_exclude_low(num_kpts))
    ndim_exclude_low = 0
    WRITE (iunit_log, *) "Writing ./dir-wan/dat.ns-nb"
    OPEN(unit=iunit_file, file="./dir-wan/dat.ns-nb", form="formatted", status="replace")
    do ik=1, num_kpts
      do n=1, num_bands
        if(lwindow(n,ik)) exit
      end do
      ndim_exclude_low(ik) = nex_lower + n - 1
      WRITE (iunit_file, '(I12, I12)') ndim_exclude_low(ik), ndimwin(ik)
    enddo
    CLOSE(iunit_file)
  end subroutine

  ! Write rotation matrices
  subroutine write_umat()
    use wan90_chk_mod, only : ndimwin, num_wann, num_kpts, u_matrix, u_matrix_opt
    implicit none
    integer ik, ib, iw, jw
    integer Mb

    ! Calculate combined rotation matrix UNT = U * Uopt
    Mb = maxval(ndimwin)
    allocate(UNT(Mb, num_wann, num_kpts)); UNT(:,:,:)=0.0d0

    do ik=1, num_kpts
      do ib=1, ndimwin(ik)
        do iw=1, num_wann
          do jw=1, num_wann
            UNT(ib, iw, ik) = UNT(ib, iw, ik) + &
                    u_matrix(jw,iw,ik) * u_matrix_opt(ib,jw,ik)
          enddo
        enddo
      enddo
    enddo !ik

    ! Write combined rotation matrix to file
    WRITE (iunit_log, *) "Writing ./dir-wan/dat.umat"
    open(iunit_file, file="./dir-wan/dat.umat", form="formatted", status="replace")
    write(iunit_file, *) num_wann
    do ik=1, num_kpts
      do ib=1, ndimwin(ik)
        write (iunit_file, *) (UNT(ib, iw, ik), iw=1, num_wann)
      end do
    end do
    close(iunit_file)
  end subroutine

  ! Write Wannier functions
  subroutine write_dat_wan()
    use wan90_chk_mod, only : ndimwin, num_wann, num_kpts
    use qe_wfn_mod,  only : C0QE, NGI, NTG, NTB, NTK, make_dat_C0
    implicit none
    integer, parameter :: iunit_file=100
    integer:: ik, ig, jw, jb
    complex(8),allocatable:: C_tilde(:,:,:), C0(:,:,:)

    ! Convert QE wavefunctions to RESPACK format
    allocate(C0(NTG,NTB,NTK));C0(:,:,:) = 0.0D0
    call make_dat_C0(C0QE, C0)

    ! Calculate Wannier functions using rotation matrices
    allocate(C_tilde(NTG,num_wann,num_kpts));C_tilde(:,:,:) = 0.0D0
    write(iunit_log, *) 'Calculating C_tilde'
    do ik=1, num_kpts
      do jw=1, num_wann
        do jb=1, ndimwin(ik)
          do ig=1, NGI(ik)
            C_tilde(ig,jw,ik) = C_tilde(ig,jw,ik) &
              + UNT(jb,jw,ik) * C0(ig,jb+ndim_exclude_low(ik),ik)
          enddo !jb
        enddo !jw
      enddo !ig
    enddo !ik

    ! Write Wannier functions to file
    write(iunit_log, '(" Writing ./dir-wan/dat.wan:  NWF, NTK = ", i5, i5)') num_wann, num_kpts
    OPEN(iunit_file, FILE='./dir-wan/dat.wan', FORM='unformatted')
    write(iunit_file) num_wann
    do ik=1, num_kpts
      write(iunit_file) ((C_tilde(ig, jw, ik), ig=1, NGI(ik)), jw=1, num_wann)
    enddo
    close(iunit_file)
    deallocate(C0, C_tilde)
  end subroutine

end program
