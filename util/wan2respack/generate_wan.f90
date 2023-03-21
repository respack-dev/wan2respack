module qe_wfn_mod
  ! Read wfn in RESPACK format
  ! All the subroutines are originally from RESPACK.
  ! Modified for full k mesh with no symmetry.
  implicit none
  integer:: NTK   ! num_kpts in wannier90
  integer:: NTG
  integer:: NTB
  real(8),allocatable::SK(:,:) !SK(3,NTK)
  integer,allocatable::NGI(:)  !NGI(NTK)
  integer,allocatable::KGI(:,:,:)  !NGI(3,NTG,NTK)
  integer,allocatable::KG0(:,:,:)  !NGI(3,NTG,NTK)
  !real(8),allocatable::E_EIGI(:,:) !E_EIGI(NTB,NTK)
  complex(8),allocatable::C0QE(:,:,:)!C0QE(NTG,NTB,NTK)
  integer,allocatable::packing(:,:,:,:)!packing(-L1:L1,-L2:L2,-L3:L3,NTK)

  real(8)::a(3,3)

  contains

  ! read dat.lattice => a
  subroutine rd_dat_lattice(file_dat_lattice)
    character(len=*), intent(in):: file_dat_lattice
    OPEN(105,FILE=file_dat_lattice)
    REWIND(105)
    READ(105,*) a(1,1),a(1,2),a(1,3)!a1 vector
    READ(105,*) a(2,1),a(2,2),a(2,3)!a2 vector
    READ(105,*) a(3,1),a(3,2),a(3,3)!a3 vector
    CLOSE(105)
  end subroutine

  ! read dat.sample_k => NTK, SK(1:NTK)
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

  ! read dat.nkm => NGI(1:NTK), NTG
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
    NTG=maxval(abs(NGI(:))) 
  end subroutine

  ! read dat.kg (g vectors for full k points (wfc for w90))
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
      do ig=1,NGI(ik)!NG_for_psi 
        read(104,*)(KGI(i,ig,ik),i=1,3) 
      enddo 
    enddo  
    close(104)

    L1=maxval(abs(KGI(1,:,:)))+1;write(6,*)'L1=',L1 
    L2=maxval(abs(KGI(2,:,:)))+1;write(6,*)'L2=',L2 
    L3=maxval(abs(KGI(3,:,:)))+1;write(6,*)'L3=',L3 
    allocate(packing(-L1:L1,-L2:L2,-L3:L3,NTK)); packing(:,:,:,:)=0 
    do ik=1,NTK 
      do ig=1,NGI(ik) 
        i1=KGI(1,ig,ik);j1=KGI(2,ig,ik);k1=KGI(3,ig,ik) 
        packing(i1,j1,k1,ik)=ig 
      enddo 
    enddo 
  end subroutine

  ! read dat.kg_respack (g vectors for irreducible k points)
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
      if(NG_for_psi/=NGI(ik)) then 
        write(6,*)'ERROR; STOP; NG_for_psi should be NGI(ik)'   
        write(6,*)'NG_for_psi=',NG_for_psi,'NG0(ik)=',NGI(ik)
        write(6,*)'ik',ik;STOP
      end if
    enddo  
    close(104)
  end subroutine

  ! read dat.eigenvalue => NTB
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

  ! read dat.wfn => C0QE
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

  ! C0QE => C0
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
  implicit none
  ! parameter
  integer, parameter :: dp = kind(1.0d0) ! double precision
  integer, parameter :: maxlen = 256 ! max length of line when reading file

  ! from commandline
  character(len=50) :: seedname

  ! from seedname.chk
  character(len=33) :: header
  integer :: num_bands
  integer :: num_exclude_bands
  integer :: num_kpts
  integer :: num_wann
  integer :: nntot
  logical :: have_disentangled
  integer :: mp_grid(3)
  real(kind=dp) :: real_lattice(3, 3), recip_lattice(3, 3)
  real(kind=dp) :: omega_invariant
  integer, allocatable :: exclude_bands(:)
  complex(kind=dp), allocatable :: u_matrix_opt(:, :, :)
  complex(kind=dp), allocatable :: u_matrix(:, :, :)
  real(kind=dp), allocatable :: wannier_centres(:, :)
  real(kind=dp), allocatable :: kpt_latt(:, :)
  character(len=20) :: checkpoint
  integer, allocatable :: ndimwin(:)
  logical, allocatable :: lwindow(:, :)

contains
  ! from wannier90: set seedname
  subroutine get_seedname()
    ! -------------------------------------
    ! Get seedname from the commandline
    ! -------------------------------------
    implicit none
    integer :: num_arg

    num_arg = COMMAND_ARGUMENT_COUNT()
    if (num_arg == 0) then
      seedname = "wannier"
    elseif (num_arg == 1) then
      call GET_COMMAND_ARGUMENT(1, seedname)
    else
      print *, "command line argument error"
    end if
  end subroutine

  ! from wannier90: read seedname.chk
  subroutine get_param_from_chk()
    implicit none
    integer :: chk_unit = 100
    integer :: i, j, k, nkp


    OPEN(unit=chk_unit, file=TRIM(seedname)//'.chk', form='unformatted')
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
    if (have_disentangled) then
      READ (chk_unit) omega_invariant
      allocate (lwindow(num_bands, num_kpts))
      READ (chk_unit) ((lwindow(i, nkp), i=1, num_bands), nkp=1, num_kpts)
      allocate (ndimwin(num_kpts))
      READ (chk_unit) (ndimwin(nkp), nkp=1, num_kpts)
      ALLOCATE(u_matrix_opt(num_bands, num_wann, num_kpts))
      READ (chk_unit) (((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann), nkp=1, num_kpts)
    else
      allocate (ndimwin(num_kpts)); ndimwin(:) = num_wann
      allocate (lwindow(num_bands, num_kpts)); lwindow(:, :) = .true.
      allocate (u_matrix_opt(num_wann, num_wann, num_kpts)); u_matrix_opt(:, :, :) = 0.0
      do i=1, num_wann
        u_matrix_opt(i, i, :) = 1.0
      enddo
    endif
    ALLOCATE(u_matrix(num_wann, num_wann, num_kpts))
    READ (chk_unit) (((u_matrix(i, j, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)
    READ (chk_unit) ! m_matrix
    ALLOCATE(wannier_centres(3, num_wann))
    READ (chk_unit) ((wannier_centres(i, j), i=1, 3), j=1, num_wann)

    CLOSE(chk_unit)

  end subroutine
end module

PROGRAM generate_wan
  use qe_wfn_mod, only : rd_dat_lattice, rd_dat_sample_k_full, rd_dat_nkm, &
                         rd_dat_kg_full, rd_dat_kg_irr, rd_dat_eigenvalue, &
                         rd_dat_wavefunction
  use wan90_chk_mod, only : get_seedname, seedname, get_param_from_chk, dp
  implicit none

  integer, parameter:: iunit_log=148, iunit_file=100
  real(kind=dp), parameter :: bohr=0.529177d0
  integer, allocatable :: ndim_exclude_low(:)
  complex(kind=dp), allocatable :: UNT(:,:,:) ! U matrix for Respack

  ! Logging file is LOG.genwan
  open(iunit_log, file='LOG.genwan', status='replace')

  call get_seedname()
  write(iunit_log, '(" Reading ", a, ".chk")') trim(seedname)
  call get_param_from_chk()

  write(iunit_log, *) 'Reading ./dir-wfn/*'
  call rd_dat_lattice("./dir-wfn/dat.lattice")
  call rd_dat_sample_k_full("./dir-wfn/dat.sample-k")
  call rd_dat_nkm("./dir-wfn/dat.nkm")
  call rd_dat_kg_full("./dir-wfn/dat.kg")
  call rd_dat_eigenvalue("./dir-wfn/dat.eigenvalue")
  call rd_dat_wavefunction("./dir-wfn/dat.wfn")

  write(iunit_log, *) 'Reading dat.kg_respack'
  call rd_dat_kg_irr("dat.kg_respack")

  call check_qe_w90()

  call write_wan_center()
  call write_nsnb()
  call write_umat()
  call write_dat_wan()

  !call SYSTEM("rm -r dir-wfn")
  !write(iunit_log, *) './dir-wfn was removed'
  close(iunit_log)

contains
  subroutine check_qe_w90()
    ! check qe(dir-wfn) and w90(seedname.chk) compatibility
    use qe_wfn_mod, only : NTK, a
    use wan90_chk_mod, only : num_kpts, real_lattice
    implicit none
    if (num_kpts /= NTK) then
      write(iunit_log, *) "Error: num of k points in dir-wfn and w90 are different." 
      write(iunit_log, *) NTK, num_kpts
      stop "Error: num of k points in dir-wfn and w90 are different." 
    end if
    if (any(abs(a(:,:)*bohr - real_lattice(:,:)) > 1e-4)) then
      write(iunit_log, *) "Error: lattice mismatch"
      write(iunit_log, *) a(:,:)
      write(iunit_log, *) real_lattice(:,:)
      stop "Error: lattice mismatch"
    end if
  end subroutine

  subroutine write_wan_center()
    use wan90_chk_mod, only : num_wann, wannier_centres
    implicit none
    ! reciprocal of length unit. Ang 1.0, Bohr 1.0/bohr
    real(kind=dp) :: lenconfac = 1/bohr
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

  subroutine write_nsnb()
    use wan90_chk_mod, only : exclude_bands, num_exclude_bands, num_bands, num_kpts, lwindow, ndimwin
    implicit none
    integer :: n, ik
    integer:: nex_lower, nex_upper

    nex_lower = 0
    do n=1, num_exclude_bands
      if(exclude_bands(n) /= n) exit
      nex_lower = n
    end do
    nex_upper = num_exclude_bands - nex_lower
    if(nex_upper > 0) then
      if(exclude_bands(nex_lower+1) /= num_bands + nex_lower + 1) then
        write(iunit_log,*) exclude_bands(nex_lower+1)
        write(iunit_log,*) num_bands, nex_lower, nex_upper
        write(iunit_log,*) (exclude_bands(n), n=1, num_exclude_bands)
        write(iunit_log,*) "Error: Excluding intermediate bands does not work."
        stop "Error: Excluding intermediate bands does not work."
      end if
    end if

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

  subroutine write_umat()
    use wan90_chk_mod, only : ndimwin, num_wann, num_kpts, u_matrix, u_matrix_opt
    implicit none
    integer ik, ib, iw, jw
    integer Mb

    Mb = maxval(ndimwin)
    allocate(UNT(Mb, num_wann, num_kpts)); UNT(:,:,:)=0.0d0

    do ik=1, num_kpts
      ! UNT^(ik) = U^(ik) * Uopt^(ik)
      do ib=1, ndimwin(ik)
        do iw=1, num_wann
          ! UNT(ib, iw) = sum_jw U(jw, iw) * Uopt(ib, jw)
          do jw=1, num_wann
            UNT(ib, iw, ik) = UNT(ib, iw, ik) + &
                    u_matrix(jw,iw,ik) * u_matrix_opt(ib,jw,ik)
          enddo
        enddo
      enddo
    enddo !ik

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

  subroutine write_dat_wan()
    use wan90_chk_mod, only : ndimwin, num_wann, num_kpts
    use qe_wfn_mod,  only : C0QE, NGI, NTG, NTB, NTK, make_dat_C0
    implicit none
    integer, parameter :: iunit_file=100
    integer:: ik, ig, jw, jb
    complex(8),allocatable:: C_tilde(:,:,:), C0(:,:,:)
    ! packing is RESPACK packing
    ! QE G vectors for Wannier90 (dat.kg) ==> RESPACK packing
    allocate(C0(NTG,NTB,NTK));C0(:,:,:) = 0.0D0
    call make_dat_C0(C0QE, C0)

    allocate(C_tilde(NTG,num_wann,num_kpts));C_tilde(:,:,:) = 0.0D0
    write(iunit_log, *) 'Calculating C_tilde'
    !write(iunit_log, *) num_kpts, NTK
    !write(iunit_log, *) NTB, maxval(ndimwin(:))
    !write(iunit_log, *) ndim_exclude_low
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
