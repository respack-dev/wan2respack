! Module for reading wavefunction data files
! This module provides functions to read various data files used in RESPACK:
! 1. dat.symmetry: Symmetry operations
! 2. dat.lattice: Lattice vectors
! 3. dat.sample-k: Sample k-points
! 4. dat.nkm: Number of G-vectors
! 5. dat.kg: G-vectors
!
! Dependencies
! -----------
! - LAPACK: Required for matrix operations (invmat)
! - Standard Fortran libraries
!
! Notes
! -----
! - All numerical operations use double precision (real(8))
! - File paths are relative to the ./dir-wfn/ directory
! - Global variables are used to maintain data across subroutines
! - Error handling uses STOP statements with error messages
module m_rd_dat_wfn
  implicit none
  ! Public subroutines for reading data files
  public::rd_dat_symmetry  ! Read symmetry operations
  public::rd_dat_lattice   ! Read lattice vectors
  public::rd_dat_sample_k  ! Read sample k-points
  public::rd_dat_nkm      ! Read number of G-vectors
  public::rd_dat_kg       ! Read G-vectors

  ! Symmetry data
  ! ------------
  integer,public::nsymq    ! Number of symmetry operations
  integer,public::nnp      ! Number of non-primitive translations
  integer,public,allocatable::rg(:,:,:)  ! Rotation matrices (3,3,nsymq)
  integer,public,allocatable::pg(:,:)    ! Translation vectors (3,nsymq)
  real(8),public,allocatable::rginv(:,:,:)  ! Inverse rotation matrices (3,3,nsymq)

  ! Band calculation parameters
  ! -------------------------
  real(8),public::Ecut_for_psi  ! Energy cutoff for wavefunctions
  real(8),public::FermiEnergy   ! Fermi energy
  real(8),public::Etot         ! Total energy

  ! Lattice vectors
  ! -------------
  real(8),public::a1(3)  ! First lattice vector
  real(8),public::a2(3)  ! Second lattice vector
  real(8),public::a3(3)  ! Third lattice vector
  real(8),public::b1(3)  ! First reciprocal lattice vector
  real(8),public::b2(3)  ! Second reciprocal lattice vector
  real(8),public::b3(3)  ! Third reciprocal lattice vector
  real(8),public::VOLUME ! Unit cell volume

  ! Sample k-points data
  ! ------------------
  integer,public::Nk_irr  ! Number of irreducible k-points
  integer,public::NTK     ! Total number of k-points
  integer,public::nkb1    ! Sampling points along b1 vector
  integer,public::nkb2    ! Sampling points along b2 vector
  integer,public::nkb3    ! Sampling points along b3 vector
  integer,public::Na1     ! Lattice translations along a1
  integer,public::Na2     ! Lattice translations along a2
  integer,public::Na3     ! Lattice translations along a3
  real(8),public,allocatable::SKI(:,:)    ! Irreducible k-points (3,Nk_irr)
  real(8),public,allocatable::SK0(:,:)    ! All k-points (3,NTK)
  integer,public,allocatable::numirr(:)   ! Irreducible k-point index (NTK)
  integer,public,allocatable::numrot(:)   ! Rotation operation index (NTK)
  integer,public,allocatable::trs(:)      ! Time-reversal symmetry flag (NTK)
  integer,public,allocatable::RW(:,:)     ! Rewind vectors (3,NTK)
  integer,public,allocatable::numMK(:)    ! Mapping to irreducible k-points (Nk_irr)

  ! Eigenvalue data
  ! -------------
  integer,public::NTB  ! Number of bands
  real(8),public,allocatable::E_EIGI(:,:)  ! Eigenvalues (NTB,Nk_irr)

  ! G-vectors data
  ! ------------
  integer,public::NTG  ! Maximum number of G-vectors
  integer,public,allocatable::NGI(:)  ! Number of G-vectors per k-point (Nk_irr)
  integer,public::L1  ! G-vector limit along b1
  integer,public::L2  ! G-vector limit along b2
  integer,public::L3  ! G-vector limit along b3
  integer,public,allocatable::KGI(:,:,:)  ! G-vectors for irreducible k-points (3,NTG,Nk_irr)
  integer,public,allocatable::KG0(:,:,:)  ! G-vectors for all k-points (3,NTG,NTK)
  integer,public,allocatable::NG0(:)      ! Number of G-vectors per k-point (NTK)
  integer,public,allocatable::packing(:,:,:,:)  ! G-vector packing array (-L1:L1,-L2:L2,-L3:L3,Nk_irr)

  ! FFT grid parameters
  ! ----------------
  integer,public::nwx2,nwy2,nwz2  ! FFT grid dimensions

  ! Wavefunction data
  ! --------------
  integer,public::ncomp  ! Number of components
  complex(4),public,allocatable::CIR(:,:,:)  ! Wavefunction coefficients (NTG,NTB,Nk_irr)
  complex(8),allocatable::CIRtmp(:)  ! Temporary array for wavefunction coefficients (NTG)

  ! Logging
  ! ------
  integer,public::iunit_log  ! Log file unit number

contains
!
subroutine rd_dat_symmetry 
  implicit none 
  integer::i,j,iop
  OPEN(100,FILE='./dir-wfn/dat.symmetry') 
  rewind(100) 
  read(100,*)nsymq 
  read(100,*)nnp 
  allocate(rg(3,3,nsymq));rg=0
  allocate(pg(3,nsymq));pg=0
  allocate(rginv(3,3,nsymq));rginv=0.0d0 
  do iop=1,nsymq
   read(100,*)((rg(i,j,iop),i=1,3),j=1,3) 
   read(100,*)(pg(i,iop),i=1,3)   
  enddo 
  close(100) 
  rginv=rg 
  do iop=1,nsymq
   call invmat(3,rginv(1,1,iop)) 
  enddo 
  do iop=1,nsymq
   write(iunit_log,*) iop
   do i=1,3
    write(iunit_log,'(3I5,1x,3F15.10)')(rg(i,j,iop),j=1,3),(rginv(i,j,iop),j=1,3)
   enddo 
  enddo 
end subroutine
!
! Read lattice vectors from dat.lattice and calculate reciprocal lattice vectors
!
! This subroutine reads the lattice vectors from dat.lattice and calculates
! the reciprocal lattice vectors and unit cell volume.
!
! File format:
! - First line: a1 vector (3 components)
! - Second line: a2 vector (3 components)
! - Third line: a3 vector (3 components)
!
! Global variables modified:
! - a1, a2, a3: Direct lattice vectors
! - b1, b2, b3: Reciprocal lattice vectors
! - VOLUME: Unit cell volume
!
! Algorithm:
! 1. Read lattice vectors from file
! 2. Calculate reciprocal lattice vectors using cross products
! 3. Calculate unit cell volume
! 4. Normalize reciprocal lattice vectors by 2π/VOLUME
!
! Notes:
! - All vectors are stored in Cartesian coordinates
! - Reciprocal lattice vectors are calculated using the formula:
!   b1 = 2π(a2 × a3)/V, b2 = 2π(a3 × a1)/V, b3 = 2π(a1 × a2)/V
!   where V is the unit cell volume
subroutine rd_dat_lattice 
  implicit none 
  real(8),parameter::tpi=2.0d0*dacos(-1.0d0)  ! 2π constant
  OPEN(105,FILE='./dir-wfn/dat.lattice') 
  REWIND(105)
  READ(105,*) a1(1),a1(2),a1(3)!a1 vector
  READ(105,*) a2(1),a2(2),a2(3)!a2 vector
  READ(105,*) a3(1),a3(2),a3(3)!a3 vector
  CLOSE(105)
  !--
  ! Calculate reciprocal lattice vectors using cross products
  call OUTER_PRODUCT(a2(1),a3(1),b1(1))  ! b1 = a2 × a3
  VOLUME=a1(1)*b1(1)+a1(2)*b1(2)+a1(3)*b1(3)  ! V = a1·(a2 × a3)
  b1(:)=b1(:)*tpi/VOLUME  ! Normalize b1 by 2π/V
  call OUTER_PRODUCT(a3(1),a1(1),b2(1))  ! b2 = a3 × a1
  b2(:)=b2(:)*tpi/VOLUME  ! Normalize b2 by 2π/V
  call OUTER_PRODUCT(a1(1),a2(1),b3(1))  ! b3 = a1 × a2
  b3(:)=b3(:)*tpi/VOLUME  ! Normalize b3 by 2π/V
end subroutine
!
! Read sample k-points from dat.sample-k and generate all k-points
!
! This subroutine reads irreducible k-points from dat.sample-k and generates
! all k-points by applying symmetry operations and time-reversal symmetry.
!
! File format:
! - First line: Nk_irr (number of irreducible k-points)
! - Next Nk_irr lines: k-point coordinates (3 components each)
!
! Global variables modified:
! - Nk_irr: Number of irreducible k-points
! - NTK: Total number of k-points
! - SKI: Irreducible k-points (3,Nk_irr)
! - SK0: All k-points (3,NTK)
! - numirr: Mapping from k-point to irreducible k-point
! - numrot: Rotation operation index for each k-point
! - trs: Time-reversal symmetry flag for each k-point
! - RW: Rewind vectors for each k-point
! - numMK: Mapping from irreducible k-point to first generated k-point
!
! Algorithm:
! 1. Read irreducible k-points from file
! 2. Estimate total number of k-points (NTK)
! 3. For each irreducible k-point:
!    a. Apply all symmetry operations
!    b. Apply time-reversal symmetry
!    c. Check for duplicates and rewind if necessary
! 4. Calculate k-point sampling parameters (nkb1, nkb2, nkb3)
!
! Notes:
! - K-points are in reduced coordinates (in units of reciprocal lattice vectors)
! - Time-reversal symmetry is applied as k → -k
! - K-points are rewound to the first Brillouin zone
! - Duplicate k-points are removed
subroutine rd_dat_sample_k 
  implicit none 
  integer::i,ik,jk,iop,iik   
  real(8)::ktmp(3) 
  integer::RWtmp(3) 
  integer::initial_flg!20180316 
  OPEN(101,FILE='./dir-wfn/dat.sample-k') 
  rewind(101) 
  read(101,*)Nk_irr 
  allocate(SKI(3,Nk_irr));SKI(:,:)=0.0D0 
  do ik=1,Nk_irr 
   read(101,*)(SKI(i,ik),i=1,3) 
  enddo 
  close(101) 
  !--
  ! Estimate total number of k-points
  call est_NTK(Nk_irr,nsymq,SKI(1,1),rg(1,1,1),NTK)
  write(iunit_log,*)'Estimated NTK=',NTK 
  !--
  ! Allocate arrays for k-point data
  allocate(SK0(3,NTK));SK0(:,:)=0.0d0
  allocate(numirr(NTK));numirr(:)=0
  allocate(numrot(NTK));numrot(:)=0
  allocate(trs(NTK));trs(:)=0
  allocate(RW(3,NTK));RW(:,:)=0
  allocate(numMK(Nk_irr));numMK=0!20180316 
  jk=0
  ! Generate all k-points
  do ik=1,Nk_irr
   initial_flg=0!20180316 
   do iop=1,nsymq
    ! Apply symmetry operation
    ktmp(:)=0.0d0; RWtmp(:)=0  
    ktmp(1)=dble(rg(1,1,iop))*SKI(1,ik)+dble(rg(1,2,iop))*SKI(2,ik)+dble(rg(1,3,iop))*SKI(3,ik)
    ktmp(2)=dble(rg(2,1,iop))*SKI(1,ik)+dble(rg(2,2,iop))*SKI(2,ik)+dble(rg(2,3,iop))*SKI(3,ik)
    ktmp(3)=dble(rg(3,1,iop))*SKI(1,ik)+dble(rg(3,2,iop))*SKI(2,ik)+dble(rg(3,3,iop))*SKI(3,ik)
    call kcheck(ktmp(1),RWtmp(1))!rewind check 
    ! Check for duplicates
    do iik=1,jk
     if(abs(SK0(1,iik)-ktmp(1))<1.0d-4.and.abs(SK0(2,iik)-ktmp(2))<1.0d-4.and.abs(SK0(3,iik)-ktmp(3))<1.0d-4) goto 1000
    enddo!iik
    ! Add new k-point
    jk=jk+1
    SK0(:,jk)=ktmp(:)
    numirr(jk)=ik
    numrot(jk)=iop
    trs(jk)=1
    RW(:,jk)=RWtmp(:)
    if(initial_flg.eq.0)then
     numMK(ik)=jk
     initial_flg=1 
    endif 
    ! Apply time-reversal symmetry
1000 ktmp(:)=0.0d0;RWtmp(:)=0  
    ktmp(1)=dble(rg(1,1,iop))*SKI(1,ik)+dble(rg(1,2,iop))*SKI(2,ik)+dble(rg(1,3,iop))*SKI(3,ik)
    ktmp(2)=dble(rg(2,1,iop))*SKI(1,ik)+dble(rg(2,2,iop))*SKI(2,ik)+dble(rg(2,3,iop))*SKI(3,ik)
    ktmp(3)=dble(rg(3,1,iop))*SKI(1,ik)+dble(rg(3,2,iop))*SKI(2,ik)+dble(rg(3,3,iop))*SKI(3,ik)
    call kcheck_trs(ktmp(1),RWtmp(1))!rewind check modified 20170316  
    ! Check for duplicates
    do iik=1,jk
     if(abs(SK0(1,iik)-(-ktmp(1)))<1.0d-4.and.abs(SK0(2,iik)-(-ktmp(2)))<1.0d-4.and.abs(SK0(3,iik)-(-ktmp(3)))<1.0d-4) goto 2000
    enddo!iik
    ! Add new k-point
    jk=jk+1
    SK0(:,jk)=-ktmp(:) 
    numirr(jk)=ik
    numrot(jk)=iop
    trs(jk)=-1
    RW(:,jk)=RWtmp(:) 
2000 enddo!iop  
  enddo!ik
  ! Calculate k-point sampling parameters
  call est_nkbi(NTK,SK0(1,1),nkb1,nkb2,nkb3)  
  Na1=nkb1/2; Na2=nkb2/2; Na3=nkb3/2
  !--
  ! Verify total number of k-points
  if(NTK/=jk) then 
   write(iunit_log,*)'ERROR;STOP;NTK should be jk'   
   write(iunit_log,*)'NTK=',NTK,'jk=',jk;STOP
  endif 
end subroutine
!--
! Read number of G-vectors for each k-point from dat.nkm
!
! This subroutine reads the number of G-vectors for each irreducible k-point
! from dat.nkm and determines the maximum number of G-vectors (NTG).
!
! File format:
! - Nk_irr lines, each containing the number of G-vectors for one k-point
!
! Global variables modified:
! - NGI: Number of G-vectors for each irreducible k-point (Nk_irr)
! - NTG: Maximum number of G-vectors across all k-points
!
! Algorithm:
! 1. Allocate array for number of G-vectors
! 2. Read number of G-vectors for each k-point
! 3. Calculate maximum number of G-vectors
!
! Notes:
! - NTG is used to allocate arrays for G-vectors in other subroutines
! - The number of G-vectors may vary for different k-points
subroutine rd_dat_nkm 
  implicit none 
  integer::ik 
  OPEN(132,FILE='./dir-wfn/dat.nkm') 
  allocate(NGI(Nk_irr));NGI(:)=0
  rewind(132)
  do ik=1,Nk_irr 
   read(132,*)NGI(ik) 
  enddo 
  close(132) 
  NTG=maxval(abs(NGI(:))) 
end subroutine
!--
! Read G-vectors from dat.kg and generate G-vectors for all k-points
!
! This subroutine reads G-vectors for irreducible k-points from dat.kg and
! generates G-vectors for all k-points by applying symmetry operations.
!
! File format:
! For each irreducible k-point:
! - First line: Number of G-vectors
! - Next lines: G-vector coordinates (3 components each)
!
! Global variables modified:
! - KGI: G-vectors for irreducible k-points (3,NTG,Nk_irr)
! - KG0: G-vectors for all k-points (3,NTG,NTK)
! - NG0: Number of G-vectors for each k-point (NTK)
! - packing: G-vector packing array (-L1:L1,-L2:L2,-L3:L3,Nk_irr)
! - nwx2, nwy2, nwz2: FFT grid dimensions
! - Ecut_for_psi: Energy cutoff for wavefunctions
!
! Algorithm:
! 1. Read G-vectors for irreducible k-points
! 2. Calculate energy cutoff from G-vectors
! 3. Create packing array for G-vectors
! 4. Calculate FFT grid dimensions
! 5. Generate G-vectors for all k-points:
!    a. Apply symmetry operations
!    b. Apply time-reversal symmetry if needed
!
! Notes:
! - G-vectors are in reduced coordinates (in units of reciprocal lattice vectors)
! - The packing array is used for efficient G-vector lookup
! - FFT grid dimensions are chosen to be compatible with 2,3,5 factors
! - G-vectors are generated to satisfy the energy cutoff condition
subroutine rd_dat_kg 
  implicit none 
  integer::i,ig,ik,jk,iop,i1,j1,k1
  integer::algn235
  real(8)::ktmp(3)
  real(8)::tmp,d1,d2,d3,qwf   
  real(8)::h1(3),h2(3),h3(3)   
  integer::NG_for_psi 
  integer,allocatable::KGtmp(:,:)!KGtmp(3,NTG) 
  real(8),allocatable::LKGI(:,:)!LKGI(NTG,Nk_irr) 
  OPEN(104,FILE='./dir-wfn/dat.kg') 
  rewind(104) 
  allocate(KGI(3,NTG,Nk_irr));KGI(:,:,:)=0 
  do ik=1,Nk_irr 
   read(104,*)NG_for_psi 
   !NGI(ik)=NG_for_psi 
   do ig=1,NGI(ik)!NG_for_psi 
    read(104,*)(KGI(i,ig,ik),i=1,3) 
   enddo 
  enddo  
  close(104)
  !--
  ! Calculate energy cutoff from G-vectors
  allocate(LKGI(NTG,Nk_irr));LKGI=0.0d0 
  do ik=1,Nk_irr 
   do ig=1,NGI(ik) 
    ktmp(1)=(SKI(1,ik)+dble(KGI(1,ig,ik)))*b1(1)+(SKI(2,ik)+dble(KGI(2,ig,ik)))*b2(1)+(SKI(3,ik)+dble(KGI(3,ig,ik)))*b3(1) 
    ktmp(2)=(SKI(1,ik)+dble(KGI(1,ig,ik)))*b1(2)+(SKI(2,ik)+dble(KGI(2,ig,ik)))*b2(2)+(SKI(3,ik)+dble(KGI(3,ig,ik)))*b3(2) 
    ktmp(3)=(SKI(1,ik)+dble(KGI(1,ig,ik)))*b1(3)+(SKI(2,ik)+dble(KGI(2,ig,ik)))*b2(3)+(SKI(3,ik)+dble(KGI(3,ig,ik)))*b3(3) 
    LKGI(ig,ik)=ktmp(1)**2+ktmp(2)**2+ktmp(3)**2
   enddo!ig 
  enddo!ik 
  Ecut_for_psi=maxval(LKGI(:,:))+1.0d-8
  deallocate(LKGI) 
  !--
  ! Create packing array for G-vectors
  L1=maxval(abs(KGI(1,:,:)))+1;write(iunit_log,*)'L1=',L1 
  L2=maxval(abs(KGI(2,:,:)))+1;write(iunit_log,*)'L2=',L2 
  L3=maxval(abs(KGI(3,:,:)))+1;write(iunit_log,*)'L3=',L3 
  allocate(packing(-L1:L1,-L2:L2,-L3:L3,Nk_irr)); packing(:,:,:,:)=0 
  do ik=1,Nk_irr 
   do ig=1,NGI(ik) 
    i1=KGI(1,ig,ik);j1=KGI(2,ig,ik);k1=KGI(3,ig,ik) 
    packing(i1,j1,k1,ik)=ig 
   enddo 
  enddo 
  ! Calculate FFT grid dimensions
  tmp=dsqrt(dot_product(a1,a1))
  h1(:)=a1(:)/tmp 
  tmp=dsqrt(dot_product(a2,a2))
  h2(:)=a2(:)/tmp 
  tmp=dsqrt(dot_product(a3,a3)) 
  h3(:)=a3(:)/tmp 
  d1=abs(dot_product(b1,h1)) 
  d2=abs(dot_product(b2,h2)) 
  d3=abs(dot_product(b3,h3)) 
  qwf=2.0d0*dsqrt(Ecut_for_psi) 
  nwx2=algn235(int(qwf/d1)+1) 
  nwy2=algn235(int(qwf/d2)+1) 
  nwz2=algn235(int(qwf/d3)+1) 
  ! Generate G-vectors for all k-points
  allocate(NG0(NTK));NG0(:)=0
  allocate(KG0(3,NTG,NTK));KG0(:,:,:)=0 
  allocate(KGtmp(3,NTG));KGtmp(:,:)=0 
  do jk=1,NTK 
   if(trs(jk)==1) then 
    ! Apply symmetry operation
    ik=numirr(jk)
    iop=numrot(jk) 
    ktmp(1)=dble(rg(1,1,iop))*SKI(1,ik)+dble(rg(1,2,iop))*SKI(2,ik)+dble(rg(1,3,iop))*SKI(3,ik)+dble(RW(1,jk)) 
    ktmp(2)=dble(rg(2,1,iop))*SKI(1,ik)+dble(rg(2,2,iop))*SKI(2,ik)+dble(rg(2,3,iop))*SKI(3,ik)+dble(RW(2,jk))  
    ktmp(3)=dble(rg(3,1,iop))*SKI(1,ik)+dble(rg(3,2,iop))*SKI(2,ik)+dble(rg(3,3,iop))*SKI(3,ik)+dble(RW(3,jk))  
    call make_KG0(NTG,b1(1),b2(1),b3(1),Ecut_for_psi,ktmp(1),ktmp(2),ktmp(3),KG0(1,1,jk),NG_for_psi)
    if(NG_for_psi/=NGI(ik)) then 
     write(iunit_log,*)'ERROR; STOP; NG_for_psi should be NGI(ik)'   
     write(iunit_log,*)'NG_for_psi=',NG_for_psi,'NG0(ik)=',NGI(ik)
     write(iunit_log,*)'ik,jk',ik,jk;STOP
    endif 
    NG0(jk)=NG_for_psi  
   elseif(trs(jk)==-1)then  
    ! Apply time-reversal symmetry
    ik=numirr(jk);iop=numrot(jk) 
    ktmp(1)=rg(1,1,iop)*SKI(1,ik)+rg(1,2,iop)*SKI(2,ik)+rg(1,3,iop)*SKI(3,ik)+dble(RW(1,jk)) 
    ktmp(2)=rg(2,1,iop)*SKI(1,ik)+rg(2,2,iop)*SKI(2,ik)+rg(2,3,iop)*SKI(3,ik)+dble(RW(2,jk))  
    ktmp(3)=rg(3,1,iop)*SKI(1,ik)+rg(3,2,iop)*SKI(2,ik)+rg(3,3,iop)*SKI(3,ik)+dble(RW(3,jk))  
    KGtmp(:,:)=0 
    call make_KG0(NTG,b1(1),b2(1),b3(1),Ecut_for_psi,ktmp(1),ktmp(2),ktmp(3),KGtmp(1,1),NG_for_psi)
    if(NG_for_psi/=NGI(ik))then 
     write(iunit_log,*)'ERROR; STOP; NG_for_psi should be NGI(ik)'   
     write(iunit_log,*)'NG_for_psi=',NG_for_psi,'NGI(ik)=',NGI(ik);STOP
    endif 
    NG0(jk)=NG_for_psi  
    KG0(:,:,jk)=-KGtmp(:,:)!notice on '-' sign 
   endif 
  enddo!jk 
  deallocate(KGtmp) 
end subroutine
!--
! Check and rewind k-point coordinates to the first Brillouin zone
!
! This subroutine checks if a k-point is outside the first Brillouin zone
! and rewinds it back if necessary. The rewind operation is performed by
! adding or subtracting integer multiples of reciprocal lattice vectors.
!
! Parameters:
! - ktmp(3): Input/output k-point coordinates (in reduced coordinates)
! - RWtmp(3): Output rewind vectors indicating how many times the k-point
!             was rewound along each direction
!
! Algorithm:
! For each component of the k-point:
! 1. If k > 1.5 + dlt_BZ, subtract 2.0 and set RWtmp = -2
! 2. If k > 0.5 + dlt_BZ, subtract 1.0 and set RWtmp = -1
! 3. If k <= -1.5 + dlt_BZ, add 2.0 and set RWtmp = 2
! 4. If k <= -0.5 + dlt_BZ, add 1.0 and set RWtmp = 1
!
! Notes:
! - The tolerance dlt_BZ is used to handle numerical errors
! - The rewind vectors are used to track the transformation
! - The k-point is rewound to the range [-0.5, 0.5] in each direction
subroutine kcheck(ktmp,RWtmp) 
  implicit none 
  real(8),intent(inout)::ktmp(3)
  integer,intent(out)::RWtmp(3) 
  real(8),parameter::dlt_BZ=1.0d-6 
  if(ktmp(1)>1.50d0+dlt_BZ)then 
   ktmp(1)=ktmp(1)-2.0d0
   RWtmp(1)=-2
  endif 
  if(ktmp(1)>0.50d0+dlt_BZ)then 
   ktmp(1)=ktmp(1)-1.0d0
   RWtmp(1)=-1
  endif 
  if(ktmp(1)<=-1.50d0+dlt_BZ)then 
   ktmp(1)=ktmp(1)+2.0d0
   RWtmp(1)=2 
  endif 
  if(ktmp(1)<=-0.50d0+dlt_BZ)then 
   ktmp(1)=ktmp(1)+1.0d0
   RWtmp(1)=1
  endif 
  !
  if(ktmp(2)>1.50d0+dlt_BZ)then 
   ktmp(2)=ktmp(2)-2.0d0
   RWtmp(2)=-2 
  endif 
  if(ktmp(2)>0.50d0+dlt_BZ)then 
   ktmp(2)=ktmp(2)-1.0d0
   RWtmp(2)=-1
  endif 
  if(ktmp(2)<=-1.50d0+dlt_BZ)then 
   ktmp(2)=ktmp(2)+2.0d0
   RWtmp(2)=2 
  endif 
  if(ktmp(2)<=-0.50d0+dlt_BZ)then 
   ktmp(2)=ktmp(2)+1.0d0
   RWtmp(2)=1
  endif 
  !
  if(ktmp(3)>1.50d0+dlt_BZ)then 
   ktmp(3)=ktmp(3)-2.0d0
   RWtmp(3)=-2 
  endif 
  if(ktmp(3)>0.50d0+dlt_BZ)then 
   ktmp(3)=ktmp(3)-1.0d0
   RWtmp(3)=-1
  endif 
  if(ktmp(3)<=-1.50d0+dlt_BZ)then 
   ktmp(3)=ktmp(3)+2.0d0
   RWtmp(3)=2 
  endif 
  if(ktmp(3)<=-0.50d0+dlt_BZ)then 
   ktmp(3)=ktmp(3)+1.0d0
   RWtmp(3)=1
  endif 
  return
end 
!
! Check and rewind k-point coordinates for time-reversal symmetry
!
! This subroutine is similar to kcheck but uses a different tolerance
! for handling time-reversal symmetry. The main difference is in the
! comparison operators (>= instead of >) and the tolerance sign.
!
! Parameters:
! - ktmp(3): Input/output k-point coordinates (in reduced coordinates)
! - RWtmp(3): Output rewind vectors indicating how many times the k-point
!             was rewound along each direction
!
! Algorithm:
! For each component of the k-point:
! 1. If k >= 1.5 + dlt_BZ, subtract 2.0 and set RWtmp = -2
! 2. If k >= 0.5 + dlt_BZ, subtract 1.0 and set RWtmp = -1
! 3. If k < -1.5 + dlt_BZ, add 2.0 and set RWtmp = 2
! 4. If k < -0.5 + dlt_BZ, add 1.0 and set RWtmp = 1
!
! Notes:
! - The tolerance dlt_BZ is negative (-1.0d-6) to handle time-reversal symmetry
! - The rewind vectors are used to track the transformation
! - The k-point is rewound to the range [-0.5, 0.5] in each direction
subroutine kcheck_trs(ktmp,RWtmp) 
  implicit none 
  real(8),intent(inout)::ktmp(3)
  integer,intent(out)::RWtmp(3) 
  real(8),parameter::dlt_BZ=-1.0d-6 
  if(ktmp(1)>=1.50d0+dlt_BZ)then 
   ktmp(1)=ktmp(1)-2.0d0
   RWtmp(1)=-2
  endif 
  if(ktmp(1)>=0.50d0+dlt_BZ)then 
   ktmp(1)=ktmp(1)-1.0d0
   RWtmp(1)=-1
  endif 
  if(ktmp(1)<-1.50d0+dlt_BZ)then 
   ktmp(1)=ktmp(1)+2.0d0
   RWtmp(1)=2 
  endif 
  if(ktmp(1)<-0.50d0+dlt_BZ)then 
   ktmp(1)=ktmp(1)+1.0d0
   RWtmp(1)=1
  endif 
  !
  if(ktmp(2)>=1.50d0+dlt_BZ)then 
   ktmp(2)=ktmp(2)-2.0d0
   RWtmp(2)=-2 
  endif 
  if(ktmp(2)>=0.50d0+dlt_BZ)then 
   ktmp(2)=ktmp(2)-1.0d0
   RWtmp(2)=-1
  endif 
  if(ktmp(2)<-1.50d0+dlt_BZ)then 
   ktmp(2)=ktmp(2)+2.0d0
   RWtmp(2)=2 
  endif 
  if(ktmp(2)<-0.50d0+dlt_BZ)then 
   ktmp(2)=ktmp(2)+1.0d0
   RWtmp(2)=1
  endif 
  !
  if(ktmp(3)>=1.50d0+dlt_BZ)then 
   ktmp(3)=ktmp(3)-2.0d0
   RWtmp(3)=-2 
  endif 
  if(ktmp(3)>=0.50d0+dlt_BZ)then 
   ktmp(3)=ktmp(3)-1.0d0
   RWtmp(3)=-1
  endif 
  if(ktmp(3)<-1.50d0+dlt_BZ)then 
   ktmp(3)=ktmp(3)+2.0d0
   RWtmp(3)=2 
  endif 
  if(ktmp(3)<-0.50d0+dlt_BZ)then 
   ktmp(3)=ktmp(3)+1.0d0
   RWtmp(3)=1
  endif 
  return
end 
!
! Generate G-vectors for a given k-point within energy cutoff
!
! This subroutine generates G-vectors for a given k-point that satisfy
! the energy cutoff condition. The G-vectors are generated by adding
! integer multiples of reciprocal lattice vectors to the k-point.
!
! Parameters:
! - NTG: Maximum number of G-vectors
! - b1, b2, b3: Reciprocal lattice vectors
! - Gcut: Energy cutoff for G-vectors
! - q1, q2, q3: K-point coordinates (in reduced coordinates)
! - KG0: Output G-vectors (3,NTG)
! - NG: Output number of G-vectors
!
! Algorithm:
! 1. Loop over integer multiples of reciprocal lattice vectors
! 2. Calculate G-vector in Cartesian coordinates
! 3. Check if G-vector satisfies energy cutoff
! 4. Store G-vector if it does
!
! Notes:
! - G-vectors are generated in a box of size NGL1×NGL2×NGL3
! - The energy cutoff condition is |k+G|² ≤ Gcut
! - G-vectors are stored in reduced coordinates
! - The maximum number of G-vectors is limited by NTG
subroutine make_KG0(NTG,b1,b2,b3,Gcut,q1,q2,q3,KG0,NG) 
  implicit none 
  integer,intent(in)::NTG
  real(8),intent(in)::b1(3),b2(3),b3(3) 
  real(8),intent(in)::Gcut,q1,q2,q3        
  integer,intent(out)::KG0(3,NTG),NG 
  integer::igL,igL1,igL2,igL3
  real(8)::qgL(3),qgL2  
  integer,parameter::NGL1=150!100
  integer,parameter::NGL2=150!100 
  integer,parameter::NGL3=150!100
  igL=0
  ! Loop over integer multiples of reciprocal lattice vectors
  do igL1=-NGL1,NGL1 
   do igL2=-NGL2,NGL2 
    do igL3=-NGL3,NGL3 
     ! Calculate G-vector in Cartesian coordinates
     qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)    
     qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
     ! Check energy cutoff condition
     if(qgL2<=Gcut)then 
      igL=igL+1 
      KG0(1,igL)=igL1
      KG0(2,igL)=igL2
      KG0(3,igL)=igL3 
     endif  
    enddo 
   enddo 
  enddo 
  NG=igL 
  RETURN 
END 
!
! Estimate k-point sampling parameters from k-point coordinates
!
! This subroutine estimates the k-point sampling parameters (nkb1, nkb2, nkb3)
! by finding the smallest non-zero component of k-points along each direction.
!
! Parameters:
! - N: Number of k-points
! - SK: K-point coordinates (3,N)
! - nkb1, nkb2, nkb3: Output sampling parameters
!
! Algorithm:
! For each direction (1,2,3):
! 1. Find the smallest non-zero component of k-points
! 2. Calculate sampling parameter as the reciprocal of the smallest component
! 3. Round to nearest integer
!
! Notes:
! - Components smaller than 1.0d-7 are considered zero
! - The sampling parameters represent the number of points along each
!   reciprocal lattice vector
! - The total number of k-points should be nkb1×nkb2×nkb3
subroutine est_nkbi(N,SK,nkb1,nkb2,nkb3)  
  implicit none 
  integer::N,nkb1,nkb2,nkb3,NTK  
  real(8)::SK(3,N) 
  integer::i 
  real(8)::x 
  ! Estimate nkb1
  x=1.0d0 
  do i=1,N
   if(abs(SK(1,i))<1.0d-7)cycle 
   if(abs(SK(1,i))<x)then 
    x=abs(SK(1,i))  
   endif 
  enddo    
  nkb1=nint(1.0d0/x)  
  ! Estimate nkb2
  x=1.0d0 
  do i=1,N
   if(abs(SK(2,i))<1.0d-7)cycle 
   if(abs(SK(2,i))<x)then 
    x=abs(SK(2,i))  
   endif 
  enddo    
  nkb2=nint(1.0d0/x)  
  ! Estimate nkb3
  x=1.0d0 
  do i=1,N
   if(abs(SK(3,i))<1.0d-7)cycle 
   if(abs(SK(3,i))<x)then 
    x=abs(SK(3,i))  
   endif 
  enddo    
  nkb3=nint(1.0d0/x)  
  !
  NTK=nkb1*nkb2*nkb3 
  !
  return 
end 
!
! Estimate total number of k-points from irreducible k-points
!
! This subroutine estimates the total number of k-points (NTK) by applying
! symmetry operations and time-reversal symmetry to irreducible k-points.
!
! Parameters:
! - Nk_irr: Number of irreducible k-points
! - Nsymq: Number of symmetry operations
! - SKI: Irreducible k-points (3,Nk_irr)
! - rg: Rotation matrices (3,3,Nsymq)
! - NTK: Output total number of k-points
!
! Algorithm:
! 1. Allocate temporary array for k-points
! 2. For each irreducible k-point:
!    a. Apply all symmetry operations
!    b. Apply time-reversal symmetry
!    c. Check for duplicates and rewind if necessary
! 3. Count unique k-points
!
! Notes:
! - The maximum number of k-points is estimated as Nk_irr×Nsymq×2
! - Duplicate k-points are removed
! - K-points are rewound to the first Brillouin zone
subroutine est_NTK(Nk_irr,Nsymq,SKI,rg,NTK)
  implicit none 
  integer::Nk_irr,Nsymq,N  
  real(8)::SKI(3,Nk_irr) 
  integer::rg(3,3,Nsymq) 
  real(8),allocatable::SK0(:,:)!SK0(3,N) 
  real(8)::ktmp(3)
  integer::RWtmp(3)
  integer::jk,ik,iop,iik 
  integer::NTK 
  N=Nk_irr*Nsymq*2
  !
  !SK0
  !
  allocate(SK0(3,N));SK0(:,:)=0.0d0
  jk=0
  do ik=1,Nk_irr
   do iop=1,Nsymq
    ktmp(:)=0.0d0; RWtmp(:)=0  
    ktmp(1)=dble(rg(1,1,iop))*SKI(1,ik)+dble(rg(1,2,iop))*SKI(2,ik)+dble(rg(1,3,iop))*SKI(3,ik)
    ktmp(2)=dble(rg(2,1,iop))*SKI(1,ik)+dble(rg(2,2,iop))*SKI(2,ik)+dble(rg(2,3,iop))*SKI(3,ik)
    ktmp(3)=dble(rg(3,1,iop))*SKI(1,ik)+dble(rg(3,2,iop))*SKI(2,ik)+dble(rg(3,3,iop))*SKI(3,ik)
    call kcheck(ktmp(1),RWtmp(1))!rewind check 
    do iik=1,jk
     if(abs(SK0(1,iik)-ktmp(1))<1.0d-4.and.abs(SK0(2,iik)-ktmp(2))<1.0d-4.and.abs(SK0(3,iik)-ktmp(3))<1.0d-4) goto 1000
    enddo!iik
    jk=jk+1
    SK0(:,jk)=ktmp(:)
1000 ktmp(:)=0.0d0; RWtmp(:)=0  
    ktmp(1)=dble(rg(1,1,iop))*SKI(1,ik)+dble(rg(1,2,iop))*SKI(2,ik)+dble(rg(1,3,iop))*SKI(3,ik)
    ktmp(2)=dble(rg(2,1,iop))*SKI(1,ik)+dble(rg(2,2,iop))*SKI(2,ik)+dble(rg(2,3,iop))*SKI(3,ik)
    ktmp(3)=dble(rg(3,1,iop))*SKI(1,ik)+dble(rg(3,2,iop))*SKI(2,ik)+dble(rg(3,3,iop))*SKI(3,ik)
    call kcheck_trs(ktmp(1),RWtmp(1))!rewind check modified 20170316  
    do iik=1,jk
     if(abs(SK0(1,iik)-(-ktmp(1)))<1.0d-4.and.abs(SK0(2,iik)-(-ktmp(2)))<1.0d-4.and.abs(SK0(3,iik)-(-ktmp(3)))<1.0d-4) goto 2000
    enddo!iik
    jk=jk+1
    SK0(:,jk)=-ktmp(:) 
2000 enddo!iop 
  enddo!ik 
  NTK=jk 
  if(NTK>N)then 
   write(iunit_log,*)'Estimated NTK is too large; stop' 
   write(iunit_log,*)'NTK, N=',NTK, N 
   stop
  endif 
  return
end subroutine 
!
! Calculate lattice parameters from lattice vectors
!
! This subroutine calculates the lattice parameters (a, b, c, α, β, γ)
! from the lattice vectors (a1, a2, a3).
!
! Parameters:
! - a1, a2, a3: Lattice vectors
! - a, b, c: Output lattice constants
! - alp, bet, gmm: Output angles (in degrees)
!
! Algorithm:
! 1. Calculate lattice constants as magnitudes of lattice vectors
! 2. Calculate angles using dot products
! 3. Convert angles from radians to degrees
!
! Notes:
! - Lattice constants are in the same units as lattice vectors
! - Angles are in degrees
! - The angles are defined as:
!   α = angle between b and c
!   β = angle between c and a
!   γ = angle between a and b
subroutine est_latparam(a1,a2,a3,a,b,c,alp,bet,gmm)   
  implicit none 
  real(8)::a1(3),a2(3),a3(3) 
  real(8)::a,b,c,alp,bet,gmm  
  real(8),parameter::pi=DACOS(-1.0d0)
  ! Calculate lattice constants
  a=0.0d0 
  b=0.0d0 
  c=0.0d0 
  a=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2) 
  b=dsqrt(a2(1)**2+a2(2)**2+a2(3)**2) 
  c=dsqrt(a3(1)**2+a3(2)**2+a3(3)**2) 
  ! Calculate angles
  alp=(a2(1)*a3(1)+a2(2)*a3(2)+a2(3)*a3(3))/b/c 
  bet=(a3(1)*a1(1)+a3(2)*a1(2)+a3(3)*a1(3))/c/a
  gmm=(a1(1)*a2(1)+a1(2)*a2(2)+a1(3)*a2(3))/a/b 
  ! Convert angles to degrees
  alp=dacos(alp)*180.0d0/pi  
  bet=dacos(bet)*180.0d0/pi  
  gmm=dacos(gmm)*180.0d0/pi  
  return 
end 
!
end module 
