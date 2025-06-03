!===============================================================================
! Utility functions for numerical computations
!===============================================================================
! This module provides functions for matrix operations and numerical utilities:
! 1. Matrix inversion (real and complex)
! 2. Matrix diagonalization
! 3. Vector operations
! 4. Number factorization utilities
!
! Dependencies:
! - LAPACK: Required for matrix operations (dgetrf, dgetri, zgetrf, zgetri, zheevd)
!===============================================================================

!-------------------------------------------------------------------------------
! Invert a real matrix using LU decomposition
!-------------------------------------------------------------------------------
! Parameters:
!   nm  : Size of the matrix
!   mat : Input/output matrix (nm, nm)
!
! Algorithm:
!   1. LU decomposition using dgetrf
!   2. Matrix inversion using dgetri
!
! Error handling:
!   - If info /= 0, prints error message and stops
!-------------------------------------------------------------------------------
subroutine invmat(nm,mat)
  implicit none 
  integer,intent(in)::nm
  real(8),intent(inout)::mat(nm,nm)
  integer::ipiv(nm)        ! Pivot indices for LU decomposition
  integer::Lwork           ! Size of work array
  real(8),allocatable::work(:)  ! Work array for dgetri
  integer::info            ! Error flag
  
  ! Set work array size (10*nm is a safe choice for most cases)
  Lwork = 10*nm
  allocate (work(Lwork))
  info = 0
  
  ! Perform LU decomposition
  call dgetrf(nm,nm,mat,nm,ipiv,info)
  
  ! Compute inverse using LU decomposition
  call dgetri(nm,mat,nm,ipiv,work,Lwork,info)
  
  ! Check for errors
  if(info /= 0) then
    write(6,*) 'info (subrouitine inv):' , info
    stop
  end if 
  
  deallocate(work)
  return 
end subroutine

!-------------------------------------------------------------------------------
! Invert a complex matrix using LU decomposition
!-------------------------------------------------------------------------------
! Parameters:
!   nm  : Size of the matrix
!   mat : Input/output complex matrix (nm, nm)
!
! Algorithm:
!   1. LU decomposition using zgetrf
!   2. Matrix inversion using zgetri
!
! Error handling:
!   - If info /= 0, prints error message and stops
!-------------------------------------------------------------------------------
subroutine invmat_complex(nm,mat)
  implicit none 
  integer,intent(in)::nm
  complex(8),intent(inout)::mat(nm,nm)
  integer::ipiv(nm)        ! Pivot indices for LU decomposition
  integer::Lwork           ! Size of work array
  complex(8),allocatable::work(:)  ! Work array for zgetri
  integer::info            ! Error flag
  
  ! Set work array size
  Lwork = 10*nm
  allocate (work(Lwork))
  info = 0
  
  ! Perform LU decomposition
  call zgetrf(nm,nm,mat,nm,ipiv,info)
  
  ! Compute inverse using LU decomposition
  call zgetri(nm,mat,nm,ipiv,work,Lwork,info)
  
  ! Check for errors
  if(info /= 0) then
    write(6,*) 'info (subrouitine inv):' , info
    stop
  end if 
  
  deallocate(work)
  return 
end subroutine

!-------------------------------------------------------------------------------
! Diagonalize a complex Hermitian matrix
!-------------------------------------------------------------------------------
! Parameters:
!   nm  : Size of the matrix
!   mat : Input/output complex Hermitian matrix (nm, nm)
!   eig : Output eigenvalues (nm)
!
! Algorithm:
!   Uses LAPACK's zheevd for Hermitian matrix diagonalization
!
! Error handling:
!   - If ind /= 0, prints error message and stops
!-------------------------------------------------------------------------------
subroutine diagV(nm,mat,eig)
  implicit none 
  integer,intent(in)::nm
  complex(8),intent(inout)::mat(nm,nm)
  real(8),intent(out)::eig(nm)
  integer::LWORK,LRWORK,LIWORK  ! Sizes of work arrays
  integer,allocatable::iwork_zheevd(:)  ! Integer work array
  real(8),allocatable::rwork_zheevd(:)  ! Real work array
  complex(8),allocatable::work_zheevd(:)  ! Complex work array
  integer::ind                 ! Error flag
  real(8)::eps                 ! Machine epsilon
  
  ! Set work array sizes
  LWORK= 2*nm+nm**2
  LRWORK=1+12*nm+3*nm**2
  LIWORK=3+10*nm 
  
  ! Allocate and initialize work arrays
  allocate(work_zheevd(LWORK));work_zheevd(:)=0.0d0
  allocate(rwork_zheevd(LRWORK));rwork_zheevd(:)=0.0d0
  allocate(iwork_zheevd(LIWORK));iwork_zheevd(:)=0
  
  eps=1.0d-18
  ind=0                 
  
  ! Perform diagonalization
  call zheevd("V","U",nm,mat,nm,eig,work_zheevd,LWORK,rwork_zheevd,LRWORK,iwork_zheevd,LIWORK,ind)
  
  ! Check for errors
  if(ind/=0)then 
   write(6,*)'ind=',ind 
   stop
  endif 
  
  deallocate(work_zheevd,rwork_zheevd,iwork_zheevd) 
  return 
end subroutine

!-------------------------------------------------------------------------------
! Calculate outer product of two vectors
!-------------------------------------------------------------------------------
! Parameters:
!   vec_x : First input vector (3)
!   vec_y : Second input vector (3)
!   vec_z : Output vector (3)
!
! Algorithm:
!   vec_z = vec_x × vec_y (cross product)
!-------------------------------------------------------------------------------
subroutine OUTER_PRODUCT(vec_x,vec_y,vec_z)
  implicit none 
  real(8)::vec_x(3),vec_y(3),vec_z(3) 
  vec_z(1)=vec_x(2)*vec_y(3)-vec_x(3)*vec_y(2)
  vec_z(2)=vec_x(3)*vec_y(1)-vec_x(1)*vec_y(3) 
  vec_z(3)=vec_x(1)*vec_y(2)-vec_x(2)*vec_y(1)
  return
end subroutine 

!-------------------------------------------------------------------------------
! Find the next number that can be factorized into powers of 2, 3, and 5
!-------------------------------------------------------------------------------
! Parameters:
!   inr : Input number
!
! Returns:
!   Next number that can be factorized into powers of 2, 3, and 5
!
! Algorithm:
!   Incrementally checks numbers until finding one that can be factorized
!   into powers of 2, 3, and 5 only
!-------------------------------------------------------------------------------
function algn235(inr)
  implicit none
  integer::algn235 
  integer,intent(in)::inr
  integer:: nr,m2,m3,m5,info
  nr=inr
  call fctck(nr,m2,m3,m5,info)
  do while (info .eq. 1)
     nr = nr + 1
     call fctck(nr,m2,m3,m5,info)
  end do
  algn235 = nr 
  return
end function algn235

!-------------------------------------------------------------------------------
! Factorize a number into powers of 2, 3, and 5
!-------------------------------------------------------------------------------
! Parameters:
!   n    : Input number to factorize
!   m2   : Output power of 2
!   m3   : Output power of 3
!   m5   : Output power of 5
!   info : Output status (0: success, 1: contains other factors)
!
! Algorithm:
!   Repeatedly divides by 2, 3, and 5 until either:
!   1. The number becomes 1 (success)
!   2. A different prime factor is found (failure)
!-------------------------------------------------------------------------------
subroutine fctck(n,m2,m3,m5,info)
  implicit none
  integer,intent(in):: n
  integer,intent(out):: m2, m3, m5, info
  integer:: i
  i=n
  m2 = 0
  m3 = 0
  m5 = 0
  info = 0
  do while (i .ne. 1)
     if (mod(i,2) .eq. 0) then
        m2 = m2 + 1
        i = i / 2
     else if (mod(i,3) .eq. 0) then
        m3 = m3 + 1
        i = i / 3
     else if (mod(i,5) .eq. 0) then
        m5 = m5 + 1
        i = i / 5
     else
        info = 1
        exit
     end if
  end do
  return
end subroutine fctck
