PROGRAM generate_sample_mk 
  use m_rd_dat_wfn
  implicit none

  integer::ik,ix,ig

  ! logging unit number
  iunit_log = 150
  open(iunit_log, FILE="LOG.mk")
  !read dir-wfn 
  call rd_dat_symmetry 
  call rd_dat_sample_k 
  call rd_dat_lattice 
  call rd_dat_nkm 
  call rd_dat_kg 
  close(iunit_log)

  OPEN(100,FILE="dat.sample_mk") 
  REWIND(100)
  write(100,'(i10)')NTK 
  do ik=1,NTK 
    write(100,'(3f15.10)')(SK0(ix,ik),ix=1,3) 
  enddo
  close(100)

  open(101, FILE="dat.kg_respack")
  REWIND(101)
  do ik=1,NTK
    write(101,*) NG0(ik)
    do ig=1, NG0(ik)
      write(101,*) KG0(:,ig,ik)
    end do
  end do
  close(101)
  stop
end
