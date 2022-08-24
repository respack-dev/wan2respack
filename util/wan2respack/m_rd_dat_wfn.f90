module m_rd_dat_wfn
  implicit none
  public::rd_dat_symmetry
  public::rd_dat_lattice 
  public::rd_dat_sample_k 
  public::rd_dat_nkm 
  public::rd_dat_kg 
  !sym(100)
  integer,public::nsymq,nnp
  integer,public,allocatable::rg(:,:,:)!3,3,nsymq)
  integer,public,allocatable::pg(:,:)!pg(3,nsymq)
  real(8),public,allocatable::rginv(:,:,:)!rginv(3,3,nsymq) 
  !bandcalc(117) 
  real(8),public::Ecut_for_psi 
  real(8),public::FermiEnergy  
  real(8),public::Etot
  !avec(105)
  real(8),public::a1(3),a2(3),a3(3)
  real(8),public::b1(3),b2(3),b3(3)
  real(8),public::VOLUME
  !sample-k(101)
  integer,public::Nk_irr,NTK 
  integer,public::nkb1!SAMPLING K POINTS ALONG b1 VECTOR
  integer,public::nkb2!SAMPLING K POINTS ALONG b2 VECTOR
  integer,public::nkb3!SAMPLING K POINTS ALONG b3 VECTOR
  integer,public::Na1!LATTICE TRANSLATIONL IN a1 20170327 
  integer,public::Na2!LATTICE TRANSLATIONL IN a2 20170327 
  integer,public::Na3!LATTICE TRANSLATIONL IN a3 20170327 
  real(8),public,allocatable::SKI(:,:)!SKI(3,Nk_irr)  
  real(8),public,allocatable::SK0(:,:)!SK0(3,NTK)  
  integer,public,allocatable::numirr(:)!numirr(NTK) 
  integer,public,allocatable::numrot(:)!numrot(NTK) 
  integer,public,allocatable::trs(:)!trs(NTK) 
  integer,public,allocatable::RW(:,:)!RW(3,NTK)
  integer,public,allocatable::numMK(:)!numMK(Nk_irr)!20180316  
  !eigenvalue(111) 
  integer,public::NTB 
  real(8),public,allocatable::E_EIGI(:,:)!E_EIGI(NTB,Nk_irr) 
  !nkm(132)  
  integer,public::NTG 
  integer,public,allocatable::NGI(:)!NG0(Nk_irr)
  !kg(104)
  integer,public::L1,L2,L3 
  integer,public,allocatable::KGI(:,:,:)!KGI(3,NTG,Nk_irr) 
  integer,public,allocatable::KG0(:,:,:)!KG0(3,NTG,NTK)
  integer,public,allocatable::NG0(:)!NG0(NTK)
  integer,public,allocatable::packing(:,:,:,:)!packing(-L1:L1,-L2:L2,-L3:L3,Nk_irr)
  !fft 
  integer,public::nwx2,nwy2,nwz2!,nfft1,nfft2,nfft3,Nl123 
  !wfn(102)
  integer,public::ncomp 
  !complex(8),public,allocatable::CIR(:,:,:)!CIR(NTG,NTB,Nk_irr) 
  complex(4),public,allocatable::CIR(:,:,:)!CIR(NTG,NTB,Nk_irr) 
  complex(8),allocatable::CIRtmp(:)!CIRtmp(NTG)
  !logging file unit
  integer,public::iunit_log
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
subroutine rd_dat_lattice 
  implicit none 
  real(8),parameter::tpi=2.0d0*dacos(-1.0d0)
  OPEN(105,FILE='./dir-wfn/dat.lattice') 
  REWIND(105)
  READ(105,*) a1(1),a1(2),a1(3)!a1 vector
  READ(105,*) a2(1),a2(2),a2(3)!a2 vector
  READ(105,*) a3(1),a3(2),a3(3)!a3 vector
  CLOSE(105)
  !--
  call OUTER_PRODUCT(a2(1),a3(1),b1(1))
  VOLUME=a1(1)*b1(1)+a1(2)*b1(2)+a1(3)*b1(3)
  b1(:)=b1(:)*tpi/VOLUME 
  call OUTER_PRODUCT(a3(1),a1(1),b2(1))
  b2(:)=b2(:)*tpi/VOLUME 
  call OUTER_PRODUCT(a1(1),a2(1),b3(1))
  b3(:)=b3(:)*tpi/VOLUME 
end subroutine
!
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
  call est_NTK(Nk_irr,nsymq,SKI(1,1),rg(1,1,1),NTK)
  write(iunit_log,*)'Estimated NTK=',NTK 
  !--
  !SK0,numirr,numrot,trs,RW
  allocate(SK0(3,NTK));SK0(:,:)=0.0d0
  allocate(numirr(NTK));numirr(:)=0
  allocate(numrot(NTK));numrot(:)=0
  allocate(trs(NTK));trs(:)=0
  allocate(RW(3,NTK));RW(:,:)=0
  allocate(numMK(Nk_irr));numMK=0!20180316 
  jk=0
  do ik=1,Nk_irr
   initial_flg=0!20180316 
   do iop=1,nsymq
    !sym
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
    numirr(jk)=ik
    numrot(jk)=iop
    trs(jk)=1
    RW(:,jk)=RWtmp(:)
    if(initial_flg.eq.0)then
     numMK(ik)=jk
     initial_flg=1 
    endif 
    !time-reversal
1000 ktmp(:)=0.0d0;RWtmp(:)=0  
    ktmp(1)=dble(rg(1,1,iop))*SKI(1,ik)+dble(rg(1,2,iop))*SKI(2,ik)+dble(rg(1,3,iop))*SKI(3,ik)
    ktmp(2)=dble(rg(2,1,iop))*SKI(1,ik)+dble(rg(2,2,iop))*SKI(2,ik)+dble(rg(2,3,iop))*SKI(3,ik)
    ktmp(3)=dble(rg(3,1,iop))*SKI(1,ik)+dble(rg(3,2,iop))*SKI(2,ik)+dble(rg(3,3,iop))*SKI(3,ik)
    call kcheck_trs(ktmp(1),RWtmp(1))!rewind check modified 20170316  
    do iik=1,jk
     if(abs(SK0(1,iik)-(-ktmp(1)))<1.0d-4.and.abs(SK0(2,iik)-(-ktmp(2)))<1.0d-4.and.abs(SK0(3,iik)-(-ktmp(3)))<1.0d-4) goto 2000
    enddo!iik
    jk=jk+1
    SK0(:,jk)=-ktmp(:) 
    numirr(jk)=ik
    numrot(jk)=iop
    trs(jk)=-1
    RW(:,jk)=RWtmp(:) 
2000 enddo!iop  
  enddo!ik
  call est_nkbi(NTK,SK0(1,1),nkb1,nkb2,nkb3)  
  Na1=nkb1/2; Na2=nkb2/2; Na3=nkb3/2
  !--
  if(NTK/=jk) then 
   write(iunit_log,*)'ERROR;STOP;NTK should be jk'   
   write(iunit_log,*)'NTK=',NTK,'jk=',jk;STOP
  endif 
end subroutine
!--
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
  allocate(LKGI(NTG,Nk_irr));LKGI=0.0d0 
  do ik=1,Nk_irr 
   do ig=1,NGI(ik) 
    ktmp(1)=(SKI(1,ik)+dble(KGI(1,ig,ik)))*b1(1)+(SKI(2,ik)+dble(KGI(2,ig,ik)))*b2(1)+(SKI(3,ik)+dble(KGI(3,ig,ik)))*b3(1) 
    ktmp(2)=(SKI(1,ik)+dble(KGI(1,ig,ik)))*b1(2)+(SKI(2,ik)+dble(KGI(2,ig,ik)))*b2(2)+(SKI(3,ik)+dble(KGI(3,ig,ik)))*b3(2) 
    ktmp(3)=(SKI(1,ik)+dble(KGI(1,ig,ik)))*b1(3)+(SKI(2,ik)+dble(KGI(2,ig,ik)))*b2(3)+(SKI(3,ik)+dble(KGI(3,ig,ik)))*b3(3) 
    LKGI(ig,ik)=ktmp(1)**2+ktmp(2)**2+ktmp(3)**2
   enddo!ig 
   !write(6,*)maxval(LKGI(:,ik))
  enddo!ik 
  !write(6,*) 
  !write(6,*)maxval(LKGI(:,:))
  Ecut_for_psi=maxval(LKGI(:,:))+1.0d-8
  deallocate(LKGI) 
  !--
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
  !fft grid
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
  !NG0,KG0 
  allocate(NG0(NTK));NG0(:)=0
  allocate(KG0(3,NTG,NTK));KG0(:,:,:)=0 
  allocate(KGtmp(3,NTG));KGtmp(:,:)=0 
  do jk=1,NTK 
   if(trs(jk)==1) then 
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
  do igL1=-NGL1,NGL1 
   do igL2=-NGL2,NGL2 
    do igL3=-NGL3,NGL3 
     qgL(:)=(q1+dble(igL1))*b1(:)+(q2+dble(igL2))*b2(:)+(q3+dble(igL3))*b3(:)    
     qgL2=qgL(1)**2+qgL(2)**2+qgL(3)**2
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
subroutine est_nkbi(N,SK,nkb1,nkb2,nkb3)  
  implicit none 
  integer::N,nkb1,nkb2,nkb3,NTK  
  real(8)::SK(3,N) 
  integer::i 
  real(8)::x 
  x=1.0d0 
  do i=1,N
   if(abs(SK(1,i))<1.0d-7)cycle 
   if(abs(SK(1,i))<x)then 
    x=abs(SK(1,i))  
   endif 
  enddo    
  nkb1=nint(1.0d0/x)  
  x=1.0d0 
  do i=1,N
   if(abs(SK(2,i))<1.0d-7)cycle 
   if(abs(SK(2,i))<x)then 
    x=abs(SK(2,i))  
   endif 
  enddo    
  nkb2=nint(1.0d0/x)  
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
subroutine est_latparam(a1,a2,a3,a,b,c,alp,bet,gmm)   
  implicit none 
  real(8)::a1(3),a2(3),a3(3) 
  real(8)::a,b,c,alp,bet,gmm  
  real(8),parameter::pi=DACOS(-1.0d0)
  a=0.0d0 
  b=0.0d0 
  c=0.0d0 
  a=dsqrt(a1(1)**2+a1(2)**2+a1(3)**2) 
  b=dsqrt(a2(1)**2+a2(2)**2+a2(3)**2) 
  c=dsqrt(a3(1)**2+a3(2)**2+a3(3)**2) 
  alp=(a2(1)*a3(1)+a2(2)*a3(2)+a2(3)*a3(3))/b/c 
  bet=(a3(1)*a1(1)+a3(2)*a1(2)+a3(3)*a1(3))/c/a
  gmm=(a1(1)*a2(1)+a1(2)*a2(2)+a1(3)*a2(3))/a/b 
  alp=dacos(alp)*180.0d0/pi  
  bet=dacos(bet)*180.0d0/pi  
  gmm=dacos(gmm)*180.0d0/pi  
  return 
end 
!
end module 
