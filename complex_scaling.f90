program gridbasis
  implicit none
  integer, parameter :: nb=2000,nx=nb,ne=nb+nb**2,lwork=3*nb-1
  real *8 pi,xo,x(nx),xst,xend,dx,rwork(2*nb),alpha
  complex*16 ho(nb,nb),t(nb,nb),vv(nx),zi,eigval(nb),work(lwork),vl(nb,nb),to(nb,nb),&
       vr(nb,nb),vi(nx),vo(nb,nb),h(nb,nb),vii(nb,nb),norm(nb,nb),y(ne),emat(nb),psi(nb,nb),cmat(nb,nb),psi1(nb,nb)
  real*8 e(nb),v(nb,nb),c(nb,nb),lambda,hmat(nb,nb),dl,lambda1,theta
  complex*16 vl1(nb,nb),vr1(nb,nb),eigval1(nb),vmat(nb,nb),v2i(nb,nb),xx(nx),v5(nb,nb)
  integer m,n,k,i,j,info,ii,itheta,path,ineta,ialpha,ix,jj
  real *8:: lamda,xl,xr,mheu,neta
  complex*16 :: cap(nb,nb),norm1(nb,nb),eigen_check(nb,nb),&
       eigen_checksq(nb,nb),hsquare(nb,nb),eigen_check2(nb,nb)
  complex*16 :: standev(nb,nb),uni1(nb,nb),uni2(nb,nb)
  complex*16 :: eigval_cmplx(nb),vr_cmplx(nb,nb),vl_cmplx(nb,nb),h_cmplx(nb,nb),wfn(nb,nb)
  integer,parameter :: neta_end=2000,diff=100
  real*8 :: x_grid=40.d0
   integer :: r1=22,r2=21,r3=23 !300
  integer :: xbox
  pi=dacos(-1.d0) 
  zi =(0.d0,1.d0)
  alpha =0.05d0
  lamda=1.5d0
  path=8
  mheu=1.d0
  xo=0.d0
  xr=+20.d0
  xl=-xr
  theta=0.0d0
  xbox=xo
  mheu=1.d0!/2.0d0
  
   do ix=1,1
      xend = x_grid+(ix-1)*5
      xst=-xend
      dx=(xend-xst)/nx
   
  
  do m=1,nb
     do n=1,nb
        if ( m .eq. n ) then
           t(m,n)=(pi**2/(3.d0*2*dx**2))
        end if
        if (  m .ne. n) then
           t(m,n)=((-1)**(m-n)/(2.d0*dx**2))*2.d0/(m-n)**2
        end if
        ! write(*,*) m,n,v(m,n)
     end do
  end do
t=t!/mheu
  
to=t
do itheta=0,100
   theta=itheta*0.01
  to=t*exp(-2.d0*zi*theta)
  do i=1,nb
     do j=1,nb
        write(88,*) t(i,j),i,j
     enddo
  enddo


  do m=1,nx
     x(m)= xst+(m-1)*dx
     xx(m)=x(m)
     xx(m)=x(m)*exp(zi*theta)
     vv(m)=(1.0d0-(1.0d0/((cosh(xx(m)))**2)))*exp(-0.0250d0*xx(m)**2)
     to(m,m)=to(m,m)+vv(m)
     write(99,*)x(m),real(vv(m)),aimag(vv(m));call flush(10)
  end do
!stop
!stop
   vo=0.d0
  !do ialpha = neta_end,neta_end
   do ineta =0,0!neta_end
      write(*,*)itheta,"*********"
        neta=ineta*0.0000005!0.00001d0
  
              do i=1,nx
                   if (x(i) .gt. 0)then
                      if (x(i) .gt. xo) vo(i,i)= -zi*(x(i)-xo)**2
                         else
                      if (x(i) .lt. -xo) vo(i,i)= -zi*(x(i)+xo)**2   
                   endif

                 if (ineta .eq. 0 ) then
                   write(11,*) x(i),real(vo(i,i)),aimag(vo(i,i));call flush(11)
                 endif
             enddo
vo=0
   
   
     
    ho=to!+neta*vo

    
     call ZGEEV( 'V', 'V', nb,ho, nb, eigval1, VL1, nb, VR1, nb, WORK, LWORK, RWORK, INFO)
     call sorting2(eigval1,vr1,nb)

     do  i = 1,nb
       !  write(20,*) real(eigval1(i)),aimag(eigval1(i)),ineta,i,neta;call flush(20)
       ! write(27,*) real(eigval1(i)),aimag(eigval1(i)),ineta,i,neta;call flush(27)
        write(33,*) real(eigval1(i)),aimag(eigval1(i)),itheta,theta;call flush(33)
     end do
     
      norm=matmul(transpose(vr1),vr1)
     do i=1,nb
        vr1(:,i)=vr1(:,i)/sqrt(norm(i,i))
        
     enddo

  
     if (ineta .eq. 0) then
        do  i = 1,nb
           !write(20,*) real(eigval1(i)),aimag(eigval1(i)),ineta,i,neta;call flush(20)
           write(4000,*) x(i),real(vr1(i,r1))+real(eigval1(r1)),real(vr1(i,r2))&
           +real(eigval1(r1)),real(vr1(i,r3))+real(eigval1(r1));call flush(4000)
          end do
       endif

   !  endif
!enddo
  !stop

     h_cmplx=matmul(transpose(vr1),matmul(to,vr1))
     
    call ZGEEV( 'V', 'V', nb,h_cmplx, nb, eigval_cmplx, VL_cmplx, nb, VR_cmplx, nb, WORK, LWORK, RWORK, INFO)
    call sorting2(eigval_cmplx,vr_cmplx,nb)
     if (ineta .eq. neta_end) then
     do  i = 1,nb 
        write(28,*) real(eigval_cmplx(i)),aimag(eigval_cmplx(i));call flush(28)
     end do
  endif
! norm=0.00
  !norm=matmul(transpose(VR_cmplx),VR_cmplx)
   !  do i=1,nb
    !    VR_cmplx(:,i)=VR_cmplx(:,i)/sqrt(norm(i,i))
      !  
    ! enddo
!wfn=VR_cmplx
     
   ! if (ineta .eq. neta_end) then
     !   do  i = 1,nb
         !    write(29,*) real(eigval_cmplx(i)),aimag(eigval_cmplx(i));call flush(29)
           !  write(2000,*) x(i),real(VR_cmplx(i,28)),real(VR_cmplx(i,28));call flush(2000)
           !   write(2001,*) x(i),real(VR_cmplx(i,r1)),real(VR_cmplx(i,r2));call flush(2001)
            ! write(3000+ineta,*) x(i),real(VR_cmplx(i,28))+real(eigval_cmplx(28)),&
          !        real(wfn(i,r2))+real(eigval_cmplx(r2)),real(wfn(i,r3))+real(eigval_cmplx(r3));call flush(3000+ineta)
        !  end do
       !endif


       wfn=matmul(vr1,vr_cmplx)
     
     if (ineta .eq. neta_end) then
        do  i = 1,nb
             write(29,*) real(eigval_cmplx(i)),aimag(eigval_cmplx(i));call flush(29)
            ! write(2000,*) x(i),real(vr_cmplx(i,r1)),real(vr_cmplx(i,r2));call flush(2000)
             write(3000+ineta,*) x(i),real(wfn(i,r1))+real(eigval_cmplx(r1)),&
                  real(wfn(i,r2))+real(eigval_cmplx(r2)),real(wfn(i,r3))+real(eigval_cmplx(r3));call flush(3000+ineta)
     end do
  endif
enddo
  

enddo                  !!!!  ix loop

enddo

 
  end program













subroutine sorting2(a,b,n)
!integer, parameter :: dp = selected_real_kind( p=15, r=307 )
  integer n,i,j
  complex*16 a(n),temp,temp2(n),b(n,n)
  do i=1,n
   do  j=1,n-i
     if (real(a(j)) .gt. real(a(j+1)))then
        temp=a(j)   ;  temp2(:)=b(:,j)
        a(j)=a(j+1) ;  b(:,j)=b(:,j+1)
        a(j+1)=temp ;  b(:,j+1)=temp2(:)
     end if
  end do
  end do
! write(6,*)'values', a
end subroutine sorting2


Subroutine driver_dsyev(hmat,cmat,emat,n)

  complex*16 a
end Subroutine driver_dsyev
