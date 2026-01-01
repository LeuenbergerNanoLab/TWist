


 Module Geometry2
    implicit none
    real(8), parameter     :: PI = 3.141592653589793238462643383279502884197169d0
    real(8)     :: a11(2),a12(2)                    ! primitive vectors for 1st layer
    real(8)     :: a21(2),a22(2)                    ! primitive vectors for 2nd layer
    real(8)     :: alat                             ! lattice constant
    integer     :: m,n
    integer     :: Ntheta                           ! number of primitive cells
    real(8)     :: theta,costh,sinth
    real(8)     :: T1(2)                            ! supercell lattice vectors
    real(8)     :: T2(2)
    real(8)     :: T1s(2)                           ! supercell lattice vectors
    real(8)     :: T2s(2)
    character(2), allocatable   :: Atoms1(:)        ! atom names for 1st layer
    character(2), allocatable   :: Atoms2(:)        ! atom names for 2nd layer
    real(8)     , allocatable   :: Pos1(:,:)        ! positions of atoms in 1st layer
    real(8)     , allocatable   :: Pos2(:,:)        ! positions of atoms in 2nd layer
    real(8)     , allocatable   :: PosX(:,:)        ! atom positions in bilayer
    character(2), allocatable   :: AtomsX(:)        ! atom names  
    real(8)     , allocatable   :: Pos(:,:)         ! temporary  
    real(8)                     :: z0               ! geometry parameter for MoS2
    real(8)                     :: dx               ! distance between MoS2 layers
 contains
 
 
 
 subroutine set_layers_a
 
  allocate(Atoms1(3*Ntheta))
  allocate(Atoms2(3*Ntheta))
  allocate(AtomsX(6*Ntheta))
  allocate(Pos1(3,3*Ntheta))
  allocate(Pos2(3,3*Ntheta))
  allocate(PosX(3,6*Ntheta))
  allocate(Pos(2,Ntheta))

  alat = 3.16d0                                   ! MoS2
  z0 = 1.585d0
  dx = 6.15d0

  a11(1) =  alat                                  ! set 1st layer
  a11(2) =  0.d0
  a12(1) = -alat*0.5d0
  a12(2) =  alat*dsqrt(3.d0)*0.5d0

  a21(1) = a11(1)*costh - a11(2)*sinth            ! set 2nd layer
  a21(2) = a11(1)*sinth + a11(2)*costh
  a22(1) = a12(1)*costh - a12(2)*sinth
  a22(2) = a12(1)*sinth + a12(2)*costh

  print *,'a11=',a11
  print *,'a12=',a12
  print *,'a21=',a21
  print *,'a22=',a22
 end subroutine set_layers_a
 
 
 
 subroutine calc_Ntheta
  Ntheta = m**2 + m*n + n**2
  print 1,Ntheta
1 format('Ntheta=',I7)  
 end subroutine calc_Ntheta
 
 
 subroutine calc_theta
  costh = dfloat(m**2+4*m*n+n**2)/dfloat(2*Ntheta)
  sinth = dsqrt(1.d0-costh**2)
  theta = dacos(costh)
  print 2,sinth**2+costh**2
  print 1,theta*180.d0/PI
1 format('theta=',F12.3)
2 format('sin2+cos2=',F22.17)
 end subroutine calc_theta
 
 
 
 subroutine calc_supercellT
  T1 =  m*a11 - n*a12
  T2 = n*a11 + (m+n)*a12
  print 1,T1,dsqrt(T1(1)**2+T1(2)**2)
  print 2,T2,dsqrt(T2(1)**2+T2(2)**2)
  
  T1s =  m*a21 - n*a22
  T2s = n*a21 + (m+n)*a22
  print 3,T1s,dsqrt(T1s(1)**2+T1s(2)**2)
  print 4,T2s,dsqrt(T2s(1)**2+T2s(2)**2)
1 format('T1=',2F20.16,'  abs =',F20.16)  
2 format('T2=',2F20.16,'  abs =',F20.16)  
3 format('T1s=',2F20.16,'  abs =',F20.16)  
4 format('T2s=',2F20.16,'  abs =',F20.16)  
 end subroutine calc_supercellT



 subroutine gen_layer1(At,a1,a2,T1a,T2a)                   ! generate atoms for supercell
  integer         :: Na
  real(8)         :: At(2)
  real(8)         :: a1(2)
  real(8)         :: a2(2)
  real(8)         :: T1a(2),T2a(2)
  integer         :: Ngen
  integer         :: i,j
  real(8)         :: R(2)
  real(8)         :: u,v
  real(8)         :: eps
  eps = 1.d-9
  Na = Ntheta
  Ngen = 0   
  Pos = 0.d0 
  do i = -Na,Na
   do j = -Na,Na
    R = At + i*a1 + j*a2
    call miller_indexes(R,T1a,T2a,u,v)
    if(-0.5d0 < u .and. u <= 0.5d0+eps .and. -0.5d0 < v .and. v <+ 0.5d0+eps) then
     print *,'u,v=',u,v
     Ngen = Ngen + 1                         ! we found atom within (T1,T2) supercell
     if(Ngen<=Ntheta) then
      Pos(1:2,Ngen) = R(1:2)
     endif
    endif
   enddo
  enddo
  if(Ngen /= Na) then
   print *,'ERROR in generation'
   print *,'Ngen =',Ngen
  endif 
 end subroutine gen_layer1



      subroutine miller_indexes(Qe,Aa,Ab,u,v)
       real(8)   :: Qe(2)
       real(8)   :: Aa(2),Ab(2)
       real(8)   :: u,v
       real(8)   :: M(2,2),Y(2)
       M(1,1:2) = Aa(1:2)
       M(2,1:2) = Ab(1:2)
       call solve2x2(M,Qe,Y)                                                  ! solve 2x2 equation M*Y=Qe
       u = Y(1)
       v = Y(2)
      end subroutine miller_indexes



     subroutine solve2x2(A,B,X)            ! solve 2x2 linear equation AX=B
      real(8)  :: A(2,2)
      real(8)  :: B(2)
      real(8)  :: X(2)
      real(8)  :: D,D1,D2
      real(8)  :: Aw(2,2)
      call det2x2(A,D)
      call make_dx(A,1,B,Aw)
      call det2x2(Aw,D1)
      call make_dx(A,2,B,Aw)
      call det2x2(Aw,D2)
      X(1) = D1/D
      X(2) = D2/D
     end subroutine solve2x2



     subroutine make_dx(A,j,B,Aw)           ! substitute column j in A
      real(8)    :: A(2,2),B(2)
      integer    :: j
      real(8)    :: Aw(2,2)
      integer    :: i
      Aw(1:2,1:2) = A(1:2,1:2)
      do i=1,2
       Aw(j,i) = B(i)
      enddo
     end subroutine make_dx
 
 
 
     subroutine det2x2(A,D)                 ! determinant of 2x2
      real(8)  :: A(2,2)
      real(8)  :: D
      D = A(1,1)*A(2,2) - A(1,2)*A(2,1)  
     end subroutine det2x2



      subroutine write_xsf(Period,Atm,Coord,Nat)
       integer             :: Nat
       character(2)        :: Atm(Nat)
       real(8)             :: Coord(3,Nat)
       real(8)             :: Period(3,3)               
       integer             :: i   !,j,k
       real(8)             :: shift(3)
       shift(1:3) = 0.5d0*(Period(1:3,1)+Period(1:3,2))
       open(unit=1,file='fileX.xsf')
        write(1,1)
        write(1,2)
        write(1,3) Period(1:3,1)
        write(1,3) Period(1:3,2)
        write(1,3) Period(1:3,3)
        write(1,4) Nat
        do i=1,Nat
         write(1,5) Atm(i),Coord(1:3,i) + shift(1:3) 
        enddo
       close(unit=1)
       print *,'writing file new.xsf'
1      format('CRYSTAL')
2      format('PRIMVEC')
3      format(3F15.10)
4      format('PRIMCOORD'/I5,'    1')
5      format(1x,A5,3F15.10)    !,2x,3F15.10)
6      format('ANIMSTEPS',I3)
7      format('ATOMS',I3)
8      format('BEGIN_BLOCK_DATAGRID_3D')
9      format('3D_PWSCF')
10     format('BEGIN_DATAGRID_3D_UNKNOWN')
11     format(3I12)
12     format(3F12.6)
13     format(6E14.6)
14     format('END_DATAGRID_3D')
15     format('END_BLOCK_DATAGRID_3D')
      end subroutine write_xsf


  subroutine fill_atoms_layer(a1,a2,Posr,Atomsr,d0,T1a,T2a)                   ! fill all Mo and S atoms for one layer
    real(8)         :: At(2)
    real(8)         :: a1(2),a2(2)
    real(8)         :: Posr(:,:)
    character(2)    :: Atomsr(:)
    real(8)         :: d0
    real(8)         :: T1a(2),T2a(2)
         
    At = 1.d0/3.d0*a1 + 2.d0/3.d0*a2                      ! Mo atoms
    print *,'generate Mo'
    call gen_layer1(At,a1,a2,T1a,T2a)
    Posr(1:2,1:Ntheta) = Pos(1:2,1:Ntheta)
    Posr(3,1:Ntheta) = d0
    Atomsr(1:Ntheta) = 'Mo'

    At = 2.d0/3.d0*a1 + 1.d0/3.d0*a2                      ! S atoms
    print *,'generate S1'
    call gen_layer1(At,a1,a2,T1a,T2a)
    Posr(1:2,Ntheta+1:2*Ntheta)   = Pos(1:2,1:Ntheta)
    Posr(3,Ntheta+1:2*Ntheta)     = d0 + z0
    Atomsr(Ntheta+1:2*Ntheta)     = 'S '

    At = 2.d0/3.d0*a1 + 1.d0/3.d0*a2                      ! S atoms
    print *,'generate S2'
    call gen_layer1(At,a1,a2,T1a,T2a)
    Posr(1:2,2*Ntheta+1:3*Ntheta) = Pos(1:2,1:Ntheta)
    Posr(3,2*Ntheta+1:3*Ntheta)   = d0 - z0
    Atomsr(2*Ntheta+1:3*Ntheta) = 'S '
  end subroutine fill_atoms_layer
  
  
  subroutine combine2
   PosX(1:3,1:3*Ntheta)          = Pos1(1:3,1:3*Ntheta)
   PosX(1:3,3*Ntheta+1:6*Ntheta) = Pos2(1:3,1:3*Ntheta)
   AtomsX(1:3*Ntheta)            = Atoms1(1:3*Ntheta)
   AtomsX(3*Ntheta+1:6*Ntheta)   = Atoms2(1:3*Ntheta)
  end subroutine combine2



  subroutine gen_bilayer(m2,n2)
    integer       :: m2,n2
    m = m2
    n = n2
    call calc_Ntheta 
    call calc_theta 
    call set_layers_a
    call calc_supercellT 
    call fill_atoms_layer(a11,a12,Pos1,Atoms1,3.d0,T1,T2)                  ! fill 1st layer 
    call fill_atoms_layer(a21,a22,Pos2,Atoms2,3.d0+dx,T1s,T2s)             ! fill 2nd layer 
    call combine2
  end subroutine gen_bilayer


 end Module Geometry2



   Program bilayer
    use Geometry2
    implicit none
    real(8)           :: Period(3,3)

 !   m = 2
 !   n = 1
    call gen_bilayer(2,1)
!    call calc_Ntheta 
!    call calc_theta 
!    call set_layers_a
!    call calc_supercellT 
!    call fill_atoms_layer(a11,a12,Pos1,Atoms1,3.d0,T1,T2)                  ! fill 1st layer 
!    call fill_atoms_layer(a21,a22,Pos2,Atoms2,3.d0+dx,T1s,T2s)             ! fill 2nd layer 
!    call combine2
    Period(1:2,1) = T1(1:2)
    Period(3,1) = 0.d0
    Period(1:2,2) = T2(1:2)
    Period(3,2) = 0.d0
    Period(1:3,3) = (/0.d0,0.d0,20.d0/)
    call write_xsf(Period,AtomsX,PosX,6*Ntheta) 
    deallocate(Atoms1,Atoms2,AtomsX)
    deallocate(Pos,Pos1,Pos2,PosX)   
   end Program bilayer



