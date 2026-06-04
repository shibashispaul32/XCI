!-----------------------------------------------------------
!       Mori et al, BioPhys 94 3684 2008
!-----------------------------------------------------------
implicit none
double precision,dimension(:),allocatable:: u1,u2,u1m,u2m,u3
double precision,dimension(:),allocatable:: p
double precision::dt,t,dx,x,l,u1_0,u2_0
double precision::pi,D_1,D_2
double precision::delta,gama,k_0,k,f
double precision::ks,ks_loc,ks_grad,ks_grad_rev
double precision::alpha1,beta1,alpha2,beta2
double precision::lap1,lap2
double precision::s,t0,t1,t2
double precision::power_stim
double precision::gama_on,gama_off
double precision::u1_max,u1_min,u1_sum
integer::i,j,k_time,n_time,n_x,n_stim,n_hill,n_wri
integer::n_new
common l,s,t0,t1,t2,power_stim,n_stim
!----------------------------------------
open(1,file='input.in',status='unknown')
!----------------------------------------
read(1,*)u1_0,u2_0                       !initial values
read(1,*)dx,n_x                          !parameters associated to space discreetization
read(1,*)dt,n_time,n_wri                 !parameters associated to time discreetization
read(1,*)D_1,D_2                         !diffusion coefficient
read(1,*)delta,k_0,k,n_hill              !System parameters
read(1,*)alpha1,beta1                    !System parameters
read(1,*)alpha2,beta2                    !System parameters
read(1,*)n_stim,power_stim               !number of stimuli
read(1,*)s,t0,t1,t2                      !stimulation parameters
read(1,*)gama_on,gama_off                !System parameters
pi=4.0d0*atan(1.0)
l=n_x*dx
!----------------------------------------
allocate(u1(n_x),u2(n_x))
allocate(u1m(n_x),u2m(n_x))
allocate(u3(n_x))
!----------------------------------------
! Assign initial values and stimulations
!----------------------------------------
u1=u1_0;u2=u2_0;t=0.0d0
u3=1
u3(1:n_x/10)=0
u3((9*n_x)/10:n_x)=0
!!----------------------------------------
!!  Assign zero flux boundary condition
!!----------------------------------------
!u1(1)=u1(2);u1(n_x)=u1(n_x-1)
!u2(1)=u2(2);u2(n_x)=u2(n_x-1)
!----------------------------------------
!  Assign periodic boundary condition
!----------------------------------------
u1(1)=u1(n_x-1);u1(n_x)=u1(2)
u2(1)=u2(n_x-1);u2(n_x)=u2(2)
!----------------------------------------
do i=1,n_x
   write(200,*)i*dx,u1(i),u2(i)
end do
!----------------------------------------
!      PDE integration starts here
!----------------------------------------
do k_time=0,n_time                 !time loop starts here
t=k_time*dt
!----------------------------------------
u1m=u1;u2m=u2
do i=2,n_x-1
x=i*dx
!----------------------------------------
! Laplacians
!----------------------------------------
  lap1=D_1*((u1m(i+1)+u1m(i-1)-(2*u1m(i)))/(2*dx**2))
  lap2=D_2*((u2m(i+1)+u2m(i-1)-(2*u2m(i)))/(2*dx**2))
!----------------------------------------
!  Kinetic terms
!----------------------------------------
  if (u3(i).eq.0)then
          gama=gama_off
  else
          gama=gama_on
  end if
  f=(u2(i)*(k_0+((gama*u1m(i)**n_hill)/(k**n_hill+u1m(i)**n_hill))))-(delta*u1m(i))
! if(t.le.200)ks=ks_grad(t,x)
! if(t.gt.200)ks=ks_grad_rev(t,x)
  ks=ks_loc(t,x)
  u1(i)=u1m(i)+(f+beta1-(alpha1*u1m(i))+lap1+(ks*u2(i)))*dt
  u2(i)=u2m(i)+(-f+beta2-(alpha2*u2m(i))+lap2-(ks*u2(i)))*dt
!----------------------------------------
end do
!!----------------------------------------
!!  Assign zero flux boundary condition
!!----------------------------------------
!u1(1)=u1(2);u1(n_x)=u1(n_x-1)
!u2(1)=u2(2);u2(n_x)=u2(n_x-1)
!----------------------------------------
!  Assign periodic boundary condition
!----------------------------------------
u1(1)=u1(n_x-1);u1(n_x)=u1(2)
u2(1)=u2(n_x-1);u2(n_x)=u2(2)
!----------------------------------------
! File writing
!----------------------------------------
if(mod(k_time,n_wri).eq.0)then
        do i=1,n_x
!        write((200+(t/t_wri)),*)(i*dx),u1(i),u2(i)
        write(600,*)t,(i*dx),u1(i)
        end do
        write(600,*)
endif
!if (t.gt.300.and.t.lt.350)then
!        u1_max=MAXVAL(u1)
!        u1_min=MINVAL(u1)
!        n_new=20
!       u1=u1_min
!       u1(n_x/2-n_new:n_x/2+n_new)=u1_max
!       u2=(((u1_0+u2_0)*n_x)-sum(u1))/n_x
!end if
!----------------------------------------
end do             !time loop ends here
end                !main program ends here






!----------------------------------------
function ks_grad(t,x)
         double precision::t,x,pi
         double precision::s,t0,t1,t2
         double precision::ks_grad
         pi=4.0d0*atan(1.0)
         s=0.07d0;t0=0.0d0;t1=20.0d0;t2=25.0d0
         ks_grad=0.0d0
         if(x.le.10)then
           if(t.le.t1)ks_grad=s*(10.0d0-x)
           if(t.gt.t1.and.t.le.t2)ks_grad=s*(1-((t-t1)/(t2-t1)))*(10-x)
         end if
         return
         end                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function ks_grad_rev(t,x)
         double precision::t,x,pi
         double precision::s,t0,t1,t2
         double precision::ks_grad_rev
         pi=4.0d0*atan(1.0)
         s=0.07d0;t0=200.0d0;t1=220.0d0;t2=225.0d0
         ks_grad_rev=0.0d0
         if(x.le.10)then
           if(t.gt.t0.and.t.le.t1)ks_grad_rev=s*x
           if(t.gt.t1.and.t.le.t2)ks_grad=s*x*(1-((t-t1)/(t2-t1)))
         end if
         return
         end
!----------------------------------------
!function gaussian_white(mu,sigma)
!        double precision::mu,sigma
!        double precision::a,b,pi
!        pi=4.0d0*atan(1.0)
!        call random_seed
!        call random_number (harvest=a)
!        call random_number (harvest=b)
!        gaussian_white=sigma*((dsqrt(-2.0d0*dlog(a))*(dcos(2*pi*b)))+mu)
!        return
!        end
!----------------------------------------
function ks_loc(t,x)
         double precision::t,x,pi
         double precision::s,t0,t1,t2
         double precision::ks_loc,l,power_stim
         integer :: n_stim
         common l,s,t0,t1,t2,power_stim,n_stim
         pi=4.0d0*atan(1.0)
         ks_loc=0.0d0
         if(t.ge.t0.and.t.le.t1)ks_loc=0.5d0*s*(dsin((n_stim*pi*x)/l)**power_stim)
         if(t.gt.t1.and.t.le.t2)ks_loc=0.25d0*s*(dcos(pi*((t-t1)/(t2-t1))))*(dsin((n_stim*pi*x)/l)**power_stim)
         return
         end        
