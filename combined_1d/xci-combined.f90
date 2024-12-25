!---------------------------------------------------
!       Mori et al, BioPhys 94 3684 2008
!---------------------------------------------------
implicit none

!----------------------------------------
!     wave pinning parameters
!----------------------------------------

double precision,dimension(:),allocatable:: u1,u2,u1m,u2m,t_mem,u1p
integer,dimension(:),allocatable:: u3,u3p
double precision::dt,t,dx,x,t_wri,u1_0,u2_0
double precision::pi,D_1,D_2
double precision::delta,gama,k_0,k,f,n_hill
double precision::ks,ks_loc,ks_grad,ks_grad_rev
double precision::lap1,lap2
double precision::gama_chro,gama_void,l
double precision::mem
double precision::alpha1,beta1,alpha2,beta2
double precision::ks_0,stim_time
integer::i,j,k_time,n_time,n_x
integer::g_i,g_t,g_stop,m

!----------------------------------------
!     xist spreading parameters
!----------------------------------------

double precision,dimension(:),allocatable:: u4,u4m
double precision::lap4,D_4                              
double precision::xic_begin,xic_end
double precision::r1,r2
double precision::beta4,beta4_0,beta4_eff,phi_1,tlim_beta4
double precision::alpha4,phi_2
double precision::dummy
double precision::thr_coeff,prob_fac
double precision::r,u4_thr
double precision::u4_0
double precision::decay,prod
integer::n_1,n_2

!----------------------------------------

open(1,file='pinning_input.in',status='unknown')
open(2,file='spreading_input.in',status='unknown')

!----------------------------------------

read(1,*)u1_0,u2_0                       !initial values
read(1,*)dx,n_x                          !parameters associated to space discreetization
read(1,*)dt,n_time,t_wri                 !parameters associated to time discreetization
read(1,*)D_1,D_2                         !diffusion coefficient
read(1,*)delta,k_0,k,n_hill              !System parameters
read(1,*)alpha1,beta1                    !System parameters
read(1,*)alpha2,beta2                    !System parameters
read(1,*)gama_chro,gama_void             !chromatin dependent positive feedback parameter
read(1,*)g_i,g_t,g_stop                  !geometry initiation and termination grid and (gi<gt<n/2)
read(1,*)m                               !compaction time step
read(1,*)ks_0,stim_time                  !stimulation parameters. Stimulation amplitude and time respectively.

!----------------------------------------

read(2,*)u4_0                        !initial value
read(2,*)r1,r2                       !switch for hill functions for prod and decay respectively
read(2,*)phi_1,phi_2                 !feedback strengths for prod and decay respectively
read(2,*)n_1,n_2                     !hill coefficients for prod and decay respectively
read(2,*)D_4                         !diffusion coefficient
read(2,*)xic_begin,xic_end           !location of X inactivation center
read(2,*)thr_coeff                   !threshold coefficient for xist tethering concentration
read(2,*)beta4,tlim_beta4            !maximal production rate
read(2,*)alpha4                      !maximal decay rate
read(2,*)prob_fac                    !probability of xist tethering
!----------------------------------------

call random_seed

!----------------------------------------

if(g_t.gt.n_x/2.or.g_i.gt.g_t)then
        print *,'incorrect geometry parameters g_t and/or g_i!'
        goto 40
end if
if(xic_end.gt.(dx*n_x))then
        write(*,*)'incorrect transcription site parameter xic_end!'
        goto 40
end if

!----------------------------------------

pi=4.0d0*atan(1.0)
beta4_0=beta4

!----------------------------------------

allocate(u1(n_x),u1m(n_x),u1p(n_x))
allocate(u3(n_x),u3p(n_x))
allocate(u2(n_x),u2m(n_x))
allocate(u4(n_x),u4m(n_x))
allocate(t_mem(n_x))

!----------------------------------------
! Assign initial values and stimulations
!----------------------------------------

l=dx*(n_x-(2*g_i))
u1=u1_0;u2=u2_0;t=0.0d0
u3=0
u3(g_i:g_t)=1;u3(n_x-g_t:n_x-g_i)=1

!----------------------------------------
!  Assign zero flux boundary condition
!----------------------------------------

u1(1)=u1(2);u1(n_x)=u1(n_x-1)
u2(1)=u2(2);u2(n_x)=u2(n_x-1)
u4(1)=u4(2);u4(n_x)=u4(n_x-1)

!----------------------------------------
!      PDE integration starts here
!----------------------------------------

do k_time=0,n_time+1                 !time loop starts here
t=k_time*dt

!----------------------------------------
!       Compaction step
!----------------------------------------

u1m=u1;u2m=u2;u4m=u4
if(mod(k_time,m).eq.0.and.g_t.le.g_stop.and.t.ne.0)then

        mem=u1m(g_t+1);u1m(g_i+1:g_t+1)=u1m(g_i:g_t);u1m(g_i)=mem
        mem=u1m(n_x-g_t-1);u1m(n_x-g_t-1:n_x-g_i-1)=u1m(n_x-g_t:n_x-g_i);u1m(n_x-g_i)=mem
        mem=u2m(g_t+1);u2m(g_i+1:g_t+1)=u2m(g_i:g_t);u2m(g_i)=mem
        mem=u2m(n_x-g_t-1);u2m(n_x-g_t-1:n_x-g_i-1)=u2m(n_x-g_t:n_x-g_i);u2m(n_x-g_i)=mem
        mem=u3(g_t+1);u3(g_i+1:g_t+1)=u3(g_i:g_t);u3(g_i)=mem
        mem=u3(n_x-g_t-1);u3(n_x-g_t-1:n_x-g_i-1)=u3(n_x-g_t:n_x-g_i);u3(n_x-g_i)=mem

        g_t=g_t+1;g_i=g_i+1

end if

!----------------------------------------

do i=2,n_x-1
x=i*dx

!----------------------------------------
! Laplacians
!----------------------------------------

  lap1=D_1*((u1m(i+1)+u1m(i-1)-(2*u1m(i)))/(dx**2))
  lap2=D_2*((u2m(i+1)+u2m(i-1)-(2*u2m(i)))/(dx**2))
  lap4=D_4*((u4m(i+1)+u4m(i-1)-(2*u4m(i)))/(dx**2))

  if(u3(i).ne.0)then
    gama=gama_chro
  else
    gama=gama_void
  end if

!----------------------------------------
!       Kinetic terms xist spreading
!----------------------------------------
  if(t.gt.tlim_beta4)then
          beta4=beta4_0*exp(-(t-tlim_beta4))
  end if

  if(x.ge.xic_begin.and.x.le.xic_end)then
        beta4_eff=beta4
  else
        beta4_eff=0.0
  end if
  if(r2.eq.0)then
        dummy=0
  else
        dummy=((u4m(i)/phi_2)**-n_2)
  end if
  decay=-(alpha4*u4m(i))/(1+dummy)
  prod=((beta4_eff)/(1+(r1*(u4m(i)/phi_1)**n_1)))
  f=prod+decay
  u4(i)=u4m(i)+(f+lap4)*dt

!----------------------------------------
!     activation based on xist conc
!----------------------------------------

  u4_thr=thr_coeff*u4(int(xic_begin/dx))
  if(u4(i).gt.u4_thr.and.u3(i).eq.1)then
       call random_number(harvest=r)
       if(r.ge.(1-prob_fac))then
               u3(i)=2
               t_mem(i)=dt*k_time
       else
               u3(i)=-1
       end if
  end if

!----------------------------------------
!      Kinetic terms wave pinning
!----------------------------------------

  f=(u2m(i)*(k_0+((gama*u1m(i)**n_hill)/(k**n_hill+u1m(i)**n_hill))))-(delta*u1m(i))
!----------------------------------------
! if(t.le.200)ks=ks_grad(t,x)
! if(t.gt.200)ks=ks_grad_rev(t,x)
!----------------------------------------
  if(u3(i).eq.2.and.(t-t_mem(i)).le.stim_time)then
        ks=ks_0
  else
        ks=0
  end if
  u1(i)=u1m(i)+(f+beta1-(alpha1*u1m(i))+lap1+(ks*u2m(i)))*dt
  u2(i)=u2m(i)+(-f+beta2-(alpha2*u2m(i))+lap2-(ks*u2m(i)))*dt
!----------------------------------------
!  if(u1(i)-u1_0.lt.dt)u1(i)=u1_0
!  if(u2(i)-u2_0.lt.dt)u2(i)=u2_0
!----------------------------------------
end do                !space loop ends here

!----------------------------------------
!     Boundary Condition (zero flux)
!----------------------------------------

u1(1)=u1(2);u1(n_x)=u1(n_x-1)
u2(1)=u2(2);u2(n_x)=u2(n_x-1)
u4(1)=u4(2);u4(n_x)=u4(n_x-1)

!----------------------------------------
!      define u1p and u3p (for plotiing) 
!----------------------------------------

do i=1,n_x
        if(u3(i).eq.0)u3p(i)=-1
        if(u3(i).eq.1)u3p(i)=0
        if(u3(i).eq.-1)u3p(i)=0
        if(u3(i).eq.2)u3p(i)=1
end do
do i=1,n_x
        if(u3(i).ne.0)then
                u1p(i)=u1(i)
        else
                u1p(i)=-1
        end if
end do

!----------------------------------------
! File writing
!----------------------------------------

if(mod(t,t_wri).eq.0)then
        do i=1,n_x
        write(600,*)t,(i*dx),u1p(i)
        end do
        write(600,*)
endif
if(mod(t,t_wri).eq.0)then
        do i=1,n_x
        write(700,*)t,(i*dx),u3p(i)
        end do
        write(700,*)
endif
if(mod(t,t_wri).eq.0)then
        do i=1,n_x
        write(800,*)t,(i*dx),u4(i)
        end do
        write(800,*)
endif

!----------------------------------------
!      conservation check
!----------------------------------------
!if(mod(t,t_wri).eq.0)write(800,*)t,((sum(u1)+sum(u2))/n_x)
!----------------------------------------

end do             !time loop ends here
40 continue
end                !main program ends here
