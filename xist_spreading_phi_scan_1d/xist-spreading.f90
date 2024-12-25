!---------------------------------------------------
!       Mori et al, BioPhys 94 3684 2008
!      Xist Spreading (phase 1 simulation)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
double precision,dimension(:),allocatable:: u4,u4m,u3
double precision::dt,t,dx,x,t_wri,u4_0
double precision::pi,D_4_hi,D_4_lo,D
double precision::lap4
double precision::f,decay
double precision::decay1,r
double precision::decay2,n_2,phi_2
double precision::beta4,beta4_eff,alpha4,phi_1,n_1
double precision::xic_begin,xic_end,u4_thr,u4_ss,thr_coeff
integer::i,k_time,n_time,n_x,n_stim,i_mem
common n_stim,n_x,dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(1,file='input.in',status='unknown')
open(2,file='phi1.in',status='unknown')
open(3,file='phi2.in',status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
read(1,*)u4_0                        !initial values
read(1,*)dx,n_x                      !parameters associated to space discreetization
read(1,*)dt,n_time,t_wri             !parameters associated to time discreetization
read(1,*)beta4,alpha4,n_1            !system parameters
read(1,*)D_4_hi,D_4_lo               !diffusion coefficient
read(1,*)r                           !decay1 parameters
read(1,*)n_2                         !decay2 parameters
read(1,*)n_stim                      !number of stimuli
read(1,*)xic_begin,xic_end           !location of X inactivation center
read(1,*)thr_coeff
read(2,*)phi_1
read(3,*)phi_2
pi=4.0d0*atan(1.0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(xic_end.gt.(dx*n_x))then
        write(*,*)'incorrect location for transcription!!'
        goto 10
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(u4(n_x),u4m(n_x),u3(n_x))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign initial values and stimulations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!u4=u4_0;t=0.0d0
!u4_ss=(sqrt(phi_1**2+((4*phi_1*beta4)/alpha4))-phi_1)/2
!u4_thr=thr_coeff*u4_ss
!print*,u4_ss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Assign zero flux boundary condition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
u4(1)=u4(2);u4(n_x)=u4(n_x-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      PDE integration starts here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do k_time=0,n_time                 !time loop starts here
t=k_time*dt
!----------------------------------------
u4m=u4
do i=2,n_x-1
x=i*dx
!----------------------------------------
! Laplacians
!----------------------------------------
  if(mod(i,20).eq.0)i_mem=i
  if(i.ge.i_mem.and.i.le.i_mem+5)then
          D=D_4_hi
  else
          D=D_4_lo
  end if
  if(t.eq.0)write(200,*)x,D
  lap4=D*((u4m(i+1)+u4m(i-1)-(2*u4m(i)))/(2*dx**2))
!----------------------------------------
!  Kinetic terms
!----------------------------------------
  if(x.ge.xic_begin.and.x.le.xic_end)then
        beta4_eff=beta4
  else
        beta4_eff=0.0
  end if
  ! decay1=-((alpha4*u4m(i))/(1+(r*exp((r*u4m(i))/u4_ss))))
  decay2=-(alpha4*u4m(i))/(1+(u4m(i)/phi_2)**-n_2)
  decay=decay2
  f=((beta4_eff)/(1+(u4m(i)/phi_1)**n_1))+decay
  u4(i)=u4m(i)+(f+lap4)*dt
  !+(ks_loc(t,x))
!----------------------------------------
end do
!----------------------------------------
!  Boundary Condition (Zero Flux)
!----------------------------------------
u4(1)=u4(2);u4(n_x)=u4(n_x-1)
!----------------------------------------
u3=0
do i=1,n_x
        if(mod(i,20).eq.0)i_mem=i
        if(i.ge.i_mem.and.i.le.i_mem+5.and.u4(i).gt.u4_thr)u3(i)=1.0
end do
!----------------------------------------
! File writing
!----------------------------------------
if(mod(t,t_wri).eq.0)then
        do i=1,n_x
!        write((200+(t/t_wri)),*)(i*dx),u3(i),u4(i)
        write(600,*)t,(i*dx),u4(i)
        write(601,*)t,(i*dx),u3(i)
        end do
        write(600,*)
        write(601,*)
endif
!----------------------------------------
end do             !time loop ends here
10 continue
end                !main program ends here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function ks_loc(t,x)
        double precision::t,x,pi,dx
        double precision::s,t0,t1,t2
        double precision::ks_loc
        integer :: n_stim,n_x
        common n_stim,n_x,dx
        pi=4.0d0*atan(1.0)
        s=0.05d0;t0=0.0d0;t1=20.0d0;t2=25.0d0
        ks_loc=0.0d0
        n=10
        if(x.ge.0.0.and.x.le.(dx*n_x))then
                if(t.ge.t0.and.t.le.t1)ks_loc=0.5d0*s*(dsin((n_stim*pi*x)/(dx*n_x))**10)
                if(t.gt.t1.and.t.le.t2)ks_loc=0.25d0*s*(dsin((n_stim*pi*x)/(dx*n_x))**10)
        end if  
        if(t.eq.0)write(200,*)x,ks_loc      
        return
        end
