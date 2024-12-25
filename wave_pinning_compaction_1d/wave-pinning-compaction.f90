!------------------------------------------------
!       Mori et al, BioPhys 94 3684 2008
!------------------------------------------------
implicit none
double precision,dimension(:),allocatable:: u1,u2,u3,u1m,u2m,u1p,u2p
double precision,dimension(:),allocatable:: p
double precision::dt,t,dx,x,t_wri,u1_0,u2_0
double precision::pi,D_1,D_2
double precision::delta,gama,k_0,k,f,n_hill
double precision::ks,ks_sin,ks_cos,ks_grad,ks_grad_rev
double precision::lap1,lap2
double precision::gama_chro,gama_void,l
double precision::mem
double precision::alpha1,beta1,alpha2,beta2
double precision::s,t0,t1,t2
double precision::r_inac,n_inac
integer::i,j,k_time,n_time,n_x
integer::g_i,g_t,g_stop,m
integer::n_stim,map,sym,boundary_con
common l,dx,s,t0,t1,t2,n_stim,g_i,g_t,n_x,sym
!------------------------------------------------
open(1,file='input.in',status='unknown')
!------------------------------------------------
read(1,*)u1_0,u2_0                       !initial values
read(1,*)dx,n_x                          !parameters associated to space discreetization
read(1,*)dt,n_time,t_wri                 !parameters associated to time discreetization
read(1,*)D_1,D_2                         !diffusion coefficient
read(1,*)delta,k_0,k,n_hill              !System parameters
read(1,*)alpha1,beta1                    !System parameters
read(1,*)alpha2,beta2                    !System parameters
read(1,*)gama_chro,gama_void             !chromatin dependent positive feedback parameter
read(1,*)g_i,g_t,g_stop                  !geometry initiation and termination grid and (gi<gt<n/2)
read(1,*)m,n_stim                        !compaction time step
read(1,*)map                             !mapping parameter
read(1,*)sym                             !spatial symmetry of the initial stimulation 1 for symmetric
read(1,*)s,t0,t1,t2
read(1,*)boundary_con                    !0 for zeroflux and periodic otherwise
if(g_t.gt.n_x/2.or.g_i.gt.g_t)then
        print *,'incorrect geometry parameters!'
        goto 40
end if
pi=4.0d0*atan(1.0)
!------------------------------------------------
allocate(u1(n_x),u2(n_x),u3(n_x))
allocate(u1m(n_x),u2m(n_x))
allocate(u1p(n_x),u2p(n_x))
!------------------------------------------------
!    Assign initial values and stimulations
!------------------------------------------------
l=dx*(n_x-(2*g_i))
u1=u1_0;u2=u2_0;t=0.0d0
u3=0
u3(g_i:g_t)=1;u3(n_x-g_t:n_x-g_i)=1
!------------------------------------------------
!          Assign boundary condition
!------------------------------------------------
if(boundary_con.eq.0) then
  u1(1)=u1(2);u1(n_x)=u1(n_x-1)
  u2(1)=u2(2);u2(n_x)=u2(n_x-1)
else
  u1(1)=u1(n_x-1);u1(n_x)=u1(2)
  u2(1)=u2(n_x-1);u2(n_x)=u2(2)
end if
!------------------------------------------------
!           PDE integration starts here
!------------------------------------------------
do k_time=0,n_time                 !time loop starts here
t=k_time*dt
!------------------------------------------------
!                Compaction step
!------------------------------------------------
u1m=u1;u2m=u2
if(mod(k_time,m).eq.0.and.g_t.lt.g_stop.and.t.ne.0)then

        if(map.eq.2)then
          mem=u1m(g_t+1);u1m(g_i+1:g_t+1)=u1m(g_i:g_t);u1m(g_i)=mem
          mem=u1m(n_x-g_t-1);u1m(n_x-g_t-1:n_x-g_i-1)=u1m(n_x-g_t:n_x-g_i);u1m(n_x-g_i)=mem
          mem=u2m(g_t+1);u2m(g_i+1:g_t+1)=u2m(g_i:g_t);u2m(g_i)=mem
          mem=u2m(n_x-g_t-1);u2m(n_x-g_t-1:n_x-g_i-1)=u2m(n_x-g_t:n_x-g_i);u2m(n_x-g_i)=mem
        end if
       
        if(map.eq.1)then
          u1m(g_i+1:g_t+1)=u1m(g_i:g_t);u1m(g_i)=u1_0
          u1m(n_x-g_t-1:n_x-g_i-1)=u1m(n_x-g_t:n_x-g_i);u1m(n_x-g_i)=u1_0
          u2m(g_i+1:g_t+1)=u2m(g_i:g_t);u2m(g_i)=u2_0
          u2m(n_x-g_t-1:n_x-g_i-1)=u2m(n_x-g_t:n_x-g_i);u2m(n_x-g_i)=u2_0
        end if

        g_t=g_t+1;g_i=g_i+1
        u3=0
        u3(g_i:g_t)=1;u3(n_x-g_t:n_x-g_i)=1

end if
!------------------------------------------------

do i=2,n_x-1
  x=i*dx
  
!------------------------------------------------
!                Laplacians
!------------------------------------------------

  lap1=D_1*((u1m(i+1)+u1m(i-1)-(2*u1m(i)))/(2*dx**2))
  lap2=D_2*((u2m(i+1)+u2m(i-1)-(2*u2m(i)))/(2*dx**2))

  if(u3(i).ne.0)then
    gama=gama_chro
  else
    gama=gama_void
  end if

!------------------------------------------------
!                Kinetic terms
!------------------------------------------------

  f=(u2m(i)*(k_0+((gama*u1m(i)**n_hill)/(k**n_hill+u1m(i)**n_hill))))-(delta*u1m(i))
! if(t.le.200)ks=ks_grad(t,x)
! if(t.gt.200)ks=ks_grad_rev(t,x)
  ks=ks_sin(t,x)
  u1(i)=u1m(i)+(f+beta1-(alpha1*u1m(i))+lap1+(ks*u2m(i)))*dt
  u2(i)=u2m(i)+(-f+beta2-(alpha2*u2m(i))+lap2-(ks*u2m(i)))*dt
!------------------------------------------------
end do
!------------------------------------------------
!           Assign boundary condition
!------------------------------------------------
if(boundary_con.eq.0) then
  u1(1)=u1(2);u1(n_x)=u1(n_x-1)
  u2(1)=u2(2);u2(n_x)=u2(n_x-1)
else
  u1(1)=u1(n_x-1);u1(n_x)=u1(2)
  u2(1)=u2(n_x-1);u2(n_x)=u2(2)
end if
!------------------------------------------------
!               File writing
!------------------------------------------------
do i=1,n_x
  if(u3(i).ne.0)then
          u1p(i)=u1(i)
          u2p(i)=u2(i)
  else
          u1p(i)=-1.0
          u2p(i)=-1.0
  end if
end do
if(mod(t,t_wri).eq.0)then
        do i=1,n_x
!        write((200+(t/t_wri)),*)(i*dx),u1(i),u2(i)
        write(600,*)t,(i*dx),u1p(i)
        write(500,*)t,(i*dx),u1(i)+u2(i)
        end do
        write(600,*)
        write(500,*)
endif
if(mod(t,t_wri).eq.0)then
        do i=1,n_x
        write(700,*)t,(i*dx),u2(i)
        end do
        write(700,*)
        n_inac=0
        do i=1,n_x
         if(u1(i).gt.1)then
            n_inac=n_inac+1
         end if
        end do
        write(901,*)t,n_inac/n_x
endif
!------------------------------------------------
!      conservation check
!------------------------------------------------
if(mod(t,t_wri).eq.0)write(800,*)t,((sum(u1)+sum(u2))/n_x)
!------------------------------------------------
end do             !time loop ends here
40 continue
end                !main program ends here
!------------------------------------------------

function ks_cos(t,x)
        double precision::t,x,pi
        double precision::s,t0,t1,t2,l,dx
        double precision::ks_cos,stim
        integer ::n_stim,g_i,g_t,n_x,sym
        common l,dx,s,t0,t1,t2,n_stim,g_i,g_t,n_x,sym
        pi=4.0d0*atan(1.0)
        ks_cos=0.0d0
        stim=0.0d0
        if(sym.eq.1)then
          if(x.ge.(g_i*dx).and.x.le.((n_x-g_i)*dx))stim=(dcos((n_stim*pi*(x-(g_i*dx)))/l))**10
        end if
        if(sym.eq.0)then
          if((x.ge.(g_i*dx)).and.(x.le.(g_t*dx)))stim=(dcos((n_stim*pi*(x-(g_i*dx)))/(3*dx*(g_t-g_i))))**10
          if((x.ge.((n_x-g_t)*dx)).and.(x.le.((n_x-g_i)*dx)))stim=(dcos((n_stim*pi*(x-((n_x-g_t)*dx)))/(3*dx*(g_t-g_i))))**10
        end if
        if(t.ge.t0.and.t.le.t1)ks_cos=0.5d0*s*stim
        if(t.gt.t1.and.t.le.t2)ks_cos=0.25d0*s*(1+dcos(pi*((t-t1)/(t2-t1))))*stim
        if(t.eq.0)write(900,*)x,stim
        return
        end        
!------------------------------------------------

function ks_sin(t,x)
        double precision::t,x,pi
        double precision::s,t0,t1,t2,l,dx
        double precision::ks_sin,stim
        integer ::n_stim,g_i,g_t,n_x,sym
        common l,dx,s,t0,t1,t2,n_stim,g_i,g_t,n_x,sym
        pi=4.0d0*atan(1.0)
        ks_sin=0.0d0
        stim=0.0d0
!------------------------
        if((x.ge.(g_i*dx)).and.(x.le.(g_t*dx)))stim=(dsin((n_stim*pi*(x-(g_i*dx)))/(dx*(g_t-g_i))))**8
        if((x.ge.((n_x-g_t)*dx)).and.(x.le.((n_x-g_i)*dx)))stim=(dsin((n_stim*pi*(x-((n_x-g_t)*dx)))/(dx*(g_t-g_i))))**8*1
!------------------------
        if(t.ge.t0.and.t.le.t1)ks_sin=0.5d0*s*stim
        if(t.gt.t1.and.t.le.t2)ks_sin=0.25d0*s*(1+dsin(pi*((t-t1)/(t2-t1))))*stim
!------------------------
        if(t.eq.0)write(900,*)x,ks_sin
        return
        end        
!------------------------------------------------
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
!------------------------------------------------
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
!------------------------------------------------
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
