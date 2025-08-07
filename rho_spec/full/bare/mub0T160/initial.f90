subroutine initial(Nflow,yflow,kappa)
!make the initialization

  implicit none
  integer Nflow
  real(8) yflow(Nflow) !sigma is the fixed expansion point at UV
  integer N_str(5) !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck,Ng
  real(8) lam0,lam1,lam2,lam3,lam4,lam5,lam6,lam7
  real(8) h
  real(8) Zphi,Zpsi,ZA
  real(8) c,kappa,kappa0
  real(8) alphaS,g,b0,MU,alphaSMU,alphaS1L
  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) lambda,nu
  real(8) k_UV,k_IR,t_UV,t_IR
  real(8) Nc,Nf
  parameter(Nc=3.,Nf=2.)
  real(8) gamma_c,gammarho_Re,gammarho_Im,omega,p,epsi


  common /strucFun/ N_str
  common /kRange/k_UV,k_IR,t_UV,t_IR
  common /gamma_c_com/gamma_c
  common /omegaps/ omega,p,epsi


  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)
  Ng=N_str(5)

  b0=11.-2./3.*Nf
  MU=91.                       !GeV
  MU=MU*1.d3/hc                !fm**(0=-1)
  alphaSMU = 0.118
  alphaS1L=alphaSMU/(1 + (b0*alphaSMU*Log(k_UV/MU))/(2.*Pi))
  alphaS=alphaS1L
  alphaS=0.235
  g=sqrt(4.*pi*alphaS)
  lambda=0.
  nu=1.d4*k_UV**2
  h=1.!expansion coefficients of Yukawa coupling
  c=3.6*(1.d3/hc)**3  !explicit chiral symmetry breaking term, in unit of fm**(-3)
  kappa0=(c/nu)**2/2.!minimal point
  lam0=(kappa0**2*lambda)/2.+nu*kappa0
  lam1=kappa0*lambda+nu
  lam2=lambda
  lam3=0.
  lam4=0.
  lam5=0.
  lam6=0.
  lam7=1.!expansion coefficients of effective potential V
  Zphi=1. !meson wave function renormalization
  Zpsi=1. !quark wave function renormalization
  ZA=1.   !gluon wave function renormalization
  gammarho_Re=-omega**2+p**2!-0.27553292193124E+12
  gammarho_Im=0.

  yflow(1)=lam1
  yflow(2)=lam2
  yflow(3)=lam3
  yflow(4)=lam4
  yflow(5)=lam5
  yflow(6)=1.2724323954768*10**9!0.12724323954768E+10
  yflow(7)=lam7
  yflow(8)=gammarho_Re
  yflow(9)=gammarho_Im
  yflow(Nv+1)=lam0
  yflow((Nv+1)+1)=h
  yflow((Nv+1)+(Nh+1)+1)=Zphi
  yflow((Nv+1)+(Nh+1)+2)=Zpsi
  yflow((Nv+1)+(Nh+1)+3)=ZA
  yflow((Nv+1)+(Nh+1)+Nz+1)=c
  yflow((Nv+1)+(Nh+1)+Nz+2)=kappa0
  yflow((Nv+1)+(Nh+1)+Nz+Nck+1)=g
  yflow((Nv+1)+(Nh+1)+Nz+Nck+2)=g   !Three gluon interaction
  yflow((Nv+1)+(Nh+1)+Nz+Nck+3)=g   !quark-gluon interaction for the strange quark

end





