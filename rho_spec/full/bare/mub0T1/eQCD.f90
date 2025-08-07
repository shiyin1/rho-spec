program QM

  implicit none

  real(8) pi,hc
  parameter(pi=3.141592653589793d0)
  parameter(hc=197.33)
  real(8) T,mu !temperature and chemical potential
  real(8) l,lb !polyakov loop
  real(8) rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,ZA,c,kappa,g,Zphi0,fpi0
  real(8) sigma_UV,kappa_UV_i,kappa_UV_i_mu,kappa_UV
  real(8) Vtotal0
  real(8) Ti,dT,mu_down,mu_up
  integer i,iTmax,j,jmumax,m,mm,i_mm,j1
  parameter(iTmax=249,jmumax=51)
  real(8) pre_res(0:jmumax,0:iTmax),T_res(0:iTmax),mu_res(0:jmumax,0:iTmax),pre_com(jmumax)
  real(8) mu_bound(0:iTmax),mu_bound_low,mu_bound_high
  real(8) T_MeV,muB_MeV
  real(8) mui,muBi,muB
  integer mmax
  parameter(mmax=10)
!order of chebyshev polynomial
  real(8) dcd0mu(mmax),dcd1mu(mmax),dcd2mu(mmax),dcd3mu(mmax),dcd4mu(mmax),chi(mmax,0:iTmax)
  real(8) factorial
  real(8) fpi_res(0:jmumax,0:iTmax),mPion_res(0:jmumax,0:iTmax),mSigma_res(0:jmumax,0:iTmax),mf_res(0:jmumax,0:iTmax)
  integer iT,iv,iTin
  real(8) l_i,lb_i,l_i_mu,lb_i_mu
  real(8) Fnf0,Fnf1,Fnf2
  external Fnf0,Fnf1,Fnf2
  real(8) nfl,nfl0,nfl1,lset
  real(8) chebev
  external chebev
  real(8) g_max,alphas_max,k_g_max,ZAm1_max,k_ZAm1_max,g3A_max,alphas3A_max,k_g3A_max
  real(8) Zphi_p0
  real(8) drho1V,drho1V_s,mass_s
  real(8) lam1,lam2,lam3,lam4,lam5,lam6,lam7
  real(8) rho_s_i,omega,epsi,p,gammarho_Re,gammarho_Im,mrho2
  real(8) hindata(250),Zindata(250),mfindata(250),lindata(250),lbindata(250),hin,mfin

  common /Tmu/ T,mu
  common /prefit/ pre_com
  common /iTiv/ iT,iv
  common /alphas_max_com/g_max,alphas_max,k_g_max,ZAm1_max,k_ZAm1_max,g3A_max,alphas3A_max,k_g3A_max
  common /Zphi_p0_com/ Zphi_p0
  common /drho1V_com/ drho1V,drho1V_s,mass_s
  common /gapPara/ lam1,lam2,lam3,lam4,lam5,lam6,lam7
  common /rho_s_i_com/ rho_s_i
  common /omegaps/ omega,p,epsi
  common /gammarho_com/ gammarho_Re,gammarho_Im,mrho2

  Ti=1./hc                 !initial temperature
  dT=1./hc                 !stepsize of temperature
  muBi=0./hc
  sigma_UV=2.8d-2/hc
  kappa_UV_i=sigma_UV**2/2.
  rho_s_i=0.203791530E+00

  !open(unit=1,file='../data/mub0/buffer/l.dat')
  !read(1,*) lindata
  !close (1)

  !open(unit=1,file='../data/mub0/buffer/lb.dat')
  !read(1,*) lbindata
  !close (1)

  !iTin=1
  l_i=1.d-10!lindata(iTin)
  lb_i=1.d-10!lbindata(iTin)

    iT=iTmax
    do i=0, 200!iT
      epsi=0.5/hc
      p=1./hc
      omega=1./hc+10.*i/hc
      T=Ti!+dT*i*2
      T_MeV=T*hc
      muB=muBi
      mu=1./3.*muB
      muB_MeV=muB*hc

      call selfEQ(kappa_UV_i,l_i,lb_i,kappa_UV,l,lb,rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,ZA,c,kappa,g)
      kappa_UV_i=kappa_UV
      l_i=l
      lb_i=lb

!goto 210
      open(unit=51,file='./buffer/Re.dat',position='append')
      write(51, "(e21.14)")(gammarho_Re)*hc**2
      close(51)

      open(unit=51,file='./buffer/Im.dat',position='append')
      write(51, "(e21.14)")(gammarho_Im)*hc**2
      close(51)

      open(unit=51,file='./buffer/mrho2bare.dat',position='append')
      write(51, "(e21.14)")(mrho2)*hc**2
      close(51)

      open(unit=51,file='./buffer/omega.dat',position='append')
      write(51, "(e21.14)")(omega)*hc
      close(51)

        open(unit=51,file='./buffer/m_rho.dat',position='append')
        write(51, "(e21.14)")sqrt(mrho2*hc**2/Zphi)
        close(51)

        open(unit=51,file='./buffer/m_rho_phy.dat',position='append')
        write(51, "(e21.14)")sqrt(mrho2*hc**2/Zphi_p0)
        close(51)

        open(unit=51,file='./buffer/Zphi_p0.dat',position='append')
        write(51, "(e21.14)")Zphi_p0
        close(51)

        open(unit=51,file='./buffer/Zphi.dat',position='append')
        write(51, "(e21.14)")Zphi
        close(51)

        open(unit=51,file='./buffer/mpion.dat',position='append')
        write(51, "(e21.14)")mpion*hc
        close(51)

        open(unit=51,file='./buffer/msigma.dat',position='append')
        write(51, "(e21.14)")msigma*hc
        close(51)

        open(unit=51,file='./buffer/mpion_phy.dat',position='append')
        write(51, "(e21.14)")mpion*hc*sqrt(Zphi/Zphi_p0)
        close(51)

        open(unit=51,file='./buffer/msigma_phy.dat',position='append')
        write(51, "(e21.14)")msigma*hc*sqrt(Zphi/Zphi_p0)
        close(51)


210 continue

      !write(*,"('j=', I4,  t25, 'i=', I4)")j, i
      write(*,"('i=', I4)")i
      write(*,"('muB_MeV=', f15.7, t25, 'T_MeV=', f15.7)")muB_MeV,T_MeV
      write(*,"('fpi=', f15.7, t25, 'mPion=', f15.7)")fpi*hc,mPion*hc*sqrt(Zphi/Zphi_p0)
      write(*,"('mf=', f15.7, t25, 'mSigma=', f15.7)")mf*hc,mSigma*hc*sqrt(Zphi/Zphi_p0)
      write(*,"('kappa_UV=', e21.14)")kappa_UV
      write(*,"('Re_rho=', e21.14)")gammarho_Re+mrho2
      write(*,"('m_rho=', f15.7)")sqrt(mrho2*hc**2/Zphi_p0)

    end do
end






