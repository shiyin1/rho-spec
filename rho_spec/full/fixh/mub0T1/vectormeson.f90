subroutine rhomeson(mp2,mf21,Zphi,etaphi,etapsi,h1,l,lb,T,mu,k,omega,p,epsi,dgammarho_Redt,dgammarho_Imdt,dmrho2dt)
!This subroutine calculate the kurtosis analytically

  implicit none

  real(8) pi,hc
  parameter(pi=3.1415926)
  parameter(hc=197.33)
  real(8) Nc,Nf
  parameter(Nc=3.,Nf=2.)
  real(8) v3
  parameter(v3=1./(2.*pi**2))
  integer n_x,i_x
  parameter(n_x=64)
  real(8) w_x(n_x),y_x(n_x) !Guass integral
  integer n_cth,i_cth
  parameter(n_cth=64)
  real(8) w_cth(n_cth),y_cth(n_cth) !Guass integral

  real(8) mp2,mf2,h,Zphi,etaphi,etapsi,g,T,mu,Zphi_p0,hp,mf2p,mf21,h1
  real(8) k,p,omega,epsi,x_q,q,costhe,xprime,xprime2,rF_plus_1,rF_plus_12
  real(8) cthint_ML_Re,cthint_QL_Re,fML1_Re,fML2_Re,fQL1_Re,fQL2_Re
  real(8) cthint_ML_Im,cthint_QL_Im,fML1_Im,fML2_Im,fQL1_Im,fQL2_Im
  real(8) dgammarho_Redt,dgammarho_Imdt,f2a,l,lb,b3pion,f3a,b2pion,dmrho2dt
  complex(8) B2B1m,B2B1p,F1F1m,F1F1p,F2F1m,F2F1p,p0

  p0=-(0,1)*(omega+(0,1)*epsi)
  mf2=mf21
  h=h1
  call gauleg(0.d0,1.d0, y_x, w_x, n_x)
  call gauleg(-1.d0,1.d0, y_cth, w_cth, n_cth)

  fQL1_Re=0.
  fQL2_Re=0.
  fQL1_Im=0.
  fQL2_Im=0.
  do i_x=1,n_x
    x_q=y_x(i_x)
    q=k*Sqrt(x_q)

    cthint_QL_Re=0.
    cthint_QL_Im=0.
 
    do i_cth=1,n_cth
      costhe=y_cth(i_cth)
      call FF(mf2,T,mu,l,lb,k,q,costhe,omega,p,epsi,F1F1m,F2F1m,F1F1p,F2F1p,f2a,f3a)
      !call BB(mp2,T,k,q,p,costhe,omega,epsi,B2B1m,B2B1p,b2pion,b3pion)

      xprime=x_q+(p**2/k**2)-2.*Sqrt(x_q)*(p/k)*costhe
      xprime2=x_q+(p**2/k**2)+2.*Sqrt(x_q)*(p/k)*costhe
      if(xprime<1.d0)then
        rF_plus_1=1./Sqrt(xprime)
      else
        rF_plus_1=1.
      end if

      if(xprime2<1.d0)then
        rF_plus_12=1./Sqrt(xprime2)
      else
        rF_plus_12=1.
      end if

      cthint_QL_Re=cthint_QL_Re+w_cth(i_cth)*real(2.*f2a*0.                         &
                   +(1.-(sqrt(x_q)*costhe**2-(p/k)*costhe)*rF_plus_1)*F1F1m         & 
                   +(1.-(sqrt(x_q)*costhe**2+(p/k)*costhe)*rF_plus_12)*F1F1p        &
                   +(2.*(x_q*costhe**2-sqrt(x_q)*(p/k)*costhe)/sqrt(x_q)*rF_plus_1  &
                    -(1.+xprime*rF_plus_1**2+(p0/k)**2))*F2F1m                      &
                   +(2.*(x_q*costhe**2+sqrt(x_q)*(p/k)*costhe)/sqrt(x_q)*rF_plus_12 &
                   -(1.+xprime2*rF_plus_12**2-(p0/k)**2))*F2F1p) !                     &
                   !-2.*(1.-costhe**2)*f2a                                           &
                   !-2.*(2.*(costhe**2)-2.)*f3a)

      cthint_QL_Im=cthint_QL_Im+w_cth(i_cth)*aimag(2.*f2a*0.                        &
                   +(1.-(sqrt(x_q)*costhe**2-(p/k)*costhe)*rF_plus_1)*F1F1m         & 
                   +(1.-(sqrt(x_q)*costhe**2+(p/k)*costhe)*rF_plus_12)*F1F1p        &
                   +(2.*(x_q*costhe**2-sqrt(x_q)*(p/k)*costhe)/sqrt(x_q)*rF_plus_1  &
                    -(1.+xprime*rF_plus_1**2+(p0/k)**2))*F2F1m                      &
                   +(2.*(x_q*costhe**2+sqrt(x_q)*(p/k)*costhe)/sqrt(x_q)*rF_plus_12 &
                    -(1.+xprime2*rF_plus_12**2-(p0/k)**2))*F2F1p)

    end do
    fQL1_Re=fQL1_Re+w_x(i_x)*Sqrt(x_q)*cthint_QL_Re
    fQL2_Re=fQL2_Re+w_x(i_x)*x_q*cthint_QL_Re
    fQL1_Im=fQL1_Im+w_x(i_x)*Sqrt(x_q)*cthint_QL_Im
    fQL2_Im=fQL2_Im+w_x(i_x)*x_q*cthint_QL_Im
  end do
  !call FF(mf2p,T,mu,l,lb,k,q,costhe,omega,p,epsi,F1F1m,F2F1m,F1F1p,F2F1p,f2a,f3a)
  h=1.
  dgammarho_Redt=-v3*h**2*k**2*((1.-etapsi)*fQL1_Re+etapsi*fQL2_Re)*Nc!*Zphi
  dgammarho_Imdt=-v3*h**2*k**2*((1.-etapsi)*fQL1_Im+etapsi*fQL2_Im)*Nc!*Zphi
  !dmrho2dt=-32./105.*v3*h**2*Zphi*k**2*(7.-etaphi)*b3pion*0.              &
  !         +8.*Nc*v3*h**2*k**2*(4.-etapsi)*(-5./36.*f2a+f3a/9.)   &
  !         +8./15.*v3*h**2*Zphi*k**2*(5.-etaphi)*b2pion *0.
  dmrho2dt=-8.*Nc*v3*h**2*k**2*(4.-etapsi)*(-5./36.*f2a+f3a/9.) 
  
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine FF(mf2,T,mu,l,lb,k,q,costhe,omega,p,epsi,F1F1mc,F2F1mc,F1F1pc,F2F1pc,f2a,f3a)
  !Calculating the right hand side of differential equations
  
    implicit none
  
    real(8) mf2,T,mu,l,lb,k,q,costhe,omega,epsi
    real(8) nff,nfd1xf,nfd2xf,nfa,nfd1xa,nfd2xa
    real(8) Eq,Eqmp,Eqpp,f2a,f3a
    real(8) nffqmp,nfaqmp,nffqpp,nfaqpp
    complex(8) p0,F1F1mc,F2F1mc,F1F1pc,F2F1pc
    real(8) k2,q2_qmp,q2_qpp,p2,p
    real(8) Fnf0,Fnf1,Fnf2
    external Fnf0,Fnf1,Fnf2
    real(8) x,nf0,nf1,nf2
    complex(8) Fdnfm1,Fdnfm2,Fdnfm1d1,Fdnfm2d1
    complex(8) Fdnfp1,Fdnfp2,Fdnfp1d1,Fdnfp2d1
  
    p0=-(0,1)*(omega+(0,1)*epsi)
    k2=k**2
    p2=p**2
    q2_qmp=q**2+p2-2.*q*p*costhe
    q2_qpp=q**2+p2+2.*q*p*costhe
    if(q2_qmp<k2)then
      q2_qmp=k2
    end if
    if(q2_qpp<k2)then
      q2_qpp=k2
    end if
    
    Eqmp=Sqrt(q2_qmp+k2*(mf2))
    Eqpp=Sqrt(q2_qpp+k2*(mf2))
    Eq=k*Sqrt(1. + mf2)

    x=Eqmp-mu
    nf0=Fnf0(x,T,l,lb)
    nf1=Fnf1(x,T,l,lb)
    nf2=Fnf2(x,T,l,lb)
    nffqmp=nf0 + lb*nf1 + l*nf2
  
    x=Eqmp+mu
    nf0=Fnf0(x,T,lb,l)
    nf1=Fnf1(x,T,lb,l)
    nf2=Fnf2(x,T,lb,l)
    nfaqmp=nf0 + l*nf1 + lb*nf2
  
    x=Eqpp-mu
    nf0=Fnf0(x,T,l,lb)
    nf1=Fnf1(x,T,l,lb)
    nf2=Fnf2(x,T,l,lb)
    nffqpp=nf0 + lb*nf1 + l*nf2
  
    x=Eqpp+mu
    nf0=Fnf0(x,T,lb,l)
    nf1=Fnf1(x,T,lb,l)
    nf2=Fnf2(x,T,lb,l)
    nfaqpp=nf0 + l*nf1 + lb*nf2
    
    x=Eq-mu
    nf0=Fnf0(x,T,l,lb)
    nf1=Fnf1(x,T,l,lb)
    nf2=Fnf2(x,T,l,lb)
  
    nff=nf0 + lb*nf1 + l*nf2
    nfd1xf=(-2*lb*nf1 - l*nf2 + 3*(nf0**2 + (lb*nf1 + l*nf2)**2 +                  &
           nf0*(-1 + 2*lb*nf1 + 2*l*nf2)))/T
    nfd2xf=(18*nf0**3 + 27*nf0**2*(-1 + 2*lb*nf1 + 2*l*nf2) +                     &
      (-1 + 3*lb*nf1 + 3*l*nf2)*(6*lb**2*nf1**2 + 4*lb*nf1*(-1 + 3*l*nf2) +       &
         l*nf2*(-1 + 6*l*nf2)) +                                                  &
      9*nf0*(1 + 6*lb**2*nf1**2 + 2*l*nf2*(-2 + 3*l*nf2) +                        &
         lb*nf1*(-5 + 12*l*nf2)))/T**2
    x=Eq+mu
    nf0=Fnf0(x,T,lb,l)
    nf1=Fnf1(x,T,lb,l)
    nf2=Fnf2(x,T,lb,l)
  
    nfa=nf0 + l*nf1 + lb*nf2
    nfd1xa=(-2*l*nf1 - lb*nf2 + 3*(nf0**2 + (l*nf1 + lb*nf2)**2 +                  &
           nf0*(-1 + 2*l*nf1 + 2*lb*nf2)))/T
    nfd2xa=(18*nf0**3 + 27*nf0**2*(-1 + 2*l*nf1 + 2*lb*nf2) +                      &
           (-1 + 3*l*nf1 + 3*lb*nf2)*(6*l**2*nf1**2 + 4*l*nf1*(-1 + 3*lb*nf2) +    &
           lb*nf2*(-1 + 6*lb*nf2)) +                                               &
           9*nf0*(1 + 6*l**2*nf1**2 + 2*lb*nf2*(-2 + 3*lb*nf2) +                   &
           l*nf1*(-5 + 12*lb*nf2)))/T**2
  
      Fdnfm1=(-nfa + nfaqmp)/(-Eq + Eqmp + (0,1)*p0)
      Fdnfm2=(nff - nffqmp)/(Eq - Eqmp + (0,1)*p0)
      Fdnfm1d1=(k**2*(-nfa + nfaqmp))/(2.*Eq*(-Eq + Eqmp + (0,1)*p0)**2) -        &
               (k**2*nfd1xa)/(2.*Eq*(-Eq + Eqmp + (0,1)*p0))
      Fdnfm2d1=(k**2*nfd1xf)/(2.*Eq*(Eq - Eqmp + (0,1)*p0)) -                     &
               (k**2*(nff - nffqmp))/(2.*Eq*(Eq - Eqmp + (0,1)*p0)**2)
  
    F1F1mc=(k**3*((-1. + nffqmp + nfa)/(-Eq - Eqmp + (0,1)*p0) +  &
           Fdnfm1 + Fdnfm2 +              &
           (1. - nff - nfaqmp)/(Eq + Eqmp + (0,1)*p0)))/(4.*Eq*Eqmp)
  
    F2F1mc=(k**5*(Fdnfm1 + Fdnfm2 + (1. - nfaqmp - nff)/(Eq + Eqmp+ (0,1)*p0) +      &
          (-1. + nfa + nffqmp)/(-Eq - Eqmp+ (0,1)*p0)))/(8.*Eq**3*Eqmp)-             &
          (k**3*(Fdnfm1d1 + Fdnfm2d1 + (k**2*nfd1xa)/(2.*Eq*(-Eq - Eqmp+ (0,1)*p0)) -&
          (k**2*nfd1xf)/(2.*Eq*(Eq + Eqmp+ (0,1)*p0)) -                              &
          (k**2*(1. - nfaqmp - nff))/(2.*Eq*(Eq + Eqmp+ (0,1)*p0)**2) +              &
          (k**2*(-1. + nfa + nffqmp))/(2.*Eq*(-Eq - Eqmp+ (0,1)*p0)**2)))/(4.*Eq*Eqmp)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Fdnfp1=(nfa - nfaqpp)/(Eq - Eqpp + (0,1)*p0)
      Fdnfp2=(-nff + nffqpp)/(-Eq + Eqpp + (0,1)*p0)
      Fdnfp1d1=-(k**2*(nfa - nfaqpp))/(2.*Eq*(Eq - Eqpp + (0,1)*p0)**2) + &
               (k**2*nfd1xa)/(2.*Eq*(Eq - Eqpp + (0,1)*p0))
      Fdnfp2d1=(k**2*(-nff + nffqpp))/(2.*Eq*(-Eq + Eqpp + (0,1)*p0)**2) - &
               (k**2*nfd1xf)/(2.*Eq*(-Eq + Eqpp + (0,1)*p0))
  
    F1F1pc=(k**3*((-1. + nfaqpp + nff)/(-Eq - Eqpp + (0,1)*p0) +  &
        Fdnfp1 +                  &
        Fdnfp2 +                &
        (1 - nfa - nffqpp)/(Eq + Eqpp + (0,1)*p0)))/(4.*Eq*Eqpp)
  
    F2F1pc=-(k**3*((k**2*(-1. + nfaqpp + nff))/(2.*Eq*(-Eq - Eqpp + (0,1)*p0)**2) + &
           (k**2*nfd1xf)/(2.*Eq*(-Eq - Eqpp + (0,1)*p0))                            &
           + Fdnfp1d1                                                               &
           + Fdnfp2d1                                                               &
           -(k**2*(1. - nfa - nffqpp))/(2.*Eq*(Eq + Eqpp + (0,1)*p0)**2) -          &
           (k**2*nfd1xa)/(2.*Eq*(Eq + Eqpp + (0,1)*p0))))/(4.*Eq*Eqpp) +            &
          (k**5*((-1. + nfaqpp + nff)/(-Eq - Eqpp + (0,1)*p0)                       &
          + Fdnfp1                                                                  &
          + Fdnfp2                                                                  &
          + (1. - nfa - nffqpp)/(Eq + Eqpp + (0,1)*p0)))/(8.*Eq**3*Eqpp)

    f2a=(1 - nfa + k*Sqrt(1 + mf2)*nfd1xa + k*Sqrt(1 + mf2)*nfd1xf - nff)/        &
          (4.*(1 + mf2)**1.5)

    f3a=(-(k**2*(1 + mf2)*(nfd2xa + nfd2xf)) +                                    &
        3*(1 - nfa + k*Sqrt(1 + mf2)*nfd1xa + k*Sqrt(1 + mf2)*nfd1xf - nff))/     &
        (16.*(1 + mf2)**2.5)
  end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine BB(mp2,T,k,q,p,costhe,omega,epsi,B2B1mc,B2B1pc,b2pion,b3pion)
    !Calculating the right hand side of differential equations
    
      implicit none
    
      real(8) mp2,T,k,q,p,costhe,omega,epsi
      real(8) k2,p2,q2_qmp,q2_qpp
      real(8) EqPion,EqmpPion,EqppPion,b3pion,b2pion
      real(8) nbqmpPion,nbqppPion,nbPion,nbd1xPion,nbd2xPion
      complex(8) B2B1mc,B2B1pc,p0
      real(8) Fnb
      external Fnb
      real(8) x,nb
    
      p0=-(0,1)*(omega+(0,1)*epsi)
      k2=k**2
      p2=p**2
      q2_qmp=q**2+p2-2.*q*p*costhe
      q2_qpp=q**2+p2+2.*q*p*costhe
      if(q2_qmp<k2)then
        q2_qmp=k2
      end if
      if(q2_qpp<k2)then
        q2_qpp=k2
      end if

      EqPion=k*Sqrt(1.+mp2)
      EqmpPion=Sqrt(q2_qmp+k2*mp2)
      EqppPion=Sqrt(q2_qpp+k2*mp2)
    
      x=EqPion
      nb=Fnb(x,T)
      nbPion=nb
      nbd1xPion=-((nb*(1 + nb))/T)
      nbd2xPion=(nb*(1 + nb)*(1 + 2*nb))/T**2

      x=EqmpPion
      nb=Fnb(x,T)
      nbqmpPion=nb

      x=EqppPion
      nb=Fnb(x,T)
      nbqppPion=nb

      B2B1mc=(k**3*((k**2*(1 + nbPion + nbqmpPion))/                  &
      (2.*EqPion*(-EqmpPion - EqPion + (0,1)*p0)**2) +                &
      (k**2*nbd1xPion)/(2.*EqPion*(-EqmpPion - EqPion + (0,1)*p0)) +  &
      (k**2*(-nbPion + nbqmpPion))/                                   &
      (2.*EqPion*(EqmpPion - EqPion + (0,1)*p0)**2) -                 &
      (k**2*nbd1xPion)/(2.*EqPion*(EqmpPion - EqPion + (0,1)*p0)) -   &
      (k**2*(nbPion - nbqmpPion))/                                    &
      (2.*EqPion*(-EqmpPion + EqPion + (0,1)*p0)**2) +                &
      (k**2*nbd1xPion)/(2.*EqPion*(-EqmpPion + EqPion + (0,1)*p0)) -  &
      (k**2*(-1 - nbPion - nbqmpPion))/                               &
      (2.*EqPion*(EqmpPion + EqPion + (0,1)*p0)**2) -                 &
      (k**2*nbd1xPion)/(2.*EqPion*(EqmpPion + EqPion + (0,1)*p0))))/  &
      (4.*EqmpPion*EqPion) - (k**5*                                   &
      ((1 + nbPion + nbqmpPion)/(-EqmpPion - EqPion + (0,1)*p0) +     &
      (-nbPion + nbqmpPion)/(EqmpPion - EqPion + (0,1)*p0) +          &
      (nbPion - nbqmpPion)/(-EqmpPion + EqPion + (0,1)*p0) +          &
      (-1 - nbPion - nbqmpPion)/(EqmpPion + EqPion + (0,1)*p0)))/     &
      (8.*EqmpPion*EqPion**3)

      B2B1pc=(k**3*((k**2*(1 + nbPion + nbqppPion))/                  & 
      (2.*EqPion*(-EqPion - EqppPion + (0,1)*p0)**2) +                &
      (k**2*nbd1xPion)/(2.*EqPion*(-EqPion - EqppPion + (0,1)*p0)) -  &
      (k**2*(nbPion - nbqppPion))/                                    &
      (2.*EqPion*(EqPion - EqppPion + (0,1)*p0)**2) +                 &
      (k**2*nbd1xPion)/(2.*EqPion*(EqPion - EqppPion + (0,1)*p0)) +   &
      (k**2*(-nbPion + nbqppPion))/                                   &
      (2.*EqPion*(-EqPion + EqppPion + (0,1)*p0)**2) -                &
      (k**2*nbd1xPion)/(2.*EqPion*(-EqPion + EqppPion + (0,1)*p0)) -  &
      (k**2*(-1 - nbPion - nbqppPion))/                               &
      (2.*EqPion*(EqPion + EqppPion + (0,1)*p0)**2) -                 &
      (k**2*nbd1xPion)/(2.*EqPion*(EqPion + EqppPion + (0,1)*p0))))/  &
      (4.*EqPion*EqppPion) - (k**5*                                   &
      ((1 + nbPion + nbqppPion)/(-EqPion - EqppPion + (0,1)*p0) +     &
      (nbPion - nbqppPion)/(EqPion - EqppPion + (0,1)*p0) +           &
      (-nbPion + nbqppPion)/(-EqPion + EqppPion + (0,1)*p0) +         &
      (-1 - nbPion - nbqppPion)/(EqPion + EqppPion + (0,1)*p0)))/     &
      (8.*EqPion**3*EqppPion)

      b3pion=((-3*k**5*nbd1xPion)/(4.*EqPion**4) + (k**5*nbd2xPion)/(4.*EqPion**3) + &
      (3*k**5*(0.5 + nbPion))/(4.*EqPion**5))/2.
  
      b2pion=-0.5*(k**3*nbd1xPion)/EqPion**2 + (k**3*(0.5 + nbPion))/(2.*EqPion**3)

    end