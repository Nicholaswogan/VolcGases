
module volc
  implicit none

  double precision, parameter :: A1 = -0.4200250000201988d0
  double precision, parameter :: A2 = -2.59560737813789d0
  double precision, parameter :: M_H2O = 18.01528d0
  double precision, parameter :: M_CO2 = 44.01d0
  double precision, parameter :: C_CO2 = 0.14d0
  double precision, parameter :: C_H2O = 0.02d0
  double precision, parameter :: x = 0.01550152865954013d0
  double precision, parameter :: a_H2O = 0.54d0
  double precision, parameter :: a_CO2 = 1.d0
  double precision, parameter :: d_H2O = 2.3d0

  double precision, parameter :: e = dexp(1.d0)

  double precision :: P
  double precision :: C1, C2, C3
  double precision F1, F2, xCO2tot, xH2Otot

  contains
    subroutine solve_gases(T,PP,f_O2,mCO2tot,mH2Otot, &
                         & P_H2O, P_H2, P_CO2, P_CO, &
                         & P_CH4, alphaG, x_CO2, x_H2O)

      implicit none
      ! Inputs
      double precision, intent(in) :: T, PP, f_O2, mCO2tot, mH2Otot

      ! Outputs
      double precision, intent(out) :: P_H2O, P_H2, P_CO2, P_CO, P_CH4
      double precision, intent(out) :: alphaG, x_CO2, x_H2O

      ! other variables
      double precision K1, K2, K3
      double precision, dimension(1000) :: find1, find2
      double precision found_high, found_low, error
      integer i,j,k

      ! for optimizations
      double precision :: out(1), fvec1(1)
      double precision, dimension(1) :: diag_s
      integer iflag,n,m,info
      integer nfev
      double precision, dimension(1,1) :: fjac
      double precision, dimension(1) :: ipvt, qtf, wa1, wa2 ,wa3, wa4, xx
      double precision, dimension(8) :: init_cond, fvec2
      double precision, dimension(8,8) :: fjac2
      double precision, dimension(8) :: diag
      double precision, dimension(8) :: ipvt2, qtf2, wa12, wa22 ,wa32, wa42
      integer njev



      P = PP

      F1 = dlog(1.d0/(M_CO2*x*1.d6))+C_CO2*P/T+A1
      F2 = dlog(1.d0/(M_H2O*x*100.d0))+C_H2O*P/T+A2

      ! calculate mol fraction of CO2 and H2O in the magma
      xCO2tot=(mCO2tot/M_CO2)/x
      xH2Otot=(mH2Otot/M_H2O)/x

      ! equilibrium constants
      ! made with Nasa thermodynamic database (Burcat database)
      K1 = dexp(-29755.11319228574d0/T+6.652127716162998d0)
      K2 = dexp(-33979.12369002451d0/T+10.418882755464773d0)
      K3 = dexp(-96444.47151911151d0/T+0.22260815074146403d0)

      ! constants
      C1 = K1/f_O2**0.5d0
      C2 = K2/f_O2**0.5d0
      C3 = K3/f_O2**2.d0

      ! now solve simple system
      call linspace(dlog10(PP), -15.d0, find1)
      find1 = 10.d0**find1
      call linspace(-15.d0, dlog10(PP), find2)
      find2 = 10.d0**find2
      do i=1,1000
        call fnc_s(1,1,(/find1(i)/),out,iflag)
        if (out(1).eq.out(1)) then
          found_high = find1(i)
          exit
        endif
      enddo

      do i=1,1000
        call fnc_s(1,1,(/find2(i)/),out,iflag)
        if (out(1).eq.out(1)) then
          found_low = find2(i)
          exit
        endif
      enddo

      ! root solves
     !  subroutine lmdif(fnc_s,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfnc_s,
     ! *                 diag,mode,factor,nprint,info,nfev,fjac,ldfjac,
     ! *                 ipvt,qtf,wa1,wa2,wa3,wa4)
      n = 1
      m = 1
      fvec1(1) = 0.d0
      xx(1) = found_high
      call lmdif(fnc_s,m,n,xx,fvec1,1.49012d-8,1.49012d-8,0.0d0,10000,0.d0, &
                 diag_s,1,100.d0,-1,info,nfev,fjac,m, &
                 ipvt,qtf,wa1,wa2,wa3,wa4)

      if ((info.ne.1) .and. (info.ne.2) .and. (info.ne.3) .and. (info.ne.4)) then
        xx(1) = found_low
        call lmdif(fnc_s,m,n,xx,fvec1,1.49012d-8,1.49012d-8,0.0d0,10000,0.d0, &
                   diag_s,1,100.d0,-1,info,nfev,fjac,m, &
                   ipvt,qtf,wa1,wa2,wa3,wa4)
        if ((info.ne.1) .and. (info.ne.2) .and. (info.ne.3) .and. (info.ne.4)) then
          print*,'Convergence issues! First optimization.'
          stop
        endif
      endif


      !
      P_CO2 = xx(1)
      P_CO = C2*P_CO2
      x_CO2 = dexp(F1+a_CO2*dlog(P_CO2))
      alphaG = P*(x_CO2-xCO2tot)/(-P_CO2+P*x_CO2)
      P_H2O = (P-P_CO2-C2*P_CO2)/(1.d0+C1+C3*P_CO2)
      P_H2 = C1*P_H2O
      P_CH4 = C3*P_CO2*P_H2O**2.d0
      x_H2O = dexp(F2+a_H2O*dlog(P_H2O))
      ! use different alphaG as inital guess
      alphaG = 0.1d0

      init_cond = (/dlog(x_H2O),dlog(x_CO2),dlog(P_H2O),dlog(P_CO2) &
                    ,alphaG,dlog(P_H2),dlog(P_CH4),dlog(P_CO)/)

      ! call fcn(8,8,init_cond,fvec2,fjac2,8,2)


      ! subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
      ! *                 maxfev,diag,mode,factor,nprint,info,nfev,njev,
      ! *                 ipvt,qtf,wa1,wa2,wa3,wa4)

      call lmder(fcn,8,8,init_cond,fvec2,fjac2,8,1.49012d-8,1.49012d-8,0.0d0, &
                 10000,diag,1,100.d0,-1,info,nfev,njev, &
                 ipvt2,qtf2,wa12,wa22,wa32,wa42)




      if ((((info.ne.1) .and. (info.ne.2) .and. (info.ne.3) .and. (info.ne.4)) .or. (init_cond(5).le.0.d0))) then
        init_cond = (/dlog(x_H2O),dlog(x_CO2),dlog(P_H2O),dlog(P_CO2) &
                      ,dlog(alphaG),dlog(P_H2),dlog(P_CH4),dlog(P_CO)/)

        call lmder1(fcn1,8,8,init_cond,fvec2,fjac2,8,1.49012d-8,1.49012d-8,0.0d0, &
                   10000,diag,1,100.d0,-1,info,nfev,njev, &
                   ipvt2,qtf2,wa12,wa22,wa32,wa42)

        call fcn1(8,8,init_cond,fvec2,fjac2,8,1)
        do i=1,8
          error = error + dabs(fvec2(i))
        enddo

        if ((((info.ne.1) .and. (info.ne.2) .and. (info.ne.3) .and. (info.ne.4)) .or.(error.ge.1d-7))) then
          x_H2O = xH2Otot
          x_CO2 = xCO2tot
          P_H2O = 0.d0
          P_CO2 = 0.d0
          alphaG = 0.d0
          P_H2 = 0.d0
          P_CH4 = 0.d0
          P_CO = 0.d0
        else
          x_H2O = dexp(init_cond(1))
          x_CO2 = dexp(init_cond(2))
          P_H2O = dexp(init_cond(3))
          P_CO2 = dexp(init_cond(4))
          alphaG = dexp(init_cond(5))
          P_H2 = dexp(init_cond(6))
          P_CH4 = dexp(init_cond(7))
          P_CO = dexp(init_cond(8))
        endif

      else
        x_H2O = dexp(init_cond(1))
        x_CO2 = dexp(init_cond(2))
        P_H2O = dexp(init_cond(3))
        P_CO2 = dexp(init_cond(4))
        alphaG = init_cond(5)
        P_H2 = dexp(init_cond(6))
        P_CH4 = dexp(init_cond(7))
        P_CO = dexp(init_cond(8))
      endif

    end subroutine

    subroutine fcn(m,n,y,out,fjac,ldfjac,iflag)
      implicit none
      integer, intent(in) :: iflag, m, n
      double precision, dimension(n), intent(in) :: y
      double precision, dimension(m), intent(out) :: out
      double precision, dimension(n,m), intent(out) :: fjac
      integer,intent(in) :: ldfjac
      double precision ln_x_H2O, ln_x_CO2, ln_H2O, ln_CO2
      double precision alphaG, ln_H2, ln_CH4, ln_CO
      double precision lnxH2O, lnxCO2, lnPH2O, lnPCO2
      double precision lnPH2, lnPCH4, lnPCO

      if (iflag.eq.1) then
        ln_x_H2O = y(1)
        ln_x_CO2 = y(2)
        ln_H2O = y(3)
        ln_CO2 = y(4)
        alphaG = y(5)
        ln_H2 = y(6)
        ln_CH4 = y(7)
        ln_CO = y(8)

        out(1) = dexp(ln_H2O)+dexp(ln_CO2)+dexp(ln_H2)+dexp(ln_CH4)+dexp(ln_CO)-P
        out(2) = -ln_x_CO2+dexp(ln_x_H2O)*d_H2O+a_CO2*ln_CO2+F1
        out(3) = -ln_x_H2O+a_H2O*ln_H2O+F2
        out(4) = -xH2Otot*P + (dexp(ln_H2O)+dexp(ln_H2) &
               +2.d0*dexp(ln_CH4))*alphaG+(1.d0-alphaG)*dexp(ln_x_H2O)*P
        out(5) = -xCO2tot*P + (dexp(ln_CO2)+dexp(ln_CO) &
               +dexp(ln_CH4))*alphaG+(1-alphaG)*dexp(ln_x_CO2)*P
        out(6) = dlog(C1)+ln_H2O-ln_H2
        out(7) = dlog(C2)+ln_CO2-ln_CO
        out(8) = dlog(C3)+ln_CO2+2.d0*ln_H2O-ln_CH4
      else if (iflag.eq.2) then
        lnxH2O = y(1)
        lnxCO2 = y(2)
        lnPH2O = y(3)
        lnPCO2 = y(4)
        alphaG = y(5)
        lnPH2 = y(6)
        lnPCH4 = y(7)
        lnPCO = y(8)


        fjac(1,:)=(/ 0.d0,0.d0,E**lnPH2O,E**lnPCO2,0.d0,E**lnPH2,E**lnPCH4,E**lnPCO /)
        fjac(2,:)=(/ d_H2O*E**lnxH2O,-1.d0,0.d0,a_CO2,0.d0,0.d0,0.d0,0.d0 /)
        fjac(3,:)=(/ -1.d0,0.d0,a_H2O,0.d0,0.d0,0.d0,0.d0,0.d0 /)

        fjac(4,:)=(/ (1.d0-alphaG)*E**lnxH2O*P,0.d0,alphaG*E**lnPH2O,0.d0,2.d0*E**lnPCH4 &
                    + E**lnPH2 + E**lnPH2O - E**lnxH2O*P, &
                    - alphaG*E**lnPH2,2.d0*alphaG*E**lnPCH4,0.d0 /)

        fjac(5,:)=(/ 0.d0,(1.d0 - alphaG)*E**lnxCO2*P,0.d0,alphaG*E**lnPCO2,E**lnPCH4 &
                    + E**lnPCO + E**lnPCO2 - E**lnxCO2*P,0.d0, &
                    - alphaG*E**lnPCH4,alphaG*E**lnPCO /)
        fjac(6,:)=(/0.d0,0.d0,1.d0,0.d0,0.d0,-1.d0,0.d0,0.d0/)
        fjac(7,:)=(/0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,-1.d0/)
        fjac(8,:)=(/0.d0,0.d0,2.d0,1.d0,0.d0,0.d0,-1.d0,0.d0/)

      else
        print*,'iflag must be 1 or 2. iflag =',iflag
        stop
      endif
    end subroutine

    subroutine fcn1(m,n,y,out,fjac,ldfjac,iflag)
      implicit none
      integer, intent(in) :: iflag, m, n
      double precision, dimension(n), intent(in) :: y
      double precision, dimension(m), intent(out) :: out
      double precision, dimension(n,m), intent(out) :: fjac
      integer,intent(in) :: ldfjac
      double precision ln_x_H2O, ln_x_CO2, ln_H2O, ln_CO2
      double precision lnalphaG, ln_H2, ln_CH4, ln_CO
      double precision lnxH2O, lnxCO2, lnPH2O, lnPCO2
      double precision lnPH2, lnPCH4, lnPCO

      if (iflag.eq.1) then
        ln_x_H2O = y(1)
        ln_x_CO2 = y(2)
        ln_H2O = y(3)
        ln_CO2 = y(4)
        lnalphaG = y(5)
        ln_H2 = y(6)
        ln_CH4 = y(7)
        ln_CO = y(8)

        out(1) = dexp(ln_H2O)+dexp(ln_CO2)+dexp(ln_H2)+dexp(ln_CH4)+dexp(ln_CO)-P
        out(2) = -ln_x_CO2+dexp(ln_x_H2O)*d_H2O+a_CO2*ln_CO2+F1
        out(1) = -ln_x_H2O+a_H2O*ln_H2O+F2
        out(4) = -xH2Otot*P + (dexp(ln_H2O)+dexp(ln_H2)+2*dexp(ln_CH4))*dexp(lnalphaG)+(1-dexp(lnalphaG))*dexp(ln_x_H2O)*P
        out(5) = -xCO2tot*P + (dexp(ln_CO2)+dexp(ln_CO)+dexp(ln_CH4))*dexp(lnalphaG)+(1-dexp(lnalphaG))*dexp(ln_x_CO2)*P
        out(6) = dlog(C1)+ln_H2O-ln_H2
        out(7) = dlog(C2)+ln_CO2-ln_CO
        out(8) = dlog(C3)+ln_CO2+2*ln_H2O-ln_CH4
      else if (iflag.eq.2) then
        lnxH2O = y(1)
        lnxCO2 = y(2)
        lnPH2O = y(3)
        lnPCO2 = y(4)
        lnalphaG = y(5)
        lnPH2 = y(6)
        lnPCH4 = y(7)
        lnPCO = y(8)


        fjac(1,:)=(/ 0.d0,0.d0,E**lnPH2O,E**lnPCO2,0.d0,E**lnPH2,E**lnPCH4,E**lnPCO /)
        fjac(2,:)=(/ d_H2O*E**lnxH2O,-1.d0,0.d0,a_CO2,0.d0,0.d0,0.d0,0.d0/)
        fjac(3,:)=(/ -1.d0,0.d0,a_H2O,0.d0,0.d0,0.d0,0.d0,0.d0 /)

        fjac(4,:)=(/ E**lnxH2O*(1.d0 - E**lnalphaG)*P,0.d0,E**(lnalphaG + lnPH2O),0.d0, &
                    - E**lnalphaG*(2.d0*E**lnPCH4 + E**lnPH2 + E**lnPH2O) - &
                    E**(lnalphaG + lnxH2O)*P,E**(lnalphaG + lnPH2), &
                    - 2.d0*E**(lnalphaG + lnPCH4),0.d0 /)

        fjac(5,:)=(/ 0.d0,E**lnxCO2*(1.d0 - E**lnalphaG)*P,0.d0,E**(lnalphaG + lnPCO2), &
                    - E**lnalphaG*(E**lnPCH4 + E**lnPCO + E**lnPCO2) &
                    - E**(lnalphaG + lnxCO2)*P,0.d0,E**(lnalphaG + lnPCH4), &
                    - E**(lnalphaG + lnPCO) /)
        fjac(6,:)=(/0.d0,0.d0,1.d0,0.d0,0.d0,-1.d0,0.d0,0.d0/)
        fjac(7,:)=(/0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,-1.d0/)
        fjac(8,:)=(/0.d0,0.d0,2.d0,1.d0,0.d0,0.d0,-1.d0,0.d0/)

      else
        print*,'iflag must be 1 or 2. iflag =',iflag
        stop
      endif
    end subroutine

    subroutine fnc_s(m,n,PCO21,out,iflag)
      implicit none
      integer, intent(in) :: m,n,iflag
      double precision, intent(in) :: PCO21(n)
      double precision, intent(out) :: out(m)
      double precision PCO2
      PCO2 = PCO21(1)
      out(1) = ( -1 * ( e )**( ( F1 + a_CO2 * dlog( PCO2 ) ) ) * ( P )**( -1 ) &
                * ( P + ( -1 * PCO2 + -1 * C2 * PCO2 ) ) * ( ( 1 + ( C1 + C3 * PCO2 ) &
                ) )**( -1 ) + ( ( P )**( -1 ) * ( P + ( -1 * PCO2 + -1 * C2 * PCO2 ) &
                ) * ( ( 1 + ( C1 + C3 * PCO2 ) ) )**( -1 ) * xCO2tot + ( -1 * ( e &
                )**( ( F2 + a_H2O * dlog( (P-PCO2-C2*PCO2)/(1+C1+C3*PCO2) ) ) ) * ( -1 * ( P )**( -1 ) * PCO2 + &
                xCO2tot ) + ( ( e )**( ( F1 + a_CO2 * dlog( PCO2 ) ) ) * xH2Otot &
                + -1 * ( P )**( -1 ) * PCO2 * xH2Otot ) ) ) )
    end subroutine


    subroutine linspace(from, to, array)
      real(8), intent(in) :: from, to
      real(8), intent(out) :: array(:)
      real(8) :: range
      integer :: n, i
      n = size(array)
      range = to - from

      if (n == 0) return

      if (n == 1) then
        array(1) = from
        return
      end if


      do i=1, n
        array(i) = from + range * (i - 1) / (n - 1)
      end do
    end subroutine

end module



program test
  use volc
  implicit none


  double precision :: T, PP, f_O2, mCO2tot, mH2Otot

  ! Outputs
  double precision :: P_H2O, P_H2, P_CO2, P_CO, P_CH4
  double precision :: alphaG, x_CO2, x_H2O
  double precision A, B, C

  double precision from, to
  double precision, dimension(10) :: array



  T = 1473.d0
  PP = 1.d0
  mCO2tot = 1000d-6
  mH2Otot = 1000d-6

  A = 25738.d0
  B = 9.d0
  C = 0.092d0

  f_O2 = 10.d0**((-A/T+B+C*(P-1.d0)/T))


  call solve_gases(T,PP,f_O2,mCO2tot,mH2Otot, &
                 & P_H2O, P_H2, P_CO2, P_CO, &
                 & P_CH4, alphaG, x_CO2, x_H2O)

  print*,P_H2O, P_H2, P_CO2, P_CO, &
  & P_CH4, alphaG, x_CO2, x_H2O
end program
