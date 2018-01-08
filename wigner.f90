!FORTRAN code to generate Clebsch-Gordon or Wigner 3j symbol coefficients
  
  
    Program Clebsch_Gordon
  !
  !  Calculate the Wigner of Clebsch-Gordon 3J symbol values for all spherical harmonics up
  !  to L1=6, L2=6 in 
  !
  !   
  !
    integer :: i, L, k, j, L1, m1, L2, m2, r, m, L_upper, L_lower, phase(0:40), Ls,nwigner,&
      ii(30), kk(30), mm1(30), mm2(30), nterms, primes(77), iii, jjj, i1, power, itxt, jtxt
    double precision :: fact(0:40), norm, CG_coeff, CGCG_coeff(30), sum, sumc
    character*1 :: neg
    character :: line1*60, line2*60
    data primes /2,3,5,7,11, 13,17,19,23,29, 31,37,41,43,47, 51, 53, 59, 61, 67,71, &
    73, 79, 83, 87, 89, 97,101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, &
    157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, &
    241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, &
    347, 349, 353, 359, 367, 373, 379/
      nwigner=0
      fact(0) = 1.d0
      phase(0) = 1
      do i = 1, 40
        fact(i) = fact(i - 1)*i
        phase(i) = -phase(i - 1)
      end do
      open(unit=7, file="CG coefficients.dat")
      open(unit=13,file="particularCG.dat")
      do L1 = 0, 6
        do L2 = 0, L1          
!
!   Calculate all Clebsch-Gordon coefficients for spherical harmonics L1 and L2 up to 5 each
!
          L_upper = L1 + L2
          L_lower = Abs(L1 - L2)
          do L = L_lower, L_upper                               
            do m = -L, L
              nterms = 0                          
              do m2 = -L2, L2
                m1 = m - m2
                if (m1 > L1 .or. m1 < -L1) cycle
                norm = 0.d0
!
!  Now follows the Wigner formula for the 3J symbol
!
                do r = 0, 20
                  if (L2 + L - r - m1 < 0) cycle
                  if (L - m - r       < 0) cycle   
                  if (L1 - m1 - r     < 0) cycle  
                  if (L2 - L + m1 + r < 0) cycle
                  norm = norm + phase(L1 + r - m1)*fact(L1 + m1 + r)*fact(L2 + L - r - m1)/ &
                  (fact(r)*fact(L - m - r)*fact(L1 - m1 - r)*fact(L2 - L + m1 + r))
                end do
                CG_coeff = norm * &
                sqrt(fact(L + m)*fact(L - m)*fact(L1 - m1)*fact(L2 - m2)*fact(L1 + L2 - L)*(2*L + 1)/ &
                  (fact(L1 + m1)*fact(L2 + m2)*fact(L1 - L2 + L)*fact(L2 - L1 + L)*fact(L1 + L2 + L + 1)))  
                if (Abs(CG_coeff) > 1.d-10) then
                  norm = 1.d0/CG_coeff**2
                  do i = 1,245025
                    k = nint(norm*i)
                    if (Abs(k - norm*i) < 1.d-8) then
                      exit
                    end if
                  end do 
                else
                  i = 0
                  k = 1  
                end if           
                nterms = nterms + 1                
                mm1(nterms) = m1
                mm2(nterms) = m2 
                ii(nterms) = i
                kk(nterms) = k
                CGCG_coeff(nterms) = CG_coeff                
              end do         
!
!  At this point all nterms coefficients for the values L, L1, L2, and M are in CGC_coeff
!  the values of m1 and m2 are in mm1 and mm2, and the rational fractions are in ii and kk.
!  These rational fractions are for use in algebraic manipulation by the user only, and are
!  not essential for the Wigner coefficients.
!  
!   Now find the common denominator, and use that for all the rational fractions.
!
              j = 1
              do i = 1, nterms
                if (kk(i) > j) j = kk(i)
              end do
              do k = 1,13
                m1 = 0
                do i = 1, nterms
                  if (mod(j*k, kk(i)) /= 0) then
                    m1 = 1
                  end if
                end do
                if (m1 == 0) exit
              end do
              j = j*k
              k = 0
              sum = 0.d0
              sumc = 0.d0
              write(7,"(3(a,i2),a,i3,a))") "L1 =",L1," L2 =",L2, " L =",L, "  m =",m,"     m1  m2    Coefficient" 
              do i = 1, nterms
                CG_coeff = CGCG_coeff(i)
                m1 = mm1(i)
                m2 = mm2(i)
                r =  nint(((ii(i)*1.d0)*j)/kk(i))
                k = k + r
                sum = sum + (ii(i)*1.d0)/kk(i)
                sumc = sumc + CG_coeff**2
                iii = ii(i)
                jjj = kk(i)
                itxt = 1
                jtxt = 1
                line1 = "1"
                line2 = "1"
                if (iii > 0) then
                  do i1 = 1, 77
                    power = 0
                    do
                      if (mod(iii, primes(i1)) == 0) then
                        iii = iii/primes(i1)
                        power = power + 1
                      else
                        exit
                      end if
                    end do
                    if (power /= 0) then
                      if (power == 1) then
                        if (primes(i1) > 99) then
                          write(line1(itxt:),"('.',i3)")primes(i1)
                          itxt = itxt + 4
                        else if (primes(i1) > 9) then
                          write(line1(itxt:),"('.',i2)")primes(i1)
                          itxt = itxt + 3
                        else
                          write(line1(itxt:),"('.',i1)")primes(i1)
                           itxt = itxt + 2
                        end if
                      else
                        if (power > 9) then
                          neg = "2"
                        else
                          neg = "1"
                        end if
                        if (primes(i1) > 99) then
                          write(line1(itxt:),"('.',i3,'^',i"//neg//")")primes(i1),power
                          itxt = itxt + 6
                        else if (primes(i1) > 9) then
                          write(line1(itxt:),"('.',i2,'^',i"//neg//")")primes(i1),power
                          itxt = itxt + 5
                        else
                          write(line1(itxt:),"('.',i1,'^',i"//neg//")")primes(i1),power
                           itxt = itxt + 4
                        end if
                        if (neg == "2") itxt = itxt + 1
                      end if
                    else
                      do
                        if (mod(jjj, primes(i1)) == 0) then
                          jjj = jjj/primes(i1)
                          power = power + 1
                        else
                          exit
                        end if
                      end do
                      if (power /= 0) then
                        if (power == 1) then
                          if (primes(i1) > 99) then
                            write(line2(jtxt:),"('.',i3)")primes(i1)
                            jtxt = jtxt + 4
                          else if (primes(i1) > 9) then
                            write(line2(jtxt:),"('.',i2)")primes(i1)
                            jtxt = jtxt + 3
                          else
                            write(line2(jtxt:),"('.',i1)")primes(i1)
                            jtxt = jtxt + 2
                          end if
                        else
                          if (primes(i1) > 99) then
                            write(line2(jtxt:),"('.',i3,'^',i1)")primes(i1),power
                            jtxt = jtxt + 6
                          else if (primes(i1) > 9) then
                            write(line2(jtxt:),"('.',i2,'^',i1)")primes(i1),power
                            jtxt = jtxt + 5
                          else
                            write(line2(jtxt:),"('.',i1,'^',i1)")primes(i1),power
                            jtxt = jtxt + 4
                          end if
                        end if
                      end if
                    end if
                  end do
                else
                  line1 = ".0"
                end if
                if (CG_coeff > -1.d-12) then
                  neg = " "
                else
                  neg = "-"
                end if
                if (line1(1:1) == ".") line1 = line1(2:)
                if (line2(1:1) == ".") line2 = line2(2:)
                write(7,"(30x,2i4,f16.10,2(a,i6,a,i7,a),a)")  m1, m2, CG_coeff, &
                      " = "//neg//"sqrt(",ii(i),"/",kk(i),")"," =  "//neg//"sqrt(",r,"/",j,") = (", &
                      trim(line1)//")/("//trim(line2)//")"
               if(m1>=0 .and. m1+m2==0)then
                nwigner=nwigner+1
                write(13,*) L1,L2,L,m1,(-1)**(L1-L2)*CG_coeff/sqrt(dble(2*L+1))
                endif
               end do
               if (k /= j) then
                 write(7,*)" An error has been detected in the Wigner coefficients just printed"
                 write(7,"(2(a,i10))")"sum of numerators:", k, "value of denominator:", j
                 stop
               end if
               if (abs(sum - 1.d0) > 1.d-10) then
                 write(7,*)" An error has been detected in the Wigner coefficients just printed"
                 write(7,"(a,f17.12)")" Sum of squares of coefficients from rational fractions", sum
                 stop
               end if
               if (abs(sumc - 1.d0) > 1.d-10) then
                 write(7,*)" An error has been detected in the Wigner coefficients just printed"
                 write(7,"(a,f17.12)")" Sum of squares of coefficients =", sum
                 stop
               end if
              write(7,"(a)")" "
            end do
          end do
        end do
      end do
 
 close(13)
 end program

