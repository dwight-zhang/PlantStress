C*******************************************************************
C  Set the seed for the random number generator 
C*******************************************************************
        subroutine set_random_seed(nrank)
        integer nrank
        integer k,i,date_time(8)
        character(len=10) big_ben(3)
        integer,allocatable :: seed(:)

        call random_seed
        call random_seed(size=k)
        allocate(seed(k))
        call date_and_time(big_ben(1),big_ben(2),big_ben(3),
     +                     date_time)
        do i = 1,k
           seed(i) = date_time(8-mod(i-1,8))+i*(nrank+1)*100
        enddo
        call random_seed(put=seed(1:k))

        end subroutine set_random_seed


C*******************************************************************
C  Generate a point on the unit sphere randomly with uniform 
C  density 
C*******************************************************************
	subroutine unif_random(pt)
        real    pt(3),crd(3),dd
        integer i

10      do i = 1,3
           call random_number(crd(i))
           crd(i) = -1.0+2.0*crd(i)
        enddo
        dd = crd(1)*crd(1)+crd(2)*crd(2)+crd(3)*crd(3)
        if (dd<=1.0) then
            dd = sqrt(dd)
            do i = 1,3
               pt(i) = crd(i)/dd;
            enddo
        else
            goto 10
        endif

	end subroutine unif_random
        

C*******************************************************************
C  Sample a point on the unit sphere randomly according to the
C  density function dens_f using the Selection Method
C*******************************************************************
        subroutine random_generator(pt,dens_max)
        real  pt(3)
        real  u,dens_max
  
20      call unif_random(pt)
        call random_number(u)
        if (u>(dens_f(pt)/dens_max)) then
           goto 20
        endif

        end subroutine random_generator
