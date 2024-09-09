      	program SCVT_By_Lloyd
      	include "scvt.m"

C*******************************************************************
C  Read the data sets for initialization                            
C*******************************************************************
      	open(15,file='scvt.in',status='unknown')
        read(15,*) 
      	read(15,*) ntemp
        read(15,*) 
C  Read the maximal number of iterations for Lloyd's method
      	read(15,*) max_iter
        read(15,*)
C  Read the tolerance for Lloyd's method         
	read(15,*) eps
        read(15,*)
C  Should we place one node at the poles ?(0--no -1--south 1--North)
        read(15,*) ip
        read(15,*)
C  Read the stride
        read(15,*) istride
        close(15)

C  The file "scvt_s.dat" must exist and then read it
       	open(16,file='scvt_s.dat',status='unknown')
C  Read the number of generators
       	read(16,*) n
        if (n.ge.nmax) then
           print *,"The number of generators must be less than ",nmax
           stop
        endif
C  Read the initial coordinates of the generators
       	do node = 1,n
           read(16,*) ntemp,x(node),y(node),z(node)
      	enddo
      	close(16)

        print *,"Number of Generators = ",n
        print *,"Number of Maximum Iterations = ",max_iter
        print *,"Maximum Tolerance = ",eps
        print *,"Stride Step = ",istride

C  Fix the south/north pole if necessary
        indx = 0   
        if (ip.ne.0) then
           print *,'The North Pole is fixed!'
           distmin = 1000000.0
           do node = 1,n
              dist = abs(x(node))+abs(y(node))+abs(z(node)-ip*1.0)           
              if (dist<distmin) then
                 distmin = dist
                 indx = node
              endif
           enddo
           x(indx) = 0.0
           y(indx) = 0.0
           z(indx) = ip*1.0
        endif
      	print *,'Initialization is done!'

C*******************************************************************
C  Start the Lloyd's algorithm 
C*******************************************************************
        print *,'Start Lloyd iteration ...'
      	ic = 0
      	do iloop = 1,max_iter
           if (mod(iloop-1,istride).eq.0) then 
              print *,"A re-triangulation is done!"
              call grid(n,x,y,z,xc,yc,zc,neigh,neisz,vortx,vorsz,ic)
           endif
           dm_LI = 0.0
           dm_L2 = 0.0
           do node = 1,n
              ns = neisz(node)
              snw(1:3) = 0.0
              vsnw  = 0.0
              v1(1) = x(node)
              v1(2) = y(node)
              v1(3) = z(node)
              dfw1 = dens_f(v1) 
              do i1 = 1,ns
                 i2 = mod(i1,ns)+1
                 if (i1.gt.1) then
                    i3 = mod(i1-1,ns)
                 else
                    i3 = ns
                 endif
                 if (mod(iloop-1,istride).eq.0) then
                    v2(1) = xc(vortx(node,i1))
                    v2(2) = yc(vortx(node,i1))
                    v2(3) = zc(vortx(node,i1))
                    v3(1) = xc(vortx(node,i2))
                    v3(2) = yc(vortx(node,i2))
                    v3(3) = zc(vortx(node,i2))
                 else
                    v4(1) = x(neigh(node,i1))-v1(1)
                    v4(2) = y(neigh(node,i1))-v1(2)
                    v4(3) = z(neigh(node,i1))-v1(3)
                    v5(1) = x(neigh(node,i2))-v1(1)
                    v5(2) = y(neigh(node,i2))-v1(2)
                    v5(3) = z(neigh(node,i2))-v1(3)
                    v6(1) = x(neigh(node,i3))-v1(1)
                    v6(2) = y(neigh(node,i3))-v1(2)
                    v6(3) = z(neigh(node,i3))-v1(3)

                    v2(1) = v6(2)*v4(3)-v6(3)*v4(2)
                    v2(2) = v6(3)*v4(1)-v6(1)*v4(3)
                    v2(3) = v6(1)*v4(2)-v6(2)*v4(1)
                    cnorm = sqrt(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))
                    v2(1:3) = v2(1:3)/cnorm

                    v3(1) = v4(2)*v5(3)-v4(3)*v5(2)
                    v3(2) = v4(3)*v5(1)-v4(1)*v5(3)
                    v3(3) = v4(1)*v5(2)-v4(2)*v5(1)
                    cnorm = sqrt(v3(1)*v3(1)+v3(2)*v3(2)+v3(3)*v3(3))
                    v3(1:3) = v3(1:3)/cnorm
                 endif
                 dfw2 = dens_f(v2)
                 dfw3 = dens_f(v3)
                 T_area = areas(v1,v2,v3)            
                 do j = 1,3 
                    crdwei = v1(j)*dfw1+v2(j)*dfw2+v3(j)*dfw3
                    snw(j) = snw(j)+T_area*crdwei
                 enddo
                 vsnw = vsnw+T_area*(dfw1+dfw2+dfw3)
              enddo  
C  Compute the centroid according to the density function
              snw(1:3) = snw(1:3)/vsnw
              st = sqrt((snw(1))**2+(snw(2))**2+(snw(3))**2)
              snw(1:3) = snw(1:3)/st
              dm = sqrt((x(node)-snw(1))**2+(y(node)-snw(2))**2
     .                  +(z(node)-snw(3))**2)
              dm_L2 = dm_L2+dm*dm  
              if (dm>dm_LI) dm_LI = dm
              x1(node) = snw(1)
              y1(node) = snw(2)
              z1(node) = snw(3)
           enddo	   
           do node= 1,n 
              if (node.ne.indx) then
                 x(node) = x1(node)
                 y(node) = y1(node)
                 z(node) = z1(node)
              endif
           enddo        
           dm_L2 = sqrt(dm_L2/(1.0*n))
           if (mod(iloop-1,istride).eq.0) then
           print *,'iloop = ',iloop, 'dm_LI = ',dm_LI, 'dm_L2 = ',dm_L2
           endif
           if (dm_LI.lt.eps) goto 10  
       	enddo


C  Write the final coornidate of nodes into file "scvt_lloyd.dat"
10      print *,"Iter = ",iloop-1," Maximal movement = ",dm_LI
        if (dm_LI.gt.eps) print *,"The maixmum iterations are reached!"
	
        open (17,file='scvt_lloyd.dat',status='unknown')
      	write(17,100) n
      	do node = 1,n             
           write(17,200) node,x(node),y(node),z(node)
      	enddo
      	close(17) 
100     format(I10)
200     format(I10,3X,F16.10,3X,F16.10,3X,F16.10)


C  Write the Delaunay triangles into the file "deltri.dat"
C  Write the Voronoi regions into the file "voronoi.dat"
      	ic = 1
      	call grid(n,x,y,z,xc,yc,zc,neigh,neisz,vortx,vorsz,ic)

      	end program SCVT_By_Lloyd
