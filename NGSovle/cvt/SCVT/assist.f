        subroutine get_angles(v1,v2,v3,angles,qual)
        implicit real(a-h,o-z)
        real  v1(3),v2(3),v3(3),angles(3),qual

	a = sqrt((v2(1)-v3(1))**2+(v2(2)-v3(2))**2+(v2(3)-v3(3))**2)
	b = sqrt((v1(1)-v3(1))**2+(v1(2)-v3(2))**2+(v1(3)-v3(3))**2)	
        c = sqrt((v1(1)-v2(1))**2+(v1(2)-v2(2))**2+(v1(3)-v2(3))**2)
        qual = (a+b-c)*(b+c-a)*(c+a-b)/(a*b*c)
	angles(1) = acos((b*b+c*c-a*a)/(2.0*b*c))
	angles(2) = acos((a*a+c*c-b*b)/(2.0*a*c))        
	angles(3) = acos((b*b+a*a-c*c)/(2.0*b*a))

        end subroutine get_angles


