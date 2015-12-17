C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C Derivatives of R(i,j) with respect to rigid body coordinates of two
C molecules x1,y1,z1,n1,l1,m1 and x2,y2,z2,l2,m2,n2.
C P(i) = (p1x,p1y,p1z), P(j)=(p2x,p2y,p2z) are the site cordinates in
C the reference geometry for the molecule at the origin and rdist is the
C actual distance between the two sites in question.
C
C This should enable rigid body systems to be coded in a
C straightforward, systematic way.
C
      DOUBLE PRECISION FUNCTION D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,
     1     l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,
     1     M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,
     1     L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,
     1     z2,l2,m2,n2,rdist

      D1=
     @    rdist*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z)
     @     - c2a2*p2x + 
     @     c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1
     @     - n2*p2y*s2 + 
     @     m2*p2z*s2 + x1 - x2)
      
      RETURN
      END

      DOUBLE PRECISION FUNCTION D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,
     1     l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,
     1     M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,
     1     L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,
     1     z2,l2,m2,n2,rdist

      D2=
     @    rdist*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z)
     @     - c2a2*p2y + 
     @     c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1
     @     + n2*p2x*s2 - 
     @     l2*p2z*s2 + y1 - y2)

      RETURN
      END

      DOUBLE PRECISION FUNCTION D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,
     1     l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,
     1     M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,
     1     L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,
     1     z2,l2,m2,n2,rdist

      D3=
     @    rdist*(c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z)
     @     - c2a2*p2z + 
     @     c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1
     @     - m2*p2x*s2 + 
     @     l2*p2y*s2 + z1 - z2)

      RETURN
      END


      DOUBLE PRECISION FUNCTION D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,
     1     l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,
     1     M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,
     1     L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,
     1     z2,l2,m2,n2,rdist

      double precision tmp

      tmp = (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @     c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1
     @     - l1*p1y*s1 - m2*p2x*s2 + l2*p2y*s2 + z1 - z2)

      D4= rdist*((c3a1*(2*l1*p1x*(-1 + l12*ralpha12) + 
     @     (m1*p1y + n1*p1z)*(-1 + 2*l12*ralpha12)) + 
     @ (l1*l12*p1x + l12*m1*p1y - l1*n1*p1y + l1*m1*p1z + l12*n1*p1z)*
     @     ralpha12*s1 + l1*(c2a1*(n1*p1y - m1*p1z)*ralpha12- p1x*s1))*
     @     (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @     c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + 
     @     (c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     @     c3a1*m1*(-p1x + 2*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @     (-(l1*p1y) + p1z + l1* (l1*m1*p1x + n1*p1x + m1**2*p1y
     @     - l1*p1z + m1*n1*p1z)*ralpha12)*
     @     s1)*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z)
     @     - c2a2*p2y + c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1
     @     + l1*p1z*s1 + n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + 
     @     (-(c2a1*l1*(-(m1*p1x) + l1*p1y)*ralpha12) + 
     @     c3a1*n1*(-p1x + 2*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @     (-p1y - l1*p1z + l1*(-(m1*p1x) + l1*n1*p1x + l1*p1y
     @     + m1*n1*p1y + n1**2*p1z)*ralpha12)*s1)*tmp)

      RETURN
      END

      DOUBLE PRECISION FUNCTION D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,
     1     l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,
     1     M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,
     1     L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,
     1     z2,l2,m2,n2,rdist

      double precision tmp

      tmp = (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @     c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1
     @     - l1*p1y*s1 - m2*p2x*s2 + l2*p2y*s2 + z1 - z2)

      D5= rdist*((-(c2a1*m1*(-(n1*p1y) + m1*p1z)*ralpha12) + 
     @     c3a1*l1*(-p1y + 2*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @     (-(m1*p1x) - p1z + m1* (l1**2*p1x + l1*m1*p1y - n1*p1y
     @     + m1*p1z + l1*n1*p1z)*ralpha12)*s1)*
     @     (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @     c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @     n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + (-(c3a1*(l1*p1x + 2*m1*p1y
     @     + n1*p1z)) + c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12
     @     + 2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     @     (-(m1*p1y) + (l1*m12*p1x + m1*n1*p1x + m1*m12*p1y
     @     - l1*m1*p1z + m12*n1*p1z)*ralpha12)*s1)*
     @     (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @     c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 
     @     + l1*p1z*s1 + n2*p2x*s2 - l2*p2z*s2 + y1 - y2)
     @     + (c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
     @     c3a1*n1*(-p1y + 2*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @     (p1x - m1*(m1 - l1*n1)*p1x*ralpha12 + 
     @     m1*((l1 + m1*n1)*p1y*ralpha12 + p1z*(-1 + n1**2*ralpha12)))
     @     *s1)*tmp)

      RETURN
      END

      DOUBLE PRECISION FUNCTION D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,
     1     l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,
     1     M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,
     1     L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,
     1     z2,l2,m2,n2,rdist

      double precision tmp

      tmp = (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @     c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1
     @     - l1*p1y*s1 - m2*p2x*s2 + l2*p2y*s2 + z1 - z2)

      D6= rdist*((c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
     @     c3a1*l1*(-p1z + 2*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12)
     @     + (-(n1*p1x) + p1y + n1* (l1**2*p1x + l1*m1*p1y - n1*p1y
     @     + m1*p1z + l1*n1*p1z)*ralpha12)*
     @     s1)*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z)
     @     - c2a2*p2x + c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)
     @     + (n1*p1y - m1*p1z)*s1 - n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + 
     @     (-(c2a1*n1*(n1*p1x - l1*p1z)*ralpha12) + 
     @     c3a1*m1*(-p1z + 2*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @     (-p1x - n1*p1y + n1*(l1*m1*p1x + n1*p1x + m1**2*p1y
     @     - l1*p1z + m1*n1*p1z)*ralpha12)*s1)*
     @     (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @     c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1
     @     + l1*p1z*s1 + n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + 
     @     (-(c3a1*(l1*p1x + m1*p1y + 2*n1*p1z)) + 
     @     c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
     @     2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     @     (-(n1*p1z) + (-(m1*n1*p1x) + l1*n12*p1x + l1*n1*p1y
     @     + m1*n12*p1y + n1*n12*p1z)*ralpha12)*s1)*tmp)

      RETURN
      END

      DOUBLE PRECISION FUNCTION D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,
     1     l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,
     1     M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,
     1     L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,
     1     z2,l2,m2,n2,rdist

      D7=
     @   -rdist*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z)
     @     - c2a2*p2x + 
     @    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1
     @     - n2*p2y*s2 + 
     @    m2*p2z*s2 + x1 - x2)

      RETURN
      END

      DOUBLE PRECISION FUNCTION D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,
     1     l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,
     1     N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,
     1     L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,
     1     z2,l2,m2,n2,rdist

      D8=
     @   -rdist*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z)
     @     - c2a2*p2y + 
     @    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1
     @     + n2*p2x*s2 - 
     @    l2*p2z*s2 + y1 - y2)

      RETURN
      END

      DOUBLE PRECISION FUNCTION D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,
     1     l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,
     1     M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,
     1     L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,
     1     z2,l2,m2,n2,rdist

      D9=
     @   -rdist*(c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z)
     @     - c2a2*p2z + 
     @    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1
     @     - m2*p2x*s2 + 
     @    l2*p2y*s2 + z1 - z2)

      RETURN
      END

      DOUBLE PRECISION FUNCTION DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,
     1     l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,
     1     M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,
     1     L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,
     1     z2,l2,m2,n2,rdist

      double precision tmp

      tmp = (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @     c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1
     @     - m2*p2x*s2 + l2*p2y*s2 + z1 - z2)


      DA= rdist*((c3a2*(2*l2*p2x + m2*p2y + n2*p2z - 
     @     2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) - 
     @     (l2*l22*p2x + l22*m2*p2y - l2*n2*p2y + l2*m2*p2z
     @     + l22*n2*p2z)*ralpha22*s2 + l2*(-(c2a2*n2*p2y*ralpha22)
     @     + c2a2*m2*p2z*ralpha22 + 
     @     p2x*s2))*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - 
     @     c2a2*p2x + c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + 
     @     (n1*p1y - m1*p1z)*s1 - n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + 
     @     (-(c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22) + 
     @     c3a2*m2*(p2x - 2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22)
     @     - (-(l2*p2y) + p2z + l2*(l2*m2*p2x + n2*p2x + m2**2*p2y
     @     - l2*p2z + m2*n2*p2z)*ralpha22)*
     @     s2)*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y
     @     + c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1
     @     + n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + 
     @     (c2a2*l2*(-(m2*p2x) + l2*p2y)*ralpha22 + 
     @     c3a2*n2*(p2x - 2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @     (p2y + l2*p2z - l2*(-(m2*p2x) + l2*n2*p2x + l2*p2y
     @     + m2*n2*p2y + n2**2*p2z)*ralpha22)*s2)*tmp)

      RETURN
      END

      DOUBLE PRECISION FUNCTION DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,
     1     l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,
     1     M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,
     1     L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,
     1     z2,l2,m2,n2,rdist

      double precision tmp

      tmp = (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @     c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1
     @     - l1*p1y*s1 - m2*p2x*s2 + l2*p2y*s2 + z1 - z2)


      DB= rdist*((c2a2*m2*(-(n2*p2y) + m2*p2z)*ralpha22 + 
     @     c3a2*l2*(p2y - 2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @     (m2*p2x + p2z - m2*(l2**2*p2x + l2*m2*p2y - n2*p2y
     @     + m2*p2z + l2*n2*p2z)*ralpha22)*s2)*
     @     (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @     c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @     n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + (c3a2*(l2*p2x
     @     + 2*m2*p2y + n2*p2z - 2*m22*(l2*p2x + m2*p2y
     @     + n2*p2z)*ralpha22) - (l2*m22*p2x + m2*n2*p2x + m2*m22*p2y
     @     - l2*m2*p2z + m22*n2*p2z)*
     @     ralpha22*s2 + m2*(c2a2*(n2*p2x - l2*p2z)*ralpha22 + p2y*s2))*
     @     (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @     c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @     n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + (-(c2a2*m2*(m2*p2x
     @     - l2*p2y)*ralpha22) + c3a2*n2*(p2y - 2*m2*(l2*p2x + m2*p2y
     @     + n2*p2z)*ralpha22) + 
     @     (-p2x + m2*p2z - m2*(-(m2*p2x) + l2*n2*p2x + l2*p2y
     @     + m2*n2*p2y + n2**2*p2z)*ralpha22)*s2)*tmp)

      RETURN
      END

      DOUBLE PRECISION FUNCTION DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,
     1     l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,
     1     M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,
     1     L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,
     1     z2,l2,m2,n2,rdist

      double precision tmp

      tmp = (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @     c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1
     @     - l1*p1y*s1 - m2*p2x*s2 + l2*p2y*s2 + z1 - z2)


      DC= rdist*((-(c2a2*n2*(n2*p2y - m2*p2z)*ralpha22) + 
     @     c3a2*l2*(p2z - 2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @     (n2*p2x - p2y - n2*(l2**2*p2x + l2*m2*p2y - n2*p2y
     @     + m2*p2z + l2*n2*p2z)*ralpha22)*s2)*
     @     (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @     c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @     n2*p2y*s2 + m2*p2z*s2 + x1 - x2)
     @     + (c2a2*n2*(n2*p2x - l2*p2z)*ralpha22 + 
     @     c3a2*m2*(p2z - 2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @     (p2x + n2*p2y - n2*(l2*m2*p2x + n2*p2x + m2**2*p2y
     @     - l2*p2z + m2*n2*p2z)*ralpha22)*s2)*
     @     (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @     c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @     n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + (c3a2*(l2*p2x + m2*p2y
     @     + 2*n2*p2z - 2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @     (n2*(m2*p2x - l2*p2y) - n22*(l2*p2x + m2*p2y
     @     + n2*p2z))*ralpha22*s2 + n2*(-(c2a2*m2*p2x*ralpha22)
     @     + c2a2*l2*p2y*ralpha22 + p2z*s2))*tmp)

       RETURN
       END
