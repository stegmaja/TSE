      module my_subs
      implicit none
      contains
        
      FUNCTION cross(a, b)
      real*8, DIMENSION(3) :: cross
      real*8, DIMENSION(3), INTENT(IN) :: a, b
      cross(1) = a(2) * b(3) - a(3) * b(2)
      cross(2) = a(3) * b(1) - a(1) * b(3)
      cross(3) = a(1) * b(2) - a(2) * b(1)
      END FUNCTION cross

      FUNCTION dotp(a, b)
      real*8 dotp
      real*8, DIMENSION(3), INTENT(IN) :: a, b
      dotp=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      END FUNCTION dotp

      end module my_subs