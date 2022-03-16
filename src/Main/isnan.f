      logical function isnanAHr(x)
      implicit none
      real x
      real a,b
      real Pinf,MInf
      integer*8 IPInf, IMinf
      data IPInf/B'01111111100000000000000000000000'/    ! +Infinity
      data IMInf/B'11111111100000000000000000000000'/    ! -Infinity


      a = transfer(IPinf,Pinf)
      b = transfer(IMinf,Minf)
!      if ((x-x).eq.0.0d0) then
!        isnanAH = .TRUE.
!      else
!        isnanAH = .FALSE.
!      endif
      isnanAHr = .true.
      If(x .ge. 0.0) then
        isnanAHr = .false.
      elseif (x .le. 0.0) then
        isnanAHr = .false.
      endif
      if (a.eq.x) isnanAHr = .true.
      if (b.eq.x) isnanAHr = .true.
      if (x.gt.1.0d200) isnanAHr = .true.
      if (x.lt.(-1.0d200)) isnanAHr = .true.
      return
      end

      logical function isnanAHd(x)
      implicit none
      real*8 x
      real a,b
      real Pinf,MInf
      integer*8 IPInf, IMinf
      data IPInf/B'01111111100000000000000000000000'/    ! +Infinity
      data IMInf/B'11111111100000000000000000000000'/    ! -Infinity


      a = transfer(IPinf,Pinf)
      b = transfer(IMinf,Minf)
!      if ((x-x).eq.0.0d0) then
!        isnanAH = .TRUE.
!      else
!        isnanAH = .FALSE.
!      endif
      isnanAHd = .true.
      If(x .ge. 0.0) then
        isnanAHd = .false.
      elseif (x .le. 0.0) then
        isnanAHd = .false.
      endif
      if (a.eq.x) isnanAHd = .true.
      if (b.eq.x) isnanAHd = .true.
      if (x.gt.1.0d200) isnanAHd = .true.
      if (x.lt.(-1.0d200)) isnanAHd = .true.
      return
      end

      ! TODO OPTY
      logical function isnan(x)
      implicit none
      real x
      real a,b
      real Pinf,MInf
      integer*8 IPInf, IMinf
      data IPInf/B'01111111100000000000000000000000'/    ! +Infinity
      data IMInf/B'11111111100000000000000000000000'/    ! -Infinity


      a = transfer(IPinf,Pinf)
      b = transfer(IMinf,Minf)
!      if ((x-x).eq.0.0d0) then
!        isnanAH = .TRUE.
!      else
!        isnanAH = .FALSE.
!      endif
      isnan = .true.
      If(x .ge. 0.0) then
        isnan = .false.
      elseif (x .le. 0.0) then
        isnan = .false.
      endif
      if (a.eq.x) isnan = .true.
      if (b.eq.x) isnan = .true.
      if (x.gt.1.0d200) isnan = .true.
      if (x.lt.(-1.0d200)) isnan = .true.
      return
      end



