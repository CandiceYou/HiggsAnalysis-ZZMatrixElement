MODULE modMisc
implicit none
public :: SetDecayModes, GetDecayModes
CONTAINS

subroutine SetDecayModes(first, second)
  implicit none
  integer, intent(in) :: first,second
  include "includeVars.F90"

  DecayMode1 = first
  DecayMode2 = second

  if( (DecayMode1.le.3) .or. ( (DecayMode1.ge.8) .and. (DecayMode1.le.9) ) ) then
    M_V = M_Z
    Ga_V= Ga_Z
  elseif( ( (DecayMode1.ge.4) .and. (DecayMode1.le.6) ) .or. ( (DecayMode1.ge.10) .and. (DecayMode1.le.11) ) ) then
    M_V = M_W
    Ga_V= Ga_W    
  elseif( DecayMode1.eq.7 ) then
    M_V = 0d0
    Ga_V= 0d0    
  endif

  write(*,*) "DecayMode1: ",DecayMode1," DecayMode2: ",DecayMode2
end subroutine

subroutine GetDecayModes(first, second)
  implicit none
  integer, intent(out) :: first,second
  include "includeVars.F90"

  first = DecayMode1
  second = DecayMode2
end subroutine


END MODULE
