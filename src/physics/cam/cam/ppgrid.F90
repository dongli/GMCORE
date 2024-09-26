
module ppgrid

  !-----------------------------------------------------------------------
  !
  ! Purpose:
  ! Initialize physics grid resolution parameters
  !  for a chunked data structure
  !
  ! Author:
  !
  !-----------------------------------------------------------------------

  implicit none

  ! Grid point resolution parameters
  integer pcols      ! number of columns (max)
  integer psubcols   ! number of sub-columns (max)
  integer pver       ! number of vertical levels
  integer pverp      ! pver + 1

  ! start, end indices for chunks owned by a given MPI task
  ! (set in phys_grid_init).
  integer :: begchunk =  0
  integer :: endchunk = -1

contains

  subroutine ppgrid_init(pcols_in, pver_in)

    integer, intent(in) :: pcols_in
    integer, intent(in) :: pver_in
    
    pcols    = pcols_in
    psubcols = 1
    pver     = pver_in
    pverp    = pver + 1

  end subroutine ppgrid_init

end module ppgrid
