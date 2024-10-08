module ioFileMod

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  !	Input/Output file manipulations. Mind file on archival system, or local
  !	disk etc.
  !
  ! Author: Mariana Vertenstein
  !
  !---------------------------------------------------------------------
 
  use shr_kind_mod  , only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use spmd_utils    , only: masterproc
  use cam_logfile   , only: iulog

  implicit none

  private

  public getfil ! Get file from archive
  public opnfil ! Open file

contains
 
  subroutine getfil(fulpath, locfn, iflag, lexist)
 
    ! --------------------------------------------------------------------
    ! Determine whether file is on local disk.
    ! . first check current working directory
    ! . next check full pathname[fulpath] on disk
    ! . by default, abort if file not found.  Setting optional iflag arg
    !   to 1 overrides this behavior, and in that case the optional lexist
    !   arg is used to return status of whether the file was found or not.
    ! --------------------------------------------------------------------
 
    character(*), intent(in)  :: fulpath     ! Full pathname on local disk
    character(*), intent(out) :: locfn       ! Local file name if found in working directory,
                                             ! Set to fulpath if not found in working dir.
    integer, optional, intent(in)  :: iflag  ! Set iflag=1 to return control to caller if
                                             ! file not found.  default is to abort.
    logical, optional, intent(out) :: lexist ! When iflag=1 then getfil will return whether the
                                             ! file is found or not.  This flag is set .true.
                                             ! if the file is found, otherwise .false.
 
    integer i, ierr
    integer klen            ! length of fulpath character string
    integer maxlen          ! length of locfn input variable
    logical lexist_in       ! true if local file exists
    logical abort_on_failure
 
    abort_on_failure = .true.
    if (present(iflag)) then
        if (iflag==1) abort_on_failure = .false.
    end if
    maxlen = len(locfn)

    ! First check if file is in current working directory.
    ! Get local file name from full name: start at end. look for first "/"
    klen = len_trim(fulpath)
    i = index(fulpath, '/', back=.true.)

    if (klen - i > maxlen) then
      if (abort_on_failure) then
        call endrun('(GETFIL): local filename variable is too short for path length')
      else
        if (masterproc) write(iulog,*) '(GETFIL): local filename variable is too short for path length', klen - i, maxlen
        if (present(lexist)) lexist = .false.
        return
      end if
    end if

    locfn = fulpath(i+1:klen)
    if (len_trim(locfn) == 0) then
      call endrun('(GETFIL): local filename has zero length')
    else if (masterproc) then
      write(iulog, *) '(GETFIL): attempting to find local file ', trim(locfn)
    end if
 
    inquire(file=locfn, exist=lexist_in)
    if (present(lexist)) lexist = lexist_in
    if (lexist_in) then
      if (masterproc) write(iulog, *) '(GETFIL): using ', trim(locfn), ' in current working directory'
      return
    end if
 
    ! Second check for full pathname on disk
    if (klen > maxlen) then
      if (abort_on_failure) then
        call endrun('(GETFIL): local filename variable is too short for path length')
      else
        if (masterproc) write(iulog, *) '(GETFIL): local filename variable is too short for path length', klen, maxlen
        if (present(lexist)) lexist = .false.
        return
      end if
    end if

    locfn = trim(fulpath)
    inquire(file=locfn, exist=lexist_in)
    if (present(lexist)) lexist = lexist_in
    if (lexist_in) then
      if (masterproc) write(iulog, *) '(GETFIL): using ', trim(fulpath)
      return
    else
      if (masterproc) write(iulog, *) '(GETFIL): all tries to get file have been unsuccessful: ', trim(fulpath)
      if (abort_on_failure) then
        call endrun ('GETFIL: FAILED to get '//trim(fulpath))
      else
        return
      end if
    end if

  end subroutine getfil
 
  subroutine opnfil(locfn, iun, form, status)
 
    !---------------------------------------------------------------------------
    ! open file locfn in unformatted or formatted form on unit iun
    !---------------------------------------------------------------------------
 
    character(*), intent(in) :: locfn
    integer, intent(in) :: iun
    character(1), intent(in) :: form   ! File format: u = unformatted. f = formatted
    character(*), optional, intent(in) :: status
 
    integer ioe
    character(11) ft, st
 
    if (len_trim(locfn) == 0) then
      call endrun('(OPNFIL): local filename has zero length')
    end if
    if (form=='u' .or. form=='U') then
      ft = 'unformatted'
    else
      ft = 'formatted'
    end if
    if (present(status)) then
      st = status
    else
      st = 'unknown'
    end if
    open (unit=iun, file=locfn, status=st, form=ft, iostat=ioe)
    if (ioe /= 0) then
      if (masterproc) write(iulog, *) '(OPNFIL): failed to open file ', trim(locfn), ' on unit ', iun, ' ierr=', ioe
      call endrun('opnfil') 
    else
      if (masterproc) write(iulog, *) '(OPNFIL): Successfully opened file ', trim(locfn), ' on unit= ', iun
    end if
 
  end subroutine opnfil
 
end module ioFileMod
