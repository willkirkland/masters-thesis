!> @file make_case_matrix_mod.f90
!! @author Will Kirkland (UTK)

module case_matrix_mod
    use kind_parameters, only: ikind, rkind, ckind, charlen
    use nestle_control_data, only: lthfdbk
    use number_constants, only: addToFahrenheitForRankine, lbmPerft3_to_gPercm3
    use rbf_interp_nd_mod
    use xml_data_mod, only: xmlfile, xml_file
    implicit none

    ! Globals
    real(kind=rkind), dimension(:,:), allocatable :: &
        weights, & ! the weight of each datapoint value calculated by radial
                   ! basis functtion fitting
        scaleFactor, & ! the scaling factor and offset to adjust each process
        scaleZero      ! variable into the desired range (normally 0 - 1)
    real(kind=rkind), target, dimension(:,:,:), allocatable :: & 
        caseMatrixRI, caseMatrixRO ! case matrices, rod in and out (see below)
    real(kind=rkind), target, dimension(:,:,:,:), allocatable :: &
        xsecRI, xsecRO ! cross sections values read from lattice physics file
    real(kind=rkind), target, dimension(:,:,:,:), allocatable :: &
        weightsRI, & ! weighting arrays for spline fits calculated from using
        weightsRO    ! rbf_interp_nd.f90
    real(kind=rkind) :: temp_factor ! converts degF or degC to absolute temp.
    real(kind=rkind) :: dens_factor ! converts density to specific gravity
    integer(kind=ikind), parameter :: ROD_OUT = 0, ROD_IN = 1
    real(kind=rkind), parameter :: r0 = 0.67 ! scaling parameter for radial
                                             ! basis function fitting

contains
  !> This subroutine creates the case matrix to be used by the cross section
  !! module get_cross_section_mod.
  subroutine make_case_matrix

  ! The case matrices created by this subroutine are arrays of rank four, as
  ! follows:
  !
  !   caseMatrixRO and caseMatrixRI(i,j,k)
  !   - Rank 1 (i) is the material number
  !   - Rank 2 (j) is the case number
  !   - Rank 3 (k) is the physical parameter number, as follows:
  !     1 - rod state (0 for rod out, 1 for rod in)
  !     2 - square root of fuel absolute temperature
  !     3 - moderator absolute temperature
  !     4 - moderator density
  !     5 - soluble poison concentration
  !
  ! The crosssection matrices xsecRO and xsecRI created by
  ! this subroutine are arrays of rank four, as follows,
  ! 
  !   xsecRO and xsecRI(i,j,k,l)
  !   - Rank 1 (i) is the material number
  !   - Rank 2 (j) is the group number
  !   - Rank 3 (k) is the case number
  !   - Rank 4 (l) is the cross section number
  !
  ! The case matrices and crosssection matrices are used by 
  ! get_cross_section_mod to determine value weighting for cross section
  ! interpolation.

    ! Locals
    integer(kind=ikind) :: num_rodstates ! number of rodstates
    integer(kind=ikind) :: rodInOut ! rod state
    integer(kind=ikind) :: i, j, i_case, i_mat, i_gh, i_gid, i_xs ! indices
    integer(kind=ikind) :: num_mat, num_tid, num_gid, num_case, &
           num_gh, num_xsec, num_xsec_likes ! loop limits
    real(kind=rkind), pointer, dimension(:,:,:) :: caseMatPtr ! Pointers
    real(kind=rkind), pointer, dimension(:,:,:,:) :: xsecPtr  !

    num_rodstates = 2 !@todo determine this somehow from the XML file
                      ! two states is typical, but not always!!!

    ! Use xml file to find temperature correction factor
    if (xmlfile%units%temperature .eq. "F") then
      temp_factor = addToFahrenheitForRankine
    elseif (xmlfile%units%temperature .eq. "C") then
      temp_factor = 273.15
    elseif (xmlfile%units%temperature .eq. "R") then
      temp_factor = 0.
    elseif (xmlfile%units%temperature .eq. "K") then
      temp_factor = 0.
    else ! Assume unspecified value is in Fahrenheit
      temp_factor = addToFahrenheitForRankine
    endif

    ! Use xml file to find density correction factor
    ! Convert anything to specific gravity; add'l 
    ! factors can be added below as needed
    if (xmlfile%units%density .eq. "lb/ft^3") then
      dens_factor = 1.
    elseif (xmlfile%units%density .eq. "kg/m^3") then
      dens_factor = 1000. * lbmPerft3_to_gPercm3
    else ! Assume unspecified value is specific gravity or g/cc
      dens_factor = lbmPerft3_to_gPercm3
    endif
    !
    ! For some reason, get_cross_section() gets passed densities in lbm/ft^3
    ! when T-H feedback is off, but g/cc when T-H Feedback is on.
    ! This will need to be removed in the future if the units inconsistency is
    ! ever resolved.
    if (lthfdbk) then
      dens_factor = dens_factor / lbmPerft3_to_gPercm3
    endif

    ! Find array sizes
    num_mat    = size(xmlfile%material)
    num_tid    = size(xmlfile%material(1)%time%timestamp)
    num_gid    = size(xmlfile%material(1)%energygroup%upbound)
    num_case   = size(xmlfile%material(1)%case)
    num_gh     = size(xmlfile%material(1)%branch(1)%grouphistory)
    num_xsec   = size(xmlfile%material(1)%branch(1)%grouphistory(1)% &
                      macrocrosssection)
    num_xsec_likes = num_xsec + 1 ! include neutrons per in addition to cross sections

    ! Allocate and initialize arrays for case, cross-section, and 
    ! weight matrices
    allocate(caseMatrixRO(num_mat,num_case/num_rodstates,5))
      caseMatrixRO = 0.0
    allocate(caseMatrixRI(num_mat,num_case/num_rodstates,5))
      caseMatrixRI = 0.0
    allocate(xsecRO(num_mat,num_gid,num_case/num_rodstates,num_xsec_likes*num_tid))
      xsecRO = 0.0
    allocate(xsecRI(num_mat,num_gid,num_case/num_rodstates,num_xsec_likes*num_tid))
      xsecRI = 0.0
    allocate(weightsRO(num_mat,num_gid,num_case/num_rodstates+5,num_xsec_likes*num_tid))
      weightsRO = 0.0
    allocate(weightsRI(num_mat,num_gid,num_case/num_rodstates+5,num_xsec_likes*num_tid))
      weightsRI = 0.0

    ! store branches in caseMatrix arrays
    do i_mat = 1,num_mat
      do i=1,num_case

        ! Associate pointers based on rod state, and set value of rod state
        ! in the applicable case matrix
        call get_case_rod_state(xmlfile, i, rodInOut, i_case)
        if (rodInOut .eq. ROD_OUT) then
          caseMatPtr => caseMatrixRO
          caseMatPtr(i_mat,i_case,1) = ROD_OUT
          xsecPtr => xsecRO
        else
          caseMatPtr => caseMatrixRI
          caseMatPtr(i_mat,i_case,1) = ROD_IN
          xsecPtr => xsecRI
        endif

        ! Set the values of parameters 2 through 5 (fuel temp, moderator
        ! temp, moderator density, and poison concentration) in the applicable
        ! case matrix
        caseMatPtr(i_mat,i_case,2) = &
          sqrt(xmlfile%material(i_mat)%case(i)%fuelcondition%history% &
          temperature%value + temp_factor)
        caseMatPtr(i_mat,i_case,3) = &
          xmlfile%material(i_mat)%case(i)%coolantcondition%history% &
          temperature%value + temp_factor
        caseMatPtr(i_mat,i_case,4) = &
          xmlfile%material(i_mat)%case(i)%coolantcondition%history% &
          density%value / dens_factor
        caseMatPtr(i_mat,i_case,5) = &
          xmlfile%material(i_mat)%case(i)%coolantcondition%history% &
          concentration%value
         
        ! store cross section values in the applicable xsec arrays for each
        ! group and burnup step
        do i_gh = 1,num_gh 
          i_gid = xmlfile%material(i_mat)%branch(i)%grouphistory(i_gh)%gid
          do i_xs=1,num_xsec
            xsecPtr(i_mat, i_gid, i_case, xmlfile%material(i_mat)%branch(i)% &
              grouphistory(i_gh)%tid + (i_xs-1)*num_tid) = &
              xmlfile%material(i_mat)%branch(i)%grouphistory(i_gh)% &
              macrocrosssection(i_xs)%value
          enddo
          xsecPtr(i_mat, i_gid, i_case, xmlfile%material(i_mat)%branch(i)% &
              grouphistory(i_gh)%tid + (i_xs-1)*num_tid) = &
              xmlfile%material(i_mat)%branch(i)%grouphistory(i_gh)% &
              neutronsper%value
        enddo
      enddo
    enddo


     ! Scale case matrices so values range between 0 and 1
    allocate(scaleFactor(num_mat,5)) 
    allocate(scaleZero(num_mat,5)) 
    scaleZero = 0.
    scaleFactor = 0.
    scaleFactor(:,1) = 1. ! No need to scale rod state

    do i_mat = 1, num_mat
      do j = 2,5
      scaleZero(i_mat,j) = min(minval(caseMatrixRO(i_mat,:,j)), &
                               minval(caseMatrixRI(i_mat,:,j)))
      scaleFactor(i_mat,j) = max(maxval(caseMatrixRO(i_mat,:,j)), &
        maxval(caseMatrixRI(i_mat,:,j))) - scaleZero(i_mat,j)
      if (scaleFactor(i_mat,j) .eq. 0.) then
        if (scaleZero(i_mat,j) .eq. 0) then
          scaleFactor(i_mat,j) = 1.
        else
          scaleFactor(i_mat,j) = scaleZero(i_mat,j)
        endif
      endif
      caseMatrixRO(i_mat,:,j) = (caseMatrixRO(i_mat,:,j) - scaleZero(i_mat,j))&
        / scaleFactor(i_mat,j)
      caseMatrixRI(i_mat,:,j) = (caseMatrixRI(i_mat,:,j) - scaleZero(i_mat,j))&
        / scaleFactor(i_mat,j)
      enddo
    enddo

    ! Perform radial basis function fitting; outputs are stored in weightsRO 
    ! and weightsRI
    do i_mat=1,num_mat
      do i_gid=1,num_gid
        call rbf_weight(4,num_case/2,num_xsec_likes*num_tid, &
          transpose(caseMatrixRO(i_mat,:,2:5)), r0,phi1, &
          xsecRO(i_mat,i_gid,:,:),weightsRO(i_mat,i_gid,:,:))
        call rbf_weight(4,num_case/2,num_xsec_likes*num_tid, &
          transpose(caseMatrixRI(i_mat,:,2:5)), r0,phi1, &
          xsecRI(i_mat,i_gid,:,:),weightsRI(i_mat,i_gid,:,:))
      enddo
    enddo
     

  end subroutine make_case_matrix

subroutine get_case_rod_state(filename, case_in, case_rod_state, case_index)

  type(xml_file), intent(in) :: filename
  integer(kind=ikind), intent(in) :: case_in
  integer(kind=ikind), intent(out) :: case_rod_state, case_index

  integer(kind=ikind) :: ii, ii_out, ii_in
  logical :: match

  ! Initialize variables
  ii = 1
  ii_out = 1
  ii_in = 1
  match = .false.

  ! Sort states as rod in or rod out
  do while (.not.match)
    if (ii .eq. case_in) then
      match = .true.
    else
      if (ii .gt. size(filename%material(1)%case)) then
        STOP "ERROR subroutine case_rod_state_sort: bad case input."
      endif
      if (filename%material(1)%case(ii)%fuelcondition%history%controlrod%state .eq. "Out" &
        .or. filename%material(1)%case(ii)%fuelcondition%history%controlrod%state .eq. "out" &
        .or. filename%material(1)%case(ii)%fuelcondition%history%controlrod%state .eq. "0" &
        ) then
        ii_out = ii_out + 1
      else
        ii_in = ii_in + 1
      endif
      ii = ii + 1
    endif
  enddo
  
  ! After match, set output variables
  if (filename%material(1)%case(ii)%fuelcondition%history%controlrod%state .eq. "Out" &
    .or. filename%material(1)%case(ii)%fuelcondition%history%controlrod%state .eq. "out" &
    .or. filename%material(1)%case(ii)%fuelcondition%history%controlrod%state .eq. "0" &
    ) then
    case_rod_state = ROD_OUT
    case_index = ii_out
  else
    case_rod_state = ROD_IN
    case_index = ii_in
  endif

  return

end subroutine get_case_rod_state  

end module case_matrix_mod
