!> @file get_cross_section_mod.f90
!> @author Will Kirkland (UTK)

module get_cross_section_mod
  use kind_parameters, only: ikind, rkind, ckind, charlen
  use case_matrix_mod, only: scaleFactor, scaleZero, caseMatrixRI, &
    caseMatrixRO, weightsRI, weightsRO, temp_factor, dens_factor, &
    ROD_OUT, ROD_IN, r0
  use ISO_C_BINDING, only: C_DOUBLE, C_SIZE_T
  use rbf_interp_nd_mod, only: rbf_interp_nd, phi1
  use xml_data_mod, only: xmlfile

  implicit none

  private
  public get_cross_section, find_group_history, find_cross_section_id

contains

    !> Returns the value of a nuclear cross-section for a given material state.
    !! Uses the case matrix constructed in case_matrix_mod.f90
    !! 
    !! @param materialID     :material number in lattice physics file
    !! @param burnupVal      :material burnup value
    !! @param groupNum       :energy group for which xsec is evaluated
    !! @param rodState       :Rod status out (=0) or in (!=0)
    !! @param tFuel          :Fuel temperature
    !! @param tMod           :Moderator temperature
    !! @param densMod        :Moderator density
    !! @param solPoisonConc  :Soluble poison concentration
    !! @param crossSectionType :String containing the xml file description of the cross-section type
    !! @param toGroupNum       :Number of destination group, for scattering xsecs
    real(kind=rkind) function get_cross_section(materialID, &
       burnupVal, groupNum, rodState, tFuel, tMod, densMod, solPoisonConc, &
       crossSectionType, toGroupNum) result(crossSection)

       ! Arguments
       integer(kind=ikind), intent(in) :: &
         materialID, & ! material number in lattice physics file
         groupNum,   & ! energy group
         rodState      ! rod out = 0, rod in != 0
       integer(kind=ikind), optional, intent(in) :: &
         toGroupNum    ! energy group scattered to, for scattering xsecs
       real(kind=rkind), intent(in) :: &
         burnupVal, &  ! fuel burnup value
         tFuel,     &  ! fuel temperature
         tMod,      &  ! moderator temperature
         densMod,   &  ! moderator density
         solPoisonConc ! soluable poison concentration
       character(kind=ckind, len=charlen), intent(in) :: crossSectionType

       ! Locals
       integer(kind=ikind) :: i, & 
         tidLo,    &  !> lowest burnup step used for spline fit (=1)
         tidHi,    &  !> highest burnup step used for spline fit (=num_tid)
         groupHistoryID, & !> the group history number for a given 
                      !> energy group and burnup step
         crossSectionID, & !> the integer index of the xsec type
         num_case, &  !> number of branch cases in the case matrix
         num_tid, &   !> number of burnup steps in the lattice physics file
         num_spline_pts  !> number of burnup steps used for spline
       integer(kind=C_SIZE_T) :: n !> num_tid cast to size_t type
       real(kind=rkind) :: tid !> the applicable burnup step
       real(kind=rkind), dimension(:), allocatable :: &
         xVals, yVals !> x and y for 1D cubic spline burnup fits
       real(kind=rkind), dimension(1:5) :: &
         inputVector  !> input state vector, projected onto case matrix
       real(kind=C_DOUBLE) :: crossSection_cdbl !> xsec from C spline function
       real(kind=rkind), pointer, dimension(:,:,:) :: caseMatPtr   ! Pointers
       real(kind=rkind), pointer, dimension(:,:,:,:) :: weightsPtr !
       real(kind=rkind), dimension(:), allocatable :: tid_list

       logical, parameter :: fourPointSplineOnly = .false. ! if true, only four
               ! points in the burnup curve are used to fit the spline. This is
               ! slightly faster but may be less accurate.

       ! Calculate tid from burnup value
       num_tid = size(xmlfile%material(materialID)%time%timestamp)
       allocate(tid_list(num_tid))
       tid_list = 0.
       do i=1,num_tid
         tid_list(i) = xmlfile%material(materialID)%time%timestamp(i)%burnup
       enddo
       i = 1
       tid = -1
       do while (tid .eq. -1)
         if (burnupVal .gt.  tid_list(i)) then
           if (i .lt. num_tid) then
             i=i+1
           else
             i = num_tid
             tid = real(i)
           endif
         elseif (i .eq. 1) then
           tid = 1.
         else
           tid = burnupVal - tid_list(i-1)
           tid = tid / (tid_list(i) - tid_list(i-1))
           tid = tid + real(i-1)
         endif
       enddo

       ! Convert input values to state vector
       if (rodState .eq. ROD_OUT) then
         caseMatPtr => caseMatrixRO
         weightsPtr => weightsRO
       else 
         caseMatPtr => caseMatrixRI
         weightsPtr => weightsRI
       endif
       inputVector(1) = rodState 
       inputVector(2) = sqrt(tFuel) ! Note that xsfdbk passes tFuel in Rankine (hence no temp_factor here) while passing tMod in Fahrenheit
       inputVector(3) = tMod + temp_factor
       inputVector(4) = densMod
       inputVector(5) = solPoisonConc
       do i=1,5
         inputVector(i) = &
           (inputVector(i) - scaleZero(materialID,i)) / scaleFactor(materialID,i)
       enddo

       num_case = size(caseMatPtr(1,:,1))

       ! Return right away if the burnup is at an evaluated point in the lattice
       ! physics file
       num_spline_pts = num_tid
       do i=1,num_tid
         if (burnupVal .eq. tid_list(i)) then
           tidLo = i
           tidHi = i
           num_spline_pts = 1
         endif
       enddo

       ! If burnup is below first burnup point or above last burnup point, use
       ! a linear extrapolation
       if (num_spline_pts .ne. 1) then
         if (burnupVal .lt. tid_list(1)) then
           tidLo = 1
           tidHi = 2
           num_spline_pts = 2
         elseif (burnupVal .gt. tid_list(num_tid)) then
           tidLo = num_tid - 1
           tidHi = num_tid
           num_spline_pts = 2

         ! Otherwise, adjust for whether four or all burnup points will be used to make the spline
         elseif (fourPointSplineOnly) then
           tidLo = 1
           tidHi = num_tid
           if (tid .lt. 2) then
             tidLo = 1
             tidHi = 4
           elseif (tid .ge. (size(xmlfile%material(materialID)%time%timestamp)-1)) then
             tidLo = size(xmlfile%material(materialID)%time%timestamp) - 3
             tidHi = tidLo + 3
           else
             tidLo = floor(tid) - 1
             tidHi = floor(tid) + 2
           endif
           num_spline_pts = 4
         else
           tidLo = 1
           tidHi = num_tid
           num_spline_pts = num_tid
         endif
       endif
       
       allocate(yVals(tidLo:tidHi))
       if (num_spline_pts .ne. 1) allocate(xVals(tidLo:tidHi))

       ! Get a vector of cross section values for tid's used for spline
       ! for the given xsec, material, and its physical state.
       do i = tidLo, tidHi
         groupHistoryID = find_group_history(materialID,1,i,groupNum)
         if (crossSectionType .eq. "scatter") then
           crossSectionID = find_cross_section_id(materialID,1, &
             groupHistoryID,crossSectionType,toGroupNum)
         else
           crossSectionID = find_cross_section_id(materialID,1, &
             groupHistoryID,crossSectionType)
         endif
         call rbf_interp_nd(4,num_case,transpose(caseMatPtr( &
           materialID,:,2:5)),r0,phi1,weightsPtr(materialID,groupNum,:, i + &
           (crossSectionID-1) * num_tid),inputVector(2:5), yVals(i))

         if (num_spline_pts .ne. 1) xVals(i) = xmlfile%material(materialID)%time%timestamp(i)%burnup
       enddo

       ! Return burnup value if it matches a lattice physics burnup point.
       if (num_spline_pts .eq. 1) then
         crossSection = yVals(1)

       ! Or use linear extrapolation if burnup point is outside of range
       elseif (num_spline_pts .eq. 2) then
         crossSection = (burnupVal - xVals(tidLo)) / (xVals(tidHi) - xVals(tidLo)) * (yVals(tidHi) - yVals(tidLo)) + yVals(tidLo)

       ! Else, fit a cubic spline to the cross-sections as a function of burnup,
       ! and read from that spline to get the xsec value.
       else
         n = tidHi - tidLo + 1
         call spline_init(xVals, yVals, n)
         call spline_read(burnupVal, crossSection_cdbl)
         crossSection = real(crossSection_cdbl,kind=rkind)
         call spline_free
       endif

    end function get_cross_section

    !> function find_group_history returns the group history index in the xml 
    !> lattice physics input file for a given burnup time and group ID.
    !! 
    !! @param materialID     :
    !! @param branchID       :The applicable branch (ie, state and burnup value)
    !! @param timeID         :The applicable burnup step
    !! @param groupID        :The applicable energy group
    integer(kind=ikind) function find_group_history(materialID, branchID, &
      timeID, groupID) result(groupHistory)

        integer(kind=ikind),intent(in) :: materialID, branchID, timeID, groupID

        integer(kind=ikind) :: i
        logical :: match

        i = 1
        match = .false.
        do while (.not.match)
          if ((timeID .eq. xmlfile%material(materialID)%branch(branchID)%&
              grouphistory(i)%tid) .and. &
              (groupID .eq. xmlfile%material(materialID)%branch(branchID)% &
              grouphistory(i)%gid)) then
            match = .true.
          else
            i = i + 1
          endif
        enddo

        groupHistory = i

    end function find_group_history

    !> function find_cross_section_id returns the cross section index in the 
    !> xml lattice physics input file for a given cross section type, passed 
    !> as a Fortran string.  The scatter cross section requires the number of a
    !> "to" group
    !! 
    !! @param materialID       :
    !! @param branchID         :The applicable branch (ie, state and burnup step)
    !! @param groupHistory     :The applicable group history index in xml file
    !! @param crossSectionType :
    !! @param toGroupNum       :Group scattered to (scatter xsecs only)
    !!
    integer(kind=ikind) function find_cross_section_id(materialID, branchID, &
      groupHistory, crossSectionType, toGroupNum) result(crossSectionID)

        integer(kind=ikind), intent(in) :: materialID, branchID, groupHistory
        integer(kind=ikind), intent(in), optional :: toGroupNum
        character(kind=ckind, len=charlen), intent(in) :: crossSectionType

        integer(kind=ikind) :: i
        logical :: match

        i = 1
        match = .false.
        do while (.not.match .and. i .le. size(xmlfile%material(materialID)%&
          branch(branchID)%grouphistory(groupHistory)%macrocrosssection))
          if (xmlfile%material(materialID)%branch(branchID)%grouphistory&
            (groupHistory)%macrocrosssection(i)%typeof .eq. crossSectionType)&
            then
            if (crossSectionType .ne. "scatter") then
              match = .true.
            elseif (xmlfile%material(materialID)%branch(branchID)%grouphistory&
              (groupHistory)%macrocrosssection(i)%togid .eq. toGroupNum) then
              match = .true.
            else
              i = i + 1
            endif
          else
            i = i + 1
          endif
        enddo
        
        if (.not.match) then
          if (crossSectionType .eq. "neutronsper") then
            match = .true.
          else
            i = i + 1
          endif
        endif

       crossSectionID = i

    end function find_cross_section_id

end module get_cross_section_mod
