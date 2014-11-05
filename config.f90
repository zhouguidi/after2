module configure
    use iso_fortran_env
    use util
    implicit none
    include "cdi.inc"

    !type definations
    type Variable
        character(100) :: name, stdName, longName, units, dependOn, dependOnInput, depended
        logical :: aveMode
    end type Variable

    type Config
        character(1000) :: iFileName, oFileName
        integer :: filetype, comptype, complev, chunks
        integer :: tType, refDate, refTime, tUnit, calendar
        logical :: silent, aveMode, tflux

        type(Variable), dimension(:), allocatable :: variables
        logical, dimension(:), allocatable :: compute, output
        integer, dimension(:), allocatable :: varOutID
        integer :: nout, nvars
    contains
        procedure :: getVarOutID => get_var_out_id
        procedure :: setVarOutID => set_var_out_id
        procedure :: getVarID => get_Var_ID
        procedure :: isVarCompute => is_Var_Compute
        procedure :: isVarOutput => is_Var_Output
        procedure :: mkDependTree => make_depended
        procedure :: disableAllVars => disable_all_vars
        procedure :: enableAllVars => enable_all_vars
        procedure :: disableVar => disable_var
        procedure :: enableVar => enable_var
        procedure :: disableVarsOn => disable_vars_on
        procedure :: disableVarsAveMode => disable_vars_ave_mode
        procedure :: varOutCount => var_out_count
        procedure :: printVarSummary => print_var_summary
        final :: clean_Config
    end type Config

    ! constructor
    interface Config
        module procedure init_Config
    end interface Config
contains
    function init_Config() result(conf)
        type(Config) :: conf

        character(1000) :: VARDEF_PATH
        integer :: stat

        call get_environment_variable("AFTER2_VARDEF_PATH", VARDEF_PATH, status = stat)
        if (stat /= 0) call errorHandler("unable to get environment variable AFTER2_VARDEF_PATH")

        ! initialize predefined variables
        conf%variables = readVarDef(trim(adjustl(VARDEF_PATH)) // "/after2.vardef")
        conf%nvars = size(conf%variables, 1)
        conf%nout = conf%nvars
        call conf%mkDependTree
        allocate(conf%compute(conf%nvars), conf%output(conf%nvars), conf%varOutID(conf%nvars))
        conf%compute = .true.
        conf%output = .true.
        conf%varOutID = -1

        ! default values
        conf%filetype = FILETYPE_NC4C
        conf%comptype = -1 ! determin from input
        conf%complev = 6
        conf%chunks = 1
        conf%silent = .false.
        conf%tflux = .false.
        conf%tType = -1 ! default is to determin from input
        conf%refDate = 19700101
        conf%refTime = 0
        conf%tUnit = TUNIT_DAY
        conf%calendar = CALENDAR_PROLEPTIC
        conf%aveMode = .false.
    end function init_Config

    subroutine clean_Config(this)
        type(Config) :: this

        deallocate(this%compute, this%output, this%variables, this%varOutID)
    end subroutine clean_Config

    subroutine print_var_summary(this)
        class(Config) :: this

        integer :: i, n
        
        n = this%nvars
        do i = 1, n
            write(*, '(A10, L2, L2, "|", A20, "|", A60, "|", A30)')trim(this%variables(i)%name), this%compute(i), this%output(i), &
                trim(this%variables(i)%dependOnInput),trim(this%variables(i)%dependOn), trim(this%variables(i)%depended)
        enddo
    end subroutine print_var_summary

    function get_Var_ID(this, varname) result(id)
        class(Config) :: this
        character(*), intent(in) :: varname
        integer :: id

        integer :: i

        id = 0
        do i = 1, this%nvars
            if (trim(adjustl(this%variables(i)%name)) == trim(adjustl(varname))) then
                id = i
                exit
            endif
        enddo
    end function get_Var_ID

    function is_Var_Compute(this, varname) result(comp)
        class(Config) :: this
        character(*), intent(in) :: varname
        logical :: comp

        integer :: ind

        ind = this%getVarID(varname)
        if (ind == 0) then
            call errorHandler("variable " // trim(varname) // " undefined")
        else
            comp = this%compute(ind)
        endif
    end function is_Var_Compute

    function is_Var_Output(this, varname) result(outp)
        class(Config) :: this
        character(*), intent(in) :: varname
        logical :: outp

        integer :: ind

        ind = this%getVarID(varname)
        if (ind == 0) then
            call errorHandler("variable " // trim(varname) // " undefined")
        else
            outp = this%output(ind)
        endif
    end function is_Var_Output

    subroutine disable_all_vars(this)
        class(Config) :: this

        this%compute = .false.
        this%output = .false.
        this%nout = 0
    end subroutine disable_all_vars

    subroutine enable_all_vars(this)
        class(Config) :: this

        this%compute = .true.
        this%output = .true.
        this%nout = count(this%output)
    end subroutine enable_all_vars

    ! enable a variable to output
    !
    ! user call this to enable a variable to output, "root" omitted, default to true
    ! internally call this to track dependOnency, "root" set to false
    ! always enable compute, if "root", enable output
    recursive subroutine enable_var(this, varname, root)
        class(Config) :: this
        character(*), intent(in) :: varname
        logical, optional, intent(in) :: root

        character(100) :: dependOn
        character(:), allocatable :: tok, rem
        integer :: ind, stat
        logical :: last, root_

        root_ = .true.
        if (present(root)) root_ = root

        ind = this%getVarID(varname)
        if (ind == 0) then
            call errorHandler("variable " // trim(varname) // " undefined")
        else
            this%compute(ind) = .true.
            if (root_) this%output(ind) = .true.

            dependOn = this%variables(ind)%dependOn
            do
                call strtok(dependOn, tok, rem, last = last, stat = stat, deli = ",", allowblk = .false.)
                if (stat /= 0) call errorHandler("invalid dependOnency list")

                if (trim(tok) /= '') call this%enableVar(trim(tok), .false.)

                dependOn = rem
                if (last) exit
            end do
        endif

        this%nout = count(this%output)
        deallocate(tok, rem)
    end subroutine enable_var

    function readVarDef(fn) result(vars)
        character(*), intent(in) :: fn
        type(Variable), dimension(:), allocatable :: vars

        integer :: stat 
        type(Variable) :: var

        open(unit = 188, file = trim(adjustl(fn)), status = 'old', action = 'read', iostat = stat)
        if (stat /= 0) then
            write(0, *)"trouble reading variable defination file"
            stop
        endif
        read(unit = 188, fmt = '(A)')
        do
            read(unit = 188, fmt = '(A)', iostat = stat)var%name
            read(unit = 188, fmt = '(A)', iostat = stat)var%dependOnInput
            read(unit = 188, fmt = '(A)', iostat = stat)
            read(unit = 188, fmt = '(L3)', iostat = stat)var%aveMode
            read(unit = 188, fmt = '(A)', iostat = stat)var%dependOn
            read(unit = 188, fmt = '(A)', iostat = stat)var%units
            read(unit = 188, fmt = '(A)', iostat = stat)var%stdName
            read(unit = 188, fmt = '(A)', iostat = stat)var%longName
            var%name = adjustl(var%name)
            var%name = var%name(2 : len_trim(var%name))
            var%dependOn = adjustl(var%dependOn)
            var%units = adjustl(var%units)
            var%stdName = adjustl(var%stdName)
            var%longName = adjustl(var%longName)
            var%depended = ""
            if (stat == -1) exit
            if (allocated(vars)) then
                vars = [vars, var]
            else
                allocate(vars(1))
                vars = [var]
            endif
        end do
        close(unit = 188)
    end function readVarDef

    subroutine make_depended(this)
        class(Config) :: this

        integer :: i, ind, stat
        character(100) :: thisvar, dependOn
        character(:), allocatable :: tok, rem
        logical :: last

        do i = 1, this%nvars
            thisvar = this%variables(i)%name
            dependOn = this%variables(i)%dependOn
            if (dependOn == "") cycle
            do
                call strtok(dependOn, tok, rem, last = last, stat = stat, deli = ',', allowblk = .false.)
                if (stat /= 0) call errorHandler("invalid dependency list")

                ind = this%getVarID(tok)
                if (ind == 0) call errorHandler("variable " // trim(tok) // " undefined")

                if (this%variables(ind)%depended == "") then
                    this%variables(ind)%depended = trim(thisvar)
                else
                    if (index(this%variables(ind)%depended, trim(thisvar)) /= 0) then
                        call errorHandler("something wrong happend")
                    else
                        this%variables(ind)%depended = trim(this%variables(ind)%depended) // "," // trim(thisvar)
                    endif
                endif
                dependOn = rem
                if (last) exit
            enddo
        end do

        deallocate(tok, rem)
    end subroutine make_depended

    ! disable a variable
    !
    ! only if this var is not required by another var which is to be computed
    ! if this var is disabled, check if the vars it depends on can also be disabled (only if that var is not output)
    recursive subroutine disable_var(this, varname)
        class(Config) :: this
        character(*), intent(in) :: varname

        integer :: thisind, ind, stat
        character(100) :: dependOn, depended
        character(:), allocatable :: tok, rem
        logical :: last, compute

        thisind = this%getVarID(varname)
        if (thisind == 0) then
            call errorHandler("variable " // trim(varname) // " undefined")
        else
            ! always disable output for this var
            this%output(thisind) = .false.

            ! if all vars depending on this var are disabled from computing, then so is this var
            compute = .false.
            depended = this%variables(thisind)%depended
            if (depended /= "") then
                do
                    call strtok(depended, tok, rem, last = last, stat = stat, deli = ",", allowblk = .false.)
                    if (stat /= 0) call errorHandler("invalid dependency list")

                    ind = this%getVarID(tok)
                    if (ind == 0) call errorHandler("variable " // trim(tok) // " undefined")

                    if (this%compute(ind)) then
                        compute = .true.
                        exit
                    endif

                    depended = rem
                    if (last) exit
                enddo
            endif
            if (.not. compute) then
                this%compute(thisind) = .false.

                ! if this var is disabled, check if its dependency can also be disabled (only when not output)
                dependOn = this%variables(thisind)%dependOn
                if (dependOn /= "") then
                    do
                        call strtok(dependOn, tok, rem, last = last, stat = stat, deli = ",", allowblk = .false.)
                        if (stat /= 0) call errorHandler("invalid dependency list")

                        ind = this%getVarID(tok)
                        if (ind == 0) call errorHandler("variable " // trim(tok) // " undefined")

                        ! the dependency var is not output
                        if (.not. this%output(ind)) then
                            call this%disableVar(tok)
                        endif

                        dependOn = rem
                        if (last) exit
                    enddo
                endif
            endif
        endif
        
        this%nout = count(this%output)
        deallocate(tok, rem)
    end subroutine disable_var

    subroutine disable_Vars_On(this, inpVarName)
        class(Config) :: this
        character(*), intent(in) :: inpVarName

        integer :: i, ind, stat
        character(100) :: dependOnInput
        character(:), allocatable :: tok, rem
        logical :: last, related

        do i = 1, this%nvars
            dependOnInput = this%variables(i)%dependOnInput
            related = .false.
            do
                call strtok(dependOnInput, tok, rem, last = last, stat = stat, deli = ",", allowblk = .false.)
                if (stat /= 0) call errorHandler("invalid dependency list")

                if (trim(tok) == trim(inpVarName)) then
                    related = .true.
                    exit
                endif

                dependOnInput = rem
                if (last) exit
            enddo

            if (related) call this%disableVar(trim(this%variables(i)%name))
        enddo

        this%nout = count(this%output)
        deallocate(tok, rem)
    end subroutine disable_Vars_on

    subroutine disable_vars_ave_mode(this, mode)
        class(Config) :: this
        logical, intent(in) :: mode

        integer :: i

        do i = 1, this%nvars
            if (this%variables(i)%aveMode .eqv. mode) call this%disableVar(trim(this%variables(i)%name))
        enddo
        this%nout = count(this%output)
    end subroutine disable_vars_ave_mode

    function var_out_count(this) result(cnt)
        class(Config) :: this
        integer :: cnt

        cnt = this%nout
    end function var_out_count

    function get_var_out_id(this, varname) result(id)
        class(Config) :: this
        character(*), intent(in) :: varname
        integer :: id

        integer :: ind

        ind = this%getVarID(varname)
        if (ind == 0) then
            call errorHandler("variable " // trim(varname) // " undefined")
        else
            id = this%varOutID(ind)
        endif
    end function get_var_out_id

    subroutine set_var_out_id(this, varname, id)
        class(Config) :: this
        character(*), intent(in) :: varname
        integer, intent(in) :: id

        integer :: ind

        ind = this%getVarID(varname)
        if (ind == 0) then
            call errorHandler("variable " // trim(varname) // " undefined")
        else
            this%varOutID(ind) = id
        endif
    end subroutine set_var_out_id
end module configure
