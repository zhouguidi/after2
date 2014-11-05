module cdiio
    use util
    use configure
    implicit none

    type FileInfoI
        character(1000) :: fileName
        integer :: streamID, vlistID, gridID, zaxisID, taxisID
        integer :: nvars, nx, ny, nz, nt
        integer :: tID, uID, vID, wID, itaID, gphID
    end type FileInfoI

    type FileInfoO
        character(1000) :: fileName
        integer :: streamID, vlistID, gridID, zaxisID, taxisID
    end type FileInfoO
contains
    function queryInput(conf) result(fi)
        type(Config), intent(inout) :: conf
        type(FileInfoI) :: fi

        integer :: ivar, stat, nout, nnout, id
        integer, dimension(6) :: ids
        character(100) :: varname

        fi%fileName = conf%iFileName
        fi%streamID = streamOpenRead(trim(fi%fileName))
        if (fi%streamID .lt. 0) call errorHandler(cdiStringError(fi%streamID))
        fi%vlistID = streamInqVlist(fi%streamID)

        ! find variables
        fi%nvars = vlistNvars(fi%vlistID)
        fi%tID = -1
        fi%uID = -1
        fi%vID = -1
        fi%wID = -1
        fi%itaID = -1
        fi%gphID = -1

        do ivar = 0, fi%nvars - 1
            call vlistInqVarName(fi%vlistID, ivar, varname)
            select case (trim(adjustl(varname)))
                case ("t")
                    fi%tID = ivar
                case ("u")
                    fi%uID = ivar
                case ("v")
                    fi%vID = ivar
                case ("omega")
                    fi%wID = ivar
                case ("svo")
                    fi%itaID = ivar
                case ("geopoth")
                    fi%gphID = ivar
            end select
        end do

        nout = conf%varOutCount()
        ! delete some output variables if required input variables are not present
        if (fi%tID .eq. -1) call conf%disableVarsOn("t")
        if (fi%uID .eq. -1) call conf%disableVarsOn("u")
        if (fi%vID .eq. -1) call conf%disableVarsOn("v")
        if (fi%wID .eq. -1) call conf%disableVarsOn("omega")
        if (fi%itaID .eq. -1) call conf%disableVarsOn("svo")
        if (fi%gphID .eq. -1) call conf%disableVarsOn("geopoth")

        ids = [fi%tID, fi%uID, fi%vID, fi%wID, fi%itaID, fi%gphID]
        id = 0
        do ivar = 1, size(ids, 1)
            if (ids(ivar) /= -1) then
                id = ivar
                exit
            endif
        enddo
            
        if (id == 0) &
            call errorHandler("input file contains no required variables")

        id = ids(id)
        ! get axis IDs and length
        fi%gridID = vlistInqVarGrid(fi%vlistID, id)
        fi%zaxisID = vlistInqVarZaxis(fi%vlistID, id)
        fi%taxisID = vlistInqTaxis(fi%vlistID)

        fi%nx = gridInqXsize(fi%gridID)
        fi%ny = gridInqYsize(fi%gridID)
        fi%nz = zaxisInqSize(fi%zaxisID)
        fi%nt = 0
        do
            stat = streamInqTimestep(fi%streamID, fi%nt)
            if (stat .eq. 0) then
                exit
            endif
            fi%nt = fi%nt + 1
        enddo

        ! if only one time step then enable average mode
        if (.not. conf%aveMode) then
            if (fi%nt .eq. 1) then
                if(.not. conf%silent) write(*, '(A)')"after2: time axis length is 1, average mode enabled..."
                conf%aveMode = .true.
            endif
        endif

        ! if average mode delete some output variable because they use turbulence data
        if (conf%aveMode) then
            conf%chunks = fi%nt
            call conf%disableVarsAveMode(.false.)
        endif

        ! if some output variables are removed from the list
        nnout = conf%varOutCount()
        if (nout /= nnout) then
            if (nnout == 0) then ! if no output vars are left
                call errorHandler("none of the desired output variables can be calculated due to lack of data in the input file")
            elseif (.not. conf%silent) then
                write(*, '(A)')"Warning(after2): some desired output variables disabled due to lack of data in the input file"
            endif
        endif
    end function queryInput

    function initOutput(conf, fi) result(fo)
        type(Config), intent(in) :: conf
        type(FileInfoI), intent(in) :: fi
        type(FileInfoO) :: fo

        character(1000) :: msg
        integer :: stat, ct, cl

        fo%fileName = conf%ofileName
        fo%vlistID = vlistCreate()
        fo%gridID = fi%gridID
        fo%zaxisID = fi%zaxisID
        call defOVars(fo, conf) ! define variables

        ! taxis setup
        if (conf%tType .eq. -1) then ! determin from input
            fo%taxisID = taxisCreate(taxisInqType(fi%taxisID))
            if (taxisInqType(fi%taxisID) .eq. TAXIS_RELATIVE) then
                call taxisDefRdate(fo%taxisID, taxisInqRdate(fi%taxisID))
                call taxisDefRtime(fo%taxisID, taxisInqRtime(fi%taxisID))
                call taxisDefCalendar(fo%taxisID, taxisInqCalendar(fi%taxisID))
            endif
        else
            fo%taxisID = taxisCreate(conf%tType)
            call taxisDefRdate(fo%taxisID, conf%refDate)
            call taxisDefRtime(fo%taxisID, conf%refTime)
            call taxisDefCalendar(fo%taxisID, conf%calendar)
            call taxisDefTunit(fo%taxisID, conf%tUnit)
        endif
        call vlistDefTaxis(fo%vlistID, fo%taxisID)

        msg = "Another Afterburner version " // trim(VERSION) // " (zhouguidi@gmail.com)"
        stat = vlistDefAttTxt(fo%vlistID, CDI_GLOBAL, "AFTER2", len_trim(msg), msg)
        msg = "Helmholtz-Center for Ocean Research"
        stat = vlistDefAttTxt(fo%vlistID, CDI_GLOBAL, "Institute", len_trim(msg), msg)

        fo%streamID = streamOpenWrite(trim(fo%fileName), conf%filetype)
        if (fo%streamID .lt. 0) call errorHandler(cdiStringError(fo%streamID))
        if (conf%comptype .ne. COMPRESS_NONE) then
            if (conf%comptype == -1) then ! determin from input
               ct = streamInqCompType(fi%streamID)
               cl = streamInqCompLevel(fi%streamID)
           else
               ct = conf%comptype
               cl = conf%complev
            endif
            call streamDefCompType(fo%streamID, ct)
            call streamDefCompLevel(fo%streamID, cl)
        endif
        call streamDefVlist(fo%streamID, fo%vlistID)
    end function initOutput

    subroutine defOVars(fo, conf)
        type(FileInfoO), intent(inout) :: fo
        type(Config), intent(in) :: conf

        integer :: varOutID, i

        do i= 1, conf%nvars
            if (.not. conf%output(i)) cycle

            varOutID = vlistDefVar(fo%vlistID, fo%gridID, fo%zaxisID, TIME_VARIABLE)
            call conf%setVarOutID(trim(conf%variables(i)%name), varOutID)

            call vlistDefVarCode(fo%vlistID, varOutID, i + 300)
            call vlistDefVarName(fo%vlistID, varOutID, trim(conf%variables(i)%name))
            call vlistDefVarlongName(fo%vlistID, varOutID, trim(conf%variables(i)%longName))
            call vlistDefVarstdName(fo%vlistID, varOutID, trim(conf%variables(i)%stdName))
            call vlistDefVarUnits(fo%vlistID, varOutID, trim(conf%variables(i)%units))
            call vlistDefVarDataType(fo%vlistID, varOutID, DATATYPE_FLT32)
        enddo
    end subroutine defOVars

    subroutine cdiFinish(fi, fo)
        type(FileInfoI), intent(in) :: fi
        type(FileInfoO), intent(in) :: fo

        call streamClose(fi%streamID)
        call streamClose(fo%streamID)
        call vlistDestroy(fo%vlistID)
        call taxisDestroy(fo%taxisID)
    end subroutine cdiFinish
end module cdiio
