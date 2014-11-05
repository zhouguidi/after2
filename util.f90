module util
    implicit none

    character(*), parameter :: VERSION = '1.3.1', VERSION_DATE = '2013-05-24'

    ! physical constants
    real(8), parameter :: pi = 3.141592653 ! pi
    real(8), parameter :: RE = 6371000 ! earth radius (m)
    real(8), parameter :: p0 = 101300 ! reference pressure (Pa)
    real(8), parameter :: R = 286.9 ! gas constant (J Kg-1 K-1)
    real(8), parameter :: cp = 1006 ! specific heat capacity (J Kg-1 K-1)
    real(8), parameter :: k = R / cp
    real(8), parameter :: radpi = pi / 180
    real(8), parameter :: REpi = radpi * RE
    integer, parameter :: spd = 86400 ! seconds per day
    real(8), parameter :: Omega = 7.2921e-5 * 2 ! 2 * earth rotate (1/s)
    real(8), parameter :: g = 9.80665 ! gravity accelaration (m/s2)

contains
    subroutine errorHandler(msg)
        character(*), intent(in) :: msg

        write(0, '(A, A)')"ERORR(after2): ", trim(adjustl(msg))
        stop
    end subroutine errorHandler

    subroutine memoryError
        call errorHandler("out of memory")
    end subroutine memoryError

    subroutine strtok(str, tok, rem, last, stat, deli, allowblk)
        character(*), intent(in) :: str
        character(:), allocatable, intent(out) :: tok, rem
        logical, optional, intent(out) :: last
        integer, optional, intent(out) :: stat
        character(1), optional, intent(in) :: deli
        logical, optional, intent(in) :: allowblk

        character(1) :: deli_
        logical :: allowblk_
        integer :: ind

        if (present(stat)) stat = 0
        rem = str ! arg is trimed str
        deli_ = " "
        if (present(deli)) deli_ = deli
        allowblk_ = .true.
        if (present(allowblk)) allowblk_ = allowblk

        if (.not. allocated(tok)) allocate(character(len_trim(str)) :: tok)
        if (.not. allocated(rem)) allocate(character(len_trim(str)) :: rem)

        ind = index(trim(rem), deli_)
        if (ind .eq. 0) then ! no commas found: stop
            tok = rem
            rem = ""
        elseif (ind .eq. 1) then ! comma at pos 1
            if (allowblk_) then ! allowed
                tok = ""
                rem = rem(2 : len_trim(rem))
            else
                if (present(stat)) then
                    stat = 1
                else
                    call errorHandler("blank token detected")
                endif
            endif
        elseif (ind .eq. len_trim(rem)) then ! comma at the end
            if (allowblk_) then
                tok = ""
                rem = ""
            else
                if (present(stat)) then
                    stat = 2
                else
                    call errorHandler("blank token detected")
                endif
            endif
        else !good comma
            tok = rem(1 : ind - 1)
            rem = rem(ind + 1 : len_trim(rem))
        endif

        if (present(last)) then
            if (len_trim(rem) .eq. 0) then
                last = .true.
            else
                last = .false.
            endif
        endif
    end subroutine strtok

    elemental function julday(year, month, day, hour, minute, second)
        integer, intent(in) :: year, month, day
        integer, optional, intent(in) :: hour, minute
        real(8), optional, intent(in) :: second
        real(8) :: julday 
        integer, parameter :: igreg = 15 + 31 * (10 + 12 * 1582)
        integer :: ja, jm, jy
        integer :: hour_, minute_
        real(8) :: second_

        hour_ = 0
        minute_ = 0
        second_ = 0
        if (present(hour)) hour_ = hour
        if (present(minute)) minute_ = minute
        if (present(second)) second_ = second

        jy = year
        if (jy .eq. 0) then
            julday = -99999
            return
        endif
        if (month .lt. 1 .or. month .gt. 12) then
            julday = -99998
            return
        endif
        if (day .lt. 1 .or. day .gt. 31) then
            julday = -99997
            return
        endif
        if (hour .lt. 0 .or. hour .gt. 23) then
            julday = -99996
            return
        endif
        if (minute .lt. 0 .or. minute .gt. 59) then
            julday = -99995
            return
        endif
        if (second .lt. 0 .or. second .ge. 60) then
            julday = -99994
            return
        endif

        if (jy .lt. 0) jy = jy + 1
        if (month .gt. 2) then
            jm = month + 1
        else
            jy = jy - 1
            jm = month + 13
        endif

        julday = int(365.25 * jy) + int(30.6001 * jm) + day + 1720995
        if (day + 31 * (month + 12 * year) .ge. igreg) then
            ja = int(0.01 * jy)
            julday = julday + 2 - ja + int(0.25 * ja)
        end if

        julday = julday + (hour + (minute + second_ / 60) / 60) / 24
    end function julday

    elemental subroutine decodeDate(date, year, month, day)
        integer, intent(in) :: date
        integer, intent(out) :: year, month, day

        year = date / 10000
        month = (date - year * 10000) / 100
        day = date - year * 10000 - month * 100
    end subroutine decodeDate

    elemental subroutine decodeTime(time, hour, minute, second)
        integer, intent(in) :: time
        integer, intent(out) :: hour, minute, second

        hour = time / 10000
        minute = (time - hour * 10000) / 100
        second = time - hour * 10000 - minute * 100
    end subroutine decodeTime

    subroutine printMatrix(varname, var, ln, cl, fm, fn)
        !printMatrix: print a Matrix to a file
        
        integer, intent(in) :: ln, cl
        real(8), dimension(ln, cl), intent(in) :: var
        character(*), intent(in) :: varname, fm
        character(*), intent(in), optional :: fn

        integer :: i, j, fid

        fid = 6
        if (present(fn)) then
            fid = 101
            open(unit = fid, file = trim(fn), status = 'replace')
        endif

        !write(unit = fid, fmt = '(A,A,I3,A,I3,A)')varname, " (", ln, "x", cl, ")"
        do i = 1, ln
            do j = 1, cl - 1
                write(unit = fid, fmt = fm, advance = 'no')var(i, j)
            enddo
            write(unit = fid, fmt = fm)var(i, cl)
        enddo

        if (present(fn)) then
            close(fid)
        endif
    end subroutine printMatrix
end module util
