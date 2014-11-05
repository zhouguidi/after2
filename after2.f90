program after2
    use iso_fortran_env
    use util
    use configure
    use cdiio
    implicit none

    ! local variables
    type(Config) :: conf
    type(FileInfoI) :: fi
    type(FileInfoO) :: fo

    ! axis variables
    real(8), dimension(:), allocatable :: lon, lat, p, pk, pk1, pki, lnp

    ! physical variables
    real(8), dimension(:, :), allocatable :: mt, mu, mv, mw, mtheta, mita, mtflxu, mtflxv, &
        mdivxtflxu, mdivytflxv, mdivtflxh, mtflxw, mdivtflxv, mtt, mtranstu, mtranstv, mtranstw, &
        mdiabh, mvflxu, mvflxv, mdivxvflxu, mdivyvflxv, mdivvflxh, mvflxw, mdivvflxv, mtransvu, &
        mtransvv, mtransvw, mdudz, mdvdz, mdtdz, mbvf, megr, mgph, mgpht, stdgph, mmflxuu, mmflxuv, &
        mmflxvv, mmflxuw, mmflxvw, mmflxww, mdivxmflxuu, mdivxmflxuv, mdivymflxuv, mdivymflxvv, &
        mdivmflxh, meke

    ! temps
    integer :: ivar, ilev, stat, nmiss
    integer :: chunklen, ichk, vdate, vtime, tsstart, tsend
    real(8) :: timebegin, timeend
    real(8) :: dlnp
    real(8), dimension(:, :), allocatable :: x, y, dx, dy, dgph, f
    real(8), dimension(:, :), allocatable :: pmtflxw, pmvflxw, pmt, pmtheta, pmita, pmu, pmv, pmgph

    ! configuration
    allocate(conf%variables(1), conf%compute(1), conf%output(1), conf%varOutID(1))
    conf = Config()
    call perseArgs(conf)
    if (conf%aveMode) then
        if (.not. conf%silent) write(*, '(A)')"after2: average mode enabled..."
    endif
    fi = queryInput(conf)

    ! start time measurement
    call cpu_time(timebegin)

    ! allocate axis variables and get value
    allocate(lon(fi%nx), &
             lat(fi%ny), &
             p(fi%nz), &
             pk(fi%nz), &
             pk1(fi%nz), &
             pki(fi%nz), &
             lnp(fi%nz), &
             f(fi%nx, fi%ny), &
             stat = stat)
    if (stat .ne. 0) call memoryError
    
    if (gridInqXVals(fi%gridID,  lon) .ne. fi%nx) call errorHandler("error reading longitude axis.")
    if (gridInqYVals(fi%gridID,  lat) .ne. fi%ny) call errorHandler("error reading latitude axis.")
    call zaxisInqLevels(fi%zaxisID,  p)

    pk = (p / p0) ** k
    pk1 = pk / p
    pki = 1 / pk
    lnp = log(p)
    f = 2 * Omega * reshape(spread(sin(radpi * lat), 1, fi%nx), [fi%nx, fi%ny])

    ! compute spatial distance and f and beta
    allocate(x(fi%nx, fi%ny), &
             y(fi%nx, fi%ny), &
             dx(fi%nx, fi%ny), &
             dy(fi%nx, fi%ny), &
             stat = stat)
    if (stat .ne. 0) call memoryError
    where (lon .lt. 0)
        lon = lon + 360
    endwhere
    x = spread(lon, 1, fi%ny)
    x = transpose(x)
    dx(2 : fi%nx - 1, :) = x(3 : fi%nx, :) - x(1 : fi%nx - 2, :)
    dx(1, :) = x(2, :) - x(fi%nx, :)
    dx(fi%nx, :) = x(1, :) - x(fi%nx - 1, :)
    where(abs(dx) > 180)
        dx = 360 + dx
    endwhere
    dx = dx * REpi * reshape(spread(cos(lat * radpi), 1, fi%nx), [fi%nx, fi%ny])
    y = reshape(spread(lat, 1, fi%nx), [fi%nx, fi%ny])
    dy(:, 2 : fi%ny - 1) = y(:, 3 : fi%ny) - y(:, 1 : fi%ny - 2)
    dy(:, 1) = y(:, 2) - y(:, 1)
    dy(:, fi%ny) = y(:, fi%ny) - y(:, fi%ny - 1)
    dy = dy  * REpi

    ! open output file and initialize
    fo = initOutput(conf, fi)

    ! allocate memory
    if (conf%isVarCompute("mt")) allocate(mt(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mu")) allocate(mu(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mv")) allocate(mv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mw")) allocate(mw(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mtheta")) allocate(mtheta(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mita")) allocate(mita(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mtflxu")) allocate(mtflxu(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mtflxv")) allocate(mtflxv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivxtflxu")) allocate(mdivxtflxu(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivytflxv")) allocate(mdivytflxv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivtflxh")) allocate(mdivtflxh(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mtflxw")) allocate(mtflxw(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivtflxv")) allocate(mdivtflxv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mtt")) allocate(mtt(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mtranstu")) allocate(mtranstu(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mtranstv")) allocate(mtranstv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mtranstw")) allocate(mtranstw(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdiabh")) allocate(mdiabh(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mvflxu")) allocate(mvflxu(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mvflxv")) allocate(mvflxv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivxvflxu")) allocate(mdivxvflxu(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivyvflxv")) allocate(mdivyvflxv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivvflxh")) allocate(mdivvflxh(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mvflxw")) allocate(mvflxw(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivvflxv")) allocate(mdivvflxv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mtransvu")) allocate(mtransvu(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mtransvv")) allocate(mtransvv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mtransvw")) allocate(mtransvw(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdudz")) allocate(mdudz(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdvdz")) allocate(mdvdz(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdtdz")) allocate(mdtdz(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mbvf")) allocate(mbvf(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("megr")) allocate(megr(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mgph")) allocate(mgph(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mgpht")) allocate(mgpht(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("stdgph")) allocate(stdgph(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivtflxv")) allocate(pmtflxw(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivvflxv")) allocate(pmvflxw(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mtranstw")) allocate(pmt(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdtdz")) allocate(pmtheta(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mtransvw")) allocate(pmita(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdudz")) allocate(pmu(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdvdz")) allocate(pmv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdudz") .or. conf%isVarCompute("mdvdz") .or. conf%isVarCompute("mdtdz")) then
        allocate(pmgph(fi%nx, fi%ny), stat = stat)
        if (stat .ne. 0) call memoryError
        allocate(dgph(fi%nx, fi%ny), stat = stat)
        if (stat .ne. 0) call memoryError
    endif
    if (conf%isVarCompute("mmflxuu")) allocate(mmflxuu(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mmflxuv")) allocate(mmflxuv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mmflxvv")) allocate(mmflxvv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mmflxuw")) allocate(mmflxuw(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mmflxvw")) allocate(mmflxvw(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mmflxww")) allocate(mmflxww(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivxmflxuu")) allocate(mdivxmflxuu(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivxmflxuv")) allocate(mdivxmflxuv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivymflxuv")) allocate(mdivymflxuv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivymflxvv")) allocate(mdivymflxvv(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("mdivmflxh")) allocate(mdivmflxh(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError
    if (conf%isVarCompute("meke")) allocate(meke(fi%nx, fi%ny), stat = stat)
    if (stat .ne. 0) call memoryError

    ! compute!
    chunklen = fi%nt / conf%chunks
    do ichk = 1, conf%chunks
        tsstart = (ichk - 1) * chunklen
        tsend = ichk * chunklen - 1
        if (ichk .eq. conf%chunks) tsend = fi%nt - 1

        stat = streamInqTimestep(fi%streamID, tsend)
        vdate = taxisInqVdate(fi%taxisID)
        vtime = taxisInqVtime(fi%taxisID)
        call taxisDefVdate(fo%taxisID, vdate)
        call taxisDefVtime(fo%taxisID, vtime)
        stat = streamDefTimestep(fo%streamID, ichk - 1)

        do ilev = fi%nz - 1, 0, -1
            ! compute space-unrelated vars
            call getLevData(conf, fi, pki, tsstart, tsend, ilev, mt, mtheta, mu, mv, mw, mtt, &
                mtflxu, mtflxv, mtflxw, mita, mvflxu, mvflxv, mvflxw, mgph, mgpht, stdgph, &
                mmflxuu, mmflxuv, mmflxvv, mmflxuw, mmflxvw, mmflxww, meke)

            ! output space-unrelated vars
            if (conf%isVarOutput("mt")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mt"), ilev, mt, 0)
            if (conf%isVarOutput("mtheta")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mtheta"), ilev, mtheta, 0)
            if (conf%isVarOutput("mu")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mu"), ilev, mu, 0)
            if (conf%isVarOutput("mv")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mv"), ilev, mv, 0)
            if (conf%isVarOutput("mw")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mw"), ilev, mw, 0)
            if (conf%isVarOutput("mtt")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mtt"), ilev, mtt, 0)
            if (conf%isVarOutput("mtflxu")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mtflxu"), ilev, mtflxu, 0)
            if (conf%isVarOutput("mtflxv")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mtflxv"), ilev, mtflxv, 0)
            if (conf%isVarOutput("mtflxw")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mtflxw"), ilev, mtflxw, 0)
            if (conf%isVarOutput("mita")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mita"), ilev, mita, 0)
            if (conf%isVarOutput("mvflxu")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mvflxu"), ilev, mvflxu, 0)
            if (conf%isVarOutput("mvflxv")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mvflxv"), ilev, mvflxv, 0)
            if (conf%isVarOutput("mvflxw")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mvflxw"), ilev, mvflxw, 0)
            if (conf%isVarOutput("mgph")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mgph"), ilev, mgph, 0)
            if (conf%isVarOutput("mgpht")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mgpht"), ilev, mgpht, 0)
            if (conf%isVarOutput("stdgph")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("stdgph"), ilev, stdgph, 0)
            if (conf%isVarOutput("mmflxuu")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mmflxuu"), ilev, mmflxuu, 0)
            if (conf%isVarOutput("mmflxuv")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mmflxuv"), ilev, mmflxuv, 0)
            if (conf%isVarOutput("mmflxvv")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mmflxvv"), ilev, mmflxvv, 0)
            if (conf%isVarOutput("mmflxuw")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mmflxuw"), ilev, mmflxuw, 0)
            if (conf%isVarOutput("mmflxvw")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mmflxvw"), ilev, mmflxvw, 0)
            if (conf%isVarOutput("mmflxww")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mmflxww"), ilev, mmflxww, 0)
            if (conf%isVarOutput("meke")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("meke"), ilev, meke, 0)

            ! compute space-related, height-unrelated vars
            if (conf%isVarCompute("mdivxmflxuu")) mdivxmflxuu = ddx(mmflxuu, dx)
            if (conf%isVarCompute("mdivxmflxuv")) mdivxmflxuv = ddx(mmflxuv, dx)
            if (conf%isVarCompute("mdivymflxuv")) mdivymflxuv = ddy(mmflxuv, dy)
            if (conf%isVarCompute("mdivymflxvv")) mdivymflxvv = ddy(mmflxvv, dy)
            if (conf%isVarCompute("mdivxtflxu")) mdivxtflxu = ddx(mtflxu, dx)
            if (conf%isVarCompute("mdivytflxv")) mdivytflxv = ddy(mtflxv, dy)
            if (conf%isVarCompute("mdivxvflxu")) mdivxvflxu = ddx(mvflxu, dx)
            if (conf%isVarCompute("mdivyvflxv")) mdivyvflxv = ddy(mvflxv, dy)
            if (conf%isVarCompute("mtranstu")) mtranstu = mu * ddx(mT, dx)
            if (conf%isVarCompute("mtranstv")) mtranstv = mv * ddy(mT, dy)
            if (conf%isVarCompute("mdivtflxh")) mdivtflxh = mdivxtflxu + mdivytflxv
            if (conf%isVarCompute("mtransvu")) mtransvu = mu * ddx(mita, dx)
            if (conf%isVarCompute("mtransvv")) mtransvv = mv * ddy(mita, dy)
            if (conf%isVarCompute("mdivvflxh")) mdivvflxh = mdivxvflxu + mdivyvflxv
            if (conf%isVarCompute("mdivmflxh")) then
                mdivmflxh = sqrt(mu ** 2 + mv ** 2)
                mdivmflxh = mu / mdivmflxh * (mdivxmflxuu + mdivymflxuv) + &
                            mv / mdivmflxh * (mdivxmflxuv + mdivymflxvv)
            endif

            ! output space-related, height-unrelated vars
            if (conf%isVarOutput("mdivxmflxuu")) &
                call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdivxmflxuu"), ilev, mdivxmflxuu, 0)
            if (conf%isVarOutput("mdivxmflxuv")) &
                call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdivxmflxuv"), ilev, mdivxmflxuv, 0)
            if (conf%isVarOutput("mdivymflxuv")) &
                call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdivymflxuv"), ilev, mdivymflxuv, 0)
            if (conf%isVarOutput("mdivymflxvv")) &
                call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdivymflxvv"), ilev, mdivymflxvv, 0)
            if (conf%isVarOutput("mdivxtflxu")) &
                call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdivxtflxu"), ilev, mdivxtflxu, 0)
            if (conf%isVarOutput("mdivytflxv")) &
                call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdivytflxv"), ilev, mdivytflxv, 0)
            if (conf%isVarOutput("mdivxvflxu")) &
                call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdivxvflxu"), ilev, mdivxvflxu, 0)
            if (conf%isVarOutput("mdivyvflxv")) &
                call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdivyvflxv"), ilev, mdivyvflxv, 0)
            if (conf%isVarOutput("mtranstu")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mtranstu"), ilev, mtranstu, 0)
            if (conf%isVarOutput("mtranstv")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mtranstv"), ilev, mtranstv, 0)
            if (conf%isVarOutput("mdivtflxh")) &
                call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdivtflxh"), ilev, mdivtflxh, 0)
            if (conf%isVarOutput("mtransvu")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mtransvu"), ilev, mtransvu, 0)
            if (conf%isVarOutput("mtransvv")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mtransvv"), ilev, mtransvv, 0)
            if (conf%isVarOutput("mdivvflxh")) &
                call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdivvflxh"), ilev, mdivvflxh, 0)
            if (conf%isVarOutput("mdivmflxh")) &
                call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdivmflxh"), ilev, mdivmflxh, 0)

            if (ilev /= fi%nz - 1) then
                ! compute pressure-related vars
                dlnp = lnp(ilev + 1) - lnp(ilev + 2)
                if (conf%isVarCompute("mdivtflxv")) mdivtflxv = ddp(mtflxw - pmtflxw, p(ilev + 1), dlnp)
                if (conf%isVarCompute("mdivvflxv")) mdivvflxv = ddp(mvflxw - pmvflxw, p(ilev + 1), dlnp)
                if (conf%isVarCompute("mtranstw")) then
                    if (conf%tflux) then
                        mtranstw = mw * ddp(mt - pmt, p(ilev + 1), dlnp)
                    else
                        mtranstw = mw * ddp(mtheta - pmt, p(ilev + 1), dlnp)
                    endif
                endif
                if (conf%isVarCompute("mtransvw")) mtransvw = mw * ddp(mita - pmita, p(ilev + 1), dlnp)

                ! output pressure-related vars
                if (conf%isVarOutput("mdivtflxv")) &
                    call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdivtflxv"), ilev, mdivtflxv, 0)
                if (conf%isVarOutput("mdivvflxv")) &
                    call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdivvflxv"), ilev, mdivvflxv, 0)
                if (conf%isVarOutput("mtranstw")) &
                    call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mtranstw"), ilev, mtranstw, 0)
                if (conf%isVarOutput("mtransvw")) &
                    call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mtransvw"), ilev, mtransvw, 0)

                ! compute diabatic heating
                if (conf%isVarCompute("mdiabh")) mdiabh = mtt + mtranstu + mtranstv + pk(ilev + 1) * mtranstw + &
                    pk(ilev + 1) * (mdivtflxh + mdivtflxv)

                ! output diabatic heating
                if (conf%isVarOutput("mdiabh")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdiabh"), ilev, mdiabh, 0) 

                ! compute altitude-related vars
                dgph = mgph - pmgph
                if (conf%isVarCompute("mdudz")) mdudz = ddz(mu - pmu, dgph)
                if (conf%isVarCompute("mdvdz")) mdvdz = ddz(mv - pmv, dgph)
                if (conf%isVarCompute("mdtdz")) mdtdz = ddz(mtheta - pmtheta, dgph)
                if (conf%isVarCompute("mbvf")) mbvf = sqrt(g * ddz(mtheta - pmtheta, dgph) / mtheta)
                if (conf%isVarCompute("megr")) megr = 0.31 * f * sqrt(mdudz ** 2 + mdvdz ** 2) / mbvf

                ! output altitude-related vars
                if (conf%isVarOutput("mdudz")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdudz"), ilev, mdudz, 0)
                if (conf%isVarOutput("mdvdz")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdvdz"), ilev, mdvdz, 0)
                if (conf%isVarOutput("mdtdz")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mdtdz"), ilev, mdtdz, 0)
                if (conf%isVarOutput("mbvf")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("mbvf"), ilev, mbvf, 0)
                if (conf%isVarOutput("megr")) call streamWriteVarSlice(fo%streamID, conf%getVarOutID("megr"), ilev, megr, 0)
            endif

            ! store this level's data to be used in the next level
            if (conf%isVarCompute("mdivtflxv")) pmtflxw = mtflxw
            if (conf%isVarCompute("mdivvflxv")) pmvflxw = mvflxw
            if (conf%isVarCompute("mtranstw")) then
                if (conf%tflux) then
                    pmt = mt
                else
                    pmt = mtheta
                endif
            endif
            if (conf%isVarCompute("mdtdz")) pmtheta = mtheta
            if (conf%isVarCompute("mtransvw")) pmita = mita
            if (conf%isVarCompute("mdudz")) pmu = mu
            if (conf%isVarCompute("mdvdz")) pmv = mv
            if (conf%isVarCompute("mdudz") .or. conf%isVarCompute("mdvdz") .or. conf%isVarCompute("mdtdz")) pmgph = mgph
        enddo !ilev
    enddo !ichk

    ! finialize
    call cdiFinish(fi, fo)
    deallocate(lon, lat, p, pk, pk1, pki, lnp, x, y, dx, dy, f)

    if (allocated(mt)) deallocate(mt)
    if (allocated(mu)) deallocate(mu)
    if (allocated(mv)) deallocate(mv)
    if (allocated(mw)) deallocate(mw)
    if (allocated(mtheta)) deallocate(mtheta)
    if (allocated(mita)) deallocate(mita)
    if (allocated(mtflxu)) deallocate(mtflxu)
    if (allocated(mtflxv)) deallocate(mtflxv)
    if (allocated(mdivxtflxu)) deallocate(mdivxtflxu)
    if (allocated(mdivytflxv)) deallocate(mdivytflxv)
    if (allocated(mdivtflxh)) deallocate(mdivtflxh)
    if (allocated(mtflxw)) deallocate(mtflxw)
    if (allocated(mdivtflxv)) deallocate(mdivtflxv)
    if (allocated(mtt)) deallocate(mtt)
    if (allocated(mtranstu)) deallocate(mtranstu)
    if (allocated(mtranstv)) deallocate(mtranstv)
    if (allocated(mtranstw)) deallocate(mtranstw)
    if (allocated(mdiabh)) deallocate(mdiabh)
    if (allocated(mvflxu)) deallocate(mvflxu)
    if (allocated(mvflxv)) deallocate(mvflxv)
    if (allocated(mdivxvflxu)) deallocate(mdivxvflxu)
    if (allocated(mdivyvflxv)) deallocate(mdivyvflxv)
    if (allocated(mdivvflxh)) deallocate(mdivvflxh)
    if (allocated(mvflxw)) deallocate(mvflxw)
    if (allocated(mdivvflxv)) deallocate(mdivvflxv)
    if (allocated(mtransvu)) deallocate(mtransvu)
    if (allocated(mtransvv)) deallocate(mtransvv)
    if (allocated(mtransvw)) deallocate(mtransvw)
    if (allocated(mdudz)) deallocate(mdudz)
    if (allocated(mdvdz)) deallocate(mdvdz)
    if (allocated(mdtdz)) deallocate(mdtdz)
    if (allocated(mbvf)) deallocate(mbvf)
    if (allocated(megr)) deallocate(megr)
    if (allocated(mgph)) deallocate(mgph)
    if (allocated(mgpht)) deallocate(mgpht)
    if (allocated(stdgph)) deallocate(stdgph)
    if (allocated(mmflxuu)) deallocate(mmflxuu)
    if (allocated(mmflxuv)) deallocate(mmflxuv)
    if (allocated(mmflxvv)) deallocate(mmflxvv)
    if (allocated(mmflxuw)) deallocate(mmflxuw)
    if (allocated(mmflxvw)) deallocate(mmflxvw)
    if (allocated(mmflxww)) deallocate(mmflxww)
    if (allocated(mdivxmflxuu)) deallocate(mdivxmflxuu)
    if (allocated(mdivxmflxuv)) deallocate(mdivxmflxuv)
    if (allocated(mdivymflxuv)) deallocate(mdivymflxuv)
    if (allocated(mdivymflxvv)) deallocate(mdivymflxvv)
    if (allocated(mdivmflxh)) deallocate(mdivmflxh)
    if (allocated(meke)) deallocate(meke)

    if (allocated(pmtflxw)) deallocate(pmtflxw)
    if (allocated(pmvflxw)) deallocate(pmvflxw)
    if (allocated(pmt)) deallocate(pmt)
    if (allocated(pmtheta)) deallocate(pmtheta)
    if (allocated(pmita)) deallocate(pmita)
    if (allocated(pmu)) deallocate(pmu)
    if (allocated(pmv)) deallocate(pmv)
    if (allocated(pmgph)) deallocate(pmgph)
    if (allocated(dgph)) deallocate(dgph)

    ! message
    call cpu_time(timeend)
    if (.not. conf%silent) call message(conf%nout, conf%chunks, chunklen, timeend - timebegin)
contains
    subroutine getLevData(conf, fi, pki, tsfirst, tslast, ilev, mt, mtheta, mu, mv, mw, mtt, &
            muptp, mvptp, mwptp, mita, mupip, mvpip, mwpip, mgph, mgpht, stdgph, &
            mmflxuu, mmflxuv, mmflxvv, mmflxuw, mmflxvw, mmflxww, meke)
        type(Config), intent(in) :: conf
        type(FileInfoI), intent(in) :: fi
        integer, intent(in) :: tsfirst, tslast, ilev
        real(8), dimension(:), intent(in) :: pki
        real(8), dimension(:, :), intent(out) :: mt, mtheta, mu, mv, mw, mtt, mita, mgph, mgpht
        real(8), dimension(:, :), intent(out) :: muptp, mvptp, mwptp, mupip, mvpip, mwpip, stdgph
        real(8), dimension(:, :), intent(out) :: mmflxuu, mmflxuv, mmflxvv, mmflxuw, mmflxvw, mmflxww
        real(8), dimension(:, :), intent(out) :: meke

        ! temporary vars
        real(8), dimension(:, :), allocatable :: t, u, v, w, pt, ita, gph, pgph

        integer :: pyear, pmonth, pday, phour, pminute, psecond
        integer :: vyear, vmonth, vday, vhour, vminute, vsecond
        integer :: itime, nmiss, stat
        real(8) :: pjul, vjul, interv

        stat = streamInqTimestep(fi%streamID, tsfirst)
        call decodeDate(taxisInqVdate(fi%taxisID), pyear, pmonth, pday)
        call decodeTime(taxisInqVtime(fi%taxisID), phour, pminute, psecond)
        pjul = julday(pyear, pmonth, pday, phour, pminute, psecond * 1.0_8)

        ! if to be computed, they must have already been allocated
        if (conf%isVarCompute("mt")) call streamReadVarSlice(fi%streamID, fi%tID, ilev, mt, nmiss)
        if (conf%isVarCompute("mu")) call streamReadVarSlice(fi%streamID, fi%uID, ilev, mu, nmiss)
        if (conf%isVarCompute("mv")) call streamReadVarSlice(fi%streamID, fi%vID, ilev, mv, nmiss)
        if (conf%isVarCompute("mw")) call streamReadVarSlice(fi%streamID, fi%wID, ilev, mw, nmiss)
        if (conf%isVarCompute("mita")) call streamReadVarSlice(fi%streamID, fi%itaID, ilev, mita, nmiss)
        if (conf%isVarCompute("mgph")) call streamReadVarSlice(fi%streamID, fi%gphID, ilev, mgph, nmiss)
        if (conf%isVarCompute("mtt")) mtt = 0
        if (conf%isVarCompute("mgpht")) mgpht = 0

        if (.not. conf%aveMode) then
            if (conf%isVarCompute("mt")) allocate(T(fi%nx, fi%ny))
            if (conf%isVarCompute("mu")) allocate(u(fi%nx, fi%ny))
            if (conf%isVarCompute("mv")) allocate(v(fi%nx, fi%ny))
            if (conf%isVarCompute("mw")) allocate(w(fi%nx, fi%ny))
            if (conf%isVarCompute("mita")) allocate(ita(fi%nx, fi%ny))
            if (conf%isVarCompute("mgph")) allocate(gph(fi%nx, fi%ny))
            if (conf%isVarCompute("mtt")) then
                allocate(pT(fi%nx, fi%ny))
                pT = 0
            endif
            if (conf%isVarCompute("mgpht")) then
                allocate(pgph(fi%nx, fi%ny))
                pgph = 0
            endif

            ! loop 1: compute time mean
            do itime = tsfirst + 1, tslast
                stat = streamInqTimestep(fi%streamID, itime)
                call decodeDate(taxisInqVdate(fi%taxisID), vyear, vmonth, vday)
                call decodeTime(taxisInqVtime(fi%taxisID), vhour, vminute, vsecond)
                vjul = julday(vyear, vmonth, vday, vhour, vminute, vsecond * 1.0_8)
                interv =  (vjul - pjul) * spd

                if (conf%isVarCompute("mt")) then
                    call streamReadVarSlice(fi%streamID, fi%tID, ilev, T, nmiss)
                    mT = mT + T
                endif
                if (conf%isVarCompute("mu")) then
                    call streamReadVarSlice(fi%streamID, fi%uID, ilev, u, nmiss)
                    mu = mu + u
                endif
                if (conf%isVarCompute("mv")) then
                    call streamReadVarSlice(fi%streamID, fi%vID, ilev, v, nmiss)
                    mv = mv + v
                endif
                if (conf%isVarCompute("mw")) then
                    call streamReadVarSlice(fi%streamID, fi%wID, ilev, w, nmiss)
                    mw = mw + w
                endif
                if (conf%isVarCompute("mita")) then
                    call streamReadVarSlice(fi%streamID, fi%itaID, ilev, ita, nmiss)
                    mita = mita + ita
                endif
                if (conf%isVarCompute("mgph")) then
                    call streamReadVarSlice(fi%streamID, fi%gphID, ilev, gph, nmiss)
                    mgph = mgph + gph
                endif

                ! dT/dt: time forward difference
                if (conf%isVarCompute("mtt")) then
                    mtt = mtt + (t - pt) / interv
                    pt = t
                endif
                if (conf%isVarCompute("mgpht")) then
                    mgpht = mgpht + (gph - pgph) / interv
                    pgph = gph
                endif

                pjul = vjul
            enddo
            if (conf%isVarCompute("mt")) mt = mt / fi%nt
            if (conf%isVarCompute("mu")) mu = mu / fi%nt
            if (conf%isVarCompute("mv")) mv = mv / fi%nt
            if (conf%isVarCompute("mw")) mw = mw / fi%nt
            if (conf%isVarCompute("mita")) mita = mita / fi%nt
            if (conf%isVarCompute("mgph")) mgph = mgph / fi%nt
            if (conf%isVarCompute("mtt")) mtt = mtt / (fi%nt - 1)
            if (conf%isVarCompute("mgpht")) mgpht = mgpht / (fi%nt - 1)
            if (conf%isVarCompute('mtheta')) mtheta = mT * pki(ilev + 1)

            ! loop 2: compute perturbence
            if (conf%isVarCompute("mtflxu")) muptp = 0
            if (conf%isVarCompute("mtflxv")) mvptp = 0
            if (conf%isVarCompute("mtflxw")) mwptp = 0
            if (conf%isVarCompute("mvflxu")) mupip = 0
            if (conf%isVarCompute("mvflxv")) mvpip = 0
            if (conf%isVarCompute("mvflxw")) mwpip = 0
            if (conf%isVarCompute("stdgph")) stdgph = 0
            if (conf%isVarCompute("mmflxuu")) mmflxuu = 0
            if (conf%isVarCompute("mmflxuv")) mmflxuv = 0
            if (conf%isVarCompute("mmflxvv")) mmflxvv = 0
            if (conf%isVarCompute("mmflxuw")) mmflxuw = 0
            if (conf%isVarCompute("mmflxvw")) mmflxvw = 0
            if (conf%isVarCompute("mmflxww")) mmflxww = 0
            if (conf%isVarCompute("meke")) meke = 0
            do itime = tsfirst, tslast
                stat = streamInqTimestep(fi%streamID, itime)
                if (conf%isVarCompute("mt")) call streamReadVarSlice(fi%streamID, fi%tID, ilev, t, nmiss)
                if (conf%isVarCompute("mita")) call streamReadVarSlice(fi%streamID, fi%itaID, ilev, ita, nmiss)
                if (conf%isVarCompute("mgph")) call streamReadVarSlice(fi%streamID, fi%gphID, ilev, gph, nmiss)
                if (conf%isVarCompute("mtflxu") .or. conf%isVarCompute("mvflxu")) &
                    call streamReadVarSlice(fi%streamID, fi%uID, ilev, u, nmiss)
                if (conf%isVarCompute("mtflxv") .or. conf%isVarCompute("mvflxv")) &
                    call streamReadVarSlice(fi%streamID, fi%vID, ilev, v, nmiss)
                if (conf%isVarCompute("mtflxw") .or. conf%isVarCompute("mvflxw")) &
                    call streamReadVarSlice(fi%streamID, fi%wID, ilev, w, nmiss)

                if (conf%isVarCompute("mt")) then
                    if (conf%tflux) then ! use temperature
                        T = T - mt
                    else  !use potential temperature
                        T = pki(ilev + 1) * T - mtheta
                    endif
                endif
                if (conf%isVarCompute("mita")) ita = ita - mita
                if (conf%isVarCompute("mu")) u = u - mu
                if (conf%isVarCompute("mv")) v = v - mv
                if (conf%isVarCompute("mw")) w = w - mw
                if (conf%isVarCompute("mtflxu")) muptp = muptp + u * T
                if (conf%isVarCompute("mtflxv")) mvptp = mvptp + v * T
                if (conf%isVarCompute("mtflxw")) mwptp = mwptp + w * T
                if (conf%isVarCompute("mvflxu")) mupip = mupip + u * ita
                if (conf%isVarCompute("mvflxv")) mvpip = mvpip + v * ita
                if (conf%isVarCompute("mvflxw")) mwpip = mwpip + w * ita
                if (conf%isVarCompute("stdgph")) stdgph = stdgph + (gph - mgph) ** 2

                if (conf%isVarCompute("mmflxuu")) mmflxuu = mmflxuu + u * u
                if (conf%isVarCompute("mmflxuv")) mmflxuv = mmflxuv + u * v
                if (conf%isVarCompute("mmflxvv")) mmflxvv = mmflxvv + v * v
                if (conf%isVarCompute("mmflxuw")) mmflxuw = mmflxuw + u * w
                if (conf%isVarCompute("mmflxvw")) mmflxvw = mmflxvw + v * w
                if (conf%isVarCompute("mmflxww")) mmflxww = mmflxww + w * w
                if (conf%isVarCompute("meke")) meke = meke + u ** 2 + v ** 2
            enddo ! itime
            if (conf%isVarCompute("mtflxu")) muptp = muptp / fi%nt
            if (conf%isVarCompute("mtflxv")) mvptp = mvptp / fi%nt
            if (conf%isVarCompute("mtflxw")) mwptp = mwptp / fi%nt
            if (conf%isVarCompute("mvflxu")) mupip = mupip / fi%nt
            if (conf%isVarCompute("mvflxv")) mvpip = mvpip / fi%nt
            if (conf%isVarCompute("mvflxw")) mwpip = mwpip / fi%nt
            if (conf%isVarCompute("stdgph")) stdgph = sqrt(stdgph / (fi%nt - 1))
            if (conf%isVarCompute("mmflxuu")) mmflxuu = mmflxuu / fi%nt
            if (conf%isVarCompute("mmflxuv")) mmflxuv = mmflxuv / fi%nt
            if (conf%isVarCompute("mmflxvv")) mmflxvv = mmflxvv / fi%nt
            if (conf%isVarCompute("mmflxuw")) mmflxuw = mmflxuw / fi%nt
            if (conf%isVarCompute("mmflxvw")) mmflxvw = mmflxvw / fi%nt
            if (conf%isVarCompute("mmflxww")) mmflxww = mmflxww / fi%nt
            if (conf%isVarCompute("meke")) meke = meke / fi%nt / 2
        endif

        if (allocated(t)) deallocate(t)
        if (allocated(u)) deallocate(u)
        if (allocated(v)) deallocate(v)
        if (allocated(w)) deallocate(w)
        if (allocated(ita)) deallocate(ita)
        if (allocated(gph)) deallocate(gph)
        if (allocated(pt)) deallocate(pt)
        if (allocated(pgph)) deallocate(pgph)
    end subroutine getLevData
    
    function ddx(var, dx) result(dvardx)
        real(8), dimension(:, :), intent(in) :: var
        real(8), dimension(:, :), intent(in) :: dx
        real(8), dimension(size(var, 1), size(var, 2)) :: dvardx

        integer :: nx, ny

        nx = size(var, 1)
        ny = size(var, 2)

        dvardx(2 : nx - 1, :) = (var(3 : nx, :) - var(1 : nx - 2, :)) / dx(2 : nx - 1, :)
        dvardx(1, :) = (var(2, :) - var(nx, :)) / dx(1, :)
        dvardx(nx, :) = (var(1, :) - var(nx - 1, :)) / dx(nx, :)
    end function ddx

    function ddy(var, dy) result(dvardy)
        real(8), dimension(:, :), intent(in) :: var
        real(8), dimension(:, :), intent(in) :: dy
        real(8), dimension(size(var, 1), size(var, 2)) :: dvardy

        integer :: nx, ny

        nx = size(var, 1)
        ny = size(var, 2)

        dvardy(:, 2 : ny - 1) = (var(:, 3 : ny) - var(:, 1 : ny - 2)) / dy(:, 2 : ny - 1)
        dvardy(:, 1) = (var(:, 2) - var(:, 1)) / dy(:, 1)
        dvardy(:, ny) = (var(:, ny) - var(:, ny - 1)) / dy(:, ny)
    end function ddy

    function ddp(dvar, p, dlnp) result(dvardp)
        real(8), dimension(:, :), intent(in) :: dvar
        real(8), intent(in) :: p, dlnp
        real(8), dimension(size(dvar, 1), size(dvar, 2)) :: dvardp

        dvardp = dvar / dlnp / p
    end function ddp

    function ddz(dvar, dlnz) result(dvardz)
        real(8), dimension(:, :), intent(in) :: dvar
        real(8), dimension(:, :), intent(in) :: dlnz
        real(8), dimension(size(dvar, 1), size(dvar, 2)) :: dvardz

        dvardz = dvar / dlnz
    end function ddz

    subroutine printUsage
        write(*,'(3A)')"AFTER2: another afterburner for ECHAM5/6 (version ", VERSION, ")"
        write(*,'(A)')"Compute additional dignostic variables"
        write(*,'(A)')"Input file must be the output of the original afterburner containing"
        write(*,'(A)')"(any subset of) the following (3D, time-varying) variables:"
        write(*,'(A)')"    t: temperature (K)"
        write(*,'(A)')"    u: u-velocity (m/s)"
        write(*,'(A)')"    v: v-velocity (m/s)"
        write(*,'(A)')"    w: vertical velocity (Pa/s)"
        write(*,'(A)')"    svo: relative vorticity (1/s)"
        write(*,'(A)')"    geopoth: geopotential height (m)"
        write(*,'(A)')"If some of the input variables are missing, the output variables"
        write(*,'(A)')"which depend on them are disabled automatically"
        write(*,'(A)')"For currently supported output variables, see the documentation"
        write(*,'(A)')"Caution: this program depends on a environment variable"
        write(*,'(A)')"    AFTER2_VARDEF_PATH which is normally the install path of the program"
        write(*,'(A)')"    You must define this path properly. For more information see the"
        write(*,'(A)')"    documentation"
        write(*,'(A)')"Usage:"
        write(*,'(A)')"   after2 [options] INPUT OUTPUT"
        write(*,'(A)')"Options:"
        write(*,'(A)')"   -v[variable list]: output variable list sepereated by commas"
        write(*,'(A)')"   -x[variable list]: excluding variables"
        write(*,'(A)')"   -ave: enable average mode"
        write(*,'(A)')"   -tflux: use temperature instead of potential temperature in heat fluxes"
        write(*,'(A)')"   -c[1..]: number of chunks to divide the input time range into"
        write(*,'(A)')"   -nc: output file type netCDF"
        write(*,'(A)')"   -nc2: output file type netCDF version 2 (64-bit)"
        write(*,'(A)')"   -nc4: output file type netCDF-4 (HDF5)"
        write(*,'(A)')"   -nc4c: output file type netCDF-4 classic"
        write(*,'(A)')"   -z: compress in ZIP format"
        write(*,'(A)')"   -gz: compress in GZIP format"
        write(*,'(A)')"   -bz2: compress in BZIP2 format"
        write(*,'(A)')"   -jpg: compress in JPEG format"
        write(*,'(A)')"   -l[0..9]: compress level (0 for none, 9 for most, default 6)"
        write(*,'(A)')"   -a: absolute time axis"
        write(*,'(A)')"   -r[YYYYMMDD[hhmmss]]: relative time axis with reference date & time"
        write(*,'(A)')"   -u[unit]: time axis unit. supported units are:"
        write(*,'(A)')"       d/day/days"
        write(*,'(A)')"       h/hour/hours"
        write(*,'(A)')"       m/min/minute/minutes"
        write(*,'(A)')"       s/sec/second/seconds"
        write(*,'(A)')"   -C[calendar]: calendar. supported calendars are:"
        write(*,'(A)')"       std/standard"
        write(*,'(A)')"       prol/proleptic/proleptic_gregorian"
        write(*,'(A)')"       360/360_days"
        write(*,'(A)')"       365/365_days"
        write(*,'(A)')"       366/366_days"
        write(*,'(A)')"   -s: silent mode"
        write(*,'(A)')"   -h: print this message and exit (must be the only argument)"
        write(*,'(A)')"   -V: print version and exit (must be the only argument)"
        write(*,'(A)')"Author:"
        write(*,'(A)')"   Guidi Zhou"
        write(*,'(3X, A)')VERSION_DATE
        write(*,'(A)')"   Helmholtz-Center for Ocean Research"
        write(*,'(A)')"   Kiel, Germany"
        write(*,'(A)')"   zhouguidi@gmail.com"
        stop
    end subroutine printUsage

    subroutine printVersion
        write(*,'(5A)')"after2 version ", VERSION, " (", VERSION_DATE, ")"
        stop
    end subroutine printVersion

    subroutine perseArgs(conf)
        type(Config), intent(inout) :: conf

        integer :: iarg, nargin, ind, i, j, stat
        character(100) :: arg
        character(:), allocatable :: tok, rem
        logical :: lst

        nargin = command_argument_count()

        select case(nargin)
        case (0) ! no conf
            call printUsage
            stop
        case (1) ! one arg
            call get_command_argument(1, arg)
            arg = adjustl(arg)
            select case(trim(arg))
            case ('-V')
                call printVersion
                stop
            case ('-h')
                call printUsage
                stop
            case default
                call errorHandler("invalid input argument.")
            end select
        case default
            ! the last two conf are file names
            call get_command_argument(nargin - 1, conf%iFileName)
            conf%iFileName = adjustl(conf%iFileName)
            call get_command_argument(nargin, conf%oFileName)
            conf%oFileName = adjustl(conf%oFileName)
            nargin = nargin - 2
            
            ! the options
            do iarg = 1, nargin
                call get_command_argument(iarg, arg)
                arg = adjustl(arg)

                if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char
                if (arg(1 : 1) .ne. "-") call errorHandler("invalid input argument") ! not starting with -

                arg = arg(2 : len_trim(arg)) ! now without -
                select case(arg(1 : 1)) ! swtich options
                case ("v")
                    if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (v)
                    arg = arg(2 : len_trim(arg)) ! now the actual variable list

                    ! given -v option, first disable every variable from output
                    call conf%disableAllVars

                    ! loop variables in the list
                    do
                        call strtok(arg, tok, rem, last = lst, stat = stat, deli = ",", allowblk = .false.)
                        if (stat /= 0) call errorHandler("invalid variable list")

                        ! now tok contains the current variable found in the list
                        call conf%enableVar(tok)

                        arg = rem
                        if (lst) exit
                    enddo
                case ("x")
                    if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (v)
                    arg = arg(2 : len_trim(arg)) ! now the actual variable list

                    ! given -x option, first enable every variable from output
                    call conf%enableAllVars

                    ! loop variables in the list
                    do
                        call strtok(arg, tok, rem, last = lst, stat = stat, deli = ",", allowblk = .false.)
                        if (stat /= 0) call errorHandler("invalid variable list")

                        ! now tok contains the current variable found in the list
                        call conf%disableVar(tok)

                        arg = rem
                        if (lst) exit
                    enddo
                case ("t")
                    if (trim(arg) == 'tflux') then
                        conf%tflux = .true.
                    else
                        call errorHandler("unknown argument")
                    endif

                case ("c")
                    if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (c)
                    tok = arg(2 : len_trim(arg))
                    read(tok, *, iostat = stat)conf%chunks
                    if (stat .ne. 0) call errorHandler("invalid chuank number") ! can't read chunk number
                case ("n")
                    select case(trim(arg))
                    case ("nc")
                        conf%filetype = FILETYPE_NC
                    case ("nc2")
                        conf%filetype = FILETYPE_NC2
                    case ("nc4")
                        conf%filetype = FILETYPE_NC4
                    case ("nc4c")
                        conf%filetype = FILETYPE_NC4C
                    case default
                        call errorHandler("invalid output file format")
                    end select
                case ("z")
                    select case(trim(arg))
                    case ("z", "zip")
                        conf%comptype = COMPRESS_ZIP
                    case default
                        call errorHandler("invalid compress type")
                    end select
                case ("g")
                    select case(trim(arg))
                    case ("g", "gz", "gzip")
                        conf%comptype = COMPRESS_GZIP
                    case default
                        call errorHandler("invalid compress type")
                    end select
                case ("b")
                    select case(trim(arg))
                    case ("b", "bz", "bz2")
                        conf%comptype = COMPRESS_BZIP2
                    case default
                        call errorHandler("invalid compress type")
                    end select
                case ("j")
                    select case(trim(arg))
                    case ("j", "jpg", "jpeg")
                        conf%comptype = COMPRESS_JPEG
                    case default
                        call errorHandler("invalid compress type")
                    end select
                case ("l")
                    if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (l)
                    tok = arg(2 : len_trim(arg))
                    read(tok, *, iostat = stat)conf%complev
                    if (stat .ne. 0) call errorHandler("invalid compress level") ! can't read compress level
                    if (conf%complev .lt. 0 .or. conf%complev .gt. 9) &
                        call errorHandler("invalid compress level") ! compress level wrong
                case ("a")
                    if (len_trim(arg) .eq. 1) then
                        conf%tType = TAXIS_ABSOLUTE
                    elseif (trim(arg) .eq. "ave") then
                        conf%aveMode = .true.
                    else
                        call errorHandler("unknown argument")
                    endif
                case ("r")
                    conf%tType = TAXIS_RELATIVE
                    if (len_trim(arg) .eq. 1) then ! default reference date & time
                        conf%refDate = 19700101
                        conf%refTime = 0
                    else
                        arg = arg(2 : len_trim(arg))
                        if (len_trim(arg) .eq. 8) then ! only ref date is given
                            read(arg, *, iostat = stat)conf%refDate
                            if (stat .ne. 0) call errorHandler("invalid reference date")
                        elseif (len_trim(arg) .eq. 14) then ! both ref date and time are givem
                            tok = arg(1 : 8)
                            read(tok, *, iostat = stat)conf%refDate
                            if (stat .ne. 0) call errorHandler("invalid reference date")
                            tok = arg(9 : 14)
                            read(tok, *, iostat = stat)conf%refTime
                            if (stat .ne. 0) call errorHandler("invalid reference time")
                        else ! incorrect ref date & time
                            call errorHandler("invalid reference date & time")
                        endif
                    endif
                case ("u")
                    if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (u)
                    arg = arg(2 : len_trim(arg))
                    select case(trim(arg))
                    case ("d", "day", "days")
                        conf%tUnit = TUNIT_DAY
                    case ("h", "hour", "hours")
                        conf%tUnit = TUNIT_HOUR
                    case ("m", "min", "minute", "minutes")
                        conf%tUnit = TUNIT_MINUTE
                    case ("s", "sec", "second", "seconds")
                        conf%tUnit = TUNIT_SECOND
                    case default
                        call errorHandler("invalid time unit")
                    end select
                case ("C")
                    if (len_trim(arg) .eq. 1) call errorHandler("invalid input argument") ! only one char (C)
                    arg = arg(2 : len_trim(arg))
                    select case(trim(arg))
                    case ("std", "standard")
                        conf%calendar = CALENDAR_STANDARD
                    case ("prol", "proleptic", "proleptic_gregorian")
                        conf%calendar = CALENDAR_PROLEPTIC
                    case ("360", "360_days")
                        conf%calendar = CALENDAR_360DAYS
                    case ("365", "365_days")
                        conf%calendar = CALENDAR_365DAYS
                    case ("366", "366_days")
                        conf%calendar = CALENDAR_366DAYS
                    case default
                        call errorHandler("invalid calendar")
                    end select
                case ("s")
                    if (len_trim(arg) .eq. 1) then
                        conf%silent = .true.
                    else
                        call errorHandler("unknown argument")
                    endif
                case ("h")
                    if (len_trim(arg) .eq. 1) then
                        call errorHandler("-h must be used as the only input argument")
                    else
                        call errorHandler("unknown argument")
                    end if
                case ("V")
                    if (len_trim(arg) .eq. 1) then
                        call errorHandler("-V must be used as the only input argument")
                    else
                        call errorHandler("unknown argument")
                    end if
                case default
                    call errorHandler("unknown argument")
                end select
            enddo
        end select
    end subroutine perseArgs

    subroutine message(nout, chunks, chunklen, time)
        integer, intent(in) :: nout, chunks, chunklen
        real(8), intent(in) :: time

        character(1000) :: msg, temp

        msg = "after2: processed"
        write(temp, '(I20)')nout
        msg = trim(msg) // " " // trim(adjustl(temp)) // " variables over"
        write(temp, '(I20)')chunks
        msg = trim(msg) // " " // trim(adjustl(temp)) // " chunks of"
        write(temp, '(I20)')chunklen
        msg = trim(msg) // " " // trim(adjustl(temp)) // " timesteps ("
        write(temp, '(F20.3)')time
        msg = trim(msg) // trim(adjustl(temp)) // "s)"
        write(*, '(A)')trim(msg)
    end subroutine message
end program after2
