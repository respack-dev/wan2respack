PROGRAM generate_sample_mk
    use m_rd_dat_wfn
    implicit none

    integer :: ik, ix, ig

    ! Initialize iunit_log in the module
    iunit_log = 150

    ! Open the log file for writing
    open(iunit_log, FILE="LOG.mk")

    ! Read data from various input files
    ! These subroutines are responsible for reading symmetry, sample k-points,
    ! lattice parameters, number of k-points, and G-vectors respectively.
    call rd_dat_symmetry
    call rd_dat_sample_k
    call rd_dat_lattice
    call rd_dat_nkm
    call rd_dat_kg

    ! Close the log file after reading operations
    close(iunit_log)

    ! Open the output file for writing sample k-points data
    open(100, FILE="dat.sample_mk")
    rewind(100)

    ! Write the number of k-points (NTK) to the file
    write(100, '(i10)') NTK

    ! Loop over each k-point and write its coordinates to the file
    do ik = 1, NTK
        ! Write the k-point coordinates (SK0) in a formatted manner
        write(100, '(3f15.10)') (SK0(ix, ik), ix = 1, 3)
    end do

    ! Close the sample k-points data file
    close(100)

    ! Open the output file for writing G-vectors data
    open(101, FILE="dat.kg_respack")
    rewind(101)

    ! Loop over each k-point to write the number of G-vectors and the G-vectors themselves
    do ik = 1, NTK
        ! Write the number of G-vectors (NG0) for the current k-point
        write(101, *) NG0(ik)

        ! Loop over each G-vector and write its components
        do ig = 1, NG0(ik)
            ! Write the G-vector components (KG0) for the current k-point
            write(101, *) KG0(:, ig, ik)
        end do
    end do

    ! Close the G-vectors data file
    close(101)

    ! Terminate the program
    stop
end PROGRAM generate_sample_mk
