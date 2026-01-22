#!/usr/bin/env julia
#
# Launcher for coupling Jexpresso with Alya
#
# This script launches Jexpresso and Alya together in a single MPI job,
# splitting ranks between the two codes.
#
# Usage:
#   mpirun -np 8 julia --project=. launch_jexpresso_alya.jl
#
# Configuration:
#   - Edit the rank distribution below
#   - Edit Jexpresso problem (CompEuler/theta by default)
#   - Implement Alya launch in the else block
#

using MPI

MPI.Init()
world_rank = MPI.Comm_rank(MPI.COMM_WORLD)
world_size = MPI.Comm_size(MPI.COMM_WORLD)

#============================================================================
# CONFIGURATION - EDIT THESE VALUES
#============================================================================

# How many ranks for each code
JEXPRESSO_RANKS = 4
ALYA_RANKS = 4

# Jexpresso problem configuration
JEXPRESSO_EQS = "CompEuler"      # Equation type
JEXPRESSO_CASE = "theta"         # Test case name

# Alya configuration (if needed)
ALYA_INPUT_FILE = "alya.dat"     # Alya input file
ALYA_EXECUTABLE = "/path/to/alya/bin/alya.x"  # Path to Alya executable

#============================================================================
# VALIDATION
#============================================================================

if world_size != JEXPRESSO_RANKS + ALYA_RANKS
    if world_rank == 0
        println("ERROR: Total MPI ranks mismatch!")
        println("  Requested: Jexpresso=$JEXPRESSO_RANKS + Alya=$ALYA_RANKS = $(JEXPRESSO_RANKS + ALYA_RANKS)")
        println("  Actual:    $world_size")
        println()
        println("Usage: mpirun -np $(JEXPRESSO_RANKS + ALYA_RANKS) julia --project=. launch_jexpresso_alya.jl")
    end
    MPI.Finalize()
    exit(1)
end

#============================================================================
# RANK ASSIGNMENT
#============================================================================

# Ranks 0 to (JEXPRESSO_RANKS-1): Jexpresso
# Ranks JEXPRESSO_RANKS to (world_size-1): Alya

if world_rank < JEXPRESSO_RANKS
    #========================================================================
    # JEXPRESSO CODE
    #========================================================================

    if world_rank == 0
        println("="^70)
        println("COUPLED RUN: Jexpresso + Alya")
        println("="^70)
        println("Total MPI ranks:     $world_size")
        println("Jexpresso ranks:     0-$(JEXPRESSO_RANKS-1) ($JEXPRESSO_RANKS total)")
        println("Alya ranks:          $JEXPRESSO_RANKS-$(world_size-1) ($ALYA_RANKS total)")
        println("Jexpresso problem:   $JEXPRESSO_EQS/$JEXPRESSO_CASE")
        println("="^70)
        flush(stdout)
    end

    # Set Jexpresso command-line arguments
    push!(empty!(ARGS),
          "--coupling",
          "--code-id", "1",
          "--n-codes", "2",
          "--code-name", "Jexpresso",
          JEXPRESSO_EQS,
          JEXPRESSO_CASE)

    # Load and run Jexpresso
    if world_rank == 0
        println("\n[Jexpresso] Loading Jexpresso module...")
        flush(stdout)
    end

    using Jexpresso

    # Jexpresso handles MPI.Finalize() internally

else
    #========================================================================
    # ALYA CODE
    #========================================================================

    if world_rank == JEXPRESSO_RANKS
        println("\n[Alya] Starting Alya on ranks $JEXPRESSO_RANKS-$(world_size-1)")
        flush(stdout)
    end

    # -----------------------------------------------------------------------
    # IMPLEMENT ALYA LAUNCH HERE
    # -----------------------------------------------------------------------

    # Option 1: If Alya has a Julia interface/wrapper
    # -----------------------------------------------------------------------
    # include("/path/to/alya/julia_wrapper.jl")
    # run_alya_coupled(
    #     input_file = ALYA_INPUT_FILE,
    #     code_id = 2,
    #     n_codes = 2,
    #     comm = MPI.COMM_WORLD
    # )

    # Option 2: If you need to call Alya's Fortran/C functions via ccall
    # -----------------------------------------------------------------------
    # This is more complex and requires proper interfacing
    # You may need to wrap Alya's main routine

    # Option 3: Placeholder for testing (just synchronize)
    # -----------------------------------------------------------------------
    # This allows testing Jexpresso side without Alya running
    if world_rank == JEXPRESSO_RANKS
        @warn """
        Alya launch not implemented yet!

        To complete this launcher:
        1. If Alya is a Julia package: uncomment Option 1 above
        2. If Alya needs ccall interface: implement Option 2
        3. If using MPMD mode: use the MPMD launch method instead

        Current behavior: Alya ranks will just wait and synchronize.
        """
    end

    # Initialize Alya coupling (pseudo-code - adapt to your Alya interface)
    # alya_coupling_init(code_id=2, n_codes=2, comm=MPI.COMM_WORLD)

    # Example: Just synchronize to keep ranks alive
    for i in 1:10
        if world_rank == JEXPRESSO_RANKS
            println("[Alya] Checkpoint $i/10")
            flush(stdout)
        end
        MPI.Barrier(MPI.COMM_WORLD)
        sleep(1)
    end

    if world_rank == JEXPRESSO_RANKS
        println("[Alya] Finished")
        flush(stdout)
    end

    # Alya should handle its own MPI.Finalize()
    # Or call it here if Alya doesn't
    MPI.Finalize()

end

# Note: Each code handles its own MPI.Finalize()
# Jexpresso does it internally, Alya should do it too
