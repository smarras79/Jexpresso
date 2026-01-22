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
    # CRITICAL: Initialize coupling on Alya side too!
    # -----------------------------------------------------------------------
    # Even if Alya isn't fully implemented yet, we need to initialize
    # coupling so Alya ranks participate in MPI collective operations

    # Load JexpressoCoupling module
    include("src/mpi/JexpressoCoupling.jl")
    using .JexpressoCoupling

    # Initialize coupling for Alya
    if world_rank == JEXPRESSO_RANKS
        println("[Alya] Initializing coupling...")
        flush(stdout)
    end

    alya_ctx = initialize_coupling(
        MPI.COMM_WORLD,
        2,  # code_id for Alya
        2;  # n_codes
        code_name="Alya"
    )

    if world_rank == JEXPRESSO_RANKS
        println("[Alya] Coupling initialized successfully")
        println("  Local size: $(alya_ctx.local_size)")
        println("  World rank: $(alya_ctx.world_rank)")
        println("  Is root: $(alya_ctx.is_root)")
        flush(stdout)
    end

    # -----------------------------------------------------------------------
    # IMPLEMENT ALYA LAUNCH HERE
    # -----------------------------------------------------------------------

    # Option 1: If Alya has a Julia interface/wrapper
    # -----------------------------------------------------------------------
    # include("/path/to/alya/julia_wrapper.jl")
    # run_alya_coupled(
    #     ctx = alya_ctx,  # Pass coupling context
    #     input_file = ALYA_INPUT_FILE
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
        Alya solver not implemented yet!

        Coupling is initialized and working.
        To complete this launcher, implement Alya's solver launch above.

        Current behavior: Alya ranks will synchronize and wait.
        """
    end

    # Simulate Alya running (replace with actual Alya code)
    for i in 1:10
        if world_rank == JEXPRESSO_RANKS
            println("[Alya] Simulation step $i/10")
            flush(stdout)
        end

        # Synchronize all codes
        synchronize_coupling(alya_ctx)

        # Simulate some work
        sleep(0.5)
    end

    if world_rank == JEXPRESSO_RANKS
        println("[Alya] Simulation finished")
        flush(stdout)
    end

    # Cleanup coupling
    finalize_coupling(alya_ctx)

    # Alya should handle its own MPI.Finalize()
    # Or call it here if Alya doesn't
    MPI.Finalize()

end

# Note: Each code handles its own MPI.Finalize()
# Jexpresso does it internally, Alya should do it too
