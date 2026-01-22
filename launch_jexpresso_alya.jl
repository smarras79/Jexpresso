#!/usr/bin/env julia
#
# Launcher for coupling Jexpresso with another code
#
# This script launches Jexpresso coupled with either:
#   - Another Jexpresso instance (for testing)
#   - Alya (for production)
#
# Usage:
#   mpirun -np 8 julia --project=. launch_jexpresso_alya.jl
#
# Configuration:
#   - Set SECOND_CODE to :Jexpresso or :Alya below
#   - Edit rank distribution
#   - Edit problem configurations
#

using MPI

MPI.Init()
world_rank = MPI.Comm_rank(MPI.COMM_WORLD)
world_size = MPI.Comm_size(MPI.COMM_WORLD)

#-----------------------------------
# CONFIGURATION - EDIT THESE VALUES
#-----------------------------------

# Choose second code: :Jexpresso (for testing) or :Alya (for production)
SECOND_CODE = :Jexpresso  # Options: :Jexpresso or :Alya

# How many ranks for each code
JEXPRESSO_RANKS = 2
SECOND_CODE_RANKS = 2

# Jexpresso problem configuration
JEXPRESSO_EQS = "CompEuler"      # Equation type
JEXPRESSO_CASE = "theta"         # Test case name

# Second Jexpresso instance configuration (if SECOND_CODE = :Jexpresso)
JEXPRESSO2_EQS = "CompEuler"     # Can be same or different
JEXPRESSO2_CASE = "theta"        # Can be same or different

# Alya configuration (if SECOND_CODE = :Alya)
ALYA_INPUT_FILE = "alya.dat"     # Alya input file
ALYA_EXECUTABLE = "/path/to/alya/bin/alya.x"  # Path to Alya executable

#-----------------------------------
# VALIDATION
#-----------------------------------

if world_size != JEXPRESSO_RANKS + SECOND_CODE_RANKS
    if world_rank == 0
        println("ERROR: Total MPI ranks mismatch!")
        println("  Requested: Jexpresso=$JEXPRESSO_RANKS + $SECOND_CODE=$SECOND_CODE_RANKS = $(JEXPRESSO_RANKS + SECOND_CODE_RANKS)")
        println("  Actual:    $world_size")
        println()
        println("Usage: mpirun -np $(JEXPRESSO_RANKS + SECOND_CODE_RANKS) julia --project=. launch_jexpresso_alya.jl")
    end
    MPI.Finalize()
    exit(1)
end

# Validate SECOND_CODE choice
if SECOND_CODE != :Jexpresso && SECOND_CODE != :Alya
    if world_rank == 0
        println("ERROR: Invalid SECOND_CODE setting!")
        println("  SECOND_CODE = $SECOND_CODE")
        println("  Valid options: :Jexpresso or :Alya")
    end
    MPI.Finalize()
    exit(1)
end

#-----------------------------------
# RANK ASSIGNMENT
#-----------------------------------

# Ranks 0 to (JEXPRESSO_RANKS-1): Jexpresso (code 1)
# Ranks JEXPRESSO_RANKS to (world_size-1): Second code (code 2)

if world_rank < JEXPRESSO_RANKS
    #-----------------------------------
    # JEXPRESSO CODE (CODE 1)
    #-----------------------------------

    if world_rank == 0
        println("="^70)
        println("COUPLED RUN: Jexpresso + $SECOND_CODE")
        println("="^70)
        println("Total MPI ranks:     $world_size")
        println("Jexpresso ranks:     0-$(JEXPRESSO_RANKS-1) ($JEXPRESSO_RANKS total)")
        println("$SECOND_CODE ranks:  $JEXPRESSO_RANKS-$(world_size-1) ($SECOND_CODE_RANKS total)")
        println("Jexpresso problem:   $JEXPRESSO_EQS/$JEXPRESSO_CASE")
        if SECOND_CODE == :Jexpresso
            println("Jexpresso-2 problem: $JEXPRESSO2_EQS/$JEXPRESSO2_CASE")
        end
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
    #-----------------------------------
    # SECOND CODE (CODE 2)
    #-----------------------------------

    if SECOND_CODE == :Jexpresso
        #-----------------------------------
        # SECOND JEXPRESSO INSTANCE
        #-----------------------------------

        if world_rank == JEXPRESSO_RANKS
            println("\n[Jexpresso-2] Starting second Jexpresso instance on ranks $JEXPRESSO_RANKS-$(world_size-1)")
            flush(stdout)
        end

        # Set command-line arguments for second Jexpresso
        push!(empty!(ARGS),
              "--coupling",
              "--code-id", "2",
              "--n-codes", "2",
              "--code-name", "Jexpresso-2",
              JEXPRESSO2_EQS,
              JEXPRESSO2_CASE)

        # Load and run second Jexpresso instance
        if world_rank == JEXPRESSO_RANKS
            println("[Jexpresso-2] Loading Jexpresso module...")
            flush(stdout)
        end

        using Jexpresso

        # Jexpresso handles MPI.Finalize() internally

    elseif SECOND_CODE == :Alya
        #-----------------------------------
        # ALYA CODE
        #-----------------------------------

        if world_rank == JEXPRESSO_RANKS
            println("\n[Alya] Starting Alya on ranks $JEXPRESSO_RANKS-$(world_size-1)")
            flush(stdout)
        end

        # -------------------------------------------------------------------
        # CRITICAL: Initialize coupling on Alya side too!
        # -------------------------------------------------------------------
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

        # -------------------------------------------------------------------
        # IMPLEMENT ALYA LAUNCH HERE
        # -------------------------------------------------------------------

        # Option 1: If Alya has a Julia interface/wrapper
        # -------------------------------------------------------------------
        # include("/path/to/alya/julia_wrapper.jl")
        # run_alya_coupled(
        #     ctx = alya_ctx,  # Pass coupling context
        #     input_file = ALYA_INPUT_FILE
        # )

        # Option 2: If you need to call Alya's Fortran/C functions via ccall
        # -------------------------------------------------------------------
        # This is more complex and requires proper interfacing
        # You may need to wrap Alya's main routine

        # Option 3: Placeholder for testing (just synchronize)
        # -------------------------------------------------------------------
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

end

# Note: Each code handles its own MPI.Finalize()
# Jexpresso does it internally, Alya should do it too
