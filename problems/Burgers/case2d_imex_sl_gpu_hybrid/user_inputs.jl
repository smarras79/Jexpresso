function user_inputs()

    #--------------------------------------------------------------------------
    # HYBRID GPU case for 2-D Burgers IMEX.
    #
    # Same problem as Burgers/case2d_imex_sl, but the per-stage IMPLICIT SOLVE
    # (I - λL) x = b is offloaded to the GPU via JACC while EVERYTHING ELSE
    # (mesh, sem_setup, the explicit rhs!, vector algebra) stays on the host.
    # This deliberately keeps `:backend => CPU()` so none of Jexpresso's
    # (incompletely ported) GPU pipeline is exercised — only the JACC BiCGSTAB
    # solver touches the device.
    #
    # IMPORTANT — JACC's backend is a COMPILE-TIME preference, so it must be set
    # to CUDA and Julia RESTARTED before JACC will use the GPU (a runtime
    # set_backend does NOT switch the current session). One-time setup:
    #
    #     using CUDA, JACC
    #     JACC.set_backend("cuda")          # writes LocalPreferences.toml
    #     # ... exit and restart Julia ...
    #
    # Then, in the fresh session:
    #
    #     using CUDA, JACC
    #     using Jexpresso
    #     Jexpresso.jacc_status()           # MUST print a CuArray type / GPU = YES
    #     Jexpresso.run_case("Burgers", "case2d_imex_sl_gpu_hybrid")
    #
    # No `backend = :cuda` kwarg is needed here: this case keeps :backend => CPU(),
    # so only JACC needs to be on its CUDA backend. If JACC is on its CPU backend
    # the solve simply runs on the CPU (still correct — and the run prints a
    # warning saying so). Requirements: CUDA.jl in the project, a functional
    # NVIDIA GPU, and the mesh ./meshes/gmsh_grids/hexa_TFI_10x10_burgers2d.msh.
    # For AMD GPUs use JACC.set_backend("amdgpu") + AMDGPU.
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # IMEX Runge-Kutta ARS(2,3,2) Butcher tableaux
    #--------------------------------------------------------------------------
    a_32 = 1.0 / 6.0 * (3.0 + 2.0 * sqrt(2.0))

    A_RK = zeros(TFloat, 3, 3)
    A_RK[2, 1] = TFloat(2.0 - sqrt(2.0))
    A_RK[3, 1] = TFloat(1.0 - a_32)
    A_RK[3, 2] = TFloat(a_32)

    b_RK = TFloat[ 1.0 / (2.0 * sqrt(2.0)),
                   1.0 / (2.0 * sqrt(2.0)),
                   1.0 - 1.0 / sqrt(2.0) ]

    c_RK = TFloat[ 0.0, 2.0 - sqrt(2.0), 1.0 ]

    A_RK_tilde = zeros(TFloat, 3, 3)
    A_RK_tilde[2, 1] = TFloat(1.0 - 1.0 / sqrt(2.0))
    A_RK_tilde[2, 2] = TFloat(1.0 - 1.0 / sqrt(2.0))
    A_RK_tilde[3, 1] = TFloat(1.0 / (2.0 * sqrt(2.0)))
    A_RK_tilde[3, 2] = TFloat(1.0 / (2.0 * sqrt(2.0)))
    A_RK_tilde[3, 3] = TFloat(1.0 - 1.0 / sqrt(2.0))

    b_RK_tilde = TFloat[ 1.0 / (2.0 * sqrt(2.0)),
                         1.0 / (2.0 * sqrt(2.0)),
                         1.0 - 1.0 / sqrt(2.0) ]

    c_RK_tilde = TFloat[ 0.0, 2.0 - sqrt(2.0), 1.0 ]

    #--------------------------------------------------------------------------
    # Explicit source S(u): the full explicit RHS via the standard rhs!.
    #--------------------------------------------------------------------------
    function S_fun!(s_j, u, time, params, sem)
        rhs!(s_j, u, params, time)
        s_j .= params.RHS
    end

    #--------------------------------------------------------------------------
    # Implicit operator  L = -μ M⁻¹ K  (host sparse; uploaded to the device CSR
    # once by the offload runner).
    #--------------------------------------------------------------------------
    L_cache = Ref{Any}(nothing)

    function _imex_L(params)
        if L_cache[] === nothing
            μ        = params.inputs[:μ]
            Minv     = Array(params.Minv)       # host vector
            L_global = params.Lap_sparse        # pre-assembled global stiffness K
            L_cache[] = -μ[1] * (Minv .* L_global)
        end
        return L_cache[]
    end

    function L_fun!(l_j, u, time, params)
        mul!(l_j, _imex_L(params), u)
    end

    function build_L(u, time, params)
        return _imex_L(params)
    end

    function lsolve(L_curr, b)
        return L_curr \ b
    end

    inputs = Dict(
        :ode_solver           => IMEX(),
        #---------------------------------------------------------------------------
        # Hybrid GPU offload: host pipeline + GPU implicit solve (see header).
        #   * :backend stays CPU (unset) so mesh/sem_setup/rhs! run on the host.
        #   * :limex_jacc + :limex_jacc_offload route the per-stage solve through
        #     the JACC BiCGSTAB on the device.
        #---------------------------------------------------------------------------
        :limex_jacc           => true,
        :limex_jacc_offload   => true,
        :tinit                => 0.0,
        :tend                 => 0.5,
        :Δt                   => 1.0e-3,
        :diagnostics_at_times => (0.1:0.1:0.5),
        :output_dir           => "./output",
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes => "lgl",
        :nop                 => 7,
        :lexact_integration  => false,
        :lsource             => false,
        #---------------------------------------------------------------------------
        # Physical parameters
        #---------------------------------------------------------------------------
        :lvisc               => true,
        :μ                   => [1.0e-2, 1.0e-2],
        #---------------------------------------------------------------------------
        # Mesh
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10_burgers2d.msh",
        #---------------------------------------------------------------------------
        # Output
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loverwrite_output   => true,
        #---------------------------------------------------------------------------
        # Matrix assembly: build the sparse Galerkin Laplacian for the implicit op.
        #---------------------------------------------------------------------------
        :ldss_laplace        => true,
        :lsparse             => true,
        :ldss_differentiation => false,
        #---------------------------------------------------------------------------
        # IMEX method configuration
        #---------------------------------------------------------------------------
        :method              => "RK",
        :delta               => 1,
        :k                   => 3,
        :coeff               => Dict(
                                    :A_RK       => A_RK,
                                    :b_RK       => b_RK,
                                    :c_RK       => c_RK,
                                    :A_RK_tilde => A_RK_tilde,
                                    :b_RK_tilde => b_RK_tilde,
                                    :c_RK_tilde => c_RK_tilde,
                                ),
        :S_fun               => S_fun!,
        :L_fun               => L_fun!,
        :build_L             => build_L,
        :upd_L               => false,
        :lsolve              => lsolve,
        :solver_precision    => TFloat,
        :nl_precision        => TFloat,
    ) #Dict

    return inputs
end
