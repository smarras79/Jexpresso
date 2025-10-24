using Trixi
using SummationByPartsOperators
using OrdinaryDiffEqSSPRK

polydeg = 3

mesh = UniformPeriodicMesh1D(xmin = 0.0, xmax = 2.0, Nx = 12)

D_operator = couple_continuously(
	legendre_derivative_operator(-1.0, 1.0, polydeg + 1),
	mesh
)

D_local = D_operator.D.basis.D

using LinearAlgebra

# Parametri
n_elements = 12
poly_deg = 3
n_dof_local = poly_deg + 1  # 4 DOF per elemento
n_dof_global = n_elements * poly_deg  # 36 DOF (non 37!)

D_local = diagm(inv.(M_loc))*D_loc
D_local[2:3,1:4] *= 10

D_global = zeros(n_dof_global, n_dof_global)
w = D_operator.D.basis.weights
function get_global_indices_periodic(element_id, poly_deg, n_dof_global)
    offset = (element_id - 1) * poly_deg
    indices = offset .+ (1:poly_deg+1)
    indices = mod1.(indices, n_dof_global)
    return indices
end

for k = 1:n_elements
    global_indices = get_global_indices_periodic(k, poly_deg, n_dof_global)
    
    for i = 1:n_dof_local
        for j = 1:n_dof_local
            i_global = global_indices[i]
            j_global = global_indices[j]
			D_global[i_global, j_global] += D_local[j, i]
        end
    end
end

#for k = 1:n_elements
#    global_indices = get_global_indices_periodic(k, poly_deg, n_dof_global)
#    D_global[global_indices, global_indices] .+= D_local[k]
#end
D = Matrix(D_operator)
M_loc = D_operator.D.basis.weights
D_loc = D_operator.D.basis.D
D_loc = diagm(inv.(M_loc)) * D_loc
D_loc[2:3,1:4] *= 10
D_loc = D_operator.D.basis.D

x = grid(D_operator) |> collect
equations = CompressibleEulerEquations1D(1004/718) # gamma = 1.4
u0 = zeros(3, length(x))
let
	for i in eachindex(x)
		if abs(x[i] - 1) < 0.25
			rho = 1.1691
			v = 0.1882 * sign(x[i] - 1)
			p = 1.245
		else
			rho = 1.0
			v = 0.0
			p = 1.0
		end
		u0[1, i] = rho
		u0[2, i] = rho * v
		u0[3, i] = 0.5 * rho * v^2 + p / (equations.gamma - 1)
	end
	u0
end
function rhs!(du, u, parameters, t)
	(; equations, D, volume_flux) = parameters

	for i in axes(du, 2)
		u_i = SVector(u[1, i], u[2, i], u[3, i])
		du_i = zero(u_i)
		
		for j in axes(du, 2)
			u_j = SVector(u[1, j], u[2, j], u[3, j])

			direction = 1 # x direction
			f = volume_flux(u_i, u_j, direction, equations)
			du_i += 2 * D[i, j] * f
		end

		du[1, i] = du_i[1]
		du[2, i] = du_i[2]
		du[3, i] = du_i[3]
	end

#	@show du[2,1]
#	throw(error)  
	return nothing
end

function rhs_local!(du, u, parameters, t)
    (; equations, D_loc, M_loc, M_inv_global, volume_flux) = parameters
        n_elements = 12
	n_nodes_per_element = 4
    
    for elem in 1:n_elements
        for i_loc in 1:n_nodes_per_element
            i_glob = (elem - 1) * (n_nodes_per_element - 1) + i_loc
            i_glob = mod1(i_glob, 36)  # Periodicità
            u_i = SVector(u[1, i_glob], u[2, i_glob], u[3, i_glob])
            du_i = zero(u_i)
            
            for j_loc in 1:n_nodes_per_element
                j_glob = (elem - 1) * (n_nodes_per_element - 1) + j_loc
                j_glob = mod1(j_glob, 36)
                
                u_j = SVector(u[1, j_glob], u[2, j_glob], u[3, j_glob])
                direction = 1
                f = volume_flux(u_i, u_j, direction, equations)
                
                du_i += 2 * D_loc[i_loc, j_loc] * f
            end
            
			du[1, i_glob] += du_i[1]/M_loc[i_loc]
			du[2, i_glob] += du_i[2]/M_loc[i_loc]
			du[3, i_glob] += du_i[3]/M_loc[i_loc]
        end
	end
	#@show du[2,1]
    return nothing
end

tspan = (0.0, 2.0)

volume_flux = flux_ranocha # flux_ranocha (EC), flux_shima_etal (KEP), flux_central

M_loc = D_operator.D.basis.weights
n_elements = 12
n_dof_global = 36
jacobian = (2/ n_elements) / 2 

M_global_diag = zeros(n_dof_global)
for elem in 1:n_elements
    for i_loc in 1:4
        i_glob = (elem - 1) * 3 + i_loc
        i_glob = mod1(i_glob, 36)
        M_global_diag[i_glob] += jacobian * M_loc[i_loc]
    end
end

M_inv_global = 1 ./ M_global_diag
parameters = (; equations, D, D_loc, volume_flux, M_loc, M_inv_global)

ode = ODEProblem(rhs_local!, u0, tspan, parameters)

sol = solve(ode, SSPRK54(); dt = 1.6781e-04, adaptive = false)

M = mass_matrix(D_operator)

let
	# total entropy over time
	U = zeros(length(sol.t)+1)
	
	for n in eachindex(sol.t)
		U_n = 0.0
		U_0 = 0.0
		k = 0
		for i in axes(sol.u[n], 2)
			k = k+1
			u_i = SVector(sol.u[n][1, i], sol.u[n][2, i], sol.u[n][3, i])
			U_n += M[i, i] * entropy(u_i, equations)
			u_i = SVector(ode.u0[1, i], ode.u0[2, i], ode.u0[3, i])
			U_0 += M[i, i] * entropy(u_i, equations)
		end
		U[n] = abs(U_n)
		if n == length(sol.t)
		U[n+1] = abs(U_0)
		end
	end
	@show U[end]
	@show U[1]
	@show U[end-1]
	@show abs(U[1]-U[end-1])
	@show maximum(U) - minimum(U)
end

