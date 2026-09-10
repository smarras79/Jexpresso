# The 1D ideal MHD equations and the Brio–Wu shock tube

## Equations

Conservative form, $\partial_t q + \partial_x F(q) = 0$ on $x \in (0, 1)$,
with $B_x$ constant (in 1D $\nabla\cdot\mathbf{B} = \partial_x B_x$, so
$B_x = \mathrm{const}$ is the solenoidal constraint and no cleaning field is
needed):

$$
q = \begin{pmatrix}\rho\\ \rho u\\ \rho v\\ \rho E\\ \rho w\\ B_x\\ B_y\\ B_z\end{pmatrix},\qquad
F = \begin{pmatrix}\rho u\\ \rho u^2 + p_T - B_x^2\\ \rho u v - B_x B_y\\ (\rho E + p_T)u - B_x(\mathbf{v}\cdot\mathbf{B})\\ \rho u w - B_x B_z\\ 0\\ B_y u - B_x v\\ B_z u - B_x w\end{pmatrix},
$$

$$
p = (\gamma - 1)\left(\rho E - \tfrac12\rho|\mathbf{v}|^2 - \tfrac12|\mathbf{B}|^2\right),\qquad
p_T = p + \tfrac12|\mathbf{B}|^2,\qquad \gamma = 2 .
$$

Units are those of Dao & Nazarov (2022): magnetic pressure $\tfrac12|\mathbf{B}|^2$
(Heaviside–Lorentz), the same convention as Jexpresso's 2D MHD cases. The
slot order is the 2D cases' $(\rho, \rho u, \rho v, \rho E, \rho w, B_x, B_y, B_z, \psi)$
without $\psi$.

Wave speeds along $x$: with $a^2 = \gamma p/\rho$, $b^2 = |\mathbf{B}|^2/\rho$,
$b_x^2 = B_x^2/\rho$,

$$
c_{f,s}^2 = \tfrac12\left(a^2 + b^2 \pm \sqrt{(a^2 + b^2)^2 - 4a^2 b_x^2}\right),\qquad c_A = |b_x| ,
$$

and the DynSGS cap uses $|u| + c_f$.

## Initial condition (Dao & Nazarov 2022, §5.2)

$(\rho, u, p, B_x, B_y) = (1, 0, 1, 0.75, 1)$ for $x < 0.5$ and
$(0.125, 0, 0.1, 0.75, -1)$ for $x \ge 0.5$; $v = w = B_z = 0$. Final time
$t = 0.1$ on $(0, 1)$, i.e. Brio & Wu's $t = 0.2$ on $(-1, 1)$, the state of
the paper's Fig. 2 (see README.md).

## Solution structure at $t = 0.1$

From left to right (Brio & Wu 1988, Fig. 2; the reference solution here):

1. a **fast rarefaction** moving left, head at $x \approx 0.32$, foot at $\approx 0.42$;
2. a **slow compound wave** at $x \approx 0.47$ — a slow shock with an attached
   slow rarefaction, the feature unique to this problem (the density spike
   to $\approx 0.72$ in Fig. 2 of the paper); it exists because the
   intermediate shock across which $B_y$ changes sign is admissible in the
   planar problem;
3. a **contact discontinuity** at $x \approx 0.56$ ($\rho$ from $\approx 0.70$ to
   $\approx 0.23$ at constant $p$, $u$, $B_y$);
4. a **slow shock** at $x \approx 0.64$;
5. a **fast rarefaction** moving right, between $x \approx 0.84$ and $0.87$.

## Regularization

DynSGS in the conserved form, $\partial_t q + \partial_x F = \partial_x(\nu\,\partial_x q)$
on every slot, with the residual-based $\nu$ of DSGS.md §4 computed by the
1D MHD kernel `compute_dsgs_viscosity!(::DSGS_MHD, ::NSD_1D)`. Departures to
expect from a CG-SEM + artificial viscosity run at 600 points against the
reference: the shocks and the contact are smeared over 2–3 nodes, the
compound wave's density peak is lower than the reference's, and small
overshoots may remain at the fast-rarefaction feet where the residual is
small.
