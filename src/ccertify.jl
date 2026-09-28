# Rigorous bounds for the complex hierarchy, in the sense of Oustry, D'Ambrosio, Liberti and Ruiz
# (PSCC 2022), transposed to complex numbers and certified NATIVELY on the Hermitian blocks.
#
# The complex SOS form (Josz & Molzahn (3.11)) is
#
#     sup  lambda   s.t.   f - lambda = sum_i sigma_i g_i,   sigma_i in Sigma_{d-k_i}[z],
#
# and a sum of squared moduli sigma_i is exactly p^H X_i p with X_i Hermitian PSD over the basis
# p. Matching coefficients of conj(z)^alpha z^beta gives one identity row per monomial pair:
#
#     C_{alpha,beta}(theta) := sum_i <X_i, G_i^{alpha,beta}> + sum_j nu_j h_j^{alpha,beta}
#                              + sum_c <mu_c, cone_c^{alpha,beta}> + lambda*[alpha=beta=empty]
#                            = f_{alpha,beta}.
#
# For ANY theta, write r = f - C(theta). Then for every feasible z,
#
#     f(z) = lambda + sum_i sigma_i(z) g_i(z) + sum_c <mu_c, ...> + sum_{a,b} r_{a,b} conj(z^a) z^b
#
# and since g_i >= 0 on the feasible set, sigma_i(z) >= lambda_min(X_i) ||p_i(z)||^2, and the cone
# multipliers are paired with points of their own (self-dual) cone,
#
#     f(z) >= lambda - sum_{a,b} |r_{a,b}| max|z^{a+b}| + sum_i rho_i min(lambda_min(X_i), 0)
#
# with rho_i >= max over the box of ||p_i(z)||^2 g_i(z). Valid at EVERY theta, so a stalled solve
# still yields a true lower bound.
#
# WHY NATIVE AND NOT ON THE REAL LIFT. A Hermitian H = A + iB lifts to [A -B; B A], whose spectrum
# is that of H with every eigenvalue REPEATED. So lambda_min is unchanged, but the trace -- and
# hence rho, which bounds it -- DOUBLES. Certifying on the lift therefore charges twice the true
# penalty for the same dual point. Working directly with the k x k Hermitian block halves the PSD
# correction, which on these instances is the dominant term.

using LinearAlgebra

struct ComplexSOSCertificate
    rows::Vector{CMono}                     # monomial pairs, in order
    f::Vector{Float64}                      # objective coefficients / scale (real part)
    mb::Vector{Float64}                     # max |z^{alpha+beta}| over the box
    blocks::Vector{Matrix{ComplexF64}}      # Hermitian Gram matrices
    rho::Vector{Float64}                    # trace bound per block
    scale::Float64
    lambda::Float64
    resid::Vector{Float64}                  # r, in the same order as rows
end

"""
    complex_certified_value(cert) -> (bound, psd_term, box_term)

Evaluate the certificate. `psd_term` uses lambda_min of the HERMITIAN block, which is the half that
the real lift would double.
"""
function complex_certified_value(cert::ComplexSOSCertificate)
    box = sum(abs.(cert.resid) .* cert.mb; init = 0.0)
    psd = 0.0
    for (k, X) in enumerate(cert.blocks)
        isempty(X) && continue
        lam = size(X, 1) == 1 ? real(X[1, 1]) : minimum(eigvals(Hermitian(X)))
        lam < 0 && (psd += cert.rho[k] * lam)
    end
    return (cert.lambda + psd - box) * cert.scale, psd * cert.scale, box * cert.scale
end

"Upper bound on |z^alpha| over the box |z_i| <= absmax_i."
_cmb(m::Vector{Int}, absmax::Vector{Float64}) = prod((absmax[v] for v in m); init = 1.0)

"""
    certify_complex(pop, rel; ...) -> (bound, info)

Build and evaluate the certificate from a solved complex moment relaxation. The Gram matrices are
the DUALS of the relaxation's Hermitian PSD constraints -- the moment form's dual IS the SOS form
-- recovered from the real lift by the structure-preserving projection

    X = (D11 + D22)/2 + i (D21 - D12)/2,

which is the adjoint of the lift and therefore the right way back.
"""
function certify_complex(pop::CPOP, rel::ComplexMomentRelaxation)
    isempty(rel.y) && return (NaN, Dict{String,Any}("status" => "no primal solution"))
    haskey(rel.info, "cert_data") ||
        return (NaN, Dict{String,Any}("status" => "relaxation was not built with certify = true"))
    cd = rel.info["cert_data"]::Dict{String,Any}
    rows = cd["rows"]::Vector{CMono}
    fvec = cd["f"]::Vector{Float64}
    resid = cd["resid"]::Vector{Float64}
    blocks = cd["blocks"]::Vector{Matrix{ComplexF64}}
    rho = cd["rho"]::Vector{Float64}
    scale = cd["scale"]::Float64
    lam = cd["lambda"]::Float64
    absmax = [isfinite(a) ? a : 1.0 for a in pop.absmax]
    mb = [_cmb(vcat(a, b), absmax) for (a, b) in rows]
    cert = ComplexSOSCertificate(rows, fvec, mb, blocks, rho, scale, lam, resid)
    bound, psd, box = complex_certified_value(cert)
    info = Dict{String,Any}("status" => "ok", "psd_correction" => psd, "box_correction" => box,
        "n_rows" => length(rows), "max_resid" => maximum(abs.(resid); init = 0.0),
        "lambda" => lam * scale)
    return bound, info
end
