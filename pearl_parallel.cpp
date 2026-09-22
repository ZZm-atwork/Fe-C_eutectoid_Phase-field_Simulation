/*
MIT License

Copyright (c) 2026 Pearl solver contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

// Pearl v0.4: documented single-file solver, with opt-in 2D filling methods.
// Existing numerical kernels, input keys and legacy initialization remain.
// The digest identifies the inherited numerical model, not this file hash.
// Opt-in filling operations add their own versioned checkpoint signature.
// Actual combined-file hashes are recorded separately in SHA256SUMS.txt.
// See README.md for build instructions and provenance.
#ifndef PEARL_SOURCE_DIGEST
#define PEARL_SOURCE_DIGEST "344eef3b2e1fdbd8a41d81378491b923f6d815d8ef917f8397729c5a57965580"
#endif

// Parallel: select the original MPI compilation path.
#ifndef PEARL_MPI
#define PEARL_MPI 1
#endif

// ============================================================================
// READING MAP — physical symbols, data flow, and units
// ============================================================================
// alpha=0: ferrite; theta=1: cementite; gamma=2: austenite.
// T is temperature [K]; temperature_slope stores derivatives per kelvin of T.
// carbon is atomic/mole fraction; mu is chemical potential [inherited J/mol convention].
// phi is the three phase fractions; chi is dc/dmu; K is D*chi mobility, not D.
// A/B/C input keys describe f(c)=A*c*c+B*c+C; descriptive members hold those values.
// The q2/q1/q0 input/output names remain compatible; they multiply mu^2, mu, and 1.
// A coefficient is not a temperature slope: only temperature_slope is multiplied by delta T.
// 1 thermodynamics -> 2 reference diffusion -> 3 capillarity -> 4 phase kinetics
// 5 serial/MPI grid -> 6 coupled evolution -> 7 input/filling -> 8 I/O -> 9 main.
// Input.in and INPUT_PARAMETERS.md describe controls; FILLING_METHODS.md describes geometry.

// BEGIN ORIGINAL FILE: include/physics.hpp

// ============================================================================
// 1. THERMODYNAMICS — phase free energies, mixture closure and driving force
// ============================================================================
// All phases share a local chemical potential. Carbon is reconstructed from
// phase fractions and susceptibility; the chemical phase force divides by Vm once.
// Three-phase binary Fe-C thermodynamics. c is carbon atomic fraction, not wt%.
// Inherited molar-code convention: f, mu J/mol; D m^2/s; chi mol/J.
// Scoped physics core; not a full replacement for MicroSim.
#include <array>
#include <cmath>
#include <stdexcept>
namespace pearl {
    constexpr int NPHASE = 3;
    enum Phase {
        Alpha = 0, Theta = 1, Gamma = 2
    };
    using Three = std::array<double, NPHASE>;
    using Matrix3 = std::array<Three, NPHASE>;
    // f(carbon) = free_energy_quadratic*carbon^2 + free_energy_linear*carbon
    //             + free_energy_constant. diffusivity is D [m^2/s], independently supplied.
    struct Parabola {
        double free_energy_quadratic, free_energy_linear, free_energy_constant, diffusivity;
        void validate() const {
        if (!std::isfinite(free_energy_quadratic) || !std::isfinite(free_energy_linear) || !std::isfinite(free_energy_constant)             || !std::isfinite(diffusivity) || free_energy_quadratic <= 0.0 || diffusivity < 0.0)
        throw std::invalid_argument("Parabola requires finite A>0, B, C, D>=0");
        }     double f(double c) const { return std::fma(std::fma(free_energy_quadratic, c, free_energy_linear), c, free_energy_constant);
        }     double mu(double c) const { return std::fma(2.0*free_energy_quadratic, c, free_energy_linear);
        }     double c_at_mu(double mu_value) const { return (mu_value-free_energy_linear)/(2.0*free_energy_quadratic);
        }     double chi() const { return 1.0/(2.0*free_energy_quadratic);
        }     double psi(double mu_value) const {
        const double d = mu_value-free_energy_linear;
        return std::fma(-0.25/free_energy_quadratic, d*d, free_energy_constant);
        // f(c_p)-mu*c_p
        }
    };
    using Thermo = std::array<Parabola, NPHASE>;
    inline void validate_phi(const Three& phi, double tolerance=1e-10) {
        double sum=0.0;
        for (double p : phi) {
            if (!std::isfinite(p) || p < -tolerance || p > 1.0+tolerance)             throw std::invalid_argument("Invalid phase fraction");
            sum += p;
        }
        if (std::abs(sum-1.0) > tolerance)         throw std::invalid_argument("Phase fractions must sum to one");
    }
    // Three-phase cubic interpolation; sum(h)=1 on the Gibbs simplex.
    // Do not replace this with three independent smoothstep functions.
    inline Three weights(const Three& phi) {
        const double triple=2.0*phi[0]*phi[1]*phi[2];
        Three h{};
        for (int p=0; p<NPHASE; ++p)         h[p]=phi[p]*phi[p]*(3.0-2.0*phi[p])+triple;
        return h;
    }
    // dh[p][a] = partial h_p / partial phi_a before constraint projection.
    inline Matrix3 weight_derivatives(const Three& phi) {
        const Three dtriple{2.0*phi[1]*phi[2],
        2.0*phi[0]*phi[2],
        2.0*phi[0]*phi[1]};
        Matrix3 dh{};
        for (int p=0; p<NPHASE; ++p)         for (int a=0; a<NPHASE; ++a)             dh[p][a]=dtriple[a]+(p==a ? 6.0*phi[p]*(1.0-phi[p]) : 0.0);
        return dh;
    }
    struct Mixture {
        double susceptibility=0.0;
        // sum(h_p chi_p)
        double c_intercept=0.0;
        // -sum(h_p chi_p B_p)
        double transport=0.0;
        // K=sum(phi_p D_p chi_p), NOT diffusivity itself
    };
    inline Mixture mixture(const Three& phi, const Thermo& phases) {
        const Three h=weights(phi);
        Mixture out;
        for (int p=0; p<NPHASE; ++p) {
            const double chi=phases[p].chi();
            out.susceptibility += h[p]*chi;
            out.c_intercept -= h[p]*chi*phases[p].free_energy_linear;
            // Retain the reviewed linear transport interpolation.
            // h weights thermodynamics; phi weights transport. They need not be identical.
            out.transport += phi[p]*phases[p].diffusivity*chi;
        }
        if (!(out.susceptibility>0.0) || !std::isfinite(out.susceptibility)         || !std::isfinite(out.c_intercept) || !std::isfinite(out.transport)         || out.transport<0.0)         throw std::runtime_error("Invalid mixture coefficients");
        return out;
    }
    inline double composition(double mu, const Mixture& mix) {
        return std::fma(mix.susceptibility, mu, mix.c_intercept);
    }
    inline double reconstruct_mu(double c_total, const Mixture& mix) {
        if (!std::isfinite(c_total)) throw std::runtime_error("Nonfinite total carbon");
        return (c_total-mix.c_intercept)/mix.susceptibility;
    }
    // Chemical derivative with respect to phi (J/m^3); add interface terms and impose the simplex constraint.
    inline Three chemical_gradient(const Three& phi, double mu,                                const Thermo& phases, double Vm) {
        if (!(Vm>0.0) || !std::isfinite(Vm)) throw std::invalid_argument("Invalid Vm");
        const Matrix3 dh=weight_derivatives(phi);
        Three g{};
        for (int a=0; a<NPHASE; ++a)         for (int p=0; p<NPHASE; ++p)             g[a] += phases[p].psi(mu)*dh[p][a]/Vm;
        return g;
    }
    // Molar-code chemical free energy for closure tests: sum(h_p f_p(c_p(mu))).
    inline double chemical_energy(const Three& phi, double c_total,                               const Thermo& phases) {
        const Three h=weights(phi);
        const double mu=reconstruct_mu(c_total, mixture(phi, phases));
        double value=0.0;
        for (int p=0; p<NPHASE; ++p)         value += h[p]*phases[p].f(phases[p].c_at_mu(mu));
        return value;
    }
    // Legacy TEST fixture: archived 2026-08-31 direct F1 coefficients at 989.76 K.
    // Diffusivities belong to the archived example, not new material-property recommendations.
    // Reconcile with the frozen Sol manifest before production use.
    inline Thermo historical_fe_c_989_76() {
        Thermo p{{
        {6.331161979528086e5,  3.374509335088470e3, -5.852000117829494e3, 2e-9},
        {1.315789473684211e7, -6.574594863157894e6,  8.165231004262522e5, 2e-9},
        {2.462118547531691e4,  2.708815789259557e3, -5.817887580401304e3, 1e-9}
        }};
        for (const auto& phase : p) phase.validate();
        return p;
    }

// -----------------------------------------------------------------------------
// TEMPERATURE PREPARATION (executed once, before initialization)
// The user supplies f(c,T_E)=A_E*c*c+B_E*c+C_E and a temperature law.
// A single-temperature parabola does NOT determine its temperature derivative.
// The legacy linear_gp law is q_k(T)=q_k(T_E)+(T-T_E)*dq_k/dT, where
// Psi(mu,T)=q2(T)*mu*mu+q1(T)*mu+q0(T). This reproduces the archived
// Ankit grand-potential interpolation; do not interpolate mu(c,T) independently.
// -----------------------------------------------------------------------------
// Grand potential Psi(mu) = mu_squared_coefficient*mu^2
//                         + mu_coefficient*mu + constant_energy.
// The same layout stores d(Psi coefficients)/dT in Config::temperature_slope.
struct GrandPotentialPolynomial {
    long double mu_squared_coefficient=0, mu_coefficient=0, constant_energy=0;
};
inline GrandPotentialPolynomial to_grand_potential(const Parabola& p) {
    p.validate();
    const long double free_energy_quadratic=p.free_energy_quadratic, free_energy_linear=p.free_energy_linear, free_energy_constant=p.free_energy_constant;
    return {-1/(4*free_energy_quadratic), free_energy_linear/(2*free_energy_quadratic), free_energy_constant-free_energy_linear*free_energy_linear/(4*free_energy_quadratic)};
}
inline Parabola from_grand_potential(const GrandPotentialPolynomial& grand_potential, double diffusivity) {
    if (!std::isfinite(grand_potential.mu_squared_coefficient)||!std::isfinite(grand_potential.mu_coefficient)||!std::isfinite(grand_potential.constant_energy)||grand_potential.mu_squared_coefficient>=0)
        throw std::invalid_argument("temperature law requires finite q2<0, q1, q0");
    const long double free_energy_quadratic=-1/(4*grand_potential.mu_squared_coefficient), free_energy_linear=-grand_potential.mu_coefficient/(2*grand_potential.mu_squared_coefficient);
    Parabola p{static_cast<double>(free_energy_quadratic),static_cast<double>(free_energy_linear),
               static_cast<double>(grand_potential.constant_energy+free_energy_linear*free_energy_linear/(4*free_energy_quadratic)),diffusivity};
    p.validate();
    return p;
}
inline Parabola at_undercooling(const Parabola& reference,
                               const GrandPotentialPolynomial& temperature_slope,
                               double undercooling) {
    if (!std::isfinite(undercooling)||undercooling<0)
        throw std::invalid_argument("undercooling must be finite and nonnegative");
    if (!std::isfinite(temperature_slope.mu_squared_coefficient)||!std::isfinite(temperature_slope.mu_coefficient)||!std::isfinite(temperature_slope.constant_energy))
        throw std::invalid_argument("nonfinite temperature derivative");
    reference.validate();
    if (undercooling==0) return reference; // preserve supplied reference exactly
    auto grand_potential=to_grand_potential(reference);
    const long double temperature_change=-static_cast<long double>(undercooling);
    grand_potential.mu_squared_coefficient+=temperature_change*temperature_slope.mu_squared_coefficient;
    grand_potential.mu_coefficient+=temperature_change*temperature_slope.mu_coefficient;
    grand_potential.constant_energy+=temperature_change*temperature_slope.constant_energy;
    return from_grand_potential(grand_potential,reference.diffusivity);
}

// -----------------------------------------------------------------------------
// PHASE-DIAGRAM TEMPERATURE LAW — reference ABC plus boundary slopes
// This is the binary fixed-curvature MicroSim construction, anchored to all
// supplied Aeq/Beq/Ceq values. It is a separate option from linear_gp.
// Slopes here mean dc_eq/dT [carbon mole fraction/K], NOT dq/dT or dT/dc.
// Zero theta slope represents a stoichiometric boundary without infinity.
// Gamma ABC stays fixed as the common affine energy-reference convention.
// The alpha/gamma and theta/gamma branch compositions move linearly with T;
// product B is linear in delta T, while C contains the required quadratic term.
// -----------------------------------------------------------------------------
struct PhaseDiagramInput {
    Three reference_composition{}; // alpha, theta, gamma contacts at T_eutectoid
    // alpha on a/g, gamma on a/g, theta on t/g, gamma on t/g
    std::array<double,4> composition_slope{};
};

inline Thermo at_phase_diagram_undercooling(const Thermo& reference,
                                          const PhaseDiagramInput& diagram,
                                          double undercooling) {
    if (!std::isfinite(undercooling)||undercooling<0)
        throw std::invalid_argument("undercooling must be finite and nonnegative");
    for (const auto& phase:reference) phase.validate();
    for (double carbon:diagram.reference_composition)
        if (!std::isfinite(carbon)||carbon<0||carbon>1)
            throw std::invalid_argument("phase-diagram reference carbon must be in [0,1]");
    for (double slope:diagram.composition_slope)
        if (!std::isfinite(slope))
            throw std::invalid_argument("phase-diagram slope must be finite dc/dT");
    const long double temperature_change=-static_cast<long double>(undercooling);
    const long double gamma_A=reference[Gamma].free_energy_quadratic;
    const long double gamma_carbon_E=diagram.reference_composition[Gamma];
    Thermo result=reference;
    for (int product=Alpha;product<=Theta;++product) {
        const long double product_A=reference[product].free_energy_quadratic;
        const long double product_carbon_E=diagram.reference_composition[product];
        const long double product_carbon_change=temperature_change*diagram.composition_slope[2*product];
        const long double gamma_carbon_change=temperature_change*diagram.composition_slope[2*product+1];
        const long double product_carbon=product_carbon_E+product_carbon_change;
        const long double gamma_carbon=gamma_carbon_E+gamma_carbon_change;
        if (!std::isfinite(product_carbon)||!std::isfinite(gamma_carbon)||
            product_carbon<0||product_carbon>1||gamma_carbon<0||gamma_carbon>1)
            throw std::invalid_argument("phase-diagram extrapolated carbon outside [0,1]");
        if (undercooling==0) continue; // preserve every supplied reference bit
        // B(T)-B_E = 2*(A_gamma*delta_c_gamma - A_product*delta_c_product).
        result[product].free_energy_linear=static_cast<double>(
            reference[product].free_energy_linear+
            2*(gamma_A*gamma_carbon_change-product_A*product_carbon_change));
        // Difference-of-squares form retains the full MicroSim C(T) law while
        // preserving the supplied reference fit's small tangent residual.
        result[product].free_energy_constant=static_cast<double>(
            reference[product].free_energy_constant+
            product_A*product_carbon_change*(2*product_carbon_E+product_carbon_change)-
            gamma_A*gamma_carbon_change*(2*gamma_carbon_E+gamma_carbon_change));
        result[product].validate();
    }
    return result;
}

}
// namespace pearl


// ============================================================================
// 2. REFERENCE DIFFUSION — conservative face fluxes and conjugate gradients
// ============================================================================
// This small serial reference implementation is retained for the original
// model tests. The production distributed counterpart is in section 6.
// Serial reference conservative backward-Euler diffusion substep.
// X periodic, Y no-flux; uniform grid; arithmetic face mobility.
// This reference substep contains no phase update, anti-trapping, MPI or I/O.
#include <algorithm>
#include <cstddef>
#include <numeric>
#include <vector>
namespace pearl {
    using Field=std::vector<double>;
    using PhaseField=std::vector<Three>;
    struct Grid {
        int nx, ny;
        double dx, dy;
        std::size_t size() const { return static_cast<std::size_t>(nx)*ny;
        }     int index(int x, int y) const { return x+nx*y;
        }     void validate() const {
        if (nx<3 || ny<2 || !(dx>0) || !(dy>0)             || !std::isfinite(dx) || !std::isfinite(dy))
        throw std::invalid_argument("Invalid 2D grid");
        }
    };
    struct Face {
        int left, right;
        double k_over_h2;
    };
    // Store each unique face once; both neighbors use the same conservative flux.
    inline std::vector<Face> make_faces(const Grid& g, const Field& k) {
        std::vector<Face> faces;
        faces.reserve(2*g.size());
        for (int y=0; y<g.ny; ++y) for (int x=0; x<g.nx; ++x) {
            const int i=g.index(x,y), east=g.index((x+1)%g.nx,y);
            faces.push_back({i,east,0.5*(k[i]+k[east])/(g.dx*g.dx)});
            if (y+1<g.ny) {
                const int north=g.index(x,y+1);
                faces.push_back({i,north,0.5*(k[i]+k[north])/(g.dy*g.dy)});
            }
        }
        return faces;
    }
    // A = diag(chi) - dt*div(K grad). Positive chi and symmetric nonnegative K give SPD A.
    inline Field apply_matrix(const Field& x, const Field& chi,                           const std::vector<Face>& faces, double dt) {
        Field y(x.size());
        for (std::size_t i=0; i<x.size(); ++i) y[i]=chi[i]*x[i];
        for (const auto& f : faces) {
            const double flux=dt*f.k_over_h2*(x[f.left]-x[f.right]);
            y[f.left] += flux;
            y[f.right] -= flux;
        }
        return y;
    }
    inline double dot(const Field& a, const Field& b) {
        return std::inner_product(a.begin(),a.end(),b.begin(),0.0);
    }
    struct CGResult {
        Field x;
        int iterations;
        double residual;
    };
    inline CGResult solve_cg(const Field& rhs, Field x, const Field& chi,                          const std::vector<Face>& faces, double dt,                          double rtol=1e-12, double atol=1e-15,                          int max_iterations=5000) {
        const std::size_t n=rhs.size();
        Field diagonal=chi;
        for (const auto& f : faces) {
            diagonal[f.left] += dt*f.k_over_h2;
            diagonal[f.right] += dt*f.k_over_h2;
        }
        const Field matrix_solution=apply_matrix(x,chi,faces,dt);
        Field residual_vector(n),preconditioned_residual(n),search_direction(n);
        for (std::size_t i=0; i<n; ++i) {
            residual_vector[i]=rhs[i]-matrix_solution[i];
            preconditioned_residual[i]=residual_vector[i]/diagonal[i];
            search_direction[i]=preconditioned_residual[i];
        }
        const double tolerance=atol+rtol*std::sqrt(dot(rhs,rhs));
        double residual_preconditioned_dot=dot(residual_vector,preconditioned_residual);
        for (int iteration=0; iteration<=max_iterations; ++iteration) {
            const double residual=std::sqrt(dot(residual_vector,residual_vector));
            if (residual<=tolerance) return {x,iteration,residual};
            if (iteration==max_iterations) break;
            const Field matrix_direction=apply_matrix(search_direction,chi,faces,dt);
            const double direction_matrix_dot=dot(search_direction,matrix_direction);
            if (!(direction_matrix_dot>0.0) || !std::isfinite(direction_matrix_dot) || !std::isfinite(residual_preconditioned_dot))             throw std::runtime_error("CG breakdown; reject time step");
            const double cg_step=residual_preconditioned_dot/direction_matrix_dot;
            for (std::size_t i=0; i<n; ++i) {
                x[i]+=cg_step*search_direction[i];
                residual_vector[i]-=cg_step*matrix_direction[i];
            }
            for (std::size_t i=0; i<n; ++i) preconditioned_residual[i]=residual_vector[i]/diagonal[i];
            const double next_residual_dot=dot(residual_vector,preconditioned_residual), direction_weight=next_residual_dot/residual_preconditioned_dot;
            for (std::size_t i=0; i<n; ++i) search_direction[i]=preconditioned_residual[i]+direction_weight*search_direction[i];
            residual_preconditioned_dot=next_residual_dot;
        }
        throw std::runtime_error("CG did not converge; reject time step");
    }
    struct DiffusionResult {
        Field carbon, mu;
        int iterations;
        double linear_residual;
        double maximum_closure_residual_before_reconstruction;
    };
    // phi_new may be unchanged for a diffusion test, or supplied by the phase substep.
    // c_old -> solve mu_new -> conservative face-flux c update -> reconstruct mu.
    inline DiffusionResult diffuse(const Grid& grid, const Field& c_old,                                const Field& mu_guess, const PhaseField& phi_new,                                const Thermo& phases, double dt) {
        grid.validate();
        if (!(dt>0) || !std::isfinite(dt)) throw std::invalid_argument("Invalid dt");
        const std::size_t n=grid.size();
        if (c_old.size()!=n || mu_guess.size()!=n || phi_new.size()!=n)         throw std::invalid_argument("Field dimensions do not match grid");
        for (const auto& phase: phases) phase.validate();
        Field chi(n), intercept(n), k(n), rhs(n);
        for (std::size_t i=0; i<n; ++i) {
            validate_phi(phi_new[i]);
            if (!std::isfinite(c_old[i]) || !std::isfinite(mu_guess[i]))             throw std::invalid_argument("Nonfinite input field");
            const auto mix=mixture(phi_new[i],phases);
            chi[i]=mix.susceptibility;
            intercept[i]=mix.c_intercept;
            k[i]=mix.transport;
            rhs[i]=c_old[i]-intercept[i];
        }
        const auto faces=make_faces(grid,k);
        auto solve=solve_cg(rhs,mu_guess,chi,faces,dt);
        Field c=c_old;
        for (const auto& f: faces) {
            const double transfer=dt*f.k_over_h2*(solve.x[f.right]-solve.x[f.left]);
            c[f.left] += transfer;
            c[f.right] -= transfer;
        }
        double closure=0.0;
        for (std::size_t i=0; i<n; ++i) {
            closure=std::max(closure,std::abs(c[i]-(intercept[i]+chi[i]*solve.x[i])));
            solve.x[i]=(c[i]-intercept[i])/chi[i];
        }
        return {std::move(c),std::move(solve.x),solve.iterations,solve.residual,closure};
    }
}
// namespace pearl


// ============================================================================
// 3. CAPILLARITY — interface energy and phase-fraction constraints
// ============================================================================
// Equal isotropic sigma, interface thickness epsilon and the triple penalty
// form the original interface functional. Projection keeps phi on the simplex.
// Isotropic three-phase interface energy and consistent discrete derivatives.
// Energy/derivative helpers; they do not by themselves validate coupled growth.
#include <functional>
namespace pearl {
    struct InterfaceParameters {
        double epsilon;
        // m; model length parameter, not the measured 10%-90% width
        double sigma;
        // J/m^2; one sigma for all three pairs in this scoped model
        double triple;
        // J/m^2; triple-phase coefficient, not a ban on triple-junction coexistence
        void validate() const {
        if (!(epsilon>0) || !(sigma>0) || triple<0             || !std::isfinite(epsilon) || !std::isfinite(sigma)             || !std::isfinite(triple))
        throw std::invalid_argument("Invalid interface parameters");
        }
    };
    struct InterfaceResult {
        double energy_sum=0.0;
        // Sum of cell/unique-face densities; multiply by dx*dy to obtain J/m
        PhaseField gradient;
        // Derivative per cell/phase, J/m^3, before simplex projection
    };
    inline InterfaceResult capillary_energy_gradient(     const Grid& grid, const PhaseField& phi, const InterfaceParameters& p) {
        grid.validate();
        p.validate();
        if (phi.size()!=grid.size()) throw std::invalid_argument("Wrong phase field size");
        // Allow off-simplex inputs here for finite-difference checks of partial derivatives.
        // The solver must validate the phase constraint before accepting a step.
        for (const auto& point: phi) for(double value:point)         if(!std::isfinite(value)) throw std::invalid_argument("Nonfinite phase field");
        InterfaceResult out;
        out.gradient.resize(grid.size(),Three{0.0,0.0,0.0});
        const double obstacle=16.0*p.sigma/(std::acos(-1.0)*std::acos(-1.0)*p.epsilon);
        const double triple=p.triple/p.epsilon;
        for(std::size_t i=0;i<grid.size();++i) {
            for(int a=0;a<3;++a) for(int b=a+1;b<3;++b) {
                out.energy_sum += obstacle*phi[i][a]*phi[i][b];
                out.gradient[i][a] += obstacle*phi[i][b];
                out.gradient[i][b] += obstacle*phi[i][a];
            }
            out.energy_sum += triple*phi[i][0]*phi[i][1]*phi[i][2];
            for(int a=0;a<3;++a)             out.gradient[i][a] += triple*phi[i][(a+1)%3]*phi[i][(a+2)%3];
        }
        // q_ab(face)=phi_a,L*phi_b,R/d - phi_a,R*phi_b,L/d。
        // Equivalent antisymmetric combination of face-averaged phi and face gradients.
        auto add_face=[&](int left,int right,double distance) {
            for(int a=0;a<3;++a) for(int b=a+1;b<3;++b) {
                const double pair_gradient=(phi[left][a]*phi[right][b]-phi[right][a]*phi[left][b])/distance;
                out.energy_sum += p.epsilon*p.sigma*pair_gradient*pair_gradient;
                const double prefactor=2.0*p.epsilon*p.sigma*pair_gradient/distance;
                out.gradient[left][a] += prefactor*phi[right][b];
                out.gradient[right][a] -= prefactor*phi[left][b];
                out.gradient[left][b] -= prefactor*phi[right][a];
                out.gradient[right][b] += prefactor*phi[left][a];
            }
        };
        for(int y=0;y<grid.ny;++y) for(int x=0;x<grid.nx;++x) {
            const int i=grid.index(x,y);
            add_face(i,grid.index((x+1)%grid.nx,y),grid.dx);
            if(y+1<grid.ny) add_face(i,grid.index(x,y+1),grid.dy);
        }
        return out;
    }
    // Euclidean projection to phi>=0 and sum(phi)=1; this does not clip carbon.
    // A three-element sort; no matrix library is needed.
    inline Three project_simplex(const Three& trial) {
        for(double value:trial)         if(!std::isfinite(value)) throw std::invalid_argument("Nonfinite projection input");
        Three sorted=trial;
        std::sort(sorted.begin(),sorted.end(),std::greater<double>());
        double sum=0.0, threshold=0.0;
        int active=0;
        for(int j=0;j<3;++j) {
            sum+=sorted[j];
            const double candidate=(sum-1.0)/(j+1);
            if(sorted[j]>candidate) {
                active=j+1;
                threshold=candidate;
            }
        }
        if(active==0) throw std::runtime_error("Simplex projection failed");
        Three result{};
        for(int a=0;a<3;++a) result[a]=std::max(0.0,trial[a]-threshold);
        return result;
    }
}
// namespace pearl


// ============================================================================
// 4. PHASE KINETICS — relaxation times and explicit phase proposal
// ============================================================================
// Pair relaxation times tau determine the phase rate. This section does not
// fit tau, change diffusivity, or add anti-trapping.
// Explicit obstacle-constrained phase step. No carbon clipping, no Laplacian gate.
namespace pearl {
    struct PhaseParameters {
        InterfaceParameters interface;
        Three tau;
        // 01, 02, 12, J s/m^4
        double Vm;
        void validate() const {
        interface.validate();
        for (double t:tau) if(!(t>0)||!std::isfinite(t))
        throw std::invalid_argument("tau must be finite and positive");
        if(!(Vm>0)||!std::isfinite(Vm)) throw std::invalid_argument("Invalid Vm");
        }
    };
    inline double local_tau(const Three& phi, const Three& tau) {
        const Three w{phi[0]*phi[1], phi[0]*phi[2], phi[1]*phi[2]};
        const double s=w[0]+w[1]+w[2];
        return s>1e-30 ? (w[0]*tau[0]+w[1]*tau[1]+w[2]*tau[2])/s                    : *std::min_element(tau.begin(),tau.end());
    }
    struct PhaseResult {
        PhaseField phi;
        double maximum_change=0;
    };
    inline PhaseResult phase_step(const Grid& g, const PhaseField& old_phi,                               const Field& mu, const Thermo& thermo,                               const PhaseParameters& p, double dt) {
        p.validate();
        if (!(dt>0)||!std::isfinite(dt)||mu.size()!=g.size())         throw std::invalid_argument("Bad phase step");
        auto cap=capillary_energy_gradient(g,old_phi,p.interface);
        PhaseResult out;
        out.phi.resize(g.size());
        for (std::size_t i=0;i<g.size();++i) {
            const Three chem=chemical_gradient(old_phi[i],mu[i],thermo,p.Vm);
            Three force{}, trial{};
            double mean=0;
            for(int a=0;a<3;++a) {
                force[a]=cap.gradient[i][a]+chem[a];
                mean+=force[a]/3;
            }
            const double dtmob=dt/(local_tau(old_phi[i],p.tau)*p.interface.epsilon);
            for(int a=0;a<3;++a) trial[a]=old_phi[i][a]-dtmob*(force[a]-mean);
            out.phi[i]=project_simplex(trial);
            for(int a=0;a<3;++a) out.maximum_change=std::max(out.maximum_change,std::abs(out.phi[i][a]-old_phi[i][a]));
        }
        return out;
    }
    inline double profile(double distance,double epsilon) {
        // Phase=1 on negative side. Compact equilibrium obstacle profile.
        const double pi=std::acos(-1.0), half=pi*pi*epsilon/8;
        if(distance<=-half) return 1;
        if(distance>=half) return 0;
        return 0.5*(1-std::sin(4*distance/(pi*epsilon)));
    }
}
// namespace pearl


// END ORIGINAL FILE: include/physics.hpp

// BEGIN ORIGINAL FILE: include/solver.hpp

// ============================================================================
// 5. SPATIAL GRID — decomposition, halo exchange, and global I/O buffers
// ============================================================================
// Global coordinates identify the same cell in serial and every MPI layout.
// X is periodic; Y has zero normal flux. Each local array includes one ghost cell.
// One-cell Cartesian halos; x periodic, y zero-flux. Same kernels in serial/MPI.
// Compile with PEARL_MPI to enable MPI. No per-rank global-field replication.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef PEARL_MPI
#include <mpi.h>
#endif
namespace pearl {
    struct Array {
        int nx=0, ny=0, nc=1;
        std::vector<double> v;
        Array()=default;
        Array(int x,int y,int n=1):nx(x),ny(y),nc(n),v(static_cast<std::size_t>(x+2)*(y+2)*n,0) {
        }
        std::size_t at(int x,int y,int c=0)const{return (static_cast<std::size_t>(y)*(nx+2)+x)*nc+c;
        }     double& operator()(int x,int y,int c=0) {
            return v[at(x,y,c)];
        }
        double operator()(int x,int y,int c=0)const{return v[at(x,y,c)];
        }
    };
    class Domain {
        public:     int gx,gy,nx,ny,x0,y0,px,py,rank=0,size=1,rx=0,ry=0;
        double dx,dy;
#ifdef PEARL_MPI
        MPI_Comm comm=MPI_COMM_NULL;
        int west,east,south,north;
#endif
        Domain(int x,int y,double hx,double hy,int xp=1,int yp=1):gx(x),gy(y),px(xp),py(yp),dx(hx),dy(hy) {
            if(x<3||y<2||hx<=0||hy<=0||!std::isfinite(hx+hy)||xp<1||yp<1)             throw std::invalid_argument("invalid grid");
#ifdef PEARL_MPI
            MPI_Comm_size(MPI_COMM_WORLD,&size);
            if(xp*yp!=size)throw std::invalid_argument("MPI size != px*py");
            int dims[2]{yp,xp},periods[2]{0,1};
            MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,0,&comm);
            MPI_Comm_rank(comm,&rank);
            int coords[2];
            MPI_Cart_coords(comm,rank,2,coords);
            ry=coords[0];
            rx=coords[1];
            MPI_Cart_shift(comm,1,1,&west,&east);
            MPI_Cart_shift(comm,0,1,&south,&north);
#else
            if(xp!=1||yp!=1)throw std::invalid_argument("serial build requires px=py=1");
#endif
            if(x<xp||y<yp)throw std::invalid_argument("empty MPI blocks are disallowed");
            nx=x/xp+(rx<x%xp);
            ny=y/yp+(ry<y%yp);
            x0=rx*(x/xp)+std::min(rx,x%xp);
            y0=ry*(y/yp)+std::min(ry,y%yp);
        }
        Domain(const Domain&)=delete;
        Domain& operator=(const Domain&)=delete;
        ~Domain() {
#ifdef PEARL_MPI
            int done=0;
            MPI_Finalized(&done);
            if(!done&&comm!=MPI_COMM_NULL)MPI_Comm_free(&comm);
#endif
        }
        double sum(double a)const{
#ifdef PEARL_MPI
        double b;
        MPI_Allreduce(&a,&b,1,MPI_DOUBLE,MPI_SUM,comm);
        return b;
#else
        return a;
#endif
        }     double maximum(double a)const{
#ifdef PEARL_MPI
        double b;
        MPI_Allreduce(&a,&b,1,MPI_DOUBLE,MPI_MAX,comm);
        return b;
#else
        return a;
#endif
        }     double minimum(double a)const{return -maximum(-a);
        }     bool all(bool a)const{return minimum(a?1.0:0.0)>0.5;
        }     void barrier()const{
#ifdef PEARL_MPI
        MPI_Barrier(comm);
#endif
        }     bool south_edge()const{return y0==0;
        }     bool north_edge()const{return y0+ny==gy;
        }     Array array(int nc=1)const{return Array(nx,ny,nc);
        }     void halo(Array& a)const{
        if(a.nx!=nx||a.ny!=ny)throw std::invalid_argument("halo shape");
#ifdef PEARL_MPI
        std::vector<double> send(static_cast<std::size_t>(ny)*a.nc),recv(send.size());
        for(int y=1;y<=ny;++y)for(int c=0;c<a.nc;++c)send[(y-1)*a.nc+c]=a(nx,y,c);
        MPI_Sendrecv(send.data(),static_cast<int>(send.size()),MPI_DOUBLE,east,100,                      recv.data(),static_cast<int>(recv.size()),MPI_DOUBLE,west,100,comm,MPI_STATUS_IGNORE);
        for(int y=1;y<=ny;++y)for(int c=0;c<a.nc;++c)a(0,y,c)=recv[(y-1)*a.nc+c];
        for(int y=1;y<=ny;++y)for(int c=0;c<a.nc;++c)send[(y-1)*a.nc+c]=a(1,y,c);
        MPI_Sendrecv(send.data(),static_cast<int>(send.size()),MPI_DOUBLE,west,101,                      recv.data(),static_cast<int>(recv.size()),MPI_DOUBLE,east,101,comm,MPI_STATUS_IGNORE);
        for(int y=1;y<=ny;++y)for(int c=0;c<a.nc;++c)a(nx+1,y,c)=recv[(y-1)*a.nc+c];
        MPI_Sendrecv(&a.v[a.at(1,ny)],nx*a.nc,MPI_DOUBLE,north,102,                      &a.v[a.at(1,0)],nx*a.nc,MPI_DOUBLE,south,102,comm,MPI_STATUS_IGNORE);
        MPI_Sendrecv(&a.v[a.at(1,1)],nx*a.nc,MPI_DOUBLE,south,103,                      &a.v[a.at(1,ny+1)],nx*a.nc,MPI_DOUBLE,north,103,comm,MPI_STATUS_IGNORE);
#else
        for(int y=1;y<=ny;++y)for(int c=0;c<a.nc;++c) {
            a(0,y,c)=a(nx,y,c);
            a(nx+1,y,c)=a(1,y,c);
        }
#endif
        if(south_edge())for(int x=1;x<=nx;++x)for(int c=0;c<a.nc;++c)a(x,0,c)=a(x,1,c);
        if(north_edge())for(int x=1;x<=nx;++x)for(int c=0;c<a.nc;++c)a(x,ny+1,c)=a(x,ny,c);
        }      // Global checkpoint/output gather on rank 0 only. Ghosts never saved.
        std::vector<double> gather(const Array& a)const{
        std::vector<double> local;
        local.reserve(static_cast<std::size_t>(nx)*ny*a.nc);
        for(int y=1;y<=ny;++y)for(int x=1;x<=nx;++x)for(int c=0;c<a.nc;++c)local.push_back(a(x,y,c));
#ifdef PEARL_MPI
        int localcount=static_cast<int>(local.size());
        std::vector<int> counts(size),offsets(size);
        MPI_Gather(&localcount,1,MPI_INT,counts.data(),1,MPI_INT,0,comm);
        int total=0;
        if(rank==0)for(int r=0;r<size;++r) {
            offsets[r]=total;
            total+=counts[r];
        }
        std::vector<double> packed(rank==0?total:0);
        MPI_Gatherv(local.data(),localcount,MPI_DOUBLE,packed.data(),counts.data(),offsets.data(),MPI_DOUBLE,0,comm);
        if(rank!=0)return {};
        std::vector<double> global(static_cast<std::size_t>(gx)*gy*a.nc);
        for(int r=0;r<size;++r) {
            int coord[2];
            MPI_Cart_coords(comm,r,2,coord);
            int yr=coord[0],xr=coord[1];
            int xn=gx/px+(xr<gx%px),yn=gy/py+(yr<gy%py),xs=xr*(gx/px)+std::min(xr,gx%px),ys=yr*(gy/py)+std::min(yr,gy%py);
            for(int y=0;y<yn;++y)for(int x=0;x<xn;++x)for(int c=0;c<a.nc;++c)                 global[(static_cast<std::size_t>(ys+y)*gx+xs+x)*a.nc+c]=packed[offsets[r]+(y*xn+x)*a.nc+c];
        }
        return global;
#else
        return local;
#endif
        }      // Only restart/rare moving-window remaps broadcast a global buffer.
        // This is intentionally simple; normal timestepping stays distributed.
        void load_global(Array& a,std::vector<double> global)const{
        const std::size_t n=static_cast<std::size_t>(gx)*gy*a.nc;
#ifdef PEARL_MPI
        if(n>static_cast<std::size_t>(std::numeric_limits<int>::max()))throw std::runtime_error("global buffer too large");
        if(rank!=0)global.resize(n);
        MPI_Bcast(global.data(),static_cast<int>(n),MPI_DOUBLE,0,comm);
#endif
        if(global.size()!=n)throw std::runtime_error("global field shape mismatch");
        for(int y=1;y<=ny;++y)for(int x=1;x<=nx;++x)for(int c=0;c<a.nc;++c)
        a(x,y,c)=global[(static_cast<std::size_t>(y0+y-1)*gx+x0+x-1)*a.nc+c];
        halo(a);
        }
    };
    inline double array_dot(const Domain& d,const Array& a,const Array& b) {
        double s=0;
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x)s+=a(x,y)*b(x,y);
        return d.sum(s);
    }
}
// namespace pearl


// ============================================================================
// 6. COUPLED EVOLUTION — initialization smoothing and accepted time steps
// ============================================================================
// Capillarity smoothing is preprocessing pseudo-time. Physical evolution uses
// a phase proposal followed by conservative implicit diffusion and acceptance checks.
// Distributed isotropic phase evolution + exactly closed conservative diffusion.
// First-order Lie splitting; no anti-trapping correction in this scoped model.
#include <limits>
namespace pearl {
    // Distributed physical fields plus accepted-step and moving-window bookkeeping.
    // mu [J/mol], carbon [atomic fraction], phi [unitless, sum=1].
    struct State {
        Array phi,carbon,mu;
        double time=0,dt=0,shift_distance=0,initial_carbon=0,carbon_added=0,carbon_removed=0;
        long long step=0,shifts=0,rejected=0;
        explicit State(const Domain& d):phi(d.array(3)),carbon(d.array()),mu(d.array()) {
        }
    };
    // Numerical acceptance limits. Their values and the original retry policy are unchanged.
    // CG is a linear solver; its step length is unrelated to alpha ferrite.
    struct Controls {
        double dt_max=1e-3,dt_min=1e-12,max_phase_change=.04;
        double cg_rtol=1e-11,cg_atol=1e-14,closure_tol=1e-10;
        double energy_rtol=1e-11,energy_atol=1e-12;
        int cg_max=2000,max_retries=24;
        bool check_energy=true,phase_only=false;
    };
    struct StepInfo {
        double accepted_dt=0,max_phase_change=0,linear_residual=0,closure=0,energy_change=0;
        int cg_iterations=0,retries=0;
    };
    inline Three cell_phi(const Array& a,int x,int y) {
        return {a(x,y,0),a(x,y,1),a(x,y,2)};
    }
    inline void set_phi(Array& a,int x,int y,const Three& p) {
        for(int i=0;i<3;++i)a(x,y,i)=p[i];
    }
    inline double total_carbon(const Domain& d,const State& s) {
        double v=0;
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x)v+=s.carbon(x,y);
        return d.sum(v)*d.dx*d.dy;
    }
    inline void close_mu(const Domain& d,State& s,const Thermo& f) {
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x)s.mu(x,y)=reconstruct_mu(s.carbon(x,y),mixture(cell_phi(s.phi,x,y),f));
        d.halo(s.mu);
    }
    inline double reduced_energy(const Domain& d,State& s,const Thermo& f,const PhaseParameters& pp,double mu_reference,bool phase_only) {
        d.halo(s.phi);
        const double eps=pp.interface.epsilon,sigma=pp.interface.sigma;
        const double obstacle=16*sigma/(std::acos(-1.)*std::acos(-1.)*eps),triple=pp.interface.triple/eps;
        const Three psi{f[0].psi(mu_reference)-f[2].psi(mu_reference),f[1].psi(mu_reference)-f[2].psi(mu_reference),0};
        double total=0;
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
            auto p=cell_phi(s.phi,x,y);
            auto h=weights(p);
            double chem=0;
            for(int a=0;a<3;++a)chem+=h[a]*psi[a];
            if(!phase_only) {
                auto m=mixture(p,f);
                double dc=s.carbon(x,y)-composition(mu_reference,m);
                chem+=.5*dc*dc/m.susceptibility;
            }
            total+=chem/pp.Vm+obstacle*(p[0]*p[1]+p[0]*p[2]+p[1]*p[2])+triple*p[0]*p[1]*p[2];
            for(int dir=0;dir<2;++dir) {
                if(dir==1 && y==d.ny && d.north_edge())continue;
                auto neighbor_phase=cell_phi(s.phi,x+(dir==0),y+(dir==1));
                double dist=dir==0?d.dx:d.dy;
                for(int a=0;a<3;++a)for(int b=a+1;b<3;++b) {
                    double z=(p[a]*neighbor_phase[b]-neighbor_phase[a]*p[b])/dist;
                    total+=eps*sigma*z*z;
                }
            }
        }
        return d.sum(total)*d.dx*d.dy;
    }
    inline double interface_energy_only(const Domain& d,State& s,const PhaseParameters& pp) {
        d.halo(s.phi);
        const double eps=pp.interface.epsilon,sigma=pp.interface.sigma;
        const double obstacle=16*sigma/(std::acos(-1.)*std::acos(-1.)*eps),triple=pp.interface.triple/eps;
        double total=0;
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
            auto p=cell_phi(s.phi,x,y);
            total+=obstacle*(p[0]*p[1]+p[0]*p[2]+p[1]*p[2])+triple*p[0]*p[1]*p[2];
            for(int dir=0;dir<2;++dir) {
                if(dir==1 && y==d.ny && d.north_edge())continue;
                auto neighbor_phase=cell_phi(s.phi,x+(dir==0),y+(dir==1));
                double dist=dir==0?d.dx:d.dy;
                for(int a=0;a<3;++a)for(int b=a+1;b<3;++b) {
                    double z=(p[a]*neighbor_phase[b]-neighbor_phase[a]*p[b])/dist;
                    total+=eps*sigma*z*z;
                }
            }
        }
        return d.sum(total)*d.dx*d.dy;
    }

    struct SmoothReport {
        int steps=0;
        bool converged=false;
        double energy_initial=0;
        double energy_final=0;
        double final_relative_energy_change=0;
        double final_max_phase_change=0;
        std::array<double,3> phase_fraction_before{0,0,0};
        std::array<double,3> phase_fraction_after{0,0,0};
    };

    inline std::array<double,3> phase_fractions(const Domain& d,const State& s) {
        std::array<double,3> local{0,0,0};
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x)
            for(int a=0;a<3;++a)local[a]+=s.phi(x,y,a);
        std::array<double,3> out{};
        const double n=static_cast<double>(d.gx)*d.gy;
        for(int a=0;a<3;++a)out[a]=d.sum(local[a])/n;
        return out;
    }

    // Numerical preprocessing only: relax a sharp phase geometry with the same
    // isotropic interface/obstacle energy, but with NO chemical driving and NO
    // carbon transport. The iteration count is pseudo-time, not physical time.
    // Physical time and solver step remain zero after this function returns.
    // ---- 6A. Capillarity-only preprocessing: no carbon diffusion or physical time ----
    inline SmoothReport smooth_interface_only(const Domain& d,State& s,
                                              const PhaseParameters& pp,
                                              int max_steps,
                                              int min_steps,
                                              double energy_rtol,
                                              double max_phase_change) {
        if(max_steps<1||min_steps<0||min_steps>max_steps||!(energy_rtol>0)||
           !(max_phase_change>0)||max_phase_change>.2)
            throw std::invalid_argument("invalid sharp-smoothing controls");
        pp.validate();
        SmoothReport report;
        report.phase_fraction_before=phase_fractions(d,s);
        report.energy_initial=interface_energy_only(d,s,pp);
        double previous=report.energy_initial;
        const double eps=pp.interface.epsilon,sigma=pp.interface.sigma;
        const double obs=16*sigma/(std::acos(-1.)*std::acos(-1.)*eps),tr=pp.interface.triple/eps;
        // Natural explicit gradient-flow scale for the discrete capillary operator.
        const double base_scale=.12*std::min(d.dx*d.dx,d.dy*d.dy)/(eps*sigma);
        Array force=d.array(3),trial_phi=d.array(3);
        for(int step=1;step<=max_steps;++step) {
            d.halo(s.phi);
            double force_max=0;
            for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
                Three p=cell_phi(s.phi,x,y),f{};
                for(int a=0;a<3;++a)
                    f[a]=obs*(p[(a+1)%3]+p[(a+2)%3])+tr*p[(a+1)%3]*p[(a+2)%3];
                const int ox[4]{-1,1,0,0},oy[4]{0,0,-1,1};
                for(int dir=0;dir<4;++dir) {
                    if((dir==2&&y==1&&d.south_edge())||(dir==3&&y==d.ny&&d.north_edge()))continue;
                    auto neighbor_phase=cell_phi(s.phi,x+ox[dir],y+oy[dir]);
                    const double dist=dir<2?d.dx:d.dy;
                    for(int a=0;a<3;++a)for(int b=a+1;b<3;++b) {
                        const double w=2*eps*sigma*(p[a]*neighbor_phase[b]-neighbor_phase[a]*p[b])/(dist*dist);
                        f[a]+=w*neighbor_phase[b];
                        f[b]-=w*neighbor_phase[a];
                    }
                }
                const double mean=(f[0]+f[1]+f[2])/3;
                for(int a=0;a<3;++a) {
                    force(x,y,a)=f[a]-mean;
                    force_max=std::max(force_max,std::abs(force(x,y,a)));
                }
            }
            force_max=d.maximum(force_max);
            if(!std::isfinite(force_max))throw std::runtime_error("nonfinite smoothing force");
            if(force_max==0) {
                report.steps=step-1;
                report.converged=true;
                report.final_relative_energy_change=0;
                report.final_max_phase_change=0;
                break;
            }
            double scale=std::min(base_scale,max_phase_change/force_max);
            bool accepted=false;
            double new_energy=previous,actual_change=0;
            for(int backtrack=0;backtrack<24&&!accepted;++backtrack) {
                actual_change=0;
                for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
                    Three p=cell_phi(s.phi,x,y),z{};
                    for(int a=0;a<3;++a)z[a]=p[a]-scale*force(x,y,a);
                    z=project_simplex(z);
                    for(int a=0;a<3;++a)actual_change=std::max(actual_change,std::abs(z[a]-p[a]));
                    set_phi(trial_phi,x,y,z);
                }
                actual_change=d.maximum(actual_change);
                State candidate=s;
                candidate.phi=trial_phi;
                new_energy=interface_energy_only(d,candidate,pp);
                const double etol=1e-14*std::max({1.0,std::abs(previous),std::abs(report.energy_initial)});
                if(std::isfinite(new_energy)&&new_energy<=previous+etol)accepted=true;
                else scale*=.5;
            }
            if(!accepted)throw std::runtime_error("sharp smoothing could not find an energy-decreasing step");
            s.phi=trial_phi;
            d.halo(s.phi);
            const double rel=std::abs(previous-new_energy)/std::max(std::abs(previous),1e-300);
            report.steps=step;
            report.energy_final=new_energy;
            report.final_relative_energy_change=rel;
            report.final_max_phase_change=actual_change;
            previous=new_energy;
            if(step>=min_steps && rel<energy_rtol && actual_change<1e-5) {
                report.converged=true;
                break;
            }
        }
        if(report.steps==0)report.energy_final=report.energy_initial;
        report.phase_fraction_after=phase_fractions(d,s);
        return report;
    }

    // ---- 6B. Explicit phase proposal: interface force + chemical driving force ----
    inline double propose_phi(const Domain& d,State& old,State& trial,const Thermo& f,const PhaseParameters& pp,double dt) {
        d.halo(old.phi);
        const double eps=pp.interface.epsilon,sigma=pp.interface.sigma;
        const double obs=16*sigma/(std::acos(-1.)*std::acos(-1.)*eps),tr=pp.interface.triple/eps;
        double change=0;
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
            Three p=cell_phi(old.phi,x,y),force{};
            // Subtract common psi_gamma: removes an irrelevant constrained force.
            const double psi0=f[2].psi(old.mu(x,y));
            auto dh=weight_derivatives(p);
            for(int a=0;a<3;++a) {
                force[a]=obs*(p[(a+1)%3]+p[(a+2)%3])+tr*p[(a+1)%3]*p[(a+2)%3];
                for(int b=0;b<3;++b)force[a]+=(f[b].psi(old.mu(x,y))-psi0)*dh[b][a]/pp.Vm;
            }
            const int ox[4]{-1,1,0,0},oy[4]{0,0,-1,1};
            for(int dir=0;dir<4;++dir) {
                if((dir==2&&y==1&&d.south_edge())||(dir==3&&y==d.ny&&d.north_edge()))continue;
                auto neighbor_phase=cell_phi(old.phi,x+ox[dir],y+oy[dir]);
                double dist=dir<2?d.dx:d.dy;
                for(int a=0;a<3;++a)for(int b=a+1;b<3;++b) {
                    double w=2*eps*sigma*(p[a]*neighbor_phase[b]-neighbor_phase[a]*p[b])/(dist*dist);
                    force[a]+=w*neighbor_phase[b];
                    force[b]-=w*neighbor_phase[a];
                }
            }
            double mean=(force[0]+force[1]+force[2])/3,scale=dt/(local_tau(p,pp.tau)*eps);
            Three z;
            for(int a=0;a<3;++a)z[a]=p[a]-scale*(force[a]-mean);
            // Model/input validation is performed before collectives; finite trial is checked below.
            bool finite=true;
            for(double a:z)finite=finite&&std::isfinite(a);
            if(!finite) {
                change=std::numeric_limits<double>::infinity();
                set_phi(trial.phi,x,y,p);
                continue;
            }
            z=project_simplex(z);
            for(int a=0;a<3;++a)change=std::max(change,std::abs(z[a]-p[a]));
            set_phi(trial.phi,x,y,z);
        }
        return d.maximum(change);
    }
    struct DiffusionInfo {
        bool ok=false;
        int it=0;
        double residual=0,closure=0;
    };
    // ---- 6C. Conservative diffusion: solve for mu-mu_reference, then close carbon ----
    // The true CG residual is recomputed for convergence; no per-step residual history is saved.
    inline DiffusionInfo implicit_diffusion(const Domain& d,const State& old,State& trial,const Thermo& f,double dt,double mu_reference,const Controls& ctrl) {
        Array chi=d.array(),intercept=d.array(),K=d.array(),rhs=d.array(),mu_offset=d.array(),preconditioner_diagonal=d.array();
        Array residual_vector=d.array(),preconditioned_residual=d.array(),search_direction=d.array(),matrix_direction=d.array();
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
            auto m=mixture(cell_phi(trial.phi,x,y),f);
            chi(x,y)=m.susceptibility;
            intercept(x,y)=composition(mu_reference,m);
            K(x,y)=m.transport;
            rhs(x,y)=old.carbon(x,y)-intercept(x,y);
            mu_offset(x,y)=old.mu(x,y)-mu_reference;
        }
        d.halo(K);
        auto weights_at=[&](int x,int y) {
            std::array<double,4>w{.5*(K(x,y)+K(x-1,y))/(d.dx*d.dx),.5*(K(x,y)+K(x+1,y))/(d.dx*d.dx),
            .5*(K(x,y)+K(x,y-1))/(d.dy*d.dy),.5*(K(x,y)+K(x,y+1))/(d.dy*d.dy)};
            if(y==1&&d.south_edge())w[2]=0;
            if(y==d.ny&&d.north_edge())w[3]=0;
            return w;
        };
        auto apply=[&](Array& in,Array& out) {
            d.halo(in);
            for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
                auto w=weights_at(x,y);
                out(x,y)=chi(x,y)*in(x,y)+dt*(w[0]*(in(x,y)-in(x-1,y))+w[1]*(in(x,y)-in(x+1,y))+w[2]*(in(x,y)-in(x,y-1))+w[3]*(in(x,y)-in(x,y+1)));
            }
        };
        apply(mu_offset,matrix_direction);
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
            auto w=weights_at(x,y);
            preconditioner_diagonal(x,y)=chi(x,y)+dt*(w[0]+w[1]+w[2]+w[3]);
            residual_vector(x,y)=rhs(x,y)-matrix_direction(x,y);
            preconditioned_residual(x,y)=residual_vector(x,y)/preconditioner_diagonal(x,y);
            search_direction(x,y)=preconditioned_residual(x,y);
        }
        const double tolerance=ctrl.cg_atol+ctrl.cg_rtol*std::sqrt(array_dot(d,rhs,rhs));
        double residual_preconditioned_dot=array_dot(d,residual_vector,preconditioned_residual);
        DiffusionInfo info;
        for(int iteration=0;iteration<=ctrl.cg_max;++iteration) {
            info.it=iteration;
            info.residual=std::sqrt(array_dot(d,residual_vector,residual_vector));
            if(!std::isfinite(info.residual))return info;
            if(info.residual<=tolerance) {
                // Recompute TRUE residual. Do not trust recursive residual alone.
                apply(mu_offset,matrix_direction);
                for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x)residual_vector(x,y)=rhs(x,y)-matrix_direction(x,y);
                info.residual=std::sqrt(array_dot(d,residual_vector,residual_vector));
                if(info.residual<=tolerance*1.05) {
                    info.ok=true;
                    break;
                }
                for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
                    preconditioned_residual(x,y)=residual_vector(x,y)/preconditioner_diagonal(x,y);
                    search_direction(x,y)=preconditioned_residual(x,y);
                }
                residual_preconditioned_dot=array_dot(d,residual_vector,preconditioned_residual);
                continue;
            }
            if(iteration==ctrl.cg_max)break;
            apply(search_direction,matrix_direction);
            double direction_matrix_dot=array_dot(d,search_direction,matrix_direction);
            if(!(direction_matrix_dot>0)||!std::isfinite(direction_matrix_dot)||!std::isfinite(residual_preconditioned_dot))return info;
            double cg_step=residual_preconditioned_dot/direction_matrix_dot;
            for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
                mu_offset(x,y)+=cg_step*search_direction(x,y);
                residual_vector(x,y)-=cg_step*matrix_direction(x,y);
                preconditioned_residual(x,y)=residual_vector(x,y)/preconditioner_diagonal(x,y);
            }
            double next_residual_dot=array_dot(d,residual_vector,preconditioned_residual),direction_weight=next_residual_dot/residual_preconditioned_dot;
            residual_preconditioned_dot=next_residual_dot;
            for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x)search_direction(x,y)=preconditioned_residual(x,y)+direction_weight*search_direction(x,y);
        }
        if(!info.ok)return info;
        d.halo(mu_offset);
        double closure=0;
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
            auto w=weights_at(x,y);
            double change=dt*(w[0]*(mu_offset(x-1,y)-mu_offset(x,y))+w[1]*(mu_offset(x+1,y)-mu_offset(x,y))+w[2]*(mu_offset(x,y-1)-mu_offset(x,y))+w[3]*(mu_offset(x,y+1)-mu_offset(x,y)));
            trial.carbon(x,y)=old.carbon(x,y)+change;
            closure=std::max(closure,std::abs(trial.carbon(x,y)-(intercept(x,y)+chi(x,y)*mu_offset(x,y))));
            trial.mu(x,y)=mu_reference+(trial.carbon(x,y)-intercept(x,y))/chi(x,y);
        }
        info.closure=d.maximum(closure);
        info.ok=info.closure<=ctrl.closure_tol;
        return info;
    }
    // ---- 6D. State checks: finite values, phase simplex, carbon and closure ----
    inline bool valid_state(const Domain& d,const State& s,const Thermo& f,bool phase_only) {
        bool good=true;
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
            double sum=0;
            auto p=cell_phi(s.phi,x,y);
            for(double a:p) {
                good=good&&std::isfinite(a)&&a>=-1e-12&&a<=1+1e-12;
                sum+=a;
            }
            good=good&&std::abs(sum-1)<1e-10;
            good=good&&std::isfinite(s.mu(x,y))&&std::isfinite(s.carbon(x,y));
            if(!phase_only) {
                good=good&&s.carbon(x,y)>=-1e-10&&s.carbon(x,y)<=1+1e-10;
                for(int a=0;a<3;++a)if(p[a]>.01) {
                    double c=f[a].c_at_mu(s.mu(x,y));
                    good=good&&c>=-1e-10&&c<=1+1e-10;
                }
            }
        }
        return d.all(good);
    }
    // ---- 6E. Accept/reject the complete coupled step, then propose the next dt ----
    // Output times can clip accepted dt even if the rejection count remains zero.
    inline StepInfo advance(const Domain& d,State& s,const Thermo& f,const PhaseParameters& pp,const Controls& ctrl,double mu_reference,double stop_time) {
        double dt=std::min({s.dt,ctrl.dt_max,stop_time-s.time});
        if(!(dt>0)||!std::isfinite(dt))throw std::runtime_error("invalid advance time");
        const double e0=ctrl.check_energy?reduced_energy(d,s,f,pp,mu_reference,ctrl.phase_only):0;
        const double mass0=ctrl.phase_only?0:total_carbon(d,s);
        std::string why;
        for(int attempt=0;attempt<=ctrl.max_retries;++attempt) {
            State trial=s;
            StepInfo out;
            out.retries=attempt;
            out.accepted_dt=dt;
            out.max_phase_change=propose_phi(d,s,trial,f,pp,dt);
            bool good=std::isfinite(out.max_phase_change)&&out.max_phase_change<=ctrl.max_phase_change;
            why="phase increment";
            if(good) {
                if(ctrl.phase_only) {
                    for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x)trial.carbon(x,y)=composition(trial.mu(x,y),mixture(cell_phi(trial.phi,x,y),f));
                }
                else {
                    auto di=implicit_diffusion(d,s,trial,f,dt,mu_reference,ctrl);
                    good=di.ok;
                    out.cg_iterations=di.it;
                    out.linear_residual=di.residual;
                    out.closure=di.closure;
                    why="linear/closure";
                }
            }
            if(good) {
                good=valid_state(d,trial,f,ctrl.phase_only);
                why="composition/phase bounds";
            }
            if(good&&!ctrl.phase_only) {
                double m=total_carbon(d,trial);
                good=std::abs(m-mass0)<=1e-11*std::max(std::abs(mass0),d.gx*d.gy*d.dx*d.dy*1e-6);
                why="mass conservation";
            }
            if(good&&ctrl.check_energy) {
                double e1=reduced_energy(d,trial,f,pp,mu_reference,ctrl.phase_only);
                out.energy_change=e1-e0;
                double tol=ctrl.energy_atol+ctrl.energy_rtol*std::max(std::abs(e0),pp.interface.sigma*std::max(d.gx*d.dx,d.gy*d.dy));
                good=std::isfinite(e1)&&e1<=e0+tol;
                why="free-energy increase";
            }
            if(good) {
                trial.time=s.time+dt;
                trial.step=s.step+1;
                trial.rejected=s.rejected+attempt;
                trial.dt=std::min(ctrl.dt_max,dt*1.1);
                s=std::move(trial);
                return out;
            }
            dt*=.5;
            if(dt<ctrl.dt_min)throw std::runtime_error("step rejected below dt_min: "+why);
        }
        throw std::runtime_error("step retry limit: "+why);
    }
}
// namespace pearl


// END ORIGINAL FILE: include/solver.hpp

// ============================================================================
// FILLING GEOMETRY — deterministic sharp phase labels in global grid coordinates
// ============================================================================
// The input commands below retain the compatible two-dimensional geometry of
// the older MicroSim filling routines. Lengths are in grid cells; angles are in
// degrees counterclockwise from +X. Their coordinates refer to integer grid
// indices, before Pearl's existing capillarity-only smoothing and carbon setup.
// Each rank evaluates the same ordered commands using global indices: this
// module contains no rank-dependent random generator and no physical evolution.
// Dependency: pearl::Three and phase IDs Alpha=0, Theta=1, Gamma=2 are defined.
#include <algorithm>
#include <cctype>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace pearl {
    struct FillingOperation {
        std::string method;
        std::vector<double> parameters;
    };

    // Parse one complete brace/comma list. Reject partially parsed numbers,
    // NaN/Inf, missing braces, empty items and unconsumed suffixes. The enclosing
    // input reader handles comments; a final semicolon remains optional.
    inline FillingOperation parse_filling_operation(const std::string& method,
                                                     const std::string& value,
                                                     int line) {
        const auto trim = [](const std::string& text) {
            const auto first = text.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) return std::string{};
            return text.substr(first, text.find_last_not_of(" \t\r\n")-first+1);
        };
        const auto malformed = [&]() -> std::invalid_argument {
            return std::invalid_argument("Invalid " + method + " filling parameters at input line "
                + std::to_string(line) + ": expected {finite_number,...}, optionally followed by ';'");
        };
        std::string list = trim(value);
        if (!list.empty() && list.back() == ';') {
            list.pop_back();
            list = trim(list);
        }
        if (list.size() < 3 || list.front() != '{' || list.back() != '}') throw malformed();
        list = list.substr(1, list.size()-2);
        FillingOperation operation{method,{}};
        std::size_t item_start = 0;
        while (true) {
            const std::size_t comma = list.find(',', item_start);
            const std::string item = trim(list.substr(item_start,
                comma == std::string::npos ? std::string::npos : comma-item_start));
            if (item.empty()) throw malformed();
            std::size_t consumed = 0;
            double parameter = 0.0;
            try { parameter = std::stod(item, &consumed); }
            catch (const std::exception&) { throw malformed(); }
            if (consumed != item.size() || !std::isfinite(parameter)) throw malformed();
            operation.parameters.push_back(parameter);
            if (comma == std::string::npos) break;
            item_start = comma+1;
        }
        return operation;
    }

    // Validate once, before any cell is initialized. Half-open boxes use a high
    // bound of nx/ny; inclusive cube bounds use nx-1/ny-1. Arms, bridge disks and
    // cylinder edges may extend beyond the domain and are clipped, as in the
    // original integer-coordinate masks. Centers must lie in the physical grid.
    inline void validate_filling_operations(const std::vector<FillingOperation>& operations,
                                            int nx, int ny) {
        if (nx <= 0 || ny <= 0) throw std::invalid_argument("Filling requires positive nx and ny");
        for (const auto& operation : operations) {
            const auto& p = operation.parameters;
            const auto fail = [&](const std::string& reason) {
                throw std::invalid_argument(operation.method + ": " + reason);
            };
            const auto count = [&](std::size_t expected) {
                if (p.size() != expected) fail("expected " + std::to_string(expected) + " parameters");
                for (double value : p) if (!std::isfinite(value)) fail("parameters must be finite");
            };
            const auto integer = [&](std::size_t index) {
                if (std::trunc(p[index]) != p[index] ||
                    p[index] < static_cast<double>(std::numeric_limits<int>::min()) ||
                    p[index] > static_cast<double>(std::numeric_limits<int>::max()))
                    fail("parameter " + std::to_string(index+1) + " must be a representable integer");
            };
            const auto phase = [&](std::size_t index) {
                integer(index);
                if (p[index] < Alpha || p[index] > Gamma)
                    fail("phase IDs must be alpha=0, theta=1 or gamma=2; gamma grain clones are unsupported");
            };
            const auto distinct_phases = [&](std::initializer_list<std::size_t> indices) {
                for (auto index : indices) phase(index);
                for (auto first = indices.begin(); first != indices.end(); ++first)
                    for (auto second = first+1; second != indices.end(); ++second)
                        if (p[*first] == p[*second]) fail("listed phase IDs must be distinct");
            };
            const auto half_open_box = [&]() {
                for (std::size_t index = 0; index < 4; ++index) integer(index);
                if (!(0 <= p[0] && p[0] < p[1] && p[1] <= nx &&
                      0 <= p[2] && p[2] < p[3] && p[3] <= ny))
                    fail("box must satisfy 0<=xlo<xhi<=nx and 0<=ylo<yhi<=ny");
            };
            const auto center = [&](std::size_t x, std::size_t y, bool integer_coordinates) {
                if (integer_coordinates) { integer(x); integer(y); }
                if (!(0 <= p[x] && p[x] < nx && 0 <= p[y] && p[y] < ny))
                    fail("center/junction must be inside the physical grid");
            };
            const auto positive = [&](std::size_t index) {
                if (!(p[index] > 0.0)) fail("parameter " + std::to_string(index+1) + " must be positive");
            };
            const auto fraction = [&](std::size_t index) {
                if (!(p[index] > 0.0 && p[index] < 1.0)) fail("cementite fraction must be strictly between 0 and 1");
            };
            const auto nonnegative = [&](std::size_t index) {
                if (p[index] < 0.0) fail("parameter " + std::to_string(index+1) + " must be nonnegative");
            };

            const auto angle = [&](std::size_t index) {
                if (!std::isfinite(p[index]*3.14159265358979323846/180.0))
                    fail("angle conversion to radians overflows");
            };
            const auto squared_length = [&](std::size_t index) {
                const double squared = p[index]*p[index];
                if (!std::isfinite(squared) || (p[index] > 0.0 && squared == 0.0))
                    fail("length squared overflows or underflows");
            };
            const auto resolved_positive = [&](double value, const std::string& name) {
                if (!(value > 0.0) || !std::isfinite(value)) fail(name + " must remain finite and positive");
            };

            if (operation.method == "FILLCUBE") {
                count(7); phase(0);
                for (std::size_t index = 1; index < 7; ++index) integer(index);
                if (p[3] != 0 || p[6] != 0) fail("2D filling requires zlo=zhi=0");
                if (!(0 <= p[1] && p[1] <= p[4] && p[4] < nx &&
                      0 <= p[2] && p[2] <= p[5] && p[5] < ny))
                    fail("inclusive bounds must satisfy 0<=xlo<=xhi<nx and 0<=ylo<=yhi<ny");
            } else if (operation.method == "FILLCENTERBOX2D") {
                count(6); half_open_box(); distinct_phases({4,5});
            } else if (operation.method == "FILLYJUNCTION2D" ||
                       operation.method == "FILLYJUNCTIONLAMELLAE2D") {
                const bool lamellar = operation.method == "FILLYJUNCTIONLAMELLAE2D";
                count(lamellar ? 16 : 13); half_open_box(); center(4,5,true);
                positive(9); positive(10); angle(6); angle(7); angle(8);
                if (lamellar) {
                    positive(11); positive(12); distinct_phases({13,14,15});
                    resolved_positive(p[11]+p[12], "lamellar period");
                }
                else distinct_phases({11,12});
            } else if (operation.method == "FILLCYLINDER") {
                count(6); phase(0); center(1,2,true); integer(3); integer(4); positive(5);
                if (p[3] != 0 || p[4] != 0) fail("2D cylinder cross-section requires zlo=zhi=0");
                squared_length(5);
            } else if (operation.method == "FILLADHEREDPEARLITEARM2D") {
                count(11); center(0,1,false); positive(3); integer(4); positive(4);
                if (p[4] > std::numeric_limits<int>::max()/2) fail("pair_count is too large");
                fraction(5); positive(6); nonnegative(7); nonnegative(8); distinct_phases({9,10});
                if (p[7] > p[6]) fail("contact_half_width must not exceed peak_half_width");
                angle(2); squared_length(3); squared_length(6); squared_length(8);
                const double period = p[3]/p[4];
                resolved_positive(period, "arm period");
                resolved_positive((1.0-p[5])*period, "ferrite segment length");
                resolved_positive(p[5]*period, "cementite segment length");
            } else if (operation.method == "FILLPEARLITEELLIPSE2D") {
                count(11); center(0,1,false); nonnegative(3); integer(4); positive(4);
                positive(5); positive(6); fraction(7); nonnegative(8); distinct_phases({9,10});
                angle(2); squared_length(8);
                const double length = p[4]*p[5];
                const double semi_major = 0.5*length, semi_minor = p[6]*semi_major;
                resolved_positive(length, "ellipse length");
                resolved_positive(semi_major, "ellipse semimajor axis");
                resolved_positive(semi_minor, "ellipse semiminor axis");
                resolved_positive(p[7]*p[5], "cementite lamellar width");
                if (!std::isfinite((p[3]+length)*(p[3]+length))) fail("ellipse endpoint distance overflows");
            } else if (operation.method == "FILLTHREEGRAINYJUNCTION2D" ||
                       operation.method == "FILLYJUNCTIONPEARLITE2D") {
                fail("historical multi-grain initialization requires separate gamma grain fields; "
                     "Pearl has one gamma field and will not merge grain IDs silently");
            } else {
                fail("unsupported filling method in the 2D three-phase Pearl solver; see FILLING_METHODS.md");
            }
        }
    }

    // One-hot labels are sharp initial phases. They are deliberately not diffuse
    // profiles: the existing sharp_smooth path smooths this map afterwards.
    inline Three filling_one_phase(int selected_phase) {
        Three phase_fraction{0.0,0.0,0.0};
        phase_fraction[selected_phase] = 1.0;
        return phase_fraction;
    }

    // Dedicated pearlite ellipses put theta in the center of every period.
    // This differs from Y-lamellae, whose periods start with a ferrite segment.
    inline int filling_centered_lamellar_phase(double distance_along_arm, double period,
                                               double cementite_fraction,
                                               int ferrite_phase, int cementite_phase) {
        double period_position = std::fmod(distance_along_arm, period);
        if (period_position < 0.0) period_position += period;
        const double cementite_width = cementite_fraction*period;
        const double cementite_start = 0.5*(period-cementite_width);
        const double cementite_end = 0.5*(period+cementite_width);
        return period_position >= cementite_start && period_position < cementite_end
            ? cementite_phase : ferrite_phase;
    }

    // Apply validated commands in input order. Whole-map recipes replace the
    // current label; overlay recipes replace only cells inside their masks.
    // Initial gamma is the original residual background for the 2D model.
    inline Three filling_phase_at(const std::vector<FillingOperation>& operations,
                                  int global_x, int global_y) {
        Three phase_fraction = filling_one_phase(Gamma);
        const double pi = 3.14159265358979323846;
        for (const auto& operation : operations) {
            const auto& p = operation.parameters;
            if (operation.method == "FILLCUBE") {
                // Legacy last-phase FILLCUBE only recomputes the residual. In
                // this one-hot map that is a no-op, not a gamma eraser rectangle.
                const int selected_phase = static_cast<int>(p[0]);
                if (selected_phase != Gamma && global_x >= p[1] && global_x <= p[4] &&
                    global_y >= p[2] && global_y <= p[5])
                    phase_fraction = filling_one_phase(selected_phase);
            } else if (operation.method == "FILLCYLINDER") {
                const double distance_x = global_x-p[1], distance_y = global_y-p[2];
                const int selected_phase = static_cast<int>(p[0]);
                if (selected_phase != Gamma && distance_x*distance_x+distance_y*distance_y <= p[5]*p[5])
                    phase_fraction = filling_one_phase(selected_phase);
            } else if (operation.method == "FILLCENTERBOX2D") {
                const bool inside = global_x >= p[0] && global_x < p[1] && global_y >= p[2] && global_y < p[3];
                phase_fraction = filling_one_phase(static_cast<int>(inside ? p[4] : p[5]));
            } else if (operation.method == "FILLYJUNCTION2D" ||
                       operation.method == "FILLYJUNCTIONLAMELLAE2D") {
                const bool lamellar = operation.method == "FILLYJUNCTIONLAMELLAE2D";
                const int ferrite_phase = static_cast<int>(p[lamellar ? 13 : 11]);
                const int austenite_phase = static_cast<int>(p[lamellar ? 15 : 12]);
                phase_fraction = filling_one_phase(austenite_phase);
                if (!(global_x >= p[0] && global_x < p[1] && global_y >= p[2] && global_y < p[3])) continue;
                const double distance_x = global_x-p[4], distance_y = global_y-p[5];
                for (int arm = 0; arm < 3; ++arm) {
                    const double angle_radians = p[6+arm]*pi/180.0;
                    const double direction_x = std::cos(angle_radians), direction_y = std::sin(angle_radians);
                    const double along = distance_x*direction_x+distance_y*direction_y;
                    const double across = -distance_x*direction_y+distance_y*direction_x;
                    if (along >= 0.0 && along <= p[9] && std::fabs(across) <= p[10]) {
                        int selected_phase = ferrite_phase;
                        if (lamellar) {
                            const double period = p[11]+p[12];
                            double period_position = std::fmod(along, period);
                            if (period_position < 0.0) period_position += period;
                            if (!(period_position < p[11])) selected_phase = static_cast<int>(p[14]);
                        }
                        phase_fraction = filling_one_phase(selected_phase);
                        break; // Original ordering: the first intersecting arm wins.
                    }
                }
            } else if (operation.method == "FILLADHEREDPEARLITEARM2D") {
                const double angle_radians = p[2]*pi/180.0;
                const double direction_x = std::cos(angle_radians), direction_y = std::sin(angle_radians);
                const double distance_x = global_x-p[0], distance_y = global_y-p[1];
                const double along = distance_x*direction_x+distance_y*direction_y;
                const double across = -distance_x*direction_y+distance_y*direction_x;
                const double arm_length = p[3];
                const int pair_count = static_cast<int>(p[4]);
                const double period = arm_length/static_cast<double>(pair_count);
                const double ferrite_width = (1.0-p[5])*period, cementite_width = p[5]*period;
                const double peak_half_width = p[6], contact_half_width = p[7], bridge_radius = p[8];
                const int ferrite_phase = static_cast<int>(p[9]), cementite_phase = static_cast<int>(p[10]);
                int selected_phase = -1;
                if (along >= 0.0 && along <= arm_length) {
                    int pair_id = along >= arm_length ? pair_count-1 : static_cast<int>(std::floor(along/period));
                    pair_id = std::max(0, std::min(pair_count-1, pair_id));
                    const double pair_start = static_cast<double>(pair_id)*period;
                    const double local_distance = along-pair_start;
                    int segment_id;
                    double segment_start, segment_length;
                    if (local_distance < ferrite_width) {
                        segment_id = 2*pair_id; segment_start = pair_start; segment_length = ferrite_width;
                        selected_phase = ferrite_phase;
                    } else {
                        segment_id = 2*pair_id+1; segment_start = pair_start+ferrite_width; segment_length = cementite_width;
                        selected_phase = cementite_phase;
                    }
                    double fraction_along_segment = (along-segment_start)/segment_length;
                    if (fraction_along_segment < 0.0) fraction_along_segment = 0.0;
                    else if (fraction_along_segment > 1.0) fraction_along_segment = 1.0;
                    const double start_half_height = segment_id == 0 ? 0.0 : contact_half_width;
                    const double end_half_height = segment_id == 2*pair_count-1 ? 0.0 : contact_half_width;
                    const double average_endpoint_height = 0.5*(start_half_height+end_half_height);
                    double half_height = (1.0-fraction_along_segment)*start_half_height+
                        fraction_along_segment*end_half_height+
                        (peak_half_width-average_endpoint_height)*std::sin(pi*fraction_along_segment);
                    if (half_height > peak_half_width) half_height = peak_half_width;
                    if (std::fabs(across) > half_height) selected_phase = -1;
                }
                if (bridge_radius > 0.0) {
                    const double bridge_radius_squared = bridge_radius*bridge_radius;
                    const double start_distance_squared = along*along+across*across;
                    const double end_distance_squared = (along-arm_length)*(along-arm_length)+across*across;
                    if (start_distance_squared <= bridge_radius_squared) selected_phase = ferrite_phase;
                    if (end_distance_squared <= bridge_radius_squared) selected_phase = cementite_phase;
                }
                if (selected_phase >= 0) phase_fraction = filling_one_phase(selected_phase);
            } else if (operation.method == "FILLPEARLITEELLIPSE2D") {
                const double angle_radians = p[2]*pi/180.0;
                const double direction_x = std::cos(angle_radians), direction_y = std::sin(angle_radians);
                const double distance_x = global_x-p[0], distance_y = global_y-p[1];
                const double along = distance_x*direction_x+distance_y*direction_y;
                const double across = -distance_x*direction_y+distance_y*direction_x;
                const double start_distance = p[3], period = p[5], length = p[4]*period;
                const double semi_major = 0.5*length, semi_minor = p[6]*semi_major;
                const double center_distance = start_distance+semi_major;
                const double normalized_along = (along-center_distance)/semi_major;
                const double normalized_across = across/semi_minor;
                const double bridge_radius = p[8];
                const int ferrite_phase = static_cast<int>(p[9]), cementite_phase = static_cast<int>(p[10]);
                int selected_phase = -1;
                if (normalized_along*normalized_along+normalized_across*normalized_across <= 1.0)
                    selected_phase = filling_centered_lamellar_phase(along-start_distance, period, p[7], ferrite_phase, cementite_phase);
                if (bridge_radius > 0.0) {
                    const double start_distance_squared = (along-start_distance)*(along-start_distance)+across*across;
                    const double end_distance_squared = (along-(start_distance+length))*(along-(start_distance+length))+across*across;
                    if (start_distance_squared <= bridge_radius*bridge_radius || end_distance_squared <= bridge_radius*bridge_radius)
                        selected_phase = ferrite_phase;
                }
                if (selected_phase >= 0) phase_fraction = filling_one_phase(selected_phase);
            } else {
                throw std::invalid_argument("Unvalidated/unsupported filling method: " + operation.method);
            }
        }
        return phase_fraction;
    }
} // namespace pearl


// BEGIN ORIGINAL FILE: include/io.hpp

// ============================================================================
// 7. INPUT PARAMETERS AND FILLING — parse, validate, then initialize
// ============================================================================
// Read physical/model/numerical controls separately from the geometry operations.
// Temperature preparation occurs once. Geometry precedes carbon initialization.
// Strict key=value parameters. Repeated FILL commands are the sole ordered-list
// exception; unknown commands and duplicate ordinary parameters remain errors.
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <iomanip>
#include <functional>
#ifndef PEARL_SOURCE_DIGEST
#define PEARL_SOURCE_DIGEST "unrecorded-direct-build"
#endif
namespace pearl {
    inline std::string trim(std::string s) {
        auto a=s.find_first_not_of(" \t\r\n");
        if(a==std::string::npos)return {};
        auto b=s.find_last_not_of(" \t\r\n");
        return s.substr(a,b-a+1);
    }
    // ---- 7A. Parsed input: physical parameters, numerical controls, filling ----
    struct Config {
        int nx=0,ny=0,seed_theta_columns=0,phase_a=0,phase_b=2;
        double dx=0,dy=0,end_time=0,output_dt=0,seed_height=0,radius=0;
        double initial_mu=0,initial_carbon=0;
        std::string composition_mode="",initialization="lamella",mode="coupled";
        int smooth_steps=100,smooth_min_steps=100;
        double smooth_energy_rtol=1e-10,smooth_max_phase_change=0.02;
        // thermo is the evaluated, fixed run-temperature thermodynamics.
        Thermo thermo{}, reference_thermo{};
        std::array<GrandPotentialPolynomial,3> temperature_slope{};
        std::string thermo_mode="direct_at_temperature";
        std::string temperature_law="linear_gp";
        PhaseDiagramInput phase_diagram;
        double T=0,T_eutectoid=0,undercooling=0,max_undercooling=0;
        PhaseParameters phase{{0,0,0},{0,0,0},0};
        Controls control;
        double dt=0,shift_trigger=0,shift_target=0,shift_gamma_carbon=0;
        long long max_steps=100000000;
        bool moving_window=false;
        // Geometry operations are opt-in; ordinary legacy input keys remain unique.
        std::vector<FillingOperation> filling_operations;
        std::map<std::string,std::string> raw;
        std::string physics_signature()const{
        std::ostringstream o;
        o<<std::setprecision(17)<<"pearlite-v0.4-reference;source="<<PEARL_SOURCE_DIGEST<<";3ph;binary;equal-sigma;Jat=0;Xperiodic;Ynoflux;";
        o<<nx<<';'<<ny<<';'<<dx<<';'<<dy<<';'<<phase.Vm<<';'<<phase.interface.epsilon<<';'<<phase.interface.sigma<<';'<<phase.interface.triple<<';';
        o<<thermo_mode<<';'<<T<<';'<<T_eutectoid<<';'<<undercooling<<';'<<max_undercooling<<';';
        if(thermo_mode=="eutectoid_reference") {
            for(auto p:reference_thermo)o<<p.free_energy_quadratic<<';'<<p.free_energy_linear<<';'<<p.free_energy_constant<<';';
            if(temperature_law=="linear_gp") {
                // Keep all legacy checkpoint signatures exactly as before.
                for(auto slope:temperature_slope)o<<slope.mu_squared_coefficient<<';'<<slope.mu_coefficient<<';'<<slope.constant_energy<<';';
            } else {
                o<<"phase-diagram-fixed-curvature-v1;";
                for(double carbon:phase_diagram.reference_composition)o<<carbon<<';';
                for(double slope:phase_diagram.composition_slope)o<<slope<<';';
            }
        }
        for(auto p:thermo)o<<p.free_energy_quadratic<<';'<<p.free_energy_linear<<';'<<p.free_energy_constant<<';'<<p.diffusivity<<';';
        for(double t:phase.tau)o<<t<<';';
        o<<mode<<';'<<initialization<<';'<<seed_height<<';'<<seed_theta_columns<<';'<<radius<<';'<<phase_a<<';'<<phase_b<<';'<<composition_mode<<';'<<initial_mu<<';'<<initial_carbon<<';';
        o<<smooth_steps<<';'<<smooth_min_steps<<';'<<smooth_energy_rtol<<';'<<smooth_max_phase_change<<';';
        o<<moving_window<<';'<<shift_trigger<<';'<<shift_target<<';'<<shift_gamma_carbon;
        // Preserve the exact old signature for every existing initialization.
        // New geometry records method order and all coordinates in its signature.
        if(initialization=="microsim") {
            o<<";filling-v1";
            for(const auto& operation:filling_operations) {
                o<<';'<<operation.method;
                for(double coordinate:operation.parameters)o<<';'<<coordinate;
            }
        }
        return o.str();
        }
    };
    // ---- 7B. Strict parameter parsing: old key names remain valid ----
    inline Config read_config(const std::string& path) {
        std::ifstream in(path);
        if(!in)throw std::runtime_error("cannot open input: "+path);
        Config config;
        std::string line;
        int ln=0;
        while(std::getline(in,line)) {
            ++ln;
            line=trim(line.substr(0,line.find('#')));
            if(line.empty())continue;
            auto pos=line.find('=');
            if(pos==std::string::npos)throw std::runtime_error("missing = on line "+std::to_string(ln));
            auto key=trim(line.substr(0,pos)),val=trim(line.substr(pos+1));
            // Repeated FILL commands form an ordered geometry recipe. They are
            // separate from the single-valued physical and numerical parameters.
            if(key.rfind("FILL",0)==0) {
                config.filling_operations.push_back(parse_filling_operation(key,val,ln));
                continue;
            }
            if(key.empty()||val.empty()||config.raw.count(key))throw std::runtime_error("empty/duplicate key: "+key);
            config.raw[key]=val;
        }
        std::set<std::string> used;
        auto str=[&](std::string k,std::string def="") {
            auto it=config.raw.find(k);
            if(it==config.raw.end())return def;
            used.insert(k);
            return it->second;
        };
        auto number=[&](std::string k,double def=std::numeric_limits<double>::quiet_NaN()) {
            auto s=str(k);
            if(s.empty()) {
                if(!std::isfinite(def))throw std::runtime_error("missing required key: "+k);
                return def;
            }
            std::size_t n=0;
            double v=std::stod(s,&n);
            if(n!=s.size()||!std::isfinite(v))throw std::runtime_error("invalid number: "+k);
            return v;
        };
        auto integer=[&](std::string k,int def=-1) {
            double a=number(k,def<0?std::numeric_limits<double>::quiet_NaN():def);
            if(a!=std::floor(a)||std::abs(a)>1e8)throw std::runtime_error("invalid integer: "+k);
            return static_cast<int>(a);
        };
        config.nx=integer("nx");
        config.ny=integer("ny");
        config.dx=number("dx");
        config.dy=number("dy",config.dx);
        config.phase.interface={number("epsilon"),number("sigma"),number("triple")};
        config.phase.Vm=number("Vm");
        config.phase.tau={number("tau_at"),number("tau_ag"),number("tau_tg")};
        std::array<std::string,3> names{"alpha","theta","gamma"};
        config.thermo_mode=str("thermo_mode");
        if(config.thermo_mode=="eutectoid_reference") {
            config.T_eutectoid=number("T_eutectoid");
            config.undercooling=number("undercooling");
            config.max_undercooling=number("max_undercooling");
            config.temperature_law=str("temperature_law");
            if(config.temperature_law!="linear_gp"&&config.temperature_law!="phase_diagram_fixed_curvature")
                throw std::runtime_error("temperature_law must be linear_gp or phase_diagram_fixed_curvature");
            config.T=config.T_eutectoid-config.undercooling;
            if(config.T_eutectoid<=0||config.T<=0||config.undercooling<0||
               config.max_undercooling<0||config.undercooling>config.max_undercooling)
                throw std::runtime_error("temperature outside declared undercooling interval");
            for(int i=0;i<3;++i) {
                const auto& n=names[i];
                config.reference_thermo[i]={number("Aeq_"+n),number("Beq_"+n),
                                       number("Ceq_"+n),number("D_"+n)};
                if(config.temperature_law=="linear_gp") {
                    config.temperature_slope[i]={number("dq2_dT_"+n),number("dq1_dT_"+n),number("dq0_dT_"+n)};
                    config.thermo[i]=at_undercooling(config.reference_thermo[i],config.temperature_slope[i],config.undercooling);
                }
            }
            if(config.temperature_law=="phase_diagram_fixed_curvature") {
                // ABC belongs to T_eutectoid. Boundary slopes supply the separate
                // undercooling correction; do not enter ABC already shifted to T.
                config.phase_diagram.reference_composition={number("ceq_alpha"),number("ceq_theta"),number("ceq_gamma")};
                config.phase_diagram.composition_slope={
                    number("slope_alpha_on_alpha_gamma"),number("slope_gamma_on_alpha_gamma"),
                    number("slope_theta_on_theta_gamma"),number("slope_gamma_on_theta_gamma")};
                config.thermo=at_phase_diagram_undercooling(config.reference_thermo,config.phase_diagram,config.undercooling);
            }
        }
        else if(config.thermo_mode=="direct_at_temperature") {
            // Explicit compatibility mode: A/B/C already belong to temperature.
            // No second undercooling correction is applied in this mode.
            config.T=number("temperature");
            if(config.T<=0)throw std::runtime_error("temperature must be positive");
            for(int i=0;i<3;++i) {
                const auto& n=names[i];
                config.thermo[i]={number("A_"+n),number("B_"+n),number("C_"+n),number("D_"+n)};
                config.thermo[i].validate();
            }
        }
        else throw std::runtime_error(
            "Set thermo_mode=eutectoid_reference or direct_at_temperature. "
            "Legacy inputs require an explicit mode and temperature; no temperature is guessed.");
        config.dt=number("dt");
        config.control.dt_max=number("dt_max",config.dt);
        config.control.dt_min=number("dt_min",config.dt/1048576.0);
        config.end_time=number("end_time");
        config.output_dt=number("output_dt");
        config.max_steps=integer("max_steps",100000000);
        config.initialization=str("initialization","lamella");
        config.mode=str("mode","coupled");
        if(config.mode!="coupled"&&config.mode!="phase_only")throw std::runtime_error("mode must be coupled or phase_only");
        config.control.phase_only=config.mode=="phase_only";
        config.seed_height=number("seed_height",config.ny*config.dy/4);
        config.seed_theta_columns=integer("theta_columns",config.nx/2);
        config.radius=number("radius",std::min(config.nx*config.dx,config.ny*config.dy)/4);
        config.phase_a=integer("phase_a",0);
        config.phase_b=integer("phase_b",2);
        if(config.initialization!="lamella"&&config.initialization!="sharp_smooth"&&config.initialization!="flat"&&config.initialization!="circle"&&config.initialization!="microsim")throw std::runtime_error("unknown initialization");
        if(config.initialization=="microsim") {
            if(config.filling_operations.empty())throw std::runtime_error("initialization=microsim requires at least one FILL operation");
            validate_filling_operations(config.filling_operations,config.nx,config.ny);
        } else if(!config.filling_operations.empty()) {
            throw std::runtime_error("FILL operations require initialization=microsim; no geometry is silently ignored");
        }
        config.smooth_steps=integer("smooth_steps",100);
        config.smooth_min_steps=integer("smooth_min_steps",config.smooth_steps);
        config.smooth_energy_rtol=number("smooth_energy_rtol",1e-10);
        config.smooth_max_phase_change=number("smooth_max_phase_change",0.02);
        if(config.phase_a<0||config.phase_a>2||config.phase_b<0||config.phase_b>2||config.phase_a==config.phase_b)throw std::runtime_error("invalid phase pair");
        int choices=config.raw.count("initial_mu")+config.raw.count("gamma_carbon")+config.raw.count("mean_carbon");
        if(choices!=1)throw std::runtime_error("specify exactly one of initial_mu, gamma_carbon, mean_carbon");
        if(config.raw.count("initial_mu")) {
            config.composition_mode="mu";
            config.initial_mu=number("initial_mu");
        }
        if(config.raw.count("gamma_carbon")) {
            config.composition_mode="gamma";
            config.initial_carbon=number("gamma_carbon");
            config.initial_mu=config.thermo[2].mu(config.initial_carbon);
        }
        if(config.raw.count("mean_carbon")) {
            config.composition_mode="mean";
            config.initial_carbon=number("mean_carbon");
        }
        config.control.max_phase_change=number("max_phase_change",.04);
        config.control.cg_rtol=number("cg_rtol",1e-11);
        config.control.cg_atol=number("cg_atol",1e-14);
        config.control.cg_max=integer("cg_max",2000);
        config.control.closure_tol=number("closure_tol",1e-10);
        config.control.energy_atol=number("energy_atol",1e-13*config.phase.interface.sigma*std::max(config.nx*config.dx,config.ny*config.dy));
        config.moving_window=integer("moving_window",0)!=0;
        config.shift_trigger=number("shift_trigger",.65*config.ny*config.dy);
        config.shift_target=number("shift_target",.35*config.ny*config.dy);
        config.shift_gamma_carbon=number("shift_gamma_carbon",config.composition_mode=="gamma"?config.initial_carbon:config.thermo[2].c_at_mu(config.initial_mu));
        for(auto& kv:config.raw)if(!used.count(kv.first))throw std::runtime_error("unknown key: "+kv.first);
        config.phase.validate();
        if(config.nx<3||config.ny<2||config.dx<=0||config.dy<=0||config.dt<=0||config.control.dt_min<=0||config.control.dt_max<config.dt||config.control.dt_min>config.dt||config.end_time<=0||config.output_dt<=0||config.max_steps<1)throw std::runtime_error("invalid geometry/time bounds");
        if(config.seed_theta_columns<1||config.seed_theta_columns>=config.nx||config.seed_height<=0||config.seed_height>=config.ny*config.dy)throw std::runtime_error("invalid seed geometry");
        if(config.control.max_phase_change<=0||config.control.max_phase_change>.2||config.control.cg_rtol<=0||config.control.cg_atol<=0||config.control.closure_tol<=0||config.control.energy_atol<0)throw std::runtime_error("invalid numerical control");
        if(config.smooth_steps<1||config.smooth_min_steps<0||config.smooth_min_steps>config.smooth_steps||config.smooth_energy_rtol<=0||config.smooth_max_phase_change<=0||config.smooth_max_phase_change>.2)throw std::runtime_error("invalid smoothing control");
        if(config.moving_window&&(config.shift_target<=0||config.shift_trigger<=config.shift_target||config.shift_trigger>=config.ny*config.dy||config.shift_gamma_carbon<0||config.shift_gamma_carbon>1))throw std::runtime_error("invalid moving window settings");
        return config;
    }
    using InitializationObserver=std::function<void(const Domain&,const State&,const char*,int)>;
    // ---- 7D. Filling -> optional smoothing -> common-mu carbon initialization ----
    inline double initialize(const Domain& d,State& s,const Config& config,SmoothReport* smooth_report=nullptr,
                             const InitializationObserver& observer={}) {
        double width=config.seed_theta_columns*d.dx,L=d.gx*d.dx,eps=config.phase.interface.epsilon;
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
            double xx=(d.x0+x-.5)*d.dx,yy=(d.y0+y-.5)*d.dy;
            Three p{};
            if(config.initialization=="lamella") {
                double dist=xx<width ? -std::min(xx,width-xx) : std::min(xx-width,L-xx);
                double th=profile(dist,eps),product=profile(yy-config.seed_height,eps);
                p={(1-th)*product,th*product,1-product};
            }
            else if(config.initialization=="microsim") {
                // Integer global indices reproduce the historical masks.
                // MPI halos are excluded; every layout sees identical labels.
                p=filling_phase_at(config.filling_operations,d.x0+x-1,d.y0+y-1);
            }
            else if(config.initialization=="sharp_smooth") {
                // Exactly sharp phase labels first. X is periodic: theta occupies
                // the first theta_columns; below seed_height the rest is alpha;
                // above seed_height everything is gamma.
                const int global_x=d.x0+x-1;
                const bool product=yy<config.seed_height;
                if(!product)p={0,0,1};
                else if(global_x<config.seed_theta_columns)p={0,1,0};
                else p={1,0,0};
            }
            else {
                double dist=config.initialization=="flat"?yy-config.seed_height:std::hypot(xx-.5*L,yy-.5*d.gy*d.dy)-config.radius;
                double u=profile(dist,eps);
                p[config.phase_a]=u;
                p[config.phase_b]=1-u;
            }
            set_phi(s.phi,x,y,p);
        }
        d.halo(s.phi);
        if(config.initialization=="sharp_smooth"||config.initialization=="microsim") {
            if(observer)observer(d,s,"init_sharp",0);
            auto report=smooth_interface_only(d,s,config.phase,config.smooth_steps,config.smooth_min_steps,
                                              config.smooth_energy_rtol,config.smooth_max_phase_change);
            if(smooth_report)*smooth_report=report;
            if(observer)observer(d,s,"init_smoothed",report.steps);
        }
        else if(smooth_report) {
            smooth_report->phase_fraction_before=phase_fractions(d,s);
            smooth_report->phase_fraction_after=smooth_report->phase_fraction_before;
            smooth_report->energy_initial=interface_energy_only(d,s,config.phase);
            smooth_report->energy_final=smooth_report->energy_initial;
            smooth_report->converged=true;
        }
        // Carbon and chemical potential are initialized ONLY after preprocessing.
        // Thus the sharp->smooth stage contains no carbon diffusion and consumes
        // no physical simulation time.
        double mu=config.initial_mu;
        if(config.composition_mode=="mean") {
            double si=0,sc=0;
            for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
                auto m=mixture(cell_phi(s.phi,x,y),config.thermo);
                si+=m.c_intercept;
                sc+=m.susceptibility;
            }
            mu=(d.gx*d.gy*config.initial_carbon-d.sum(si))/d.sum(sc);
        }
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
            s.mu(x,y)=mu;
            s.carbon(x,y)=composition(mu,mixture(cell_phi(s.phi,x,y),config.thermo));
        }
        s.time=0;
        s.step=0;
        s.rejected=0;
        s.dt=config.dt;
        s.initial_carbon=total_carbon(d,s);
        d.halo(s.phi);
        d.halo(s.carbon);
        d.halo(s.mu);
        if(!valid_state(d,s,config.thermo,config.control.phase_only))throw std::runtime_error("initial state out of model range");
        return mu;
    }

}
// namespace pearl


// ============================================================================
// 8. OUTPUT AND RESTART — HDF5 fields, provenance, and moving window
// ============================================================================
// Rank zero saves five physical fields without halos. Restart preserves
// physical time, carbon bookkeeping, and the physics/filling signature.
// Rank-0 serial HDF5 I/O is deliberate: works with ASU's serial HDF5 module.
// Five physical datasets, no ghosts; checkpoints are repartitionable.
#include <hdf5.h>
#include <filesystem>
#include <functional>
namespace pearl {
    inline void hcheck(herr_t s,const char* name) {
        if(s<0)throw std::runtime_error(std::string("HDF5: ")+name);
    }
    struct HObject {
        hid_t id=-1;
        herr_t (*close)(hid_t)=nullptr;
        HObject(hid_t a,herr_t(*b)(hid_t)):id(a),close(b) {
            if(a<0)throw std::runtime_error("HDF5 open/create failure");
        }
        ~HObject() {
            if(id>=0&&close)close(id);
        }
        HObject(const HObject&)=delete;
        operator hid_t()const{return id;
        }
    };
    inline void attr_double(hid_t file,const char* name,double value) {
        HObject space(H5Screate(H5S_SCALAR),H5Sclose);
        HObject a(H5Acreate2(file,name,H5T_NATIVE_DOUBLE,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
        hcheck(H5Awrite(a,H5T_NATIVE_DOUBLE,&value),"write scalar");
    }
    inline double read_double(hid_t file,const char* name) {
        HObject a(H5Aopen(file,name,H5P_DEFAULT),H5Aclose);
        double x;
        hcheck(H5Aread(a,H5T_NATIVE_DOUBLE,&x),"read scalar");
        return x;
    }
    inline void attr_string(hid_t file,const char* name,const std::string& value) {
        HObject type(H5Tcopy(H5T_C_S1),H5Tclose);
        hcheck(H5Tset_size(type,value.size()+1),"string size");
        HObject sp(H5Screate(H5S_SCALAR),H5Sclose);
        HObject a(H5Acreate2(file,name,type,sp,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
        hcheck(H5Awrite(a,type,value.c_str()),"string");
    }
    inline std::string read_string(hid_t file,const char* name) {
        HObject a(H5Aopen(file,name,H5P_DEFAULT),H5Aclose);
        HObject t(H5Aget_type(a),H5Tclose);
        std::vector<char> s(H5Tget_size(t)+1,0);
        hcheck(H5Aread(a,t,s.data()),"read string");
        return s.data();
    }
    inline void write_dataset(hid_t file,const std::string& name,const std::vector<double>& v,int nx,int ny) {
        hsize_t dims[2]{static_cast<hsize_t>(ny),static_cast<hsize_t>(nx)},chunk[2]{static_cast<hsize_t>(std::min(ny,64)),static_cast<hsize_t>(std::min(nx,64))};
        HObject sp(H5Screate_simple(2,dims,nullptr),H5Sclose),prop(H5Pcreate(H5P_DATASET_CREATE),H5Pclose);
        hcheck(H5Pset_chunk(prop,2,chunk),"chunk");
        hcheck(H5Pset_fletcher32(prop),"checksum");
        HObject set(H5Dcreate2(file,name.c_str(),H5T_IEEE_F64LE,sp,H5P_DEFAULT,prop,H5P_DEFAULT),H5Dclose);
        hcheck(H5Dwrite(set,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,v.data()),"dataset write");
    }
    inline std::vector<double> read_dataset(hid_t file,const char* name,int nx,int ny) {
        HObject set(H5Dopen2(file,name,H5P_DEFAULT),H5Dclose),sp(H5Dget_space(set),H5Sclose);
        hsize_t dims[2]{};
        if(H5Sget_simple_extent_ndims(sp)!=2)throw std::runtime_error("dataset is not 2D");
        H5Sget_simple_extent_dims(sp,dims,nullptr);
        if(dims[0]!=static_cast<hsize_t>(ny)||dims[1]!=static_cast<hsize_t>(nx))throw std::runtime_error("checkpoint dimensions mismatch");
        std::vector<double> v(static_cast<std::size_t>(nx)*ny);
        hcheck(H5Dread(set,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,v.data()),"dataset read");
        return v;
    }
    // ---- 8A. Collective gather; rank-zero field and metadata output ----
    inline void write_frame(const Domain& d,const State& s,const Config& config,double mu_reference,const std::filesystem::path& path) {
        auto phi=d.gather(s.phi),carbon=d.gather(s.carbon),mu=d.gather(s.mu);
        if(d.rank==0) {
            if(std::filesystem::exists(path))throw std::runtime_error("refuse overwriting frame "+path.string());
            auto temp=path;
            temp+=".part";
            {
            HObject f(H5Fcreate(temp.c_str(),H5F_ACC_EXCL,H5P_DEFAULT,H5P_DEFAULT),H5Fclose);
            std::array<std::string,3> names{"alpha","theta","gamma"};
            for(int a=0;a<3;++a) {
                std::vector<double> v(carbon.size());
                for(std::size_t i=0;i<v.size();++i)v[i]=phi[3*i+a];
                write_dataset(f,names[a],v,d.gx,d.gy);
            }
            write_dataset(f,"carbon",carbon,d.gx,d.gy);
            write_dataset(f,"mu",mu,d.gx,d.gy);
            attr_string(f,"physics_signature",config.physics_signature());
            attr_string(f,"format","pearlite-v0.4-5fields");
            attr_string(f,"thermo_mode",config.thermo_mode);
            attr_double(f,"temperature_K",config.T);
            if(config.thermo_mode=="eutectoid_reference") {
                attr_double(f,"T_eutectoid_K",config.T_eutectoid);
                attr_double(f,"undercooling_K",config.undercooling);
            }
            attr_double(f,"time",s.time);
            attr_double(f,"dt_next",s.dt);
            attr_double(f,"step",static_cast<double>(s.step));
            attr_double(f,"dx",d.dx);
            attr_double(f,"dy",d.dy);
            attr_double(f,"nx",d.gx);
            attr_double(f,"ny",d.gy);
            attr_double(f,"shift_distance",s.shift_distance);
            attr_double(f,"shifts",static_cast<double>(s.shifts));
            attr_double(f,"initial_carbon",s.initial_carbon);
            attr_double(f,"carbon_added",s.carbon_added);
            attr_double(f,"carbon_removed",s.carbon_removed);
            attr_double(f,"mu_reference",mu_reference);
            attr_double(f,"rejected_steps",static_cast<double>(s.rejected));
            hcheck(H5Fflush(f,H5F_SCOPE_GLOBAL),"flush");
            }std::filesystem::rename(temp,path);
            auto xmf=path;
            xmf.replace_extension(".xmf");
            auto xt=xmf;
            xt+=".part";
            std::ofstream out(xt);
            if(!out)throw std::runtime_error("cannot write XMF");
            out<<std::setprecision(17)<<"<?xml version=\"1.0\"?>\n<Xdmf Version=\"3.0\"><Domain><Grid Name=\"pearlite\" GridType=\"Uniform\">\n<Time Value=\""<<s.time<<"\"/>\n<Topology TopologyType=\"2DCoRectMesh\" Dimensions=\""<<d.gy<<" "<<d.gx<<"\"/>\n<Geometry GeometryType=\"ORIGIN_DXDY\"><DataItem Dimensions=\"2\" Format=\"XML\">"<<.5*d.dy+s.shift_distance<<" "<<.5*d.dx<<"</DataItem><DataItem Dimensions=\"2\" Format=\"XML\">"<<d.dy<<" "<<d.dx<<"</DataItem></Geometry>\n";
            for(auto name:{"alpha","theta","gamma","carbon","mu"})out<<"<Attribute Name=\""<<name<<"\" AttributeType=\"Scalar\" Center=\"Node\"><DataItem Dimensions=\""<<d.gy<<" "<<d.gx<<"\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">"<<path.filename().string()<<":/"<<name<<"</DataItem></Attribute>\n";
            out<<"</Grid></Domain></Xdmf>\n";
            out.close();
            if(!out)throw std::runtime_error("XMF write failure");
            std::filesystem::rename(xt,xmf);
        }
        d.barrier();
    }
    // ---- 8B. Restart: validate the complete signature before loading fields ----
    inline double read_checkpoint(const Domain& d,State& s,const Config& config,const std::string& path) {
        std::vector<double> phi,carbon,mu;
        std::array<double,10> meta{};
        bool okay=true;
        std::string error;
        if(d.rank==0) {
            try {
                HObject f(H5Fopen(path.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose);
                if(read_string(f,"format")!="pearlite-v0.4-5fields"||read_string(f,"physics_signature")!=config.physics_signature())throw std::runtime_error("checkpoint physics signature mismatch");
                auto a=read_dataset(f,"alpha",d.gx,d.gy),t=read_dataset(f,"theta",d.gx,d.gy),g=read_dataset(f,"gamma",d.gx,d.gy);
                carbon=read_dataset(f,"carbon",d.gx,d.gy);
                mu=read_dataset(f,"mu",d.gx,d.gy);
                phi.resize(3*carbon.size());
                for(std::size_t i=0;i<carbon.size();++i) {
                    phi[3*i]=a[i];
                    phi[3*i+1]=t[i];
                    phi[3*i+2]=g[i];
                }
                meta={read_double(f,"time"),read_double(f,"dt_next"),read_double(f,"step"),read_double(f,"shift_distance"),read_double(f,"shifts"),read_double(f,"initial_carbon"),read_double(f,"carbon_added"),read_double(f,"carbon_removed"),read_double(f,"mu_reference"),read_double(f,"rejected_steps")};
                for(double v:meta)if(!std::isfinite(v))throw std::runtime_error("nonfinite checkpoint metadata");
                if(meta[0]<0||meta[1]<=0||meta[2]<0||meta[2]!=std::floor(meta[2])||meta[2]>9e15)throw std::runtime_error("bad checkpoint time/step");
            }
            catch(const std::exception& e) {
                okay=false;
                error=e.what();
            }
        }
        if(!d.all(okay))throw std::runtime_error(d.rank==0?error:"rank-0 checkpoint read failed");
#ifdef PEARL_MPI
        MPI_Bcast(meta.data(),static_cast<int>(meta.size()),MPI_DOUBLE,0,d.comm);
#endif
        d.load_global(s.phi,std::move(phi));
        d.load_global(s.carbon,std::move(carbon));
        d.load_global(s.mu,std::move(mu));
        s.time=meta[0];
        s.dt=std::min(config.control.dt_max,meta[1]);
        s.step=static_cast<long long>(meta[2]);
        s.shift_distance=meta[3];
        s.shifts=static_cast<long long>(meta[4]);
        s.initial_carbon=meta[5];
        s.carbon_added=meta[6];
        s.carbon_removed=meta[7];
        s.rejected=static_cast<long long>(meta[9]);
        if(!valid_state(d,s,config.thermo,config.control.phase_only))throw std::runtime_error("invalid checkpoint state");
        double closure=0;
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x)closure=std::max(closure,std::abs(s.carbon(x,y)-composition(s.mu(x,y),mixture(cell_phi(s.phi,x,y),config.thermo))));
        if(d.maximum(closure)>config.control.closure_tol)throw std::runtime_error("checkpoint closure failed");
        double ledger=total_carbon(d,s)+s.carbon_removed-s.carbon_added-s.initial_carbon;
        if(!config.control.phase_only&&std::abs(ledger)>1e-10*d.gx*d.gy*d.dx*d.dy)throw std::runtime_error("checkpoint mass ledger failed");
        return meta[8];
    }
    // Positive shift removes bottom rows, inserts fresh gamma at top. Carbon ledger explicit.
    // ---- 8C. Optional moving window and conservative carbon bookkeeping ----
    inline void shift_rows(const Domain& d,State& s,const Thermo& f,int rows,double c_gamma) {
        if(rows<=0||rows>=d.gy||c_gamma<0||c_gamma>1)throw std::invalid_argument("invalid shift");
        auto phi=d.gather(s.phi),carbon=d.gather(s.carbon);
        double removed=0,added=rows*d.gx*c_gamma*d.dx*d.dy;
        if(d.rank==0) {
            for(int y=0;y<rows;++y)for(int x=0;x<d.gx;++x)removed+=carbon[y*d.gx+x]*d.dx*d.dy;
            for(int y=0;y<d.gy;++y)for(int x=0;x<d.gx;++x) {
                std::size_t i=static_cast<std::size_t>(y)*d.gx+x;
                if(y+rows<d.gy) {
                    std::size_t j=static_cast<std::size_t>(y+rows)*d.gx+x;
                    carbon[i]=carbon[j];
                    for(int a=0;a<3;++a)phi[3*i+a]=phi[3*j+a];
                }
                else {
                    carbon[i]=c_gamma;
                    phi[3*i]=0;
                    phi[3*i+1]=0;
                    phi[3*i+2]=1;
                }
            }
        }
        removed=d.sum(removed);
        d.load_global(s.phi,std::move(phi));
        d.load_global(s.carbon,std::move(carbon));
        close_mu(d,s,f);
        s.shift_distance+=rows*d.dy;
        s.shifts++;
        s.carbon_added+=added;
        s.carbon_removed+=removed;
    }
    inline bool maybe_shift(const Domain& d,State& s,const Config& config) {
        if(!config.moving_window)return false;
        d.halo(s.phi);
        double local=-1;
        for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
            if(y==d.ny&&d.north_edge())continue;
            double a=1-2*s.phi(x,y,2),b=1-2*s.phi(x,y+1,2);
            if(a>=0&&b<0)local=std::max(local,(d.y0+y-.5+a/(a-b))*d.dy);
        }
        double front=d.maximum(local);
        if(front<config.shift_trigger)return false;
        int rows=static_cast<int>(std::floor((front-config.shift_target)/d.dy));
        // Do not erase mother phase. This does not by itself prove diffusion-boundary independence.
        double gamma=0;
        for(int y=1;y<=d.ny;++y)if(d.y0+y-1<rows)for(int x=1;x<=d.nx;++x)gamma=std::max(gamma,s.phi(x,y,2));
        if(d.maximum(gamma)>1e-6)throw std::runtime_error("unsafe shift: removed strip contains gamma");
        shift_rows(d,s,config.thermo,rows,config.shift_gamma_carbon);
        return true;
    }

    // This is not a checkpoint: preprocessing phase fields have no carbon/mu yet.
    // All ranks call this function; rank-0 failure is propagated collectively.
    inline void write_initial_phases(const Domain& d,const State& s,
                                     const std::filesystem::path& path,
                                     const std::string& stage,int pseudo_step) {
        auto phi=d.gather(s.phi);
        bool okay=true;
        std::string error;
        if(d.rank==0) {
            try {
                if(std::filesystem::exists(path))throw std::runtime_error("initial snapshot exists");
                auto temp=path; temp+=".part";
                {
                    HObject f(H5Fcreate(temp.c_str(),H5F_ACC_EXCL,H5P_DEFAULT,H5P_DEFAULT),H5Fclose);
                    std::array<std::string,3> names{"alpha","theta","gamma"};
                    for(int a=0;a<3;++a) {
                        std::vector<double> v(static_cast<std::size_t>(d.gx)*d.gy);
                        for(std::size_t i=0;i<v.size();++i)v[i]=phi[3*i+a];
                        write_dataset(f,names[a],v,d.gx,d.gy);
                    }
                    attr_string(f,"format","pearlite-preprocessing-phases-only");
                    attr_string(f,"stage",stage);
                    attr_double(f,"pseudo_step",pseudo_step);
                    attr_double(f,"dx",d.dx); attr_double(f,"dy",d.dy);
                    hcheck(H5Fflush(f,H5F_SCOPE_GLOBAL),"flush initial snapshot");
                }
                std::filesystem::rename(temp,path);
                auto xmf=path; xmf.replace_extension(".xmf");
                std::ofstream out(xmf);
                out<<std::setprecision(17)<<"<?xml version=\"1.0\"?>\n<Xdmf Version=\"3.0\"><Domain><Grid Name=\""<<stage<<"\" GridType=\"Uniform\">\n"
                   <<"<Topology TopologyType=\"2DCoRectMesh\" Dimensions=\""<<d.gy<<" "<<d.gx<<"\"/>\n"
                   <<"<Geometry GeometryType=\"ORIGIN_DXDY\"><DataItem Dimensions=\"2\" Format=\"XML\">"<<.5*d.dy<<" "<<.5*d.dx
                   <<"</DataItem><DataItem Dimensions=\"2\" Format=\"XML\">"<<d.dy<<" "<<d.dx<<"</DataItem></Geometry>\n";
                for(const char* name:{"alpha","theta","gamma"})
                    out<<"<Attribute Name=\""<<name<<"\" AttributeType=\"Scalar\" Center=\"Node\"><DataItem Dimensions=\""<<d.gy<<" "<<d.gx
                       <<"\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">"<<path.filename().string()<<":/"<<name<<"</DataItem></Attribute>\n";
                out<<"</Grid></Domain></Xdmf>\n";
                out.close();
                if(!out)throw std::runtime_error("initial XMF write failure");
            } catch(const std::exception& e) {okay=false;error=e.what();}
        }
        if(!d.all(okay))throw std::runtime_error(d.rank==0?error:"rank-0 initial snapshot failed");
    }

    // ---- 8D. Record evaluated thermodynamics and the inherited unit convention ----
    inline void write_thermodynamic_report(const Config& config,const std::filesystem::path& path) {
        std::ofstream out(path);
        if(!out)throw std::runtime_error("cannot write thermodynamics report");
        out<<std::setprecision(17);
        out<<"thermo_mode="<<config.thermo_mode<<"\nrun_temperature_K="<<config.T<<"\nVm_m3_per_mol="<<config.phase.Vm<<"\n";
        out<<"energy_convention=legacy v0.3 molar-code convention; chemical term Psi/Vm, c=atomic fraction\n";
        out<<"unit_scope=this preserves the inherited conversion; it does not certify the material's physical molar volume\n";
        out<<"thermal_history=isothermal at the run temperature; no finite-rate cooling segment\n";
        if(config.thermo_mode=="eutectoid_reference") {
            out<<"T_eutectoid_K="<<config.T_eutectoid<<"\nundercooling_K="<<config.undercooling
               <<"\nmax_undercooling_K="<<config.max_undercooling<<"\ntemperature_law="<<config.temperature_law<<"\n";
            if(config.temperature_law=="phase_diagram_fixed_curvature") {
                out<<"phase_diagram_slope_units=carbon mole fraction/K (dc_eq/dT)\n";
                out<<"phase_diagram_curvature=fixed at supplied reference Aeq\n";
                out<<"phase_diagram_gamma_reference=ABC held at reference values\n";
                out<<"phase_diagram_reference_c_alpha_theta_gamma="<<config.phase_diagram.reference_composition[0]<<','<<config.phase_diagram.reference_composition[1]<<','<<config.phase_diagram.reference_composition[2]<<"\n";
                out<<"phase_diagram_slopes_alpha_ag_gamma_ag_theta_tg_gamma_tg=";
                for(int j=0;j<4;++j)out<<(j?",":"")<<config.phase_diagram.composition_slope[j];
                out<<"\nphase_diagram_reference_residual=preserved, not retuned\n";
            }
        }
        const std::array<std::string,3> names{"alpha","theta","gamma"};
        for(int i=0;i<3;++i) {
            const auto& n=names[i];
            if(config.thermo_mode=="eutectoid_reference") {
                auto p=config.reference_thermo[i]; auto slope=config.temperature_slope[i];
                out<<"reference_"<<n<<"_ABC="<<p.free_energy_quadratic<<','<<p.free_energy_linear<<','<<p.free_energy_constant<<"\n";
                if(config.temperature_law=="linear_gp") out<<"dq_dT_"<<n<<"_q2_q1_q0="<<slope.mu_squared_coefficient<<','<<slope.mu_coefficient<<','<<slope.constant_energy<<"\n";
            }
            auto p=config.thermo[i];auto grand_potential=to_grand_potential(p);
            out<<"effective_"<<n<<"_ABC="<<p.free_energy_quadratic<<','<<p.free_energy_linear<<','<<p.free_energy_constant<<"\n";
            out<<"effective_"<<n<<"_q2_q1_q0="<<grand_potential.mu_squared_coefficient<<','<<grand_potential.mu_coefficient<<','<<grand_potential.constant_energy<<"\n";
            out<<"effective_"<<n<<"_chi="<<p.chi()<<"\nD_"<<n<<"="<<p.diffusivity<<"\n";
        }
        out<<"initial_composition_mode="<<config.composition_mode<<"\n";
        if(config.composition_mode!="mean") out<<"initial_mu_at_run_temperature="<<config.initial_mu<<"\n";
        else out<<"initial_mu_at_run_temperature=solved after smoothing from specified mean_carbon\n";
        out<<"initialization_note=one local mu field; tip samples at distinct positions need not be equal\n";
        out<<"tau_at_ag_tg="<<config.phase.tau[0]<<','<<config.phase.tau[1]<<','<<config.phase.tau[2]<<"\n";
        out<<"tau_note=explicit input, not recalibrated by changing temperature; physical applicability remains a separate check\n";
        out.close();
        if(!out)throw std::runtime_error("thermodynamics report write failed");
    }
}
// namespace pearl


// END ORIGINAL FILE: include/io.hpp

// ============================================================================
// 9. MAIN PROGRAM — command line, run setup, integration, and clean shutdown
// ============================================================================
// Read one immutable input snapshot; create a new output directory; initialize
// or restart; advance to output/end times; save an explicit final checkpoint.
// Solver completion is not a certificate of cooperative growth.

// BEGIN ORIGINAL FILE: main.cpp
// Main loop: initialize/restart -> phase proposal -> implicit diffusion -> checks -> output.
#include <chrono>
#include <iostream>
#include <csignal>
using namespace pearl;
volatile std::sig_atomic_t interrupted=0;
void on_signal(int) {
    interrupted=1;
}
int main(int argc,char** argv) {
#ifdef PEARL_MPI
    MPI_Init(&argc,&argv);
#endif
    int status=0;
    try {
        if(argc<3)throw std::invalid_argument("Usage: pearl input.in NEW_OUTPUT_DIR [--restart frame.h5] [--px N --py N] [--max-steps N] [--thermo-only]");
        std::string input=argv[1],output=argv[2],restart;
        bool thermo_only=false;
        int px=1,py=1;
        long long stop_steps=-1;
        for(int k=3;k<argc;++k) {
            std::string key=argv[k];
            if(key=="--thermo-only") {thermo_only=true;continue;}
            if(k+1==argc)throw std::invalid_argument("missing CLI value");
            std::string val=argv[++k];
            if(key=="--restart")restart=val;
            else if(key=="--px")px=std::stoi(val);
            else if(key=="--py")py=std::stoi(val);
            else if(key=="--max-steps")stop_steps=std::stoll(val);
            else throw std::invalid_argument("unknown CLI flag: "+key);
        }
        if(thermo_only&&!restart.empty())throw std::invalid_argument("--thermo-only cannot be combined with --restart");
        Config config=read_config(input);
        if(stop_steps==0||stop_steps< -1)throw std::invalid_argument("max-steps must be positive");
        Domain d(config.nx,config.ny,config.dx,config.dy,px,py);
        bool outok=true;
        if(d.rank==0) {
            try {
                if(std::filesystem::exists(output))throw std::runtime_error("output path already exists; refuse overwrite");
                std::filesystem::create_directories(output);
                std::filesystem::copy_file(input,std::filesystem::path(output)/"Input.in");
                write_thermodynamic_report(config,std::filesystem::path(output)/"thermodynamics.txt");
            }
            catch(...) {
                outok=false;
            }
        }
        if(!d.all(outok))throw std::runtime_error("cannot create NEW output directory");
        if(thermo_only) {
            if(d.rank==0)std::cout<<"THERMO_ONLY: report written; no initialization or evolution performed\n";
        } else {
        State state(d);
        SmoothReport smooth_report;
        InitializationObserver initial_writer=[&](const Domain& dom,const State& st,const char* stage,int pseudo_step) {
            write_initial_phases(dom,st,std::filesystem::path(output)/(std::string(stage)+".h5"),stage,pseudo_step);
        };
        double mu_reference=restart.empty()?initialize(d,state,config,&smooth_report,initial_writer):read_checkpoint(d,state,config,restart);
        if(restart.empty()&&d.rank==0) {
            std::ofstream initlog(std::filesystem::path(output)/"initialization_report.txt");
            initlog<<std::setprecision(17);
            initlog<<"initialization="<<config.initialization<<"\n";
            initlog<<"physical_time_after_preprocess="<<state.time<<"\n";
            initlog<<"solver_step_after_preprocess="<<state.step<<"\n";
            initlog<<"smooth_steps_completed="<<smooth_report.steps<<"\n";
            initlog<<"smooth_energy_converged="<<(smooth_report.converged?1:0)<<"\n";
            initlog<<"smooth_energy_initial="<<smooth_report.energy_initial<<"\n";
            initlog<<"smooth_energy_final="<<smooth_report.energy_final<<"\n";
            initlog<<"smooth_final_relative_energy_change="<<smooth_report.final_relative_energy_change<<"\n";
            initlog<<"smooth_final_max_phase_change="<<smooth_report.final_max_phase_change<<"\n";
            initlog<<"phase_fraction_before="<<smooth_report.phase_fraction_before[0]<<','<<smooth_report.phase_fraction_before[1]<<','<<smooth_report.phase_fraction_before[2]<<"\n";
            initlog<<"phase_fraction_after="<<smooth_report.phase_fraction_after[0]<<','<<smooth_report.phase_fraction_after[1]<<','<<smooth_report.phase_fraction_after[2]<<"\n";
            initlog<<"mu_reference="<<mu_reference<<"\n";
            initlog<<"NOTE=smoothing is capillarity-only preprocessing pseudo-time; carbon diffusion is disabled and physical time is reset to zero before coupled evolution.\n";
        }
        if(state.time>=config.end_time)throw std::runtime_error("checkpoint time >= end_time");
        std::signal(SIGINT,on_signal);
        std::signal(SIGTERM,on_signal);
        const long long first_step=state.step;
        long long frames=0;
        double next_output=(std::floor((state.time+1e-12*config.output_dt)/config.output_dt)+1.0)*config.output_dt;
        StepInfo last;
        long long cg_total=0;
        int cg_peak=0;
        const auto start=std::chrono::steady_clock::now();
        std::ofstream log;
        if(d.rank==0) {
            log.open(std::filesystem::path(output)/"diagnostics.csv");
            if(!log)throw std::runtime_error("cannot open diagnostics");
            log<<"time,step,dt_next,mean_carbon,mass_ledger_error,reduced_energy,phi_sum_error,c_min,c_max,mu_min,mu_max,cg_last,cg_peak,closure_last,rejected,shift_distance,wall_s\n";
            std::cout<<"resources: ranks="<<d.size<<" layout="<<px<<"x"<<py<<" grid="<<config.nx<<"x"<<config.ny<<"\n";
        }
        auto snapshot=[&]() {
            double minc=1e300,maxc=-1e300,minmu=1e300,maxmu=-1e300,sumerr=0;
            for(int y=1;y<=d.ny;++y)for(int x=1;x<=d.nx;++x) {
                minc=std::min(minc,state.carbon(x,y));
                maxc=std::max(maxc,state.carbon(x,y));
                minmu=std::min(minmu,state.mu(x,y));
                maxmu=std::max(maxmu,state.mu(x,y));
                auto p=cell_phi(state.phi,x,y);
                sumerr=std::max(sumerr,std::abs(p[0]+p[1]+p[2]-1));
            }
            minc=d.minimum(minc);
            maxc=d.maximum(maxc);
            minmu=d.minimum(minmu);
            maxmu=d.maximum(maxmu);
            sumerr=d.maximum(sumerr);
            double mass=total_carbon(d,state),ledger=mass+state.carbon_removed-state.carbon_added-state.initial_carbon;
            double energy=reduced_energy(d,state,config.thermo,config.phase,mu_reference,config.control.phase_only);
            double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            std::ostringstream fn;
            fn<<"frame_"<<std::setfill('0')<<std::setw(6)<<frames++<<".h5";
            write_frame(d,state,config,mu_reference,std::filesystem::path(output)/fn.str());
            if(d.rank==0) {
                log<<std::setprecision(17)<<state.time<<','<<state.step<<','<<state.dt<<','<<mass/(d.gx*d.gy*d.dx*d.dy)<<','<<ledger<<','<<energy<<','<<sumerr<<','<<minc<<','<<maxc<<','<<minmu<<','<<maxmu<<','<<last.cg_iterations<<','<<cg_peak<<','<<last.closure<<','<<state.rejected<<','<<state.shift_distance<<','<<wall<<'\n';
                log.flush();
                std::cout<<"t="<<state.time<<" step="<<state.step<<" dt_next="<<state.dt<<" mass_error="<<ledger<<" CG="<<last.cg_iterations<<" rejected="<<state.rejected<<"\n";
            }
        };
        snapshot();
        bool stopped=false;
        double saved_time=state.time;
        while(state.time<config.end_time) {
            bool stop=d.maximum(interrupted?1.:0.)>.5||(stop_steps>0&&state.step-first_step>=stop_steps);
            if(stop) {
                stopped=true;
                break;
            }
            if(state.step>=config.max_steps)throw std::runtime_error("configured max_steps reached before end_time");
            last=advance(d,state,config.thermo,config.phase,config.control,mu_reference,std::min(config.end_time,next_output));
            cg_total+=last.cg_iterations;
            cg_peak=std::max(cg_peak,last.cg_iterations);
            double eps=1e-12*std::max(config.end_time,config.output_dt);
            if(state.time>=next_output-eps||state.time>=config.end_time-eps) {
                // Snap only the clock within rounding tolerance, not any physical field.
                if(std::abs(state.time-next_output)<=eps)state.time=next_output;
                if(std::abs(state.time-config.end_time)<=eps)state.time=config.end_time;
                maybe_shift(d,state,config);
                snapshot();
                saved_time=state.time;
                next_output+=config.output_dt;
            }
        }
        if(state.time!=saved_time)snapshot();
        write_frame(d,state,config,mu_reference,std::filesystem::path(output)/"checkpoint.h5");
        if(d.rank==0) {
            std::ofstream fin(std::filesystem::path(output)/(stopped?"STOPPED.txt":"COMPLETED.txt"));
            fin<<std::setprecision(17)<<"status="<<(stopped?"USER_STOPPED":"HORIZON_COMPLETED")<<"\ntime="<<state.time<<"\nstep="<<state.step<<"\nCG_total="<<cg_total<<"\nrejected="<<state.rejected<<"\nThis is a solver status, not certification of cooperative growth.\n";
        }
        } // evolution (not --thermo-only)
    }
    catch(const std::exception& e) {
        std::cerr<<"ERROR: "<<e.what()<<'\n';
        status=1;
#ifdef PEARL_MPI
        MPI_Abort(MPI_COMM_WORLD,1);
#endif
    }
#ifdef PEARL_MPI
    MPI_Finalize();
#endif
    return status;
}

// END ORIGINAL FILE: main.cpp
