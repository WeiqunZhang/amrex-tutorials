#include <AMReX.H>
#include <AMReX_ForkJoin.H>
#include <AMReX_Math.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParReduce.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <string>
#include <type_traits>

using namespace amrex;

namespace {

// Parameters of the variable-coefficient problem  a*u - b*div(beta*grad u) = f,
// where beta jumps smoothly from 1 to sigma across the sphere r = 0.25.
struct Problem
{
    Real a = 1.e-3_rt;
    Real b = 1._rt;
    Real sigma = 10._rt;               // beta outside the sphere
    Real theta = 0.5_rt*std::log(3._rt)/0.05_rt; // steepness of the jump
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real radius (GpuArray<Real,AMREX_SPACEDIM> const& x) noexcept
{
    Real r2 = 0._rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { r2 += (x[d]-0.5_rt)*(x[d]-0.5_rt); }
    return std::sqrt(r2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real beta_coef (GpuArray<Real,AMREX_SPACEDIM> const& x, Problem const& p) noexcept
{
    return 0.5_rt*(p.sigma-1._rt)*std::tanh(p.theta*(radius(x)-0.25_rt))
        +  0.5_rt*(p.sigma+1._rt);
}

// Exact solution u = prod sin(2 pi x_d) + 1/4 prod sin(4 pi x_d), zero on the boundary
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real exact_solution (GpuArray<Real,AMREX_SPACEDIM> const& x) noexcept
{
    constexpr Real pi = Math::pi<Real>();
    Real s2 = 1._rt, s4 = 1._rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        s2 *= std::sin(2._rt*pi*x[d]);
        s4 *= std::sin(4._rt*pi*x[d]);
    }
    return s2 + 0.25_rt*s4;
}

// f = a*u - b*(beta*lap(u) + dbeta/dr * (x-0.5)/r . grad(u))
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real rhs_value (GpuArray<Real,AMREX_SPACEDIM> const& x, Problem const& p) noexcept
{
    constexpr Real pi = Math::pi<Real>();
    constexpr Real tpi = 2._rt*pi;
    constexpr Real fpi = 4._rt*pi;

    Real s2 = 1._rt, s4 = 1._rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        s2 *= std::sin(tpi*x[d]);
        s4 *= std::sin(fpi*x[d]);
    }
    Real neg_lap = AMREX_SPACEDIM*tpi*tpi*(s2 + s4);

    Real r = radius(x);
    Real dbdr_grad = 0._rt;
    if (r > 0._rt) {
        Real ch = std::cosh(p.theta*(r-0.25_rt));
        Real dbdr = 0.5_rt*(p.sigma-1._rt)*p.theta/(ch*ch);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            Real g2 = tpi*std::cos(tpi*x[d]);
            Real g4 = pi *std::cos(fpi*x[d]);
            for (int e = 0; e < AMREX_SPACEDIM; ++e) {
                if (e != d) {
                    g2 *= std::sin(tpi*x[e]);
                    g4 *= std::sin(fpi*x[e]);
                }
            }
            dbdr_grad += dbdr*(x[d]-0.5_rt)/r*(g2 + g4);
        }
    }
    return p.a*(s2 + 0.25_rt*s4) + p.b*(beta_coef(x,p)*neg_lap - dbdr_grad);
}

// Physical location of index iv, cell-centered in all directions except
// the nodal direction inode (if any).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuArray<Real,AMREX_SPACEDIM>
location (IntVect const& iv, int inode, GpuArray<Real,AMREX_SPACEDIM> const& problo,
          GpuArray<Real,AMREX_SPACEDIM> const& dx) noexcept
{
    GpuArray<Real,AMREX_SPACEDIM> x;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real offset = (d == inode) ? 0._rt : 0.5_rt;
        x[d] = problo[d] + (static_cast<Real>(iv[d]) + offset)*dx[d];
    }
    return x;
}

// Solve a single-component problem with MLMG. Inside a ForkJoin task, MLMG
// uses the task's MPI communicator (ParallelContext) automatically.
void solve (MultiFab& phi, MultiFab const& rhs, MultiFab const& alpha,
            Array<MultiFab const*,AMREX_SPACEDIM> const& beta,
            Geometry const& geom, Problem const& prob, int verbose)
{
    MLABecLaplacian mlabec({geom}, {phi.boxArray()}, {phi.DistributionMap()});
    mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                     LinOpBCType::Dirichlet,
                                     LinOpBCType::Dirichlet)},
                       {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                     LinOpBCType::Dirichlet,
                                     LinOpBCType::Dirichlet)});
    mlabec.setLevelBC(0, nullptr); // homogeneous Dirichlet
    mlabec.setScalars(prob.a, prob.b);
    mlabec.setACoeffs(0, alpha);
    mlabec.setBCoeffs(0, beta);

    constexpr Real tol_rel = std::is_same_v<Real,float> ? 1.e-5_rt : 1.e-10_rt;
    MLMG mlmg(mlabec);
    mlmg.setVerbose(verbose);
    mlmg.solve({&phi}, {&rhs}, tol_rel, 0._rt);
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int ncomp = 4;       // number of independent problems
        int ntasks = 2;      // number of ForkJoin tasks
        int fj_verbose = 0;
        int mlmg_verbose = 0;
        std::string task_output_dir; // if set, each task writes its output there
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("ncomp", ncomp);
            pp.query("ntasks", ntasks);
            pp.query("fj_verbose", fj_verbose);
            pp.query("mlmg_verbose", mlmg_verbose);
            pp.query("task_output_dir", task_output_dir);
        }
        // Every task needs at least one rank and one component.
        ntasks = std::clamp(ntasks, 1, std::min(ParallelDescriptor::NProcs(), ncomp));

        Box domain(IntVect(0), IntVect(n_cell-1));
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);
        Geometry geom(domain, RealBox(AMREX_D_DECL(0._rt,0._rt,0._rt),
                                      AMREX_D_DECL(1._rt,1._rt,1._rt)),
                      CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});
        auto const problo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();
        Problem const prob{};

        // Component n holds problem n, whose solution is (n+1)*u.
        MultiFab phi(ba, dm, ncomp, 0);
        MultiFab rhs(ba, dm, ncomp, 0);
        auto const& rhs_ma = rhs.arrays();
        ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k)
        {
            auto x = location(IntVect(AMREX_D_DECL(i,j,k)), -1, problo, dx);
            Real f = rhs_value(x, prob);
            for (int n = 0; n < ncomp; ++n) { rhs_ma[bi](i,j,k,n) = Real(n+1)*f; }
        });

        // Coefficients shared by all problems: alpha (cells) and beta (faces).
        MultiFab alpha(ba, dm, 1, 0);
        alpha.setVal(1._rt);
        Array<MultiFab,AMREX_SPACEDIM> beta;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            beta[idim].define(amrex::convert(ba, IntVect::TheDimensionVector(idim)), dm, 1, 0);
            auto const& beta_ma = beta[idim].arrays();
            ParallelFor(beta[idim], [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k)
            {
                auto x = location(IntVect(AMREX_D_DECL(i,j,k)), idim, problo, dx);
                beta_ma[bi](i,j,k) = beta_coef(x, prob);
            });
        }
        Gpu::streamSynchronize();

        // Split the ranks evenly into ntasks tasks.
        ForkJoin fj(ntasks);
        fj.SetVerbose(fj_verbose);
        fj.set_task_output_dir(task_output_dir);

        // split:     each task gets a subset of the components
        // duplicate: each task gets a full copy
        // (single:   only one task gets a copy; not used here)
        // Intent controls the copies before (in) and after (out) the fork.
        fj.reg_mf(rhs, "rhs", ForkJoin::Strategy::split, ForkJoin::Intent::in);
        fj.reg_mf(phi, "phi", ForkJoin::Strategy::split, ForkJoin::Intent::out);
        fj.reg_mf(alpha, "alpha", ForkJoin::Strategy::duplicate, ForkJoin::Intent::in);
        fj.reg_mf_vec(Vector<MultiFab const*>{AMREX_D_DECL(&beta[0],&beta[1],&beta[2])},
                      "beta", ForkJoin::Strategy::duplicate, ForkJoin::Intent::in);

        amrex::Print() << "Solving " << ncomp << " problems with " << fj.NTasks()
                       << " ForkJoin tasks on " << ParallelDescriptor::NProcs() << " ranks\n";

        fj.fork_join([&] (ForkJoin& f)
        {
            // Here, ParallelContext refers to the task's ranks only, and
            // amrex::Print prints from the task's first rank.
            MultiFab& my_phi = f.get_mf("phi");
            MultiFab const& my_rhs = f.get_mf("rhs");
            MultiFab const& my_alpha = f.get_mf("alpha");
            auto my_beta = f.get_mf_vec("beta");
            auto comps = f.ComponentBounds("phi");

            amrex::Print() << "  Task " << f.MyTask() << ": "
                           << ParallelContext::NProcsSub() << " rank(s), components ["
                           << comps.lo << ", " << comps.hi << ")\n";

            for (int n = 0; n < my_phi.nComp(); ++n) {
                MultiFab phi_n(my_phi, amrex::make_alias, n, 1);
                MultiFab rhs_n(my_rhs, amrex::make_alias, n, 1);
                phi_n.setVal(0._rt);
                solve(phi_n, rhs_n, my_alpha,
                      {AMREX_D_DECL(my_beta[0],my_beta[1],my_beta[2])},
                      geom, prob, mlmg_verbose);
            }
        });
        // After the join, phi holds the solutions from all tasks.

        Real const tol = 1.e-2_rt;
        auto const& phi_ma = phi.const_arrays();
        for (int n = 0; n < ncomp; ++n) {
            Real err = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, phi, IntVect(0),
                [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) -> GpuTuple<Real>
                {
                    auto x = location(IntVect(AMREX_D_DECL(i,j,k)), -1, problo, dx);
                    return { std::abs(phi_ma[bi](i,j,k,n) - Real(n+1)*exact_solution(x)) };
                });
            ParallelDescriptor::ReduceRealMax(err);
            err /= Real(n+1);
            amrex::Print() << "Problem " << n << ": max error of phi/(n+1) = " << err << "\n";
            if (!(err < tol)) { status = 1; }
        }
        if (status != 0) { amrex::Print() << "FAILED\n"; }
    }
    amrex::Finalize();
    return status;
}
