import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import openmdao.api as om
from openmdao.utils.spline_distributions import cell_centered
from julia.OpenMDAO import make_component
from julia import CCBladeLoadingExample
from paropt.paropt_driver import ParOptDriver

from ccblade_openmdao_examples.structural_group import StructuralGroup
from ccblade_openmdao_examples.gxbeam_openmdao_component import MassComp
from ccblade_openmdao_examples.ccblade_openmdao_component import BEMTRotorCAComp


def get_problem():
    B = 3  # Number of blades.
    D = 24.0*0.0254  # Diameter in meters.
    P = 16.0*0.0254  # Pitch in meters (used in the twist distribution).
    c = 1.5*0.0254   # Constant chord in meters.
    Rtip = 0.5*D
    Rhub = 0.2*Rtip  # Just guessing on the hub diameter.

    # Not sure about the atmospheric conditions, so I'll just use the [ICAO standard
    # atmosphere at sealevel.](https://www.engineeringtoolbox.com/international-standard-atmosphere-d_985.html)
    p0 = 101325.0
    T0 = 273.15 + 15.0
    gam = 1.4
    speedofsound = np.sqrt(gam*287.058*T0)
    rho0 = gam*p0/speedofsound**2
    mu = rho0*1.461e-5

    # Operating parameters for this case.
    rpm = 7200.0
    M_infty = 0.11
    v = M_infty*speedofsound  # axial velocity in m/sec.
    omega = 2*np.pi/60.0*rpm  # propeller rotation rate in rad/sec.
    pitch = 0.0

    # Target thrust value in Newtons.
    thrust_target = 97.246

    # Get the initial blade geometry.
    num_cp = 8
    radii_cp0 = np.linspace(Rhub, Rtip, num_cp)
    chord_cp0 = c*np.ones(num_cp)
    theta_cp0 = np.arctan(P/(np.pi*D*radii_cp0/Rtip))

    num_operating_points = 1
    num_radial = 30
    nnodes = num_radial
    nelems = nnodes-1
    num_stress_eval_points = 20

    prob = om.Problem()

    ivc = om.IndepVarComp()
    ivc.add_output("Rhub", val=Rhub, units="m")
    ivc.add_output("Rtip", val=Rtip, units="m")
    ivc.add_output("radii_cp", val=radii_cp0, units="inch")
    ivc.add_output("chord_cp", val=chord_cp0, units="inch")
    ivc.add_output("theta_cp", val=theta_cp0, units="rad")
    ivc.add_output("v", val=v, shape=num_operating_points, units="m/s")
    ivc.add_output("omega", val=omega, shape=num_operating_points, units="rad/s")
    ivc.add_output("pitch", val=pitch, shape=num_operating_points, units="rad")
    prob.model.add_subsystem("ivc", ivc, promotes=["*"])

    x_cp = np.linspace(0.0, 1.0, num_cp)
    x_interp = cell_centered(nelems, 0.0, 1.0)
    interp_options = {"delta_x": 0.1}
    comp = om.SplineComp(method="akima", interp_options=interp_options, x_cp_val=x_cp, x_interp_val=x_interp)
    comp.add_spline(y_cp_name="radii_cp", y_interp_name="radii", y_units="inch")
    comp.add_spline(y_cp_name="chord_cp", y_interp_name="chord", y_units="inch")
    comp.add_spline(y_cp_name="theta_cp", y_interp_name="theta", y_units="rad")
    prob.model.add_subsystem("akima_comp", comp,
                             promotes_inputs=["radii_cp", "chord_cp", "theta_cp"],
                             promotes_outputs=["radii", "chord", "theta"])

    struc_group = StructuralGroup(nnodes=nnodes, num_stress_eval_points=num_stress_eval_points)
    prob.model.add_subsystem("structural_group", struc_group,
                             promotes_inputs=["omega", "Tp", "Np", "chord", "theta"],
                             promotes_outputs=["sigma1", "m"])

    af_fname = "../data/xf-n0012-il-500000.dat"
    comp = make_component(
        BEMTRotorCAComp(
            af_fname=af_fname, cr75=c/Rtip, Re_exp=0.6,
            num_operating_points=num_operating_points, num_blades=B,
            num_radial=num_radial, rho=rho0, mu=mu, speedofsound=speedofsound))
    prob.model.add_subsystem("bemt_rotor_comp", comp, promotes_inputs=["Rhub", "Rtip", "radii", "chord", "theta", "v", "omega", "pitch"], promotes_outputs=["thrust", "torque", "efficiency"])

    # TODO: connect the aero forces to the structural group

    # Aggregate the stress
    prob.model.add_subsystem('ks', om.KSComp(width=nelems*num_stress_eval_points*2,
                                             add_constraint=False, ref=1.0,
                                             units=None))
    prob.model.connect('sigma1', 'ks.g')

    prob.model.linear_solver = om.DirectSolver()
    prob.driver = ParOptDriver()

    # Set the optimization options
    options = {'algorithm': 'tr',
               'tr_max_size': 0.1,
               #'tr_min_size': 1e-4,
               'tr_adaptive_gamma_update': True,
               'qn_diag_type': 'yts_over_sts',
               'tr_use_soc': True,
               'tr_max_iterations': 1000,
               'penalty_gamma': 5.0,
               'tr_penalty_gamma_max': 10.0,
               'qn_subspace_size': 20,
               'qn_type': 'bfgs',
               'abs_res_tol': 1e-8,
               'starting_point_strategy': 'affine_step',
               'barrier_strategy': 'mehrotra_predictor_corrector',
               'use_line_search': False,
               'max_major_iters': 200}
    for key in options:
        prob.driver.options[key] = options[key]

    # Lower and upper limits on the chord design variable, in inches.
    chord_lower = 0.5
    chord_upper = 2.0

    # Lower and upper limits on the twist design variable, radians.
    theta_lower =  0.0*np.pi/180.0
    theta_upper =  85.0*np.pi/180.0

    # Define the mass threshold using a percent of the upper chord limit
    mass_frac = 0.40
    m_max = 2600.0 * ((chord_upper/100.0)**2) * 821.8 * 12.0*0.0254 # maximum possible mass = rho * (c_max/c_ref)^2 * A_ref * span

    # Set the optimization problem formulation
    prob.model.add_design_var("chord_cp", lower=chord_lower, upper=chord_upper, ref=1e0, units="inch")
    prob.model.add_design_var("twist_cp", lower=theta_lower, upper=theta_upper, ref=1e0, units="rad")

    prob.model.add_objective("efficiency", ref=-1e0)
    prob.model.add_constraint("thrust", lower=thrust_target, upper=thrust_target, units="N", ref=1e2)
    prob.model.add_constraint("ks.KS", upper=1.0, ref=1e0)

    prob.setup()

    return x, prob, Np, Tp


if __name__ == "__main__":
    x, p, Np, Tp = get_problem()
    p.run_driver()

    chord = p.get_val("chord", units="inch")[0]
    theta = p.get_val("twist", units="deg")[0]
    print("mass = ", p.get_val("m", units="kg"))
    print("KS(sigma1)/sigma_y = ", p.get_val("ks.KS"))
    print("max(sigma1) = ", np.amax(p.get_val("sigma1")))

    span = 12.0*0.0254
    nelems = len(x)-1
    dl = span/nelems
    xe = (np.array(x[0:-1]) + dl/2)/0.0254 # x-location of each element midpoint

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
    ax0.plot(xe, chord, color='tab:blue')
    ax1.plot(xe, theta, color='tab:blue')

    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    ax0.tick_params(axis='x', direction='out', length=0.0, width=0.0)
    #ax0.xaxis.set_visible(False)
    ax0.yaxis.set_ticks_position('left')
    ax0.xaxis.set_ticks_position('bottom')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax0.grid(True)
    ax1.grid(True)

    # ticks_loc = ax1.get_xticks().tolist()
    # ax1.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # ax1.set_xticklabels([label_format.format(x) for x in ticks_loc])
    # ax0.set_xticklabels([])

    ax0.set_ylabel('Chord (in)')
    ax1.set_xlabel(r'$x_1$ (in)')
    ax1.set_ylabel('Twist (deg.)')

    plt.savefig("chord_theta.pdf", transparent=True)

    df = pd.DataFrame({'chord':chord, 'theta':theta})
    df.to_csv('chord_theta.csv', index=False)

    Tpe = np.interp(xe, np.array(x)/0.0254, Tp)
    Npe = np.interp(xe, np.array(x)/0.0254, Np)

    # Plot everything
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, sharex=True, constrained_layout=True)

    E = 72.4e9
    Hyy = E*p.get_val('structural_group.Iyy')
    Hzz = E*p.get_val('structural_group.Izz')
    Hyz = E*p.get_val('structural_group.Iyz')
    Fx = p.get_val('structural_group.solver_comp.Fx')
    My = p.get_val('structural_group.solver_comp.My')
    Mz = p.get_val('structural_group.solver_comp.Mz')
    sigma1 = p.get_val('sigma1')
    sigma_x1 = np.zeros(nelems)
    num_stress_eval_points = 20
    for i in range(nelems):
        sigma_x1[i] = np.amax(sigma1[2*i*num_stress_eval_points:2*(i+1)*num_stress_eval_points])

    ax0.plot(xe, Tpe, label=r'$T_p$')
    ax0.plot(xe, Npe, label=r'$N_p$')

    ax1.plot(xe, Hyy, label=r'$H_{yy}$')
    ax1.plot(xe, Hzz, label=r'$H_{zz}$')
    ax1.plot(xe, Hyz, label=r'$H_{yz}$')

    ax2r = ax2.twinx()
    ax2r.plot(xe, Fx, c='C0', label=r'$N_{x}$')
    ax2.plot(xe, My, c='C1', label=r'$M_{y}$')
    ax2.plot(xe, Mz, c='C2', label=r'$M_{z}$')

    ax3.plot(xe, sigma_x1, label=r'$\sigma_1$')

    for ax in [ax0, ax1, ax2]:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='x', direction='out', length=0.0, width=0.0)
        #ax.xaxis.set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.grid(True)

    ax2r.spines['top'].set_visible(False)
    ax2r.spines['bottom'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)

    ax2r.grid(None)
    ax2.grid(True)
    ax3.grid(True)

    # ticks_loc = ax1.get_xticks().tolist()
    # ax1.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # ax1.set_xticklabels([label_format.format(x) for x in ticks_loc])
    # ax0.set_xticklabels([])

    fs = 6
    ax0.set_ylabel('Aero forces (N/m)', fontsize=fs)
    ax1.set_ylabel(r'Bending stiffness (N-m$^2$)', fontsize=fs)
    ax2.set_ylabel('Bending moment (N-m)', fontsize=fs)
    ax2r.set_ylabel('Axial force (N)', fontsize=fs)
    ax3.set_ylabel(r'$\sigma_1(x_1)/\sigma_y$')
    ax3.set_xlabel(r'$x_1$ (in)')

    leg0 = ax0.legend()
    leg1 = ax1.legend()
    leg3 = ax1.legend()

    h1, l1 = ax2r.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    leg2 = ax2.legend([h1[0]]+h2, [l1[0]]+l2)

    plt.savefig("opt_vals.pdf", transparent=True)