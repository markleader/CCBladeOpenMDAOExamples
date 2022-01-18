import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import openmdao.api as om
from openmdao.utils.spline_distributions import cell_centered
from julia.OpenMDAO import make_component
from julia import CCBladeLoadingExample
from paropt.paropt_driver import ParOptDriver

from ccblade_openmdao_examples.structural_group import StructuralGroup
from ccblade_openmdao_examples.gxbeam_openmdao_component import MassComp


def get_problem(optimizer="SNOPT"):

    # Propeller dimensions
    D = 24.0*0.0254  # Diameter in meters.
    Rtip = 0.5*D
    Rhub = 0.2*Rtip  # Just guessing on the hub diameter.

    # Compute the aerodynamic loads
    x, Np, Tp = CCBladeLoadingExample.run_ccblade()
    nelems = len(x)
    num_stress_eval_points = 20

    # Initialize the design variables
    chord = 1.0*np.ones(nelems)  # (inch)
    theta = (45.0*np.pi/180.0)*np.ones(nelems)  # (rad)

    # Define the angular rotation
    rpm = 7110.0
    omega = rpm*2*np.pi/60

    prob = om.Problem()

    ivc = om.IndepVarComp()
    ivc.add_output("omega", omega, units="rad/s")
    ivc.add_output("Tp", Tp, units="N/m")
    ivc.add_output("Np", Np, units="N/m")
    ivc.add_output("chord", chord, units="inch")
    ivc.add_output("theta", theta, units="rad")
    prob.model.add_subsystem("ivc", ivc, promotes=["*"])

    struc_group = StructuralGroup(nelems=nelems, num_stress_eval_points=num_stress_eval_points, Rhub=Rhub, span=Rtip)
    prob.model.add_subsystem("structural_group", struc_group,
                             promotes_inputs=["omega", "Tp", "Np", "chord", "theta"],
                             promotes_outputs=["sigma_vm", "m"])

    # Aggregate the stress
    prob.model.add_subsystem("ks", om.KSComp(width=nelems*num_stress_eval_points*2,
                                             add_constraint=False, ref=1.0,
                                             units=None))
    prob.model.connect("sigma_vm", "ks.g")

    if optimizer == "SNOPT":
        prob.driver = om.pyOptSparseDriver(optimizer="SNOPT")
    else:  # use ParOpt
        prob.driver = ParOptDriver()

        # Set the optimization options
        options = {"algorithm": "tr",
                "tr_max_size": 0.1,
                #"tr_min_size": 1e-4,
                "tr_adaptive_gamma_update": True,
                "qn_diag_type": "yts_over_sts",
                "tr_use_soc": True,
                "tr_max_iterations": 1000,
                "penalty_gamma": 5.0,
                "tr_penalty_gamma_max": 10.0,
                "qn_subspace_size": 20,
                "qn_type": "bfgs",
                "abs_res_tol": 1e-8,
                "starting_point_strategy": "affine_step",
                "barrier_strategy": "mehrotra_predictor_corrector",
                "use_line_search": False,
                "max_major_iters": 200}
        for key in options:
            prob.driver.options[key] = options[key]

    # Lower and upper limits on the chord design variable, in inches.
    chord_lower = 0.5
    chord_upper = 2.0

    # Lower and upper limits on the twist design variable, radians.
    theta_lower = 5.0*np.pi/180.0
    theta_upper = 85.0*np.pi/180.0

    prob.model.add_design_var("chord", lower=chord_lower, upper=chord_upper, units="inch")
    prob.model.add_design_var("theta", lower=theta_lower, upper=theta_upper, units="rad")

    # Stress-constrained mass minimization
    prob.model.add_objective("m", ref=1e-2)
    prob.model.add_constraint("ks.KS", upper=0.5)

    prob.setup()
    om.n2(prob, show_browser=False, outfile='struc_opt.html')

    return x, prob, Np, Tp

def get_problem_w_splines(optimizer="SNOPT"):

    # Propeller dimensions
    D = 24.0*0.0254  # Diameter in meters.
    Rtip = 0.5*D
    Rhub = 0.2*Rtip  # Just guessing on the hub diameter.

    # Compute the aerodynamic loads
    x, Np, Tp = CCBladeLoadingExample.run_ccblade()
    nelems = len(x)
    num_stress_eval_points = 20

    # Initialize the design variables
    num_cp = 8
    chord_cp0 = 1.0*np.ones(num_cp)  # (inch)
    theta_cp0 = (45.0*np.pi/180.0)*np.ones(num_cp)  # (rad)

    # Define the angular rotation
    rpm = 7110.0
    omega = rpm*2*np.pi/60

    prob = om.Problem()

    ivc = om.IndepVarComp()
    ivc.add_output("omega", omega, units="rad/s")
    ivc.add_output("Tp", Tp, units="N/m")
    ivc.add_output("Np", Np, units="N/m")
    ivc.add_output("chord_cp", chord_cp0, units="inch")
    ivc.add_output("theta_cp", theta_cp0, units="rad")
    prob.model.add_subsystem("ivc", ivc, promotes=["*"])

    x_cp = np.linspace(0.0, 1.0, num_cp)
    x_interp = cell_centered(nelems, 0.0, 1.0)
    interp_options = {"delta_x": 0.1}
    comp = om.SplineComp(method="akima", interp_options=interp_options, x_cp_val=x_cp, x_interp_val=x_interp)
    comp.add_spline(y_cp_name="chord_cp", y_interp_name="chord", y_units="inch")
    comp.add_spline(y_cp_name="theta_cp", y_interp_name="theta", y_units="rad")
    prob.model.add_subsystem("akima_comp", comp,
                             promotes_inputs=["chord_cp", "theta_cp"],
                             promotes_outputs=["chord", "theta"])

    struc_group = StructuralGroup(nelems=nelems, num_stress_eval_points=num_stress_eval_points, Rhub=Rhub, span=Rtip)
    prob.model.add_subsystem("structural_group", struc_group,
                             promotes_inputs=["omega", "Tp", "Np", "chord", "theta"],
                             promotes_outputs=["sigma_vm", "m"])

    # Aggregate the stress
    prob.model.add_subsystem("ks", om.KSComp(width=nelems*num_stress_eval_points*2,
                                             add_constraint=False, ref=1.0,
                                             units=None))
    prob.model.connect("sigma_vm", "ks.g")

    if optimizer == "SNOPT":
        prob.driver = om.pyOptSparseDriver(optimizer="SNOPT")
    else:  # use ParOpt
        prob.driver = ParOptDriver()

        # Set the optimization options
        options = {"algorithm": "tr",
                "tr_max_size": 0.1,
                #"tr_min_size": 1e-4,
                "tr_adaptive_gamma_update": True,
                "qn_diag_type": "yts_over_sts",
                "tr_use_soc": True,
                "tr_max_iterations": 1000,
                "penalty_gamma": 5.0,
                "tr_penalty_gamma_max": 10.0,
                "qn_subspace_size": 20,
                "qn_type": "bfgs",
                "abs_res_tol": 1e-8,
                "starting_point_strategy": "affine_step",
                "barrier_strategy": "mehrotra_predictor_corrector",
                "use_line_search": False,
                "max_major_iters": 200}
        for key in options:
            prob.driver.options[key] = options[key]

    # Lower and upper limits on the chord design variable, in inches.
    chord_lower = 0.5
    chord_upper = 2.0

    # Lower and upper limits on the twist design variable, radians.
    theta_lower = 0.0*np.pi/180.0
    theta_upper = 85.0*np.pi/180.0

    prob.model.add_design_var("chord_cp", lower=chord_lower, upper=chord_upper, units="inch")
    prob.model.add_design_var("theta_cp", lower=theta_lower, upper=theta_upper, units="rad")

    # Stress-constrained mass minimization
    prob.model.add_objective("m", ref=1e-2)
    prob.model.add_constraint("ks.KS", upper=1.0)

    prob.setup()
    om.n2(prob, show_browser=False, outfile='struc_opt_w_splines.html')

    return x, prob, Np, Tp

def run_optimization(use_splines=False, optimizer='SNOPT'):
    # Run the structural optimization problem and plot the outputs

    if use_splines:
        xe, p, Np, Tp = get_problem_w_splines(optimizer=optimizer)
    else:
        xe, p, Np, Tp = get_problem(optimizer=optimizer)
    p.run_driver()

    print("mass = ", p.get_val("m", units="kg"))
    print("KS(sigma_vm)/sigma_y = ", p.get_val("ks.KS"))
    print("max(sigma_vm) = ", np.amax(p.get_val("sigma_vm")))

    # Save the chord and twist distribution to a csv
    if use_splines:
        fname = "chord_theta_{0}_w_splines.csv".format(optimizer)
        chord = p.get_val("chord", units="inch")[0]
        theta = p.get_val("theta", units="deg")[0]
    else:
        fname = "chord_theta_{0}.csv".format(optimizer)
        chord = p.get_val("chord", units="inch")
        theta = p.get_val("theta", units="deg")
    df = pd.DataFrame({"chord":chord, "theta":theta})
    df.to_csv(fname, index=False)

    # Plot the chord and twist distribution
    xe = np.array(xe)/0.0254
    plot_chord_theta(xe, chord, theta)

    # Plot the other values of interest
    plot_extras(p, xe, Tp, Np)

    return


class ConvertScalarToVec(om.ExplicitComponent):

    def initialize(self):
        self.options.declare("nelems", types=int)
        self.options.declare("units")

        return

    def setup(self):

        self.nelems = self.options["nelems"]
        units = self.options["units"]
        self.add_input("val", 0.0, units=units)
        self.add_output("val_vec", np.zeros(self.nelems), units=units)
        self.declare_partials(of="val_vec", wrt="val")

        return

    def compute(self, inputs, outputs):

        val = inputs["val"]
        outputs["val_vec"] = val*np.ones(self.nelems)

        return

    def compute_partials(self, inputs, partials):

        partials["val_vec", "val"] = np.ones(self.nelems)

        return

def get_1d_problem(chord_val=None, theta_val=None, optimizer="SNOPT"):

    # If chord_val is None, treat chord as a design variable
    # If theta_val is None, treat theta as a design variable

    # Propeller dimensions
    D = 24.0*0.0254  # Diameter in meters.
    Rtip = 0.5*D
    Rhub = 0.2*Rtip  # Just guessing on the hub diameter.

    # Compute the aerodynamic loads
    x, Np, Tp = CCBladeLoadingExample.run_ccblade()
    nelems = len(x)
    num_stress_eval_points = 20

    # Initialize the design variables
    if not chord_val:
        init_chord_val = 1.0
    else:
        init_chord_val = chord_val

    if not theta_val:
        init_theta_val = 45.0*np.pi/180.0  # (rad)
    else:
        init_theta_val = theta_val

    # Define the angular rotation
    rpm = 7110.0
    omega = rpm*2*np.pi/60

    # Set dummy loads
    Np = np.zeros(nelems)
    Tp = np.zeros(nelems)

    prob = om.Problem()

    ivc = om.IndepVarComp()
    ivc.add_output("omega", omega, units="rad/s")
    ivc.add_output("Tp", Tp, units="N/m")
    ivc.add_output("Np", Np, units="N/m")
    ivc.add_output("chord_val", init_chord_val, units="inch")
    ivc.add_output("theta_val", init_theta_val, units="rad")
    prob.model.add_subsystem("ivc", ivc, promotes=["*"])

    chord_comp = ConvertScalarToVec(nelems=nelems, units="inch")
    prob.model.add_subsystem("chord_comp", chord_comp)
    prob.model.connect("chord_val", "chord_comp.val")

    theta_comp = ConvertScalarToVec(nelems=nelems, units="rad")
    prob.model.add_subsystem("theta_comp", theta_comp)
    prob.model.connect("theta_val", "theta_comp.val")

    struc_group = StructuralGroup(nelems=nelems, num_stress_eval_points=num_stress_eval_points, Rhub=Rhub, span=Rtip)
    prob.model.add_subsystem("structural_group", struc_group,
                             promotes_inputs=["omega", "Tp", "Np", "chord", "theta"],
                             promotes_outputs=["sigma_vm", "m"])

    prob.model.connect("chord_comp.val_vec", "chord")
    prob.model.connect("theta_comp.val_vec", "theta")

    # Aggregate the stress
    prob.model.add_subsystem("ks", om.KSComp(width=nelems*num_stress_eval_points*2,
                                             add_constraint=False, ref=1.0,
                                             units=None))
    prob.model.connect("sigma_vm", "ks.g")
    #prob.model.connect("structural_group.stress_comp.sigma1", "ks.g")

    if optimizer == "SNOPT":
        prob.driver = om.pyOptSparseDriver(optimizer="SNOPT")
    else:
        prob.driver = ParOptDriver()

        # Set the optimization options
        options = {"algorithm": "tr",
                "tr_max_size": 0.1,
                #"tr_min_size": 1e-4,
                "tr_adaptive_gamma_update": True,
                "qn_diag_type": "yts_over_sts",
                "tr_use_soc": True,
                "tr_max_iterations": 1000,
                "penalty_gamma": 5.0,
                "tr_penalty_gamma_max": 10.0,
                "qn_subspace_size": 20,
                "qn_type": "bfgs",
                "abs_res_tol": 1e-8,
                "starting_point_strategy": "affine_step",
                "barrier_strategy": "mehrotra_predictor_corrector",
                "use_line_search": False,
                "max_major_iters": 200}
        for key in options:
            prob.driver.options[key] = options[key]

    # Lower and upper limits on the chord design variable, in inches.
    chord_lower = 0.5
    chord_upper = 2.0

    # Lower and upper limits on the twist design variable, radians.
    theta_lower = 5.0*np.pi/180.0
    theta_upper = 85.0*np.pi/180.0

    if chord_val is None:
        prob.model.add_design_var("chord_val", lower=chord_lower, upper=chord_upper, ref=1e0, units="inch")
    if theta_val is None:
        prob.model.add_design_var("theta_val", lower=theta_lower, upper=theta_upper, ref=1e0, units="rad")

    # Stress minimization
    prob.model.add_objective("ks.KS", ref=1e0)

    prob.setup()
    om.n2(prob, show_browser=False, outfile='struc_1d.html')

    return x, prob, Np, Tp

def check_theta_opt_val():
    # Run a sweep of uniform twist values and plot the KS stress at each point
    # Then solve an optimization with uniform twist as the only design variable and check that point

    theta_vals = (np.pi/180)*np.linspace(5.0, 85.0, 100)
    ks_vals = np.zeros(len(theta_vals))
    Iyy = np.zeros(len(theta_vals))
    Izz = np.zeros(len(theta_vals))
    Iyz = np.zeros(len(theta_vals))
    Fx = np.zeros(len(theta_vals))
    My = np.zeros(len(theta_vals))
    Mz = np.zeros(len(theta_vals))
    _, p, _, _ = get_1d_problem(theta_val=0.0)

    for i in range(len(theta_vals)):
        p.set_val("theta_val", theta_vals[i])
        p.run_model()
        ks_vals[i] = p.get_val("ks.KS")
        #print(f"Theta = {theta_vals[i]}; KS = {ks_vals[i]}")
        Iyy[i] = p.get_val("structural_group.Iyy", units="inch**4")[0]
        Izz[i] = p.get_val("structural_group.Izz", units="inch**4")[0]
        Iyz[i] = p.get_val("structural_group.Iyz", units="inch**4")[0]
        Fx[i] = np.amax(p.get_val("structural_group.solver_comp.Fx"))
        My[i] = np.amax(p.get_val("structural_group.solver_comp.My"))
        Mz[i] = np.amax(p.get_val("structural_group.solver_comp.Mz"))

    # Reset the starting twist value for optimization
    # _, p, _, _ = get_1d_problem(chord_val=1.0)
    # p.run_driver()
    # theta_opt = p.get_val("theta_val")
    # ks_opt = p.get_val("ks.KS")

    # Hard-code the optimal
    theta_opt = 1.48352986*180.0/np.pi
    ks_opt = 0.34065576

#     fig = plt.figure()
#     ax = plt.subplot(1, 1, 1)
#     ax.plot(theta_vals*180.0/np.pi, ks_vals, label="Value sweep")
#     ax.scatter(theta_opt, ks_opt, color="tab:red", label="Optimal", zorder=100)
#
#     ax.spines["right"].set_visible(False)
#     ax.spines["top"].set_visible(False)
#     ax.set_xlabel("Feather angle (deg.)")
#     ax.set_ylabel(r"KS($\sigma$)")
#     ax.legend(frameon=False)

    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True, constrained_layout=True)
    ax0.plot((180.0/np.pi)*theta_vals, ks_vals, label="Value sweep")
    ax0.scatter(theta_opt, ks_opt, color="tab:red", label="Optimal", zorder=100)

    ax0.spines["right"].set_visible(False)
    ax0.spines["top"].set_visible(False)
    ax0.set_ylabel(r"KS($\sigma$)")
    # ax0.legend(frameon=False)

    ax1.plot((180.0/np.pi)*theta_vals, Iyy, label=r"$I_{yy}$")
    ax1.plot((180.0/np.pi)*theta_vals, Izz, label=r"$I_{zz}$")
    ax1.plot((180.0/np.pi)*theta_vals, Iyz, label=r"$I_{yz}$")

    ax2r = ax2.twinx()
    ax2r.plot((180.0/np.pi)*theta_vals, Fx, c="C0", label=r"$F_{x}$")
    ax2.plot((180.0/np.pi)*theta_vals, My, c="C1", label=r"$M_{y}$")
    ax2.plot((180.0/np.pi)*theta_vals, Mz, c="C2", label=r"$M_{z}$")

    for ax in [ax0, ax1]:
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", direction="out", length=0.0, width=0.0)
        ax.yaxis.set_ticks_position("left")
        ax.xaxis.set_ticks_position("bottom")
        ax.grid(True)

    ax2r.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    ax2.set_xlabel("Feather angle (deg.)")

    ax1.legend()
    h1, l1 = ax2r.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax2.legend([h1[0]]+h2, [l1[0]]+l2)

    plt.savefig("theta_opt_sweep.pdf")

    return

def check_chord_opt_val():
    # Run a sweep of uniform chord values and plot the KS stress at each point
    # Then solve an optimization with uniform chord as the only design variable and check that point

    chord_vals = np.linspace(0.5, 2.0, 16)
    ks_vals = np.zeros(len(chord_vals))
    _, p, _, _ = get_1d_problem(chord_val=1.0)

    for i in range(len(chord_vals)):
        p.set_val("chord_val", chord_vals[i])
        p.run_model()
        ks_vals[i] = p.get_val("ks.KS")

    # # Reset the starting twist value for optimization
    # _, p, _, _ = get_1d_problem(theta_val=0.0)
    # p.run_driver()
    # chord_opt = p.get_val("chord_val")
    # ks_opt = p.get_val("ks.KS")
    # print("Chord = ", chord_opt)
    # print("KS(sigma) = ", ks_opt)
    chord_opt = 0.9471346
    ks_opt = 0.44174087

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    ax.plot(chord_vals, ks_vals, label="Value sweep")
    ax.scatter(chord_opt, ks_opt, color="tab:red", label="Optimal", zorder=100)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xlabel("Chord (inches)")
    ax.set_ylabel(r"KS($\sigma$)")
    ax.legend(frameon=False)
    plt.savefig("chord_opt_sweep.pdf")

    return

def plot_stress_distribution(theta=0.0):
    # Check the stress distribution using the 1D problem

    x, p, _, _ = get_1d_problem(theta_val=theta)
    p.run_model()
    ks = p.get_val("ks.KS")
    #print(f"Theta = {theta}; KS = {ks}")
    sigma_vec = p.get_val("sigma_vm")
    nelems = len(x)
    num_stress_eval_points = int(len(sigma_vec)/(2*nelems))
    sigma = np.zeros((2*num_stress_eval_points, nelems))
    sigma_max = np.zeros(nelems)
    for i in range(nelems):
        sigma_max[i] = np.amax(sigma_vec[2*i*num_stress_eval_points:2*(i+1)*num_stress_eval_points])
        for j in range(2*num_stress_eval_points):
            sigma[j, i] = sigma_vec[2*i*num_stress_eval_points + j]

    #sigma = np.reshape(sigma, (2*num_stress_eval_points, nelems))

    # Plot the stress distribution of each eval point along the span
    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    ax.plot(x, sigma_max, color="C0")
    for i in range(2*num_stress_eval_points):
        ax.plot(x, sigma[i, :], color="C0", alpha=0.2)

    ax.set_xlabel("Span (inches)")
    ax.set_ylabel("Stress")
    plt.savefig(f"stress_distribution_theta_{int(theta*180.0/np.pi)}.pdf")

    return

def check_theta_area_props():

    theta_vals = (np.pi/180)*np.linspace(0.0, 90.0, 19)
    A_vals = np.zeros(len(theta_vals))
    Iyy_vals = np.zeros(len(theta_vals))
    Izz_vals = np.zeros(len(theta_vals))
    Iyz_vals = np.zeros(len(theta_vals))
    _, p, _, _ = get_1d_problem(theta_val=0.0)

    for i in range(len(theta_vals)):
        p.set_val("theta_val", theta_vals[i])
        p.run_model()
        Iyy_vals[i] = p.get_val("structural_group.Iyy", units="inch**4")[0]
        Izz_vals[i] = p.get_val("structural_group.Izz", units="inch**4")[0]
        Iyz_vals[i] = p.get_val("structural_group.Iyz", units="inch**4")[0]

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    ax.plot((180.0/np.pi)*theta_vals, Iyy_vals, label=r"$I_{yy}$")
    ax.plot((180.0/np.pi)*theta_vals, Izz_vals, label=r"$I_{zz}$")
    ax.plot((180.0/np.pi)*theta_vals, Iyz_vals, label=r"$I_{yz}$")

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xlabel("Twist (deg.)")
    ax.set_ylabel(r"$I$ (inches$^4$)")
    ax.legend(frameon=False)

    plt.savefig("theta_sweep_area_props.pdf")

    return

def plot_chord_theta(xe, chord, theta):

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
    ax0.plot(xe, chord, color="tab:blue")
    ax1.plot(xe, theta, color="tab:blue")

    ax0.spines["right"].set_visible(False)
    ax0.spines["top"].set_visible(False)
    ax0.spines["bottom"].set_visible(False)
    ax0.tick_params(axis="x", direction="out", length=0.0, width=0.0)
    #ax0.xaxis.set_visible(False)
    ax0.yaxis.set_ticks_position("left")
    ax0.xaxis.set_ticks_position("bottom")
    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)
    ax0.grid(True)
    ax1.grid(True)

    # ticks_loc = ax1.get_xticks().tolist()
    # ax1.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # ax1.set_xticklabels([label_format.format(x) for x in ticks_loc])
    # ax0.set_xticklabels([])

    ax0.set_ylabel("Chord (in)")
    ax1.set_xlabel(r"$x_1$ (in)")
    ax1.set_ylabel("Twist (deg.)")

    plt.savefig("chord_theta.pdf", transparent=False)

    return

def plot_extras(p, xe, Tp, Np):

    nelems = len(Tp)

    # Plot everything
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, sharex=True, constrained_layout=True)

    E = 72.4e9
    Hyy = E*p.get_val("structural_group.Iyy")
    Hzz = E*p.get_val("structural_group.Izz")
    Hyz = E*p.get_val("structural_group.Iyz")
    Fx = p.get_val("structural_group.solver_comp.Fx")
    My = p.get_val("structural_group.solver_comp.My")
    Mz = p.get_val("structural_group.solver_comp.Mz")
    sigma_vm = p.get_val("sigma_vm")
    sigma_x1 = np.zeros(nelems)
    num_stress_eval_points = 20
    for i in range(nelems):
        sigma_x1[i] = np.amax(sigma_vm[2*i*num_stress_eval_points:2*(i+1)*num_stress_eval_points])

    ax0.plot(xe, Tp, label=r"$T_p$")
    ax0.plot(xe, Np, label=r"$N_p$")

    ax1.plot(xe, Hyy, label=r"$H_{yy}$")
    ax1.plot(xe, Hzz, label=r"$H_{zz}$")
    ax1.plot(xe, Hyz, label=r"$H_{yz}$")

    ax2r = ax2.twinx()
    ax2r.plot(xe, Fx, c="C0", label=r"$N_{x}$")
    ax2.plot(xe, My, c="C1", label=r"$M_{y}$")
    ax2.plot(xe, Mz, c="C2", label=r"$M_{z}$")

    ax3.plot(xe, sigma_x1, label=r"$\sigma_1$")

    for ax in [ax0, ax1, ax2]:
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", direction="out", length=0.0, width=0.0)
        #ax.xaxis.set_visible(False)
        ax.yaxis.set_ticks_position("left")
        ax.xaxis.set_ticks_position("bottom")
        ax.grid(True)

    ax2r.spines["top"].set_visible(False)
    ax2r.spines["bottom"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.spines["top"].set_visible(False)

    ax2r.grid(None)
    ax2.grid(True)
    ax3.grid(True)

    # ticks_loc = ax1.get_xticks().tolist()
    # ax1.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    # ax1.set_xticklabels([label_format.format(x) for x in ticks_loc])
    # ax0.set_xticklabels([])

    fs = 6
    ax0.set_ylabel("Aero forces (N/m)", fontsize=fs)
    ax1.set_ylabel(r"Bending stiffness (N-m$^2$)", fontsize=fs)
    ax2.set_ylabel("Bending moment (N-m)", fontsize=fs)
    ax2r.set_ylabel("Axial force (N)", fontsize=fs)
    ax3.set_ylabel(r"$\sigma_1(x_1)/\sigma_y$")
    ax3.set_xlabel(r"$x_1$ (in)")

    leg0 = ax0.legend()
    leg1 = ax1.legend()
    leg3 = ax1.legend()

    h1, l1 = ax2r.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    leg2 = ax2.legend([h1[0]]+h2, [l1[0]]+l2)

    plt.savefig("opt_vals.pdf", transparent=True)

    return

if __name__ == "__main__":

    #check_theta_opt_val()
    #check_chord_opt_val()
    # _, p, _, _ = get_1d_problem(chord_val=1.0)
    # p.run_driver()
    # print("chord = ", p.get_val("chord_val"))
    # print("theta = ", p.get_val("theta_val"))
    # print("KS = ", p.get_val("ks.KS"))

    #plot_stress_distribution(theta=89.0*np.pi/180.0)

    # xe, p, Np, Tp = get_1d_problem()
    # p.run_model()

    run_optimization(use_splines=False, optimizer="SNOPT")