import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import openmdao.api as om
from openmdao.devtools import iprofile
from julia.OpenMDAO import make_component
from julia import CCBladeLoadingExample
from paropt.paropt_driver import ParOptDriver

from ccblade_openmdao_examples.structural_group import StructuralGroup
from ccblade_openmdao_examples.gxbeam_openmdao_component import MassComp


def get_problem():

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
    twist = (45.0*np.pi/180.0)*np.ones(nelems)  # (rad)

    # Define the angular rotation
    rpm = 7110.0
    omega = rpm*2*np.pi/60

    prob = om.Problem()

    ivc = om.IndepVarComp()
    ivc.add_output("omega", omega, units="rad/s")
    ivc.add_output("Tp", Tp, units="N/m")
    ivc.add_output("Np", Np, units="N/m")
    ivc.add_output("chord", chord, units="inch")
    ivc.add_output("twist", twist, units="rad")
    prob.model.add_subsystem("ivc", ivc, promotes=["*"])

    struc_group = StructuralGroup(nelems=nelems, num_stress_eval_points=num_stress_eval_points)
    prob.model.add_subsystem("structural_group", struc_group,
                             promotes_inputs=["omega", "Tp", "Np", "chord", "twist"],
                             promotes_outputs=["sigma1", "m"])

    # Aggregate the stress
    prob.model.add_subsystem("ks", om.KSComp(width=nelems*num_stress_eval_points*2,
                                             add_constraint=False, ref=1.0,
                                             units=None))
    prob.model.connect("sigma1", "ks.g")

    #prob.driver = ParOptDriver()

    # # Set the optimization options
    # options = {"algorithm": "tr",
    #            "tr_max_size": 0.1,
    #            #"tr_min_size": 1e-4,
    #            "tr_adaptive_gamma_update": True,
    #            "qn_diag_type": "yts_over_sts",
    #            "tr_use_soc": True,
    #            "tr_max_iterations": 1000,
    #            "penalty_gamma": 5.0,
    #            "tr_penalty_gamma_max": 10.0,
    #            "qn_subspace_size": 20,
    #            "qn_type": "bfgs",
    #            "abs_res_tol": 1e-8,
    #            "starting_point_strategy": "affine_step",
    #            "barrier_strategy": "mehrotra_predictor_corrector",
    #            "use_line_search": False,
    #            "max_major_iters": 200}
    # for key in options:
    #     prob.driver.options[key] = options[key]

    prob.driver = om.pyOptSparseDriver(optimizer="SNOPT")

    # Lower and upper limits on the chord design variable, in inches.
    chord_lower = 0.5
    chord_upper = 2.0

    # Lower and upper limits on the twist design variable, radians.
    theta_lower =  0.0*np.pi/180.0
    theta_upper =  90.0*np.pi/180.0

    prob.model.add_design_var("chord", lower=chord_lower, upper=chord_upper, ref=1e0, units="inch")
    prob.model.add_design_var("twist", lower=theta_lower, upper=theta_upper, ref=1e0, units="rad")

    # Stress-constrained mass minimization
    prob.model.add_objective("m", ref=1e-2)
    prob.model.add_constraint("ks.KS", upper=1.0, ref=1e0)

    prob.setup()

    return x, prob, Np, Tp

def run_optimization():
    # Run the structural optimization problem and plot the outputs

    xe, p, Np, Tp = get_problem()
    p.run_driver()

    print("mass = ", p.get_val("m", units="kg"))
    print("KS(sigma1)/sigma_y = ", p.get_val("ks.KS"))
    print("max(sigma1) = ", np.amax(p.get_val("sigma1")))

    # Save the chord and twist distribution to a csv
    chord = p.get_val("chord", units="inch")
    theta = p.get_val("twist", units="deg")
    df = pd.DataFrame({"chord":chord, "theta":theta})
    df.to_csv("chord_theta.csv", index=False)

    # Plot the chord and twist distribution
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

def get_1d_problem(chord_val=None, twist_val=None):

    # If chord_val is None, treat chord as a design variable
    # If twist_val is None, treat twist as a design variable

    # Compute the aerodynamic loads
    x, Np, Tp = CCBladeLoadingExample.run_ccblade()
    nelems = len(x)
    num_stress_eval_points = 20

    # Initialize the design variables
    if not chord_val:
        init_chord_val = 1.0
    else:
        init_chord_val = chord_val

    if not twist_val:
        init_twist_val = 45.0*np.pi/180.0  # (rad)
    else:
        init_twist_val = twist_val

    # Define the angular rotation
    rpm = 7110.0
    omega = rpm*2*np.pi/60

    prob = om.Problem()

    ivc = om.IndepVarComp()
    ivc.add_output("omega", omega, units="rad/s")
    ivc.add_output("Tp", Tp, units="N/m")
    ivc.add_output("Np", Np, units="N/m")
    ivc.add_output("chord_val", init_chord_val, units="inch")
    ivc.add_output("twist_val", init_twist_val, units="rad")
    prob.model.add_subsystem("ivc", ivc, promotes=["*"])

    chord_comp = ConvertScalarToVec(nelems=nelems, units="inch")
    prob.model.add_subsystem("chord_comp", chord_comp)
    prob.model.connect("chord_val", "chord_comp.val")

    twist_comp = ConvertScalarToVec(nelems=nelems, units="rad")
    prob.model.add_subsystem("twist_comp", twist_comp)
    prob.model.connect("twist_val", "twist_comp.val")

    struc_group = StructuralGroup(nelems=nelems, num_stress_eval_points=num_stress_eval_points)
    prob.model.add_subsystem("structural_group", struc_group,
                             promotes_inputs=["omega", "Tp", "Np", "chord", "twist"],
                             promotes_outputs=["sigma1", "m"])

    prob.model.connect("chord_comp.val_vec", "chord")
    prob.model.connect("twist_comp.val_vec", "twist")

    # Aggregate the stress
    prob.model.add_subsystem("ks", om.KSComp(width=nelems*num_stress_eval_points*2,
                                             add_constraint=False, ref=1.0,
                                             units=None))
    prob.model.connect("sigma1", "ks.g")

#     prob.driver = ParOptDriver()
#
#     # Set the optimization options
#     options = {"algorithm": "tr",
#                "tr_max_size": 0.1,
#                #"tr_min_size": 1e-4,
#                "tr_adaptive_gamma_update": True,
#                "qn_diag_type": "yts_over_sts",
#                "tr_use_soc": True,
#                "tr_max_iterations": 1000,
#                "penalty_gamma": 5.0,
#                "tr_penalty_gamma_max": 10.0,
#                "qn_subspace_size": 20,
#                "qn_type": "bfgs",
#                "abs_res_tol": 1e-8,
#                "starting_point_strategy": "affine_step",
#                "barrier_strategy": "mehrotra_predictor_corrector",
#                "use_line_search": False,
#                "max_major_iters": 200}
#     for key in options:
#         prob.driver.options[key] = options[key]

    prob.driver = om.pyOptSparseDriver(optimizer="SNOPT")

    # Lower and upper limits on the chord design variable, in inches.
    chord_lower = 0.5
    chord_upper = 2.0

    # Lower and upper limits on the twist design variable, radians.
    theta_lower =  0.0*np.pi/180.0
    theta_upper =  90.0*np.pi/180.0

    if chord_val is None:
        prob.model.add_design_var("chord_val", lower=chord_lower, upper=chord_upper, ref=1e0, units="inch")
    if twist_val is None:
        prob.model.add_design_var("twist_val", lower=theta_lower, upper=theta_upper, ref=1e0, units="rad")

    # Stress minimization
    prob.model.add_objective("ks.KS", ref=1e0)

    prob.setup()
    om.n2(prob, show_browser=False, outfile='struc_1d.html')

    return x, prob, Np, Tp

# def run_1d_optimization():
#
#     return
#
# def run_2d_optimization():
#
#     return

def check_twist_opt_val():
    # Run a sweep of uniform twist values and plot the KS stress at each point
    # Then solve an optimization with uniform twist as the only design variable and check that point

    theta_vals = (np.pi/180)*np.linspace(0.0, 90.0, 19)
    ks_vals = np.zeros(len(theta_vals))
    _, p, _, _ = get_1d_problem(twist_val=0.0)

    for i in range(len(theta_vals)):
        p.set_val("twist_val", theta_vals[i])
        p.run_model()
        ks_vals[i] = p.get_val("ks.KS")

    # Reset the starting twist value for optimization
    # _, p, _, _ = get_1d_problem(chord_val=1.0)
    # p.run_driver()
    # theta_opt = p.get_val("twist_val")
    # ks_opt = p.get_val("ks.KS")

    # Hard-code the optimal
    theta_opt = 1.198250*180.0/np.pi
    ks_opt = 0.9297127

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    ax.plot((180.0/np.pi)*theta_vals, ks_vals, label="Value sweep")
    ax.scatter(theta_opt, ks_opt, color="tab:red", label="Optimal", zorder=100)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xlabel("Twist (deg.)")
    ax.set_ylabel(r"KS($\sigma$)")
    ax.legend(frameon=False)
    plt.savefig("twist_opt_sweep.pdf")

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
    # _, p, _, _ = get_1d_problem(twist_val=0.0)
    # p.run_driver()
    # chord_opt = p.get_val("chord_val")
    # ks_opt = p.get_val("ks.KS")
    # print("Chord = ", chord_opt)
    # print("KS(sigma) = ", ks_opt)
    chord_opt = 1.245781
    ks_opt = 0.9851451

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

def check_twist_area_props():

    theta_vals = (np.pi/180)*np.linspace(0.0, 90.0, 19)
    A_vals = np.zeros(len(theta_vals))
    Iyy_vals = np.zeros(len(theta_vals))
    Izz_vals = np.zeros(len(theta_vals))
    Iyz_vals = np.zeros(len(theta_vals))
    _, p, _, _ = get_1d_problem(twist_val=0.0)

    for i in range(len(theta_vals)):
        p.set_val("twist_val", theta_vals[i])
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

    plt.savefig("twist_sweep_area_props.pdf")

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

    plt.savefig("chord_theta.png", transparent=True)

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
    sigma1 = p.get_val("sigma1")
    sigma_x1 = np.zeros(nelems)
    num_stress_eval_points = 20
    for i in range(nelems):
        sigma_x1[i] = np.amax(sigma1[2*i*num_stress_eval_points:2*(i+1)*num_stress_eval_points])

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

    #check_twist_opt_val()
    #check_chord_opt_val()
    # _, p, _, _ = get_1d_problem(chord_val=1.0)
    # p.run_driver()
    # print(p.get_val("twist_val"))
    # print(p.get_val("ks.KS"))

    # xe, p, Np, Tp = get_1d_problem()
    # p.run_model()

    run_optimization()