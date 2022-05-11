import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import openmdao.api as om
from julia.OpenMDAO import make_component

from ccblade_openmdao_examples.ccblade_openmdao_component import BEMTRotorCAComp
from ccblade_openmdao_examples.gxbeam_openmdao_component import DiffBSplineComp
from ccblade_openmdao_examples.structural_group import StructuralGroup


def get_problem(sf=1.0, use_curvature_constraint=False, use_theta_dv=True, use_omega_dv=False, nonlinear=False):

    B = 3  # Number of blades.
    D = 24.0*0.0254  # Diameter in inches.
    P = 16.0*0.0254  # Pitch in inches (used in the twist distribution).
    c = 1.5*0.0254   # Constant chord in inches.
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

    # Lower and upper limits on the chord design variable, inches.
    chord_lower = 0.1*0.0254
    chord_upper = 5.0*0.0254

    # Lower and upper limits on the twist design variable, radians.
    theta_lower = 5.0*np.pi/180.0
    theta_upper = 85.0*np.pi/180.0

    # Target thrust value in Newtons.
    thrust_target = 97.246

    # Lower and upper limits on omega, rad/sec.
    omega_lower = 0.25*omega
    omega_upper = speedofsound/Rtip

    # Get the initial blade geometry.
    num_cp = 8
    radii_cp0 = np.linspace(Rhub, Rtip, num_cp)
    chord_cp0 = c*np.ones(num_cp)
    theta_cp0 = np.arctan(P/(np.pi*D*radii_cp0/Rtip))

    num_operating_points = 1
    num_radial = 30
    num_stress_eval_points = 20
    nelems = num_radial

    # Define the nodal coordinates of the finite element mesh
    x0 = np.linspace(Rhub, Rtip, nelems+1)
    u1_dv = np.zeros(nelems+1)
    u2_dv = np.zeros(nelems+1)
    u3_dv = np.zeros(nelems+1)

    prob = om.Problem()

    ivc = om.IndepVarComp()
    ivc.add_output("Rhub", val=Rhub, units="m")
    ivc.add_output("Rtip", val=Rtip, units="m")
    ivc.add_output("radii_cp", val=radii_cp0, units="m")
    ivc.add_output("chord_cp", val=chord_cp0, units="m")
    ivc.add_output("theta_cp", val=theta_cp0, units="rad")
    ivc.add_output("v", val=v, shape=num_operating_points, units="m/s")
    ivc.add_output("omega", val=omega, shape=num_operating_points, units="rad/s")
    ivc.add_output("pitch", val=pitch, shape=num_operating_points, units="rad")
    ivc.add_output("collective", val=0.0, units="rad")
    ivc.add_output("x0", val=x0, units="m")
    ivc.add_output("u1_dv", val=u1_dv, units="m")
    ivc.add_output("u2_dv", val=u2_dv, units="m")
    ivc.add_output("u3_dv", val=u3_dv, units="m")
    prob.model.add_subsystem("ivc", ivc, promotes_outputs=["*"])

    prob.model.add_subsystem("theta_adder", om.ExecComp("total_theta_cp=theta_cp+collective", total_theta_cp=theta_cp0, theta_cp=theta_cp0, collective=0.0, units="rad"), promotes=["*"])
    prob.model.add_subsystem("x_adder", om.ExecComp("xvals=x0+u1_dv", xvals=x0, x0=x0, u1_dv=u1_dv, units="m"), promotes=["*"])

    spline_comp = make_component(DiffBSplineComp(ncp=num_cp, nelems=nelems))
    prob.model.add_subsystem("spline_comp", spline_comp,
                             promotes_inputs=["radii_cp", "chord_cp"],
                             promotes_outputs=["radii", "chord", "theta", "d2c_dr2", "d2t_dr2"])
    prob.model.connect("total_theta_cp", "spline_comp.theta_cp")

    af_fname = "../data/xf-n0012-il-500000.dat"
    comp = make_component(
        BEMTRotorCAComp(
            af_fname=af_fname, cr75=c/Rtip, Re_exp=0.6,
            num_operating_points=num_operating_points, num_blades=B,
            num_radial=num_radial, rho=rho0, mu=mu, speedofsound=speedofsound))
    prob.model.add_subsystem("bemt_rotor_comp", comp, promotes_inputs=["Rhub", "Rtip", "radii", "chord", "theta", "v", "omega", "pitch"], promotes_outputs=["thrust", "torque", "efficiency", "Tp", "Np"])

    struc_group = StructuralGroup(nelems=nelems, num_stress_eval_points=num_stress_eval_points)
    prob.model.add_subsystem("structural_group", struc_group,
                             promotes_inputs=["omega", "Tp", "Np", "chord", "theta"],
                             promotes_outputs=["sigma_vm", "m", "u1", "u2", "u3"])

    prob.model.connect("xvals", "structural_group.xvals")
    prob.model.connect("u2_dv", "structural_group.yvals")
    prob.model.connect("u3_dv", "structural_group.zvals")
    if nonlinear:
        prob.model.add_subsystem("u1_con", om.ExecComp("u1_eq=u1_dv-u1", u1_eq=u1_dv, u1_dv=u1_dv, u1=np.zeros(nelems+1), units="m"), promotes=["*"])
        prob.model.add_subsystem("u2_con", om.ExecComp("u2_eq=u2_dv-u2", u2_eq=u2_dv, u2_dv=u2_dv, u2=np.zeros(nelems+1), units="m"), promotes=["*"])
        prob.model.add_subsystem("u3_con", om.ExecComp("u3_eq=u3_dv-u3", u3_eq=u3_dv, u3_dv=u3_dv, u3=np.zeros(nelems+1), units="m"), promotes=["*"])

    # Set the optimizer
    prob.model.linear_solver = om.DirectSolver()
    prob.driver = om.pyOptSparseDriver(optimizer="SNOPT")
    #prob.driver = om.ScipyOptimizeDriver(maxiter=200)
    #prob.driver.options["optimizer"] = "SLSQP"

    # Define the optimization problem
    prob.model.add_design_var("chord_cp", ref=1e-2)#, lower=chord_lower, upper=chord_upper, ref=1e-2)
    prob.model.add_constraint("chord", lower=chord_lower, upper=chord_upper, ref=1e-2)

    if use_theta_dv:
        prob.model.add_design_var("theta_cp")#, lower=theta_lower, upper=theta_upper, ref=1e0)
        prob.model.add_constraint("theta", lower=theta_lower, upper=theta_upper, ref=1e0)
    else:
        prob.model.add_design_var("collective", lower=theta_lower-np.amin(theta_cp0), upper=theta_upper-np.amax(theta_cp0), ref=1e0)

    if use_omega_dv:
        prob.model.add_design_var("omega", lower=omega_lower, upper=omega_upper, ref=omega)

    prob.model.add_objective("efficiency", ref=-1e0)
    prob.model.add_constraint("thrust", lower=thrust_target, upper=thrust_target, units="N", ref=1e2)

    if sf > 0.0:
        prob.model.add_constraint("sigma_vm", upper=1.0/sf)

    if use_curvature_constraint:
        prob.model.add_constraint("d2c_dr2", upper=0.0)
        if use_theta_dv:
            prob.model.add_constraint("d2t_dr2", lower=0.0)

    if nonlinear and sf > 0.0:
        prob.model.add_design_var("u1_dv", ref=1e-3)
        prob.model.add_design_var("u2_dv", ref=1e-3)
        prob.model.add_design_var("u3_dv", ref=1e-3)
        prob.model.add_constraint("u1_eq", lower=0.0, upper=0.0, ref=1e-3)
        prob.model.add_constraint("u2_eq", lower=0.0, upper=0.0, ref=1e-3)
        prob.model.add_constraint("u3_eq", lower=0.0, upper=0.0, ref=1e-3)

    prob.setup(check=True)
    om.n2(prob, show_browser=False, outfile='struc_aero_opt.html')

    return prob

def run_optimization(sf=1.0, use_curvature_constraint=False, use_theta_dv=True, use_omega_dv=False, nonlinear=False):
    # Run the coupled aero-structural optimization problem and plot the outputs

    if use_omega_dv and use_theta_dv:
        dv_fname = "dvs_using_theta_omega_sf_{0:.1f}.csv".format(sf)
        dv_cp_fname = "dvs_cp_using_theta_omega_sf_{0:.1f}.csv".format(sf)
        force_fname = "aero_forces_using_theta_omega_sf_{0:.1f}.csv".format(sf)
        stress_fname = "stress_using_theta_omega_sf_{0:.1f}.csv".format(sf)
        disp_fname = "disp_using_theta_omega_sf_{0:.1f}.csv".format(sf)
    elif use_omega_dv:
        dv_fname = "dvs_using_omega_sf_{0:.1f}.csv".format(sf)
        dv_cp_fname = "dvs_cp_using_omega_sf_{0:.1f}.csv".format(sf)
        force_fname = "aero_forces_using_omega_sf_{0:.1f}.csv".format(sf)
        stress_fname = "stress_using_omega_sf_{0:.1f}.csv".format(sf)
        disp_fname = "disp_using_omega_sf_{0:.1f}.csv".format(sf)
    elif use_theta_dv:
        dv_fname = "dvs_using_theta_sf_{0:.1f}.csv".format(sf)
        dv_cp_fname = "dvs_cp_using_theta_sf_{0:.1f}.csv".format(sf)
        force_fname = "aero_forces_using_theta_sf_{0:.1f}.csv".format(sf)
        stress_fname = "stress_using_theta_sf_{0:.1f}.csv".format(sf)
        disp_fname = "disp_using_theta_sf_{0:.1f}.csv".format(sf)
    else:
        dv_fname = "dvs_using_chord_sf_{0:.1f}.csv".format(sf)
        dv_cp_fname = "dvs_cp_using_chord_sf_{0:.1f}.csv".format(sf)
        force_fname = "aero_forces_using_chord_sf_{0:.1f}.csv".format(sf)
        stress_fname = "stress_using_chord_sf_{0:.1f}.csv".format(sf)
        disp_fname = "disp_using_chord_sf_{0:.1f}.csv".format(sf)

    if nonlinear:
        dv_fname = "nl_" + dv_fname
        dv_cp_fname = "nl_" + dv_cp_fname
        force_fname = "nl_" + force_fname
        stress_fname = "nl_" + stress_fname
        disp_fname = "nl_" + disp_fname

    p = get_problem(sf=sf, use_curvature_constraint=use_curvature_constraint, use_theta_dv=use_theta_dv, use_omega_dv=use_omega_dv, nonlinear=nonlinear)
    p.run_driver()

    xe = p.get_val("radii", units="inch")
    x_cp = p.get_val("radii_cp", units="inch")
    print("mass = ", p.get_val("m", units="kg"))
    print("max(sigma) = ", np.amax(p.get_val("sigma_vm")))
    print("efficiency = ", p.get_val("efficiency"))
    print("omega = ", p.get_val("omega"))
    print("collective = ", p.get_val("collective")[0]*(180.0/np.pi))
    print("thrust_eq rel = ", (p.get_val("thrust")[0]-97.246)/97.246)

    # Save the chord and twist distribution to a csv
    chord = p.get_val("chord", units="inch")
    theta = p.get_val("theta", units="deg")
    chord_cp = p.get_val("chord_cp", units="inch")
    theta_cp = p.get_val("theta_cp", units="deg")
    df = pd.DataFrame({"radii":xe,
                       "chord":chord,
                       "theta":theta})
    df.to_csv(dv_fname, index=False)
    df = pd.DataFrame({"radii_cp":x_cp,
                       "chord_cp":chord_cp,
                       "theta_cp":theta_cp})
    df.to_csv(dv_cp_fname, index=False)

    # Plot the chord and twist distribution
    plot_chord_theta(xe, chord, theta)

    # Plot the other values of interest
    Tp = p.get_val("Tp")[0]
    Np = p.get_val("Np")[0]
    plot_extras(p, xe, Tp, Np)
    df = pd.DataFrame({"Np":Np, "Tp":Tp})
    df.to_csv(force_fname, index=False)

    # Write out the element-wise max stress to a csv file
    sigma1 = p.get_val("sigma_vm")
    sigma_x1 = np.zeros(len(xe))
    num_stress_eval_points = 20
    for i in range(len(xe)):
        sigma_x1[i] = np.amax(sigma1[2*i*num_stress_eval_points:2*(i+1)*num_stress_eval_points])
    df = pd.DataFrame({"sigma":sigma_x1})
    df.to_csv(stress_fname, index=False)

    # Write out the displacement to a csv file
    x0 = p.get_val("x0", units="inch")
    u1 = p.get_val("u1", units="inch")
    u2 = p.get_val("u2", units="inch")
    u3 = p.get_val("u3", units="inch")
    u1_dv = p.get_val("u1_dv", units="inch")
    u2_dv = p.get_val("u2_dv", units="inch")
    u3_dv = p.get_val("u3_dv", units="inch")
    # u1e = p.get_val("structural_group.u1e", units="inch")
    # u2e = p.get_val("structural_group.u2e", units="inch")
    # u3e = p.get_val("structural_group.u3e", units="inch")
    # ue = np.sqrt((u1e**2) + (u2e**2) + (u3e**2))
    df = pd.DataFrame({"x0":x0, "u1":u1, "u2":u2, "u3":u3,
                       "u1_dv":u1_dv, "u2_dv":u2_dv, "u3_dv":u3_dv})
    df.to_csv(disp_fname, index=False)

    # Print the nonlinear displacement convergence
    if nonlinear:
        print("max rel u1 diff = ", np.amax(np.abs(u1_dv-u1)/(np.abs(u1)+1e-10)))
        print("max rel u2 diff = ", np.amax(np.abs(u2_dv-u2)/(np.abs(u2)+1e-10)))
        print("max rel u3 diff = ", np.amax(np.abs(u3_dv-u3)/(np.abs(u3)+1e-10)))

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
    ax1.set_xlabel(r"$x$ (in)")
    ax1.set_ylabel("Twist (deg.)")

    plt.savefig("chord_theta.pdf", transparent=False)

    return

def plot_extras(p, xe, Tp, Np):

    nelems = len(Tp)

    # Plot everything
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, sharex=True, constrained_layout=True)

    E = 72.4e9
    Iyy = E*p.get_val("structural_group.Iyy", units="inch**4")
    Izz = E*p.get_val("structural_group.Izz", units="inch**4")
    Iyz = E*p.get_val("structural_group.Iyz", units="inch**4")
    Fx = p.get_val("structural_group.solver_comp.Fx")
    My = p.get_val("structural_group.solver_comp.My")
    Mz = p.get_val("structural_group.solver_comp.Mz")
    sigma1 = p.get_val("sigma_vm")
    sigma_x1 = np.zeros(nelems)
    num_stress_eval_points = 20
    for i in range(nelems):
        sigma_x1[i] = np.amax(sigma1[2*i*num_stress_eval_points:2*(i+1)*num_stress_eval_points])

    ax0.plot(xe, Tp, label=r"$T_p$")
    ax0.plot(xe, Np, label=r"$N_p$")

    ax1.plot(xe, Iyy, label=r"$I_{yy}$")
    ax1.plot(xe, Izz, label=r"$I_{zz}$")
    ax1.plot(xe, Iyz, label=r"$I_{yz}$")

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
    ax1.set_ylabel(r"Moment of inertia (in$^4$)", fontsize=fs)
    ax2.set_ylabel("Bending moment (N-m)", fontsize=fs)
    ax2r.set_ylabel("Axial force (N)", fontsize=fs)
    ax3.set_ylabel(r"$\sigma_vM(x)/\sigma_y$")
    ax3.set_xlabel(r"$x$ (in)")

    leg0 = ax0.legend()
    leg1 = ax1.legend()
    leg3 = ax1.legend()

    h1, l1 = ax2r.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    leg2 = ax2.legend([h1[0]]+h2, [l1[0]]+l2)

    plt.savefig("opt_vals.pdf", transparent=True)

    return

if __name__ == "__main__":

    run_optimization(sf=4.0, use_curvature_constraint=True, use_theta_dv=True, use_omega_dv=True, nonlinear=True)