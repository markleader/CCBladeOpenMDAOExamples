using BSplineKit: BSplineKit
using FLOWMath: ksmax
using ComponentArrays
using ConcreteStructs
using ForwardDiff
using OpenMDAO: AbstractExplicitComp, VarData, PartialsData

@concrete struct DiffBSplineComp <: AbstractExplicitComp
    ncp
    nelems
    apply_nonlinear_forwarddiffable!
    x
    y
    J
    forwarddiff_config
end

function DiffBSplineComp(; ncp, nelems)

    function apply_nonlinear_forwarddiffable!(y, x)
        T = eltype(x)

        # Unpack the inputs.
        r_cp = x[:r_cp]
        chord_cp = x[:chord_cp]
        theta_cp = x[:theta_cp]

        x_cp = collect(range(0.0, 1.0, length=ncp))
        dx = 1.0/nelems
        x_interp = collect(range(dx/2, 1.0-dx/2, length=nelems))

        radii_interp = BSplineKit.interpolate(x_cp, r_cp, BSplineKit.BSplineOrder(4))

        chord_interp = BSplineKit.interpolate(x_cp, chord_cp, BSplineKit.BSplineOrder(4))
        d2c_dx2_func = BSplineKit.diff(chord_interp, BSplineKit.Derivative(2))

        theta_interp = BSplineKit.interpolate(x_cp, theta_cp, BSplineKit.BSplineOrder(4))
        d2t_dx2_func = BSplineKit.diff(theta_interp, BSplineKit.Derivative(2))

        y[:r] .= radii_interp.(x_interp)
        y[:chord] .= chord_interp.(x_interp)
        y[:theta] .= theta_interp.(x_interp)
        y[:d2c_dr2] .= d2c_dx2_func.(x_interp)
        y[:d2t_dr2] .= d2t_dx2_func.(x_interp)

        return nothing
    end

    # Initialize the input and output vectors needed by ForwardDiff.jl.
    X = ComponentArray(
        r_cp=zeros(Float64, ncp), chord_cp=zeros(Float64, ncp),
        theta_cp=zeros(Float64, ncp))
    Y = ComponentArray(
        r=zeros(Float64, nelems), chord=zeros(Float64, nelems), theta=zeros(Float64, nelems),
        d2c_dr2=zeros(Float64, nelems), d2t_dr2=zeros(Float64, nelems))
    J = Y.*X'

    # Get the JacobianConfig object, which we'll reuse each time when calling
    # the ForwardDiff.jacobian! function (apparently good for efficiency).
    config = ForwardDiff.JacobianConfig(apply_nonlinear_forwarddiffable!, Y, X)

    return DiffBSplineComp(ncp, nelems, apply_nonlinear_forwarddiffable!, X, Y, J, config)
end

function OpenMDAO.setup(self::DiffBSplineComp)

    # Declare the OpenMDAO inputs.
    input_data = Vector{VarData}()
    push!(input_data, VarData("radii_cp", shape=self.ncp, val=collect(range(2.4*.0254, 12.0*.0254, length=self.ncp)), units="m"))
    push!(input_data, VarData("chord_cp", shape=self.ncp, units="m"))
    push!(input_data, VarData("theta_cp", shape=self.ncp, units="rad"))

    # Declare the OpenMDAO outputs.
    output_data = Vector{VarData}()
    push!(output_data, VarData("radii", shape=self.nelems, units="m"))
    push!(output_data, VarData("chord", shape=self.nelems, units="m"))
    push!(output_data, VarData("theta", shape=self.nelems, units="rad"))
    push!(output_data, VarData("d2c_dr2", shape=self.nelems))
    push!(output_data, VarData("d2t_dr2", shape=self.nelems))

    # Declare the OpenMDAO partial derivatives.
    partials_data = Vector{PartialsData}()
    push!(partials_data, PartialsData("radii", "radii_cp"))
    push!(partials_data, PartialsData("chord", "chord_cp"))
    push!(partials_data, PartialsData("theta", "theta_cp"))
    push!(partials_data, PartialsData("d2c_dr2", "chord_cp"))
    push!(partials_data, PartialsData("d2t_dr2", "theta_cp"))

    return input_data, output_data, partials_data
end

function OpenMDAO.compute!(self::DiffBSplineComp, inputs, outputs)

    # Unpack the inputs.
    r_cp = inputs["radii_cp"]
    chord_cp = inputs["chord_cp"]
    theta_cp = inputs["theta_cp"]

    # Unpack the outputs.
    r = outputs["radii"]
    chord = outputs["chord"]
    theta = outputs["theta"]
    d2c_dx2 = outputs["d2c_dr2"]
    d2t_dx2 = outputs["d2t_dr2"]

    x_cp = collect(range(0.0, 1.0, length=self.ncp))
    dx = 1.0/self.nelems
    x_interp = collect(range(dx/2, 1.0-dx/2, length=self.nelems))

    radii_interp = BSplineKit.interpolate(x_cp, r_cp, BSplineKit.BSplineOrder(4))
    r .= radii_interp.(x_interp)

    chord_interp = BSplineKit.interpolate(x_cp, chord_cp, BSplineKit.BSplineOrder(4))
    chord .= chord_interp.(x_interp)
    d2c_dx2_func = BSplineKit.diff(chord_interp, BSplineKit.Derivative(2))
    d2c_dx2 .= d2c_dx2_func.(x_interp)

    theta_interp = BSplineKit.interpolate(x_cp, theta_cp, BSplineKit.BSplineOrder(4))
    theta .= theta_interp.(x_interp)
    d2t_dx2_func = BSplineKit.diff(theta_interp, BSplineKit.Derivative(2))
    d2t_dx2 .= d2t_dx2_func.(x_interp)

end

function OpenMDAO.compute_partials!(self::DiffBSplineComp, inputs, partials)

    # Unpack the inputs.
    r_cp = inputs["radii_cp"]
    chord_cp = inputs["chord_cp"]
    theta_cp = inputs["theta_cp"]

    # Working arrays and configuration for ForwardDiff's Jacobian routine.
    x = self.x
    y = self.y
    J = self.J
    config = self.forwarddiff_config

    x[:r_cp] .= r_cp
    x[:chord_cp] .= chord_cp
    x[:theta_cp] .= theta_cp

    dr_dr_cp = partials["radii", "radii_cp"]
    dchord_dchord_cp = partials["chord", "chord_cp"]
    dd2c_dr2_dchord_cp = partials["d2c_dr2", "chord_cp"]
    dtheta_dtheta_cp = partials["theta", "theta_cp"]
    dd2t_dr2_dtheta_cp = partials["d2t_dr2", "theta_cp"]

    # Get the Jacobian.
    ForwardDiff.jacobian!(J, self.apply_nonlinear_forwarddiffable!, y, x, config)

    dr_dr_cp .= J[:r, :r_cp]
    dchord_dchord_cp .= J[:chord, :chord_cp]
    dd2c_dr2_dchord_cp .= J[:d2c_dr2, :chord_cp]
    dtheta_dtheta_cp .= J[:theta, :theta_cp]
    dd2t_dr2_dtheta_cp .= J[:d2t_dr2, :theta_cp]

end