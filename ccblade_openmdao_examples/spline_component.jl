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
        r = x[:r]

        chord_interp = BSplineKit.interpolate(r_cp, chord_cp, BSplineKit.BSplineOrder(4))
        d2c_dr2_func = BSplineKit.diff(chord_interp, BSplineKit.Derivative(2))

        theta_interp = BSplineKit.interpolate(r_cp, theta_cp, BSplineKit.BSplineOrder(4))
        d2t_dr2_func = BSplineKit.diff(theta_interp, BSplineKit.Derivative(2))

        y[:chord] .= chord_interp.(r)
        y[:theta] .= theta_interp.(r)
        y[:d2c_dr2] .= ksmax(d2c_dr2_func.(r))
        y[:d2t_dr2] .= ksmax(d2t_dr2_func.(r))

        return nothing
    end

    # Initialize the input and output vectors needed by ForwardDiff.jl.
    X = ComponentArray(
        r_cp=zeros(Float64, ncp), chord_cp=zeros(Float64, ncp),
        theta_cp=zeros(Float64, ncp), r=zeros(Float64, nelems))
    Y = ComponentArray(
        chord=zeros(Float64, nelems), theta=zeros(Float64, nelems),
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
    push!(input_data, VarData("radii_cp", shape=self.ncp, val=collect(range(2.4, 12.0, length=self.ncp)), units="inch"))
    push!(input_data, VarData("chord_cp", shape=self.ncp, units="inch"))
    push!(input_data, VarData("theta_cp", shape=self.ncp, units="rad"))
    push!(input_data, VarData("radii", shape=self.nelems, val=collect(range(2.4, 12.0, length=self.nelems)), units="inch"))

    # Declare the OpenMDAO outputs.
    output_data = Vector{VarData}()
    push!(output_data, VarData("chord", shape=self.nelems, units="inch"))
    push!(output_data, VarData("theta", shape=self.nelems, units="rad"))
    push!(output_data, VarData("d2c_dr2", shape=1))#self.nelems))
    push!(output_data, VarData("d2t_dr2", shape=1))#self.nelems))

    # Declare the OpenMDAO partial derivatives.
    partials_data = Vector{PartialsData}()
    push!(partials_data, PartialsData("chord", "chord_cp"))
    push!(partials_data, PartialsData("chord", "radii_cp"))
    push!(partials_data, PartialsData("chord", "radii"))

    push!(partials_data, PartialsData("theta", "theta_cp"))
    push!(partials_data, PartialsData("theta", "radii_cp"))
    push!(partials_data, PartialsData("theta", "radii"))

    push!(partials_data, PartialsData("d2c_dr2", "chord_cp"))
    push!(partials_data, PartialsData("d2c_dr2", "radii_cp"))
    push!(partials_data, PartialsData("d2c_dr2", "radii"))

    push!(partials_data, PartialsData("d2t_dr2", "theta_cp"))
    push!(partials_data, PartialsData("d2t_dr2", "radii_cp"))
    push!(partials_data, PartialsData("d2t_dr2", "radii"))


    return input_data, output_data, partials_data
end

function OpenMDAO.compute!(self::DiffBSplineComp, inputs, outputs)

    # Unpack the inputs.
    r_cp = inputs["radii_cp"]
    chord_cp = inputs["chord_cp"]
    theta_cp = inputs["theta_cp"]
    r = inputs["radii"]

    # Unpack the outputs.
    chord = outputs["chord"]
    theta = outputs["theta"]
    d2c_dr2 = outputs["d2c_dr2"]
    d2t_dr2 = outputs["d2t_dr2"]

    chord_interp = BSplineKit.interpolate(r_cp, chord_cp, BSplineKit.BSplineOrder(4))
    chord .= chord_interp.(r)
    d2c_dr2_func = BSplineKit.diff(chord_interp, BSplineKit.Derivative(2))
    d2c_dr2 .= ksmax(d2c_dr2_func.(r))

    theta_interp = BSplineKit.interpolate(r_cp, theta_cp, BSplineKit.BSplineOrder(4))
    theta .= theta_interp.(r)
    d2t_dr2_func = BSplineKit.diff(theta_interp, BSplineKit.Derivative(2))
    d2t_dr2 .= ksmax(d2t_dr2_func.(r))

end

function OpenMDAO.compute_partials!(self::DiffBSplineComp, inputs, partials)

    # Unpack the inputs.
    r_cp = inputs["radii_cp"]
    chord_cp = inputs["chord_cp"]
    theta_cp = inputs["theta_cp"]
    r = inputs["radii"]

    # Working arrays and configuration for ForwardDiff's Jacobian routine.
    x = self.x
    y = self.y
    J = self.J
    config = self.forwarddiff_config

    x[:r_cp] .= r_cp
    x[:chord_cp] .= chord_cp
    x[:theta_cp] .= theta_cp
    x[:r] .= r

    dchord_dchord_cp = partials["chord", "chord_cp"]
    dchord_dr_cp = partials["chord", "radii_cp"]
    dchord_dr = partials["chord", "radii"]

    dd2c_dr2_dchord_cp = partials["d2c_dr2", "chord_cp"]
    dd2c_dr2_dr_cp = partials["d2c_dr2", "radii_cp"]
    dd2c_dr2_dr = partials["d2c_dr2", "radii"]

    dtheta_dtheta_cp = partials["theta", "theta_cp"]
    dtheta_dr_cp = partials["theta", "radii_cp"]
    dtheta_dr = partials["theta", "radii"]

    dd2t_dr2_dtheta_cp = partials["d2t_dr2", "theta_cp"]
    dd2t_dr2_dr_cp = partials["d2t_dr2", "radii_cp"]
    dd2t_dr2_dr = partials["d2t_dr2", "radii"]

    # Get the Jacobian.
    ForwardDiff.jacobian!(J, self.apply_nonlinear_forwarddiffable!, y, x, config)

    dchord_dchord_cp .= J[:chord, :chord_cp]
    dchord_dr_cp .= J[:chord, :r_cp]
    dchord_dr .= J[:chord, :r]

    dd2c_dr2_dchord_cp .= J[:d2c_dr2, :chord_cp]
    dd2c_dr2_dr_cp .= J[:d2c_dr2, :r_cp]
    dd2c_dr2_dr .= J[:d2c_dr2, :r]

    dtheta_dtheta_cp .= J[:theta, :theta_cp]
    dtheta_dr_cp .= J[:theta, :r_cp]
    dtheta_dr .= J[:theta, :r]

    dd2t_dr2_dtheta_cp .= J[:d2t_dr2, :theta_cp]
    dd2t_dr2_dr_cp .= J[:d2t_dr2, :r_cp]
    dd2t_dr2_dr .= J[:d2t_dr2, :r]

end