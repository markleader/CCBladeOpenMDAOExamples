using ConcreteStructs
using OpenMDAO: AbstractExplicitComp, VarData, PartialsData, get_rows_cols

@concrete struct AreaComp <: AbstractExplicitComp
    nelems
    A_ref
    Iyy_ref
    Izz_ref
end

function AreaComp(; nelems, A_ref, Iyy_ref, Izz_ref)

    return AreaComp(nelems, A_ref, Iyy_ref, Izz_ref)
end

function OpenMDAO.setup(self::AreaComp)

    # Declare the OpenMDAO inputs.
    input_data = Vector{VarData}()
    push!(input_data, VarData("chord", shape=self.nelems, units="m"))
    push!(input_data, VarData("theta", shape=self.nelems, units="rad"))

    # Declare the OpenMDAO outputs.
    output_data = Vector{VarData}()
    push!(output_data, VarData("A", shape=self.nelems, units="m**2"))
    push!(output_data, VarData("Iyy", shape=self.nelems, units="m**4"))
    push!(output_data, VarData("Izz", shape=self.nelems, units="m**4"))
    push!(output_data, VarData("Iyz", shape=self.nelems, units="m**4"))

    # Declare the OpenMDAO partial derivatives.
    partials_data = Vector{PartialsData}()
    rows = collect(range(0, self.nelems-1, length=self.nelems))
    cols = rows
    push!(partials_data, PartialsData("A", "chord", rows=rows, cols=cols))
    push!(partials_data, PartialsData("Iyy", "chord", rows=rows, cols=cols))
    push!(partials_data, PartialsData("Iyy", "theta", rows=rows, cols=cols))
    push!(partials_data, PartialsData("Izz", "chord", rows=rows, cols=cols))
    push!(partials_data, PartialsData("Izz", "theta", rows=rows, cols=cols))
    push!(partials_data, PartialsData("Iyz", "chord", rows=rows, cols=cols))
    push!(partials_data, PartialsData("Iyz", "theta", rows=rows, cols=cols))

    return input_data, output_data, partials_data
end

function OpenMDAO.compute!(self::AreaComp, inputs, outputs)

    # Unpack the inputs.
    chord = inputs["chord"]
    theta = inputs["theta"]

    # Unpack the outputs.
    A = outputs["A"]
    Iyy = outputs["Iyy"]
    Izz = outputs["Izz"]
    Iyz = outputs["Iyz"]

    s = sin.(theta)
    c = cos.(theta)
    s2 = sin.(theta).^2
    c2 = cos.(theta).^2

    k = chord/100  # scale factor
    A[1:end] = (k.^2)*self.A_ref
    Iyy0 = (k.^4)*self.Iyy_ref
    Izz0 = (k.^4)*self.Izz_ref
    Iyy[1:end] = @. c2 *Iyy0 + s2 *Izz0
    Izz[1:end] = @. s2 *Iyy0 + c2 *Izz0
    Iyz[1:end] = @. (Izz0 - Iyy0)*s*c

end

function OpenMDAO.compute_partials!(self::AreaComp, inputs, partials)

    # Unpack the inputs.
    chord = inputs["chord"]
    theta = inputs["theta"]

    dA_dc = partials["A", "chord"]
    dIyy_dc = partials["Iyy", "chord"]
    dIyy_dtheta = partials["Iyy", "theta"]
    dIzz_dc = partials["Izz", "chord"]
    dIzz_dtheta = partials["Izz", "theta"]
    dIyz_dc = partials["Iyz", "chord"]
    dIyz_dtheta = partials["Iyz", "theta"]

    s = sin.(theta)
    c = cos.(theta)
    s2 = s.^2
    c2 = c.^2

    k = chord/100  # scale factor
    A = (k.^2)*self.A_ref
    Iyy0 = (k.^4)*self.Iyy_ref
    Izz0 = (k.^4)*self.Izz_ref
    Iyy = @. c2 *Iyy0 + s2 *Izz0
    Izz = @. s2 *Iyy0 + c2 *Izz0
    Iyz = @. (Izz0 - Iyy0)*s*c

    ds2 = @. 2 *s*c
    dc2 = @. -2 *s*c
    dcs = @. c2 - s2

#     dA_dc .= @. (2.0/chord)*A
#     dIyy_dc .= (4.0/chord)*Iyy
#     dIzz_dc .= (4.0/chord)*Izz
#     dIyz_dc .= (4.0/chord)*Iyz
#
#     dIyy_dtheta .= @. dc2*Iyy0 + ds2*Izz0
#     dIzz_dtheta .= @. ds2*Iyy0 + dc2*Izz0
#     dIyz_dtheta .= @. (Iyy0 - Izz0)*dcs

    for i = 1:length(chord)

        dA_dc[i]   = (2.0/chord[i])*A[i]
        dIyy_dc[i] = (4.0/chord[i])*Iyy[i]
        dIzz_dc[i] = (4.0/chord[i])*Izz[i]
        dIyz_dc[i] = (4.0/chord[i])*Iyz[i]

        dIyy_dtheta[i] = dc2[i]*Iyy0[i] + ds2[i]*Izz0[i]
        dIzz_dtheta[i] = ds2[i]*Iyy0[i] + dc2[i]*Izz0[i]
        dIyz_dtheta[i] = (Izz0[i] - Iyy0[i])*dcs[i]
    end

end