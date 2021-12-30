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
    push!(input_data, VarData("twist", shape=self.nelems, units="rad"))

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
    push!(partials_data, PartialsData("Iyy", "twist", rows=rows, cols=cols))
    push!(partials_data, PartialsData("Izz", "chord", rows=rows, cols=cols))
    push!(partials_data, PartialsData("Izz", "twist", rows=rows, cols=cols))
    push!(partials_data, PartialsData("Iyz", "chord", rows=rows, cols=cols))
    push!(partials_data, PartialsData("Iyz", "twist", rows=rows, cols=cols))

    return input_data, output_data, partials_data
end

function OpenMDAO.compute!(self::AreaComp, inputs, outputs)

    # Unpack the inputs.
    chord = inputs["chord"]
    twist = inputs["twist"]

    # Unpack the outputs.
    A = outputs["A"]
    Iyy = outputs["Iyy"]
    Izz = outputs["Izz"]
    Iyz = outputs["Iyz"]

    s = sin.(twist)
    c = cos.(twist)
    s2 = sin.(twist).^2
    c2 = cos.(twist).^2

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
    twist = inputs["twist"]

    dA_dc = partials["A", "chord"]
    dIyy_dc = partials["Iyy", "chord"]
    dIyy_dtwist = partials["Iyy", "twist"]
    dIzz_dc = partials["Izz", "chord"]
    dIzz_dtwist = partials["Izz", "twist"]
    dIyz_dc = partials["Iyz", "chord"]
    dIyz_dtwist = partials["Iyz", "twist"]

    s = sin.(twist)
    c = cos.(twist)
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
#     dIyy_dtwist .= @. dc2*Iyy0 + ds2*Izz0
#     dIzz_dtwist .= @. ds2*Iyy0 + dc2*Izz0
#     dIyz_dtwist .= @. (Iyy0 - Izz0)*dcs

    for i = 1:length(chord)

        dA_dc[i]   = (2.0/chord[i])*A[i]
        dIyy_dc[i] = (4.0/chord[i])*Iyy[i]
        dIzz_dc[i] = (4.0/chord[i])*Izz[i]
        dIyz_dc[i] = (4.0/chord[i])*Iyz[i]

        dIyy_dtwist[i] = dc2[i]*Iyy0[i] + ds2[i]*Izz0[i]
        dIzz_dtwist[i] = ds2[i]*Iyy0[i] + dc2[i]*Izz0[i]
        dIyz_dtwist[i] = (Izz0[i] - Iyy0[i])*dcs[i]
    end

end