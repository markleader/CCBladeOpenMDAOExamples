using ConcreteStructs
using OpenMDAO: AbstractExplicitComp, VarData, PartialsData

@concrete struct AreaComp <: AbstractExplicitComp
    num_elem
    A_ref
    Iyy_ref
    Izz_ref
end

function AreaComp(; num_elem, A_ref, Iyy_ref, Izz_ref)

    return AreaComp(num_elem, A_ref, Iyy_ref, Izz_ref)
end

function OpenMDAO.setup(self::AreaComp)

    # Declare the OpenMDAO inputs.
    input_data = Vector{VarData}()
    push!(input_data, VarData("chord", shape=self.num_elem, units="m"))
    push!(input_data, VarData("twist", shape=self.num_elem, units="rad"))

    # Declare the OpenMDAO outputs.
    output_data = Vector{VarData}()
    push!(output_data, VarData("A", shape=self.num_elem, units="m**2"))
    push!(output_data, VarData("Iyy", shape=self.num_elem, units="m**4"))
    push!(output_data, VarData("Izz", shape=self.num_elem, units="m**4"))
    push!(output_data, VarData("Iyz", shape=self.num_elem, units="m**4"))

    # Declare the OpenMDAO partial derivatives.
    partials_data = Vector{PartialsData}()
    push!(partials_data, PartialsData("A", "chord"))
    push!(partials_data, PartialsData("Iyy", "chord"))
    push!(partials_data, PartialsData("Iyy", "twist"))
    push!(partials_data, PartialsData("Izz", "chord"))
    push!(partials_data, PartialsData("Izz", "twist"))
    push!(partials_data, PartialsData("Iyz", "chord"))
    push!(partials_data, PartialsData("Iyz", "twist"))

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
    Iyz[1:end] = @. (Iyy0 - Izz0)*s*c

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
    Iyz = @. (Iyy0 - Izz0)*s*c

#     ds2 = @. 2 *s*c
#     dc2 = @. -2 *s*c
#     dcs = @. c2 - s2
#
#     dA_dc = @. (2 /chord)*A
#     dIyy_dc = @. (4 /chord)*Iyy
#     dIzz_dc = @. (4 /chord)*Izz
#     dIyz_dc = @. (4 /chord)*Iyz
#
#     dIyy_dtwist = @. dc2 *Iyy0 + ds2 *Izz0
#     dIzz_dtwist = @. ds2 *Iyy0 + dc2 *Izz0
#     dIyz_dtwist = @. (Iyy0 - Izz0)*dcs

    for i = 1:length(chord)

        ds2 = 2*s[i]*c[i]
        dc2 = -2*s[i]*c[i]
        dcs = c2[i] - s2[i]

        dA_dc[i, i]   = (2.0/chord[i])*A[i]
        dIyy_dc[i, i] = (4.0/chord[i])*Iyy[i]
        dIzz_dc[i, i] = (4.0/chord[i])*Izz[i]
        dIyz_dc[i, i] = (4.0/chord[i])*Iyz[i]

        dIyy_dtwist[i, i] = dc2*Iyy0[i] + ds2*Izz0[i]
        dIzz_dtwist[i, i] = ds2*Iyy0[i] + dc2*Izz0[i]
        dIyz_dtwist[i, i] = (Iyy0[i] - Izz0[i])*dcs
    end

end