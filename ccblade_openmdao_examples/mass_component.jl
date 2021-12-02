using ConcreteStructs
using OpenMDAO: AbstractExplicitComp, VarData, PartialsData

@concrete struct MassComp <: AbstractExplicitComp
    rho
    span
    nelems
end

function MassComp(; rho, span, nelems)

    return MassComp(rho, span, nelems)
end

function OpenMDAO.setup(self::MassComp)

    # Declare the OpenMDAO inputs.
    input_data = Vector{VarData}()
    push!(input_data, VarData("A", shape=self.nelems, units="m**2"))

    # Declare the OpenMDAO outputs.
    output_data = Vector{VarData}()
    push!(output_data, VarData("m", val=0.0, units="kg"))

    # Declare the OpenMDAO partial derivatives.
    partials_data = Vector{PartialsData}()
    push!(partials_data, PartialsData("m", "A"))

    return input_data, output_data, partials_data
end

function OpenMDAO.compute!(self::MassComp, inputs, outputs)

    # Unpack the inputs.
    A = inputs["A"]

    # Unpack the outputs.
    m = outputs["m"]

    m .= (sum(A)/length(A))*self.span*self.rho

end

function OpenMDAO.compute_partials!(self::MassComp, inputs, partials)

    # Unpack the inputs.
    A = inputs["A"]
    dm_dA = partials["m", "A"]
    dm_dA[1, :] .= ones(length(A))*self.span*self.rho/length(A)

end