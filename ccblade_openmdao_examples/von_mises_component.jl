using ConcreteStructs
using OpenMDAO: AbstractExplicitComp, VarData, PartialsData

@concrete struct VonMisesComp <: AbstractExplicitComp
    nelems
    num_stress_eval_points
end

function VonMisesComp(; nelems, num_stress_eval_points)

    return VonMisesComp(nelems, num_stress_eval_points)
end

function OpenMDAO.setup(self::VonMisesComp)

    # Declare the OpenMDAO inputs.
    input_data = Vector{VarData}()
    push!(input_data, VarData("sigma1", shape=2*self.nelems*self.num_stress_eval_points))

    # Declare the OpenMDAO outputs.
    output_data = Vector{VarData}()
    push!(output_data, VarData("sigma_vm", shape=2*self.nelems*self.num_stress_eval_points))

    # Declare the OpenMDAO partial derivatives.
    partials_data = Vector{PartialsData}()
    rows = collect(range(0, 2*self.nelems*self.num_stress_eval_points-1,
                         length=2*self.nelems*self.num_stress_eval_points))
    cols = rows
    push!(partials_data, PartialsData("sigma_vm", "sigma1", rows=rows, cols=cols))

    return input_data, output_data, partials_data
end

function OpenMDAO.compute!(self::VonMisesComp, inputs, outputs)

    # Unpack the inputs.
    sigma1 = inputs["sigma1"]

    # Unpack the outputs.
    sigma_vm = outputs["sigma_vm"]

    sigma_vm .= (sigma1.^2).^0.5

end

function OpenMDAO.compute_partials!(self::VonMisesComp, inputs, partials)

    # Unpack the inputs.
    sigma1 = inputs["sigma1"]
    dvm_ds = partials["sigma_vm", "sigma1"]
    sigma_vm = (sigma1.^2).^0.5

    for i = 1:length(sigma1)

        dvm_ds[i] = sigma1[i]/sigma_vm[i]

    end

end