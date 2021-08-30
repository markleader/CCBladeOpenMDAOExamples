using ConcreteStructs
using OpenMDAO: AbstractExplicitComp

function NACA0012(z)
    # For relative chord dimension z in [0, 1], return the +/- relative
    # thickness coordinates of a NACA 0012 airfoil

    y = 0.594689181*[0.298222773 *(z.^0.5) - 0.127125232 .*z - 0.357907906 .*(z.^2) + 0.291984971 .*(z.^3) - 0.105174606 .*(z.^4)]

    return y[1], -1 .*y[1]
end

@concrete struct StressComp <: AbstractExplicitComp
    num_elem
    num_stress_eval_points
end

function StressComp(; num_elem, num_stress_eval_points)

    return StressComp(num_elem, num_stress_eval_points)
end

function OpenMDAO.setup(self::StressComp)

    # Declare the OpenMDAO inputs.
    input_data = Vector{VarData}()
    push!(input_data, VarData("Nx", shape=self.num_elem, units="N"))
    push!(input_data, VarData("My", shape=self.num_elem, units="N*m"))
    push!(input_data, VarData("Mz", shape=self.num_elem, units="N*m"))
    push!(input_data, VarData("chord", shape=self.num_elem, units="m"))
    push!(input_data, VarData("twist", shape=self.num_elem, units="rad"))
    push!(input_data, VarData("A", shape=self.num_elem, units="m**2"))
    push!(input_data, VarData("Iyy", shape=self.num_elem, units="m**4"))
    push!(input_data, VarData("Izz", shape=self.num_elem, units="m**4"))
    push!(input_data, VarData("Iyz", shape=self.num_elem, units="m**4"))

    # Declare the OpenMDAO outputs.
    output_data = Vector{VarData}()
    push!(output_data, VarData("sigma1", shape=2*self.num_elem*self.num_stress_eval_points, units="N/m**2"))

    # Declare the OpenMDAO partial derivatives.
    partials_data = Vector{PartialsData}()
    push!(partials_data, PartialsData("sigma1", "Nx"))
    push!(partials_data, PartialsData("sigma1", "My"))
    push!(partials_data, PartialsData("sigma1", "Mz"))
    push!(partials_data, PartialsData("sigma1", "chord"))
    push!(partials_data, PartialsData("sigma1", "twist"))
    push!(partials_data, PartialsData("sigma1", "A"))
    push!(partials_data, PartialsData("sigma1", "Iyy"))
    push!(partials_data, PartialsData("sigma1", "Izz"))
    push!(partials_data, PartialsData("sigma1", "Iyz"))

    return input_data, output_data, partials_data
end

function OpenMDAO.compute!(self::StressComp, inputs, outputs)

    # Unpack the inputs.
    Nx = inputs["Nx"]
    My = inputs["My"]
    Mz = inputs["Mz"]
    chord = inputs["chord"]
    twist = inputs["twist"]
    A = inputs["A"]
    Iyy = inputs["Iyy"]
    Izz = inputs["Izz"]
    Iyz = inputs["Iyz"]

    # Unpack the outputs.
    sigma1 = outputs["sigma1"]

    s = sin.(twist)
    c = cos.(twist)
    dI = @. Iyy*Izz - Iyz*Iyz
    z_rel = range(0.0, 1.0, length=self.num_stress_eval_points) .- 0.421
    y1_rel, y2_rel = NACA0012(z_rel .+ 0.421)

    for i = 1:length(chord)
        for j = 1:self.num_stress_eval_points
            z1 = chord[i]*(z_rel[j]*c[i] + y1_rel[j]*s[i])
            y1 = chord[i]*(-z_rel[j]*s[i] + y1_rel[j]*c[i])
            z2 = chord[i]*(z_rel[j]*c[i] + y2_rel[j]*s[i])
            y2 = chord[i]*(-z_rel[j]*s[i] + y2_rel[j]*c[i])

            idx = 2*(i-1)*self.num_stress_eval_points + 2*j - 1
            sigma1[idx]   = Nx[i]/A[i] - My[i]*(y1*Iyz[i] - z1*Izz[i])/dI[i] - Mz[i]*(y1*Iyy[i] - z1*Iyz[i])/dI[i]
            sigma1[idx+1] = Nx[i]/A[i] - My[i]*(y2*Iyz[i] - z2*Izz[i])/dI[i] - Mz[i]*(y2*Iyy[i] - z2*Iyz[i])/dI[i]
        end
    end

end

function OpenMDAO.compute_partials!(self::StressComp, inputs, partials)

    # Unpack the inputs.
    Nx = inputs["Nx"]
    My = inputs["My"]
    Mz = inputs["Mz"]
    chord = inputs["chord"]
    twist = inputs["twist"]
    A = inputs["A"]
    Iyy = inputs["Iyy"]
    Izz = inputs["Izz"]
    Iyz = inputs["Iyz"]

    ds_dNx = partials["sigma1", "Nx"]
    ds_dMy = partials["sigma1", "My"]
    ds_dMz = partials["sigma1", "Mz"]
    ds_dc = partials["sigma1", "chord"]
    ds_dtheta = partials["sigma1", "twist"]
    ds_dA = partials["sigma1", "A"]
    ds_dIyy = partials["sigma1", "Iyy"]
    ds_dIzz = partials["sigma1", "Izz"]
    ds_dIyz = partials["sigma1", "Iyz"]

    s = sin.(twist)
    c = cos.(twist)
    dI = @. Iyy*Izz - Iyz*Iyz
    z_rel = range(0.0, 1.0, length=self.num_stress_eval_points) .- 0.421
    y1_rel, y2_rel = NACA0012(z_rel .+ 0.421)

    for i = 1:length(chord)
        for j = 1:self.num_stress_eval_points

            z1 = chord[i]*( z_rel[j]*c[i] + y1_rel[j]*s[i])
            y1 = chord[i]*(-z_rel[j]*s[i] + y1_rel[j]*c[i])
            z2 = chord[i]*( z_rel[j]*c[i] + y2_rel[j]*s[i])
            y2 = chord[i]*(-z_rel[j]*s[i] + y2_rel[j]*c[i])

            idx = 2*(i-1)*self.num_stress_eval_points + 2*j - 1

            ds_dNx[idx, i]   = 1/A[i]
            ds_dNx[idx+1, i] = 1/A[i]
            ds_dA[idx, i]    = -Nx[i]/A[i]^2
            ds_dA[idx+1, i]  = -Nx[i]/A[i]^2

            ds_dMy[idx, i]   = -(y1*Iyz[i] - z1*Izz[i])/dI[i]
            ds_dMy[idx+1, i] = -(y2*Iyz[i] - z2*Izz[i])/dI[i]
            ds_dMz[idx, i]   = -(y1*Iyy[i] - z1*Iyz[i])/dI[i]
            ds_dMz[idx+1, i] = -(y2*Iyy[i] - z2*Iyz[i])/dI[i]

            ds_dy = -(My[i]*Iyz[i] + Mz[i]*Iyy[i])/dI[i]
            ds_dz =  (My[i]*Izz[i] + Mz[i]*Iyz[i])/dI[i]

            ds_dc[idx, i]   = ds_dy*(-z_rel[j]*s[i] + y1_rel[j]*c[i]) + ds_dz*(z_rel[j]*c[i] + y1_rel[j]*s[i])
            ds_dc[idx+1, i] = ds_dy*(-z_rel[j]*s[i] + y2_rel[j]*c[i]) + ds_dz*(z_rel[j]*c[i] + y2_rel[j]*s[i])
            ds_dtheta[idx, i]   = chord[i]*(ds_dy*(-z_rel[j]*c[i] - y1_rel[j]*s[i]) + ds_dz*(-z_rel[j]*s[i] + y1_rel[j]*c[i]))
            ds_dtheta[idx+1, i] = chord[i]*(ds_dy*(-z_rel[j]*c[i] - y2_rel[j]*s[i]) + ds_dz*(-z_rel[j]*s[i] + y2_rel[j]*c[i]))

            ds_dIyy[idx, i]   = -(My[i]*Izz[i] + Mz[i]*Iyz[i])*(Izz[i]*z1 - Iyz[i]*y1)/dI[i]^2
            ds_dIyy[idx+1, i] = -(My[i]*Izz[i] + Mz[i]*Iyz[i])*(Izz[i]*z2 - Iyz[i]*y2)/dI[i]^2
            ds_dIzz[idx, i]   = (My[i]*Iyz[i] + Mz[i]*Iyy[i])*(Iyy[i]*y1 - Iyz[i]*z1)/dI[i]^2
            ds_dIzz[idx+1, i] = (My[i]*Iyz[i] + Mz[i]*Iyy[i])*(Iyy[i]*y2 - Iyz[i]*z2)/dI[i]^2
            ds_dIyz[idx, i]   = -((y1*(Iyy[i]*Izz[i] + Iyz[i]^2) - 2*z1*Izz[i]*Iyz[i])*My[i] + (2*y1*Iyy[i]*Iyz[i] - z1*(Iyy[i]*Izz[i] + Iyz[i]^2))*Mz[i])/(dI[i]^2)
            ds_dIyz[idx+1, i] = -((y2*(Iyy[i]*Izz[i] + Iyz[i]^2) - 2*z2*Izz[i]*Iyz[i])*My[i] + (2*y2*Iyy[i]*Iyz[i] - z2*(Iyy[i]*Izz[i] + Iyz[i]^2))*Mz[i])/(dI[i]^2)

        end
    end

end