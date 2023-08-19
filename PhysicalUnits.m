classdef PhysicalUnits
    %stores and converts the units of SI to FD
    %Use with FDTD method
    properties
        % FDTD units
        speedOfLight_FD = 1.0;
        electronCharge_FD = 1.0;
        length_FD = 1.0;
        vacuumPermittivity_FD = 1.0;
        vacuumPermeability_FD = 1.0;

        % SI Units
        speedOfLight_SI = 3e8;
        electronCharge_SI = 1.602e-19;
        electronMass_SI = 3.109e-31;
        vacuumPermittivity_SI = 8.854e-12;
        vacuumPermeability_SI = 1.257e-6;
        length_SI = 10e-6;     % FDTD unit length in SI units, it can be set based on a specific wavelength for example

        % conversion factors
        % _SI_to_FD : SI unit/ FD unit ratio
        length_SI_to_FD;     % L_SI/L_FD
        area_SI_to_FD;
        volume_SI_to_FD;
        time_SI_to_FD;
        frequency_SI_to_FD;
        electricPotential_SI_to_FD;
        electricField_SI_to_FD;
        electricCurrent_SI_to_FD;
        energy_SI_to_FD;
        c_SI_to_FD;
        mass_SI_to_FD;
        magneticField_SI_to_FD;
    end
    
    methods
        function obj = PhysicalUnits(length_si)
            obj.length_SI = length_si;
           % set conversion factors
            obj.length_SI_to_FD = obj.length_SI / obj.length_FD;
            obj.area_SI_to_FD = obj.length_SI_to_FD * obj.length_SI_to_FD;
            obj.volume_SI_to_FD = obj.area_SI_to_FD * obj.length_SI_to_FD;
            obj.time_SI_to_FD = (obj.length_SI / obj.speedOfLight_SI) / (obj.length_FD / obj.speedOfLight_FD);
            obj.frequency_SI_to_FD = 1.0 / obj.time_SI_to_FD;
            obj.electricPotential_SI_to_FD = (obj.electronCharge_SI / (obj.vacuumPermittivity_SI * obj.length_SI)) / (obj.electronCharge_FD / (obj.vacuumPermittivity_FD * obj.length_FD));
            obj.electricField_SI_to_FD = (obj.electronCharge_SI / (obj.vacuumPermittivity_SI * obj.length_SI * obj.length_SI)) / (obj.electronCharge_FD / (obj.vacuumPermittivity_FD * obj.length_FD * obj.length_FD));
            obj.electricCurrent_SI_to_FD = (obj.electronCharge_SI / obj.electronCharge_FD) / obj.time_SI_to_FD;
            obj.energy_SI_to_FD = (obj.electronCharge_SI / obj.electronCharge_FD) * obj.electricPotential_SI_to_FD;
            obj.c_SI_to_FD = obj.speedOfLight_SI / obj.speedOfLight_FD;
            obj.mass_SI_to_FD = obj.energy_SI_to_FD /  (obj.c_SI_to_FD * obj.c_SI_to_FD);
            obj.magneticField_SI_to_FD = obj.electricField_SI_to_FD * obj.c_SI_to_FD * (obj.vacuumPermittivity_SI/obj.vacuumPermittivity_FD);
        end
        
        function out = ConvertFDLengthToSIUnits(obj,l_fd)
            out = l_fd * obj.length_SI_to_FD;
        end
        
        function out = ConvertSILengthToFDUnits(obj,l_si)
            out = l_si / obj.length_SI_to_FD;
        end

        function out = ConvertFDAreaToSIUnits(obj,a_fd)
            out = a_fd * obj.area_SI_to_FD;
        end

        function out = ConvertSIAreaToFDUnits(obj,a_si)
            out = a_si / obj.area_SI_to_FD;
        end
        
        function out = ConvertFDVolumeToSIUnits(obj,v_fd)
            out = v_fd * obj.volume_SI_to_FD;
        end

        function out = ConvertSIVolumeToFDUnits(obj,v_si)
            out = v_si / obj.volume_SI_to_FD;
        end

        function out = ConvertFDTimeToSIUnits(obj,t_fd)
            out = t_fd * obj.time_SI_to_FD;
        end

        function out = ConvertSITimeToFDUnits(obj,t_si)
            out = t_si / obj.time_SI_to_FD;
        end

        function out = ConvertFDFrequencyToSIUnits(obj,f_fd)
            out = f_fd * obj.frequency_SI_to_FD;
        end

        function out = ConvertSIFrequencyToFDUnits(obj,f_si)
            out = f_si / obj.frequency_SI_to_FD;
        end

        function out = ConvertFDElectricPotentialToSIUnits(obj,v_fd)
            out = v_fd * obj.electricPotential_SI_to_FD;
        end

        function out = ConvertSIElectricPotentialToFDUnits(obj,v_si)
            out = v_si / obj.electricPotential_SI_to_FD;
        end

        function out = ConvertFDElectricFieldToSIUnits(obj,eField_fd)
            out = eField_fd * obj.electricField_SI_to_FD;
        end

        function out = ConvertSIElectricFieldToFDUnits(obj,eField_si)
            out = eField_si / obj.electricField_SI_to_FD;
        end

        function out = ConvertSIEnergyToFDUnits(obj,energy_si)
            out = energy_si / obj.energy_SI_to_FD;
        end

        function out = ConvertFDEnergyToSIUnits(obj,energy_fd)
            out = energy_fd * obj.energy_SI_to_FD;
        end
        
        function out = ConvertSIMassToFDUnits(obj,mass_si)
            out = mass_si / obj.mass_SI_to_FD;
        end

        function out = ConvertFDMassToSIUnits(obj,mass_fd)
            out = mass_fd * obj.mass_SI_to_FD;
        end
        
        function out = ConvertSIMagneticFieldToFDUnits(obj,hField_si)
            out = hField_si / obj.magneticField_SI_to_FD;
        end
        
        function out = ConvertFDMagneticFieldToSIUnits(obj,hField_fd)
            out = hField_fd * obj.magneticField_SI_to_FD;
        end
    end
end