from peritheos.eos import ideal_gas, van_der_waals

def main():
    # Example usage of the equation of state functions
    pressure = ideal_gas.calculate_pressure(temperature=298.15, volume=0.024, moles=1.0)
    print(f"Ideal gas pressure: {pressure:.2f} Pa")
    
    # Van der Waals example
    a = 0.1364  # (Pa·m^6/mol^2) for nitrogen
    b = 3.87e-5  # (m^3/mol) for nitrogen
    pressure_vdw = van_der_waals.calculate_pressure(temperature=298.15, volume=0.024, moles=1.0, a=a, b=b)
    print(f"Van der Waals pressure: {pressure_vdw:.2f} Pa")


if __name__ == "__main__":
    main()
