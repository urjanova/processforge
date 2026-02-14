import pandas as pd

def generate_validation_report(input_csv, output_excel):
    # 1. Load your simulation data
    df = pd.read_csv(input_csv)
    
    # 2. Conservation of Mass Calculation
    # We check if (Inflow - Outflow) is within a 0.5% tolerance
    df['Mass_Delta'] = df['Inlet_Mass_Flow'] - df['Outlet_Mass_Flow']
    df['Mass_Error_Percent'] = (df['Mass_Delta'] / df['Inlet_Mass_Flow']) * 100
    
    # 3. Energy Consistency Check
    # Simple logic: If it's an exothermic process, Outlet Temp MUST be > Inlet Temp
    df['Energy_Logic_Check'] = df['Outlet_Temp_K'] > df['Inlet_Temp_K']
    
    # 4. Human-Readable Status
    df['Validation_Status'] = df.apply(
        lambda x: '✅ PASS' if abs(x['Mass_Error_Percent']) < 0.5 and x['Energy_Logic_Check'] 
        else '❌ FAIL', axis=1
    )
    
    # 5. Export to Excel for the non-coder
    df.to_excel(output_excel, index=False)
    print(f"Validation report ready: {output_excel}")


if __name__ == "__main__":
    # Example usage
    generate_validation_report('outputs/closed-loop-chain_timeseries.csv', 'validation_report.xlsx')