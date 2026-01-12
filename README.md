# Process Forge

A Python-based process simulation framework for chemical process engineering applications.

## Features

- Steady-state and dynamic process simulations
- Thermodynamic property calculations
- Unit operations: Flash, Heater, etc.
- Flowsheet modeling and solving with closed-loop (recycle) support
- Automatic tear stream detection and convergence for recycle loops
- Results export to CSV and visualization

## Installation

1. Clone the repository:
   ```
   git clone <repository-url>
   cd processforge
   ```

2. Install the package using uv:
   ```
   uv sync
   ```

   This will create a virtual environment (`.venv`) and install all dependencies.

   **Note:** uv is the recommended package manager for this project. If you prefer pip, you can use `pip install -r requirements.txt` for basic dependencies, but uv provides better dependency resolution and virtual environment management.

3. Activate the virtual environment:
   ```
   source .venv/bin/activate  # On Unix/macOS
   # or
   .venv\Scripts\activate     # On Windows
   ```

## Usage

Run simulations using the `processforge` command:

```bash
processforge flowsheets/example_flash.json
```

Example flowsheet files are available in the `flowsheets/` directory.

## Project Structure

- `processforge/`: Core source code
  - `flowsheet.py`: Flowsheet modeling with closed-loop handling
  - `thermo.py`: Thermodynamic calculations
  - `units/`: Unit operations (flash, heater, solver)
  - `simulate.py`: Main simulation script
- `flowsheets/`: Example flowsheet configurations

## Dependencies

- numpy
- scipy
- coolprop
- networkx
- h5py
- casadi

## License

This project is proprietary software. See the [LICENSE](LICENSE) file for details.

For licensing inquiries, please [contact the development](https://forms.gle/wUweVnoSqA9VeD7m9) team.