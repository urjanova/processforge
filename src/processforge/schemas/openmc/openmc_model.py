from enum import Enum
from typing import Dict, List, Optional, Union

from pydantic import BaseModel, Field, field_validator




# Material definition
class NuclideFraction(BaseModel):
    nuclide: str = Field(examples=["U235", "H1", "O16"], description="Nuclide symbol")  # e.g., 'U235', 'H1', 'O16'
    fraction: float = Field(examples=[0.5, 0.3, 0.2], description="atom or weight fraction")


class ElementFraction(BaseModel):
    element: str = Field(examples=["Fe", "C", "O"], description="Element symbol")  # e.g., 'Fe', 'C', 'O'
    percent: float = Field(..., description="atom or weight fraction")
    percent_type: str = Field(
        examples=["ao", "wo"],
        description=" ‘ao’ for atom percent and ‘wo’ for weight percent.",
    )  # 'atom' or 'weight'


class Material(BaseModel):
    name: str = Field(
        examples=["divertor_upper", "blanket_rear_wall"],
        description="Descriptive name of the material.",
    )
    friendly_material_id: int = Field(examples=[42], description="Unique identifier for the material as an integer.")
    density: float = Field(examples=[10.5], description="Density of the material.")
    density_units: str = Field(examples=["g/cm3", "atoms/b-cm"], description="Units for the density.")  # 'g/cm3' or 'atoms/b-cm'
    nuclides: List[NuclideFraction] = Field(description="List of nuclides with their fractions in the material.")
    temperature: Optional[float] = Field(description="Temperature of the material in Kelvin.")
    elements: Optional[List[ElementFraction]] = None



class MaterialLibrary(BaseModel):
    materials: List[Material] = Field(description="List of materials in the library.")


# Geometry definitions
class Surface(BaseModel):
    surface_id: int = Field(
        examples=[1, 2, 3],
        description="Unique identifier for the surface as an integer.",
    )
    type: str = Field(examples=["plane", "sphere", "cylinder"], description="Type of the surface.")
    params: Dict[str, float] = Field(description="Parameters defining the surface.")
    boundary_type: Optional[str] = Field(examples=["vacuum", "reflective"], description="Boundary type of the surface.")


# Specialized plane definitions for OpenMC
class XPlane(BaseModel):
    surface_id: int = Field(
        examples=[1, 2, 3],
        description="Unique identifier for the surface as an integer.",
    )
    x0: float = Field(description="X-coordinate of the plane.")
    boundary_type: Optional[str] = Field(examples=["vacuum", "reflective"], description="Boundary type of the surface.")


class YPlane(BaseModel):
    surface_id: int = Field(
        examples=[1, 2, 3],
        description="Unique identifier for the surface as an integer.",
    )
    y0: float = Field(description="Y-coordinate of the plane.")
    boundary_type: Optional[str] = Field(examples=["vacuum", "reflective"], description="Boundary type of the surface.")


class ZPlane(BaseModel):
    surface_id: int = Field(
        examples=[1, 2, 3],
        description="Unique identifier for the surface as an integer.",
    )
    z0: float = Field(description="Z-coordinate of the plane.")
    boundary_type: Optional[str] = Field(examples=["vacuum", "reflective"], description="Boundary type of the surface.")


class Cell(BaseModel):
    cell_id: int = Field(examples=[1, 2, 3], description="Unique identifier for the cell as an integer.")
    name: Optional[str] = Field(
        default=None,
        examples=["Cell 1", "Cell 2"],
        description="Descriptive name of the cell.",
    )
    fill: Optional[Union[int, str]] = Field(
        default=None,
        examples=["Material 1", "Material 2"],
        description="Material ID or Universe ID.",
    )
    region: Optional[str] = Field(
        default=None,
        examples=["Region 1", "Region 2"],
        description="Boolean expression of surfaces.",
    )
    temperature: Optional[float] = Field(default=None, examples=[300.0, 600.0], description="Temperature of the cell.")
    material: Optional[int] = Field(default=None, examples=[1, 2], description="Material ID.")
    universe: Optional[int] = Field(default=None, examples=[1, 2], description="Universe ID.")
    metadata: Optional[Dict[str, Union[str, float, int]]] = Field(
        default=None,
        examples=[{"key": "value"}, {"key": 42}],
        description="Additional metadata for the cell.",
    )


class Universe(BaseModel):
    universe_id: int = Field(
        examples=[1, 2, 3],
        description="Unique identifier for the universe as an integer.",
    )
    cells: List[int]  # cell_ids
    name: Optional[str] = Field(
        default=None,
        examples=["Universe 1", "Universe 2"],
        description="Descriptive name of the universe.",
    )


class Lattice(BaseModel):
    lattice_id: int = Field(
        examples=[1, 2, 3],
        description="Unique identifier for the lattice as an integer.",
    )
    name: Optional[str] = Field(
        default=None,
        examples=["Lattice 1", "Lattice 2"],
        description="Descriptive name of the lattice.",
    )
    universes: List[List[int]] = Field(examples=[[1, 2], [3]], description="2D or 3D array of universe_ids.")
    pitch: List[float] = Field(description="[x, y, (z)]")
    lower_left: List[float] = Field(description="[x, y, (z)]")


# Tally definitions
class TallyFilter(BaseModel):
    type: str = Field(description="Type of the tally filter (e.g., 'energy', 'cell', 'material', etc.)")
    bins: List[Union[int, float, str]] |List[List[Union[int, float, str]]] = Field(description="List of bins for the tally filter, can be nested for multi-dimensional filters.")


class Tally(BaseModel):
    friendly_tally_id: int = Field(description="Unique identifier for the tally.")
    name: Optional[str] = Field(default=None, description="Descriptive name of the tally.")
    filters: List[TallyFilter] = Field(description="List of filters applied to the tally.")
    scores: List[str] = Field(description="List of scores to be tallied (e.g., ['flux', 'fission', 'absorption']).")
    nuclides: Optional[List[str]] = Field(default=None, description="List of nuclides to be tallied.")
    estimator: Optional[str] = Field(
        default=None,
        description="Estimation method used for the tally (e.g., 'tracklength', 'analog').",
    )


class BoundingBox(BaseModel):
    x_min: float = Field(examples=[0.0, 1.0], description="Minimum x-coordinate of the bounding box.")
    y_min: float = Field(examples=[0.0, 1.0], description="Minimum y-coordinate of the bounding box.")
    z_min: float = Field(examples=[0.0, 1.0], description="Minimum z-coordinate of the bounding box.")
    x_max: float = Field(examples=[1.0, 2.0], description="Maximum x-coordinate of the bounding box.")
    y_max: float = Field(examples=[1.0, 2.0], description="Maximum y-coordinate of the bounding box.")
    z_max: float = Field(examples=[1.0, 2.0], description="Maximum z-coordinate of the bounding box.")


# Top-level OpenMC model
class Geometry(BaseModel):
    root: int = Field(examples=[0], description="Universe ID of the root universe.")
    merge_surfaces: Optional[bool] = Field(default=False, description="Flag to indicate if surfaces should be merged.")
    surface_precision: Optional[float] = Field(default=1e-6, description="Precision for merging surfaces.")
    bounding_box: Optional[BoundingBox] = Field(default=None, description="Bounding box for the geometry.")


# OpenMC settings schema
class SourceSettings(BaseModel):
    space: Optional[Dict[str, float]] = Field(default=None)  # e.g., {'x': 0.0, 'y': 0.0, 'z': 0.0}
    energy: Optional[float] = Field(default=None)  # e.g., 14e6 for 14 MeV
    angle: Optional[Dict[str, float]] = Field(default=None)  # e.g., {'mu': 1.0}


class SourcePointSettings(BaseModel):
    batches: Optional[List[int]] = Field(
        default=None,
        examples=[[10, 20]],
        description="List of batches at which to write source.",
    )
    overwrite: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Whether to overwrite existing source files.",
    )
    separate: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Whether the source should be written as a separate file.",
    )
    write: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Whether or not to write the source.",
    )
    mcpl: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Whether to write the source as an MCPL file.",
    )


class SurfSourceWriteSettings(BaseModel):
    surface_ids: Optional[List[int]] = Field(
        default=None,
        description="List of surface ids at which crossing particles are to be banked.",
    )
    max_particles: Optional[int] = Field(
        default=None,
        description="Maximum number of particles to be banked on surfaces per process.",
    )
    max_source_files: Optional[int] = Field(
        default=None,
        description="Maximum number of surface source files to be created.",
    )
    mcpl: Optional[bool] = Field(default=None, description="Output in the form of an MCPL-file.")
    cell: Optional[int] = Field(
        default=None,
        description="Cell ID for banking particles crossing identified surfaces (from or to this cell).",
    )
    cellfrom: Optional[int] = Field(
        default=None,
        description="Cell ID for banking particles coming from this declared cell.",
    )
    cellto: Optional[int] = Field(
        default=None,
        description="Cell ID for banking particles going to this declared cell.",
    )


class RandomRaySettings(BaseModel):
    distance_inactive: Optional[float] = Field(default=None, description="Inactive distance for random ray solver.")
    distance_active: Optional[float] = Field(default=None, description="Active distance for random ray solver.")
    ray_source: Optional[Dict] = Field(default=None, description="Source settings for random ray solver.")
    volume_estimator: Optional[str] = Field(
        default=None,
        examples=["hybrid", "analog"],
        description="Type of volume estimator to use.",
    )
    source_shape: Optional[str] = Field(
        default=None,
        examples=["flat", "sphere"],
        description="Shape of the source region.",
    )
    volume_normalized_flux_tallies: Optional[bool] = Field(default=None, description="Whether to use volume-normalized flux tallies.")
    adjoint: Optional[bool] = Field(default=None, description="Whether to use adjoint random ray solver.")
    sample_method: Optional[str] = Field(
        default=None,
        examples=["prng", "quasirandom"],
        description="Sampling method for random rays.",
    )
    source_region_meshes: Optional[List] = Field(default=None, description="List of meshes for source region.")


class ResonanceScatteringSettings(BaseModel):
    enable: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Enable resonance elastic scattering.",
    )
    method: Optional[str] = Field(
        default=None,
        examples=["rvs", "dbrc"],
        description="Method for resonance elastic scattering (e.g., 'rvs', 'dbrc').",
    )
    energy_min: Optional[float] = Field(
        default=None,
        examples=[0.01],
        description="Minimum energy for resonance scattering (in eV).",
    )
    energy_max: Optional[float] = Field(
        default=None,
        examples=[1000.0],
        description="Maximum energy for resonance scattering (in eV).",
    )
    nuclides: Optional[List[str]] = Field(
        default=None,
        examples=[["U238", "U235"]],
        description="List of nuclides for which resonance scattering is enabled.",
    )

    class Config:
        json_schema_extra = {
            "examples": [
                {
                    "enable": True,
                    "method": "rvs",
                    "energy_min": 0.01,
                    "energy_max": 1000.0,
                    "nuclides": ["U238", "U235"],
                }
            ]
        }


class StatepointSettings(BaseModel):
    batches: Optional[List[int]] = Field(
        default=None,
        examples=[[10, 20, 30]],
        description="List of batches at which to write state points.",
    )


class RunMode(str, Enum):
    fixed_source = "fixed source"
    eigenvalue = "eigenvalue"
    plot = "plot"
    volume = "volume"
    particle_restart = "particle restart"


class OutputSettings(BaseModel):
    path: Optional[str] = Field(default=None, examples=["output/"], description="Path to the output directory.")
    summary: Optional[bool] = Field(default=None, examples=[True, False], description="Whether to output summary files.")
    tallies: Optional[bool] = Field(default=None, examples=[True, False], description="Whether to output tallies files.")


class VolumeCalculationSettings(BaseModel):
    domains: Optional[List[int]] = Field(
        default=None,
        description="IDs of domains to find volumes of (cells, materials, or universes).",
    )
    domain_type: Optional[str] = Field(
        default=None,
        description="Type of each domain: 'cell', 'material', or 'universe'.",
    )
    samples: Optional[int] = Field(default=None, description="Number of samples used to generate volume estimates.")
    lower_left: Optional[List[float]] = Field(
        default=None,
        description="Lower-left coordinates of bounding box used to sample points.",
    )
    upper_right: Optional[List[float]] = Field(
        default=None,
        description="Upper-right coordinates of bounding box used to sample points.",
    )
    volumes: Optional[Dict[int, float]] = Field(
        default=None,
        description="Dictionary mapping unique IDs of domains to estimated volumes in cm^3.",
    )
    atoms: Optional[Dict[int, Dict[str, float]]] = Field(
        default=None,
        description="Dictionary mapping domain IDs to a mapping of nuclides to total number of atoms.",
    )
    threshold: Optional[float] = Field(
        default=None,
        description="Threshold for the maximum standard deviation of volumes.",
    )
    trigger_type: Optional[str] = Field(
        default=None,
        description="Value type used to halt volume calculation: 'variance', 'std_dev', or 'rel_err'.",
    )
    iterations: Optional[int] = Field(
        default=None,
        description="Number of iterations over samples (for calculations with a trigger).",
    )


class TemperatureSettings(BaseModel):
    default: Optional[float] = Field(
        default=None, examples=[293.6, 600.0], description="Default temperature in Kelvin for materials without specified temperatures."
    )
    method: Optional[str] = Field(
        default=None,
        examples=["nearest", "interpolation"],
        description="Method for treating intermediate temperatures between cross section data points.",
    )
    tolerance: Optional[float] = Field(default=None, examples=[10.0], description="Temperature tolerance in Kelvin for matching cross section data.")
    multipole: Optional[bool] = Field(
        default=None, examples=[True, False], description="Whether to use windowed multipole method for temperature treatment."
    )


class KeffTrigger(BaseModel):
    type: str = Field(examples=["std_dev", "variance", "rel_err"], description="Type of trigger: 'std_dev', 'variance', or 'rel_err'")
    threshold: float = Field(examples=[0.001, 0.0005], description="Threshold value for the trigger")


class OpenMCSetting(BaseModel):
    friendly_setting_id: int = Field(examples=[1], description="Unique identifier for the setting as an integer.")
    batches: int = Field(examples=[100, 200], description="Number of batches to simulate.")
    inactive: int = Field(examples=[10, 20], description="Number of inactive batches.")
    particles: int = Field(examples=[1000, 2000], description="Number of particles per generation.")
    source: Optional[SourceSettings] = Field(
        default=None,
        examples=[{"space": {"x": 0.0, "y": 0.0, "z": 0.0}, "energy": 14e6}],
        description="Distribution of source sites in space, angle, and energy.",
    )
    run_mode: Optional[RunMode] = Field(
        default=RunMode.eigenvalue,
        examples=[
            None,
            "fixed source",
            "eigenvalue",
            "plot",
            "volume",
            "particle restart",
        ],
        description="Type of calculation to perform.",
    )
    # In OpenMCSettings:
    output: Optional[OutputSettings] = Field(default=None, description="Settings for output files and directories.")

    temperature: Optional[TemperatureSettings] = Field(
        default=None,
        examples=[{"default": 293.6, "method": "nearest", "tolerance": 10.0}],
        description="Temperature settings for cross section data and material treatment.",
    )
    confidence_intervals: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="If True, uncertainties on tally results will be reported as the half-width of the 95% two-sided confidence interval.",
    )
    create_fission_neutrons: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Indicate whether fission neutrons should be created or not.",
    )
    cutoff: Optional[Dict[str, Union[float, bool]]] = Field(
        default=None,
        examples=[{"weight": 0.001, "energy": 1e-5, "time": 1e-9}],
        description="Dictionary defining weight, energy, and time cutoffs.",
    )
    delayed_photon_scaling: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Scale the fission photon yield by (EGP + EGD)/EGP.",
    )
    electron_treatment: Optional[str] = Field(
        default=None,
        examples=["led", "ttb"],
        description="How to treat electrons: 'led' or 'ttb'.",
    )
    energy_mode: Optional[str] = Field(
        default=None,
        examples=["continuous-energy", "multi-group"],
        description="Calculation mode: continuous-energy or multi-group.",
    )
    entropy_mesh: Optional[Dict[str, Union[int, float]]] = Field(
        default=None,
        examples=[
            {
                "lower_left": [0.0, 0.0],
                "upper_right": [10.0, 10.0],
                "dimension": [10, 10],
            }
        ],
        description="Mesh to be used to calculate Shannon entropy.",
    )
    event_based: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Use event-based parallelism instead of history-based.",
    )
    generations_per_batch: Optional[int] = Field(default=None, examples=[1, 2], description="Number of generations per batch.")
    max_lost_particles: Optional[int] = Field(
        default=None,
        examples=[100, 1000],
        description="Maximum number of lost particles.",
    )
    rel_max_lost_particles: Optional[float] = Field(
        default=None,
        examples=[0.1, 0.01],
        description="Maximum number of lost particles, relative to total.",
    )
    keff_trigger: Optional[KeffTrigger] = Field(
        default=None,
        examples=[{"type": "std_dev", "threshold": 0.001}],
        description="Trigger on eigenvalue.",
    )
    log_grid_bins: Optional[int] = Field(
        default=None,
        examples=[8000, 16000],
        description="Number of bins for logarithmic energy grid search.",
    )
    material_cell_offsets: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Generate an offset table for material cells.",
    )
    max_particles_in_flight: Optional[int] = Field(
        default=None,
        examples=[10000, 20000],
        description="Number of neutrons to run concurrently with event-based parallelism.",
    )
    max_particle_events: Optional[int] = Field(
        default=None,
        examples=[100000, 200000],
        description="Maximum number of allowed particle events per source particle.",
    )
    max_order: Optional[int] = Field(
        default=None,
        examples=[5, 10],
        description="Maximum scattering order to apply globally in multi-group mode.",
    )
    max_history_splits: Optional[int] = Field(
        default=None,
        examples=[100, 200],
        description="Maximum number of times a particle can split during a history.",
    )
    max_tracks: Optional[int] = Field(
        default=None,
        examples=[1000, 2000],
        description="Maximum number of tracks written to a track file per MPI process.",
    )
    max_write_lost_particles: Optional[int] = Field(
        default=None,
        examples=[10, 20],
        description="Maximum number of particle restart files to write for lost particles per MPI process.",
    )
    no_reduce: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Do not reduce tallies across processes in parallel calculation.",
    )
    photon_transport: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Whether to use photon transport.",
    )
    plot_seed: Optional[int] = Field(
        default=None,
        examples=[12345, 67890],
        description="Initial seed for randomly generated plot colors.",
    )
    ptables: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Determine whether probability tables are used.",
    )

    random_ray: Optional[RandomRaySettings] = Field(
        default=None,
        examples=[
            {
                "distance_inactive": 10.0,
                "distance_active": 20.0,
                "ray_source": {},
                "volume_estimator": "hybrid",
                "source_shape": "flat",
                "volume_normalized_flux_tallies": False,
                "adjoint": False,
                "sample_method": "prng",
                "source_region_meshes": [],
            }
        ],
        description="Options for configuring the random ray solver.",
    )
    resonance_scattering: Optional[ResonanceScatteringSettings] = Field(default=None, description="Settings for resonance elastic scattering.")
    seed: Optional[int] = Field(
        default=None,
        examples=[123456789],
        description="Seed for the pseudorandom number generator.",
    )
    stride: Optional[int] = Field(
        default=None,
        examples=[100, 200],
        description="Number of random numbers allocated for each source particle history.",
    )
    sourcepoint: Optional[SourcePointSettings] = Field(
        default=None,
        examples=[
            {
                "batches": [10, 20],
                "write": True,
                "overwrite": False,
                "separate": True,
                "mcpl": False,
            }
        ],
        description="Options for writing source points.",
    )
    statepoint: Optional[StatepointSettings] = Field(
        default=None,
        examples=[{"batches": [10, 20, 30]}],
        description="Options for writing state points.",
    )
    surf_source_read: Optional[Dict[str, str]] = Field(
        default=None,
        examples=[{"path": "surf_source.h5"}],
        description="Options for reading surface source points.",
    )
    surf_source_write: Optional[SurfSourceWriteSettings] = Field(default=None, description="Options for writing surface source points.")
    survival_biasing: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Indicate whether survival biasing is to be used.",
    )
    tabular_legendre: Optional[Dict[str, Union[bool, int]]] = Field(
        default=None,
        examples=[{"enable": True, "n_points": 8}],
        description="Multi-group scattering moment kernel tabular distribution settings.",
    )
    trace: Optional[Union[List[int], tuple]] = Field(
        default=None,
        examples=[[1, 2, 3], (1, 2, 3)],
        description="Show detailed information about a single particle (batch, generation, particle).",
    )
    track: Optional[Union[List[int], tuple]] = Field(
        default=None,
        examples=[[1, 2, 3], (1, 2, 3)],
        description="Specify particles for which track files should be written.",
    )
    trigger_active: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Indicate whether tally triggers are used.",
    )
    trigger_batch_interval: Optional[int] = Field(
        default=None,
        examples=[5, 10],
        description="Number of batches between convergence checks.",
    )
    trigger_max_batches: Optional[int] = Field(
        default=None,
        examples=[100, 200],
        description="Maximum number of batches simulated.",
    )
    uniform_source_sampling: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Sample among multiple sources uniformly.",
    )
    ufs_mesh: Optional[Dict[str, Union[int, float]]] = Field(
        default=None,
        examples=[
            {
                "dimension": [10, 10],
                "lower_left": [0.0, 0.0],
                "upper_right": [10.0, 10.0],
            }
        ],
        description="Mesh for uniform fission site method.",
    )
    use_decay_photons: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Produce decay photons from neutron reactions instead of prompt.",
    )
    verbosity: Optional[int] = Field(
        default=None,
        examples=[1, 5, 10],
        description="Verbosity during simulation (1-10).",
    )
    volume_calculations: Optional[Union[VolumeCalculationSettings, List[VolumeCalculationSettings]]] = Field(
        default=None,
        examples=[
            {
                "domains": [1, 2],
                "domain_type": "cell",
                "samples": 1000,
                "lower_left": [0.0, 0.0, 0.0],
                "upper_right": [10.0, 10.0, 10.0],
                "volumes": {1: 100.0, 2: 200.0},
                "threshold": 0.01,
                "trigger_type": "std_dev",
                "iterations": 5,
                "atoms": {1: {"U235": 1.0e22, "U238": 5.0e22}},
            }
        ],
        description="Stochastic volume calculation specifications and results.",
    )

    # atoms_dataframe is omitted as it is a runtime pandas.DataFrame, not a serializable schema field.
    weight_windows: Optional[Union[Dict, List[Dict]]] = Field(
        default=None,
        examples=[
            {
                "energy_bins": [0.0, 1.0e6, 20.0e6],
                "lower_weight": 0.1,
                "upper_weight": 10.0,
            }
        ],
        description="Weight windows for variance reduction.",
    )
    weight_window_checkpoints: Optional[Dict[str, bool]] = Field(
        default=None,
        examples=[{"enable": True}],
        description="Checkpoints for weight window split/roulettes.",
    )
    weight_window_generators: Optional[Union[Dict, List[Dict]]] = Field(
        default=None,
        examples=[{"method": "auto", "target": 0.5}],
        description="Weight window generation parameters.",
    )
    create_delayed_neutrons: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Whether delayed neutrons are created in fission.",
    )
    weight_windows_on: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Whether weight windows are enabled.",
    )
    weight_windows_file: Optional[str] = Field(
        default=None,
        examples=["weight_windows.h5"],
        description="Path to a weight window file to load during simulation initialization.",
    )
    write_initial_source: Optional[bool] = Field(
        default=None,
        examples=[True, False],
        description="Write the initial source distribution to file.",
    )



class OpenMCModel(BaseModel):
    materials: List[Material] = Field(description="List of materials in the model.")
    geometry: Geometry = Field(description="Geometry definition for the model.")
    tallies: Optional[List[Tally]] = Field(default=None, description="List of tallies in the model.")
