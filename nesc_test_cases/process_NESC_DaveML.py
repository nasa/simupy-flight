import os
import pathlib
from simupy_flight.parse_daveml import ProcessDaveML

base_dir = os.path.join(pathlib.Path(__file__).parent.resolve(),
                        '..', 'NESC_data', 'All_models')

filenames = [
    os.path.join("F16_package", "F16_S119_source", "F16_control.dml"),
    os.path.join("F16_package", "F16_S119_source", "F16_aero.dml"),
    os.path.join("F16_package", "F16_S119_source", "F16_prop.dml"),
    os.path.join("F16_package", "F16_S119_source", "F16_gnc.dml"),
    os.path.join("F16_package", "F16_S119_source", "F16_inertia.dml"),
    os.path.join("two-stage_package", "twostage_aero.dml"),
    os.path.join("two-stage_package", "twostage_inertia.dml"),
    os.path.join("two-stage_package", "twostage_prop.dml"),
    "brick_aero.dml",
    "brick_inertia.dml",
    "cannonball_aero.dml",
    "cannonball_inertia.dml",
]

for filename in filenames:
    ProcessDaveML(os.path.join(base_dir, filename))
