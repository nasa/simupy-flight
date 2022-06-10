import os
import pathlib
from simupy_flight.parse_daveml import ProcessDaveML

here = os.path.abspath(os.path.dirname(__file__))

base_dir = os.path.join(here, "..", "NESC_data", "All_models")
out_dir = here

filenames = [
    os.path.join("F16_package", "F16_S119_source", "F16_control.dml"),
    os.path.join("F16_package", "F16_S119_source", "F16_aero.dml"),
    os.path.join("F16_package", "F16_S119_source", "F16_prop.dml"),
    os.path.join("F16_package", "F16_S119_source", "F16_gnc.dml"),
    os.path.join("F16_package", "F16_S119_source", "F16_inertia.dml"),
]

for filename in filenames:
    ProcessDaveML(os.path.join(base_dir, filename), outdir=out_dir)
