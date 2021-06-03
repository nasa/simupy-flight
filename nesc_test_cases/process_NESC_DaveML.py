import os
from parse_daveml import ProcessDaveML

filenames = [
    os.path.join('F16_package', 'F16_package', 'F16_S119_source', 'F16_control.dml'),
    os.path.join('F16_package', 'F16_package', 'F16_S119_source', 'F16_aero.dml'),
    os.path.join('F16_package', 'F16_package', 'F16_S119_source', 'F16_prop.dml'),

    os.path.join('F16_package', 'F16_package', 'F16_S119_source', 'F16_gnc.dml'),
    os.path.join('F16_package', 'F16_package', 'F16_S119_source', 'F16_inertia.dml'),

    os.path.join('two-stage_package', 'two-stage_package', 'twostage_aero.dml'),
    os.path.join('two-stage_package', 'two-stage_package', 'twostage_inertia.dml'),
    os.path.join('two-stage_package', 'two-stage_package', 'twostage_prop.dml'),



    'brick_aero.dml',
    'brick_inertia.dml',

    'cannonball_aero.dml',
    'cannonball_inertia.dml',

    'orbital_station_inertia.dml',
    'orbital_cylinder_inertia.dml',
    'orbital_sphere_inertia.dml',

]

for filename in filenames:
    ProcessDaveML(filename)