# Fortran source files: sw
rrtmg_sw_sources = files(
  'rrtmg_sw_v4.0/gcm_model/modules/parkind.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/parrrsw.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_aer.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_cld.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_con.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg16.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg17.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg18.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg19.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg20.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg21.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg22.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg23.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg24.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg25.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg26.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg27.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg28.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg29.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_ncpar.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_ref.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_tbl.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_vsn.f90',
  'rrtmg_sw_v4.0/gcm_model/modules/rrsw_wvn.f90',
  'rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_k_g.f90',
  'rrtmg_sw_v4.0/gcm_model/src/mcica_random_numbers.f90',
  'rrtmg_sw_v4.0/gcm_model/src/mcica_subcol_gen_sw.f90',
  'rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_vrtqdr.f90',
  'rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_reftra.f90',
  'rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_taumol.f90',
  'rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_spcvmc.f90',
  'rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_setcoef.f90',
  'rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_init.f90',
  'rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_cldprmc.f90',
  'sourcemods/rrtmg_sw_rad.f90',
  'Driver.f90',
)
#unoptimized_src = ['rrtmg_sw_k_g.f90']

python_sources = [
   '__init__.py',
]

py.install_sources(
   python_sources,
   pure: false,
   subdir: 'climlab_rrtmg/rrtmg_sw'
)
 
 incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

rrtmg_sw_library = static_library(
  'rrtmg_sw',
  rrtmg_sw_sources,
  include_directories: inc_np
)

py.extension_module('_rrtmg_sw',
  sources: [
    '_rrtmg_swmodule.c',
    incdir_f2py / 'fortranobject.c'
  ],
  include_directories: inc_np,
  link_with: rrtmg_sw_library,
  dependencies: py_dep,
  install: true,
  subdir: 'climlab_rrtmg/rrtmg_sw'
)
