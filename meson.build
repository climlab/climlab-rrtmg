project('climlab_rrtmg', 'c', 
  version : '0.4',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2'],
)

add_languages('fortran')
# Add Fortran compiler flags globally (you already have this)
add_project_arguments(
    '-fPIC', '-fno-range-check', '-w',
    language: 'fortran'
)

# Fortran source files: sw
rrtmg_sw_sources = files(
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/parkind.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/parrrsw.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_aer.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_cld.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_con.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg16.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg17.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg18.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg19.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg20.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg21.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg22.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg23.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg24.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg25.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg26.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg27.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg28.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_kg29.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_ncpar.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_ref.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_tbl.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_vsn.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/modules/rrsw_wvn.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_k_g.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/src/mcica_random_numbers.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/src/mcica_subcol_gen_sw.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_vrtqdr.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_reftra.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_taumol.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_spcvmc.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_setcoef.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_init.f90',
  'climlab_rrtmg/rrtmg_sw/rrtmg_sw_v4.0/gcm_model/src/rrtmg_sw_cldprmc.f90',
  'climlab_rrtmg/rrtmg_sw/sourcemods/rrtmg_sw_rad.f90',
  'climlab_rrtmg/rrtmg_sw/Driver.f90',
)
#unoptimized_src = ['rrtmg_sw_k_g.f90']

# Fortran source files: lw
rrtmg_lw_sources = files(
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/parkind.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/parrrtm.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_cld.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_con.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg01.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg02.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg03.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg04.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg05.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg06.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg07.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg08.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg09.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg10.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg11.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg12.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg13.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg14.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg15.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg16.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_ncpar.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_ref.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_tbl.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_vsn.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/modules/rrlw_wvn.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/src/rrtmg_lw_k_g.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/src/rrtmg_lw_taumol.f90',
  'climlab_rrtmg/rrtmg_lw/sourcemods/rrtmg_lw_setcoef.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/src/rrtmg_lw_rtrnmc.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/src/rrtmg_lw_cldprmc.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/src/mcica_random_numbers.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/src/mcica_subcol_gen_lw.f90',
  'climlab_rrtmg/rrtmg_lw/rrtmg_lw_v4.85/gcm_model/src/rrtmg_lw_init.f90',
  'climlab_rrtmg/rrtmg_lw/sourcemods/rrtmg_lw_rad.f90',
  'climlab_rrtmg/rrtmg_lw/Driver.f90',
)
#unoptimized_src = ['rrtmg_lw_k_g.f90']

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

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
  install: true
)

rrtmg_lw_library = static_library(
  'rrtmg_lw',
  rrtmg_lw_sources,
  include_directories: inc_np
)

py.extension_module('_rrtmg_lw',
  sources: [
    '_rrtmg_lwmodule.c',
    incdir_f2py / 'fortranobject.c'
  ],
  include_directories: inc_np,
  link_with: rrtmg_lw_library,
  dependencies: py_dep,
  install: true
)