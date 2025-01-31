python_sources = [
  '__init__.py',
]

py.install_sources(
  python_sources,
  pure: false,
  subdir: 'climlab_rrtmg/rrtmg_lw'
)

# Fortran source files: lw
rrtmg_lw_sources = files(
  'rrtmg_lw_v4.85/gcm_model/modules/parkind.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/parrrtm.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_cld.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_con.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg01.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg02.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg03.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg04.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg05.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg06.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg07.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg08.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg09.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg10.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg11.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg12.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg13.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg14.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg15.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_kg16.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_ncpar.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_ref.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_tbl.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_vsn.f90',
  'rrtmg_lw_v4.85/gcm_model/modules/rrlw_wvn.f90',
  'rrtmg_lw_v4.85/gcm_model/src/rrtmg_lw_k_g.f90',
  'rrtmg_lw_v4.85/gcm_model/src/rrtmg_lw_taumol.f90',
  'sourcemods/rrtmg_lw_setcoef.f90',
  'rrtmg_lw_v4.85/gcm_model/src/rrtmg_lw_rtrnmc.f90',
  'rrtmg_lw_v4.85/gcm_model/src/rrtmg_lw_cldprmc.f90',
  'rrtmg_lw_v4.85/gcm_model/src/mcica_random_numbers.f90',
  'rrtmg_lw_v4.85/gcm_model/src/mcica_subcol_gen_lw.f90',
  'rrtmg_lw_v4.85/gcm_model/src/rrtmg_lw_init.f90',
  'sourcemods/rrtmg_lw_rad.f90',
  'Driver.f90',
)
#unoptimized_src = ['rrtmg_lw_k_g.f90']



incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

lwmodule = custom_target('_rrtmg_lw.c',
   input : 'Driver.f90',  
   output : ['_rrtmg_lwmodule.c', 
            #  '_rrtmg_lw-f2pywrappers2.f90',
            #  '_rrtmg_lw-f2pywrappers.f',
             ],
   command : [py, '-m', 'numpy.f2py', '@INPUT@', 
              '-m', '_rrtmg_lw', '--lower',
              '--build-dir', '@OUTDIR@']
)

rrtmg_lw_library = static_library(
  'rrtmg_lw',
  rrtmg_lw_sources,
  include_directories: inc_np
)

py.extension_module('_rrtmg_lw',
  sources: [
    lwmodule,
    incdir_f2py / 'fortranobject.c'
  ],
  include_directories: inc_np,
  link_with: rrtmg_lw_library,
  dependencies: py_dep,
  install: true,
  subdir: 'climlab_rrtmg/rrtmg_lw'
)