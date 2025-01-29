project('climlab_rrtmg', 'c', 
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2'],
)

add_languages('fortran')
# Add Fortran compiler flags globally (you already have this)
add_project_arguments(
    '-O0', '-fPIC', '-ftree-vectorize', '-fno-range-check', '-fvisibility=hidden', 
    '-fstack-protector-strong', '-ffunction-sections', 
    '-pipe', '-Wno-unused-variable', language: 'fortran'
)

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

subdir('climlab_rrtmg')