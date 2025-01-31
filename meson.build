project('climlab_rrtmg', 'c',)

add_languages('fortran')
# Add Fortran compiler flags globally (you already have this)
add_project_arguments(
    '-O0', '-fPIC', '-fno-range-check', '-w', '-fno-var-tracking-assignments',
    language: 'fortran'
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

subdir('climlab_rrtmg')
