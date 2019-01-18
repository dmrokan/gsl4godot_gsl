import os

excluded_dir_names = [
    'autom4te.cache'
    , 'build'
    , 'doc'
    , 'doc_texinfo'
    , 'test'
    , '.git'
]

def source_dirs(root_dir, dir_list, cwd='./'):
    files = os.listdir(root_dir)

    for f in files:
        file_name = cwd + f

        if os.path.isdir(file_name):
            if f in excluded_dir_names:
                continue

            dir_list.append(file_name)
            source_dirs(file_name, dir_list, file_name + '/')

def scons_paths(fh, dirs):
    fh.write('scons_paths = [\n')
    fh.write('    \'.\'\n')
    fh.write('    , \'build/\'\n')

    for d in dirs:
        fh.write('    , \'' + d + '/\'\n')

    fh.write(']\n\n')

def scons_sources(fh, dirs):
    fh.write('scons_sources = [\n')


    fh.write('    [x for x in Glob(\'' + '.' + '/*.c\') ')
    fh.write(' if not os.path.basename(str(x)).endswith(\'_source.c\') ')
    fh.write('and os.path.basename(str(x)) not in excluded_files ')
    fh.write('and not os.path.basename(str(x)).startswith(\'demo\') ')
    fh.write('and not os.path.basename(str(x)).startswith(\'c_pass_\') ')
    fh.write('and not os.path.basename(str(x)).startswith(\'hc_pass_\') ')
    fh.write('and not os.path.basename(str(x)).startswith(\'fp-\') ')
    fh.write('and not os.path.basename(str(x)).startswith(\'real_pass_\') ')
    fh.write('and not os.path.basename(str(x)).startswith(\'test_\')]\n')

    for d in dirs:
        fh.write('    , [x for x in Glob(\'' + d + '/*.c\') ')
        fh.write(' if not os.path.basename(str(x)).endswith(\'_source.c\') ')
        fh.write('and os.path.basename(str(x)) not in excluded_files ')
        fh.write('and not os.path.basename(str(x)).startswith(\'demo\') ')
        fh.write('and not os.path.basename(str(x)).startswith(\'c_pass_\') ')
        fh.write('and not os.path.basename(str(x)).startswith(\'hc_pass_\') ')
        fh.write('and not os.path.basename(str(x)).startswith(\'fp-\') ')
        fh.write('and not os.path.basename(str(x)).startswith(\'real_pass_\') ')
        fh.write('and not os.path.basename(str(x)).startswith(\'test_\')]\n')

    fh.write('    , \'godot_gsl.cpp\'\n')
    fh.write('    , \'register_types.cpp\'\n')
    fh.write(']\n\n')

def excluded_files(fh):
    fh.write('excluded_files = [\n')
    fh.write('    \'common.c\'\n')
    fh.write('    , \'apply_givens.c\'\n')
    fh.write('    , \'cholesky_common.c\'\n')
    fh.write('    , \'svdstep.c\'\n')
    fh.write('    , \'lmiterate.c\'\n')
    fh.write('    , \'lmmisc.c\'\n')
    fh.write('    , \'lmpar.c\'\n')
    fh.write('    , \'lmset.c\'\n')
    fh.write('    , \'linear_common.c\'\n')
    fh.write('    , \'qrsolv.c\'\n')
    fh.write('    , \'hypot.c\'\n')
    fh.write('    , \'tests.c\'\n')
    fh.write('    , \'urand.c\'\n')
    fh.write('    , \'c_init.c\'\n')
    fh.write('    , \'c_main.c\'\n')
    fh.write('    , \'bitreverse.c\'\n')
    fh.write('    , \'c_radix2.c\'\n')
    fh.write('    , \'bitreverse.c\'\n')
    fh.write('    , \'factorize.c\'\n')
    fh.write('    , \'hc_init.c\'\n')
    fh.write('    , \'hc_main.c\'\n')
    fh.write('    , \'hc_radix2.c\'\n')
    fh.write('    , \'hc_unpack.c\'\n')
    fh.write('    , \'real_init.c\'\n')
    fh.write('    , \'real_main.c\'\n')
    fh.write('    , \'real_radix2.c\'\n')
    fh.write('    , \'real_unpack.c\'\n')
    fh.write('    , \'linear_minimize.c\'\n')
    fh.write('    , \'linear_wrapper.c\'\n')
    fh.write('    , \'directional_minimize.c\'\n')
    fh.write('    , \'find.c\'\n')
    fh.write('    , \'find2d.c\'\n')
    fh.write('    , \'qrstep.c\'\n')
    fh.write('    , \'beta_inc.c\'\n')
    fh.write('    , \'testint.c\'\n')
    fh.write('    , \'ringbuf.c\'\n')
    fh.write('    , \'deque.c\'\n')
    fh.write('    , \'cquad_const.c\'\n')
    fh.write('    , \'initialise.c\'\n')
    fh.write('    , \'append.c\'\n')
    fh.write('    , \'set_initial.c\'\n')
    fh.write('    , \'qpsrt.c\'\n')
    fh.write('    , \'util.c\'\n')
    fh.write('    , \'qc25s.c\'\n')
    fh.write('    , \'qc25c.c\'\n')
    fh.write('    , \'qelg.c\'\n')
    fh.write('    , \'err.c\'\n')
    fh.write('    , \'reset.c\'\n')
    fh.write('    , \'qpsrt2.c\'\n')
    fh.write('    , \'positivity.c\'\n')
    fh.write('    , \'qc25f.c\'\n')
    fh.write('    , \'ptsort.c\'\n')
    fh.write('    , \'companion.c\'\n')
    fh.write('    , \'balance.c\'\n')
    fh.write('    , \'qr.c\'\n')
    fh.write('    , \'demof.c\'\n')
    fh.write('    , \'nielsen.c\'\n')
    fh.write('    , \'step_utils.c\'\n')
    fh.write('    , \'control_utils.c\'\n')
    fh.write('    , \'rksubs.c\'\n')
    fh.write('    , \'modnewton1.c\'\n')
    fh.write('    , \'cheb_eval.c\'\n')
    fh.write('    , \'cheb_eval_mode.c\'\n')
    fh.write('    , \'enorm.c\'\n')
    fh.write('    , \'dogleg.c\'\n')
    fh.write('    , \'test.c\'\n')
    fh.write('    , \'endian.c\'\n')
    fh.write('    , \'standardize.c\'\n')
    fh.write('    , \'lmutil.c\'\n')
    fh.write('    , \'hess.c\'\n')
    fh.write('    , \'wigner.c\'\n')
    fh.write('    , \'benchmark.c\'\n')
    fh.write('    , \'errs.c\'\n')
    fh.write('    , \'poly.c\'\n')
    fh.write('    , \'testgen.c\'\n')
    fh.write('    , \'norm.c\'\n')
    fh.write(']\n\n')

def SCsub(root_dir):
    with open(root_dir + '/' + 'SCsub', 'w') as sc:
        sc.write('#!/usr/bin/env python\n')
        sc.write('# SCsub\n')
        sc.write('import re\n')
        sc.write('import os\n')

        sc.write('Import(\'env\')\n\n')

        dir_list = []
        source_dirs(root_dir, dir_list)
        excluded_files(sc)
        scons_paths(sc, dir_list)
        scons_sources(sc, dir_list)

        sc.write('env_godot_gsl = env\n')
        sc.write('env_godot_gsl.Append(CPPPATH=scons_paths)\n')
        sc.write('env_godot_gsl.add_source_files(env.modules_sources, scons_sources)')

def main():
    root_dir = os.path.dirname(os.path.realpath(__file__))
    dir_list = []

    SCsub(root_dir)
    # source_dirs(root_dir, dir_list)

    # print(dir_list)


if __name__ == "__main__":
    main()
