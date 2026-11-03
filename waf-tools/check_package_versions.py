#! /usr/bin/env python
# encoding: utf-8

# Detect packages shadowed by another installation of the same package.
#
# The build merges include and library paths of all packages in the order of the 'use' list.
# If e.g. --elemental-dir also contains FLINT, its headers and library can shadow those from --flint-dir.
# No order works for all setups, so we detect the conflict instead.
# For each package, we print its header and library versions twice:
# with the package's own paths first ("expected"), and in the build order ("actual").

import glob, os, re

from waflib import Configure

# package -> (#include directives, header version expression, library version expression or None)
PROBES = {
    'gmpxx': ('#include <gmp.h>',
              'STR(__GNU_MP_VERSION) "." STR(__GNU_MP_VERSION_MINOR) "." STR(__GNU_MP_VERSION_PATCHLEVEL)',
              'gmp_version'),
    'mpfr': ('#include <mpfr.h>', 'MPFR_VERSION_STRING', 'mpfr_get_version()'),
    'boost': ('#include <boost/version.hpp>', 'BOOST_LIB_VERSION', None),
    'elemental': ('#include <El/config.h>', 'EL_VERSION_MAJOR "." EL_VERSION_MINOR " " EL_GIT_SHA1', None),
    'libxml2': ('#include <libxml/parser.h>', 'LIBXML_DOTTED_VERSION', 'xmlParserVersion'),
    'rapidjson': ('#include <rapidjson/rapidjson.h>', 'RAPIDJSON_VERSION_STRING', None),
    'libarchive': ('#include <archive.h>', 'ARCHIVE_VERSION_STRING', 'archive_version_string()'),
    'flint': ('#include <flint/flint.h>', 'FLINT_VERSION', 'flint_version'),
    # mps/version.h declares mps_get_version() without extern "C"
    'mpsolve': ('extern "C" {\n#include <mps/version.h>\n}',
                'STR(MPS_MAJOR_VERSION) "." STR(MPS_MINOR_VERSION) "." STR(MPS_PATCH_VERSION)',
                'mps_get_version()'),
}


# Returns {'header': version, 'library': version}, or None if the program fails to build or run
def _probe(conf, package, use):
    include, header, library = PROBES[package]
    lines = [include,
             '#include <cstdio>',
             '#define STR_(x) #x',
             '#define STR(x) STR_(x)',
             'int main()',
             '{',
             '  std::printf("header %%s\\n", %s);' % header]
    if library:
        lines.append('  std::printf("library %%s\\n", %s);' % library)
    lines.append('}')
    # Repeat all libraries at the end of the link command.
    # A static library (e.g. libflint.a) needs its dependencies (e.g. -lgmp) after it,
    # and our program alone does not keep them linked.
    ldflags = ['-l' + lib for p in use for lib in conf.env['LIB_' + p]]
    output = conf.check_cxx(fragment='\n'.join(lines), use=use, ldflags=ldflags, execute=True, define_ret=True,
                            mandatory=False)
    if not output:
        return None
    return dict(line.split(' ', 1) for line in output.splitlines())


# Files of the package found in directories of packages preceding it in the build order
def _shadowing_files(conf, package, use):
    headers = re.findall(r'<([^>]+)>', PROBES[package][0])
    libs = conf.env['LIB_' + package]
    own_dirs = conf.env['INCLUDES_' + package] + conf.env['LIBPATH_' + package]
    result = []
    for p in use[:use.index(package)]:
        for d in conf.env['INCLUDES_' + p]:
            if d not in own_dirs:
                result += [os.path.join(d, h) for h in headers if os.path.isfile(os.path.join(d, h))]
        for d in conf.env['LIBPATH_' + p]:
            if d not in own_dirs:
                for lib in libs:
                    result += sorted(glob.glob(os.path.join(d, 'lib' + lib + '.*')))
    return list(dict.fromkeys(result))


# use: packages in the same order as in the build
@Configure.conf
def check_package_versions(conf, use):
    errors = []
    for package in [p for p in use if p in PROBES]:
        conf.start_msg('Checking ' + package + ' version')
        expected = _probe(conf, package, [package] + [p for p in use if p != package])
        actual = _probe(conf, package, use)
        if not expected or not actual:
            conf.end_msg('skipped: failed to build or run', 'YELLOW')
            continue
        mismatches = ['  %s: %s, expected %s' % (kind, actual.get(kind), version)
                      for kind, version in sorted(expected.items()) if actual.get(kind) != version]
        if not mismatches:
            conf.end_msg(expected['header'])
            continue
        conf.end_msg('shadowed', 'RED')

        errors.append('\n'.join(
            ['The build sees a different %s than the one found in its own directories:' % package]
            + mismatches
            + ['Another %s installation shadows it:' % package]
            + ['  ' + path for path in _shadowing_files(conf, package, use)]
            + ['Install %s in a separate directory.' % package]))

    if errors:
        conf.fatal('\n\n'.join(errors))
