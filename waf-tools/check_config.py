#! /usr/bin/env python
# encoding: utf-8

def check_pkg_config(conf, package, uselib_store, fragment, use=None, msg_prefix=None, **kwargs):
    msg_prefix = msg_prefix or ''
    if not conf.check_cfg(msg=msg_prefix + 'Searching for ' + package + ' in pkg-config',
                          path='pkg-config',
                          args='--libs --cflags',
                          package=package,
                          uselib_store=uselib_store,
                          **kwargs):
        return False

    # check_cfg() does not support 'fragment' argument,
    # so we have to check compilation separately via check_cxx()

    # We also set rpath=libpath, same as we do in check_config() below
    # You can find LIBPATH values (set by conf.check_cfg)
    # in build/c4che/_cache.py, e.g.
    # LIBPATH_cblas = ['/usr/lib/x86_64-linux-gnu/openblas-pthread/']
    libpath = getattr(conf.env, 'LIBPATH_' + uselib_store, None)
    return conf.check_cxx(msg=msg_prefix + 'Checking ' + package + ' compilation',
                          fragment=fragment,
                          rpath=libpath,
                          use=(use or []) + [uselib_store],
                          uselib_store=uselib_store,
                          **kwargs)


def check_config(conf, fragment, uselib_store, includes, libpath, lib, msg=None, packages=None, use=None,
                 mandatory=True,
                 **kwargs):
    # pkg-config package names
    if packages is None:
        packages = [uselib_store]
    # extra dependencies, e.g. cxx17
    if use is None:
        use = []

    if conf.check_cxx(msg=msg or 'Checking for ' + uselib_store,
                      fragment=fragment,
                      includes=includes,
                      uselib_store=uselib_store,
                      libpath=libpath,
                      rpath=libpath,
                      lib=lib,
                      use=use,
                      mandatory=False,
                      **kwargs):
        return True

    for package in packages:
        if check_pkg_config(conf,
                            fragment=fragment,
                            package=package,
                            uselib_store=uselib_store,
                            use=use,
                            mandatory=False,
                            msg_prefix='  ',
                            **kwargs):
            return True

    if mandatory:
        conf.fatal('Could not find ' + uselib_store)
    return False
