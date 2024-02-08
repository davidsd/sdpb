#! /usr/bin/env python
# encoding: utf-8

# -fsanitize=undefined checks for undefined behaviour in runtime
# This could prevent some weird UB errors, e.g. caused by -O3 compiler optimizations.
def configure(conf):
    cxx_fragment = '''
    int main() {
        int x = 1;
        return 0;
    }
    '''
    flag = '-fsanitize=undefined'
    try:
        conf.check_cxx(msg=f'Checking for {flag}',
                       fragment=cxx_fragment,
                       cxxflags=flag, linkflags=flag, uselib_store='fsanitize')
    except conf.errors.ConfigurationError as e:
        conf.to_log(f'{flag} not supported')
