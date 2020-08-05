import os, subprocess

def options(opt):
    opt.load(['compiler_cxx','gnu_dirs','cxx17','boost','gmpxx','mpfr',
              'elemental','libxml2', 'rapidjson'])

def configure(conf):
    if not 'CXX' in os.environ or os.environ['CXX']=='g++' or os.environ['CXX']=='icpc':
        conf.environ['CXX']='mpicxx'

    conf.load(['compiler_cxx','gnu_dirs','cxx17','boost','gmpxx','mpfr',
               'elemental','libxml2', 'rapidjson'])

    conf.env.git_version=subprocess.check_output('git describe --dirty', universal_newlines=True, shell=True).rstrip()
    
def build(bld):
    default_flags=['-Wall', '-Wextra', '-O3', '-D SDPB_VERSION_STRING="' + bld.env.git_version + '"']
    # default_flags=['-Wall', '-Wextra', '-g', '-D SDPB_VERSION_STRING="' + bld.env.git_version + '"']
    use_packages=['cxx17','boost','gmpxx','mpfr','elemental','libxml2', 'rapidjson']
    
    sdp_read_sources=['src/sdp_read/read_input/read_input.cxx',
                      'src/sdp_read/read_input/read_json/read_json.cxx',
                      'src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_key.cxx',
                      'src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_string.cxx',
                      'src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_start_array.cxx',
                      'src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_end_array.cxx',
                      'src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_start_object.cxx',
                      'src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_end_object.cxx',
                      'src/sdp_read/read_input/read_json/Damped_Rational_State/json_key.cxx',
                      'src/sdp_read/read_input/read_json/Damped_Rational_State/json_string.cxx',
                      'src/sdp_read/read_input/read_json/Damped_Rational_State/json_start_array.cxx',
                      'src/sdp_read/read_input/read_json/Damped_Rational_State/json_end_array.cxx',
                      'src/sdp_read/read_input/read_json/Damped_Rational_State/json_start_object.cxx',
                      'src/sdp_read/read_input/read_json/Damped_Rational_State/json_end_object.cxx',
                      'src/sdp_read/read_input/read_json/JSON_Parser/Key.cxx',
                      'src/sdp_read/read_input/read_json/JSON_Parser/String.cxx',
                      'src/sdp_read/read_input/read_json/JSON_Parser/StartArray.cxx',
                      'src/sdp_read/read_input/read_json/JSON_Parser/EndArray.cxx',
                      'src/sdp_read/read_input/read_json/JSON_Parser/StartObject.cxx',
                      'src/sdp_read/read_input/read_json/JSON_Parser/EndObject.cxx',
                      'src/sdp_read/read_input/read_mathematica/read_mathematica.cxx',
                      'src/sdp_read/read_input/read_mathematica/parse_SDP/parse_SDP.cxx',
                      'src/sdp_read/read_input/read_mathematica/parse_SDP/parse_matrices.cxx',
                      'src/sdp_read/read_input/read_mathematica/parse_SDP/parse_number.cxx',
                      'src/sdp_read/read_input/read_mathematica/parse_SDP/parse_polynomial.cxx',
                      'src/sdp_read/read_input/read_mathematica/parse_SDP/parse_matrix/parse_matrix.cxx',
                      'src/sdp_read/read_input/read_mathematica/parse_SDP/parse_matrix/parse_damped_rational.cxx']

    bld.stlib(source=sdp_read_sources,
              target='sdp_read',
              cxxflags=default_flags,
              use=use_packages)

    bld.program(source=['src/outer/main.cxx',
                        'src/outer/is_feasible.cxx',
                        'src/outer/compute_optimal.cxx',
                        'src/outer/solve_LP.cxx',
                        'src/outer/get_new_points.cxx',
                        'src/outer/load_vector.cxx',
                        'src/outer/Mesh/Mesh.cxx',
                        'src/outer/Mesh/ostream.cxx',
                        'src/outer/Functional/Functional.cxx',
                        'src/outer/Functional/prefactor.cxx'],
                target='outer',
                cxxflags=default_flags,
                use=use_packages
                )
