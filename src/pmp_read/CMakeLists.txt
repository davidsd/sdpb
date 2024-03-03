add_library(pmp_read SHARED)
target_sources(pmp_read
                PRIVATE
                collect_files_expanding_nsv.cxx
                PMP_File_Parse_Result.cxx
                read_nsv_file_list.cxx
                read_polynomial_matrix_program.cxx
                read_json/read_json.cxx
                read_json/Json_PMP_Parser.cxx
                read_mathematica/read_mathematica.cxx
                read_mathematica/parse_SDP/parse_SDP.cxx
                read_mathematica/parse_SDP/parse_matrices.cxx
                read_mathematica/parse_SDP/parse_number.cxx
                read_mathematica/parse_SDP/parse_polynomial.cxx
                read_mathematica/parse_SDP/parse_matrix/parse_matrix.cxx
                read_mathematica/parse_SDP/parse_matrix/parse_damped_rational.cxx
                read_xml/read_xml.cxx
                read_xml/Xml_Parser/on_start_element.cxx
                read_xml/Xml_Parser/on_end_element.cxx
                read_xml/Xml_Parser/on_characters.cxx)
target_link_libraries(pmp_read PUBLIC ${USE_PACKAGES} pmp pmp2sdp_lib)