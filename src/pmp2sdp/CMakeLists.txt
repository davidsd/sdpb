add_library(pmp2sdp_lib SHARED)
target_sources(pmp2sdp_lib
                PRIVATE
                Dual_Constraint_Group/Dual_Constraint_Group.cxx
                Dual_Constraint_Group/sample_bilinear_basis.cxx
                Output_SDP/Output_SDP.cxx
                write_block_data.cxx
                write_block_info_json.cxx
                write_objectives_json.cxx
                write_normalization_json.cxx
                write_sdp.cxx
                write_control_json.cxx
                Archive_Writer/Archive_Writer.cxx
                Archive_Writer/write_entry.cxx
                Archive_Entry.cxx)
target_link_libraries(pmp2sdp_lib PUBLIC ${USE_PACKAGES} pmp)