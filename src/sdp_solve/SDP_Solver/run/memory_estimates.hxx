#pragma once

#include "sdp_solve/Block_Diagonal_Matrix.hxx"
#include "sdp_solve/Block_Info.hxx"
#include "sdp_solve/SDP.hxx"
#include "sdpb_util/Environment.hxx"

size_t
get_max_shared_memory_bytes(size_t default_max_shared_memory_bytes,
                            size_t nonshared_memory_required_per_node_bytes,
                            const Environment &env, bool debug);

size_t get_matrix_size_local(const Block_Diagonal_Matrix &X);
size_t get_A_X_size_local(const Block_Info &block_info, const SDP &sdp);
size_t get_schur_complement_size_local(const Block_Info &block_info);
size_t get_B_size_local(const SDP &sdp);
size_t get_Q_size_local(const SDP &sdp);
size_t get_SDP_size_local(const SDP &sdp);

size_t bigfloat_bytes();

