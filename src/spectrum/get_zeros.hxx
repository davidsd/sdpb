#pragma once

#include "sdpb_util/Mesh.hxx"

std::vector<El::BigFloat>
get_zeros(const Mesh &mesh, const El::BigFloat &threshold);
