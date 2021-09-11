#pragma once

#include "../Mesh.hxx"

#include <deque>

std::deque<El::BigFloat>
get_zeros(const Mesh &mesh, const El::BigFloat &threshold);

