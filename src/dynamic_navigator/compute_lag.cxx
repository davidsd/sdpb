#include

// SDP_solver 
// dual_residues[p] = c[p] - A[p,a,b] Y[a,b] - B[p,a] y[a]
// Lagrangian = dual_residues.x + b.y + Tr(X,Y) - mu log det (X) 


void compute_dual_residues_and_error(
  const Block_Info &block_info, const SDP &sdp, const Block_Vector &y,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  Block_Vector &dual_residues, El::BigFloat &dual_error, Timers &timers)



frobenius_product_symmetric(X, Y)


  // b.y
    objective
    = El::Dot(sdp.dual_objective_b, y.blocks.at(0)) + sdp.objective_const;

// temp = dB.y
        El::DistMatrix<El::BigFloat> temp(x.blocks.at(block));
        El::Zero(temp);
        El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
                 d_sdp.free_var_matrix.blocks[block], y.blocks.at(0),
                 El::BigFloat(0.0), temp);
 // -x.dB.y/
local_linear -= El::Dotu(temp, x.blocks.at(block));



