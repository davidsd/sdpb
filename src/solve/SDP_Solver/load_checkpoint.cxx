#include "../SDP_Solver.hxx"

#include <boost/archive/text_iarchive.hpp>

void SDP_Solver::load_checkpoint(const boost::filesystem::path &checkpoint_file)
{
  boost::filesystem::ifstream ifs(checkpoint_file);
  boost::archive::text_iarchive ar(ifs);
  if(El::mpi::Rank() == 0)
    {
      std::cout << "Loading checkpoint from : " << checkpoint_file << '\n';
    }
  // boost::serialization::serialize_SDP_solver_state(ar, x_elemental, X,
  //                                                  y_elemental, Y);
}
