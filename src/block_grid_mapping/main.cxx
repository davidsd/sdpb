#include "../compute_block_grid_mapping.hxx"

#include <algorithm>
#include <iostream>

int main(int argc, char *argv[])
{
  if(argc < 4)
    {
      std::cerr
        << "Need at least 3 arguments: procs_per_node, num_nodes, costs...\n";
      exit(1);
    }
  size_t num_procs(std::stoi(argv[1])), procs_per_node(std::stoi(argv[2]));
  std::vector<Block_Cost> costs;
  for(int ii = 3; ii < argc; ++ii)
    {
      costs.emplace_back(std::stoi(argv[ii]), ii - 3);
    }
  std::sort(costs.rbegin(), costs.rend());
  std::vector<std::vector<Block_Map>> mapping(
    compute_block_grid_mapping(num_procs, procs_per_node, costs));

  for(size_t node = 0; node < mapping.size(); ++node)
    {
      for(auto &m : mapping[node])
        {
          std::cout << node << " " << m.num_procs << ": "
                    << m.cost / static_cast<double>(m.num_procs) << ", {";
          for(size_t ii = 0; ii < m.block_indices.size(); ++ii)
            {
              if(ii != 0)
                {
                  std::cout << ",";
                }
              std::cout << m.block_indices[ii];
            }
          std::cout << "}\n";
        }
      std::cout << "\n";
    }
}
