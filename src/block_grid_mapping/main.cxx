#include "../compute_block_grid_mapping.hxx"

#include <algorithm>
#include <iostream>

int main(int argc, char *argv[])
{
  if(argc < 3)
    {
      std::cerr << "Need 3 arguments\n";
      exit(1);
    }
  size_t num_procs(std::stoi(argv[1]));
  std::vector<Block_Cost> costs;
  for(int ii = 2; ii < argc; ++ii)
    {
      costs.emplace_back(std::stoi(argv[ii]), ii-2);
    }
  std::sort(costs.rbegin(),costs.rend());
  std::vector<Block_Map> mapping(compute_block_grid_mapping(num_procs, costs));

  for(auto &m: mapping)
    {
      std::cout << m.num_procs << ": " << m.cost/static_cast<double>(m.num_procs) << ", {";
      for(auto &index : m.block_indices)
        {
          std::cout << index << ",";
        }
      std::cout << "}\n";
    }
}
