#pragma once

#include "sdpb_util/assert.hxx"

#include <memory>
#include <string>
#include <vector>

// Create tree of nodes describing memory allocations and their lifetimes,
// track peak memory usage.
// Usage example:
// Memory_Tracker tracker("Total memory estimate");
// {
//   // Allocate 100 bytes
//   Allocation buffer("buffer", 100, tracker);
//   {
//     // Short-living allocations inside foo() scope
//     Scope foo("foo()", tracker);
//     Allocation x("x", 200, tracker);
//     Allocation y("y", 300, tracker);
//   }
//   {
//     // Short-living allocations inside bar() scope
//     Scope bar("bar()", tracker);
//     Allocation z("z", 400, tracker);
//     Allocation w("w", 500, tracker);
//   }
// }
// bool also_print_exact_bytes = false;
// // Will print peak memory usage at each step
// tracker.print(std::cout, also_print_exact_bytes);
// // Will return 100 + 400 + 500 = 1000
// size_t max_mem = tracker.peak_memory();

class Memory_Tracker
{
  // Helper classes

  struct Node
  {
    Node *parent;
    std::string name;
    size_t curr_bytes = 0;
    size_t peak_bytes = 0;
    std::vector<std::unique_ptr<Node>> children;

    Node(Node *parent, const std::string &name, size_t bytes);

    Node *add_child(const std::string &child_name, size_t bytes);

    void allocate(size_t bytes);
    void free(size_t bytes);
    void exit();
  };

public:
  // RAII interface for new scope
  struct [[nodiscard]] Scope
  {
    Scope(const std::string &name, Memory_Tracker &tracker);
    ~Scope();

  private:
    Memory_Tracker &tracker;
    Node *node;
  };

  // RAII interface for memory allocation within current scope
  struct [[nodiscard]] Allocation
  {
    Allocation(const std::string &name, size_t bytes, Memory_Tracker &tracker);
    Allocation(const std::string &name, size_t bytes, Memory_Tracker &tracker,
               const Allocation &parent);
    ~Allocation();

  private:
    Memory_Tracker &tracker;
    Node *node;
  };

  // Memory_Tracker members

private:
  Node root;
  Node *current_scope;
  std::vector<Node *> active_nodes;

  // Memory_Tracker public API

public:
  explicit Memory_Tracker(const std::string &name);

  // Peak memory usage, bytes
  size_t peak_memory() const;
  bool is_finished() const;
  void print(std::ostream &os, bool also_print_exact_bytes) const;

  // Memory_Tracker implementation

private:
  static void print(std::ostream &os, const Node &node, size_t level,
                    bool also_print_exact_bytes);

  // Enter new named scope
  Node *enter_scope(const std::string &name);
  // Exit named scope
  void exit_scope(Node *scope);
  // Make an allocation in the current scope
  Node *allocate(const std::string &name, size_t bytes, Node *parent);
  // Free an allocation in the current scope
  void free(Node *node);

  // implementation for enter_scope() and allocate()
  Node *
  open(const std::string &name, size_t bytes, Node *parent, bool enter_scope);
  // implementation for exit_scope() and free()
  void close(Node *node, bool exit_scope);
};
