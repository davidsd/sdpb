#pragma once

#include "sdpb_util/assert.hxx"

#include <boost/noncopyable.hpp>

#include <functional>
#include <memory>
#include <string>
#include <vector>

// Create tree of nodes describing memory allocations and their lifetimes,
// track peak memory usage.
// Usage example:
// Memory_Tracker tracker("Total memory estimate");
// {
//   // Allocate 100 bytes
//   Allocation buffer(tracker, "buffer", 100);
//   {
//     // Short-living allocations inside foo() scope
//     Scope foo(tracker, "foo()");
//     Allocation x(tracker, "x", 200);
//     Allocation y(tracker, "y", 300);
//   }
//   {
//     // Short-living allocations inside bar() scope
//     Scope bar(tracker, "bar()");
//     Group z_and_w(tracker, "z and w");
//     Allocation z(tracker, z_and_w, "z", 400);
//     Allocation w(tracker, z_and_w, "w", 500);
//   }
// }
// std::string initial_indent = "";
// std::string indent = "\t";
// bool also_print_exact_bytes = false;
// // Will print peak memory usage at each step
// tracker.print(std::cout, initial_indent, indent_string, also_print_exact_bytes);
// // Will return 100 + 400 + 500 = 1000
// size_t max_mem = tracker.peak_memory();

class Memory_Tracker
{
public:
  Memory_Tracker(const Memory_Tracker &other) = delete;
  Memory_Tracker(Memory_Tracker &&other) noexcept;
  Memory_Tracker &operator=(const Memory_Tracker &other) = delete;
  Memory_Tracker &operator=(Memory_Tracker &&other) noexcept;

private:
  // Helper classes

  struct Node
  {
    Node *parent;
    std::string name;
    size_t curr_bytes = 0;
    size_t peak_bytes = 0;
    std::vector<std::unique_ptr<Node>> children;
    std::function<void(size_t)> on_peak_changed = {};

    Node(Node *parent, const std::string &name, size_t bytes,
         const std::function<void(size_t)> &on_peak_changed);

    Node *add_child(const std::string &child_name, size_t bytes,
                    const std::function<void(size_t)> &on_child_peak_changed);

    void allocate(size_t bytes);
    void free(size_t bytes);
    void exit();
  };

public:
  // RAII interface for new scope
  struct [[nodiscard]] Scope : boost::noncopyable
  {
    Scope(Memory_Tracker &tracker, const std::string &name,
          const std::function<void(size_t)> &on_peak_changed = {});
    ~Scope();

  private:
    Memory_Tracker &tracker;
    Node *node;
  };

  // RAII interface for memory allocation within current scope
  struct [[nodiscard]] Allocation : boost::noncopyable
  {
    Allocation(Memory_Tracker &tracker, const std::string &name, size_t bytes,
               const std::function<void(size_t)> &on_peak_changed = {});
    Allocation(Memory_Tracker &tracker, const Allocation &parent,
               const std::string &name, size_t bytes,
               const std::function<void(size_t)> &on_peak_changed = {});
    ~Allocation();

  private:
    Memory_Tracker &tracker;
    Node *node;
  };
  // This is just an empty Allocation (own bytes=0),
  // intended to be used as a parent for other allocations.
  struct [[nodiscard]] Group : Allocation
  {
    Group(Memory_Tracker &tracker, const std::string &name,
          const std::function<void(size_t)> &on_peak_changed = {});
    Group(Memory_Tracker &tracker, const Allocation &parent,
          const std::string &name,
          const std::function<void(size_t)> &on_peak_changed = {});
  };

  // Memory_Tracker members

private:
  Node root;
  Node *current_scope;
  std::vector<Node *> active_nodes;

  // Memory_Tracker public API

public:
  explicit Memory_Tracker(
    const std::string &name,
    const std::function<void(size_t)> &on_root_peak_changed = {});

  // Peak memory usage, bytes
  [[nodiscard]] size_t peak_memory() const;
  [[nodiscard]] bool is_finished() const;
  void
  print(std::ostream &os, const std::string &initial_indent,
        const std::string &indent_string, bool also_print_exact_bytes) const;
  [[nodiscard]] std::string
  to_string(const std::string &initial_indent,const std::string &indent_string,
            bool also_print_exact_bytes) const;
  [[nodiscard]] std::string
  to_string(const std::string &indent_string = "  ",
            bool also_print_exact_bytes = true) const;
  [[nodiscard]] std::string
  to_string(bool also_print_exact_bytes) const;

  // Memory_Tracker implementation

private:
  static void
  print(std::ostream &os, const Node &node, const std::string &initial_indent,
        const std::string &indent_string, bool also_print_exact_bytes,
        size_t indent_level);

  // Enter new named scope
  Node *enter_scope(const std::string &name,
                    const std::function<void(size_t)> &on_peak_changed);
  // Exit named scope
  void exit_scope(Node *scope);
  // Make an allocation in the current scope
  Node *allocate(const std::string &name, size_t bytes, Node *parent,
                 const std::function<void(size_t)> &on_peak_changed);
  // Free an allocation in the current scope
  void free(Node *node);

  // implementation for enter_scope() and allocate()
  Node *
  open(const std::string &name, size_t bytes, Node *parent, bool enter_scope,
       const std::function<void(size_t)> &on_peak_changed);
  // implementation for exit_scope() and free()
  void close(Node *node, bool exit_scope);
};
