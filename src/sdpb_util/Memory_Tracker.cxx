#include "Memory_Tracker.hxx"

#include "ostream/pretty_print_bytes.hxx"
#include "sdpb_util/assert.hxx"

Memory_Tracker::Node::Node(Node *parent, const std::string &name,
                           const size_t bytes)
    : parent(parent), name(name)
{
  allocate(bytes);
}
Memory_Tracker::Node *
Memory_Tracker::Node::add_child(const std::string &child_name,
                                const size_t bytes)
{
  return children.emplace_back(std::make_unique<Node>(this, child_name, bytes))
    .get();
}
void Memory_Tracker::Node::allocate(const size_t bytes)
{
  if(bytes == 0)
    return;
  for(auto *node = this; node != nullptr; node = node->parent)
    {
      node->curr_bytes += bytes;
      node->peak_bytes = std::max(node->curr_bytes, node->peak_bytes);
    }
}
void Memory_Tracker::Node::free(const size_t bytes)
{
  if(bytes == 0)
    return;
  for(auto *node = this; node != nullptr; node = node->parent)
    {
      ASSERT(node->curr_bytes >= bytes, "Cannot deallocate: too many bytes",
             DEBUG_STRING(node->name), DEBUG_STRING(node->curr_bytes),
             DEBUG_STRING(bytes));
      node->curr_bytes -= bytes;
    }
}
void Memory_Tracker::Node::exit()
{
  for(const auto &child : children)
    {
      ASSERT_EQUAL(child->curr_bytes, 0,
                   "Cannot exit node: child is still alive",
                   DEBUG_STRING(name), DEBUG_STRING(child->name));
    }
  free(curr_bytes);
}
Memory_Tracker::Memory_Tracker(const std::string &name)
    : root(nullptr, name, 0), current_scope(&root), active_nodes{current_scope}
{}
size_t Memory_Tracker::peak_memory() const
{
  return root.peak_bytes;
}
bool Memory_Tracker::is_finished() const
{
  if(current_scope != &root)
    return false;
  if(active_nodes.size() > 1)
    return false;
  ASSERT_EQUAL(active_nodes.size(), 1);
  ASSERT_EQUAL(active_nodes.front(), &root,
               DEBUG_STRING(active_nodes.front()->name));
  return true;
}
void Memory_Tracker::print(std::ostream &os,
                           const bool also_print_exact_bytes) const
{
  constexpr size_t level = 0;
  print(os, root, level, also_print_exact_bytes);
}
void Memory_Tracker::print(std::ostream &os, const Node &node,
                           const size_t level,
                           const bool also_print_exact_bytes)
{
  // indent = 2 spaces
  const std::string indent(level * 2, ' ');
  os << indent << node.name << ": "
     << pretty_print_bytes(node.peak_bytes, also_print_exact_bytes)
     << std::endl;
  for(const auto &child : node.children)
    {
      print(os, *child, level + 1, also_print_exact_bytes);
    }
}
Memory_Tracker::Node *Memory_Tracker::enter_scope(const std::string &name)
{
  return open(name, 0, true);
}
void Memory_Tracker::exit_scope(Node *scope)
{
  close(scope, true);
}
Memory_Tracker::Node *
Memory_Tracker::allocate(const std::string &name, size_t bytes)
{
  return open(name, bytes, false);
}
void Memory_Tracker::free(Node *node)
{
  return close(node, false);
}
Memory_Tracker::Node *
Memory_Tracker::open(const std::string &name, const size_t bytes,
                     const bool enter_scope)
{
  auto *child = current_scope->add_child(name, bytes);
  active_nodes.push_back(child);
  if(enter_scope)
    current_scope = child;
  return child;
}
void Memory_Tracker::close(Node *node, const bool exit_scope)
{
  ASSERT(node != nullptr);
  ASSERT(node != &root);
  ASSERT(!active_nodes.empty());
  ASSERT_EQUAL(node, active_nodes.back(), "Closing nodes in wrong order!",
               DEBUG_STRING(node->name),
               DEBUG_STRING(active_nodes.back()->name));
  node->exit();
  active_nodes.pop_back();
  if(exit_scope)
    {
      ASSERT(current_scope != &root);
      ASSERT(current_scope->parent != nullptr);
      current_scope->exit();
      current_scope = current_scope->parent;
    }
}
Memory_Tracker::Scope::Scope(const std::string &name, Memory_Tracker &tracker)
    : tracker(tracker), node(tracker.enter_scope(name))
{}
Memory_Tracker::Scope::~Scope()
{
  tracker.exit_scope(node);
}
Memory_Tracker::Allocation::Allocation(const std::string &name,
                                       const size_t bytes,
                                       Memory_Tracker &tracker)
    : tracker(tracker), node(tracker.allocate(name, bytes))
{}
Memory_Tracker::Allocation::~Allocation()
{
  tracker.free(node);
}
