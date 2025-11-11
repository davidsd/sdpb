#include "Memory_Tracker.hxx"

#include "ostream/pretty_print_bytes.hxx"
#include "sdpb_util/assert.hxx"

Memory_Tracker::Memory_Tracker(Memory_Tracker &&other) noexcept
    : root(std::move(other.root)),
      current_scope(other.current_scope),
      active_nodes(std::move(other.active_nodes))
{}
Memory_Tracker &Memory_Tracker::operator=(Memory_Tracker &&other) noexcept
{
  if(this == &other)
    return *this;
  root = std::move(other.root);
  current_scope = other.current_scope;
  active_nodes = std::move(other.active_nodes);
  return *this;
}
Memory_Tracker::Node::Node(Node *parent, const std::string &name,
                           const size_t bytes,
                           const std::function<void(size_t)> &on_peak_changed)
    : parent(parent), name(name), on_peak_changed(on_peak_changed)
{
  allocate(bytes);
}
Memory_Tracker::Node *Memory_Tracker::Node::add_child(
  const std::string &child_name, const size_t bytes,
  const std::function<void(size_t)> &on_child_peak_changed)
{
  return children
    .emplace_back(
      std::make_unique<Node>(this, child_name, bytes, on_child_peak_changed))
    .get();
}
void Memory_Tracker::Node::allocate(const size_t bytes)
{
  if(bytes == 0)
    return;
  for(auto *node = this; node != nullptr; node = node->parent)
    {
      node->curr_bytes += bytes;
      if(node->curr_bytes > node->peak_bytes)
        {
          node->peak_bytes = node->curr_bytes;
          if(node->on_peak_changed)
            node->on_peak_changed(node->peak_bytes);
        }
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
Memory_Tracker::Memory_Tracker(
  const std::string &name,
  const std::function<void(size_t)> &on_root_peak_changed)
    : root(nullptr, name, 0, on_root_peak_changed),
      current_scope(&root),
      active_nodes{current_scope}
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
void Memory_Tracker::print(std::ostream &os, const std::string &initial_indent,
                           const std::string &indent_string,
                           const bool also_print_exact_bytes) const
{
  constexpr size_t level = 0;
  print(os, root, initial_indent, indent_string, also_print_exact_bytes,
        level);
}
std::string Memory_Tracker::to_string(const std::string &initial_indent,
                                      const std::string &indent_string,
                                      bool also_print_exact_bytes) const
{
  std::stringstream ss;
  print(ss, initial_indent, indent_string, also_print_exact_bytes);
  return ss.str();
}
std::string Memory_Tracker::to_string(const std::string &indent_string,
                                      const bool also_print_exact_bytes) const
{
  return to_string("", indent_string, also_print_exact_bytes);
}
std::string Memory_Tracker::to_string(bool also_print_exact_bytes) const
{
  return to_string("  ", also_print_exact_bytes);
}
void Memory_Tracker::print(std::ostream &os, const Node &node,
                           const std::string &initial_indent,
                           const std::string &indent_string,
                           const bool also_print_exact_bytes,
                           const size_t indent_level)
{
  os << initial_indent;
  for(size_t level = 0; level < indent_level; ++level)
    os << indent_string;
  os << node.name << ": "
     << pretty_print_bytes(node.peak_bytes, also_print_exact_bytes)
     << std::endl;
  for(const auto &child : node.children)
    {
      print(os, *child, initial_indent, indent_string, also_print_exact_bytes,
            indent_level + 1);
    }
}
Memory_Tracker::Node *
Memory_Tracker::enter_scope(const std::string &name,
                            const std::function<void(size_t)> &on_peak_changed)
{
  return open(name, 0, current_scope, true, on_peak_changed);
}
void Memory_Tracker::exit_scope(Node *scope)
{
  close(scope, true);
}
Memory_Tracker::Node *
Memory_Tracker::allocate(const std::string &name, const size_t bytes,
                         Node *parent,
                         const std::function<void(size_t)> &on_peak_changed)
{
  return open(name, bytes, parent, false, on_peak_changed);
}
void Memory_Tracker::free(Node *node)
{
  return close(node, false);
}
Memory_Tracker::Node *
Memory_Tracker::open(const std::string &name, const size_t bytes, Node *parent,
                     const bool enter_scope,
                     const std::function<void(size_t)> &on_peak_changed)
{
  // if(!parent)
  //   parent = current_scope;
  ASSERT(parent != nullptr);
  auto *child = parent->add_child(name, bytes, on_peak_changed);
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

Memory_Tracker::Scope::Scope(Memory_Tracker &tracker, const std::string &name,
                             const std::function<void(size_t)> &on_peak_changed)
    : tracker(tracker), node(tracker.enter_scope(name, on_peak_changed))
{}
Memory_Tracker::Scope::~Scope() noexcept
{
  try
    {
      tracker.exit_scope(node);
    }
  catch(const std::exception &e)
    {
      PRINT_WARNING("Memory_Tracker::Scope::~Scope() exception: ", e.what());
    }
  catch(...)
    {
      PRINT_WARNING("Memory_Tracker::Scope::~Scope(): unknown exception");
    }
}

Memory_Tracker::Allocation::Allocation(
  Memory_Tracker &tracker, const std::string &name, const size_t bytes,
  const std::function<void(size_t)> &on_peak_changed)
    : tracker(tracker),
      node(
        tracker.allocate(name, bytes, tracker.current_scope, on_peak_changed))
{}
Memory_Tracker::Allocation::Allocation(
  Memory_Tracker &tracker, const Allocation &parent, const std::string &name,
  const size_t bytes, const std::function<void(size_t)> &on_peak_changed)
    : tracker(tracker),
      node(tracker.allocate(name, bytes, parent.node, on_peak_changed))
{}
Memory_Tracker::Allocation::~Allocation() noexcept
{
  try
    {
      tracker.free(node);
    }
  catch(const std::exception &e)
    {
      PRINT_WARNING("Memory_Tracker::Allocation::~Allocation() exception: ",
                    e.what());
    }
  catch(...)
    {
      PRINT_WARNING(
        "Memory_Tracker::Allocation::~Allocation(): unknown exception");
    }
}
Memory_Tracker::Group::Group(Memory_Tracker &tracker, const std::string &name,
                             const std::function<void(size_t)> &on_peak_changed)
    : Allocation(tracker, name, 0, on_peak_changed)
{}
Memory_Tracker::Group::Group(Memory_Tracker &tracker, const Allocation &parent,
                             const std::string &name,
                             const std::function<void(size_t)> &on_peak_changed)
    : Allocation(tracker, parent, name, 0, on_peak_changed)
{}
