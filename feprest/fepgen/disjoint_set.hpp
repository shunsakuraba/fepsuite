#include <vector>

class disjoint_set
{
public:
  disjoint_set(int n)
  : link(n)
  {
    for(int i = 0; i < n; ++i) {
      link[i] = i;
    }
  }

  void join(int i, int j)
  {
    link[find_parent(i)] = find_parent(j);
  }

  int find_parent(int i)
  {
    if(link[i] == i) return i;
    return link[i] = find_parent(link[i]);
  }

private:
  std::vector<int> link;
};
