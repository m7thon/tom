#include "tom.h"

int main() {
  auto oom = tom::Oom(3,5);
  std::cout << oom.toJSON() << std::endl;
  return 1;
}
