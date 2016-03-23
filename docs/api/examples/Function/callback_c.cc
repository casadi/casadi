#include "casadi/casadi.hpp"
using namespace casadi;

class MyCallback : public Callback {
public:
  // Creator function, creates an owning reference
  static Function create(const std::string& name, const Dict& opts=Dict()) {
    return Callback::create(name, new MyCallback(), opts);
  }

  // Number of inputs and outputs
  virtual int get_n_in() { return 1;}
  virtual int get_n_out() { return 1;}

  // Evaluate numerically
  virtual std::vector<DM> eval(const std::vector<DM>& arg) {
    DM x = arg.at(0);
    return {sin(x)};
  }
};

int main() {
  Function f = MyCallback::create("f");
  std::vector<DM> arg = {2};
  std::vector<DM> res = f(arg);
  std::cout << res << std::endl;
  return 0;
}
