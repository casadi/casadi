#include "casadi/casadi.hpp"

using namespace casadi;
using namespace std;

class MyCallback : public Callback {
public:
  // Creator function, creates an owning reference
  static Function fun(const string& name, const Dict& opts=Dict()) {
    return Callback::fun(name, new MyCallback(), opts);
  }

  // Number of inputs and outputs
  virtual int get_n_in() { return 1;}
  virtual int get_n_out() { return 1;}

  // Evaluate numerically
  virtual vector<DM> eval(const vector<DM>& arg) {
    vector<DM> ret(1);
    ret.at(0) = sin(arg.at(0));
    return ret;
  }
};

int main(int argc, char *argv[]) {

  Function f = MyCallback::fun("f");

  DM arg = DM(2);
  DM res = f(vector<DM>(1,arg)).at(0);

  cout << "out:" << res << endl;

  return 0;
}
