#ifdef SWIGPYTHON
#ifdef WITH_PYTHON_INTERRUPTS
%{
#include <pythonrun.h>

void SigIntHandler(int) {
  std::cerr << "Keyboard Interrupt" << std::endl;
	signal(SIGINT, SIG_DFL);
	kill(getpid(), SIGINT);
}
%}

%init %{

PyOS_setsig(SIGINT, SigIntHandler);
%}
#endif // WITH_PYTHON_INTERRUPTS
#endif // SWIGPYTHON

