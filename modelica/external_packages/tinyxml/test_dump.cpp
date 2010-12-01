#include "dump.hpp"

// ----------------------------------------------------------------------
// main() for printing files named on the command line
// ----------------------------------------------------------------------
int main(int argc, char* argv[])
{
	for (int i=1; i<argc; i++)
	{
		dump_to_stdout(argv[i]);
	}
	return 0;
}
