//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

//
// Brancher class finds suitable candidates for branching. e.g. integer
// variable, 'non convex' variable.
//

#include "MinotaurConfig.h"
#include "Brancher.h"
#include "Node.h"

using namespace Minotaur;


Brancher::Brancher()
  : logger_(LoggerPtr()) // NULL
{
}


Brancher::~Brancher()
{
}


void Brancher::updateAfterLP(NodePtr , ConstSolutionPtr )
{
}


// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
