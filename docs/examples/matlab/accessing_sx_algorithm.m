%
%     This file is part of CasADi.
%
%     CasADi -- A symbolic framework for dynamic optimization.
%     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
%                             K.U. Leuven. All rights reserved.
%     Copyright (C) 2011-2014 Greg Horn
%
%     CasADi is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 3 of the License, or (at your option) any later version.
%
%     CasADi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     Lesser General Public License for more details.
%
%     You should have received a copy of the GNU Lesser General Public
%     License along with CasADi; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
% -*- coding: utf-8 -*-

% Demonstration on how the algorithm of an MX function can be accessed and its operations can be transversed.

import casadi.*

% Create a function
a = SX.sym('a');
b = SX.sym('b',2);
f = Function('f', {a,b}, {2*a + b}, {'a', 'b'}, {'r'});

% Input values of the same dimensions as the above
input_val = {[2.0],
             [3.0;4.0]};

% Output values to be calculated of the same dimensions as the above
output_val = {zeros(2, 1)};

% Work vector
work = cell(1, f.sz_w());

% For debugging
instr = f.instructions_sx();

% Loop over the algorithm
for k=0:f.n_instructions()-1

  % Get the atomic operation
  op = f.instruction_id(k);

  o = f.instruction_output(k)+1;
  i = f.instruction_input(k)+1;

  if op==OP_CONST
    v = f.instruction_constant(k);
    work{o(1)} = v;
    disp(['work{' num2str(o(1)) '} = ' num2str(v) ';'])
  elseif op==OP_INPUT
    v = input_val{i(1)};
    work{o(1)} = v(i(2));
    disp(['work{' num2str(o(1)) '} = input{' num2str(i(1)) '}{' num2str(i(2)) '};     ---> ' num2str(v(i(2)))]);
  elseif op==OP_OUTPUT
    v = output_val{o(1)};
    v(o(2)) = work{i(1)};
    output_val{o(1)} = v;
    disp(['output{' num2str(o(1)) '}{' num2str(o(2)) '} = work{' num2str(i(1)) '};     ---> ' num2str(work{i(1)})]);
  else
    if op==OP_ADD
      work{o(1)} = work{i(1)} + work{i(2)};
      disp(['work{' num2str(o(1)) '} = work{' num2str(i(1)) '} + work{' num2str(i(2)) '};     ---> ' num2str(work{o(1)})]);
    elseif op==OP_MUL
      work{o(1)} = work{i(1)} * work{i(2)};
     disp(['work{' num2str(o(1)) '} = work{' num2str(i(1)) '} * work{'  num2str(i(2)) '};     ---> ' num2str(work{o(1)})]);
    else
      disp_in = {};
      for a=i
         disp_in{end+1} = ['work{' num2str(a) '}'];
      end
      debug_str = print_operator(instr{k+1},disp_in);
      error(['Unknown operation: '  num2str(op) ' -- ' debug_str]);
    end
  end
end
disp('------')
disp(['Evaluated ' str(f)])
disp('Expected: ')
celldisp(f.call(input_val))
disp('Got:      ')
celldisp(output_val)
