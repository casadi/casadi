#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
# 
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
# 
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
# 
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# 
# 
# -*- coding: utf-8 -*-

from casadi_main import *
from casadi_primitive import *
from casadi_primitive_tools import *
from casadi_noncore import *

interfaces = []
failed_interfaces = {}

try:
  from casadi_interface_ipopt import *
  interfaces.append("ipopt")
except Exception as e:
  failed_interfaces["ipopt"] = str(e)


try:
  from casadi_interface_sundials import *
  interfaces.append("sundials")
except Exception as e:
  failed_interfaces["sundials"] = str(e)


try:
  from casadi_interface_qpoases import *
  interfaces.append("qpoases")
except Exception as e:
  failed_interfaces["qpoases"] = str(e)


try:
  from casadi_interface_dsdp import *
  interfaces.append("dsdp")
except Exception as e:
  failed_interfaces["dsdp"] = str(e)


try:
  from casadi_interface_csparse import *
  interfaces.append("csparse")
except Exception as e:
  failed_interfaces["csparse"] = str(e)


try:
  from casadi_interface_knitro import *
  interfaces.append("knitro")
except Exception as e:
  failed_interfaces["knitro"] = str(e)


try:
  from casadi_interface_cplex import *
  interfaces.append("cplex")
except Exception as e:
  failed_interfaces["cplex"] = str(e)


try:
  from casadi_interface_ooqp import *
  interfaces.append("ooqp")
except Exception as e:
  failed_interfaces["ooqp"] = str(e)


try:
  from casadi_interface_slicot import *
  interfaces.append("slicot")
except Exception as e:
  failed_interfaces["slicot"] = str(e)


try:
  from casadi_interface_worhp import *
  interfaces.append("worhp")
except Exception as e:
  failed_interfaces["worhp"] = str(e)


try:
  from casadi_interface_snopt import *
  interfaces.append("snopt")
except Exception as e:
  failed_interfaces["snopt"] = str(e)


try:
  from casadi_interface_lapack import *
  interfaces.append("lapack")
except Exception as e:
  failed_interfaces["lapack"] = str(e)


