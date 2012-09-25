/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2011 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file include/qpOASES/extras/SolutionAnalysis.hpp
 *	\author Boris Houska, Hans Joachim Ferreau
 *	\version 3.0beta
 *	\date 2008-2009
 *
 *	Declaration of the SolutionAnalysis class designed to perform
 *	additional analysis after solving a QP with qpOASES.
 */


#ifndef QPOASES_SOLUTIONANALYSIS_HPP
#define QPOASES_SOLUTIONANALYSIS_HPP


#include <qpOASES/SQProblem.hpp>


BEGIN_NAMESPACE_QPOASES


/** 
 *	\brief Provides additional tools for analysing QP solutions.
 *
 *	This class is intended to provide additional tools for analysing
 *  a QP solution obtained with qpOASES.
 *
 *	\author Boris Houska, Hans Joachim Ferreau
 *	\version 3.0beta
 *	\date 2007-2011
 */
class SolutionAnalysis
{
	/*
	 *	PUBLIC MEMBER FUNCTIONS
	 */
	public:
		/** Default constructor. */
		SolutionAnalysis( );

		/** Copy constructor (deep copy). */
		SolutionAnalysis(	const SolutionAnalysis& rhs		/**< Rhs object. */
							);

		/** Destructor. */
		~SolutionAnalysis( );

		/** Copy asingment operator (deep copy). */
		SolutionAnalysis& operator=(	const SolutionAnalysis& rhs		/**< Rhs object. */
										);


		/** Determines the maximum violation of the KKT optimality conditions
		 *  of the current iterate within the QProblemB object.
		 *	\return SUCCESSFUL_RETURN \n
					RET_UNABLE_TO_ANALYSE_QPROBLEM */
		returnValue getMaxKKTviolation(	QProblemB* qp,			/**< QProblemB to be analysed. */
										real_t& maxKKTviolation	/**< OUTPUT: maximum violation of the KKT conditions. */
										) const;

		/** Determines the maximum violation of the KKT optimality conditions
		 *  of the current iterate within the QProblem object.
		 *	\return SUCCESSFUL_RETURN \n
					RET_UNABLE_TO_ANALYSE_QPROBLEM */
		returnValue getMaxKKTviolation(	QProblem* qp,			/**< QProblem to be analysed. */
										real_t& maxKKTviolation	/**< OUTPUT: maximum violation of the KKT conditions. */
										) const;

		/** Determines the maximum violation of the KKT optimality conditions
		 *  of the current iterate within the SQProblem object.
		 *	\return SUCCESSFUL_RETURN \n
					RET_UNABLE_TO_ANALYSE_QPROBLEM */
		returnValue getMaxKKTviolation(	SQProblem* qp,			/**< SQProblem to be analysed. */
										real_t& maxKKTviolation	/**< OUTPUT: maximum violation of the KKT conditions. */
										) const;


		/** Computes the variance-covariance matrix of the QP output for uncertain	\n
			inputs.
		 *	\return SUCCESSFUL_RETURN \n
					RET_HOTSTART_FAILED \n
		 			RET_STEPDIRECTION_FAILED_TQ \n
					RET_STEPDIRECTION_FAILED_CHOLESKY */
		returnValue getVarianceCovariance(	QProblemB* qp,			/**< QProblemB to be analysed. */
											real_t* g_b_bA_VAR,		/**< INPUT : Variance-covariance of g, the bounds lb and ub, and lbA and ubA respectively. Dimension: 2nV x 2nV */
											real_t* Primal_Dual_VAR	/**< OUTPUT: The result for the variance-covariance of the primal and dual variables. Dimension: 2nV x 2nV */
											) const;

		/** Computes the variance-covariance matrix of the QP output for uncertain	\n
			inputs.
		 *	\return SUCCESSFUL_RETURN \n
					RET_HOTSTART_FAILED \n
		 			RET_STEPDIRECTION_FAILED_TQ \n
					RET_STEPDIRECTION_FAILED_CHOLESKY */
		returnValue getVarianceCovariance(	QProblem* qp,			/**< QProblem to be analysed. */
											real_t* g_b_bA_VAR,		/**< INPUT : Variance-covariance of g, the bounds lb and ub, and lbA and ubA respectively. Dimension:  (2nV+nC) x (2nV+nC) */
											real_t* Primal_Dual_VAR	/**< OUTPUT: The result for the variance-covariance of the primal and dual variables. Dimension:  (2nV+nC) x (2nV+nC) */
											) const;

		/** Computes the variance-covariance matrix of the QP output for uncertain	\n
			inputs.
		 *	\return SUCCESSFUL_RETURN \n
					RET_HOTSTART_FAILED \n
		 			RET_STEPDIRECTION_FAILED_TQ \n
					RET_STEPDIRECTION_FAILED_CHOLESKY */
		returnValue getVarianceCovariance(	SQProblem* qp,			/**< SQProblem to be analysed. */
											real_t* g_b_bA_VAR,		/**< INPUT : Variance-covariance of g, the bounds lb and ub, and lbA and ubA respectively. Dimension:  (2nV+nC) x (2nV+nC) */
											real_t* Primal_Dual_VAR	/**< OUTPUT: The result for the variance-covariance of the primal and dual variables. Dimension:  (2nV+nC) x (2nV+nC) */
											) const;


	/*
	 *	PROTECTED MEMBER VARIABLES
	 */
	protected:

};


END_NAMESPACE_QPOASES

#include <qpOASES/extras/SolutionAnalysis.ipp>

#endif	/* QPOASES_SOLUTIONANALYSIS_HPP */


/*
 *	end of file
 */
