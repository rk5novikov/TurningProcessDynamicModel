	Summary:
Turning machining dynamic simulation based on 2DOF tool model

	Desription:
This thesis project is devoted to the study of the dynamics of cutting tool in the process of turning based on a single-mass model with two degrees of freedom. As part of the work, a program was written that allows obtain time realizations of the process of tool vibrations and cutting forces during turning, as well as to build a geometric model of the workpiece surface after machining. The developed program makes it possible to integrate the equations of tool vibrations for different machining parameters, with different cutting edge geometry, and with different models of cutting forces. In the course of the work, the influence of the cutting depth on tool vibrations was studied, and the main vibration modes were obtained with simulation and discussed.


	Files:
calc_*** - functions for solving differential equation
geom_*** - functions for creating and working with geom. models.
cut - material cutting function.
main_*** - programs for process modeling.
main_const - program for modeling with forces constant at a step by time.
main_lin - program for simulating linear forces at a step by time (without iterative refinement).
main_iter - a program for simulating linear forces at a step by time (with iterative refinement).
main_iter_DOF1 - the same, but with only one degree of freedom.
main_iter_puankaracii - a function that returns the resulting simulation 
	dynamic displacements, which are then used to construct punching lines.
main_iter_puankaracii_8 - the same, but with a cut depth of cut = 8 feeds.
main_iter_puankaracii_1DOF - the same, but with one degree of freedom.
main_iter_puankaracii_1DOF_8 - the same, but with a depth of cut of 8 feeds.
puank - is the main program for building Puankare plots.
