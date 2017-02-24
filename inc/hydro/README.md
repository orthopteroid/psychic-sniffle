# Hydropower Optimization with GA: Theory of operation
 
  Starting from some initial conditions, the task is to find a sequence of changes to the operation of a
  river system's units that meet the criteria of feasibility (ie. realistic operations) and provide either
  some maximum benefit (ie. [profit] = [power] x [price]) or some minimum cost (ie. environmental impact,
  or cost of extra purchased electricity from fossil sources).
 
  The states of the units and the plants are separate from the "operations" performed on the states of those
  physical objects. In the current implementation, the GA finds operations that change the state of the units.
  Unit states are prescribed by the StateType enum class and unit operations are prescribed by the OpType
  enum class. The GA determines sequences of OpType operations that a simulation model uses to operate the units
  and have them change into different StateType states.
 
## An example
 
  Initially, a unit may be in the STOPPED state. The GA may prescribe a WARMUP_GEN operation for the next timestep.
  The simulation calls CalNextState with these params and the new state of the unit is now WARMUP_G. At this point
  ths unit is warming up and the unit simulation code for warming up discharges water from the pool but generates
  no power (ie the unit is not at-speed and is not yet frequency locked). The next timestep the GA prescribes
  CONTINUE. CalcNextState(WARMUP_G, CONTINUE) then results in GENERATE - and the unit is now generating power
  at the current head. But what is the power setting of the unit?
 
  First, some background on the implementation: Not only does the OpType enum class prescribe which state units
  should transition to but to also can provide a numerical value to that new state: OpTypes are BYTES with 2 bits
  reserved for the operation (4 operations) and 6 bits reserved for the numerical parameter (64 values).
 
  For any particular head there are lower and upper bounds on the operation of the unit (and over this range
  of operation the efficiency changes, see m_PHEArr). When a unit should be generating, the simulation code
  determines the min and range of feasible power (see how m_FeasZoneArr is used) and uses the OpType's numerical parameter
  to select operations in one 64 increments across that range (see how CalcSpan and getFrac are used).
 
## Summary of the workings
 
  With plants and units starting in an initial state, the GA provides changes to the unit states so that repeated
  simulations of the system will improve the performance, measured by some objective result of the simulation. The
  'objective result' if the value of the objective function and it is fed back to the GA to help it decide how
  to recombine previous simulation trials to create new and better simulation trials.
 
## Brief analysis of model & sniffle GA

Coming from a background where these kinds of problems are solved with LP/QP methods it can be pretty frustrating
trying to debug the code as the results are always different each time the darn thing runs. I suppose the prng could
be cooked to seed the same way each time, but that would be beside the point of having an effective GA in my opinion.

The real promise I had hoped to find with a GA method is to overcome many of the linear
limitations on constraints in LP/QP models (ie. piecewise constraints, operating points as a function of time,
 linear convolution for river routing).

In looking at how this particular problem seems to be solved by the sniffle GA I'm not sure it's the right fit;
the sniffle GA keeps distributions of the OpTypes for specific genes but it doesn't hold any auto/correlation
information that would allow it to synthesize patterns of operation. The emergence of patterns are left to
come out in the wash of randomness. It may be that better performance could be achieved by having a
heuristic first pass (it may not even need to be greatly feasible) and using sniffle to do all the tidying
 that makes it feasible and valuable.

  Thanks for reading!
