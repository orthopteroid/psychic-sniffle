# psychic-sniffle

A minimalistic but many featured genetic algorithm (maximizer) designed specifically to make use of phenotype-bytes' expected-values. The idea was to track and manipulate these distributions in order to control saturation and test distribution-related termination conditions (like signal/noise ratio).
 
 This GA seems to give OK performance. It is small and should be easy to follow and hack on. It's got:
* <300 lines of open-mp friendly code,
* flat (open-mp friendly) arrays for phenotypes,
* an open-mp friendly version of the Tausme88 PRNG,
* no solver-loop (you have that in your problem, along with your termination conditions),
* phenotype-byte distribution tracking in order to:
  1. avoid phenotype saturation
  2. perform 'jumping mutation' on phenotype-bytes
  3. construct random phenotypes using weighted high-value phenotype-bytes 

Populations are double-buffered and partitioned into 7 groups, not necessarily of equal size. Each new population is built from the previous generation through the sequential assembly of the 7 groups (the reason for all this silliness with groups is just to reduce indeterminate branching in the code so that I can keep the cache and the pipeline on my old machine happy):
 1. picking the maximally-best,
 2. weighted selection of highest value
 3. weighted selection of highest value with phenotype-byte jump mutations
 4. splicing maximally-best with a weighted high-value parent
 5. splicing a weighted high-value parent with maximally-best
 6. splicing pairs of weighted high-value parents
 7. assembling weighted high-value phenotype-bytes into random phenotypes

Categories 2 & 3 are considered 'elite' in that they are picked from groups 2-7 of the previous generation. Overall, the weighted selection process only includes groups 2-7, group 1 (the maximally-best) is never included in order to avoid saturation.

Current problems solved with this GA include:
* The Schwefel Function, over [-500,+500]. A 16-bit fixed-point 20 dimensional Schwefel, with
 a population of 400, can be solved in 168 iterations (so, 240000 evaluations). State-of-the-art
  (2016) is about 40000 evaluations (reference?).

* A 12 timestep, 2 Plant, 3 Unit hydropower nonlinear optimization problem. The plant reservoirs and tailwater
 curves are assumed to be linear but the unit performance curves are interpolated from a sampling-point
 cloud over the surface of a real unit performance curve (http://encyclopedia2.thefreedictionary.com/Hydroturbine
 (fig 6)). Some of the implementation tricks include: polygons to represent feasible and roughzone regions,
 64 discretization levels for power at each head, generate vs spin mode (for ancillary benifits: ie hz support).
  The objective function is an as-hoc mishmash of things polluted with weighting coefficients (ah well);
   1. efficiency is maximized,
   2. deviation from demand is minimzed,
   3. start-stops are minimized,
   4. roughzone operations are minimzed.
