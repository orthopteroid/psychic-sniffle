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

Populations are double-buffered and partitioned into 7 groups, not necessarily of equal size. Each new population is built from the previous generation through the sequential assembly of the 7 groups:
 1. picking the maximally-best,
 2. weighted selection of highest value
 3. weighted selection of highest value with phenotype-byte jump mutations
 4. splicing maximally-best with a weighted high-value parent
 5. splicing a weighted high-value parent with maximally-best
 6. splicing pairs of weighted high-value parents
 7. assembling weighted high-value phenotype-bytes into random phenotypes

Categories 2 & 3 are considered 'elite' in that they are picked from groups 2-7 of the previous generation. Overall, the weighted selection process only includes groups 2-7, group 1 (the maximally-best) is never included in order to avoid saturation.

Current problems solved with this GA include:
* The Schwefel Function, over [-500,+500]. A 16-bit fixed-point 20 dimensional Schwefel, with a population of 400, can be solved in 168 iterations (so, 240000 evaluations).
