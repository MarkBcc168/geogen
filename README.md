# geogen

`geogen` is a fast C++20 geometry configuration explorer and lightweight
Olympiad-geometry prover. It builds a configuration from a small whitespace DSL,
finds collinear and concyclic point sets, runs a fixed-point angle chase, and
prints only coincidences which that chase does not already make routine.

This is an experimental problem-discovery tool, not a replacement for a formal
proof assistant. Initial figures are symbolic: the program samples several
generic realizations to *discover* candidate statements, while symbolic facts are
used to prove and filter them.

## Build and run

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/geogen examples/midpoint_perpendicular.geogen
ctest --test-dir build --output-on-failure
```

The executable reads standard input when no filename is supplied.

## Input language

One command appears on each line. Tokens are separated by whitespace and `#`
starts a comment. Objects must be declared before use. No coordinates are entered
for the initial triangle or quadrilateral.

### Modes and initial configurations

```text
mode generate
mode prove
option show_easy 0             # print/suppress coincidences proved by the chase
option circle_budget 25000000  # skip general scan above this many triples
option trials 5                # independent random realizations in generate mode
option seed 20260828            # reproducible pseudorandom seed
option max_points 30            # automatically expand to at most 30 points
option line_circle_intersections 1 # generate quadratic-free second intersections

triangle A B C
quadrilateral A B C D
cyclic_quad A B C D
```

For every trial, `triangle` and `quadrilateral` generate a fresh nondegenerate
generic realization. `cyclic_quad` generates four random points on a circle and
immediately adds the full directed-angle relation for a cyclic quadrilateral.
The seed is printed in the report so a run can be reproduced.

In generation mode, a coincidence is emitted only if the same named point set is
found in every trial. This makes accidental isosceles, right-angle, parallel, or
concyclic behavior overwhelmingly unlikely. It remains randomized evidence, not
a proof; `mode prove` uses the symbolic fact system for proof goals.

The advanced `point P x y` command remains available for diagnostics, but a
fixed-coordinate point changes the problem and should not normally be used in a
universal triangle configuration.

### Point generation and listing

Without `max_points`, the configuration contains exactly the initial and explicit
construction commands. Setting `option max_points N` enables automatic expansion.
The generator repeatedly tries circumcenters, orthocenters, and midpoints in a
stable symbolic order until it reaches `N` points (or exhausts applicable
constructions). The limit includes initial and explicitly constructed points and
may be set from 1 through 5000.

`option line_circle_intersections 1` additionally scans every non-segment line and
circle with a symbolically known common point. It constructs the other
intersection before the ordinary point expansion. This pass also works without
`max_points`; with a maximum, it respects the remaining point capacity.

Every report lists the resulting points without trial-specific coordinates:

```text
points=6
POINT A [random initial triangle]
POINT B [random initial triangle]
POINT C [random initial triangle]
POINT O(A,B,C) [circumcenter]
POINT H(A,B,C) [orthocenter]
POINT M(A,B) [midpoint]
```

The symbolic names are identical across random trials. In generation mode, only
points successfully constructed in every trial are retained in the final list.

Before insertion, every constructed point is compared with all existing points
using a scale-aware numerical tolerance. A coincident result is discarded: it is
not counted, printed, or used as a seed for automatic expansion. Its requested
name is retained only as an internal alias of the canonical point, so subsequent
explicit DSL commands remain valid. For example, constructing the orthocenter of
`A,B,H(A,B,C)` aliases that result to `C` instead of creating a recursive duplicate.

### Point, line, and circle constructions

```text
line AB A B
midpoint M A B
perp_bisector p A B
parallel q P AB
perpendicular r P AB
angle_bisector s A B C       # internal bisector of angle ABC
reflection_line X P AB
reflection_point X P O
foot H P AB
intersection_ll X p q

circumcenter O A B C
orthocenter H A B C
incenter I A B C
circle omega O A
circumcircle gamma A B C
incircle inc I A B C
```

The quadratic-free intersection commands require one already known root and
return the other root:

```text
intersection_lc_known X line_name circle_name K
intersection_cc_known Y circle_one circle_two K
```

`K` is checked to lie on both input objects. Tangency is rejected because it has
no distinct second point. Thus configuration construction never invokes a
general quadratic solver.

For example, the following automatically creates a point named
`X(secant,omega,A)` at the second intersection:

```text
mode generate
option line_circle_intersections 1
triangle A B C
incenter I A B C
line secant A I
circumcircle omega A B C
```

### Independent prover

Use `mode prove`, describe the configuration with the same commands, then add
one or more goals:

```text
prove_collinear A B C D
prove_concyclic A B C D E
prove_parallel A B C D
prove_perpendicular A B C D
prove_equal_distance A B C D
```

The last three mean `AB || CD`, `AB perpendicular CD`, and `AB = CD`. A
successful angular proof lists the construction/theorem facts used. See
[`examples/orthocenter.geogen`](examples/orthocenter.geogen).

## Search and filtering design

- Collinear sets are found by grouping normalized directions about each pivot.
  This takes `O(n^2 log n)` time and keeps only actual candidate groups.
- Declared circles are scanned in `O(n)`. General concyclicity uses a fixed-anchor
  circumcircle hash: `O(n^3)` time, `O(n^2)` peak memory, and no `O(n^4)` scan.
  `circle_budget` lets very large runs retain only declared-circle detection.
- Directed line angles live modulo 180 degrees. An exact sparse integer-lattice
  basis tracks equations such as
  `angle(AB)+angle(CD)=angle(AD)+angle(BC)` without ever dividing a geometric
  relation. The 90-degree constant has order two, so two perpendicular relations
  correctly add to 180 degrees, while `2*x=0` never incorrectly implies `x=0`.
- The fixed-point chase adds cyclic converse facts and indexed kite consequences;
  orthocenter, circumcenter, incenter, midpoint, incidence, parallel,
  perpendicular, reflection, and angle-bisector constructions seed their standard
  relations. It stops when a pass adds no facts.
- Midpoint closure includes the triangle midline theorem and the right-triangle
  hypotenuse-midpoint theorem. Perpendicular-bisector incidences add equal
  distances, kites add their full reflection-angle relation, and the standard six
  side-midpoint/altitude-foot configuration invokes the nine-point-circle rule.
- A numerical coincidence which follows from this fact base is labeled `EASY`
  and suppressed by default. Remaining statements are printed as `NONTRIVIAL`.

The optional Pappus, radical-center, and Newton--Gauss rules are intentionally
not enabled in this first version: applying them blindly creates many auxiliary
objects and can dominate both runtime and output. They fit behind the same fact
closure interface if later benchmarks justify targeted versions.

## Example from the prompt

```text
mode generate
option trials 7
triangle A B C
midpoint M B C
perp_bisector p A M
line AB A B
line AC A C
intersection_ll E p AC
intersection_ll F p AB
```

This is available as
[`examples/midpoint_perpendicular.geogen`](examples/midpoint_perpendicular.geogen).

## Numerical limits

Discovery uses internally sampled `long double` coordinates with scale-aware
validation after quantized hashing. A reported `NONTRIVIAL` statement has survived
all configured random trials, but is still a conjecture: prove it using the
independent prover, a stronger prover, or by hand.
