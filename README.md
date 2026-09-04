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
option line_circle_intersections 1 # default: enable known-root circle intersections
option angle_coefficient_limit 10000 # maximum accepted proof coefficient
option proof_scope ancestry     # prove each candidate from its definitions only

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

The symbolic construction and theorem closure are identical across generic
random realizations, so generation mode classifies coincidences during the first
trial only. Later trials perform the numerical construction and coincidence
scans, then intersect their findings with that classified set. This preserves
the output while avoiding repeated angle-lattice proofs.

By default, `option proof_scope global` lets the filter use every constructed
point, matching the traditional full-configuration prover. The conservative
`option proof_scope ancestry` rebuilds a separate subconfiguration for each
detected coincidence. It recursively retains only the points in the statement,
their defining point/line/circle constructions, and the initial configuration
commands needed by those definitions. Unrelated generated auxiliary points are
therefore unavailable to its proof. This mode is intentionally conservative and
may leave more coincidences unproved; it never promotes an unproved global
candidate to an easy one merely from numerical incidence.

The advanced `point P x y` command remains available for diagnostics, but a
fixed-coordinate point changes the problem and should not normally be used in a
universal triangle configuration.

### Point generation and listing

Without `max_points`, the configuration contains the initial and explicit
construction commands, plus the finite set of known-root line-circle
intersections already present in that configuration. Setting `option max_points N`
enables seeded random expansion until the configuration reaches `N` points (or
the available construction search is exhausted). The limit includes initial and
explicitly constructed points and may be set from 1 through 5000.

The random point-producing mix enables every implemented construction that does
not require finding unknown quadratic roots:

- midpoint, point reflection, line reflection, and perpendicular foot;
- circumcenter, orthocenter, and incenter;
- intersection of two nonparallel lines;
- the second line-circle or circle-circle intersection when one intersection is
  already known.

Supporting objects are generated randomly as well. Lines may be defined by two
points, as a perpendicular bisector or angle bisector, or as a parallel or
perpendicular through a point. Circles may be center-radius circles,
circumcircles, or incircles. Circumcircles through three existing points are the
dominant random circle choice. Known-root line-circle and circle-circle
intersections first reuse existing circles; when they need a supporting circle,
they create a circumcircle through the known intersection whenever possible.
These auxiliary objects are interleaved with point generation and bounded
indirectly by the point cap.

Input choices use weights proportional to `4^(-depth)`, normalized to the
shallowest currently available object. Each additional construction level is
therefore four times less likely to be selected. Initial and shallow points,
lines, and circles are strongly preferred over deeply nested ones.
The symbolic construction RNG uses `option seed` and is shared by all numerical
trials; only the initial coordinates vary. Consequently every trial tests the
same named construction graph. Set `option line_circle_intersections 0` to
exclude both known-root circle-intersection methods.

Every report lists the resulting points as parseable symbolic assignments,
without trial-specific coordinates. Supporting line and circle constructions are
expanded recursively, so each record is self-contained:

```text
points=7
POINT A = initial(A)
POINT B = initial(B)
POINT C = initial(C)
POINT M = midpoint(B,C)
POINT E = intersect(perpendicular_bisector(A,midpoint(B,C)),line(A,C))
POINT X = reflect(A,line(B,C))
POINT Y = other_intersection(line(A,C),circle(B,A),A)
```

The identifier before `=` is the stable name used by coincidence statements.
The expression after `=` records the complete definition. This keeps generated
IDs compact while showing exactly which earlier points and auxiliary objects
produce them. Points explicitly constructed in the input are treated as atomic
inside later generated definitions. Thus the prompt's `E` and `F` remain `E` and
`F`, rather than being repeatedly expanded into their line intersections; their
own `POINT E = ...` and `POINT F = ...` records still show their definitions.

The symbolic names are identical across random trials. In generation mode, only
points successfully constructed in every trial are retained in the final list.

Before insertion, every constructed point is compared with all existing points
using a scale-aware numerical tolerance. A coincident result is discarded: it is
not counted, printed, or used as a seed for automatic expansion. Its requested
name is retained only as an internal alias of the canonical point, so subsequent
explicit DSL commands remain valid. For example, constructing the orthocenter of
`A,B,H(A,B,C)` aliases that result to `C` instead of creating a recursive duplicate.

Lines are canonicalized in the same way. A construction that reproduces an
existing geometric line remains usable by its requested name, but it does not
create another line or another expansion seed. Point definitions are printed
with the canonical line expression. For example,
`perpendicular(C,perpendicular_bisector(B,C))` is rendered as `line(B,C)`.

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

The prover includes the general angle-defined kite congruence: if `AC`
bisects both endpoint angles of quadrilateral `ABCD`, then it derives
`AB = AD` and `CB = CD`. In particular, equal tangent lengths and equal
inradii for feet from an incenter are consequences of this kite rule; they are
not registered as a separate incenter-specific metric shortcut.

The quadratic-free intersection commands require one already known root and
return the other root:

```text
intersection_lc_known X line_name circle_name K
intersection_cc_known Y circle_one circle_two K
```

`K` is checked to lie on both input objects. Tangency is rejected because it has
no distinct second point. Thus configuration construction never invokes a
general quadratic solver.

When both circles have named centers, a known-root circle-circle intersection
also registers the standard common-chord theorem: the two intersection points
are reflections across the center line, and their chord is perpendicular to it.

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

- For configurations initialized by one free triangle, a small affine prover
  tracks barycentric coordinates through midpoints, point reflections, ordinary
  and parallel lines, and eligible line-line intersections. It registers
  collinear groups, parallel carriers, and midpoint identities before the angle
  closure runs, so those facts can trigger midline, parallelogram, cyclic, and
  other downstream rules. Identities must agree over two independent prime
  fields; metric constructions are intentionally left to the metric/angle layer.
- A complementary two-prime orthogonal-component prover works inside certified
  parallel/perpendicular direction classes. It combines line incidences,
  midpoints, feet, and reflections as scalar linear equations. In particular,
  it recognizes an isosceles trapezoid when the two leg vectors are certified
  mirror images across the direction perpendicular to its bases, then registers
  equal legs and cyclicity without relying on numerical equality.

- Collinear sets are found by grouping normalized directions about each pivot.
  This takes `O(n^2 log n)` time and keeps only actual candidate groups.
- Declared circles are scanned in `O(n)`. General concyclicity uses a fixed-anchor
  circumcircle hash: `O(n^3)` time, `O(n^2)` peak memory, and no `O(n^4)` scan.
  `circle_budget` lets very large runs retain only declared-circle detection.
- Directed line angles live modulo 180 degrees. For ordinary fact bases, an exact
  sparse integer-lattice basis tracks equations such as
  `angle(AB)+angle(CD)=angle(AD)+angle(BC)` without ever dividing a geometric
  relation. The 90-degree constant has order two, so two perpendicular relations
  correctly add to 180 degrees, while `2*x=0` never incorrectly implies `x=0`.
- Large generated configurations automatically switch to sparse elimination over
  two large prime fields before exact Hermite coefficients can swell. The prover
  tracks the linear combination of original theorem facts, converts every field
  coefficient to its balanced signed representative, and accepts the proof only
  when both primes reconstruct the same coefficients and every absolute value is
  at most `angle_coefficient_limit` (10,000 by default). Modular division such as
  `1/2` reconstructs near half the prime and is therefore rejected.
- Identical angle facts are stored once, and successful validated queries are
  cached. Prime-field membership is tested using coefficient-only sparse rows;
  the more expensive original-fact combination is reconstructed only after a
  candidate is known to lie in both field spans.
- Profiling the 100-point generator attributes roughly 30% of samples to the
  modular reducer as a whole, chiefly sparse-map traversal and certificate
  updates rather than isolated multiplication. Residues therefore remain in
  ordinary representation: converting the entire proof basis and certificates
  to Montgomery form did not have a sufficiently large measured target to
  justify the added representation and decoding complexity.
- The prover runs a global fixed point across affine incidence, metric closure,
  center and perpendicular-bisector loci, reflection propagation, directed-angle
  chasing, cyclic converses, indexed kite consequences, and the orthogonal
  component prover. A theorem discovered by a late phase is therefore available
  to every earlier phase on the next pass. Closure stops only when a complete
  pass adds no stored incidence, angle, cyclic, midpoint, or equal-length fact;
  there is no fixed round limit.
- Orthogonal projection preserves certified midpoints: if `M` is the midpoint
  of `AB`, their projections onto one line are `A1`, `M1`, and `B1`, then `M1`
  is registered as the midpoint of `A1B1`. A second formal-affine replay treats
  otherwise metric points as independent vectors, allowing midpoint and
  parallel-intersection identities to cancel those unknowns soundly.
- Orthocenters are recognized from their defining property, not only from an
  `orthocenter` command. On a relevant four-point set, any two certified
  altitudes register the third altitude and feed it back into angle closure.
- Midpoint closure includes the triangle midline theorem and the right-triangle
  hypotenuse-midpoint theorem. Perpendicular-bisector incidences add equal
  distances, kites add their full reflection-angle relation, and the standard six
  side-midpoint/altitude-foot configuration invokes the nine-point-circle rule.
- A point reflection `C = reflect(A,B)` registers `B` as the midpoint of `AC`.
  Whenever one point is the midpoint of both `AB` and `CD`, the parallelogram
  closure adds `AC || BD` and `AD || BC`. Conversely, the certified intersection
  of `AB` with its perpendicular bisector is registered as the midpoint of `AB`;
  this includes a circumcenter's perpendicular foot to a chord.
- Every known equality `XA=XB` adds the undivided isosceles relation
  `angle(XA)+angle(XB)=2*angle(AB)`. This retains integer coefficients and lets
  circumcenter-radius facts participate in the full angle chase. Conversely, a
  proved base-angle relation of this form registers `XA=XB`.
- The incenter of a certified isosceles triangle inherits the triangle's
  reflection symmetry: it is equidistant from the base endpoints, and its line
  to the apex is perpendicular to the base. This also lets nested centers reuse
  that symmetry axis.
- The perpendicular-bisector theorem is bidirectional: every point on the
  locus gains equal endpoint distances, and every known equal-distance point
  is attached to an existing canonical locus for those endpoints.
- A line reflection records its full undivided direction identity at every
  certified point of the mirror. In particular, reflection across an angle
  bisector maps one sideline to the other without requiring angle division.
  When the mirror is a perpendicular bisector, its defining endpoints are
  treated as a reflected pair, so distances and line directions are transported
  through the mirror without constructing additional image points.
  Reflection symmetry is also propagated to points attached to the mirror only
  during later theorem closure, such as a circumcenter on a perpendicular
  bisector.
- If a source point is reflected in a line and half-turned about a point on that
  line, the two images determine a line parallel to the mirror. This general
  composition rule is reused by cyclic and right-angle deductions.
- For a circle with a named center, every known line through the center and
  perpendicular to a certified chord is registered as the reflection axis of
  the chord endpoints. Subsequent reflection, kite, and isosceles-trapezoid
  rules then apply without constructing the chord midpoint.
- Definition-certified line carriers sharing two points, or sharing one point
  with a known parallel direction, are merged. Circles sharing three certified
  points are merged likewise. Thus every line-line, known-root line-circle, and
  known-root circle-circle intersection immediately inherits all corresponding
  collinear and cyclic incidences, even in a large maximal set.
- Feet from an incenter to adjacent sidelines obtain equal inradii and equal
  tangent lengths through the general angle-defined kite rule. Equal-length
  classes are closed transitively before all orientations of the isosceles and
  metric kite rules are applied.
- A constructed incenter contributes the three undivided internal-bisector
  cross-angle sums. They resolve the internal/external branch that doubled
  equations modulo 180 degrees deliberately cannot select. Consequently, the
  ordinary angle chase proves the right-angles-on-an-incircle-chord result: if
  `X = AI intersect DE` in the standard contact-triangle notation, it discovers
  `B,I,X,D,F` on the circle with diameter `BI` and derives
  `AX perpendicular BX` without constructing an auxiliary point.
- If `A,B,C,D` are cyclic and a point `O` has three known equal distances to
  members of the quadruple, the fourth radius is registered as equal. Named
  centers of constructed circles also contribute equal-radius facts for every
  point constructed on that circle.
- A constructed circumcenter also adds all cyclic forms of
  `angle(ACB)+angle(OAB)=90 degrees`. This resolves the factor-of-two ambiguity
  that equal-radius isosceles relations alone cannot remove modulo 180 degrees.
- Circumcenters through a fixed pair are registered on the pair's canonical
  perpendicular-bisector locus. Orthocenters are likewise registered on the
  altitude determined by a vertex and any known carrier of the opposite side.
  Consequently, collinear families such as `O(A,B,*)` and `H(A,P,Q)` for
  collinear `P,Q` are treated as construction-level facts. The chord midpoint
  and nested centers such as `H(A,B,O(A,B,C))` are attached to that same locus.
  A parity union-find matches known parallel and perpendicular carriers without
  repeatedly querying the full angle lattice.
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
