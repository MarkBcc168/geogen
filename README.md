# geogen

`geogen` is a fast C++20 geometry configuration explorer and lightweight
Olympiad-geometry prover. It builds a configuration from a small whitespace DSL,
finds collinear and concyclic point sets, runs a fixed-point angle chase, and
prints only coincidences which that chase does not already make routine.

This is an experimental problem-discovery tool, not a replacement for a formal
proof assistant. Coordinates are used to *discover* candidate statements;
symbolic facts are used to prove and filter them.

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
starts a comment. Objects must be declared before use. Initial coordinates should
be generic (avoid an isosceles or right triangle unless that is intentional),
because accidental coordinate symmetry creates accidental conjectures.

### Modes and initial configurations

```text
mode generate
mode prove
option show_easy 0             # print/suppress coincidences proved by the chase
option circle_budget 25000000  # skip general scan above this many triples

triangle A 0 0 B 6 0 C 1 4
quadrilateral A 0 0 B 5 0 C 6 3 D 1 4
cyclic_quad A 1 0 B 0 1 C -1 0 D 0 -1
point P 2.5 1.25
```

`cyclic_quad` validates its coordinates and immediately adds the full directed
angle relation for a cyclic quadrilateral.

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
- Directed line angles live modulo 180 degrees. Sparse Gaussian bases over two
  large primes track equations such as
  `angle(AB)+angle(CD)=angle(AD)+angle(BC)` without ever dividing a geometric
  relation by two. Two moduli make hash-field false proofs negligibly unlikely.
- The fixed-point chase adds cyclic converse facts and indexed kite consequences;
  orthocenter, circumcenter, incenter, midpoint, incidence, parallel,
  perpendicular, reflection, and angle-bisector constructions seed their standard
  relations. It stops when a pass adds no facts.
- A numerical coincidence which follows from this fact base is labeled `EASY`
  and suppressed by default. Remaining statements are printed as `NONTRIVIAL`.

The optional Pappus, radical-center, and Newton--Gauss rules are intentionally
not enabled in this first version: applying them blindly creates many auxiliary
objects and can dominate both runtime and output. They fit behind the same fact
closure interface if later benchmarks justify targeted versions.

## Example from the prompt

```text
mode generate
triangle A 0 0 B 6 0 C 1 4
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

Discovery uses `long double` with scale-aware validation after quantized hashing.
Do not use coordinates with wildly different magnitudes. A reported
`NONTRIVIAL` statement is a conjecture: verify it with several generic initial
coordinate choices, then prove it using a stronger prover or by hand.
