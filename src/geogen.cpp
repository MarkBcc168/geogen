#include <algorithm>
#include <array>
#include <cmath>
#include <compare>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace geogen {

struct AffineScalar {
  static constexpr std::int64_t p1=1000000007,p2=1000000009;
  std::int64_t a=0,b=0;
  AffineScalar()=default;
  AffineScalar(std::int64_t x):a((x%p1+p1)%p1),b((x%p2+p2)%p2){}
  static std::int64_t power(std::int64_t x,std::int64_t n,std::int64_t p){std::int64_t r=1;while(n){if(n&1)r=r*x%p;x=x*x%p;n>>=1;}return r;}
  AffineScalar inverse()const{if(a==0||b==0)throw std::runtime_error("singular affine field element");AffineScalar r;r.a=power(a,p1-2,p1);r.b=power(b,p2-2,p2);return r;}
  friend AffineScalar operator+(AffineScalar x,AffineScalar y){AffineScalar r;r.a=(x.a+y.a)%p1;r.b=(x.b+y.b)%p2;return r;}
  friend AffineScalar operator-(AffineScalar x,AffineScalar y){AffineScalar r;r.a=(x.a-y.a+p1)%p1;r.b=(x.b-y.b+p2)%p2;return r;}
  friend AffineScalar operator*(AffineScalar x,AffineScalar y){AffineScalar r;r.a=x.a*y.a%p1;r.b=x.b*y.b%p2;return r;}
  friend AffineScalar operator/(AffineScalar x,AffineScalar y){return x*y.inverse();}
  auto operator<=>(const AffineScalar&)const=default;
};
using AffineVector = std::array<AffineScalar,3>;

AffineVector affine_cross(const AffineVector&a,const AffineVector&b){
  return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}
bool affine_zero(const AffineVector&a){return a[0]==AffineScalar(0)&&a[1]==AffineScalar(0)&&a[2]==AffineScalar(0);}
AffineVector affine_normalize_projective(AffineVector a){
  for(const auto&x:a)if(x.a!=0&&x.b!=0){AffineScalar pivot=x;for(auto&y:a)y=y/pivot;return a;}
  return {AffineScalar(0),AffineScalar(0),AffineScalar(0)};
}
std::optional<AffineVector> affine_normalize_point(AffineVector a){
  AffineScalar sum=a[0]+a[1]+a[2];if(sum.a==0||sum.b==0)return std::nullopt;
  for(auto&x:a)x=x/sum;
  return a;
}

// Homogeneous scalar linear equations checked simultaneously over two large
// prime fields. Unlike the directed-angle lattice, ordinary coordinate
// equations may be divided during Gaussian elimination.
class PairedLinearSystem {
  using Row=std::map<int,AffineScalar>;
  std::map<int,Row> basis_;
  static Row make_row(std::initializer_list<std::pair<int,int>> terms){
    Row row;for(auto [variable,coefficient]:terms){auto value=row[variable]+AffineScalar(coefficient);
      if(value==AffineScalar(0))row.erase(variable);else row[variable]=value;}return row;
  }
  void reduce(Row&row)const{while(!row.empty()){auto found=basis_.find(row.begin()->first);if(found==basis_.end())break;
      AffineScalar factor=AffineScalar(0)-row.begin()->second;for(const auto&[variable,coefficient]:found->second){
        auto value=row[variable]+factor*coefficient;if(value==AffineScalar(0))row.erase(variable);else row[variable]=value;}}}
 public:
  void add(std::initializer_list<std::pair<int,int>> terms){Row row=make_row(terms);reduce(row);if(row.empty())return;
    AffineScalar inverse=row.begin()->second.inverse();for(auto&[_,coefficient]:row)coefficient=coefficient*inverse;
    basis_[row.begin()->first]=std::move(row);
  }
  bool proves(std::initializer_list<std::pair<int,int>> terms)const{Row row=make_row(terms);reduce(row);return row.empty();}
};

constexpr long double EPS = 1e-9L;
constexpr long double PI = 3.141592653589793238462643383279502884L;

struct Point {
  std::string name;
  long double x{}, y{};
  std::string origin;
};

struct Line {
  std::string name;
  long double a{}, b{}, c{}; // ax + by + c = 0, a^2+b^2=1
  std::string origin;
};

struct Circle {
  std::string name;
  Point center;
  long double r2{};
  std::string origin;
};

long double sq(long double x) { return x * x; }
long double dist2(const Point& a, const Point& b) {
  return sq(a.x - b.x) + sq(a.y - b.y);
}
long double cross(const Point& a, const Point& b, const Point& c) {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}
long double dot(long double ax, long double ay, long double bx, long double by) {
  return ax * bx + ay * by;
}
long double scale(const Point& a, const Point& b, const Point& c) {
  return 1.0L + std::sqrt(dist2(a, b)) + std::sqrt(dist2(a, c));
}
bool near(long double a, long double b, long double s = 1.0L) {
  return std::fabs(a - b) <= EPS * (s + std::fabs(a) + std::fabs(b));
}

Line through(std::string name, const Point& p, const Point& q, std::string why) {
  long double dx = q.x - p.x, dy = q.y - p.y;
  long double z = std::hypotl(dx, dy);
  if (z <= EPS) throw std::runtime_error("cannot define a line through coincident points");
  long double a = -dy / z, b = dx / z, c = -(a * p.x + b * p.y);
  if (a < -EPS || (std::fabs(a) <= EPS && b < 0)) { a = -a; b = -b; c = -c; }
  return {std::move(name), a, b, c, std::move(why)};
}

Point intersect(const Line& l, const Line& m, std::string name, std::string why) {
  long double d = l.a * m.b - m.a * l.b;
  if (std::fabs(d) <= EPS) throw std::runtime_error("parallel lines have no finite intersection");
  return {std::move(name), (l.b * m.c - m.b * l.c) / d,
          (l.c * m.a - m.c * l.a) / d, std::move(why)};
}

Point circumcenter(const Point& a, const Point& b, const Point& c, std::string name,
                   std::string why) {
  long double d = 2 * cross(a, b, c);
  if (std::fabs(d) <= EPS * scale(a, b, c))
    throw std::runtime_error("circumcenter needs three non-collinear points");
  long double aa = a.x * a.x + a.y * a.y;
  long double bb = b.x * b.x + b.y * b.y;
  long double cc = c.x * c.x + c.y * c.y;
  return {std::move(name),
          (aa * (b.y-c.y) + bb * (c.y-a.y) + cc * (a.y-b.y)) / d,
          (aa * (c.x-b.x) + bb * (a.x-c.x) + cc * (b.x-a.x)) / d,
          std::move(why)};
}

// Exact integer-lattice elimination for directed angles in R/(pi Z). Relations
// may be added and subtracted, but never divided: 2*x=0 must not imply x=0.
class AngleSystem {
 public:
  using Integer=__int128_t;
  using Coeff=std::vector<std::pair<int,Integer>>;
  struct Equation {Coeff c;std::set<int> reasons;};

 private:
  using ModularCoeff=std::vector<std::pair<int,std::int64_t>>;
  static constexpr int HALF_TURN_COLUMN=1000000000;
  std::vector<Equation> rows_;
  std::set<Coeff> row_keys_;
  mutable std::vector<Equation> basis_;
  mutable std::map<Coeff,std::set<int>> proven_cache_;
  mutable bool dirty_=true;
  mutable bool modular_mode_=false;
  struct ModularEquation {
    ModularCoeff c;
    ModularCoeff combination; // coefficients of original fact rows
  };
  static constexpr std::array<std::int64_t,2> MODULI{1000000007LL,1000000009LL};
  mutable std::vector<std::unordered_map<int,ModularEquation>> modular_basis_;
  std::int64_t coefficient_limit_=10000;
  std::vector<std::string> reason_text_;

  template<class Value>
  static Value coefficient(const std::vector<std::pair<int,Value>>&row,int column){
    auto it=std::lower_bound(row.begin(),row.end(),column,
      [](const auto&entry,int key){return entry.first<key;});
    return it==row.end()||it->first!=column?Value(0):it->second;
  }
  static Integer checked_multiply(Integer a,Integer b){
    Integer result;
    if(__builtin_mul_overflow(a,b,&result))
      throw std::overflow_error("angle lattice coefficient overflow");
    return result;
  }
  static Integer checked_add(Integer a,Integer b){
    Integer result;
    if(__builtin_add_overflow(a,b,&result))
      throw std::overflow_error("angle lattice coefficient overflow");
    return result;
  }
  static void add_scaled(Equation&x,const Equation&y,const Integer&f){
    if(f==0)return;
    Coeff merged;merged.reserve(x.c.size()+y.c.size());std::size_t i=0,j=0;
    while(i<x.c.size()||j<y.c.size()){
      if(j==y.c.size()||(i<x.c.size()&&x.c[i].first<y.c[j].first))merged.push_back(x.c[i++]);
      else if(i==x.c.size()||y.c[j].first<x.c[i].first){Integer value=checked_multiply(f,y.c[j].second);if(value)merged.push_back({y.c[j].first,value});++j;}
      else {Integer value=checked_add(x.c[i].second,checked_multiply(f,y.c[j].second));if(value)merged.push_back({x.c[i].first,value});++i;++j;}
    }
    x.c=std::move(merged);
    x.reasons.insert(y.reasons.begin(), y.reasons.end());
  }
  static std::int64_t mod_norm(std::int64_t x,std::int64_t p){x%=p;if(x<0)x+=p;return x;}
  static std::int64_t mod_mul(std::int64_t a,std::int64_t b,std::int64_t p){return (a*b)%p;}
  static std::int64_t mod_power(std::int64_t a,std::int64_t e,std::int64_t p){std::int64_t r=1;while(e){if(e&1)r=mod_mul(r,a,p);a=mod_mul(a,a,p);e>>=1;}return r;}
  static ModularEquation modularized(const Equation&e,std::int64_t p){
    ModularEquation out;for(const auto&[v,a]:e.c){
      std::int64_t value=mod_norm(static_cast<std::int64_t>(a%p),p);if(value)out.c.push_back({v,value});}
    return out;
  }
  static void modular_add_scaled(ModularCoeff&x,const ModularCoeff&y,std::int64_t f,std::int64_t p){
    if(!f)return;
    ModularCoeff merged;merged.reserve(x.size()+y.size());std::size_t i=0,j=0;
    while(i<x.size()||j<y.size()){
      if(j==y.size()||(i<x.size()&&x[i].first<y[j].first))merged.push_back(x[i++]);
      else if(i==x.size()||y[j].first<x[i].first){auto value=mod_mul(f,y[j].second,p);if(value)merged.push_back({y[j].first,value});++j;}
      else {auto value=x[i].second+mod_mul(f,y[j].second,p);if(value>=p)value-=p;if(value)merged.push_back({x[i].first,value});++i;++j;}
    }
    x=std::move(merged);
  }
  static void modular_add_scaled(ModularEquation&x,const ModularEquation&y,std::int64_t f,std::int64_t p){
    modular_add_scaled(x.c,y.c,f,p);
    modular_add_scaled(x.combination,y.combination,f,p);
  }
  static void modular_reduce(ModularEquation&e,const std::unordered_map<int,ModularEquation>&basis,std::int64_t p){
    while(!e.c.empty()){int pivot=e.c.front().first;auto it=basis.find(pivot);if(it==basis.end())break;modular_add_scaled(e,it->second,p-e.c.front().second,p);}
  }
  static bool modular_member(const Equation&e,const std::unordered_map<int,ModularEquation>&basis,std::int64_t p){
    auto row=modularized(e,p);
    while(!row.c.empty()){
      int pivot=row.c.front().first;auto it=basis.find(pivot);if(it==basis.end())return false;
      modular_add_scaled(row.c,it->second.c,p-row.c.front().second,p);
    }
    return true;
  }
  void modular_insert(const Equation&e,std::size_t mi,int fact_id)const{
    auto p=MODULI[mi];auto row=modularized(e,p);row.combination.push_back({fact_id,1});
    modular_reduce(row,modular_basis_[mi],p);if(row.c.empty())return;
    int pivot=row.c.front().first;auto inv=mod_power(row.c.front().second,p-2,p);
    for(auto&[_,a]:row.c)a=mod_mul(a,inv,p);
    for(auto&[_,a]:row.combination)a=mod_mul(a,inv,p);
    modular_basis_[mi][pivot]=std::move(row);
  }
  void rebuild_modular()const{
    modular_basis_.assign(MODULI.size(),{});Equation period;period.c.push_back({HALF_TURN_COLUMN,2});
    for(std::size_t mi=0;mi<MODULI.size();++mi){modular_insert(period,mi,-1);for(std::size_t i=0;i<rows_.size();++i)modular_insert(rows_[i],mi,static_cast<int>(i));}
    basis_.clear();modular_mode_=true;dirty_=false;
  }
  static Equation converted(const Coeff& c,std::int64_t numerator,
                            std::int64_t denominator){
    if(denominator==0||(2*numerator)%denominator!=0)
      throw std::runtime_error("angle constant is not an integral multiple of pi/2");
    Equation e;e.c=c;
    // h=pi/2 is deliberately ordered after all line variables. Eliminating this
    // dense torsion column first causes catastrophic coefficient swell.
    Integer half_turns=(2*numerator)/denominator;
    if(half_turns!=0)e.c.push_back({HALF_TURN_COLUMN,-half_turns});
    return e;
  }

  void rebuild()const{
    if(!dirty_)return;
    // Exact HNF is preferred for small theorem bases. Large generated figures
    // use a multi-modulus lattice certificate to avoid unbounded coefficient
    // swell; modulus 2 specifically preserves the no-angle-halving invariant.
    if(rows_.size()>120){rebuild_modular();return;}
    modular_mode_=false;
    std::vector<Equation> work=rows_;
    Equation period;period.c.push_back({HALF_TURN_COLUMN,2});work.push_back(std::move(period)); // pi=0
    std::size_t rank=0;std::set<int> occupied;
    for(const auto&row:work)for(const auto&[column,_]:row.c)occupied.insert(column);
    for(int col:occupied){if(rank>=work.size())break;
      std::size_t chosen=rank;while(chosen<work.size()&&coefficient(work[chosen].c,col)==0)++chosen;
      if(chosen==work.size())continue;
      std::swap(work[rank],work[chosen]);
      for(std::size_t i=rank+1;i<work.size();++i){
        while(coefficient(work[i].c,col)!=0){
          Integer a=coefficient(work[rank].c,col),b=coefficient(work[i].c,col);
          add_scaled(work[rank],work[i],-(a/b));
          std::swap(work[rank],work[i]);
        }
      }
      if(coefficient(work[rank].c,col)<0)
        for(auto& [_,value]:work[rank].c)value=-value;
      ++rank;
    }
    basis_.assign(work.begin(),work.begin()+static_cast<std::ptrdiff_t>(rank));dirty_=false;
  }

 public:
  std::size_t fact_count()const{return rows_.size();}
  void set_coefficient_limit(std::int64_t limit){
    if(limit<1)throw std::runtime_error("angle coefficient limit must be positive");
    coefficient_limit_=limit;proven_cache_.clear();
  }
  int add_reason(std::string s) {
    reason_text_.push_back(std::move(s)); return static_cast<int>(reason_text_.size()) - 1;
  }

  bool add(const Coeff& c, std::int64_t num, std::int64_t den, const std::string& why) {
    Equation e=converted(c,num,den);
    if(!row_keys_.insert(e.c).second)return false;
    e.reasons.insert(add_reason(why));
    rows_.push_back(e);
    if(modular_mode_&&!dirty_)for(std::size_t mi=0;mi<MODULI.size();++mi)modular_insert(e,mi,static_cast<int>(rows_.size()-1));
    else dirty_=true;
    return true;
  }

  bool proves(const Coeff& c, std::int64_t num, std::int64_t den,
              std::set<int>* reasons = nullptr) const {
    Equation e=converted(c,num,den);
    Coeff query=e.c;
    if(auto it=proven_cache_.find(query);it!=proven_cache_.end()){
      if(reasons)*reasons=it->second;
      return true;
    }
    rebuild();
    if(modular_mode_){
      // Most theorem probes fail. Test field membership without carrying the
      // much larger proof-combination maps, then reconstruct and validate a
      // small integer certificate only for successful candidates.
      for(std::size_t mi=0;mi<MODULI.size();++mi)
        if(!modular_member(e,modular_basis_[mi],MODULI[mi]))return false;
      ModularCoeff expected;for(std::size_t mi=0;mi<MODULI.size();++mi){
      auto p=MODULI[mi];auto row=modularized(e,p);modular_reduce(row,modular_basis_[mi],p);
      ModularCoeff signed_coefficients;
      for(auto[id,value]:row.combination){if(value>p/2)value-=p;if(std::llabs(value)>coefficient_limit_)return false;if(value)signed_coefficients.push_back({id,value});}
      if(mi==0)expected=std::move(signed_coefficients);else if(signed_coefficients!=expected)return false;
    }
      std::set<int> proof_reasons;
      for(const auto&[fact,coefficient]:expected)if(fact>=0&&coefficient)
        proof_reasons.insert(rows_[static_cast<std::size_t>(fact)].reasons.begin(),rows_[static_cast<std::size_t>(fact)].reasons.end());
      proven_cache_[query]=proof_reasons;if(reasons)*reasons=std::move(proof_reasons);
      return true;
    }
    for(const auto& row:basis_){
      if(row.c.empty())continue;
      int pivot=row.c.front().first;
      Integer value=coefficient(e.c,pivot);if(value==0)continue;
      Integer divisor=row.c.front().second;
      if(value%divisor!=0)return false;
      add_scaled(e,row,-(value/divisor));
    }
    if(!e.c.empty())return false;
    proven_cache_[query]=e.reasons;if(reasons)*reasons=std::move(e.reasons);
    return true;
  }

  std::vector<std::string> explain(const std::set<int>& ids) const {
    std::vector<std::string> out;
    for (int id : ids) if (id >= 0 && static_cast<std::size_t>(id) < reason_text_.size())
      out.push_back(reason_text_[static_cast<std::size_t>(id)]);
    return out;
  }
};

struct Goal { std::string kind; std::vector<std::string> args; };
struct Candidate { std::string kind; std::vector<int> points; std::string source; };
struct MidpointFact {int midpoint,a,b;};
struct PerpendicularBisectorFact {int line,a,b;};
struct FootFact {int foot,source,line;};
struct CircumcenterFact {int center,a,b,c;};
struct OrthocenterFact {int center,a,b,c;};
struct IncenterFact {int center,a,b,c;};
struct AngleBisectorLineFact {int line,vertex;};
struct LineReflectionFact {int image,source,line;};
enum class ObjectKind {point,line,circle};
struct ObjectRef {
  ObjectKind kind;std::string name;
  auto operator<=>(const ObjectRef&)const=default;
};

class Engine {
  std::vector<Point> points_;
  std::vector<Line> lines_;
  std::vector<Circle> circles_;
  std::vector<int> point_depth_,line_depth_,circle_depth_;
  std::vector<std::vector<int>> line_points_;
  std::vector<std::vector<int>> circle_points_;
  std::vector<int> circle_center_ids_;
  std::vector<int> direction_parent_,direction_rank_,direction_parity_;
  std::unordered_map<std::string, int> point_id_, line_id_, circle_id_;
  std::unordered_map<std::string,std::string> line_canonical_name_;
  std::map<std::pair<int,int>, int> segment_line_;
  AngleSystem angles_;
  std::set<std::pair<std::pair<int,int>, std::pair<int,int>>> equal_lengths_;
  std::vector<std::array<int,4>> cyclic_facts_;
  std::vector<Goal> goals_;
  std::vector<Candidate> circle_cache_;
  std::vector<MidpointFact> midpoint_facts_;
  std::vector<PerpendicularBisectorFact> perpendicular_bisectors_;
  std::vector<FootFact> foot_facts_;
  std::vector<CircumcenterFact> circumcenter_facts_;
  std::vector<OrthocenterFact> orthocenter_facts_;
  std::vector<IncenterFact> incenter_facts_;
  std::vector<AngleBisectorLineFact> angle_bisector_line_facts_;
  std::vector<std::pair<int,int>> angle_bisector_diagonals_;
  std::vector<LineReflectionFact> line_reflection_facts_;
  std::map<std::pair<int,int>,int> perpendicular_bisector_loci_;
  std::vector<std::vector<std::string>> construction_commands_;
  std::set<std::string> input_point_names_;
  std::set<std::string> initial_point_names_;
  bool record_commands_=true,ancestry_scope_=false;
  std::set<std::string> disabled_constructions_;
  bool symmetry_enabled_=false,symmetry_ready_=false,symmetric_only_=false;
  std::string symmetry_first_,symmetry_second_;
  std::map<ObjectRef,ObjectRef> symmetry_dual_;
  std::set<int> symmetry_primary_points_,symmetry_primary_lines_,symmetry_primary_circles_;
  std::int64_t angle_coefficient_limit_=10000;
  std::vector<std::array<int,3>> initial_triangles_;
  bool prove_mode_ = false, show_easy_ = false;
  std::size_t circle_budget_ = 25000000;
  std::size_t max_points_=0;
  std::size_t automatic_serial_=0;
  std::string automatic_namespace_;
  std::uint64_t seed_,generation_seed_;
  std::mt19937_64 rng_,generation_rng_;

  static bool is_construction(const std::string&op){
    static const std::set<std::string> operations{
      "triangle","quadrilateral","cyclic_quad","point","line","midpoint",
      "perp_bisector","parallel","perpendicular","angle_bisector",
      "reflection_line","reflection_point","foot","intersection_ll",
      "circumcenter","orthocenter","incenter","circle","circumcircle",
      "incircle","intersection_lc_known","intersection_cc_known"};
    return operations.count(op)!=0;
  }
  bool construction_enabled(const std::string&op)const{return !disabled_constructions_.count(op);}
  static const std::set<std::string>& automatic_construction_types(){
    static const std::set<std::string> operations{
      "midpoint","reflection_point","reflection_line","foot","circumcenter",
      "orthocenter","incenter","intersection_ll","intersection_lc_known",
      "intersection_cc_known","line","perp_bisector","parallel",
      "perpendicular","angle_bisector","circle","circumcircle","incircle"};
    return operations;
  }
  static std::string canonical_construction_option(std::string op){
    static const std::map<std::string,std::string> aliases{
      {"point_reflection","reflection_point"},{"line_reflection","reflection_line"},
      {"perpendicular_bisector","perp_bisector"},{"line_line_intersection","intersection_ll"},
      {"line_circle_intersection","intersection_lc_known"},
      {"circle_circle_intersection","intersection_cc_known"}};
    if(auto found=aliases.find(op);found!=aliases.end())return found->second;
    return op;
  }
  static std::vector<ObjectRef> command_outputs(const std::vector<std::string>&t){
    const auto&op=t[0];std::vector<ObjectRef> out;
    if(op=="triangle")for(int i=1;i<=3;++i)out.push_back({ObjectKind::point,t[static_cast<std::size_t>(i)]});
    else if(op=="quadrilateral"||op=="cyclic_quad")for(int i=1;i<=4;++i)out.push_back({ObjectKind::point,t[static_cast<std::size_t>(i)]});
    else if(op=="line"||op=="perp_bisector"||op=="parallel"||op=="perpendicular"||op=="angle_bisector")out.push_back({ObjectKind::line,t[1]});
    else if(op=="circle"||op=="circumcircle"||op=="incircle")out.push_back({ObjectKind::circle,t[1]});
    else out.push_back({ObjectKind::point,t[1]});
    return out;
  }
  static std::vector<ObjectRef> command_inputs(const std::vector<std::string>&t){
    const auto&op=t[0];std::vector<ObjectRef> in;
    auto point=[&](std::size_t i){in.push_back({ObjectKind::point,t[i]});};
    auto line=[&](std::size_t i){in.push_back({ObjectKind::line,t[i]});};
    auto circle=[&](std::size_t i){in.push_back({ObjectKind::circle,t[i]});};
    if(op=="line"||op=="midpoint"||op=="perp_bisector"||op=="reflection_point"){point(2);point(3);}
    else if(op=="parallel"||op=="perpendicular"||op=="reflection_line"||op=="foot"){point(2);line(3);}
    else if(op=="angle_bisector"||op=="circumcenter"||op=="orthocenter"||op=="incenter"||op=="circumcircle"){point(2);point(3);point(4);}
    else if(op=="intersection_ll"){line(2);line(3);}
    else if(op=="circle"){point(2);point(3);}
    else if(op=="incircle"){point(2);point(3);point(4);point(5);}
    else if(op=="intersection_lc_known"){line(2);circle(3);point(4);}
    else if(op=="intersection_cc_known"){circle(2);circle(3);point(4);}
    return in;
  }

  ObjectRef canonical_object(ObjectRef object)const{
    if(object.kind==ObjectKind::point)object.name=points_[static_cast<std::size_t>(pid(object.name))].name;
    else if(object.kind==ObjectKind::line)object.name=lines_[static_cast<std::size_t>(lid(object.name))].name;
    return object;
  }
  void mark_primary_output(const ObjectRef&output){
    auto object=canonical_object(output);
    if(object.kind==ObjectKind::point)symmetry_primary_points_.insert(pid(object.name));
    else if(object.kind==ObjectKind::line)symmetry_primary_lines_.insert(lid(object.name));
    else symmetry_primary_circles_.insert(cid(object.name));
  }
  void mark_primary_commands(std::size_t begin,std::size_t end){
    if(!symmetry_enabled_)return;
    for(std::size_t i=begin;i<end;++i)
      for(const auto&output:command_outputs(construction_commands_[i]))mark_primary_output(output);
  }
  std::string dual_name(const ObjectRef&object)const{
    ObjectRef canonical=canonical_object(object);auto found=symmetry_dual_.find(canonical);
    if(found==symmetry_dual_.end())throw std::runtime_error("symmetry has no dual for "+object.name);
    return found->second.name;
  }
  std::vector<std::string> make_dual_command(const std::vector<std::string>&source){
    const auto&op=source[0];auto dual=source;
    if(op=="triangle"||op=="quadrilateral"||op=="cyclic_quad"||op=="point")
      throw std::runtime_error("cannot replay an initial declaration as a symmetric construction");
    dual[1]="S$"+source[1];
    auto point=[&](std::size_t i){dual[i]=dual_name({ObjectKind::point,source[i]});};
    auto line=[&](std::size_t i){dual[i]=dual_name({ObjectKind::line,source[i]});};
    auto circle=[&](std::size_t i){dual[i]=dual_name({ObjectKind::circle,source[i]});};
    if(op=="line"||op=="midpoint"||op=="perp_bisector"||op=="reflection_point"){point(2);point(3);}
    else if(op=="parallel"||op=="perpendicular"||op=="reflection_line"||op=="foot"){point(2);line(3);}
    else if(op=="angle_bisector"||op=="circumcenter"||op=="orthocenter"||op=="incenter"||op=="circumcircle"){point(2);point(3);point(4);}
    else if(op=="intersection_ll"){line(2);line(3);}
    else if(op=="circle"){point(2);point(3);}
    else if(op=="incircle"){point(2);point(3);point(4);point(5);}
    else if(op=="intersection_lc_known"){line(2);circle(3);point(4);}
    else if(op=="intersection_cc_known"){circle(2);circle(3);point(4);}
    else throw std::runtime_error("symmetry does not support construction "+op);
    return dual;
  }
  void register_dual_outputs(const std::vector<std::string>&source,const std::vector<std::string>&dual){
    auto source_outputs=command_outputs(source),dual_outputs=command_outputs(dual);
    for(std::size_t i=0;i<source_outputs.size();++i){
      auto a=canonical_object(source_outputs[i]),b=canonical_object(dual_outputs[i]);
      symmetry_dual_[a]=b;symmetry_dual_[b]=a;
      if(a.kind==ObjectKind::point){int ai=pid(a.name),bi=pid(b.name);int depth=point_depth_[static_cast<std::size_t>(ai)];
        if(points_[static_cast<std::size_t>(bi)].name==dual_outputs[i].name)point_depth_[static_cast<std::size_t>(bi)]=depth;
        else point_depth_[static_cast<std::size_t>(bi)]=std::min(point_depth_[static_cast<std::size_t>(bi)],depth);
      } else if(a.kind==ObjectKind::line){int ai=lid(a.name),bi=lid(b.name);int depth=line_depth_[static_cast<std::size_t>(ai)];
        if(lines_[static_cast<std::size_t>(bi)].name==dual_outputs[i].name)line_depth_[static_cast<std::size_t>(bi)]=depth;
        else line_depth_[static_cast<std::size_t>(bi)]=std::min(line_depth_[static_cast<std::size_t>(bi)],depth);
      } else circle_depth_[static_cast<std::size_t>(cid(b.name))]=circle_depth_[static_cast<std::size_t>(cid(a.name))];
    }
  }
  void replay_symmetric_commands(std::size_t begin,std::size_t end){
    for(std::size_t i=begin;i<end;++i){
      auto source=construction_commands_[i];
      if(source[0]=="triangle"||source[0]=="quadrilateral"||source[0]=="cyclic_quad")continue;
      if(source[0]=="point")throw std::runtime_error("option symmetry does not support coordinate point declarations");
      auto dual=make_dual_command(source);execute(dual,0);register_dual_outputs(source,dual);
    }
  }
  void establish_symmetry(){
    if(!symmetry_enabled_||symmetry_ready_)return;
    if(!initial_point_names_.count(symmetry_first_)||!initial_point_names_.count(symmetry_second_)||symmetry_first_==symmetry_second_)
      throw std::runtime_error("option symmetry requires two distinct initial point names");
    for(const auto&name:initial_point_names_)symmetry_dual_[{ObjectKind::point,name}]={ObjectKind::point,name};
    symmetry_dual_[{ObjectKind::point,symmetry_first_}]={ObjectKind::point,symmetry_second_};
    symmetry_dual_[{ObjectKind::point,symmetry_second_}]={ObjectKind::point,symmetry_first_};
    std::size_t explicit_end=construction_commands_.size();
    mark_primary_commands(0,explicit_end);
    replay_symmetric_commands(0,explicit_end);symmetry_ready_=true;
  }

  long double random_real(long double lo, long double hi) {
    return std::uniform_real_distribution<long double>(lo, hi)(rng_);
  }

  int pid(const std::string& s) const {
    auto it = point_id_.find(s); if (it == point_id_.end()) throw std::runtime_error("unknown point: " + s);
    return it->second;
  }
  int lid(const std::string& s) const {
    auto it = line_id_.find(s); if (it == line_id_.end()) throw std::runtime_error("unknown line: " + s);
    return it->second;
  }
  int cid(const std::string& s) const {
    auto it = circle_id_.find(s); if (it == circle_id_.end()) throw std::runtime_error("unknown circle: " + s);
    return it->second;
  }
  bool add_point(Point p) {
    if (point_id_.count(p.name)) throw std::runtime_error("duplicate point: " + p.name);
    for(std::size_t i=0;i<points_.size();++i){const auto&q=points_[i];
      long double magnitude=1+std::hypotl(p.x,p.y)+std::hypotl(q.x,q.y);
      if(std::sqrt(dist2(p,q))<=10*EPS*magnitude){
        // Preserve the discarded name as an alias for later DSL references, but
        // never expose it as a point or use it as an automatic construction seed.
        point_id_[p.name]=static_cast<int>(i);return false;
      }
    }
    if(max_points_&&points_.size()>=max_points_)
      throw std::runtime_error("configuration exceeds option max_points");
    point_id_[p.name] = static_cast<int>(points_.size()); points_.push_back(std::move(p));
    point_depth_.push_back(0);
    return true;
  }
  int add_line(Line l) {
    if (line_id_.count(l.name)) throw std::runtime_error("duplicate line: " + l.name);
    long double norm=std::hypotl(l.a,l.b);if(norm<=EPS)throw std::runtime_error("cannot add a degenerate line");
    l.a/=norm;l.b/=norm;l.c/=norm;
    if(l.a < -EPS || (std::fabs(l.a)<=EPS&&l.b<0)){l.a=-l.a;l.b=-l.b;l.c=-l.c;}
    for(std::size_t i=0;i<lines_.size();++i){const auto&q=lines_[i];
      long double tolerance=10*EPS*(1+std::fabs(l.c)+std::fabs(q.c));
      if(std::fabs(l.a-q.a)<=10*EPS&&std::fabs(l.b-q.b)<=10*EPS&&std::fabs(l.c-q.c)<=tolerance){
        line_id_[l.name]=static_cast<int>(i);line_canonical_name_[l.name]=q.name;return static_cast<int>(i);
      }
    }
    int id = static_cast<int>(lines_.size()); line_id_[l.name] = id; lines_.push_back(std::move(l));
    line_canonical_name_[lines_.back().name]=lines_.back().name;
    line_depth_.push_back(0);
    line_points_.emplace_back();direction_parent_.push_back(id);
    direction_rank_.push_back(0);direction_parity_.push_back(0);return id;
  }
  void add_circle(Circle c,int center_id=-1) {
    if (circle_id_.count(c.name)) throw std::runtime_error("duplicate circle: " + c.name);
    circle_id_[c.name] = static_cast<int>(circles_.size()); circles_.push_back(std::move(c));
    circle_depth_.push_back(0);
    circle_points_.emplace_back();
    circle_center_ids_.push_back(center_id);
  }
  int segment(int a, int b) {
    if (a == b) throw std::runtime_error("zero segment has no angle");
    if (a > b) std::swap(a, b);
    auto key = std::make_pair(a, b); auto it = segment_line_.find(key);
    if (it != segment_line_.end()) return it->second;
    std::string n = "@" + points_[a].name + points_[b].name;
    int id = add_line(through(n, points_[a], points_[b], "segment"));
    auto&on=line_points_[static_cast<std::size_t>(id)];
    if(std::find(on.begin(),on.end(),a)==on.end())on.push_back(a);
    if(std::find(on.begin(),on.end(),b)==on.end())on.push_back(b);
    segment_line_[key] = id; return id;
  }
  static AngleSystem::Coeff equation(std::initializer_list<std::pair<int,int>> xs) {
    std::vector<std::pair<int,int>> terms(xs);std::sort(terms.begin(),terms.end());
    AngleSystem::Coeff c;
    for(auto [variable,coefficient]:terms){
      if(!c.empty()&&c.back().first==variable)c.back().second+=coefficient;
      else c.push_back({variable,coefficient});
      if(c.back().second==0)c.pop_back();
    }
    return c;
  }
  std::pair<int,int> direction_find(int x){
    if(direction_parent_[static_cast<std::size_t>(x)]==x)return {x,0};
    auto [root,parity]=direction_find(direction_parent_[static_cast<std::size_t>(x)]);
    direction_parity_[static_cast<std::size_t>(x)]^=parity;
    direction_parent_[static_cast<std::size_t>(x)]=root;
    return {root,direction_parity_[static_cast<std::size_t>(x)]};
  }
  void direction_union(int x,int y,int parity){
    auto [rx,px]=direction_find(x);auto [ry,py]=direction_find(y);if(rx==ry)return;
    if(direction_rank_[static_cast<std::size_t>(rx)]<direction_rank_[static_cast<std::size_t>(ry)]){
      std::swap(rx,ry);std::swap(px,py);
    }
    direction_parent_[static_cast<std::size_t>(ry)]=rx;
    direction_parity_[static_cast<std::size_t>(ry)]=px^py^parity;
    if(direction_rank_[static_cast<std::size_t>(rx)]==direction_rank_[static_cast<std::size_t>(ry)])
      ++direction_rank_[static_cast<std::size_t>(rx)];
  }
  bool direction_known(int x,int y,int parity){
    auto [rx,px]=direction_find(x);auto [ry,py]=direction_find(y);
    return rx==ry&&((px^py)==parity);
  }
  bool parallel_fact(int x, int y, const std::string& why) {
    bool known=direction_known(x,y,0);
    direction_union(x,y,0);
    bool added=angles_.add(equation({{x,1},{y,-1}}), 0, 1, why);
    return !known||added;
  }
  bool perpendicular_fact(int x, int y, const std::string& why) {
    bool known=direction_known(x,y,1);
    direction_union(x,y,1);
    bool added=angles_.add(equation({{x,1},{y,-1}}), 1, 2, why);
    return !known||added;
  }
  void circumcenter_angle_fact(int o,int a,int b,int c,const std::string& why){
    // angle(ACB)-angle(OAB)=90 degrees in the line-angle convention used
    // here. Written without dividing a relation.
    angles_.add(equation({{segment(b,c),1},{segment(a,c),-1},
                          {segment(a,b),-1},{segment(a,o),1}}),1,2,why);
  }
  void incidence(int p, int l, const std::string& why) {
    auto known=line_points_[static_cast<std::size_t>(l)];
    for (int q : known)
      if (p != q) parallel_fact(segment(p, q), l, why);
    auto& on = line_points_[static_cast<std::size_t>(l)];
    if (std::find(on.begin(), on.end(), p) == on.end()) on.push_back(p);
  }
  void inherit_collinearity(int p,int a,int b,const std::string& why){
    segment(a,b);
    std::size_t count=line_points_.size();
    for(std::size_t l=0;l<count;++l){const auto&on=line_points_[l];
      if(std::find(on.begin(),on.end(),a)!=on.end()&&
         std::find(on.begin(),on.end(),b)!=on.end())incidence(p,static_cast<int>(l),why);
    }
  }
  int perpendicular_bisector_locus(int a,int b){
    auto key=lenkey(a,b);auto known=perpendicular_bisector_loci_.find(key);
    if(known!=perpendicular_bisector_loci_.end())return known->second;
    for(const auto&pb:perpendicular_bisectors_)
      if(lenkey(pb.a,pb.b)==key){perpendicular_bisector_loci_[key]=pb.line;return pb.line;}
    Point m{"",(points_[a].x+points_[b].x)/2,(points_[a].y+points_[b].y)/2,""};
    Line base=through("",points_[a],points_[b],"");
    std::string name="@perp_bisector("+points_[key.first].name+","+points_[key.second].name+")";
    int line=add_line({name,base.b,-base.a,-(base.b*m.x-base.a*m.y),"circumcenter locus"});
    perpendicular_fact(line,segment(a,b),"circumcenter perpendicular-bisector locus");
    perpendicular_bisectors_.push_back({line,a,b});
    perpendicular_bisector_loci_[key]=line;return line;
  }
  void register_center_loci(){
    // Circumcenters sharing two defining points lie on one canonical
    // perpendicular bisector. Registering incidence avoids rediscovering this
    // definition-level locus from each pair of centers.
    for(const auto&f:circumcenter_facts_)for(auto [a,b]:
        {std::pair{f.a,f.b},std::pair{f.a,f.c},std::pair{f.b,f.c}}){
      int line=perpendicular_bisector_locus(a,b);
      const auto&on=line_points_[static_cast<std::size_t>(line)];
      if(std::find(on.begin(),on.end(),f.center)==on.end())
        incidence(f.center,line,"circumcenter lies on perpendicular bisector");
    }
    // The midpoint of a chord is also on its perpendicular bisector. Keeping
    // this as explicit incidence is important for maximal detected lines such
    // as M(A,B), O(A,B,*), H(A,B,O(A,B,*)); otherwise a large configuration
    // asks the modular angle backend to reconstruct a needlessly indirect proof.
    for(const auto&f:midpoint_facts_){
      auto it=perpendicular_bisector_loci_.find(lenkey(f.a,f.b));
      if(it==perpendicular_bisector_loci_.end())continue;
      const auto&on=line_points_[static_cast<std::size_t>(it->second)];
      if(std::find(on.begin(),on.end(),f.midpoint)==on.end())
        incidence(f.midpoint,it->second,"midpoint lies on perpendicular bisector");
    }

    // Snapshot point-to-line incidence after the circumcenter loci above. For
    // H(A,P,Q), every known carrier of P,Q determines the same altitude through
    // A. The (A, carrier) key makes all such orthocenters share that altitude.
    std::vector<std::vector<int>> point_lines(points_.size());
    for(std::size_t line=0;line<line_points_.size();++line)
      for(int p:line_points_[line])point_lines[static_cast<std::size_t>(p)].push_back(static_cast<int>(line));
    std::map<std::pair<int,int>,int> altitude_loci;
    for(const auto&f:orthocenter_facts_)for(auto [apex,b,c]:
        {std::array{f.a,f.b,f.c},std::array{f.b,f.a,f.c},std::array{f.c,f.a,f.b}}){
      for(int carrier:point_lines[static_cast<std::size_t>(b)]){
        const auto&base_points=line_points_[static_cast<std::size_t>(carrier)];
        if(std::find(base_points.begin(),base_points.end(),c)==base_points.end())continue;
        auto key=std::make_pair(apex,carrier);int altitude=-1;
        if(auto it=altitude_loci.find(key);it!=altitude_loci.end())altitude=it->second;
        else {
          int best_score=-1;
          for(int candidate:point_lines[static_cast<std::size_t>(apex)])
            if(candidate!=carrier&&direction_known(candidate,carrier,1)){
              // Prefer a named/derived carrier over the private two-point
              // segment from this very orthocenter construction. This merges
              // nested centers into an already known locus such as the
              // perpendicular bisector containing O(A,B,*) and M(A,B).
              int score=static_cast<int>(line_points_[static_cast<std::size_t>(candidate)].size());
              if(lines_[static_cast<std::size_t>(candidate)].origin!="segment")score+=1000000;
              if(score>best_score){best_score=score;altitude=candidate;}
            }
          if(altitude<0){
            const Line&base=lines_[static_cast<std::size_t>(carrier)];
            std::string name="@alt("+points_[apex].name+","+std::to_string(carrier)+")";
            altitude=add_line({name,base.b,-base.a,
              -base.b*points_[apex].x+base.a*points_[apex].y,"orthocenter locus"});
            incidence(apex,altitude,"orthocenter altitude through vertex");
            perpendicular_fact(altitude,carrier,"orthocenter altitude locus");
          }
          altitude_loci[key]=altitude;
        }
        const auto&on=line_points_[static_cast<std::size_t>(altitude)];
        if(std::find(on.begin(),on.end(),f.center)==on.end())
          incidence(f.center,altitude,"orthocenter lies on altitude");
      }
    }
  }
  void register_incenter_loci(){
    // Internal angle bisectors depend on rays, not just their unoriented carrier
    // lines. Build a small directed-ray union-find from certified midpoint
    // betweenness: if M is the midpoint of AB, rays AM and AB coincide (and so
    // do BM and BA). This distinguishes the internal and external bisectors
    // without using numerical side tests or illegally halving an angle relation.
    std::map<std::pair<int,int>,int> ray_node;std::vector<int> parent;
    auto get=[&](int a,int b){auto [it,inserted]=ray_node.emplace(std::pair{a,b},static_cast<int>(parent.size()));if(inserted)parent.push_back(it->second);return it->second;};
    auto find=[&](int x){while(parent[static_cast<std::size_t>(x)]!=x){parent[static_cast<std::size_t>(x)]=parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];x=parent[static_cast<std::size_t>(x)];}return x;};
    auto unite=[&](int a,int b){a=find(a);b=find(b);if(a!=b)parent[static_cast<std::size_t>(b)]=a;};
    for(const auto&f:midpoint_facts_){unite(get(f.a,f.midpoint),get(f.a,f.b));unite(get(f.b,f.midpoint),get(f.b,f.a));}
    std::map<std::tuple<int,int,int>,std::vector<int>> loci;
    auto add_locus=[&](int vertex,int side1,int side2,int center){int r1=find(get(vertex,side1)),r2=find(get(vertex,side2));if(r2<r1)std::swap(r1,r2);loci[{vertex,r1,r2}].push_back(center);};
    for(const auto&f:incenter_facts_){add_locus(f.a,f.b,f.c,f.center);add_locus(f.b,f.a,f.c,f.center);add_locus(f.c,f.a,f.b,f.center);}
    for(const auto&[key,centers]:loci)if(centers.size()>1){int vertex=std::get<0>(key),carrier=segment(vertex,centers[0]);for(std::size_t i=1;i<centers.size();++i)incidence(centers[i],carrier,"incenters on the same certified internal angle-bisector ray locus");}
  }
  void circle_incidence(int p, int c) {
    auto& on = circle_points_[static_cast<std::size_t>(c)];
    if (std::find(on.begin(), on.end(), p) == on.end()) on.push_back(p);
  }
  static std::pair<int,int> lenkey(int a, int b) { if (a>b) std::swap(a,b); return std::make_pair(a,b); }
  bool equal_length(int a, int b, int c, int d, const std::string&) {
    auto x=lenkey(a,b), y=lenkey(c,d); if (y<x) std::swap(x,y);
    return equal_lengths_.insert({x,y}).second;
  }
  bool register_line_reflection(int image,int source,int line,const std::string&why){
    if(image==source)return false;
    bool changed=perpendicular_fact(segment(source,image),line,why);
    for(int axis_point:line_points_[static_cast<std::size_t>(line)])
      if(axis_point!=source&&axis_point!=image){
        changed|=equal_length(axis_point,source,axis_point,image,why+" axis distances");
        changed|=angles_.add(equation({{segment(axis_point,source),1},
                                      {segment(axis_point,image),1},{line,-2}}),
                             0,1,why+" angle symmetry");
      }
    for(const auto&f:line_reflection_facts_)
      if(f.line==line&&((f.image==image&&f.source==source)||
                       (f.image==source&&f.source==image)))return changed;
    line_reflection_facts_.push_back({image,source,line});
    return true;
  }
  bool length_equal(int a,int b,int c,int d) const {
    auto x=lenkey(a,b), y=lenkey(c,d); if (y<x) std::swap(x,y);
    if(x==y||equal_lengths_.count({x,y}))return true;
    std::map<std::pair<int,int>,std::vector<std::pair<int,int>>> graph;
    for(const auto&f:equal_lengths_){graph[f.first].push_back(f.second);graph[f.second].push_back(f.first);}
    std::set<std::pair<int,int>> seen{x};std::vector<std::pair<int,int>> todo{x};
    while(!todo.empty()){auto u=todo.back();todo.pop_back();for(auto v:graph[u]){
      if(v==y)return true;
      if(seen.insert(v).second)todo.push_back(v);}}
    return false;
  }
  bool prove_equal_distance(int a,int b,int c,int d){
    if(length_equal(a,b,c,d))return true;
    int vertex=-1,x=-1,y=-1;
    if(a==c){vertex=a;x=b;y=d;}else if(a==d){vertex=a;x=b;y=c;}
    else if(b==c){vertex=b;x=a;y=d;}else if(b==d){vertex=b;x=a;y=c;}
    if(vertex<0||x==y)return false;
    auto relation=equation({{segment(vertex,x),1},{segment(vertex,y),1},{segment(x,y),-2}});
    if(!angles_.proves(relation,0,1))return false;
    equal_length(vertex,x,vertex,y,"converse isosceles triangle theorem");return true;
  }
  void register_midpoint_fact(int midpoint,int a,int b,const std::string&why){
    for(const auto&f:midpoint_facts_)if(f.midpoint==midpoint&&lenkey(f.a,f.b)==lenkey(a,b))return;
    // Two midpoints on sides from a shared vertex give the opposite-side
    // parallel. Keeping this here lets constructed midpoints and point
    // reflections participate in exactly the same theorem closure.
    for(const auto&f:midpoint_facts_){int shared=-1;
      if(a==f.a||a==f.b)shared=a;else if(b==f.a||b==f.b)shared=b;if(shared<0)continue;
      int p=a==shared?b:a,q=f.a==shared?f.b:f.a;
      if(p!=q)parallel_fact(segment(midpoint,f.midpoint),segment(p,q),"triangle midline theorem "+why+","+points_[f.midpoint].name);
    }
    midpoint_facts_.push_back({midpoint,a,b});
  }
  void add_cyclic(int a, int b, int c, int d, const std::string& why) {
    std::array<int,4> q{a,b,c,d}; auto sorted=q; std::sort(sorted.begin(),sorted.end());
    for (auto old:cyclic_facts_) { std::sort(old.begin(),old.end()); if(old==sorted) return; }
    cyclic_facts_.push_back(q);
    int ab=segment(a,b), cd=segment(c,d), ad=segment(a,d), bc=segment(b,c);
    int ac=segment(a,c), bd=segment(b,d);
    angles_.add(equation({{ab,1},{cd,1},{ad,-1},{bc,-1}}),0,1,why+" [cyclic pair-sum 1]");
    angles_.add(equation({{ab,1},{cd,1},{ac,-1},{bd,-1}}),0,1,why+" [cyclic pair-sum 2]");
  }
  bool proves_collinear(const std::vector<int>& p, std::set<int>* why=nullptr) {
    if (p.size()<3) return false;
    for(const auto&on:line_points_){bool all=true;
      for(int x:p)if(std::find(on.begin(),on.end(),x)==on.end()){all=false;break;}
      if(all)return true;
    }
    int base=segment(p[0],p[1]);
    std::set<int> all;
    for(std::size_t i=2;i<p.size();++i){ std::set<int> w;
      if(!angles_.proves(equation({{base,1},{segment(p[0],p[i]),-1}}),0,1,&w)) return false;
      all.insert(w.begin(),w.end()); }
    if(why)*why=std::move(all);
    return true;
  }
  bool proves_cyclic(const std::vector<int>& p, std::set<int>* why=nullptr) {
    if(p.size()<4)return false;
    std::set<int> all;
    for(std::size_t i=3;i<p.size();++i){
      int a=p[0],b=p[1],c=p[2],d=p[i]; std::set<int>w;
      int ab=segment(a,b),cd=segment(c,d),ac=segment(a,c),bd=segment(b,d);
      int ad=segment(a,d),bc=segment(b,c);
      auto e12=equation({{ab,1},{cd,1},{ac,-1},{bd,-1}});
      auto e13=equation({{ab,1},{cd,1},{ad,-1},{bc,-1}});
      auto e23=equation({{ac,1},{bd,1},{ad,-1},{bc,-1}});
      if(!angles_.proves(e12,0,1,&w)&&!angles_.proves(e13,0,1,&w)&&
         !angles_.proves(e23,0,1,&w)) return false;
      all.insert(w.begin(),w.end());
    }
    if(why)*why=std::move(all);
    return true;
  }
  std::vector<int> names_to_points(const std::vector<std::string>& a) const {
    std::vector<int> v; for(auto&s:a)v.push_back(pid(s)); return v;
  }

  void initial_triangle(const std::vector<std::string>& t) {
    if(t.size()!=4) throw std::runtime_error("triangle expects three point names");
    std::array<Point,3> p;
    do {
      for(std::size_t i=0;i<3;++i)
        p[i]={t[i+1],random_real(-10,10),random_real(-10,10),"random initial triangle"};
    } while(std::fabs(cross(p[0],p[1],p[2]))<5.0L);
    for(auto& q:p) add_point(std::move(q));
    initial_point_names_.insert(t[1]);initial_point_names_.insert(t[2]);initial_point_names_.insert(t[3]);
    initial_triangles_.push_back({pid(t[1]),pid(t[2]),pid(t[3])});
  }

  void initial_quadrilateral(const std::vector<std::string>& t, bool cyclic) {
    if(t.size()!=5) throw std::runtime_error(t[0]+" expects four point names");
    std::array<Point,4> p;
    for(;;) {
      std::array<long double,4> angle;
      for(auto& a:angle) a=random_real(0,2*PI);
      std::sort(angle.begin(),angle.end());
      long double cx=random_real(-3,3),cy=random_real(-3,3),base=random_real(4,8);
      for(std::size_t i=0;i<4;++i){
        long double r=cyclic?base:base*random_real(0.80L,1.20L);
        p[i]={t[i+1],cx+r*std::cos(angle[i]),cy+r*std::sin(angle[i]),
              cyclic?"random initial cyclic quadrilateral":"random initial quadrilateral"};
      }
      long double sign=0;bool good=true;
      for(int i=0;i<4;++i){long double z=cross(p[i],p[(i+1)%4],p[(i+2)%4]);
        if(std::fabs(z)<1.0L){good=false;break;}if(i==0)sign=z;else if(z*sign<=0){good=false;break;}}
      if(good)break;
    }
    for(auto& q:p) add_point(std::move(q));
    for(std::size_t i=1;i<5;++i)initial_point_names_.insert(t[i]);
    if(cyclic) add_cyclic(pid(t[1]),pid(t[2]),pid(t[3]),pid(t[4]),
                          "initial cyclic quadrilateral");
  }

  void execute(const std::vector<std::string>& t, int line_no) {
    auto need=[&](std::size_t n){if(t.size()!=n)throw std::runtime_error(t[0]+" expects "+std::to_string(n-1)+" arguments");};
    const auto& op=t[0];
    if(record_commands_&&is_construction(op))construction_commands_.push_back(t);
    if(op=="mode"){need(2);prove_mode_=t[1]=="prove";if(t[1]!="prove"&&t[1]!="generate")throw std::runtime_error("mode is generate or prove");}
    else if(op=="option"){
      if(t.size()<3)throw std::runtime_error("option expects a name and value");
      if(t[1]=="symmetry"){
        need(4);symmetry_enabled_=true;symmetry_first_=t[2];symmetry_second_=t[3];
      } else if(t[1]=="construction"){
        need(4);auto construction=canonical_construction_option(t[2]);bool enabled=std::stoi(t[3])!=0;
        if(construction=="all"){
          if(enabled)disabled_constructions_.clear();else disabled_constructions_=automatic_construction_types();
        } else {
          if(!automatic_construction_types().count(construction))throw std::runtime_error("unknown automatic construction "+t[2]);
          if(enabled)disabled_constructions_.erase(construction);else disabled_constructions_.insert(construction);
        }
      } else {
        need(3);if(t[1]=="show_easy")show_easy_=std::stoi(t[2])!=0;else if(t[1]=="show_all_points"){(void)std::stoi(t[2]);}else if(t[1]=="symmetric_coincidences_only")symmetric_only_=std::stoi(t[2])!=0;else if(t[1]=="circle_budget")circle_budget_=std::stoull(t[2]);else if(t[1]=="angle_coefficient_limit"){angle_coefficient_limit_=std::stoll(t[2]);angles_.set_coefficient_limit(angle_coefficient_limit_);}else if(t[1]=="proof_scope"){if(t[2]=="global")ancestry_scope_=false;else if(t[2]=="ancestry")ancestry_scope_=true;else throw std::runtime_error("option proof_scope is global or ancestry");}else if(t[1]=="max_points"){max_points_=std::stoull(t[2]);if(max_points_>5000)throw std::runtime_error("option max_points cannot exceed 5000");if(max_points_&&points_.size()>max_points_)throw std::runtime_error("option max_points is below the number of points already declared");}else if(t[1]=="trials"||t[1]=="seed"){}else throw std::runtime_error("unknown option "+t[1]);
      }
    }
    else if(op=="triangle") initial_triangle(t);
    else if(op=="quadrilateral") initial_quadrilateral(t,false);
    else if(op=="cyclic_quad") initial_quadrilateral(t,true);
    else if(op=="point"){need(4);add_point({t[1],std::stold(t[2]),std::stold(t[3]),"declared"});}
    else if(op=="line"){need(4);int a=pid(t[2]),b=pid(t[3]);int id=add_line(through(t[1],points_[a],points_[b],op));parallel_fact(id,segment(a,b),"definition line("+t[2]+","+t[3]+")");incidence(a,id,"line incidence");incidence(b,id,"line incidence");}
    else if(op=="midpoint"){need(4);int a=pid(t[2]),b=pid(t[3]);if(!add_point({t[1],(points_[a].x+points_[b].x)/2,(points_[a].y+points_[b].y)/2,op}))return;int m=pid(t[1]);parallel_fact(segment(a,m),segment(a,b),"midpoint collinearity "+t[1]);parallel_fact(segment(b,m),segment(a,b),"midpoint collinearity "+t[1]);inherit_collinearity(m,a,b,"midpoint incidence "+t[1]);equal_length(a,m,m,b,"midpoint lengths");register_midpoint_fact(m,a,b,t[1]);}
    else if(op=="perp_bisector"){need(4);int a=pid(t[2]),b=pid(t[3]);Point m{"",(points_[a].x+points_[b].x)/2,(points_[a].y+points_[b].y)/2,""};Line base=through("",points_[a],points_[b],"");Line l{t[1],base.b,-base.a,-(base.b*m.x-base.a*m.y),op};int id=add_line(l);perpendicular_fact(id,segment(a,b),"perpendicular bisector "+t[1]);perpendicular_bisectors_.push_back({id,a,b});}
    else if(op=="parallel"||op=="perpendicular"){need(4);int p=pid(t[2]),base=lid(t[3]);const Line& q=lines_[base];Line l;if(op=="parallel")l={t[1],q.a,q.b,-q.a*points_[p].x-q.b*points_[p].y,op};else l={t[1],q.b,-q.a,-q.b*points_[p].x+q.a*points_[p].y,op};int id=add_line(l);incidence(p,id,op+" through point");if(op=="parallel")parallel_fact(id,base,"parallel construction "+t[1]);else perpendicular_fact(id,base,"perpendicular construction "+t[1]);}
    else if(op=="angle_bisector"){need(5);int a=pid(t[2]),b=pid(t[3]),c=pid(t[4]);long double ux=points_[a].x-points_[b].x,uy=points_[a].y-points_[b].y,vx=points_[c].x-points_[b].x,vy=points_[c].y-points_[b].y;long double un=std::hypotl(ux,uy),vn=std::hypotl(vx,vy);Point q{"",points_[b].x+ux/un+vx/vn,points_[b].y+uy/un+vy/vn,""};int id=add_line(through(t[1],points_[b],q,op));incidence(b,id,"angle bisector through vertex");angles_.add(equation({{id,2},{segment(a,b),-1},{segment(b,c),-1}}),0,1,"angle bisector "+t[1]);angle_bisector_line_facts_.push_back({id,b});}
    else if(op=="reflection_line"){need(4);int p=pid(t[2]),l=lid(t[3]);auto&q=lines_[l];long double d=q.a*points_[p].x+q.b*points_[p].y+q.c;if(!add_point({t[1],points_[p].x-2*q.a*d,points_[p].y-2*q.b*d,op}))return;int x=pid(t[1]);register_line_reflection(x,p,l,"line reflection "+t[1]);}
    else if(op=="reflection_point"){need(4);int p=pid(t[2]),o=pid(t[3]);if(!add_point({t[1],2*points_[o].x-points_[p].x,2*points_[o].y-points_[p].y,op}))return;int x=pid(t[1]);parallel_fact(segment(p,o),segment(p,x),"point reflection collinearity "+t[1]);inherit_collinearity(x,p,o,"point reflection incidence "+t[1]);equal_length(p,o,o,x,"point reflection lengths");register_midpoint_fact(o,p,x,"point reflection "+t[1]);}
    else if(op=="foot"){need(4);int p=pid(t[2]),l=lid(t[3]);auto&q=lines_[l];long double d=q.a*points_[p].x+q.b*points_[p].y+q.c;if(!add_point({t[1],points_[p].x-q.a*d,points_[p].y-q.b*d,op}))return;int x=pid(t[1]);incidence(x,l,"foot incidence "+t[1]);perpendicular_fact(segment(p,x),l,"foot "+t[1]);foot_facts_.push_back({x,p,l});}
    else if(op=="intersection_ll"){need(4);int a=lid(t[2]),b=lid(t[3]);if(!add_point(intersect(lines_[a],lines_[b],t[1],op)))return;int x=pid(t[1]);incidence(x,a,"intersection incidence "+t[1]+" on "+t[2]);incidence(x,b,"intersection incidence "+t[1]+" on "+t[3]);}
    else if(op=="circumcenter"){need(5);int a=pid(t[2]),b=pid(t[3]),c=pid(t[4]);if(!add_point(circumcenter(points_[a],points_[b],points_[c],t[1],op)))return;int o=pid(t[1]);equal_length(o,a,o,b,"circumcenter radii");equal_length(o,a,o,c,"circumcenter radii");equal_length(o,b,o,c,"circumcenter radii");circumcenter_angle_fact(o,a,b,c,"circumcenter angle theorem at "+t[1]);circumcenter_angle_fact(o,b,c,a,"circumcenter angle theorem at "+t[1]);circumcenter_angle_fact(o,c,a,b,"circumcenter angle theorem at "+t[1]);circumcenter_facts_.push_back({o,a,b,c});}
    else if(op=="orthocenter"){need(5);int a=pid(t[2]),b=pid(t[3]),c=pid(t[4]);Line bc=through("",points_[b],points_[c],""),ac=through("",points_[a],points_[c],"");Line ha{"",bc.b,-bc.a,-bc.b*points_[a].x+bc.a*points_[a].y,""},hb{"",ac.b,-ac.a,-ac.b*points_[b].x+ac.a*points_[b].y,""};if(!add_point(intersect(ha,hb,t[1],op)))return;int h=pid(t[1]);perpendicular_fact(segment(a,h),segment(b,c),"orthocenter altitude 1 "+t[1]);perpendicular_fact(segment(b,h),segment(a,c),"orthocenter altitude 2 "+t[1]);perpendicular_fact(segment(c,h),segment(a,b),"orthocenter closure "+t[1]);orthocenter_facts_.push_back({h,a,b,c});}
    else if(op=="incenter"){need(5);int a=pid(t[2]),b=pid(t[3]),c=pid(t[4]);long double la=std::sqrt(dist2(points_[b],points_[c])),lb=std::sqrt(dist2(points_[a],points_[c])),lc=std::sqrt(dist2(points_[a],points_[b])),s=la+lb+lc;if(!add_point({t[1],(la*points_[a].x+lb*points_[b].x+lc*points_[c].x)/s,(la*points_[a].y+lb*points_[b].y+lc*points_[c].y)/s,op}))return;int i=pid(t[1]);angles_.add(equation({{segment(a,i),2},{segment(a,b),-1},{segment(a,c),-1}}),0,1,"incenter bisector at "+t[2]);angles_.add(equation({{segment(b,i),2},{segment(a,b),-1},{segment(b,c),-1}}),0,1,"incenter bisector at "+t[3]);angles_.add(equation({{segment(c,i),2},{segment(a,c),-1},{segment(b,c),-1}}),0,1,"incenter closure at "+t[4]);
      // The doubled bisector equations do not distinguish the internal and
      // external branches modulo 180 degrees. `incenter` certifies the internal
      // choice, whose three undivided cross-angle sums are 90 degrees. These
      // facts resolve that branch without ever dividing an arbitrary relation.
      angles_.add(equation({{segment(b,i),1},{segment(c,i),1},{segment(b,c),-1},{segment(a,i),-1}}),1,2,"internal incenter angle sum opposite "+t[2]);
      angles_.add(equation({{segment(a,i),1},{segment(c,i),1},{segment(a,c),-1},{segment(b,i),-1}}),1,2,"internal incenter angle sum opposite "+t[3]);
      angles_.add(equation({{segment(a,i),1},{segment(b,i),1},{segment(a,b),-1},{segment(c,i),-1}}),1,2,"internal incenter angle sum opposite "+t[4]);
      incenter_facts_.push_back({i,a,b,c});angle_bisector_diagonals_.push_back({a,i});angle_bisector_diagonals_.push_back({b,i});angle_bisector_diagonals_.push_back({c,i});}
    else if(op=="circle"){need(4);int o=pid(t[2]),p=pid(t[3]);add_circle({t[1],points_[o],dist2(points_[o],points_[p]),op},o);circle_incidence(p,cid(t[1]));}
    else if(op=="circumcircle"){need(5);int a=pid(t[2]),b=pid(t[3]),c=pid(t[4]);Point o=circumcenter(points_[a],points_[b],points_[c],"@"+t[1],op);add_circle({t[1],o,dist2(o,points_[a]),op});int z=cid(t[1]);circle_incidence(a,z);circle_incidence(b,z);circle_incidence(c,z);}
    else if(op=="incircle"){need(6);int i=pid(t[2]),a=pid(t[3]),b=pid(t[4]),c=pid(t[5]);Line ab=through("",points_[a],points_[b],"");long double r=ab.a*points_[i].x+ab.b*points_[i].y+ab.c;add_circle({t[1],points_[i],r*r,op},i);(void)c;}
    else if(op=="intersection_lc_known"){need(5);int l=lid(t[2]),c=cid(t[3]),k=pid(t[4]);auto&ln=lines_[l];auto&cc=circles_[c];if(std::fabs(ln.a*points_[k].x+ln.b*points_[k].y+ln.c)>EPS*10||!near(dist2(cc.center,points_[k]),cc.r2,10))throw std::runtime_error(t[4]+" is not a known intersection of "+t[2]+" and "+t[3]);long double dx=ln.b,dy=-ln.a;long double vx=points_[k].x-cc.center.x,vy=points_[k].y-cc.center.y;long double tt=-2*dot(vx,vy,dx,dy);if(std::fabs(tt)<=EPS)throw std::runtime_error("known line-circle intersection is tangent; no distinct second point");if(!add_point({t[1],points_[k].x+tt*dx,points_[k].y+tt*dy,op}))return;int x=pid(t[1]);incidence(k,l,"known line-circle incidence");incidence(x,l,"line-circle second incidence");circle_incidence(k,c);circle_incidence(x,c);parallel_fact(segment(k,x),l,"line-circle known-root incidence");}
    else if(op=="intersection_cc_known"){need(5);int c1=cid(t[2]),c2=cid(t[3]),k=pid(t[4]);auto&a=circles_[c1];auto&b=circles_[c2];if(!near(dist2(a.center,points_[k]),a.r2,10)||!near(dist2(b.center,points_[k]),b.r2,10))throw std::runtime_error(t[4]+" is not on both circles");Line axis=through("",a.center,b.center,"");long double d=axis.a*points_[k].x+axis.b*points_[k].y+axis.c;if(std::fabs(d)<=EPS)throw std::runtime_error("known circle-circle intersection is tangent; no distinct second point");if(!add_point({t[1],points_[k].x-2*axis.a*d,points_[k].y-2*axis.b*d,op}))return;int x=pid(t[1]);circle_incidence(k,c1);circle_incidence(x,c1);circle_incidence(k,c2);circle_incidence(x,c2);
      int o1=circle_center_ids_[static_cast<std::size_t>(c1)],o2=circle_center_ids_[static_cast<std::size_t>(c2)];
      if(o1>=0&&o2>=0&&o1!=o2){int centers=segment(o1,o2);equal_length(o1,k,o1,x,"circle-circle equal radii");equal_length(o2,k,o2,x,"circle-circle equal radii");perpendicular_fact(segment(k,x),centers,"circle-circle common chord perpendicular to center line");angles_.add(equation({{segment(o1,k),1},{segment(o1,x),1},{centers,-2}}),0,1,"circle-circle reflection across center line");angles_.add(equation({{segment(o2,k),1},{segment(o2,x),1},{centers,-2}}),0,1,"circle-circle reflection across center line");}
    }
    else if(op.rfind("prove_",0)==0){goals_.push_back({op.substr(6),std::vector<std::string>(t.begin()+1,t.end())});}
    else throw std::runtime_error("line "+std::to_string(line_no)+": unknown command "+op);
  }

  void normalize_definition_incidences(){
    // Two certified points determine a unique line. Merge every carrier that
    // contains the same point pair; also merge known-parallel carriers through
    // a common point. This turns long maximal collinear sets assembled through
    // line(), foot(), reflection(), and intersections into a direct incidence
    // fact rather than a large angle-lattice certificate.
    std::vector<int> parent(lines_.size());std::iota(parent.begin(),parent.end(),0);
    auto find=[&](int x){while(parent[static_cast<std::size_t>(x)]!=x){parent[static_cast<std::size_t>(x)]=parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];x=parent[static_cast<std::size_t>(x)];}return x;};
    auto unite=[&](int a,int b){a=find(a);b=find(b);if(a!=b)parent[static_cast<std::size_t>(b)]=a;};
    std::map<std::pair<int,int>,int> pair_owner;
    for(std::size_t l=0;l<line_points_.size();++l){auto on=line_points_[l];std::sort(on.begin(),on.end());on.erase(std::unique(on.begin(),on.end()),on.end());
      for(std::size_t i=0;i<on.size();++i)for(std::size_t j=i+1;j<on.size();++j){auto key=std::pair{on[i],on[j]};auto [it,inserted]=pair_owner.emplace(key,static_cast<int>(l));if(!inserted)unite(static_cast<int>(l),it->second);}}
    std::map<std::tuple<int,int,int>,int> direction_point_owner;
    for(std::size_t l=0;l<line_points_.size();++l){auto [root,parity]=direction_find(static_cast<int>(l));for(int p:line_points_[l]){
      auto key=std::tuple{root,parity,p};auto [it,inserted]=direction_point_owner.emplace(key,static_cast<int>(l));if(!inserted)unite(static_cast<int>(l),it->second);}}
    std::map<int,std::set<int>> component_points;std::map<int,std::vector<int>> component_lines;
    for(std::size_t l=0;l<line_points_.size();++l){int root=find(static_cast<int>(l));component_lines[root].push_back(static_cast<int>(l));component_points[root].insert(line_points_[l].begin(),line_points_[l].end());}
    for(const auto&[root,members]:component_lines){std::vector<int> on(component_points[root].begin(),component_points[root].end());
      for(int line:members)line_points_[static_cast<std::size_t>(line)]=on;}

    // Three certified points determine a unique circle. Equivalent declared or
    // constructed circles therefore share every known incidence.
    std::vector<int> circle_parent(circles_.size());std::iota(circle_parent.begin(),circle_parent.end(),0);
    auto circle_find=[&](int x){while(circle_parent[static_cast<std::size_t>(x)]!=x){circle_parent[static_cast<std::size_t>(x)]=circle_parent[static_cast<std::size_t>(circle_parent[static_cast<std::size_t>(x)])];x=circle_parent[static_cast<std::size_t>(x)];}return x;};
    auto circle_unite=[&](int a,int b){a=circle_find(a);b=circle_find(b);if(a!=b)circle_parent[static_cast<std::size_t>(b)]=a;};
    std::map<std::tuple<int,int,int>,int> triple_owner;
    for(std::size_t c=0;c<circle_points_.size();++c){auto on=circle_points_[c];std::sort(on.begin(),on.end());on.erase(std::unique(on.begin(),on.end()),on.end());
      for(std::size_t i=0;i<on.size();++i)for(std::size_t j=i+1;j<on.size();++j)for(std::size_t k=j+1;k<on.size();++k){auto key=std::tuple{on[i],on[j],on[k]};auto [it,inserted]=triple_owner.emplace(key,static_cast<int>(c));if(!inserted)circle_unite(static_cast<int>(c),it->second);}}
    std::map<int,std::set<int>> component_circle_points;std::map<int,int> component_center;
    for(std::size_t c=0;c<circle_points_.size();++c){int root=circle_find(static_cast<int>(c));component_circle_points[root].insert(circle_points_[c].begin(),circle_points_[c].end());if(circle_center_ids_[c]>=0)component_center.emplace(root,circle_center_ids_[c]);}
    for(std::size_t c=0;c<circle_points_.size();++c){int root=circle_find(static_cast<int>(c));circle_points_[c].assign(component_circle_points[root].begin(),component_circle_points[root].end());if(auto it=component_center.find(root);it!=component_center.end())circle_center_ids_[c]=it->second;}
  }

  void register_affine_facts(){
    // Use two-field barycentric coordinates only when the initial configuration is
    // one free triangle. A generic quadrilateral has a fourth independent
    // coordinate choice which this deliberately small affine layer does not
    // attempt to encode.
    int triangle_commands=0;std::vector<std::string> initial;
    for(const auto&t:construction_commands_){if(t[0]=="triangle"){++triangle_commands;initial=t;}
      else if(t[0]=="quadrilateral"||t[0]=="cyclic_quad"||t[0]=="point")return;}
    if(triangle_commands!=1||initial.size()!=4)return;

    std::vector<std::optional<AffineVector>> point(points_.size());
    std::vector<std::optional<AffineVector>> line(lines_.size());
    point[static_cast<std::size_t>(pid(initial[1]))]=AffineVector{AffineScalar(1),AffineScalar(0),AffineScalar(0)};
    point[static_cast<std::size_t>(pid(initial[2]))]=AffineVector{AffineScalar(0),AffineScalar(1),AffineScalar(0)};
    point[static_cast<std::size_t>(pid(initial[3]))]=AffineVector{AffineScalar(0),AffineScalar(0),AffineScalar(1)};
    auto assign_point=[&](const std::string&name,const std::optional<AffineVector>&value){if(!value)return;int id=pid(name);if(!point[static_cast<std::size_t>(id)])point[static_cast<std::size_t>(id)]=*value;};
    auto assign_line=[&](const std::string&name,const std::optional<AffineVector>&value){if(!value)return;int id=lid(name);if(static_cast<std::size_t>(id)>=line.size())line.resize(lines_.size());if(!line[static_cast<std::size_t>(id)])line[static_cast<std::size_t>(id)]=affine_normalize_projective(*value);};
    auto point_value=[&](const std::string&name)->std::optional<AffineVector>{return point[static_cast<std::size_t>(pid(name))];};
    auto line_value=[&](const std::string&name)->std::optional<AffineVector>{int id=lid(name);if(static_cast<std::size_t>(id)>=line.size())return std::nullopt;return line[static_cast<std::size_t>(id)];};

    // Metric closure can discover an affine fact after the first affine pass;
    // notably, the intersection of a segment and its perpendicular bisector is
    // a midpoint even when that point was originally named as a foot. Seed such
    // learned midpoints before replaying subsequent affine constructions.
    bool seeded=true;while(seeded){seeded=false;for(const auto&f:midpoint_facts_)
      if(!point[static_cast<std::size_t>(f.midpoint)]&&point[static_cast<std::size_t>(f.a)]&&point[static_cast<std::size_t>(f.b)]){
        AffineVector m;for(int i=0;i<3;++i)m[static_cast<std::size_t>(i)]=
          ((*point[static_cast<std::size_t>(f.a)])[static_cast<std::size_t>(i)]+(*point[static_cast<std::size_t>(f.b)])[static_cast<std::size_t>(i)])/2;
        point[static_cast<std::size_t>(f.midpoint)]=m;seeded=true;
      }}

    const AffineVector infinity{AffineScalar(1),AffineScalar(1),AffineScalar(1)};
    for(const auto&t:construction_commands_){const auto&op=t[0];
      if(op=="line"){auto a=point_value(t[2]),b=point_value(t[3]);if(a&&b){auto l=affine_cross(*a,*b);if(!affine_zero(l))assign_line(t[1],l);}}
      else if(op=="midpoint"){auto a=point_value(t[2]),b=point_value(t[3]);if(a&&b){AffineVector m;for(int i=0;i<3;++i)m[static_cast<std::size_t>(i)]=((*a)[static_cast<std::size_t>(i)]+(*b)[static_cast<std::size_t>(i)])/2;assign_point(t[1],m);}}
      else if(op=="reflection_point"){auto p=point_value(t[2]),o=point_value(t[3]);if(p&&o){AffineVector x;for(int i=0;i<3;++i)x[static_cast<std::size_t>(i)]=2*(*o)[static_cast<std::size_t>(i)]-(*p)[static_cast<std::size_t>(i)];assign_point(t[1],x);}}
      else if(op=="parallel"){auto p=point_value(t[2]),base=line_value(t[3]);if(p&&base){auto direction=affine_cross(*base,infinity);auto l=affine_cross(*p,direction);if(!affine_zero(l))assign_line(t[1],l);}}
      else if(op=="intersection_ll"){auto a=line_value(t[2]),b=line_value(t[3]);if(a&&b)assign_point(t[1],affine_normalize_point(affine_cross(*a,*b)));}
    }

    std::map<AffineVector,std::vector<int>> at;
    std::vector<int> known_points;for(std::size_t i=0;i<point.size();++i)if(point[i]){known_points.push_back(static_cast<int>(i));at[*point[i]].push_back(static_cast<int>(i));}
    // Discover every already constructed point which is the certified midpoint
    // of two affine-known points. This includes nested midpoint/reflection
    // identities not visible from the local construction that named the point.
    for(std::size_t i=0;i<known_points.size();++i)for(std::size_t j=i+1;j<known_points.size();++j){int a=known_points[i],b=known_points[j];AffineVector m;
      for(int k=0;k<3;++k)m[static_cast<std::size_t>(k)]=((*point[static_cast<std::size_t>(a)])[static_cast<std::size_t>(k)]+(*point[static_cast<std::size_t>(b)])[static_cast<std::size_t>(k)])/2;
      auto found=at.find(m);if(found==at.end())continue;for(int center:found->second)if(center!=a&&center!=b){equal_length(a,center,center,b,"affine midpoint identity");register_midpoint_fact(center,a,b,"affine midpoint identity");inherit_collinearity(center,a,b,"affine midpoint incidence");}}

    // Group affine-known points by certified carrier. Only groups of at least three
    // add information; pair-only carriers are left unmaterialized.
    std::map<AffineVector,std::set<int>> carriers;
    for(std::size_t i=0;i<known_points.size();++i)for(std::size_t j=i+1;j<known_points.size();++j){int a=known_points[i],b=known_points[j];auto l=affine_cross(*point[static_cast<std::size_t>(a)],*point[static_cast<std::size_t>(b)]);if(!affine_zero(l)){l=affine_normalize_projective(std::move(l));carriers[l].insert(a);carriers[l].insert(b);}}
    for(const auto&[_,on]:carriers)if(on.size()>=3){auto it=on.begin();int a=*it++,b=*it++;int carrier=segment(a,b);for(;it!=on.end();++it)incidence(*it,carrier,"affine collinearity certificate");}

    // Recover equations for every existing carrier containing two known
    // affine points, then register direction classes with the angle engine.
    line.resize(lines_.size());
    for(std::size_t l=0;l<line_points_.size();++l)if(!line[l]){std::vector<int> on;for(int p:line_points_[l])if(point[static_cast<std::size_t>(p)])on.push_back(p);if(on.size()>=2){auto value=affine_cross(*point[static_cast<std::size_t>(on[0])],*point[static_cast<std::size_t>(on[1])]);if(!affine_zero(value))line[l]=affine_normalize_projective(std::move(value));}}
    std::map<AffineVector,std::vector<int>> directions;
    for(std::size_t l=0;l<line.size();++l)if(line[l]){auto direction=affine_normalize_projective(affine_cross(*line[l],infinity));if(!affine_zero(direction))directions[direction].push_back(static_cast<int>(l));}
    for(const auto&[_,same]:directions)if(same.size()>=2)for(std::size_t i=1;i<same.size();++i)parallel_fact(same[0],same[i],"affine parallelism certificate");
  }

  void register_formal_affine_facts(){
    // Unlike the three-coordinate affine replay above, this version permits a
    // metrically constructed point (such as a foot) to remain an independent
    // formal vector. Midpoints and parallel-line intersections can still
    // cancel that unknown vector and prove universal affine consequences.
    using Expr=std::map<int,AffineScalar>;
    struct FormalLine {Expr point,direction;};
    auto add=[](Expr a,const Expr&b,AffineScalar factor=AffineScalar(1)){
      for(const auto&[variable,coefficient]:b){auto value=a[variable]+factor*coefficient;
        if(value==AffineScalar(0))a.erase(variable);else a[variable]=value;}return a;};
    auto scaled=[](Expr a,AffineScalar factor){for(auto&[_,coefficient]:a)coefficient=coefficient*factor;return a;};
    auto direction_key=[&](Expr value)->std::optional<Expr>{
      if(value.empty())return std::nullopt;
      AffineScalar pivot=value.begin()->second;
      if(pivot.a==0||pivot.b==0)return std::nullopt;
      return scaled(std::move(value),pivot.inverse());
    };
    auto difference=[&](const Expr&a,const Expr&b){return add(a,b,AffineScalar(-1));};

    std::vector<Expr> point(points_.size());
    for(std::size_t i=0;i<point.size();++i)point[i][static_cast<int>(i)]=AffineScalar(1);
    std::vector<std::optional<FormalLine>> formal_line(lines_.size());
    auto solve_intersection=[&](const FormalLine&a,const FormalLine&b)->std::optional<Expr>{
      Expr rhs=difference(b.point,a.point);std::set<int> variables;
      for(const auto&[v,_]:a.direction)variables.insert(v);
      for(const auto&[v,_]:b.direction)variables.insert(v);
      for(const auto&[v,_]:rhs)variables.insert(v);
      auto coefficient=[](const Expr&e,int v){auto it=e.find(v);return it==e.end()?AffineScalar(0):it->second;};
      std::vector<int> ids(variables.begin(),variables.end());
      for(std::size_t i=0;i<ids.size();++i)for(std::size_t j=i+1;j<ids.size();++j){
        AffineScalar a1=coefficient(a.direction,ids[i]),b1=AffineScalar(0)-coefficient(b.direction,ids[i]);
        AffineScalar a2=coefficient(a.direction,ids[j]),b2=AffineScalar(0)-coefficient(b.direction,ids[j]);
        AffineScalar det=a1*b2-a2*b1;if(det.a==0||det.b==0)continue;
        AffineScalar r1=coefficient(rhs,ids[i]),r2=coefficient(rhs,ids[j]);
        AffineScalar t=(r1*b2-r2*b1)/det,s=(a1*r2-a2*r1)/det;bool valid=true;
        for(int v:ids)if(coefficient(a.direction,v)*t-coefficient(b.direction,v)*s!=coefficient(rhs,v)){valid=false;break;}
        if(valid)return add(a.point,scaled(a.direction,t));
      }
      return std::nullopt;
    };

    // Learned midpoint identities may concern points whose original command was
    // metric. Apply them before and after replaying parallel intersections.
    for(int round=0;round<4;++round){
      for(const auto&mf:midpoint_facts_){Expr value=scaled(add(point[static_cast<std::size_t>(mf.a)],point[static_cast<std::size_t>(mf.b)]),AffineScalar(1)/AffineScalar(2));point[static_cast<std::size_t>(mf.midpoint)]=std::move(value);}
      formal_line.assign(lines_.size(),std::nullopt);
      for(const auto&t:construction_commands_){const auto&op=t[0];
        if(op=="line"){int a=pid(t[2]),b=pid(t[3]),l=lid(t[1]);formal_line[static_cast<std::size_t>(l)]=FormalLine{point[static_cast<std::size_t>(a)],difference(point[static_cast<std::size_t>(b)],point[static_cast<std::size_t>(a)])};}
        else if(op=="parallel"){int p=pid(t[2]),l=lid(t[1]),base=lid(t[3]);if(formal_line[static_cast<std::size_t>(base)])formal_line[static_cast<std::size_t>(l)]=FormalLine{point[static_cast<std::size_t>(p)],formal_line[static_cast<std::size_t>(base)]->direction};}
        else if(op=="intersection_ll"){int x=pid(t[1]),a=lid(t[2]),b=lid(t[3]);if(formal_line[static_cast<std::size_t>(a)]&&formal_line[static_cast<std::size_t>(b)])if(auto value=solve_intersection(*formal_line[static_cast<std::size_t>(a)],*formal_line[static_cast<std::size_t>(b)]))point[static_cast<std::size_t>(x)]=std::move(*value);}
      }
    }

    // Reconstruct every existing carrier from any two formal points on it and
    // merge carriers with proportional formal direction vectors.
    std::map<Expr,std::vector<int>> directions;
    for(std::size_t l=0;l<line_points_.size();++l){const auto&on=line_points_[l];if(on.size()<2)continue;
      if(auto key=direction_key(difference(point[static_cast<std::size_t>(on[1])],point[static_cast<std::size_t>(on[0])])) )directions[*key].push_back(static_cast<int>(l));}
    for(const auto&[_,same]:directions)if(same.size()>=2)
      for(std::size_t i=1;i<same.size();++i)parallel_fact(same[0],same[i],"formal affine parallelism certificate");
  }

  bool register_orthocenter_closure(){
    // In any triangle, two concurrent altitudes determine the orthocenter, so
    // the third vertex-to-center line is perpendicular to the opposite side.
    // Restrict candidate quadruples to proof goals and numerically detected
    // circles; this avoids an O(n^4) scan while still feeding discovered facts
    // back into the ordinary fixed-point angle chase.
    std::set<std::array<int,4>> candidates;
    for(const auto&g:goals_)if(g.args.size()==4){std::array<int,4> q;
      for(int i=0;i<4;++i)q[static_cast<std::size_t>(i)]=pid(g.args[static_cast<std::size_t>(i)]);
      std::sort(q.begin(),q.end());if(std::adjacent_find(q.begin(),q.end())==q.end())candidates.insert(q);}
    for(const auto&candidate:circle_cache_)if(candidate.points.size()>=4)
      for(std::size_t i=3;i<candidate.points.size();++i){std::array<int,4> q{candidate.points[0],candidate.points[1],candidate.points[2],candidate.points[i]};std::sort(q.begin(),q.end());candidates.insert(q);}

    bool changed=false;
    for(const auto&q:candidates)for(int hi=0;hi<4;++hi){int h=q[static_cast<std::size_t>(hi)];std::array<int,3> v;int at=0;
      for(int i=0;i<4;++i)if(i!=hi)v[static_cast<std::size_t>(at++)]=q[static_cast<std::size_t>(i)];
      std::array<std::pair<int,int>,3> altitudes{{
        {segment(h,v[0]),segment(v[1],v[2])},
        {segment(h,v[1]),segment(v[0],v[2])},
        {segment(h,v[2]),segment(v[0],v[1])}}};
      int known=0;for(auto [a,b]:altitudes)known+=direction_known(a,b,1)?1:0;
      if(known<2)continue;
      for(auto [a,b]:altitudes)if(!direction_known(a,b,1))
        changed|=perpendicular_fact(a,b,"third-altitude orthocenter closure");
    }
    return changed;
  }

  void register_perpendicular_bisector_midpoints(){
    // The unique intersection of a segment carrier and its perpendicular
    // bisector is the segment midpoint. This includes the familiar fact that
    // the perpendicular foot from a circumcenter to a chord bisects the chord.
    for(const auto&pb:perpendicular_bisectors_)for(int x:line_points_[static_cast<std::size_t>(pb.line)]){
      if(x==pb.a||x==pb.b)continue;
      bool on_base=false;
      for(const auto&carrier:line_points_)if(std::find(carrier.begin(),carrier.end(),pb.a)!=carrier.end()&&
          std::find(carrier.begin(),carrier.end(),pb.b)!=carrier.end()&&std::find(carrier.begin(),carrier.end(),x)!=carrier.end()){on_base=true;break;}
      if(on_base){equal_length(pb.a,x,x,pb.b,"perpendicular-bisector intersection midpoint lengths");register_midpoint_fact(x,pb.a,pb.b,"perpendicular-bisector intersection");}
    }
  }

  void register_derived_perpendicular_bisectors(){
    // Any known line through a segment midpoint and perpendicular to the
    // segment is its perpendicular bisector. Register the locus so reflection
    // transport and equal-distance closure can use it just like an explicitly
    // constructed perpendicular_bisector.
    std::size_t midpoint_count=midpoint_facts_.size();
    for(std::size_t i=0;i<midpoint_count;++i){const auto f=midpoint_facts_[i];int base=segment(f.a,f.b);std::size_t line_count=lines_.size();
      for(std::size_t l=0;l<line_count;++l){const auto&on=line_points_[l];
        if(std::find(on.begin(),on.end(),f.midpoint)==on.end()||!direction_known(static_cast<int>(l),base,1))continue;
        auto key=lenkey(f.a,f.b);perpendicular_bisector_loci_[key]=static_cast<int>(l);
        bool known=false;for(const auto&pb:perpendicular_bisectors_)if(pb.line==static_cast<int>(l)&&lenkey(pb.a,pb.b)==key){known=true;break;}
        if(!known)perpendicular_bisectors_.push_back({static_cast<int>(l),f.a,f.b});
      }
    }
  }

  bool register_projection_midpoints(){
    // Orthogonal projection onto a fixed line is an affine map. Hence it sends
    // the midpoint of AB to the midpoint of the projections of A and B. A
    // point already on the target line is its own projection.
    bool changed=false;std::set<int> target_lines;
    for(const auto&f:foot_facts_)target_lines.insert(f.line);
    auto projection=[&](int source,int target){
      const auto&on=line_points_[static_cast<std::size_t>(target)];
      if(std::find(on.begin(),on.end(),source)!=on.end())return source;
      for(const auto&f:foot_facts_){if(f.source!=source)continue;
        const auto&foot_on=line_points_[static_cast<std::size_t>(target)];
        if(direction_known(f.line,target,0)&&
           std::find(foot_on.begin(),foot_on.end(),f.foot)!=foot_on.end())return f.foot;
      }
      return -1;
    };
    std::size_t midpoint_count=midpoint_facts_.size();
    for(std::size_t index=0;index<midpoint_count;++index){auto mf=midpoint_facts_[index];for(int line:target_lines){
      int a=projection(mf.a,line),m=projection(mf.midpoint,line),b=projection(mf.b,line);
      if(a<0||m<0||b<0||a==b||m==a||m==b)continue;
      std::size_t before=midpoint_facts_.size();
      register_midpoint_fact(m,a,b,"orthogonal projection of midpoint");
      changed|=midpoint_facts_.size()!=before;
    }}
    return changed;
  }

  bool register_circle_chord_reflections(){
    // The perpendicular through a circle's named center bisects every chord
    // and is its reflection axis. Register the reflected endpoint pair so all
    // ordinary reflection, kite, and isosceles-trapezoid rules can reuse it.
    bool changed=false;
    for(std::size_t c=0;c<circle_points_.size();++c){int center=circle_center_ids_[c];
      if(center<0)continue;
      const auto on=circle_points_[c];
      for(std::size_t i=0;i<on.size();++i)for(std::size_t j=i+1;j<on.size();++j){
        int chord=segment(on[i],on[j]);
        for(std::size_t l=0;l<lines_.size();++l){const auto&axis=line_points_[l];
          if(std::find(axis.begin(),axis.end(),center)==axis.end()||
             !direction_known(static_cast<int>(l),chord,1))continue;
          changed|=register_line_reflection(on[j],on[i],static_cast<int>(l),
                                            "circle center perpendicular-to-chord reflection");
        }
      }
    }
    return changed;
  }


  bool register_orthogonal_trapezoids(){
    // Materialize the candidate chords before building direction-component
    // models; segment() can create a previously unseen canonical carrier.
    for(const auto&candidate:circle_cache_)if(candidate.points.size()>=4)
      for(std::size_t i=0;i<candidate.points.size();++i)for(std::size_t j=i+1;j<candidate.points.size();++j)
        segment(candidate.points[i],candidate.points[j]);

    struct Model {PairedLinearSystem along,across;};
    constexpr std::array<std::array<int,4>,3> partitions{{{{0,1,2,3}},{{0,2,1,3}},{{0,3,1,2}}}};
    std::set<int> relevant_roots;
    for(const auto&candidate:circle_cache_)if(candidate.points.size()>=4)
      for(std::size_t i=3;i<candidate.points.size();++i){std::array<int,4> q{candidate.points[0],candidate.points[1],candidate.points[2],candidate.points[i]};
        for(const auto&z:partitions){int first=segment(q[static_cast<std::size_t>(z[0])],q[static_cast<std::size_t>(z[1])]);int second=segment(q[static_cast<std::size_t>(z[2])],q[static_cast<std::size_t>(z[3])]);
          auto [root1,parity1]=direction_find(first);auto [root2,parity2]=direction_find(second);if(root1==root2&&parity1==parity2)relevant_roots.insert(root1);}}
    std::map<int,Model> models;for(int root:relevant_roots)models.try_emplace(root);

    // Midpoint and point-reflection constructions are vector equations in both
    // orthogonal components, independent of the chosen direction component.
    for(auto&[_,model]:models)for(const auto&f:midpoint_facts_){
      model.along.add({{f.midpoint,2},{f.a,-1},{f.b,-1}});
      model.across.add({{f.midpoint,2},{f.a,-1},{f.b,-1}});
    }

    // Points on a carrier share the coordinate perpendicular to its direction.
    for(std::size_t l=0;l<lines_.size();++l){auto [root,parity]=direction_find(static_cast<int>(l));auto found=models.find(root);if(found==models.end())continue;auto&system=parity==0?found->second.across:found->second.along;
      const auto&on=line_points_[l];if(on.size()<2)continue;int representative=on[0];
      for(std::size_t i=1;i<on.size();++i)system.add({{on[i],1},{representative,-1}});
    }

    // Reflection preserves the along-axis component and reverses the normal
    // component about any certified point of the mirror. A perpendicular
    // bisector supplies the same normal equation through its endpoint average
    // even when its midpoint was not explicitly constructed.
    for(const auto&f:line_reflection_facts_){auto [root,parity]=direction_find(f.line);auto found=models.find(root);if(found==models.end())continue;auto&model=found->second;
      auto&along=parity==0?model.along:model.across;auto&normal=parity==0?model.across:model.along;
      along.add({{f.source,1},{f.image,-1}});
      for(int axis_point:line_points_[static_cast<std::size_t>(f.line)])
        normal.add({{f.source,1},{f.image,1},{axis_point,-2}});
      for(const auto&pb:perpendicular_bisectors_)if(pb.line==f.line)
        normal.add({{f.source,1},{f.image,1},{pb.a,-1},{pb.b,-1}});
    }

    auto mirrored=[&](const Model&model,int base_parity,int a,int c,int b,int d){
      const auto&along=base_parity==0?model.along:model.across;
      const auto&normal=base_parity==0?model.across:model.along;
      // Vector AC is the reflection of vector BD across the direction normal
      // to the parallel bases: along components are opposite, normal components
      // agree. This is the non-parallelogram branch of the equal-leg condition.
      return along.proves({{c,1},{a,-1},{d,1},{b,-1}})&&
             normal.proves({{c,1},{a,-1},{d,-1},{b,1}});
    };

    bool changed=false;
    for(const auto&candidate:circle_cache_)if(candidate.points.size()>=4)
      for(std::size_t i=3;i<candidate.points.size();++i){std::array<int,4> q{candidate.points[0],candidate.points[1],candidate.points[2],candidate.points[i]};bool proved=false;
        for(const auto&z:partitions){int a=q[static_cast<std::size_t>(z[0])],b=q[static_cast<std::size_t>(z[1])],c=q[static_cast<std::size_t>(z[2])],d=q[static_cast<std::size_t>(z[3])];
          int first=segment(a,b),second=segment(c,d);auto [root1,parity1]=direction_find(first);auto [root2,parity2]=direction_find(second);
          if(root1!=root2||parity1!=parity2)continue;
          const auto&model=models[root1];
          if(mirrored(model,parity1,a,c,b,d)){changed|=equal_length(a,c,b,d,"orthogonal-component isosceles trapezoid");auto before=cyclic_facts_.size();add_cyclic(a,b,c,d,"orthogonal-component isosceles trapezoid theorem");changed|=cyclic_facts_.size()!=before;proved=true;break;}
          if(mirrored(model,parity1,a,d,b,c)){changed|=equal_length(a,d,b,c,"orthogonal-component isosceles trapezoid");auto before=cyclic_facts_.size();add_cyclic(a,b,c,d,"orthogonal-component isosceles trapezoid theorem");changed|=cyclic_facts_.size()!=before;proved=true;break;}
        }
        if(proved)continue;
      }
    return changed;
  }

  void geometry_closure() {
    // Only declared construction incidences enter the proof layer. Numerical
    // discoveries remain conjectures and therefore cannot prove themselves.
    auto fact_count=[&]{
      std::size_t total=angles_.fact_count()+equal_lengths_.size()+cyclic_facts_.size()+
                        midpoint_facts_.size()+perpendicular_bisectors_.size()+
                        line_reflection_facts_.size();
      for(const auto&on:line_points_)total+=on.size();
      return total;
    };
    bool candidates_detected=false;
    for(;;){std::size_t facts_before=fact_count();
    register_affine_facts();
    normalize_definition_incidences();
    register_center_loci();
    normalize_definition_incidences();
    register_incenter_loci();
    normalize_definition_incidences();
    register_perpendicular_bisector_midpoints();
    register_projection_midpoints();
    register_affine_facts();
    register_formal_affine_facts();
    normalize_definition_incidences();
    register_derived_perpendicular_bisectors();
    // A mirror may acquire additional certified points after its reflection was
    // constructed (notably a circumcenter on a perpendicular-bisector mirror).
    // Propagate reflection symmetry to those late incidences as well.
    for(const auto&f:line_reflection_facts_)for(int axis_point:line_points_[static_cast<std::size_t>(f.line)])
      if(axis_point!=f.source&&axis_point!=f.image){
        equal_length(axis_point,f.source,axis_point,f.image,"late reflection-axis equal distances");
        angles_.add(equation({{segment(axis_point,f.source),1},{segment(axis_point,f.image),1},{f.line,-2}}),0,1,"late line-reflection angle symmetry");
      }
    // A perpendicular bisector is also the reflection axis exchanging its two
    // defining endpoints. Transport both distances and line directions through
    // that mirror. This is the ordinary reflection-isometry rule; it lets the
    // angle chase handle a reflected point without constructing the reflected
    // copies of every other point first.
    for(const auto&f:line_reflection_facts_)for(const auto&pb:perpendicular_bisectors_)
      if(f.line==pb.line){
        if(f.image!=pb.a&&f.source!=pb.b)
          equal_length(f.image,pb.a,f.source,pb.b,"reflection transports distance across perpendicular bisector");
        if(f.image!=pb.b&&f.source!=pb.a)
          equal_length(f.image,pb.b,f.source,pb.a,"reflection transports distance across perpendicular bisector");
        if(f.source!=pb.a&&f.image!=pb.b)
          angles_.add(equation({{segment(f.source,pb.a),1},{segment(f.image,pb.b),1},{f.line,-2}}),0,1,
                      "reflection transports line across perpendicular bisector");
        if(f.source!=pb.b&&f.image!=pb.a)
          angles_.add(equation({{segment(f.source,pb.b),1},{segment(f.image,pb.a),1},{f.line,-2}}),0,1,
                      "reflection transports line across perpendicular bisector");
      }
    // A half-turn about a point on a mirror, composed with reflection in that
    // mirror, sends the source to two images whose connector is parallel to the
    // mirror. This is the elementary composition of two plane isometries.
    for(const auto&reflection:line_reflection_facts_)for(const auto&mf:midpoint_facts_){
      if(mf.midpoint!=reflection.source&&mf.a!=reflection.source&&mf.b!=reflection.source)continue;
      int half_turn_image=-1;
      if(mf.a==reflection.source)half_turn_image=mf.b;
      else if(mf.b==reflection.source)half_turn_image=mf.a;
      else continue;
      const auto&axis=line_points_[static_cast<std::size_t>(reflection.line)];
      if(std::find(axis.begin(),axis.end(),mf.midpoint)==axis.end()||
         half_turn_image==reflection.image)continue;
      parallel_fact(segment(reflection.image,half_turn_image),reflection.line,
                    "line-reflection and half-turn composition");
    }
    // If two segments have the same midpoint, their endpoints form a
    // parallelogram in the crossed order. Both opposite-side parallels are
    // direct affine consequences and use no angle division.
    for(std::size_t i=0;i<midpoint_facts_.size();++i)for(std::size_t j=i+1;j<midpoint_facts_.size();++j){
      const auto&x=midpoint_facts_[i];const auto&y=midpoint_facts_[j];if(x.midpoint!=y.midpoint)continue;
      std::set<int> endpoints{x.a,x.b,y.a,y.b};if(endpoints.size()!=4)continue;
      std::string why="shared-midpoint parallelogram theorem at "+points_[x.midpoint].name;
      parallel_fact(segment(x.a,y.a),segment(x.b,y.b),why);
      parallel_fact(segment(x.a,y.b),segment(x.b,y.a),why);
    }
    // Reflecting an orthocenter across a sideline places its image on both the
    // corresponding altitude and the circumcircle. Register this as a direct
    // construction theorem, rather than forcing the angle lattice to recover
    // the reflection orientation by cancelling a factor of two.
    for(const auto&reflection:line_reflection_facts_)for(const auto&center:orthocenter_facts_)
      if(reflection.source==center.center)for(auto [apex,b,c]:
          {std::array{center.a,center.b,center.c},std::array{center.b,center.a,center.c},std::array{center.c,center.a,center.b}}){
        const auto&axis=line_points_[static_cast<std::size_t>(reflection.line)];
        if(std::find(axis.begin(),axis.end(),b)==axis.end()||std::find(axis.begin(),axis.end(),c)==axis.end())continue;
        inherit_collinearity(reflection.image,apex,center.center,"orthocenter reflection lies on altitude");
        if(apex!=reflection.image&&b!=reflection.image&&c!=reflection.image)
          add_cyclic(apex,b,c,reflection.image,"orthocenter reflection across sideline theorem");
      }
    for (std::size_t c=0;c<circle_points_.size();++c) {
      const auto& on=circle_points_[c];
      int center=circle_center_ids_[c];
      if(center>=0&&!on.empty())for(std::size_t i=1;i<on.size();++i)
        equal_length(center,on[0],center,on[i],"equal radii of circle "+circles_[c].name);
      if(on.size()>=4) for(std::size_t i=3;i<on.size();++i)
        add_cyclic(on[0],on[1],on[2],on[i],"points constructed on circle "+circles_[c].name);
    }
    register_circle_chord_reflections();
    if(!candidates_detected){circle_cache_=detect_circles(true);candidates_detected=true;}
    std::vector<std::pair<std::size_t,int>> pending_midpoint_carriers;
    for(std::size_t i=0;i<midpoint_facts_.size();++i){const auto&f=midpoint_facts_[i];
      for(int p=0;p<(int)points_.size();++p){if(p==f.midpoint||p==f.a||p==f.b)continue;
        if(std::fabs(cross(points_[p],points_[f.a],points_[f.b]))<=EPS*scale(points_[p],points_[f.a],points_[f.b])*10)
          pending_midpoint_carriers.push_back({i,p});}}
    bool changed=true;
    while(changed){changed=false;
      // A midpoint inherits every subsequently proved carrier of its endpoints.
      // Some carriers become available only after theorem closure (for example,
      // reflecting a point across an angle bisector). A numerical test merely
      // prunes candidates; the collinearity premise is proved symbolically.
      std::vector<std::pair<std::size_t,int>> still_pending;still_pending.reserve(pending_midpoint_carriers.size());
      for(auto [index,p]:pending_midpoint_carriers){const auto&f=midpoint_facts_[index];
        if(!proves_collinear({p,f.a,f.b})){still_pending.push_back({index,p});continue;}
        auto target=equation({{segment(p,f.midpoint),1},{segment(p,f.b),-1}});
        if(!angles_.proves(target,0,1))changed|=parallel_fact(segment(p,f.midpoint),segment(p,f.b),"midpoint inherits proved endpoint carrier");
      }
      pending_midpoint_carriers=std::move(still_pending);
      // Angle-defined kite congruence. If AC bisects both endpoint angles of
      // quadrilateral ABCD, then AB=AD and CB=CD. The two hypotheses are
      // checked as undivided equations, so this rule is sound for directed
      // angles modulo 180 degrees and never cancels a factor of two.
      //
      // Candidate diagonals come only from constructed angle bisectors and
      // incenter bisectors. This avoids an O(n^4) scan over all quadrilaterals.
      std::set<std::pair<int,int>> bisector_diagonals(angle_bisector_diagonals_.begin(),angle_bisector_diagonals_.end());
      for(const auto&f:angle_bisector_line_facts_)for(int p:line_points_[static_cast<std::size_t>(f.line)])
        if(p!=f.vertex)bisector_diagonals.insert({f.vertex,p});
      std::vector<std::set<int>> neighbors(points_.size());
      for(const auto&[ends,_]:segment_line_){neighbors[static_cast<std::size_t>(ends.first)].insert(ends.second);neighbors[static_cast<std::size_t>(ends.second)].insert(ends.first);}
      for(auto [a,c]:bisector_diagonals){if(a==c)continue;std::vector<int> common;
        std::set_intersection(neighbors[static_cast<std::size_t>(a)].begin(),neighbors[static_cast<std::size_t>(a)].end(),
          neighbors[static_cast<std::size_t>(c)].begin(),neighbors[static_cast<std::size_t>(c)].end(),std::back_inserter(common));
        for(std::size_t i=0;i<common.size();++i)for(std::size_t j=i+1;j<common.size();++j){int b=common[i],d=common[j];
          if(b==a||b==c||d==a||d==c)continue;
          // Numerical equality is only a candidate filter; both symbolic angle
          // hypotheses below are still required before any fact is registered.
          // This avoids expensive lattice probes for overwhelmingly impossible
          // pairs in large generated configurations.
          if(!near(dist2(points_[a],points_[b]),dist2(points_[a],points_[d]),10)||
             !near(dist2(points_[c],points_[b]),dist2(points_[c],points_[d]),10))continue;
          auto at_a=equation({{segment(a,c),2},{segment(a,b),-1},{segment(a,d),-1}});
          if(!angles_.proves(at_a,0,1))continue;
          auto at_c=equation({{segment(a,c),2},{segment(c,b),-1},{segment(c,d),-1}});
          if(!angles_.proves(at_c,0,1))continue;
          std::string why="angle-defined kite congruence "+points_[a].name+points_[b].name+points_[c].name+points_[d].name;
          changed|=equal_length(a,b,a,d,why);changed|=equal_length(c,b,c,d,why);
        }
      }
      // Every point of a perpendicular bisector is equidistant from its endpoints.
      for(const auto&pb:perpendicular_bisectors_)for(int x:line_points_[pb.line])
        changed|=equal_length(x,pb.a,x,pb.b,"perpendicular bisector distance theorem");

      // In a right triangle, the midpoint of the hypotenuse is equidistant from
      // all three vertices. This is deliberately derived without halving angles.
      for(const auto&mf:midpoint_facts_)for(int b=0;b<(int)points_.size();++b)
        if(b!=mf.a&&b!=mf.b&&b!=mf.midpoint){
          if(angles_.proves(equation({{segment(mf.a,b),1},{segment(b,mf.b),-1}}),1,2)){
            changed|=equal_length(mf.midpoint,mf.a,mf.midpoint,b,"right triangle midpoint theorem");
            changed|=equal_length(mf.midpoint,mf.a,mf.midpoint,mf.b,"right triangle midpoint theorem");
          }
        }

      // Targeted nine-point-circle theorem. The component facts (three side
      // midpoints and three altitude feet) are all construction-certified.
      auto midpoint_of=[&](int a,int b){
        for(const auto&f:midpoint_facts_)
          if((f.a==a&&f.b==b)||(f.a==b&&f.b==a))return f.midpoint;
        return -1;
      };
      auto foot_on=[&](int source,int a,int b){for(const auto&f:foot_facts_)if(f.source==source){
        const auto&on=line_points_[f.line];if(std::find(on.begin(),on.end(),a)!=on.end()&&
          std::find(on.begin(),on.end(),b)!=on.end())return f.foot;}return -1;};
      for(const auto&t:initial_triangles_){int a=t[0],b=t[1],c=t[2];
        std::vector<int> nine{midpoint_of(a,b),midpoint_of(b,c),midpoint_of(c,a),
          foot_on(a,b,c),foot_on(b,c,a),foot_on(c,a,b)};
        if(std::find(nine.begin(),nine.end(),-1)==nine.end())for(std::size_t i=3;i<nine.size();++i){
          auto before=cyclic_facts_.size();add_cyclic(nine[0],nine[1],nine[2],nine[i],"nine-point circle theorem");changed|=cyclic_facts_.size()!=before;}
      }

      // Equal-radius components immediately supply cyclic point sets.
      using SegmentKey=std::pair<int,int>;
      std::map<SegmentKey,int> node;std::vector<int> parent;
      auto get_node=[&](SegmentKey s){auto [it,inserted]=node.emplace(s,(int)parent.size());if(inserted)parent.push_back(it->second);return it->second;};
      auto find_root=[&](int x){while(parent[x]!=x){parent[x]=parent[parent[x]];x=parent[x];}return x;};
      auto unite=[&](int a,int b){a=find_root(a);b=find_root(b);if(a!=b)parent[b]=a;};
      for(const auto&f:equal_lengths_)unite(get_node(f.first),get_node(f.second));
      // Every triangle center fixed by the symmetry of an isosceles triangle
      // lies on that symmetry axis. For an incenter, this first registers the
      // equal distances to the base endpoints; the ordinary kite closure then
      // derives that the apex-incenter line is perpendicular to the base.
      auto lengths_same=[&](int a,int b,int c,int d){auto x=node.find(lenkey(a,b)),y=node.find(lenkey(c,d));return x!=node.end()&&y!=node.end()&&find_root(x->second)==find_root(y->second);};
      for(const auto&f:incenter_facts_)for(auto [apex,u,v]:
          {std::array{f.a,f.b,f.c},std::array{f.b,f.a,f.c},std::array{f.c,f.a,f.b}})
        if(lengths_same(apex,u,apex,v)){changed|=equal_length(f.center,u,f.center,v,"incenter symmetry in isosceles triangle");
          int axis=segment(apex,f.center),base=segment(u,v);if(!direction_known(axis,base,1))changed|=perpendicular_fact(axis,base,"isosceles-triangle incenter symmetry axis");}
      for(int o=0;o<(int)points_.size();++o){std::map<int,std::vector<int>> groups;
        for(int p=0;p<(int)points_.size();++p)if(p!=o){auto it=node.find(lenkey(o,p));if(it!=node.end())groups[find_root(it->second)].push_back(p);}
        for(const auto&[_,on]:groups)if(on.size()>=4)for(std::size_t i=3;i<on.size();++i){auto before=cyclic_facts_.size();add_cyclic(on[0],on[1],on[2],on[i],"equal radii about "+points_[o].name);changed|=cyclic_facts_.size()!=before;}
      }

      // Three equal radii to a known cyclic quadruple identify its circle
      // center, so the radius to the fourth point is equal as well. Try every
      // omitted vertex to cover all orderings of the cyclic fact.
      for(const auto&q:cyclic_facts_)for(int o=0;o<(int)points_.size();++o){
        if(std::find(q.begin(),q.end(),o)!=q.end())continue;
        for(int omitted=0;omitted<4;++omitted){int representative=-1,root=-1;bool known=true;
          for(int i=0;i<4;++i)if(i!=omitted){auto it=node.find(lenkey(o,q[static_cast<std::size_t>(i)]));if(it==node.end()){known=false;break;}int current=find_root(it->second);if(root<0){root=current;representative=q[static_cast<std::size_t>(i)];}else if(root!=current){known=false;break;}}
          if(known)changed|=equal_length(o,representative,o,q[static_cast<std::size_t>(omitted)],"cyclic circle-center radius completion");
        }
      }

      // A detected numerical circle can contain both elementary and genuinely
      // difficult points. Certify its members independently: requiring one
      // proof for the entire maximal set lets a single hard member suppress
      // useful cyclic facts for every easy member on that circle.
      for(const auto&x:circle_cache_)if(x.points.size()>=4)
        for(std::size_t i=3;i<x.points.size();++i){
          std::vector<int> quadruple{x.points[0],x.points[1],x.points[2],x.points[i]};
          std::set<int>w;if(proves_cyclic(quadruple,&w)){
            auto before=cyclic_facts_.size();
            add_cyclic(x.points[0],x.points[1],x.points[2],x.points[i],
                       "converse cyclic angle theorem "+point_list(quadruple));
            changed|=cyclic_facts_.size()!=before;
          }
        }
      // Kite: two equidistant vertices give both perpendicular diagonals and the
      // full reflection angle relation between the four corresponding sides.
      std::map<int,std::vector<SegmentKey>> length_classes;
      for(const auto&[segment_key,id]:node)length_classes[find_root(id)].push_back(segment_key);
      std::map<std::pair<int,int>,std::set<int>> wings;
      for(const auto&[_,segments]:length_classes)for(std::size_t i=0;i<segments.size();++i)for(std::size_t j=i+1;j<segments.size();++j){
        auto x=segments[i],y=segments[j];int vertex=-1,b=-1,c=-1;
        if(x.first==y.first){vertex=x.first;b=x.second;c=y.second;}
        else if(x.first==y.second){vertex=x.first;b=x.second;c=y.first;}
        else if(x.second==y.first){vertex=x.second;b=x.first;c=y.second;}
        else if(x.second==y.second){vertex=x.second;b=x.first;c=y.first;}
        if(vertex<0||b==c)continue;
        if(b>c)std::swap(b,c);
        wings[{b,c}].insert(vertex);
        auto isosceles=equation({{segment(vertex,b),1},{segment(vertex,c),1},{segment(b,c),-2}});
        if(!angles_.proves(isosceles,0,1))changed|=angles_.add(isosceles,0,1,"isosceles triangle theorem "+points_[vertex].name+points_[b].name+points_[c].name);
      }
      // Converse perpendicular-bisector theorem. Reuse an existing canonical
      // locus instead of creating a new line for every isolated equality.
      for(const auto&[base,vertices]:wings){auto locus=perpendicular_bisector_loci_.find(base);if(locus==perpendicular_bisector_loci_.end())continue;
        for(int vertex:vertices){auto&on=line_points_[static_cast<std::size_t>(locus->second)];if(std::find(on.begin(),on.end(),vertex)==on.end()){incidence(vertex,locus->second,"equal distances imply perpendicular-bisector incidence");changed=true;}}}
      for(const auto& [base,vertices]:wings)for(auto ai=vertices.begin();ai!=vertices.end();++ai)for(auto di=std::next(ai);di!=vertices.end();++di){
        std::string why="kite theorem "+points_[*ai].name+points_[base.first].name+points_[*di].name+points_[base.second].name;
        int ad=segment(*ai,*di),bc=segment(base.first,base.second);if(!angles_.proves(equation({{ad,1},{bc,-1}}),1,2))changed|=perpendicular_fact(ad,bc,why);
        auto symmetry=equation({{segment(*ai,base.first),1},{segment(*ai,base.second),1},{segment(*di,base.first),-1},{segment(*di,base.second),-1}});
        if(!angles_.proves(symmetry,0,1))changed|=angles_.add(symmetry,0,1,why+" [reflection angles]");
      }

      changed|=register_orthocenter_closure();

    }
    // Component models are comparatively expensive and depend on the completed
    // incidence/direction closure. Build them once, after the ordinary fixed
    // point, rather than once per theorem round.
    register_orthogonal_trapezoids();
    // Theorems found during closure can identify previously separate carriers
    // (for example, two proved-parallel lines through the same point). Merge
    // them once more so definition-level incidence queries use the completed
    // fact base rather than only the pre-closure construction graph.
    normalize_definition_incidences();
    if(fact_count()==facts_before)break;
    }
  }

  std::string automatic_name(const std::string&prefix){
    return automatic_namespace_+prefix+"$"+std::to_string(automatic_serial_++);
  }
  template<class Depths>
  int depth_weighted_pick(const std::vector<int>&items,const Depths&depths){
    if(items.empty())return -1;
    int minimum=std::numeric_limits<int>::max();
    for(int id:items)minimum=std::min(minimum,depths[static_cast<std::size_t>(id)]);
    std::vector<double> weights;weights.reserve(items.size());
    for(int id:items){int relative=depths[static_cast<std::size_t>(id)]-minimum;
      weights.push_back(relative>500?0.0:std::ldexp(1.0,-2*relative));}
    std::discrete_distribution<std::size_t> choose(weights.begin(),weights.end());
    return items[choose(generation_rng_)];
  }
  int random_point(){
    std::vector<int> ids;
    if(symmetry_enabled_)ids.assign(symmetry_primary_points_.begin(),symmetry_primary_points_.end());
    else {ids.resize(points_.size());std::iota(ids.begin(),ids.end(),0);}
    return depth_weighted_pick(ids,point_depth_);
  }
  std::vector<int> random_distinct_points(std::size_t count,bool noncollinear=false){
    std::vector<int> primary;
    if(symmetry_enabled_)primary.assign(symmetry_primary_points_.begin(),symmetry_primary_points_.end());
    else {primary.resize(points_.size());std::iota(primary.begin(),primary.end(),0);}
    if(primary.size()<count)return {};
    for(int attempt=0;attempt<64;++attempt){
      std::vector<int> available=primary,out;
      while(out.size()<count){int p=depth_weighted_pick(available,point_depth_);out.push_back(p);available.erase(std::find(available.begin(),available.end(),p));}
      if(!noncollinear||count<3||std::fabs(cross(points_[out[0]],points_[out[1]],points_[out[2]]))>
          EPS*scale(points_[out[0]],points_[out[1]],points_[out[2]]))return out;
    }
    return {};
  }
  std::vector<int> named_lines()const{
    // A public line(P,Q) normally canonicalizes to the already existing
    // internal segment carrier @PQ.  Eligibility belongs to the public alias,
    // not to the canonical carrier's original "segment" label.
    std::vector<bool> exposed(lines_.size(),false);
    for(const auto&[name,id]:line_id_)if(!name.empty()&&name[0]!='@')
      exposed[static_cast<std::size_t>(id)]=true;
    std::vector<int> out;for(std::size_t i=0;i<lines_.size();++i)
      if(exposed[i]&&(!symmetry_enabled_||symmetry_primary_lines_.count(static_cast<int>(i))))
        out.push_back(static_cast<int>(i));
    return out;
  }
  std::string public_line_name(int id)const{
    std::string best;
    for(const auto&[name,candidate]:line_id_)if(candidate==id&&!name.empty()&&name[0]!='@'&&
        (best.empty()||name<best))best=name;
    return best.empty()?lines_[static_cast<std::size_t>(id)].name:best;
  }
  std::vector<int> named_circles()const{
    std::vector<int> out;for(std::size_t i=0;i<circles_.size();++i)
      if(!circles_[i].name.empty()&&circles_[i].name[0]!='@'&&
         (!symmetry_enabled_||symmetry_primary_circles_.count(static_cast<int>(i))))out.push_back(static_cast<int>(i));
    return out;
  }
  int create_random_line(){
    if(points_.size()<2)return -1;
    for(int attempt=0;attempt<64;++attempt){
      int kind=std::uniform_int_distribution<int>(0,4)(generation_rng_);
      const char*selected=kind==0?"line":kind==1?"perp_bisector":kind==2?"parallel":kind==3?"perpendicular":"angle_bisector";
      if(!construction_enabled(selected))continue;
      auto pair=random_distinct_points(2);if(pair.empty())return -1;
      std::string name;int depth=1+std::max(point_depth_[pair[0]],point_depth_[pair[1]]);
      if(kind==0){name=automatic_name("L");execute({"line",name,points_[pair[0]].name,points_[pair[1]].name},0);}
      else if(kind==1){name=automatic_name("PB");execute({"perp_bisector",name,points_[pair[0]].name,points_[pair[1]].name},0);}
      else if(kind==2||kind==3){
        auto pool=named_lines();if(pool.empty()){if(!construction_enabled("line"))continue;kind=0;name=automatic_name("L");execute({"line",name,points_[pair[0]].name,points_[pair[1]].name},0);}
        else {int base=depth_weighted_pick(pool,line_depth_);int p=random_point();depth=1+std::max(point_depth_[p],line_depth_[base]);
          name=automatic_name(kind==2?"Par":"Perp");execute({kind==2?"parallel":"perpendicular",name,points_[p].name,public_line_name(base)},0);}
      } else {
        auto triple=random_distinct_points(3,true);if(triple.empty())continue;
        depth=1+std::max({point_depth_[triple[0]],point_depth_[triple[1]],point_depth_[triple[2]]});
        name=automatic_name("Bis");execute({"angle_bisector",name,points_[triple[0]].name,points_[triple[1]].name,points_[triple[2]].name},0);
      }
      int id=lid(name);
      if(lines_[static_cast<std::size_t>(id)].name==name)line_depth_[static_cast<std::size_t>(id)]=depth;
      else line_depth_[static_cast<std::size_t>(id)]=std::min(line_depth_[static_cast<std::size_t>(id)],depth);
      if(symmetry_enabled_)symmetry_primary_lines_.insert(id);
      return id;
    }
    return -1;
  }
  int create_random_circle(){
    if(points_.size()<2)return -1;
    for(int attempt=0;attempt<64;++attempt){
      // Circumcircles are the standard Olympiad circle construction, so give
      // them most of the automatic-circle probability. Center--point circles
      // and incircles remain available for configurations that need them.
      int kind=std::uniform_int_distribution<int>(0,5)(generation_rng_);
      std::string name;int depth=0;
      if(kind<=3){
        if(!construction_enabled("circumcircle"))continue;
        auto triple=random_distinct_points(3,true);if(triple.empty())continue;
        depth=1+std::max({point_depth_[triple[0]],point_depth_[triple[1]],point_depth_[triple[2]]});
        name=automatic_name("Circ");execute({"circumcircle",name,points_[triple[0]].name,points_[triple[1]].name,points_[triple[2]].name},0);
      } else if(kind==4||incenter_facts_.empty()) {
        if(!construction_enabled("circle"))continue;
        auto pair=random_distinct_points(2);int o=pair[0],p=pair[1];depth=1+std::max(point_depth_[o],point_depth_[p]);
        name=automatic_name("Circle");execute({"circle",name,points_[o].name,points_[p].name},0);
      } else {
        if(!construction_enabled("incircle"))continue;
        const auto&f=incenter_facts_[std::uniform_int_distribution<std::size_t>(0,incenter_facts_.size()-1)(generation_rng_)];
        depth=1+std::max({point_depth_[f.center],point_depth_[f.a],point_depth_[f.b],point_depth_[f.c]});
        name=automatic_name("Inc");execute({"incircle",name,points_[f.center].name,points_[f.a].name,points_[f.b].name,points_[f.c].name},0);
      }
      int id=cid(name);circle_depth_[static_cast<std::size_t>(id)]=depth;
      if(symmetry_enabled_)symmetry_primary_circles_.insert(id);
      return id;
    }
    return -1;
  }
  bool create_random_point(){
    if(points_.size()<2)return false;
    int kind=std::uniform_int_distribution<int>(0,9)(generation_rng_);
    const char*selected=kind==0?"midpoint":kind==1?"reflection_point":kind==2?"circumcenter":
      kind==3?"orthocenter":kind==4?"incenter":kind==5?"reflection_line":kind==6?"foot":
      kind==7?"intersection_ll":kind==8?"intersection_lc_known":"intersection_cc_known";
    if(!construction_enabled(selected))return false;
    std::size_t before=points_.size();std::string name;
    auto finish=[&](int depth){if(points_.size()>before)point_depth_[static_cast<std::size_t>(pid(name))]=depth;return points_.size()>before;};
    if(kind==0||kind==1){
      auto pair=random_distinct_points(2);if(pair.empty())return false;int depth=1+std::max(point_depth_[pair[0]],point_depth_[pair[1]]);
      name=automatic_name(kind==0?"M":"Rp");execute({kind==0?"midpoint":"reflection_point",name,points_[pair[0]].name,points_[pair[1]].name},0);return finish(depth);
    }
    if(kind>=2&&kind<=4){
      auto triple=random_distinct_points(3,true);if(triple.empty())return false;
      int depth=1+std::max({point_depth_[triple[0]],point_depth_[triple[1]],point_depth_[triple[2]]});
      const char*op=kind==2?"circumcenter":kind==3?"orthocenter":"incenter";
      name=automatic_name(kind==2?"O":kind==3?"H":"I");execute({op,name,points_[triple[0]].name,points_[triple[1]].name,points_[triple[2]].name},0);return finish(depth);
    }
    if(kind==5||kind==6){
      auto pool=named_lines();if(pool.empty())create_random_line();pool=named_lines();if(pool.empty())return false;
      int p=random_point(),line=depth_weighted_pick(pool,line_depth_);int depth=1+std::max(point_depth_[p],line_depth_[line]);
      name=automatic_name(kind==5?"Rl":"Foot");execute({kind==5?"reflection_line":"foot",name,points_[p].name,public_line_name(line)},0);return finish(depth);
    }
    if(kind==7){
      if(named_lines().size()<2){create_random_line();create_random_line();}
      for(int attempt=0;attempt<32;++attempt){auto pool=named_lines();if(pool.size()<2)return false;
        int a=depth_weighted_pick(pool,line_depth_),b=depth_weighted_pick(pool,line_depth_);if(a==b)continue;
        if(std::fabs(lines_[a].a*lines_[b].b-lines_[b].a*lines_[a].b)<=EPS)continue;
        int depth=1+std::max(line_depth_[a],line_depth_[b]);name=automatic_name("Xll");
        execute({"intersection_ll",name,public_line_name(a),public_line_name(b)},0);return finish(depth);
      }return false;
    }
    if(kind==8){
      auto line_pool=named_lines();if(line_pool.empty()){create_random_line();line_pool=named_lines();}if(line_pool.empty())return false;
      for(int attempt=0;attempt<32;++attempt){int line=depth_weighted_pick(line_pool,line_depth_),circle=-1,k=-1;auto circle_pool=named_circles();
        // Prefer an existing circle--especially one of the automatically
        // generated circumcircles--when it shares a certified point with the
        // selected line.
        if(!circle_pool.empty()){int candidate=depth_weighted_pick(circle_pool,circle_depth_);std::vector<int> common;
          for(int p:line_points_[static_cast<std::size_t>(line)])if(std::find(circle_points_[static_cast<std::size_t>(candidate)].begin(),circle_points_[static_cast<std::size_t>(candidate)].end(),p)!=circle_points_[static_cast<std::size_t>(candidate)].end())common.push_back(p);
          if(!common.empty()){circle=candidate;k=common[std::uniform_int_distribution<std::size_t>(0,common.size()-1)(generation_rng_)];}}
        // If no existing pair has a known root, create a circumcircle through a
        // point already on the line and two additional low-depth points.
        if(circle<0){if(!construction_enabled("circumcircle"))return false;const auto&on=line_points_[static_cast<std::size_t>(line)];if(on.empty())continue;k=on[std::uniform_int_distribution<std::size_t>(0,on.size()-1)(generation_rng_)];
          auto pair=random_distinct_points(2);if(pair.empty()||pair[0]==k||pair[1]==k||std::fabs(cross(points_[k],points_[pair[0]],points_[pair[1]]))<=EPS*scale(points_[k],points_[pair[0]],points_[pair[1]]))continue;
          std::string circle_name=automatic_name("Circ");int circle_depth=1+std::max({point_depth_[k],point_depth_[pair[0]],point_depth_[pair[1]]});
          execute({"circumcircle",circle_name,points_[k].name,points_[pair[0]].name,points_[pair[1]].name},0);circle=cid(circle_name);circle_depth_[static_cast<std::size_t>(circle)]=circle_depth;}
        const auto&ln=lines_[line];const auto&cc=circles_[circle];long double dx=ln.b,dy=-ln.a;
        if(std::fabs(-2*dot(points_[k].x-cc.center.x,points_[k].y-cc.center.y,dx,dy))<=EPS)continue;
        int depth=1+std::max(line_depth_[line],circle_depth_[circle]);name=automatic_name("Xlc");
        execute({"intersection_lc_known",name,public_line_name(line),circles_[circle].name,points_[k].name},0);return finish(depth);
      }return false;
    }
    for(int attempt=0;attempt<32;++attempt){
      int c1=-1,c2=-1,k=-1;auto pool=named_circles();
      if(pool.size()>=2){c1=depth_weighted_pick(pool,circle_depth_);c2=depth_weighted_pick(pool,circle_depth_);if(c1==c2)continue;std::vector<int> common;
        for(int p:circle_points_[static_cast<std::size_t>(c1)])if(std::find(circle_points_[static_cast<std::size_t>(c2)].begin(),circle_points_[static_cast<std::size_t>(c2)].end(),p)!=circle_points_[static_cast<std::size_t>(c2)].end())common.push_back(p);
        if(!common.empty())k=common[std::uniform_int_distribution<std::size_t>(0,common.size()-1)(generation_rng_)];}
      if(k<0){if(!construction_enabled("circumcircle"))return false;k=random_point();auto first=random_distinct_points(2),second=random_distinct_points(2);if(first.empty()||second.empty())continue;
        if(first[0]==k||first[1]==k||second[0]==k||second[1]==k||std::fabs(cross(points_[k],points_[first[0]],points_[first[1]]))<=EPS*scale(points_[k],points_[first[0]],points_[first[1]])||std::fabs(cross(points_[k],points_[second[0]],points_[second[1]]))<=EPS*scale(points_[k],points_[second[0]],points_[second[1]]))continue;
        std::string c1_name=automatic_name("Circ"),c2_name=automatic_name("Circ");
        execute({"circumcircle",c1_name,points_[k].name,points_[first[0]].name,points_[first[1]].name},0);execute({"circumcircle",c2_name,points_[k].name,points_[second[0]].name,points_[second[1]].name},0);
        c1=cid(c1_name);c2=cid(c2_name);circle_depth_[static_cast<std::size_t>(c1)]=1+std::max({point_depth_[k],point_depth_[first[0]],point_depth_[first[1]]});circle_depth_[static_cast<std::size_t>(c2)]=1+std::max({point_depth_[k],point_depth_[second[0]],point_depth_[second[1]]});}
      if(c1<0||c2<0||k<0)continue;
      Point center1=circles_[static_cast<std::size_t>(c1)].center;
      Point center2=circles_[static_cast<std::size_t>(c2)].center;
      if(dist2(center1,center2)<=EPS*EPS)continue;
      Line axis=through("",center1,center2,"");
      if(std::fabs(axis.a*points_[k].x+axis.b*points_[k].y+axis.c)<=EPS)continue;
      int depth=1+std::max(circle_depth_[c1],circle_depth_[c2]);name=automatic_name("Xcc");
      execute({"intersection_cc_known",name,circles_[c1].name,circles_[c2].name,points_[k].name},0);return finish(depth);
    }
    return false;
  }
  bool create_self_dual_point(){
    // When exactly one cap slot remains, an ordinary construction/dual pair
    // cannot fit. The midpoint of a swapped point pair is fixed by the swap,
    // so its dual construction aliases the same new point and consumes one
    // slot while preserving construction-level symmetry.
    std::vector<std::tuple<int,int,int>> pairs;
    for(std::size_t p=0;p<points_.size();++p){
      ObjectRef object{ObjectKind::point,points_[p].name};auto found=symmetry_dual_.find(object);
      if(found==symmetry_dual_.end())continue;
      int q=pid(found->second.name);
      if(static_cast<int>(p)>=q)continue;
      pairs.push_back({std::max(point_depth_[p],point_depth_[static_cast<std::size_t>(q)]),static_cast<int>(p),q});
    }
    std::sort(pairs.begin(),pairs.end());
    for(const auto&[input_depth,a,b]:pairs){
      std::size_t before=points_.size(),begin=construction_commands_.size();
      std::string name=automatic_name("SymM");execute({"midpoint",name,points_[a].name,points_[b].name},0);
      std::size_t end=construction_commands_.size();mark_primary_commands(begin,end);replay_symmetric_commands(begin,end);
      if(points_.size()>before){point_depth_[static_cast<std::size_t>(pid(name))]=1+input_depth;return true;}
    }
    return false;
  }
  void expand_points(){
    establish_symmetry();
    auto room=[&]{return !max_points_||points_.size()<max_points_;};
    // Without a point cap, retain the finite known-root scan for explicitly
    // supplied lines and circles. With a cap, known-root intersections join the
    // randomized construction mix below.
    if(construction_enabled("intersection_lc_known")&&!max_points_){
      std::map<int,std::string> scan_lines;
      for(const auto&[name,id]:line_id_)if(!name.empty()&&name[0]!='@'&&
          (!symmetry_enabled_||symmetry_primary_lines_.count(id)))scan_lines.emplace(id,name);
      std::size_t circle_count=circles_.size();
      for(const auto&[line_id,line_name]:scan_lines){std::size_t l=static_cast<std::size_t>(line_id);
        for(std::size_t c=0;c<circle_count;++c){if(symmetry_enabled_&&!symmetry_primary_circles_.count(static_cast<int>(c)))continue;auto on_line=line_points_[l],on_circle=circle_points_[c];
          for(int k:on_line)if(std::find(on_circle.begin(),on_circle.end(),k)!=on_circle.end()){
            const auto&ln=lines_[l];const auto&cc=circles_[c];long double dx=ln.b,dy=-ln.a;
            if(std::fabs(-2*dot(points_[k].x-cc.center.x,points_[k].y-cc.center.y,dx,dy))<=EPS)continue;
            std::string name="X("+line_name+","+circles_[c].name+","+points_[k].name+")";
            if(!point_id_.count(name)){
              std::size_t begin=construction_commands_.size();
              execute({"intersection_lc_known",name,line_name,circles_[c].name,points_[k].name},0);
              std::size_t end=construction_commands_.size();mark_primary_commands(begin,end);if(symmetry_enabled_)replay_symmetric_commands(begin,end);
            }
          }}
      }return;
    }
    if(!max_points_||!room())return;
    std::size_t failed=0,max_failed=std::max<std::size_t>(2000,100*max_points_);
    while(room()&&failed<max_failed){
      if(symmetry_enabled_&&max_points_-points_.size()<2){if(create_self_dual_point())continue;break;}
      std::size_t batch_begin=construction_commands_.size(),point_count=points_.size();
      // Auxiliary constructions are deliberately interleaved. They make line
      // intersections and known-root circle intersections available without
      // letting the supporting-object count grow independently of the point cap.
      if(named_lines().size()<2||std::uniform_int_distribution<int>(0,3)(generation_rng_)==0)create_random_line();
      if(named_circles().empty()||std::uniform_int_distribution<int>(0,7)(generation_rng_)==0)create_random_circle();
      bool created=create_random_point();std::size_t batch_end=construction_commands_.size();
      mark_primary_commands(batch_begin,batch_end);
      if(symmetry_enabled_)replay_symmetric_commands(batch_begin,batch_end);
      if(created||points_.size()>point_count)failed=0;else ++failed;
    }
    if(room()){
      if(symmetry_enabled_&&max_points_-points_.size()<2)
        std::cerr<<"warning: symmetry pairing stopped at "<<points_.size()<<" of "<<max_points_<<" points\n";
      else std::cerr<<"warning: random construction search exhausted at "<<points_.size()<<" of "<<max_points_<<" points\n";
    }
  }

  static long long quant(long double x, long double step=1e-8L) { return std::llround(x/step); }
  std::vector<Candidate> detect_lines() const {
    std::map<std::vector<int>,Candidate> uniq;int n=(int)points_.size();
    for(int a=0;a<n;++a){std::map<long long,std::vector<int>> groups;
      for(int b=0;b<n;++b)if(a!=b){long double ang=std::atan2(points_[b].y-points_[a].y,points_[b].x-points_[a].x);while(ang<0)ang+=PI;while(ang>=PI)ang-=PI;groups[quant(ang)].push_back(b);}
      for(auto&[_,g]:groups)if(g.size()>=2){g.push_back(a);std::sort(g.begin(),g.end());g.erase(std::unique(g.begin(),g.end()),g.end());bool ok=true;for(int x:g)if(std::fabs(cross(points_[g[0]],points_[g[1]],points_[x]))>EPS*scale(points_[g[0]],points_[g[1]],points_[x])*10)ok=false;if(ok)uniq[g]={"collinear",g,"direction hash"};}
    }std::vector<Candidate> out;for(auto&[_,v]:uniq)out.push_back(v);return out;
  }
  std::vector<Candidate> detect_circles(bool respect_budget=true) const {
    std::map<std::vector<int>,Candidate> uniq;int n=(int)points_.size();
    // Every declared circle is cheap to scan.
    for(auto&c:circles_){std::vector<int> on;for(int i=0;i<n;++i)if(near(dist2(c.center,points_[i]),c.r2,10))on.push_back(i);if(on.size()>=4)uniq[on]={"concyclic",on,"declared circle "+c.name};}
    auto dual_of=[&](int p)->int{auto found=symmetry_dual_.find({ObjectKind::point,points_[static_cast<std::size_t>(p)].name});
      return found==symmetry_dual_.end()?-1:pid(found->second.name);};
    std::uint64_t triples=(std::uint64_t)n*(n-1)*(n-2)/6;
    if(symmetric_only_){std::uint64_t pairs=0,self=0;for(int p=0;p<n;++p){int q=dual_of(p);pairs+=q>p;self+=q==p;}
      triples=pairs*static_cast<std::uint64_t>(std::max(0,n-2))+(self<3?0:self*(self-1)*(self-2)/6);}
    if(respect_budget&&triples>circle_budget_){std::cerr<<"warning: general circle scan skipped ("<<triples<<" triples > circle_budget)\n";std::vector<Candidate>out;for(auto&[_,v]:uniq)out.push_back(v);return out;}
    // O(n^3) time and O(n^2) peak memory. Fixing one anchor retains enough
    // information to discover every circle while allowing each hash table to be freed.
    for(int a=0;a<n;++a){
      std::map<std::tuple<long long,long long,long long>,std::set<int>> bins;
      auto add_triple=[&](int b,int c){
        if(std::fabs(cross(points_[a],points_[b],points_[c]))<=EPS*scale(points_[a],points_[b],points_[c]))return;
        Point o=circumcenter(points_[a],points_[b],points_[c],"","scan");long double r2=dist2(o,points_[a]);
        auto key=std::make_tuple(quant(o.x),quant(o.y),quant(r2));auto&s=bins[key];s.insert(a);s.insert(b);s.insert(c);
      };
      if(symmetric_only_){int dual=dual_of(a);
        if(dual>a)for(int c=0;c<n;++c)if(c!=a&&c!=dual)add_triple(dual,c);
        if(dual==a)for(int b=a+1;b<n;++b)if(dual_of(b)==b)
          for(int c=b+1;c<n;++c)if(dual_of(c)==c)add_triple(b,c);
      } else for(int b=0;b<n;++b)if(b!=a)for(int c=b+1;c<n;++c)if(c!=a)add_triple(b,c);
      for(auto&[_,s]:bins)if(s.size()>=4){std::vector<int>v(s.begin(),s.end());Point o=circumcenter(points_[v[0]],points_[v[1]],points_[v[2]],"","verify");bool ok=true;for(int x:v)if(!near(dist2(o,points_[x]),dist2(o,points_[v[0]]),10))ok=false;if(ok)uniq[v]={"concyclic",v,"circle hash"};}
    }
    std::vector<Candidate>out;for(auto&[_,v]:uniq)out.push_back(v);return out;
  }

  using ProducerMap=std::map<ObjectRef,std::size_t>;
  std::string construction_expression(const ObjectRef&object,const ProducerMap&producer,
      std::map<ObjectRef,std::string>&memo,std::set<ObjectRef>&active,bool force=false)const{
    ObjectRef canonical=object;
    if(canonical.kind==ObjectKind::line)if(auto alias=line_canonical_name_.find(canonical.name);alias!=line_canonical_name_.end())canonical.name=alias->second;
    // Points supplied by the input configuration are named lemmas from the
    // user's perspective. Expand their own POINT assignment, but keep them
    // atomic whenever a later generated definition refers to them.
    if(!force&&canonical.kind==ObjectKind::point&&input_point_names_.count(canonical.name))return canonical.name;
    if(auto known=memo.find(canonical);known!=memo.end())return known->second;
    if(!active.insert(canonical).second)return canonical.name;
    auto done=[&](std::string value){active.erase(canonical);memo[canonical]=value;return value;};
    auto found=producer.find(canonical);if(found==producer.end()){
      if(canonical.kind==ObjectKind::line){int id=lid(canonical.name);const auto&on=line_points_[static_cast<std::size_t>(id)];if(on.size()>=2)
        return done("line("+construction_expression({ObjectKind::point,points_[on[0]].name},producer,memo,active)+","+construction_expression({ObjectKind::point,points_[on[1]].name},producer,memo,active)+")");}
      return done(canonical.name);
    }
    const auto&t=construction_commands_[found->second];const auto&op=t[0];
    auto point=[&](std::size_t i){return construction_expression({ObjectKind::point,t[i]},producer,memo,active);};
    auto line=[&](std::size_t i){return construction_expression({ObjectKind::line,t[i]},producer,memo,active);};
    auto circle=[&](std::size_t i){return construction_expression({ObjectKind::circle,t[i]},producer,memo,active);};
    if(op=="triangle"||op=="quadrilateral"||op=="cyclic_quad")return done(object.name);
    if(op=="point")return done("point("+t[2]+","+t[3]+")");
    if(op=="line")return done("line("+point(2)+","+point(3)+")");
    if(op=="midpoint")return done("midpoint("+point(2)+","+point(3)+")");
    if(op=="perp_bisector")return done("perpendicular_bisector("+point(2)+","+point(3)+")");
    if(op=="parallel")return done("parallel("+point(2)+","+line(3)+")");
    if(op=="perpendicular")return done("perpendicular("+point(2)+","+line(3)+")");
    if(op=="angle_bisector")return done("angle_bisector("+point(2)+","+point(3)+","+point(4)+")");
    if(op=="reflection_line")return done("reflect("+point(2)+","+line(3)+")");
    if(op=="reflection_point")return done("reflect("+point(2)+","+point(3)+")");
    if(op=="foot")return done("foot("+point(2)+","+line(3)+")");
    if(op=="intersection_ll")return done("intersect("+line(2)+","+line(3)+")");
    if(op=="circumcenter"||op=="orthocenter"||op=="incenter")
      return done(op+"("+point(2)+","+point(3)+","+point(4)+")");
    if(op=="circle")return done("circle("+point(2)+","+point(3)+")");
    if(op=="circumcircle")return done("circumcircle("+point(2)+","+point(3)+","+point(4)+")");
    if(op=="incircle")return done("incircle("+point(2)+","+point(3)+","+point(4)+","+point(5)+")");
    if(op=="intersection_lc_known")return done("other_intersection("+line(2)+","+circle(3)+","+point(4)+")");
    if(op=="intersection_cc_known")return done("other_intersection("+circle(2)+","+circle(3)+","+point(4)+")");
    return done(object.name);
  }
  std::vector<std::string> point_definitions()const{
    ProducerMap producer;
    for(std::size_t i=0;i<construction_commands_.size();++i)
      for(auto output:command_outputs(construction_commands_[i]))producer.emplace(std::move(output),i);
    std::map<ObjectRef,std::string> memo;std::set<ObjectRef> active;std::vector<std::string> out;out.reserve(points_.size());
    for(const auto&p:points_){
      std::string expression;
      if(p.origin.rfind("random initial",0)==0)expression="initial("+p.name+")";
      else expression=construction_expression({ObjectKind::point,p.name},producer,memo,active,true);
      out.push_back("POINT "+p.name+" = "+expression);
    }
    return out;
  }
  std::string point_list(const std::vector<int>& p) const {std::string s;for(std::size_t i=0;i<p.size();++i){if(i)s+=",";s+=points_[p[i]].name;}return s;}
  std::set<std::string> required_initial_points(const std::vector<int>& candidate)const{
    ProducerMap producer;
    for(std::size_t i=0;i<construction_commands_.size();++i)
      for(auto output:command_outputs(construction_commands_[i]))producer.emplace(std::move(output),i);
    std::vector<ObjectRef> todo;todo.reserve(candidate.size());
    for(int p:candidate)todo.push_back({ObjectKind::point,points_[static_cast<std::size_t>(p)].name});
    std::set<ObjectRef> seen;std::set<std::string> required;
    while(!todo.empty()){
      ObjectRef object=std::move(todo.back());todo.pop_back();
      if(!seen.insert(object).second)continue;
      if(object.kind==ObjectKind::point&&input_point_names_.count(object.name))required.insert(object.name);
      auto found=producer.find(object);if(found==producer.end())continue;
      for(auto input:command_inputs(construction_commands_[found->second]))todo.push_back(std::move(input));
    }
    return required;
  }
  std::string unused_initial_annotation(const std::vector<int>& candidate)const{
    auto required=required_initial_points(candidate);std::string names;
    for(const auto&name:input_point_names_)if(!required.count(name)){
      if(!names.empty())names+=",";
      names+=name;
    }
    return names.empty()?"":" [unused_initial="+names+"]";
  }
  bool numerically_parallel(int a,int b,int c,int d)const{
    long double ux=points_[b].x-points_[a].x,uy=points_[b].y-points_[a].y;
    long double vx=points_[d].x-points_[c].x,vy=points_[d].y-points_[c].y;
    return std::fabs(ux*vy-uy*vx)<=10*EPS*(1+std::hypotl(ux,uy)*std::hypotl(vx,vy));
  }
  bool numerically_right_angle(int vertex,int a,int b)const{
    long double ux=points_[a].x-points_[vertex].x,uy=points_[a].y-points_[vertex].y;
    long double vx=points_[b].x-points_[vertex].x,vy=points_[b].y-points_[vertex].y;
    return std::fabs(ux*vx+uy*vy)<=10*EPS*(1+std::hypotl(ux,uy)*std::hypotl(vx,vy));
  }
  std::string cyclic_shape_annotation(const Candidate&candidate)const{
    if(candidate.kind!="concyclic"||candidate.points.size()<4)return "";
    constexpr std::size_t annotation_limit=16;
    std::set<std::vector<int>> trapezoid_sets,right_sets;
    std::map<long long,std::vector<std::pair<int,int>>> chords;
    for(std::size_t i=0;i<candidate.points.size();++i)for(std::size_t j=i+1;j<candidate.points.size();++j){
      int a=candidate.points[i],b=candidate.points[j];long double angle=std::atan2(points_[b].y-points_[a].y,points_[b].x-points_[a].x);
      while(angle<0)angle+=PI;
      while(angle>=PI)angle-=PI;
      chords[quant(angle)].push_back({a,b});
    }
    // On one circle, disjoint parallel chords are opposite sides of an
    // isosceles trapezoid. Direction buckets avoid enumerating all quadruples.
    for(const auto&[_,same_direction]:chords)for(std::size_t i=0;i<same_direction.size();++i)
      for(std::size_t j=i+1;j<same_direction.size();++j){auto [a,b]=same_direction[i];auto [c,d]=same_direction[j];
        if(a==c||a==d||b==c||b==d||!numerically_parallel(a,b,c,d))continue;
        std::vector<int> q{a,b,c,d};std::sort(q.begin(),q.end());trapezoid_sets.insert(std::move(q));
      }
    // A pair is a diameter exactly when every other circle point sees it at a
    // right angle. Scan endpoint pairs and combine the certified right vertices.
    for(std::size_t i=0;i<candidate.points.size();++i)for(std::size_t j=i+1;j<candidate.points.size();++j){
      int a=candidate.points[i],b=candidate.points[j];std::vector<int> right;
      for(int vertex:candidate.points)if(vertex!=a&&vertex!=b&&numerically_right_angle(vertex,a,b))right.push_back(vertex);
      for(std::size_t u=0;u<right.size();++u)for(std::size_t v=u+1;v<right.size();++v){
        std::vector<int> q{a,b,right[u],right[v]};std::sort(q.begin(),q.end());right_sets.insert(std::move(q));
      }
    }
    auto describe=[&](const std::set<std::vector<int>>&sets){std::vector<std::string> values;std::size_t count=0;
      for(const auto&q:sets){if(count++==annotation_limit){values.push_back("...");break;}values.push_back(point_list(q));}return values;};
    auto trapezoids=describe(trapezoid_sets),right_pairs=describe(right_sets);
    auto join=[](const std::vector<std::string>&values){std::string out;for(const auto&value:values){if(!out.empty())out+=";";out+=value;}return out;};
    std::string result;
    if(!trapezoids.empty())result+=" [isosceles_trapezoid="+join(trapezoids)+"]";
    if(!right_pairs.empty())result+=" [two_right_angles="+join(right_pairs)+"]";
    return result;
  }
  std::string candidate_core(const Candidate&candidate)const{
    return candidate.kind+"("+point_list(candidate.points)+")";
  }
  std::string symmetric_core(const Candidate&candidate)const{
    if(!symmetry_enabled_)return "";
    std::vector<int> dual;dual.reserve(candidate.points.size());
    for(int p:candidate.points){ObjectRef object{ObjectKind::point,points_[static_cast<std::size_t>(p)].name};
      auto found=symmetry_dual_.find(object);if(found==symmetry_dual_.end())return "";
      dual.push_back(pid(found->second.name));
    }
    std::sort(dual.begin(),dual.end());dual.erase(std::unique(dual.begin(),dual.end()),dual.end());
    return candidate.kind+"("+point_list(dual)+")";
  }
  std::string candidate_statement(const Candidate&candidate,const std::set<std::string>&all_candidates)const{
    auto core=candidate_core(candidate),dual=symmetric_core(candidate);std::string symmetry;
    if(!dual.empty()&&dual==core)symmetry=" [symmetry=self]";
    else if(!dual.empty()&&all_candidates.count(dual))
      symmetry=" [symmetry=asymmetric] [symmetric_partner="+dual+"]";
    else if(symmetry_enabled_)
      symmetry=" [symmetry=asymmetric] [symmetric_partner=not_detected]";
    return candidate.kind+"("+point_list(candidate.points)+")"+
           symmetry+unused_initial_annotation(candidate.points)+cyclic_shape_annotation(candidate);
  }
  void print_proof(const std::string& label,const std::set<int>& w)const{std::cout<<"PROVED "<<label<<"\n";int step=1;for(auto&s:angles_.explain(w))std::cout<<"  "<<step++<<". "<<s<<"\n";if(step==1)std::cout<<"  1. direct known fact\n";}
  bool ancestry_proves(const std::string&kind,const std::vector<int>&candidate)const{
    std::map<ObjectRef,std::size_t> producer;
    for(std::size_t i=0;i<construction_commands_.size();++i)
      for(auto output:command_outputs(construction_commands_[i]))producer.emplace(std::move(output),i);
    std::vector<ObjectRef> todo;todo.reserve(candidate.size());
    for(int p:candidate)todo.push_back({ObjectKind::point,points_[static_cast<std::size_t>(p)].name});
    std::set<ObjectRef> seen;std::set<std::size_t> selected;
    while(!todo.empty()){
      ObjectRef object=std::move(todo.back());todo.pop_back();if(!seen.insert(object).second)continue;
      auto it=producer.find(object);if(it==producer.end())
        throw std::runtime_error("missing construction definition for "+object.name);
      if(selected.insert(it->second).second)
        for(auto input:command_inputs(construction_commands_[it->second]))todo.push_back(std::move(input));
    }
    Engine sub(seed_,generation_seed_);sub.record_commands_=true;sub.angles_.set_coefficient_limit(angle_coefficient_limit_);
    for(std::size_t command:selected)sub.execute(construction_commands_[command],0);
    sub.geometry_closure();std::vector<int> local;local.reserve(candidate.size());
    for(int p:candidate)local.push_back(sub.pid(points_[static_cast<std::size_t>(p)].name));
    if(kind=="collinear")return sub.proves_collinear(local);
    if(kind=="concyclic")return sub.proves_cyclic(local);
    throw std::runtime_error("ancestry scope supports collinear and concyclic candidates");
  }
  void run_goals(){for(auto&g:goals_){try{std::set<int>w;bool ok=false;std::string label=g.kind+"(";for(std::size_t i=0;i<g.args.size();++i){if(i)label+=",";label+=g.args[i];}label+=")";
      if(g.kind=="collinear")ok=proves_collinear(names_to_points(g.args),&w);
      else if(g.kind=="concyclic")ok=proves_cyclic(names_to_points(g.args),&w);
      else if((g.kind=="parallel"||g.kind=="perpendicular")&&g.args.size()==4){int a=pid(g.args[0]),b=pid(g.args[1]),c=pid(g.args[2]),d=pid(g.args[3]);ok=angles_.proves(equation({{segment(a,b),1},{segment(c,d),-1}}),g.kind=="perpendicular"?1:0,g.kind=="perpendicular"?2:1,&w);}
      else if(g.kind=="equal_distance"&&g.args.size()==4)ok=prove_equal_distance(pid(g.args[0]),pid(g.args[1]),pid(g.args[2]),pid(g.args[3]));
      else throw std::runtime_error("bad or unsupported proof goal");
      if(ok) print_proof(label,w); else std::cout<<"UNPROVED "<<label<<"\n";
    }catch(const std::exception&e){std::cout<<"ERROR goal: "<<e.what()<<"\n";}}}

  void materialize_goal_segments(){
    // Later closure phases (notably affine direction recovery and component
    // reasoning) can only attach facts to segment carriers that already exist.
    // Proof goals are parsed before closure, so expose all of their point pairs
    // up front instead of first creating them while printing the final verdict.
    for(const auto&g:goals_){std::vector<int> p;
      for(const auto&name:g.args)p.push_back(pid(name));
      for(std::size_t i=0;i<p.size();++i)for(std::size_t j=i+1;j<p.size();++j)
        if(p[i]!=p[j])segment(p[i],p[j]);
    }
  }

 public:
  explicit Engine(std::uint64_t seed,std::uint64_t generation_seed,std::string automatic_namespace={}):
    automatic_namespace_(std::move(automatic_namespace)),seed_(seed),generation_seed_(generation_seed),
    rng_(seed),generation_rng_(generation_seed){}
  void parse(std::istream& in) {std::string line;int no=0;while(std::getline(in,line)){++no;auto hash=line.find('#');if(hash!=std::string::npos)line.resize(hash);std::istringstream ss(line);std::vector<std::string>t;std::string x;while(ss>>x)t.push_back(x);if(t.empty())continue;try{execute(t,no);}catch(const std::exception&e){throw std::runtime_error("line "+std::to_string(no)+": "+e.what());}}}
  void report(bool classify=true){
    input_point_names_.clear();
    for(const auto&command:construction_commands_)
      for(const auto&output:command_outputs(command))if(output.kind==ObjectKind::point)input_point_names_.insert(output.name);
    expand_points();
    if(prove_mode_||!goals_.empty())materialize_goal_segments();
    if(classify||prove_mode_||!goals_.empty())geometry_closure();
    std::cout<<"GEOGEN REPORT\npoints="<<points_.size()<<" lines="<<lines_.size()<<" circles="<<circles_.size()<<"\n";
    for(const auto&definition:point_definitions())std::cout<<definition<<'\n';
    if(prove_mode_||!goals_.empty()){run_goals();return;}
    auto ls=detect_lines();
    if(!classify){
      auto cs=detect_circles();
      std::set<std::string> keys;for(const auto&x:ls)keys.insert(candidate_core(x));for(const auto&x:cs)keys.insert(candidate_core(x));
      for(const auto&x:ls){auto statement=candidate_statement(x,keys);
        std::cout<<"NONTRIVIAL "<<statement<<'\n';if(show_easy_)std::cout<<"EASY "<<statement<<'\n';}
      for(const auto&x:cs){auto statement=candidate_statement(x,keys);
        std::cout<<"NONTRIVIAL "<<statement<<'\n';if(show_easy_)std::cout<<"EASY "<<statement<<'\n';}
      return;
    }
    std::size_t easy=0,hard=0;
    std::set<std::string> keys;for(const auto&x:ls)keys.insert(candidate_core(x));for(const auto&x:circle_cache_)keys.insert(candidate_core(x));
    for(auto&x:ls){std::set<int>w;bool e=ancestry_scope_?ancestry_proves("collinear",x.points):proves_collinear(x.points,&w);auto statement=candidate_statement(x,keys);if(e)++easy;else{++hard;std::cout<<"NONTRIVIAL "<<statement<<"\n";}if(e&&show_easy_)std::cout<<"EASY "<<statement<<"\n";}
    for(auto&x:circle_cache_){std::set<int>w;bool e=ancestry_scope_?ancestry_proves("concyclic",x.points):proves_cyclic(x.points,&w);auto statement=candidate_statement(x,keys);if(e)++easy;else{++hard;std::cout<<"NONTRIVIAL "<<statement<<"\n";}if(e&&show_easy_)std::cout<<"EASY "<<statement<<"\n";}
    std::cout<<"summary nontrivial="<<hard<<" filtered_easy="<<easy<<"\n";
  }
};

} // namespace geogen

namespace {

struct RunSettings {
  bool prove=false,show_easy=false,show_all_points=false,symmetry=false,symmetric_only=false;
  int trials=5;
  std::uint64_t seed=0x47454f47454eULL;
};

RunSettings read_settings(const std::string& input) {
  RunSettings s;std::istringstream in(input);std::string line;
  while(std::getline(in,line)){
    auto hash=line.find('#');if(hash!=std::string::npos)line.resize(hash);
    std::istringstream row(line);std::string a,b,c;if(!(row>>a))continue;
    if(a=="mode"&&row>>b)s.prove=(b=="prove");
    if(a=="option"&&row>>b>>c){
      if(b=="trials")s.trials=std::stoi(c);
      else if(b=="seed")s.seed=std::stoull(c);
      else if(b=="show_easy")s.show_easy=std::stoi(c)!=0;
      else if(b=="show_all_points")s.show_all_points=std::stoi(c)!=0;
      else if(b=="symmetric_coincidences_only")s.symmetric_only=std::stoi(c)!=0;
      else if(b=="symmetry")s.symmetry=true;
    }
  }
  if(s.trials<1||s.trials>100)throw std::runtime_error("option trials must be between 1 and 100");
  if(s.symmetric_only&&!s.symmetry)throw std::runtime_error("option symmetric_coincidences_only requires option symmetry");
  return s;
}

std::string execute_once(const std::string& input,std::uint64_t seed,
                         std::uint64_t generation_seed,bool classify=true,
                         std::string automatic_namespace={}) {
  geogen::Engine e(seed,generation_seed,std::move(automatic_namespace));std::istringstream in(input);e.parse(in);
  std::ostringstream captured;auto* old=std::cout.rdbuf(captured.rdbuf());
  try{e.report(classify);std::cout.rdbuf(old);}catch(...){std::cout.rdbuf(old);throw;}
  return captured.str();
}

std::set<std::string> findings(const std::string& report,bool show_easy,bool symmetric_only) {
  std::set<std::string> out;std::istringstream in(report);std::string line;
  while(std::getline(in,line))if((line.rfind("NONTRIVIAL ",0)==0||
      (show_easy&&line.rfind("EASY ",0)==0))&&
      (!symmetric_only||line.find("[symmetry=self]")!=std::string::npos))out.insert(line);
  return out;
}

std::vector<std::string> point_listing(const std::string& report){
  std::vector<std::string> out;std::istringstream in(report);std::string line;
  while(std::getline(in,line))if(line.rfind("POINT ",0)==0)out.push_back(line);
  return out;
}

std::set<std::string> directly_reported_points(const std::set<std::string>& report_lines){
  std::set<std::string> out;
  for(const auto&line:report_lines){
    auto open=line.find('('),close=line.find(')',open==std::string::npos?0:open+1);
    if(open==std::string::npos||close==std::string::npos)continue;
    std::istringstream names(line.substr(open+1,close-open-1));std::string name;
    while(std::getline(names,name,','))if(!name.empty())out.insert(name);
  }
  return out;
}

std::string listed_point_name(const std::string&definition){
  constexpr std::size_t prefix=6;auto end=definition.find(" = ",prefix);
  return end==std::string::npos?"":definition.substr(prefix,end-prefix);
}

std::string symmetry_group_key(const std::string&line){
  constexpr std::string_view marker="[symmetric_partner=";auto marked=line.find(marker);
  if(marked==std::string::npos)return "~"+line;
  auto partner_begin=marked+marker.size(),partner_end=line.find(']',partner_begin);
  auto prefix_end=line.find(' '),annotation=line.find(" [",prefix_end+1);
  if(partner_end==std::string::npos||prefix_end==std::string::npos)return "~"+line;
  std::string own=line.substr(prefix_end+1,(annotation==std::string::npos?line.size():annotation)-(prefix_end+1));
  std::string partner=line.substr(partner_begin,partner_end-partner_begin);
  if(partner=="not_detected")return "~"+line;
  return line.substr(0,prefix_end)+" "+std::min(own,partner);
}

int symmetry_sort_rank(const std::string&line){
  if(line.find("[symmetry=self]")!=std::string::npos)return 0;
  if(line.find("[symmetry=asymmetric]")!=std::string::npos&&
     line.find("[symmetric_partner=not_detected]")==std::string::npos)return 1;
  if(line.find("[symmetry=asymmetric]")!=std::string::npos)return 2;
  return 3;
}

int shape_sort_rank(const std::string&line){
  return line.find("[isosceles_trapezoid=")!=std::string::npos||
         line.find("[two_right_angles=")!=std::string::npos;
}

} // namespace

int main(int argc,char**argv){try{
  std::ios::sync_with_stdio(false);std::cin.tie(nullptr);
  if(argc>2){std::cerr<<"usage: geogen [input.geogen]\n";return 2;}
  std::string input;
  if(argc==2){std::ifstream f(argv[1]);if(!f)throw std::runtime_error("cannot open input file");input.assign(std::istreambuf_iterator<char>(f),{});}
  else input.assign(std::istreambuf_iterator<char>(std::cin),{});
  RunSettings settings=read_settings(input);
  if(settings.prove){std::cout<<execute_once(input,settings.seed,settings.seed);return 0;}
  std::set<std::string> combined;
  std::vector<std::string> combined_points;
  std::map<std::string,std::string> point_by_name;
  for(int trial=0;trial<settings.trials;++trial){
    std::uint64_t trial_seed=settings.seed+0x9e3779b97f4a7c15ULL*static_cast<std::uint64_t>(trial+1);
    std::uint64_t generation_seed=settings.seed+0xd1b54a32d192ed03ULL*static_cast<std::uint64_t>(trial);
    std::string automatic_namespace=settings.trials==1?"":"T"+std::to_string(trial+1)+"$";
    std::string report=execute_once(input,trial_seed,generation_seed,true,std::move(automatic_namespace));
    auto current=findings(report,settings.show_easy,settings.symmetric_only);auto current_points=point_listing(report);
    combined.insert(current.begin(),current.end());
    for(auto&definition:current_points){auto name=listed_point_name(definition);auto [where,inserted]=point_by_name.emplace(name,definition);
      if(inserted)combined_points.push_back(std::move(definition));
      else if(where->second!=definition)throw std::runtime_error("conflicting definitions for point "+name);
    }
  }
  std::size_t total_points=combined_points.size();
  if(settings.symmetric_only||!settings.show_all_points){auto reported_points=directly_reported_points(combined);
    combined_points.erase(std::remove_if(combined_points.begin(),combined_points.end(),[&](const std::string&definition){
      return !reported_points.count(listed_point_name(definition));
    }),combined_points.end());
  }
  std::cout<<"GEOGEN REPORT\nrandom_trials="<<settings.trials<<" seed="<<settings.seed
           <<"\ngenerated_points="<<total_points<<" reported_points="<<combined_points.size()
           <<" hidden_points="<<(total_points-combined_points.size())<<"\n";
  for(const auto& point:combined_points)std::cout<<point<<'\n';
  std::vector<std::string> ordered(combined.begin(),combined.end());
  std::sort(ordered.begin(),ordered.end(),[](const std::string&a,const std::string&b){
    auto ra=symmetry_sort_rank(a),rb=symmetry_sort_rank(b);if(ra!=rb)return ra<rb;
    auto sa=shape_sort_rank(a),sb=shape_sort_rank(b);if(sa!=sb)return sa<sb;
    auto ga=symmetry_group_key(a),gb=symmetry_group_key(b);return ga==gb?a<b:ga<gb;
  });
  for(const auto& line:ordered)std::cout<<line<<'\n';
  std::cout<<"summary combined_coincidences="<<combined.size();
  if(settings.symmetry){std::size_t self=0,asymmetric=0,unpaired=0;
    for(const auto&line:combined){self+=line.find("[symmetry=self]")!=std::string::npos;
      asymmetric+=line.find("[symmetry=asymmetric]")!=std::string::npos;
      unpaired+=line.find("[symmetric_partner=not_detected]")!=std::string::npos;}
    std::cout<<" symmetric_self="<<self<<" asymmetric_statements="<<asymmetric
             <<" unpaired_asymmetric="<<unpaired;
  }
  std::cout<<"\n";
  return 0;
}catch(const std::exception&e){std::cerr<<"geogen: "<<e.what()<<'\n';return 1;}}
