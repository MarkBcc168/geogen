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
  using Coeff=std::map<int,Integer>;
  struct Equation {Coeff c;std::set<int> reasons;};

 private:
  static constexpr int HALF_TURN_COLUMN=1000000000;
  std::vector<Equation> rows_;
  std::set<Coeff> row_keys_;
  mutable std::vector<Equation> basis_;
  mutable std::map<Coeff,std::set<int>> proven_cache_;
  mutable bool dirty_=true;
  mutable bool modular_mode_=false;
  struct ModularEquation {
    std::map<int,std::int64_t> c;
    std::map<int,std::int64_t> combination; // coefficients of original fact rows
  };
  static constexpr std::array<std::int64_t,2> MODULI{1000000007LL,1000000009LL};
  mutable std::vector<std::unordered_map<int,ModularEquation>> modular_basis_;
  std::int64_t coefficient_limit_=10000;
  std::vector<std::string> reason_text_;

  static Integer coefficient(const Equation& e,int column){
    auto it=e.c.find(column);return it==e.c.end()?Integer(0):it->second;
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
  static void add_scaled(Equation& x,const Equation& y,const Integer& f){
    if(f==0)return;
    for (auto [v, a] : y.c) {
      Integer value=checked_add(x.c[v],checked_multiply(f,a));
      if(value==0)x.c.erase(v);else x.c[v]=std::move(value);
    }
    x.reasons.insert(y.reasons.begin(), y.reasons.end());
  }
  static std::int64_t mod_norm(std::int64_t x,std::int64_t p){x%=p;if(x<0)x+=p;return x;}
  static std::int64_t mod_mul(std::int64_t a,std::int64_t b,std::int64_t p){return (a*b)%p;}
  static std::int64_t mod_power(std::int64_t a,std::int64_t e,std::int64_t p){std::int64_t r=1;while(e){if(e&1)r=mod_mul(r,a,p);a=mod_mul(a,a,p);e>>=1;}return r;}
  static ModularEquation modularized(const Equation&e,std::int64_t p){
    ModularEquation out;for(const auto&[v,a]:e.c){
      std::int64_t value=mod_norm(static_cast<std::int64_t>(a%p),p);if(value)out.c[v]=value;}
    return out;
  }
  static void modular_add_scaled(ModularEquation&x,const ModularEquation&y,std::int64_t f,std::int64_t p){
    if(!f)return;
    for(auto[v,a]:y.c){auto it=x.c.find(v);auto old=it==x.c.end()?0:it->second;
      auto value=mod_norm(old+mod_mul(f,a,p),p);if(value)x.c[v]=value;else if(it!=x.c.end())x.c.erase(it);}
    for(auto[v,a]:y.combination){auto it=x.combination.find(v);auto old=it==x.combination.end()?0:it->second;
      auto value=mod_norm(old+mod_mul(f,a,p),p);if(value)x.combination[v]=value;else if(it!=x.combination.end())x.combination.erase(it);}
  }
  static void modular_reduce(ModularEquation&e,const std::unordered_map<int,ModularEquation>&basis,std::int64_t p){
    while(!e.c.empty()){int pivot=e.c.begin()->first;auto it=basis.find(pivot);if(it==basis.end())break;modular_add_scaled(e,it->second,mod_norm(-e.c.begin()->second,p),p);}
  }
  static bool modular_member(const Equation&e,const std::unordered_map<int,ModularEquation>&basis,std::int64_t p){
    auto row=modularized(e,p);
    while(!row.c.empty()){
      int pivot=row.c.begin()->first;auto it=basis.find(pivot);if(it==basis.end())return false;
      auto factor=mod_norm(-row.c.begin()->second,p);
      for(auto[v,a]:it->second.c){auto jt=row.c.find(v);auto old=jt==row.c.end()?0:jt->second;
        auto value=mod_norm(old+mod_mul(factor,a,p),p);
        if(value)row.c[v]=value;else if(jt!=row.c.end())row.c.erase(jt);
      }
    }
    return true;
  }
  void modular_insert(const Equation&e,std::size_t mi,int fact_id)const{
    auto p=MODULI[mi];auto row=modularized(e,p);row.combination[fact_id]=1;
    modular_reduce(row,modular_basis_[mi],p);if(row.c.empty())return;
    int pivot=row.c.begin()->first;auto inv=mod_power(row.c.begin()->second,p-2,p);
    for(auto&[_,a]:row.c)a=mod_mul(a,inv,p);
    for(auto&[_,a]:row.combination)a=mod_mul(a,inv,p);
    modular_basis_[mi][pivot]=std::move(row);
  }
  void rebuild_modular()const{
    modular_basis_.assign(MODULI.size(),{});Equation period;period.c[HALF_TURN_COLUMN]=2;
    for(std::size_t mi=0;mi<MODULI.size();++mi){modular_insert(period,mi,-1);for(std::size_t i=0;i<rows_.size();++i)modular_insert(rows_[i],mi,static_cast<int>(i));}
    basis_.clear();modular_mode_=true;dirty_=false;
  }
  static Equation converted(const Coeff& c,std::int64_t numerator,
                            std::int64_t denominator){
    if(denominator==0||(2*numerator)%denominator!=0)
      throw std::runtime_error("angle constant is not an integral multiple of pi/2");
    Equation e;
    // h=pi/2 is deliberately ordered after all line variables. Eliminating this
    // dense torsion column first causes catastrophic coefficient swell.
    Integer half_turns=(2*numerator)/denominator;
    if(half_turns!=0)e.c[HALF_TURN_COLUMN]=-half_turns;
    for(auto [v,a]:c)if(a!=0)e.c[v]=a;
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
    Equation period;period.c[HALF_TURN_COLUMN]=2;work.push_back(std::move(period)); // pi=0
    std::size_t rank=0;std::set<int> occupied;
    for(const auto&row:work)for(const auto&[column,_]:row.c)occupied.insert(column);
    for(int col:occupied){if(rank>=work.size())break;
      std::size_t chosen=rank;while(chosen<work.size()&&coefficient(work[chosen],col)==0)++chosen;
      if(chosen==work.size())continue;
      std::swap(work[rank],work[chosen]);
      for(std::size_t i=rank+1;i<work.size();++i){
        while(coefficient(work[i],col)!=0){
          Integer a=coefficient(work[rank],col),b=coefficient(work[i],col);
          add_scaled(work[rank],work[i],-(a/b));
          std::swap(work[rank],work[i]);
        }
      }
      if(coefficient(work[rank],col)<0)
        for(auto& [_,value]:work[rank].c)value=-value;
      ++rank;
    }
    basis_.assign(work.begin(),work.begin()+static_cast<std::ptrdiff_t>(rank));dirty_=false;
  }

 public:
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
    if(auto it=proven_cache_.find(e.c);it!=proven_cache_.end()){
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
      std::map<int,std::int64_t> expected;for(std::size_t mi=0;mi<MODULI.size();++mi){
      auto p=MODULI[mi];auto row=modularized(e,p);modular_reduce(row,modular_basis_[mi],p);
      std::map<int,std::int64_t> signed_coefficients;
      for(auto[id,value]:row.combination){if(value>p/2)value-=p;if(std::llabs(value)>coefficient_limit_)return false;if(value)signed_coefficients[id]=value;}
      if(mi==0)expected=std::move(signed_coefficients);else if(signed_coefficients!=expected)return false;
    }
      std::set<int> proof_reasons;
      for(const auto&[fact,coefficient]:expected)if(fact>=0&&coefficient)
        proof_reasons.insert(rows_[static_cast<std::size_t>(fact)].reasons.begin(),rows_[static_cast<std::size_t>(fact)].reasons.end());
      proven_cache_[e.c]=proof_reasons;if(reasons)*reasons=std::move(proof_reasons);
      return true;
    }
    for(const auto& row:basis_){
      if(row.c.empty())continue;
      int pivot=row.c.begin()->first;
      Integer value=coefficient(e,pivot);if(value==0)continue;
      Integer divisor=row.c.begin()->second;
      if(value%divisor!=0)return false;
      add_scaled(e,row,-(value/divisor));
    }
    if(!e.c.empty())return false;
    proven_cache_[e.c]=e.reasons;if(reasons)*reasons=std::move(e.reasons);
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
  std::vector<LineReflectionFact> line_reflection_facts_;
  std::map<std::pair<int,int>,int> perpendicular_bisector_loci_;
  std::vector<std::vector<std::string>> construction_commands_;
  std::set<std::string> input_point_names_;
  bool record_commands_=true,ancestry_scope_=false;
  std::int64_t angle_coefficient_limit_=10000;
  std::vector<std::array<int,3>> initial_triangles_;
  bool prove_mode_ = false, show_easy_ = false;
  bool auto_line_circle_=true;
  std::size_t circle_budget_ = 25000000;
  std::size_t max_points_=0;
  std::size_t automatic_serial_=0;
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
    int id = static_cast<int>(lines_.size()); line_id_[l.name] = id; lines_.push_back(std::move(l));
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
    line_points_[static_cast<std::size_t>(id)]={a,b};
    segment_line_[key] = id; return id;
  }
  static AngleSystem::Coeff equation(std::initializer_list<std::pair<int,int>> xs) {
    AngleSystem::Coeff c; for (auto [v,a] : xs) { c[v] += a; if (!c[v]) c.erase(v); } return c;
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
  void parallel_fact(int x, int y, const std::string& why) {
    direction_union(x,y,0);
    angles_.add(equation({{x,1},{y,-1}}), 0, 1, why);
  }
  void perpendicular_fact(int x, int y, const std::string& why) {
    direction_union(x,y,1);
    angles_.add(equation({{x,1},{y,-1}}), 1, 2, why);
  }
  void circumcenter_angle_fact(int o,int a,int b,int c,const std::string& why){
    // angle(ACB)+angle(OAB)=90 degrees, written without dividing a relation.
    angles_.add(equation({{segment(b,c),1},{segment(a,c),-1},
                          {segment(a,b),1},{segment(a,o),-1}}),1,2,why);
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
    if(cyclic) add_cyclic(pid(t[1]),pid(t[2]),pid(t[3]),pid(t[4]),
                          "initial cyclic quadrilateral");
  }

  void execute(const std::vector<std::string>& t, int line_no) {
    auto need=[&](std::size_t n){if(t.size()!=n)throw std::runtime_error(t[0]+" expects "+std::to_string(n-1)+" arguments");};
    const auto& op=t[0];
    if(record_commands_&&is_construction(op))construction_commands_.push_back(t);
    if(op=="mode"){need(2);prove_mode_=t[1]=="prove";if(t[1]!="prove"&&t[1]!="generate")throw std::runtime_error("mode is generate or prove");}
    else if(op=="option"){need(3);if(t[1]=="show_easy")show_easy_=std::stoi(t[2])!=0;else if(t[1]=="circle_budget")circle_budget_=std::stoull(t[2]);else if(t[1]=="angle_coefficient_limit"){angle_coefficient_limit_=std::stoll(t[2]);angles_.set_coefficient_limit(angle_coefficient_limit_);}else if(t[1]=="proof_scope"){if(t[2]=="global")ancestry_scope_=false;else if(t[2]=="ancestry")ancestry_scope_=true;else throw std::runtime_error("option proof_scope is global or ancestry");}else if(t[1]=="line_circle_intersections")auto_line_circle_=std::stoi(t[2])!=0;else if(t[1]=="max_points"){max_points_=std::stoull(t[2]);if(max_points_>5000)throw std::runtime_error("option max_points cannot exceed 5000");if(max_points_&&points_.size()>max_points_)throw std::runtime_error("option max_points is below the number of points already declared");}else if(t[1]=="trials"||t[1]=="seed"){}else throw std::runtime_error("unknown option "+t[1]);}
    else if(op=="triangle") initial_triangle(t);
    else if(op=="quadrilateral") initial_quadrilateral(t,false);
    else if(op=="cyclic_quad") initial_quadrilateral(t,true);
    else if(op=="point"){need(4);add_point({t[1],std::stold(t[2]),std::stold(t[3]),"declared"});}
    else if(op=="line"){need(4);int a=pid(t[2]),b=pid(t[3]);int id=add_line(through(t[1],points_[a],points_[b],op));parallel_fact(id,segment(a,b),"definition line("+t[2]+","+t[3]+")");incidence(a,id,"line incidence");incidence(b,id,"line incidence");}
    else if(op=="midpoint"){need(4);int a=pid(t[2]),b=pid(t[3]);if(!add_point({t[1],(points_[a].x+points_[b].x)/2,(points_[a].y+points_[b].y)/2,op}))return;int m=pid(t[1]);parallel_fact(segment(a,m),segment(a,b),"midpoint collinearity "+t[1]);parallel_fact(segment(b,m),segment(a,b),"midpoint collinearity "+t[1]);inherit_collinearity(m,a,b,"midpoint incidence "+t[1]);equal_length(a,m,m,b,"midpoint lengths");register_midpoint_fact(m,a,b,t[1]);}
    else if(op=="perp_bisector"){need(4);int a=pid(t[2]),b=pid(t[3]);Point m{"",(points_[a].x+points_[b].x)/2,(points_[a].y+points_[b].y)/2,""};Line base=through("",points_[a],points_[b],"");Line l{t[1],base.b,-base.a,-(base.b*m.x-base.a*m.y),op};int id=add_line(l);perpendicular_fact(id,segment(a,b),"perpendicular bisector "+t[1]);perpendicular_bisectors_.push_back({id,a,b});}
    else if(op=="parallel"||op=="perpendicular"){need(4);int p=pid(t[2]),base=lid(t[3]);const Line& q=lines_[base];Line l;if(op=="parallel")l={t[1],q.a,q.b,-q.a*points_[p].x-q.b*points_[p].y,op};else l={t[1],q.b,-q.a,-q.b*points_[p].x+q.a*points_[p].y,op};int id=add_line(l);incidence(p,id,op+" through point");if(op=="parallel")parallel_fact(id,base,"parallel construction "+t[1]);else perpendicular_fact(id,base,"perpendicular construction "+t[1]);}
    else if(op=="angle_bisector"){need(5);int a=pid(t[2]),b=pid(t[3]),c=pid(t[4]);long double ux=points_[a].x-points_[b].x,uy=points_[a].y-points_[b].y,vx=points_[c].x-points_[b].x,vy=points_[c].y-points_[b].y;long double un=std::hypotl(ux,uy),vn=std::hypotl(vx,vy);Point q{"",points_[b].x+ux/un+vx/vn,points_[b].y+uy/un+vy/vn,""};int id=add_line(through(t[1],points_[b],q,op));incidence(b,id,"angle bisector through vertex");angles_.add(equation({{id,2},{segment(a,b),-1},{segment(b,c),-1}}),0,1,"angle bisector "+t[1]);}
    else if(op=="reflection_line"){need(4);int p=pid(t[2]),l=lid(t[3]);auto&q=lines_[l];long double d=q.a*points_[p].x+q.b*points_[p].y+q.c;if(!add_point({t[1],points_[p].x-2*q.a*d,points_[p].y-2*q.b*d,op}))return;int x=pid(t[1]);perpendicular_fact(segment(p,x),l,"line reflection "+t[1]);for(int axis_point:line_points_[static_cast<std::size_t>(l)])if(axis_point!=p&&axis_point!=x)equal_length(axis_point,p,axis_point,x,"reflection-axis equal distances");line_reflection_facts_.push_back({x,p,l});}
    else if(op=="reflection_point"){need(4);int p=pid(t[2]),o=pid(t[3]);if(!add_point({t[1],2*points_[o].x-points_[p].x,2*points_[o].y-points_[p].y,op}))return;int x=pid(t[1]);parallel_fact(segment(p,o),segment(p,x),"point reflection collinearity "+t[1]);inherit_collinearity(x,p,o,"point reflection incidence "+t[1]);equal_length(p,o,o,x,"point reflection lengths");register_midpoint_fact(o,p,x,"point reflection "+t[1]);}
    else if(op=="foot"){need(4);int p=pid(t[2]),l=lid(t[3]);auto&q=lines_[l];long double d=q.a*points_[p].x+q.b*points_[p].y+q.c;if(!add_point({t[1],points_[p].x-q.a*d,points_[p].y-q.b*d,op}))return;int x=pid(t[1]);incidence(x,l,"foot incidence "+t[1]);perpendicular_fact(segment(p,x),l,"foot "+t[1]);foot_facts_.push_back({x,p,l});}
    else if(op=="intersection_ll"){need(4);int a=lid(t[2]),b=lid(t[3]);if(!add_point(intersect(lines_[a],lines_[b],t[1],op)))return;int x=pid(t[1]);incidence(x,a,"intersection incidence "+t[1]+" on "+t[2]);incidence(x,b,"intersection incidence "+t[1]+" on "+t[3]);}
    else if(op=="circumcenter"){need(5);int a=pid(t[2]),b=pid(t[3]),c=pid(t[4]);if(!add_point(circumcenter(points_[a],points_[b],points_[c],t[1],op)))return;int o=pid(t[1]);equal_length(o,a,o,b,"circumcenter radii");equal_length(o,a,o,c,"circumcenter radii");equal_length(o,b,o,c,"circumcenter radii");circumcenter_angle_fact(o,a,b,c,"circumcenter angle theorem at "+t[1]);circumcenter_angle_fact(o,b,c,a,"circumcenter angle theorem at "+t[1]);circumcenter_angle_fact(o,c,a,b,"circumcenter angle theorem at "+t[1]);circumcenter_facts_.push_back({o,a,b,c});}
    else if(op=="orthocenter"){need(5);int a=pid(t[2]),b=pid(t[3]),c=pid(t[4]);Line bc=through("",points_[b],points_[c],""),ac=through("",points_[a],points_[c],"");Line ha{"",bc.b,-bc.a,-bc.b*points_[a].x+bc.a*points_[a].y,""},hb{"",ac.b,-ac.a,-ac.b*points_[b].x+ac.a*points_[b].y,""};if(!add_point(intersect(ha,hb,t[1],op)))return;int h=pid(t[1]);perpendicular_fact(segment(a,h),segment(b,c),"orthocenter altitude 1 "+t[1]);perpendicular_fact(segment(b,h),segment(a,c),"orthocenter altitude 2 "+t[1]);perpendicular_fact(segment(c,h),segment(a,b),"orthocenter closure "+t[1]);orthocenter_facts_.push_back({h,a,b,c});}
    else if(op=="incenter"){need(5);int a=pid(t[2]),b=pid(t[3]),c=pid(t[4]);long double la=std::sqrt(dist2(points_[b],points_[c])),lb=std::sqrt(dist2(points_[a],points_[c])),lc=std::sqrt(dist2(points_[a],points_[b])),s=la+lb+lc;if(!add_point({t[1],(la*points_[a].x+lb*points_[b].x+lc*points_[c].x)/s,(la*points_[a].y+lb*points_[b].y+lc*points_[c].y)/s,op}))return;int i=pid(t[1]);angles_.add(equation({{segment(a,i),2},{segment(a,b),-1},{segment(a,c),-1}}),0,1,"incenter bisector at "+t[2]);angles_.add(equation({{segment(b,i),2},{segment(a,b),-1},{segment(b,c),-1}}),0,1,"incenter bisector at "+t[3]);angles_.add(equation({{segment(c,i),2},{segment(a,c),-1},{segment(b,c),-1}}),0,1,"incenter closure at "+t[4]);incenter_facts_.push_back({i,a,b,c});}
    else if(op=="circle"){need(4);int o=pid(t[2]),p=pid(t[3]);add_circle({t[1],points_[o],dist2(points_[o],points_[p]),op},o);circle_incidence(p,cid(t[1]));}
    else if(op=="circumcircle"){need(5);int a=pid(t[2]),b=pid(t[3]),c=pid(t[4]);Point o=circumcenter(points_[a],points_[b],points_[c],"@"+t[1],op);add_circle({t[1],o,dist2(o,points_[a]),op});int z=cid(t[1]);circle_incidence(a,z);circle_incidence(b,z);circle_incidence(c,z);}
    else if(op=="incircle"){need(6);int i=pid(t[2]),a=pid(t[3]),b=pid(t[4]),c=pid(t[5]);Line ab=through("",points_[a],points_[b],"");long double r=ab.a*points_[i].x+ab.b*points_[i].y+ab.c;add_circle({t[1],points_[i],r*r,op},i);(void)c;}
    else if(op=="intersection_lc_known"){need(5);int l=lid(t[2]),c=cid(t[3]),k=pid(t[4]);auto&ln=lines_[l];auto&cc=circles_[c];if(std::fabs(ln.a*points_[k].x+ln.b*points_[k].y+ln.c)>EPS*10||!near(dist2(cc.center,points_[k]),cc.r2,10))throw std::runtime_error(t[4]+" is not a known intersection of "+t[2]+" and "+t[3]);long double dx=ln.b,dy=-ln.a;long double vx=points_[k].x-cc.center.x,vy=points_[k].y-cc.center.y;long double tt=-2*dot(vx,vy,dx,dy);if(std::fabs(tt)<=EPS)throw std::runtime_error("known line-circle intersection is tangent; no distinct second point");if(!add_point({t[1],points_[k].x+tt*dx,points_[k].y+tt*dy,op}))return;int x=pid(t[1]);incidence(k,l,"known line-circle incidence");incidence(x,l,"line-circle second incidence");circle_incidence(k,c);circle_incidence(x,c);parallel_fact(segment(k,x),l,"line-circle known-root incidence");}
    else if(op=="intersection_cc_known"){need(5);int c1=cid(t[2]),c2=cid(t[3]),k=pid(t[4]);auto&a=circles_[c1];auto&b=circles_[c2];if(!near(dist2(a.center,points_[k]),a.r2,10)||!near(dist2(b.center,points_[k]),b.r2,10))throw std::runtime_error(t[4]+" is not on both circles");Line axis=through("",a.center,b.center,"");long double d=axis.a*points_[k].x+axis.b*points_[k].y+axis.c;if(std::fabs(d)<=EPS)throw std::runtime_error("known circle-circle intersection is tangent; no distinct second point");if(!add_point({t[1],points_[k].x-2*axis.a*d,points_[k].y-2*axis.b*d,op}))return;int x=pid(t[1]);circle_incidence(k,c1);circle_incidence(x,c1);circle_incidence(k,c2);circle_incidence(x,c2);}
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

  void geometry_closure() {
    // Only declared construction incidences enter the proof layer. Numerical
    // discoveries remain conjectures and therefore cannot prove themselves.
    normalize_definition_incidences();
    register_center_loci();
    normalize_definition_incidences();
    register_incenter_loci();
    normalize_definition_incidences();
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
    // Feet from an incenter to two adjacent sidelines are tangency points from
    // the same vertex. Register both equal tangent lengths and equal inradii;
    // these are the two equal-side pairs of the corresponding kite.
    auto foot_on_side=[&](int source,int a,int b){for(const auto&f:foot_facts_)if(f.source==source){
      const auto&on=line_points_[static_cast<std::size_t>(f.line)];if(std::find(on.begin(),on.end(),a)!=on.end()&&std::find(on.begin(),on.end(),b)!=on.end())return f.foot;}return -1;};
    for(const auto&f:incenter_facts_){int ab=foot_on_side(f.center,f.a,f.b),ac=foot_on_side(f.center,f.a,f.c),bc=foot_on_side(f.center,f.b,f.c);
      auto tangent_pair=[&](int vertex,int x,int y){if(x<0||y<0||x==y)return;equal_length(vertex,x,vertex,y,"equal tangents from vertex to incircle");equal_length(f.center,x,f.center,y,"equal inradii to tangent points");};
      tangent_pair(f.a,ab,ac);tangent_pair(f.b,ab,bc);tangent_pair(f.c,ac,bc);
    }
    circle_cache_ = detect_circles(true);
    bool changed=true;int rounds=0;
    while(changed&&rounds++<8){changed=false;
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

      for(const auto&x:circle_cache_)if(x.points.size()>=4){std::set<int>w;if(proves_cyclic(x.points,&w)){auto before=cyclic_facts_.size();add_cyclic(x.points[0],x.points[1],x.points[2],x.points[3],"converse cyclic angle theorem");changed|=cyclic_facts_.size()!=before;}}
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
        if(!angles_.proves(isosceles,0,1)){angles_.add(isosceles,0,1,"isosceles triangle theorem "+points_[vertex].name+points_[b].name+points_[c].name);changed=true;}
      }
      for(const auto& [base,vertices]:wings)for(auto ai=vertices.begin();ai!=vertices.end();++ai)for(auto di=std::next(ai);di!=vertices.end();++di){
        std::string why="kite theorem "+points_[*ai].name+points_[base.first].name+points_[*di].name+points_[base.second].name;
        int ad=segment(*ai,*di),bc=segment(base.first,base.second);if(!angles_.proves(equation({{ad,1},{bc,-1}}),1,2)){perpendicular_fact(ad,bc,why);changed=true;}
        auto symmetry=equation({{segment(*ai,base.first),1},{segment(*ai,base.second),1},{segment(*di,base.first),-1},{segment(*di,base.second),-1}});
        if(!angles_.proves(symmetry,0,1)){angles_.add(symmetry,0,1,why+" [reflection angles]");changed=true;}
      }

    }
  }

  std::string automatic_name(const std::string&prefix){
    return prefix+"$"+std::to_string(automatic_serial_++);
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
    std::vector<int> ids(points_.size());std::iota(ids.begin(),ids.end(),0);
    return depth_weighted_pick(ids,point_depth_);
  }
  std::vector<int> random_distinct_points(std::size_t count,bool noncollinear=false){
    if(points_.size()<count)return {};
    for(int attempt=0;attempt<64;++attempt){
      std::vector<int> available(points_.size());std::iota(available.begin(),available.end(),0);std::vector<int> out;
      while(out.size()<count){int p=depth_weighted_pick(available,point_depth_);out.push_back(p);available.erase(std::find(available.begin(),available.end(),p));}
      if(!noncollinear||count<3||std::fabs(cross(points_[out[0]],points_[out[1]],points_[out[2]]))>
          EPS*scale(points_[out[0]],points_[out[1]],points_[out[2]]))return out;
    }
    return {};
  }
  std::vector<int> named_lines()const{
    std::vector<int> out;for(std::size_t i=0;i<lines_.size();++i)
      if(lines_[i].origin!="segment"&&!lines_[i].name.empty()&&lines_[i].name[0]!='@')out.push_back(static_cast<int>(i));
    return out;
  }
  std::vector<int> named_circles()const{
    std::vector<int> out;for(std::size_t i=0;i<circles_.size();++i)
      if(!circles_[i].name.empty()&&circles_[i].name[0]!='@')out.push_back(static_cast<int>(i));
    return out;
  }
  int create_random_line(){
    if(points_.size()<2)return -1;
    for(int attempt=0;attempt<64;++attempt){
      int kind=std::uniform_int_distribution<int>(0,4)(generation_rng_);
      auto pair=random_distinct_points(2);if(pair.empty())return -1;
      std::string name;int depth=1+std::max(point_depth_[pair[0]],point_depth_[pair[1]]);
      if(kind==0){name=automatic_name("L");execute({"line",name,points_[pair[0]].name,points_[pair[1]].name},0);}
      else if(kind==1){name=automatic_name("PB");execute({"perp_bisector",name,points_[pair[0]].name,points_[pair[1]].name},0);}
      else if(kind==2||kind==3){
        auto pool=named_lines();if(pool.empty()){kind=0;name=automatic_name("L");execute({"line",name,points_[pair[0]].name,points_[pair[1]].name},0);}
        else {int base=depth_weighted_pick(pool,line_depth_);int p=random_point();depth=1+std::max(point_depth_[p],line_depth_[base]);
          name=automatic_name(kind==2?"Par":"Perp");execute({kind==2?"parallel":"perpendicular",name,points_[p].name,lines_[base].name},0);}
      } else {
        auto triple=random_distinct_points(3,true);if(triple.empty())continue;
        depth=1+std::max({point_depth_[triple[0]],point_depth_[triple[1]],point_depth_[triple[2]]});
        name=automatic_name("Bis");execute({"angle_bisector",name,points_[triple[0]].name,points_[triple[1]].name,points_[triple[2]].name},0);
      }
      int id=lid(name);line_depth_[static_cast<std::size_t>(id)]=depth;return id;
    }
    return -1;
  }
  int create_random_circle(){
    if(points_.size()<2)return -1;
    for(int attempt=0;attempt<64;++attempt){
      int kind=std::uniform_int_distribution<int>(0,2)(generation_rng_);
      std::string name;int depth=0;
      if(kind==1){
        auto triple=random_distinct_points(3,true);if(triple.empty())continue;
        depth=1+std::max({point_depth_[triple[0]],point_depth_[triple[1]],point_depth_[triple[2]]});
        name=automatic_name("Circ");execute({"circumcircle",name,points_[triple[0]].name,points_[triple[1]].name,points_[triple[2]].name},0);
      } else if(kind==0||incenter_facts_.empty()) {
        auto pair=random_distinct_points(2);int o=pair[0],p=pair[1];depth=1+std::max(point_depth_[o],point_depth_[p]);
        name=automatic_name("Circle");execute({"circle",name,points_[o].name,points_[p].name},0);
      } else {
        const auto&f=incenter_facts_[std::uniform_int_distribution<std::size_t>(0,incenter_facts_.size()-1)(generation_rng_)];
        depth=1+std::max({point_depth_[f.center],point_depth_[f.a],point_depth_[f.b],point_depth_[f.c]});
        name=automatic_name("Inc");execute({"incircle",name,points_[f.center].name,points_[f.a].name,points_[f.b].name,points_[f.c].name},0);
      }
      int id=cid(name);circle_depth_[static_cast<std::size_t>(id)]=depth;return id;
    }
    return -1;
  }
  bool create_random_point(){
    if(points_.size()<2)return false;
    int kind=std::uniform_int_distribution<int>(0,9)(generation_rng_);
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
      name=automatic_name(kind==5?"Rl":"Foot");execute({kind==5?"reflection_line":"foot",name,points_[p].name,lines_[line].name},0);return finish(depth);
    }
    if(kind==7){
      if(named_lines().size()<2){create_random_line();create_random_line();}
      for(int attempt=0;attempt<32;++attempt){auto pool=named_lines();if(pool.size()<2)return false;
        int a=depth_weighted_pick(pool,line_depth_),b=depth_weighted_pick(pool,line_depth_);if(a==b)continue;
        if(std::fabs(lines_[a].a*lines_[b].b-lines_[b].a*lines_[a].b)<=EPS)continue;
        int depth=1+std::max(line_depth_[a],line_depth_[b]);name=automatic_name("Xll");
        execute({"intersection_ll",name,lines_[a].name,lines_[b].name},0);return finish(depth);
      }return false;
    }
    if(kind==8){
      if(!auto_line_circle_)return false;
      auto line_pool=named_lines();if(line_pool.empty()){create_random_line();line_pool=named_lines();}if(line_pool.empty())return false;
      for(int attempt=0;attempt<32;++attempt){int line=depth_weighted_pick(line_pool,line_depth_);const auto&on=line_points_[line];if(on.empty())continue;
        int k=on[std::uniform_int_distribution<std::size_t>(0,on.size()-1)(generation_rng_)],o=random_point();if(o==k)continue;
        std::string circle_name=automatic_name("Circle");int circle_depth=1+std::max(point_depth_[o],point_depth_[k]);
        execute({"circle",circle_name,points_[o].name,points_[k].name},0);int circle=cid(circle_name);circle_depth_[circle]=circle_depth;
        const auto&ln=lines_[line];const auto&cc=circles_[circle];long double dx=ln.b,dy=-ln.a;
        if(std::fabs(-2*dot(points_[k].x-cc.center.x,points_[k].y-cc.center.y,dx,dy))<=EPS)continue;
        int depth=1+std::max(line_depth_[line],circle_depth_[circle]);name=automatic_name("Xlc");
        execute({"intersection_lc_known",name,lines_[line].name,circles_[circle].name,points_[k].name},0);return finish(depth);
      }return false;
    }
    if(!auto_line_circle_)return false;
    for(int attempt=0;attempt<32;++attempt){
      int k=random_point(),o1=random_point(),o2=random_point();if(k==o1||k==o2||o1==o2)continue;
      if(std::fabs(cross(points_[o1],points_[k],points_[o2]))<=EPS*scale(points_[o1],points_[k],points_[o2]))continue;
      std::string c1_name=automatic_name("Circle"),c2_name=automatic_name("Circle");
      execute({"circle",c1_name,points_[o1].name,points_[k].name},0);execute({"circle",c2_name,points_[o2].name,points_[k].name},0);
      int c1=cid(c1_name),c2=cid(c2_name);circle_depth_[c1]=1+std::max(point_depth_[o1],point_depth_[k]);circle_depth_[c2]=1+std::max(point_depth_[o2],point_depth_[k]);
      int depth=1+std::max(circle_depth_[c1],circle_depth_[c2]);name=automatic_name("Xcc");
      execute({"intersection_cc_known",name,circles_[c1].name,circles_[c2].name,points_[k].name},0);return finish(depth);
    }
    return false;
  }
  void expand_points(){
    auto room=[&]{return !max_points_||points_.size()<max_points_;};
    // Without a point cap, retain the finite known-root scan for explicitly
    // supplied lines and circles. With a cap, known-root intersections join the
    // randomized construction mix below.
    if(auto_line_circle_&&!max_points_){
      std::size_t line_count=lines_.size(),circle_count=circles_.size();
      for(std::size_t l=0;l<line_count;++l){if(lines_[l].origin=="segment")continue;
        for(std::size_t c=0;c<circle_count;++c){auto on_line=line_points_[l],on_circle=circle_points_[c];
          for(int k:on_line)if(std::find(on_circle.begin(),on_circle.end(),k)!=on_circle.end()){
            const auto&ln=lines_[l];const auto&cc=circles_[c];long double dx=ln.b,dy=-ln.a;
            if(std::fabs(-2*dot(points_[k].x-cc.center.x,points_[k].y-cc.center.y,dx,dy))<=EPS)continue;
            std::string name="X("+lines_[l].name+","+circles_[c].name+","+points_[k].name+")";
            if(!point_id_.count(name))execute({"intersection_lc_known",name,lines_[l].name,circles_[c].name,points_[k].name},0);
          }}
      }return;
    }
    if(!max_points_||!room())return;
    std::size_t failed=0,max_failed=std::max<std::size_t>(2000,100*max_points_);
    while(room()&&failed<max_failed){
      // Auxiliary constructions are deliberately interleaved. They make line
      // intersections and known-root circle intersections available without
      // letting the supporting-object count grow independently of the point cap.
      if(named_lines().size()<2||std::uniform_int_distribution<int>(0,3)(generation_rng_)==0)create_random_line();
      if(named_circles().empty()||std::uniform_int_distribution<int>(0,7)(generation_rng_)==0)create_random_circle();
      if(create_random_point())failed=0;else ++failed;
    }
    if(room())std::cerr<<"warning: random construction search exhausted at "<<points_.size()<<" of "<<max_points_<<" points\n";
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
    std::uint64_t triples=(std::uint64_t)n*(n-1)*(n-2)/6;
    if(respect_budget&&triples>circle_budget_){std::cerr<<"warning: general circle scan skipped ("<<triples<<" triples > circle_budget)\n";std::vector<Candidate>out;for(auto&[_,v]:uniq)out.push_back(v);return out;}
    // O(n^3) time and O(n^2) peak memory. Fixing one anchor retains enough
    // information to discover every circle while allowing each hash table to be freed.
    for(int a=0;a<n;++a){
      std::map<std::tuple<long long,long long,long long>,std::set<int>> bins;
      for(int b=0;b<n;++b)if(b!=a)for(int c=b+1;c<n;++c)if(c!=a){
        if(std::fabs(cross(points_[a],points_[b],points_[c]))<=EPS*scale(points_[a],points_[b],points_[c]))continue;
        Point o=circumcenter(points_[a],points_[b],points_[c],"","scan");long double r2=dist2(o,points_[a]);
        auto key=std::make_tuple(quant(o.x),quant(o.y),quant(r2));auto&s=bins[key];s.insert(a);s.insert(b);s.insert(c);
      }
      for(auto&[_,s]:bins)if(s.size()>=4){std::vector<int>v(s.begin(),s.end());Point o=circumcenter(points_[v[0]],points_[v[1]],points_[v[2]],"","verify");bool ok=true;for(int x:v)if(!near(dist2(o,points_[x]),dist2(o,points_[v[0]]),10))ok=false;if(ok)uniq[v]={"concyclic",v,"circle hash"};}
    }
    std::vector<Candidate>out;for(auto&[_,v]:uniq)out.push_back(v);return out;
  }

  using ProducerMap=std::map<ObjectRef,std::size_t>;
  std::string construction_expression(const ObjectRef&object,const ProducerMap&producer,
      std::map<ObjectRef,std::string>&memo,std::set<ObjectRef>&active,bool force=false)const{
    // Points supplied by the input configuration are named lemmas from the
    // user's perspective. Expand their own POINT assignment, but keep them
    // atomic whenever a later generated definition refers to them.
    if(!force&&object.kind==ObjectKind::point&&input_point_names_.count(object.name))return object.name;
    if(auto known=memo.find(object);known!=memo.end())return known->second;
    if(!active.insert(object).second)return object.name;
    auto done=[&](std::string value){active.erase(object);memo[object]=value;return value;};
    auto found=producer.find(object);if(found==producer.end())return done(object.name);
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
    Engine sub(seed_,generation_seed_);sub.record_commands_=false;sub.angles_.set_coefficient_limit(angle_coefficient_limit_);
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

 public:
  explicit Engine(std::uint64_t seed,std::uint64_t generation_seed):
    seed_(seed),generation_seed_(generation_seed),rng_(seed),generation_rng_(generation_seed){}
  void parse(std::istream& in) {std::string line;int no=0;while(std::getline(in,line)){++no;auto hash=line.find('#');if(hash!=std::string::npos)line.resize(hash);std::istringstream ss(line);std::vector<std::string>t;std::string x;while(ss>>x)t.push_back(x);if(t.empty())continue;try{execute(t,no);}catch(const std::exception&e){throw std::runtime_error("line "+std::to_string(no)+": "+e.what());}}}
  void report(bool classify=true){
    input_point_names_.clear();for(const auto&p:points_)input_point_names_.insert(p.name);
    expand_points();
    if(classify||prove_mode_||!goals_.empty())geometry_closure();
    std::cout<<"GEOGEN REPORT\npoints="<<points_.size()<<" lines="<<lines_.size()<<" circles="<<circles_.size()<<"\n";
    for(const auto&definition:point_definitions())std::cout<<definition<<'\n';
    if(prove_mode_||!goals_.empty()){run_goals();return;}
    auto ls=detect_lines();
    if(!classify){
      auto cs=detect_circles();
      for(const auto&x:ls){auto statement="collinear("+point_list(x.points)+")";
        std::cout<<"NONTRIVIAL "<<statement<<'\n';if(show_easy_)std::cout<<"EASY "<<statement<<'\n';}
      for(const auto&x:cs){auto statement="concyclic("+point_list(x.points)+")";
        std::cout<<"NONTRIVIAL "<<statement<<'\n';if(show_easy_)std::cout<<"EASY "<<statement<<'\n';}
      return;
    }
    std::size_t easy=0,hard=0;
    for(auto&x:ls){std::set<int>w;bool e=ancestry_scope_?ancestry_proves("collinear",x.points):proves_collinear(x.points,&w);if(e)++easy;else{++hard;std::cout<<"NONTRIVIAL collinear("<<point_list(x.points)<<")\n";}if(e&&show_easy_)std::cout<<"EASY collinear("<<point_list(x.points)<<")\n";}
    for(auto&x:circle_cache_){std::set<int>w;bool e=ancestry_scope_?ancestry_proves("concyclic",x.points):proves_cyclic(x.points,&w);if(e)++easy;else{++hard;std::cout<<"NONTRIVIAL concyclic("<<point_list(x.points)<<")\n";}if(e&&show_easy_)std::cout<<"EASY concyclic("<<point_list(x.points)<<")\n";}
    std::cout<<"summary nontrivial="<<hard<<" filtered_easy="<<easy<<"\n";
  }
};

} // namespace geogen

namespace {

struct RunSettings {
  bool prove=false,show_easy=false;
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
    }
  }
  if(s.trials<1||s.trials>100)throw std::runtime_error("option trials must be between 1 and 100");
  return s;
}

std::string execute_once(const std::string& input,std::uint64_t seed,
                         std::uint64_t generation_seed,bool classify=true) {
  geogen::Engine e(seed,generation_seed);std::istringstream in(input);e.parse(in);
  std::ostringstream captured;auto* old=std::cout.rdbuf(captured.rdbuf());
  try{e.report(classify);std::cout.rdbuf(old);}catch(...){std::cout.rdbuf(old);throw;}
  return captured.str();
}

std::set<std::string> findings(const std::string& report,bool show_easy) {
  std::set<std::string> out;std::istringstream in(report);std::string line;
  while(std::getline(in,line))if(line.rfind("NONTRIVIAL ",0)==0||
      (show_easy&&line.rfind("EASY ",0)==0))out.insert(line);
  return out;
}

std::vector<std::string> point_listing(const std::string& report){
  std::vector<std::string> out;std::istringstream in(report);std::string line;
  while(std::getline(in,line))if(line.rfind("POINT ",0)==0)out.push_back(line);
  return out;
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
  std::set<std::string> common;
  std::vector<std::string> common_points;
  for(int trial=0;trial<settings.trials;++trial){
    std::uint64_t trial_seed=settings.seed+0x9e3779b97f4a7c15ULL*static_cast<std::uint64_t>(trial+1);
    std::string report=execute_once(input,trial_seed,settings.seed,trial==0);
    auto current=findings(report,settings.show_easy);auto current_points=point_listing(report);
    if(trial==0){common=std::move(current);common_points=std::move(current_points);}
    else {
      std::set<std::string> both;std::set_intersection(common.begin(),common.end(),current.begin(),current.end(),std::inserter(both,both.end()));common=std::move(both);
      std::set<std::string> present(current_points.begin(),current_points.end());
      common_points.erase(std::remove_if(common_points.begin(),common_points.end(),[&](const std::string& p){return !present.count(p);}),common_points.end());
    }
  }
  std::cout<<"GEOGEN REPORT\nrandom_trials="<<settings.trials<<" seed="<<settings.seed<<"\npoints="<<common_points.size()<<"\n";
  for(const auto& point:common_points)std::cout<<point<<'\n';
  for(const auto& line:common)std::cout<<line<<'\n';
  std::cout<<"summary stable_coincidences="<<common.size()<<"\n";
  return 0;
}catch(const std::exception&e){std::cerr<<"geogen: "<<e.what()<<'\n';return 1;}}
