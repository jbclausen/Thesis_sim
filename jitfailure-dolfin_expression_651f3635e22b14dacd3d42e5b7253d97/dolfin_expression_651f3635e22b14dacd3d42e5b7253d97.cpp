
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif

#include <dolfin/function/Expression.h>
#include <dolfin/math/basic.h>
#include <Eigen/Dense>


// cmath functions
using std::cos;
using std::sin;
using std::tan;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cosh;
using std::sinh;
using std::tanh;
using std::exp;
using std::frexp;
using std::ldexp;
using std::log;
using std::log10;
using std::modf;
using std::pow;
using std::sqrt;
using std::ceil;
using std::fabs;
using std::floor;
using std::fmod;
using std::max;
using std::min;

const double pi = DOLFIN_PI;


namespace dolfin
{
  class dolfin_expression_651f3635e22b14dacd3d42e5b7253d97 : public Expression
  {
     public:
       std::shared_ptr<dolfin::GenericFunction> generic_function_p1;
std::shared_ptr<dolfin::GenericFunction> generic_function_p2;


       dolfin_expression_651f3635e22b14dacd3d42e5b7253d97()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          double p2;
            generic_function_p2->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&p2), x);
          double p1;
            generic_function_p1->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&p1), x);
          values[0] = x[1] < 0 + DOLFIN_EPS ? p1 : 0);

       }

       void set_property(std::string name, double _value) override
       {

       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {

       throw std::runtime_error("No such property");
       return 0.0;
       }

       void set_generic_function(std::string name, std::shared_ptr<dolfin::GenericFunction> _value) override
       {
          if (name == "p1") { generic_function_p1 = _value; return; }          if (name == "p2") { generic_function_p2 = _value; return; }
       throw std::runtime_error("No such property");
       }

       std::shared_ptr<dolfin::GenericFunction> get_generic_function(std::string name) const override
       {
          if (name == "p1") return generic_function_p1;          if (name == "p2") return generic_function_p2;
       throw std::runtime_error("No such property");
       }

  };
}

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_651f3635e22b14dacd3d42e5b7253d97()
{
  return new dolfin::dolfin_expression_651f3635e22b14dacd3d42e5b7253d97;
}

