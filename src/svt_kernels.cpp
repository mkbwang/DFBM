// Elementwise kernels over the mask stack for bmfsvt() and the fold-in.
//
// Each function reproduces the R expression it replaces operation for
// operation, so results agree with the R versions to the last bit or to
// summation roundoff:
//   - log1exp(x) is log(1 + exp(x)) for x < 20 and x otherwise, as in R;
//   - the sigmoid is 1 / (1 + exp(-x)), which is what stats::plogis() computes;
//   - sums accumulate in long double, as R's sum() does, are cast to double
//     per mask, and are added across masks in ascending t, as the R loops do.
// The R wrappers live in R/svt_utils.R; see there for the mathematics.
//
// Arguments arrive as Rcpp types, which check and, if needed, coerce them. The
// loops then run over raw pointers from begin(): unlike operator[], a pointer
// read stays cheap when devtools::load_all() compiles without optimization.

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

static inline double log1exp1(double x) {
  return x < 20 ? std::log(1 + std::exp(x)) : x;
}

static inline double sigmoid1(double x) {
  return 1 / (1 + std::exp(-x));
}

// A list of equally sized matrices, each converted to a NumericMatrix once.
// The matrices are held, not only their pointers, so that a coerced copy (an
// integer matrix, say) lives as long as the pointers into it are used.
struct MatrixStack {
  std::vector<NumericMatrix> mats;// vector of Numeric Matrices
  std::vector<const double*> p;// vector of pointers
  int nr, nc; // number of rows and number of columns
  R_xlen_t n; // total number of entries

  MatrixStack(const List& L, const char* what) {
    int Tn = L.size();
    if (Tn < 1) stop("`%s` is empty.", what);
    mats.reserve(Tn);
    p.reserve(Tn);
    for (int t = 0; t < Tn; ++t) {
      NumericMatrix m = L[t];
      if (t > 0 && (m.nrow() != nr || m.ncol() != nc)) {
        stop("every element of `%s` must have the same dimensions.", what);
      }
      nr = m.nrow();
      nc = m.ncol();
      mats.push_back(m);
      p.push_back(m.begin());
    }
    n = (R_xlen_t) nr * nc;
  }

  int size() const { return (int) p.size(); }
};

// The training mask, or an empty vector when every entry is used. `w` is set
// to its data, or to NULL. The caller keeps the returned vector alive while
// it uses `w`.
static NumericVector train_vector(const Nullable<NumericMatrix>& train,
                                  R_xlen_t n, const double*& w) {
  if (train.isNull()) {
    w = NULL;
    return NumericVector(0);
  }
  NumericVector tr(train.get());
  if (tr.size() != n) stop("`train` must have the dimensions of the masks.");
  w = tr.begin(); // pointer to the start of the training mask
  return tr;
}

// [[Rcpp::export]]
double logistic_ce_cpp(NumericVector X, NumericVector Y,
                       Nullable<NumericMatrix> train) {
  R_xlen_t n = X.size();
  if (Y.size() != n) stop("`X` and `Y` must have the same size.");
  const double* x = X.begin();
  const double* y = Y.begin();
  const double* w;
  NumericVector tr = train_vector(train, n, w);
  long double s = 0;
  if (w == NULL) {
    for (R_xlen_t i = 0; i < n; ++i) s += log1exp1(x[i]) - y[i] * x[i];
  } else {
    for (R_xlen_t i = 0; i < n; ++i) s += w[i] * (log1exp1(x[i]) - y[i] * x[i]);
  }
  return (double) s;
}

// T x ncol table of offsets, off[t * nc + j]: mu0 for t = 0, nu_t after.
// Column offsets (`mu0` of length > 1 or a matrix `nu`) vary with j; scalar
// ones are repeated, which is what rep(x, each = n) and recycling give in R.
static std::vector<double> offset_table(const NumericVector& mu0,
                                        const NumericVector& nu, int Tn,
                                        int nc) {
  bool nu_mat = nu.hasAttribute("dim");
  bool column = nu_mat || mu0.size() > 1;
  int nu_rows = 0, nu_cols = 0;
  if (nu_mat) {
    IntegerVector dim = nu.attr("dim");
    if (dim.size() != 2) stop("`nu` must be a vector or a matrix.");
    nu_rows = dim[0];
    nu_cols = dim[1];
  }
  if (Tn > 1 && (nu_mat ? (nu_rows < Tn - 1 || nu_cols != nc)
                        : nu.size() < Tn - 1)) {
    stop("`nu` must have T-1 elements, or T-1 rows and ncol(Z) columns.");
  }
  if (mu0.size() < 1 || (column && mu0.size() != 1 && mu0.size() != nc)) {
    stop("`mu0` must have length 1 or ncol(Z).");
  }
  std::vector<double> off((size_t) Tn * nc);
  for (int j = 0; j < nc; ++j) {
    off[j] = mu0.size() == 1 ? mu0[0] : mu0[j];
  }
  const double* nup = nu.begin();
  for (int t = 1; t < Tn; ++t) {
    for (int j = 0; j < nc; ++j) {
      off[(size_t) t * nc + j] =
        nu_mat ? nup[(t - 1) + (R_xlen_t) nu_rows * j] : nup[t - 1];
    }
  }
  return off;
}

// [[Rcpp::export]]
List stack_forward_cpp(List Z, NumericVector mu0, NumericVector nu) {
  MatrixStack zs(Z, "Z");
  int Tn = zs.size(), nr = zs.nr, nc = zs.nc;
  std::vector<double> off = offset_table(mu0, nu, Tn, nc);
  List X(Tn);
  const double* prev = NULL;
  for (int t = 0; t < Tn; ++t) {
    const double* z = zs.p[t];
    NumericMatrix xt(nr, nc);
    double* o = xt.begin();
    for (int j = 0; j < nc; ++j) {
      double v = off[(size_t) t * nc + j];
      for (int i = 0; i < nr; ++i) {
        R_xlen_t k = (R_xlen_t) j * nr + i;
        o[k] = t == 0 ? z[k] + v : prev[k] + z[k] + v;
      }
    }
    X[t] = xt;
    prev = o;
  }
  return X;
}

// Reverse cumulative residuals of X against Y.
// [[Rcpp::export]]
List grad_backward_cpp(List X, List Y, Nullable<NumericMatrix> train) {
  MatrixStack xs(X, "X"), ys(Y, "Y");
  int Tn = xs.size();
  if (ys.size() != Tn || ys.n != xs.n) {
    stop("`X` and `Y` must have the same length and dimensions.");
  }
  R_xlen_t n = xs.n;
  const double* w;
  NumericVector tr = train_vector(train, n, w);
  List Psi(Tn);
  const double* run = NULL;
  for (int t = Tn - 1; t >= 0; --t) {
    const double* x = xs.p[t];
    const double* y = ys.p[t];
    NumericMatrix pt(xs.nr, xs.nc);
    double* o = pt.begin();
    for (R_xlen_t i = 0; i < n; ++i) {
      double r = sigmoid1(x[i]) - y[i];
      if (w != NULL) r = r * w[i];
      o[i] = run == NULL ? r : run[i] + r;
    }
    Psi[t] = pt;
    run = o;
  }
  return Psi;
}

// Scratch space for the natural parameter stack. It lives as long as the
// process and is reused across calls, so the stack built inside the kernels
// below never goes through R's allocator or garbage collector. R calls these
// from one thread; each forked worker gets its own copy.
static std::vector<double>& scratch(size_t size) {
  static std::vector<double> buf;
  if (buf.size() < size) buf.resize(size);
  return buf;
}

// X^(t) = mu_t + sum_{t' <= t} Z^(t') into x + t * n, in the order of
// operations of stack_forward().
static void forward_raw(const MatrixStack& zs, const std::vector<double>& off,
                        double* x) {
  int Tn = zs.size(), nr = zs.nr, nc = zs.nc;
  R_xlen_t n = zs.n;
  for (int t = 0; t < Tn; ++t) {
    double* o = x + (size_t) t * n;
    const double* prev = o - n;
    const double* z = zs.p[t];
    for (int j = 0; j < nc; ++j) {
      double v = off[(size_t) t * nc + j];
      R_xlen_t k0 = (R_xlen_t) j * nr;
      if (t == 0) {
        for (R_xlen_t k = k0; k < k0 + nr; ++k) o[k] = z[k] + v;
      } else {
        for (R_xlen_t k = k0; k < k0 + nr; ++k) o[k] = prev[k] + z[k] + v;
      }
    }
  }
}

// Gradient Psi of the stacked cross entropy at the stack built from the
// increments Z, i.e. grad_backward(stack_forward(Z, mu0, nu), Y, train), and
// with `ce` also stack_ce() at the same point. The stack itself is built in
// scratch space and never returned.
// [[Rcpp::export]]
List stack_grad_cpp(List Z, NumericVector mu0, NumericVector nu, List Y,
                    Nullable<NumericMatrix> train, bool ce) {
  MatrixStack zs(Z, "Z"), ys(Y, "Y");
  int Tn = zs.size(), nr = zs.nr, nc = zs.nc;
  if (ys.size() != Tn || ys.n != zs.n) {
    stop("`Z` and `Y` must have the same length and dimensions.");
  }
  R_xlen_t n = zs.n;
  std::vector<double> off = offset_table(mu0, nu, Tn, nc);
  const double* w;
  NumericVector tr = train_vector(train, n, w);
  double* xs = scratch((size_t) Tn * n).data();
  forward_raw(zs, off, xs);

  List Psi(Tn);
  double f = 0;
  std::vector<double> cet(Tn);
  const double* run = NULL;
  for (int t = Tn - 1; t >= 0; --t) {
    const double* x = xs + (size_t) t * n;
    const double* y = ys.p[t];
    NumericMatrix pt(nr, nc);
    double* o = pt.begin();
    if (ce) {
      // One exp serves both terms below the log1exp cutoff, where exp(x) is
      // finite; e / (1 + e) is the sigmoid to within an ulp or two.
      long double s = 0;
      for (R_xlen_t i = 0; i < n; ++i) {
        double xi = x[i], sg, l;
        if (xi < 20) {
          double e = std::exp(xi);
          l = std::log(1 + e);
          sg = e / (1 + e);
        } else {
          l = xi;
          sg = sigmoid1(xi);
        }
        double r = sg - y[i];
        double v = l - y[i] * xi;
        if (w != NULL) {
          r = r * w[i];
          v = w[i] * v;
        }
        o[i] = run == NULL ? r : run[i] + r;
        s += v;
      }
      cet[t] = (double) s;
    } else {
      for (R_xlen_t i = 0; i < n; ++i) {
        double r = sigmoid1(x[i]) - y[i];
        if (w != NULL) r = r * w[i];
        o[i] = run == NULL ? r : run[i] + r;
      }
    }
    Psi[t] = pt;
    run = o;
  }
  if (ce) {
    for (int t = 0; t < Tn; ++t) f = f + cet[t];
  } else {
    f = NA_REAL;
  }
  return List::create(_["Psi"] = Psi, _["f"] = f); // gradient matricces and sum of cross entropy losses
}

// Cross entropy f at the stack built from Zp and the backtracking bound
// Q = f_W + sum_t <Zp_t - W_t, Psi_t> + L_t / 2 ||Zp_t - W_t||^2, in one pass
// per mask. The running stack is kept in one scratch layer.
// [[Rcpp::export]]
List step_eval_cpp(List Zp, List W, List Psi, NumericVector L,
                   NumericVector mu0, NumericVector nu, List Y,
                   Nullable<NumericMatrix> train, double f_W) {
  MatrixStack zs(Zp, "Zp"), ws(W, "W"), ps(Psi, "Psi"), ys(Y, "Y");
  int Tn = zs.size(), nr = zs.nr, nc = zs.nc;
  if (ws.size() != Tn || ps.size() != Tn || ys.size() != Tn || L.size() != Tn) {
    stop("`Zp`, `W`, `Psi`, `Y` and `L` must have the same length.");
  }
  if (ws.n != zs.n || ps.n != zs.n || ys.n != zs.n) {
    stop("`Zp`, `W`, `Psi` and `Y` must have the same dimensions.");
  }
  R_xlen_t n = zs.n;
  std::vector<double> off = offset_table(mu0, nu, Tn, nc);
  const double* w;
  NumericVector tr = train_vector(train, n, w);
  double* x = scratch((size_t) n).data();

  double f = 0;
  double quad = f_W;
  for (int t = 0; t < Tn; ++t) {
    const double* z = zs.p[t];
    const double* wt = ws.p[t];
    const double* pt = ps.p[t];
    const double* y = ys.p[t];
    long double sce = 0, s1 = 0, s2 = 0;
    for (int j = 0; j < nc; ++j) {
      double v0 = off[(size_t) t * nc + j];
      R_xlen_t k0 = (R_xlen_t) j * nr;
      for (R_xlen_t k = k0; k < k0 + nr; ++k) {
        double xk = t == 0 ? z[k] + v0 : x[k] + z[k] + v0;
        x[k] = xk;
        double v = log1exp1(xk) - y[k] * xk;
        sce += w == NULL ? v : w[k] * v;
        double d = z[k] - wt[k];
        s1 += d * pt[k];
        s2 += d * d;
      }
    }
    f = f + (double) sce;
    quad = quad + (double) s1 + L[t] / 2 * (double) s2;
  }
  return List::create(_["f"] = f, _["quad"] = quad);
}

// [[Rcpp::export]]
NumericVector fista_stats_cpp(List Zp, List Z, List W) {
  MatrixStack as(Zp, "Zp"), bs(Z, "Z"), cs(W, "W");
  int Tn = as.size();
  if (bs.size() != Tn || cs.size() != Tn) {
    stop("`Zp`, `Z` and `W` must have the same length.");
  }
  if (bs.n != as.n || cs.n != as.n) {
    stop("`Zp`, `Z` and `W` must have the same dimensions.");
  }
  R_xlen_t n = as.n;
  double num = 0, den = 0, ip = 0;
  for (int t = 0; t < Tn; ++t) {
    const double* a = as.p[t];
    const double* b = bs.p[t];
    const double* c = cs.p[t];
    long double sn = 0, sd = 0, si = 0;
    for (R_xlen_t i = 0; i < n; ++i) {
      double d = a[i] - b[i];
      sn += d * d;
      sd += b[i] * b[i];
      si += (c[i] - a[i]) * d;
    }
    num = num + (double) sn;// Frobenius norm of Zp-Z
    den = den + (double) sd;// Fronenious norm of Z
    ip = ip + (double) si;// used to decide whether to reset step size
  }
  return NumericVector::create(_["num"] = num, _["den"] = den, _["ip"] = ip);
}

// [[Rcpp::export]]
List clip_block_cpp(NumericMatrix mat, NumericVector level) {
  int nr = mat.nrow(), nc = mat.ncol();
  if (level.size() != 1 && level.size() != nc) {
    stop("`level` must have length 1 or ncol(mat).");
  }
  NumericMatrix out(nr, nc);
  const double* a = mat.begin();
  double* o = out.begin();
  double count = 0;
  for (int j = 0; j < nc; ++j) {
    double lv = level.size() == 1 ? level[0] : level[j];
    for (int i = 0; i < nr; ++i) {
      R_xlen_t k = (R_xlen_t) j * nr + i;
      if (a[k] > lv) {
        o[k] = lv;
        count += 1;
      } else {
        o[k] = a[k];
      }
    }
  }
  return List::create(_["mat"] = out, _["n_clipped"] = count);
}

// [[Rcpp::export]]
List extrapolate_cpp(List Zp, List Z, double fac) {
  MatrixStack as(Zp, "Zp"), bs(Z, "Z");
  int Tn = as.size();
  if (bs.size() != Tn || bs.n != as.n) {
    stop("`Zp` and `Z` must have the same length and dimensions.");
  }
  R_xlen_t n = as.n;
  List W(Tn);
  for (int t = 0; t < Tn; ++t) {
    const double* a = as.p[t];
    const double* b = bs.p[t];
    NumericMatrix wt(as.nr, as.nc);
    double* o = wt.begin();
    for (R_xlen_t i = 0; i < n; ++i) o[i] = a[i] + fac * (a[i] - b[i]);
    W[t] = wt;
  }
  return W;
}

// [[Rcpp::export]]
NumericMatrix prox_point_cpp(NumericMatrix W, NumericMatrix Psi, double L) {
  R_xlen_t n = W.size();
  if (Psi.size() != n) stop("`W` and `Psi` must have the same size.");
  NumericMatrix out(W.nrow(), W.ncol());
  const double* a = W.begin();
  const double* b = Psi.begin();
  double* o = out.begin();
  for (R_xlen_t i = 0; i < n; ++i) o[i] = a[i] - b[i] / L;
  return out;
}
