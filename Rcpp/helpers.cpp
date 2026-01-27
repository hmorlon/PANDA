#include <Rcpp.h>
using namespace Rcpp;

//
// grow_init_cpp(mode, initial_capacity)
// [[Rcpp::export]]
List grow_init_cpp(std::string mode, int initial_capacity) {
  SEXP buf;

  if (mode == "integer") {
    buf = IntegerVector(initial_capacity);
  } else if (mode == "numeric") {
    buf = NumericVector(initial_capacity);
  } else if (mode == "logical") {
    buf = LogicalVector(initial_capacity);
  } else if (mode == "character") {
    buf = CharacterVector(initial_capacity);
  } else {
    stop("Unsupported mode");
  }

  return List::create(
    _["buf"]  = buf,
    _["len"]  = 0,
    _["cap"]  = initial_capacity,
    _["mode"] = mode
  );
}

//
// grow_append_cpp(obj, values)
// [[Rcpp::export]]
List grow_append_cpp(List obj, SEXP values) {
  int len = obj["len"];
  int cap = obj["cap"];
  std::string mode = as<std::string>(obj["mode"]);
  SEXP buf = obj["buf"];

  int n = Rf_length(values);
  if (n == 0) return obj;

  int new_len = len + n;

  if (new_len > cap) {
    int new_cap = std::max(cap * 2, new_len);

    if (mode == "integer") {
      IntegerVector oldv(buf);
      IntegerVector newv(new_cap);
      std::copy(oldv.begin(), oldv.begin() + len, newv.begin());
      buf = newv;
    } else if (mode == "numeric") {
      NumericVector oldv(buf);
      NumericVector newv(new_cap);
      std::copy(oldv.begin(), oldv.begin() + len, newv.begin());
      buf = newv;
    } else if (mode == "logical") {
      LogicalVector oldv(buf);
      LogicalVector newv(new_cap);
      std::copy(oldv.begin(), oldv.begin() + len, newv.begin());
      buf = newv;
    } else if (mode == "character") {
      CharacterVector oldv(buf);
      CharacterVector newv(new_cap);
      for (int i=0; i<len; i++) newv[i] = oldv[i];
      buf = newv;
    }

    cap = new_cap;
  }

  // append
  if (mode == "integer") {
    IntegerVector b(buf);
    IntegerVector v(values);
    std::copy(v.begin(), v.end(), b.begin() + len);
  } else if (mode == "numeric") {
    NumericVector b(buf);
    NumericVector v(values);
    std::copy(v.begin(), v.end(), b.begin() + len);
  } else if (mode == "logical") {
    LogicalVector b(buf);
    LogicalVector v(values);
    std::copy(v.begin(), v.end(), b.begin() + len);
  } else if (mode == "character") {
    CharacterVector b(buf);
    CharacterVector v(values);
    for (int i=0; i<n; i++) b[len + i] = v[i];
  }

  obj["buf"] = buf;
  obj["len"] = new_len;
  obj["cap"] = cap;

  return obj;
}

//
// grow_finalize_cpp(obj)
// [[Rcpp::export]]
SEXP grow_finalize_cpp(List obj) {
  int len = obj["len"];
  std::string mode = as<std::string>(obj["mode"]);
  SEXP buf = obj["buf"];

  if (mode == "integer") {
    IntegerVector b(buf);
    IntegerVector out(len);
    std::copy(b.begin(), b.begin() + len, out.begin());
    return out;
  } else if (mode == "numeric") {
    NumericVector b(buf);
    NumericVector out(len);
    std::copy(b.begin(), b.begin() + len, out.begin());
    return out;
  } else if (mode == "logical") {
    LogicalVector b(buf);
    LogicalVector out(len);
    std::copy(b.begin(), b.begin() + len, out.begin());
    return out;
  } else if (mode == "character") {
    CharacterVector b(buf);
    CharacterVector out(len);
    for (int i=0; i<len; i++) out[i] = b[i];
    return out;
  }

  stop("Unsupported mode");
}


//
// reset_grow_cpp(obj)
// [[Rcpp::export]]
List reset_grow_cpp(List obj) {
  obj["len"] = 0;
  return obj;
}

#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
void inplace_cpp_num(NumericMatrix mat,
                     IntegerVector rows,
                     IntegerVector cols,
                     NumericVector vals) {

    const int nr = rows.size();
    const int nc = cols.size();

    // Determine total operations based on broadcast rules
    int n;

    bool broadcast_rows_cols = false;
    bool broadcast_vals = (vals.size() == 1);

    if (nr == nc) {
        n = nr;   // pairwise scatter
    } else if (nr > 1 && nc == 1) {
        n = nr;   // broadcast single col
        broadcast_rows_cols = true;
    } else if (nr == 1 && nc > 1) {
        n = nc;   // broadcast single row
        broadcast_rows_cols = true;
    } else {
        stop("Unsupported input: lengths(rows) and lengths(cols) must match, or one must be length 1");
    }

    if (!broadcast_vals && vals.size() != n) {
        stop("vals must be length 1 (broadcast) or length equal to the number of assignments");
    }

    for (int k = 0; k < n; k++) {
        int i, j;

        if (broadcast_rows_cols) {
            if (nc == 1) {
                // rows vector, single col
                i = rows[k] - 1;
                j = cols[0] - 1;
            } else {
                // single row, cols vector
                i = rows[0] - 1;
                j = cols[k] - 1;
            }
        } else {
            // pairwise
            i = rows[k] - 1;
            j = cols[k] - 1;
        }

        if (i < 0 || i >= mat.nrow() || j < 0 || j >= mat.ncol()) {
            stop("row/col index out of bounds");
        }

        double v = broadcast_vals ? vals[0] : vals[k];
        mat(i, j) = v;  // in-place
    }
}


// Same function as above but for integer-filled matrices
// [[Rcpp::export]]
void inplace_cpp_int(IntegerMatrix mat,
                     IntegerVector rows,
                     IntegerVector cols,
                     IntegerVector vals) {

    const int nr = rows.size();
    const int nc = cols.size();

    // Determine total operations based on broadcast rules
    int n;

    bool broadcast_rows_cols = false;
    bool broadcast_vals = (vals.size() == 1);

    if (nr == nc) {
        n = nr;   // pairwise scatter
    } else if (nr > 1 && nc == 1) {
        n = nr;   // broadcast single col
        broadcast_rows_cols = true;
    } else if (nr == 1 && nc > 1) {
        n = nc;   // broadcast single row
        broadcast_rows_cols = true;
    } else {
        stop("Unsupported input: lengths(rows) and lengths(cols) must match, or one must be length 1");
    }

    if (!broadcast_vals && vals.size() != n) {
        stop("vals must be length 1 (broadcast) or length equal to the number of assignments");
    }

    for (int k = 0; k < n; k++) {
        int i, j;

        if (broadcast_rows_cols) {
            if (nc == 1) {
                // rows vector, single col
                i = rows[k] - 1;
                j = cols[0] - 1;
            } else {
                // single row, cols vector
                i = rows[0] - 1;
                j = cols[k] - 1;
            }
        } else {
            // pairwise
            i = rows[k] - 1;
            j = cols[k] - 1;
        }

        if (i < 0 || i >= mat.nrow() || j < 0 || j >= mat.ncol()) {
            stop("row/col index out of bounds");
        }

        double v = broadcast_vals ? vals[0] : vals[k];
        mat(i, j) = v;  // in-place
    }
}
