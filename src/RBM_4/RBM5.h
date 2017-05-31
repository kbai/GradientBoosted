class RBM {

public:
  int N;
  int n_visible;
  int n_hidden;
  double **W,**vW;
  double **dW,**sW;
  double *hbias,*vhbias;
  double *dh,*sh;
  double *vbias,*vvbias;
  double *dv,*sv;
  double *bif,*but,*bu,*vbif,*vbut,*vbu;
  RBM(int, int, int, double**, double*, double*,vector<vector<float>>&);
  RBM(const RBM&A);
  ~RBM();
  void operator=(const RBM&A);
  void contrastive_divergence(const vector<int>& ,const vector<int>&,const double,const double,const int);
  void mpigath();
  void sample_h_given_v(const int*,const vector<int>&, double*, int*) const;
  void sample_v_given_h(const int*,const vector<int>&, double*, int*) const;
  double propup(const int*,const vector<int> &imovie, double*, double) const;
  void propdown(const int*,const int,const  double*,const  double*,const double*,const double*, double*) const;
  void gibbs_hvh(const int*,const vector<int>&, double*, int*, double*, int*) const;
  void reconstruct(const vector<int>&, const vector<int>&, vector<double>&) const;
  void reconstruct(const vector<int>&, const vector<int>&, const vector<int>&, vector<double>&) const;
  void reconstruct(const vector<int>&, const vector<int>&, const vector<int>&, const vector<int>&, vector<double>&) const;
  void output();


};

