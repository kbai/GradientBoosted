class RBM {

public:
  int N;
  int n_visible;
  int n_hidden;
  double **W;
  double **dW,**sW;
  double *hbias;
  double *dh,*sh;
  double *vbias;
  double *dv,*sv;
  RBM(int, int, int, double**, double*, double*);
  ~RBM();
  void contrastive_divergence(vector<int>& ,vector<int>&, double, int);
  void mpigath();
  void sample_h_given_v(int*,vector<int>&, double*, int*);
  void sample_v_given_h(int*,vector<int>&, double*, int*);
  double propup(int*,vector<int> &imovie, double*, double);
  void propdown(int*, int, double*, double*);
  void gibbs_hvh(int*,vector<int>&, double*, int*, double*, int*);
  void reconstruct(vector<int>&, vector<int>&, vector<double>&);
  void reconstruct(vector<int>&, vector<int>&, vector<int>&, vector<double>&);
  void reconstruct(vector<int>&, vector<int>&, vector<int>&,vector<int>&, vector<double>&);


};

