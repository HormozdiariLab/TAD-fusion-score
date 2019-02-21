//#: Date         :05/05/2017
//#: Author       :Linh Huynh
//#: Version      :1.0.0 

#include <ilcplex/cplex.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

using namespace std;

typedef vector<int> IdList;
typedef vector<double> ValueList;
typedef vector<IdList> IdMatrix;
typedef vector<ValueList> ValueMatrix; 
typedef vector<string> StringList;
typedef vector<StringList> StringMatrix;

typedef map<string,int> IdMap;

#define STR_EQ        0
#define ROW_DELIM     '\n'
#define FIELD_DELIM   '\t'
#define COMMA_DELIM   ','

#define ZERO          1e-1

struct HiCDataModel {
  ValueList alpha, beta_1, beta_2, insulation;
};

struct ModelFittingParameter {
  int resolution, start, end;
  int min_length_in_obj, mid_length_in_obj, max_length_in_obj;
  double min_beta;
  string fitting_method;   // full or segmentation
  int seg_length, min_seg_overlap;
  double zero_const;
  int print_debug_info;
  string output_filename, output_comment;
  ModelFittingParameter() {
    resolution = start = end = -1;
    min_length_in_obj = mid_length_in_obj = max_length_in_obj = -1;
    min_beta = -2;
    fitting_method = "";
    seg_length = min_seg_overlap = -1;
    zero_const = 1e-1;
    print_debug_info = 0;
    output_filename = output_comment = "";
  }
};

struct HicFileInfo {
  string filename, file_format;
  int header_line_num;
  HicFileInfo() {
    filename = file_format = "";
    header_line_num = 0;
  }
};

// Auxiliary functions
string num_to_string (int);
double max(double, double);
double min(double, double);
StringList split(const string&s, char delim);
// Common data structure functions
void resize_and_fill(ValueList&, int size, int value);                  // size, value
void resize_and_fill(ValueMatrix&, int size, int value);                // Squre matrix only
void read_a_table_file (string filename, char row_delim_char, char field_delim_char, 
	char comment_char, int skipping_header_line_num, 
	const IdList& col_list, StringMatrix& str_mat);                 // index of 1st col = 0
void read_a_value_list (string filename, int col, ValueList& val_list); // index of 1st col = 0
void write_a_value_list (string filename, ValueList&);
void sort_by_value (const ValueList& l, ValueList& sorted_index, ValueList& rank, bool is_asc);
// Read HiC data
void read_full_hic_data_file (string filename, int header_row_num, int header_col_num, ValueMatrix& hi_c_matrix); // Full matrix format (Bing Ren data)
void read_sparse_hic_data_file (string filename, int resolution, ValueMatrix& hi_c_matrix);  // 3-col format (Rao data)
void read_5C_file (string filename, int header_line_num, int matrix_size,  ValueMatrix& hi_c_matrix);
int read_Mundlos_file(string filename, int& resolution, int& min_start, int& max_end, ValueMatrix& hi_c_matrix);
// Fit the model without partitioning
void full_fitting_hic_data_model (const ValueMatrix& hic_mat, const ModelFittingParameter& model_fitting_parameter, HiCDataModel& hic_data_model);
// Fit the model with/without partitioning and print the model to file
void fit_hic_data_model (const ValueMatrix& hic_mat, const ModelFittingParameter& model_fitting_parameter);

// Parse the parameter
void parse_parameter(int argc, char** argv, HicFileInfo& hic_file_info, ModelFittingParameter& model_fitting_parameter) {
  if (argc % 2 == 0) {
    cout << "ERROR: Number of parameters must be even" << endl;
    return;
  }
  string summary_str = "";
  for (int i = 1; i < argc; i = i + 2) {
    string option(argv[i]);
    cout << option << "\t" << argv[i+1] << endl;
    bool is_summary_added = true;
    if (option.compare("-fn") == STR_EQ) {
      hic_file_info.filename = argv[i+1];
      is_summary_added = false;
    }
    else if (option.compare("-ff") == STR_EQ) {
      hic_file_info.file_format = argv[i+1];
      is_summary_added = false;
    }
    else if (option.compare("-fhl") == STR_EQ) {
      hic_file_info.header_line_num = atoi(argv[i+1]);
      is_summary_added = false;
    }
    else if (option.compare("-res") == STR_EQ) {
      model_fitting_parameter.resolution = atoi(argv[i+1]);
      is_summary_added = false;
    }
    else if (option.compare("-start") == STR_EQ) {
      model_fitting_parameter.start = atoi(argv[i+1]);
      is_summary_added = false;
    }
    else if (option.compare("-end") == STR_EQ) {
      model_fitting_parameter.end = atoi(argv[i+1]);
      is_summary_added = false;
    }
    else if (option.compare("-mn") == STR_EQ) 
      model_fitting_parameter.min_length_in_obj = atoi(argv[i+1]);
    else if (option.compare("-md") == STR_EQ)
      model_fitting_parameter.mid_length_in_obj = atoi(argv[i+1]);
    else if (option.compare("-mx") == STR_EQ)
      model_fitting_parameter.max_length_in_obj = atoi(argv[i+1]);
    else if (option.compare("-mbt") == STR_EQ)
      model_fitting_parameter.min_beta = atof(argv[i+1]);
    else if (option.compare("-method") == STR_EQ)
      model_fitting_parameter.fitting_method = argv[i+1];
    else if (option.compare("-sg") == STR_EQ)
      model_fitting_parameter.seg_length = atoi(argv[i+1]);
    else if (option.compare("-mso") == STR_EQ)
      model_fitting_parameter.min_seg_overlap = atoi(argv[i+1]);
    else if (option.compare("-zero") == STR_EQ)
      model_fitting_parameter.zero_const = atof(argv[i+1]);
    else if (option.compare("-debug") == STR_EQ) {
      model_fitting_parameter.print_debug_info = atoi(argv[i+1]);
      is_summary_added = false;
    }
    else if (option.compare("-cm") == STR_EQ) {
      model_fitting_parameter.output_comment = argv[i+1];
      is_summary_added = false;
    }
    else if (option.compare("-of") == STR_EQ) {
      model_fitting_parameter.output_filename = argv[i+1];
      is_summary_added = false;
    }
    else {
      cout << "ERROR: Can not recognize the option " << option << endl;
      return;
    }
    if (is_summary_added)
      summary_str += ("_" + option + ":" + (::string(argv[i+1])));
  }
  model_fitting_parameter.output_comment += summary_str;
}


int main (int argc, char** argv) {
  HicFileInfo hic_file_info;
  ModelFittingParameter model_fitting_parameter;
  parse_parameter(argc, argv, hic_file_info, model_fitting_parameter);
  // Read Hi-C data
  ValueMatrix hi_c_mat;
  if (hic_file_info.file_format.compare("full_matrix_format") == STR_EQ)
    read_full_hic_data_file (hic_file_info.filename, hic_file_info.header_line_num, 0, hi_c_mat);
  else if (hic_file_info.file_format.compare("sparse_matrix_format") == STR_EQ)
    read_sparse_hic_data_file (hic_file_info.filename, model_fitting_parameter.resolution, hi_c_mat);
  else {
    cout << "Can not find the format " << hic_file_info.file_format << endl;
    return -1;
  }
  fit_hic_data_model(hi_c_mat, model_fitting_parameter);
  cout << "Complete!" << endl;
}

// Only interactions that have the length <= max_length_in_obj will contribute to the obj function
//void fit_hi_c_data_model (const ValueMatrix& hi_c_mat, int left, int right, int min_length_in_obj, int mid_length_in_obj, int max_length_in_obj, double zero_weight, HiCDataModel& hi_c_data_model, bool print_solver_status) {
void full_fitting_hic_data_model (const ValueMatrix& hic_mat, const ModelFittingParameter& model_fitting_parameter, HiCDataModel& hic_data_model) {
  int start = model_fitting_parameter.start,
      end = model_fitting_parameter.end,
      min_length_in_obj = model_fitting_parameter.min_length_in_obj,
      mid_length_in_obj = model_fitting_parameter.mid_length_in_obj,
      max_length_in_obj = model_fitting_parameter.max_length_in_obj,
      print_solver_status = model_fitting_parameter.print_debug_info;
  double zero_const = model_fitting_parameter.zero_const;
  
  cout << hic_mat.size() << "\t" << start << "\t" << end << "\t" << zero_const << endl;
  CPXENVptr env = NULL;
  CPXLPptr lp = NULL;

  int length = end - start + 1;
  int alpha_size = length;
  int beta_size = 2;
  int e_size = ((length > max_length_in_obj)? 
                 ((length - max_length_in_obj)*max_length_in_obj + max_length_in_obj*(max_length_in_obj - 1)/2):
                 (length*(length - 1)/2));	// slack
      e_size = e_size - (length - min_length_in_obj + 1)*(min_length_in_obj - 1) - (min_length_in_obj - 1)*(min_length_in_obj - 2)/2;
  int r_size = length;	// insulator

  int var_num = alpha_size + beta_size + e_size + r_size; 
  int status = 0;
  env = CPXopenCPLEX (&status);
  // Turn on output to the screen 
  status = CPXsetintparam (env, CPXPARAM_ScreenOutput, (print_solver_status? CPX_ON:CPX_OFF));
  // Turn on data checking 
  status = CPXsetintparam (env, CPXPARAM_Read_DataCheck, CPX_ON);
  //status = CPXsetintparam (env, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
  // Create the problem 
  lp = CPXcreateprob (env, &status, "linear_model");
  // Populate the problem
  // Problem is minimization
  status = CPXchgobjsen (env, lp, CPX_MIN);
  // For obj and bounds
  double* obj = (double*) malloc(var_num*sizeof(double));
  double* lb = (double*) malloc(var_num*sizeof(double));
  double* ub = (double*) malloc(var_num*sizeof(double));
  char** colname = (char**) malloc(var_num*sizeof(char*));
	
  int  current_id = 0;
  // alpha
  for (int i = 0; i < alpha_size; i++) {
    lb[current_id] = -CPX_INFBOUND;
    ub[current_id] = CPX_INFBOUND;
    colname[current_id] = (char*) malloc(50*sizeof(char));
    sprintf(colname[current_id], "a_%i\0", i);
    obj[current_id] = 0;
    current_id++;
  }

  for (int i = 0; i < beta_size; i++) {
    lb[current_id] = model_fitting_parameter.min_beta;
    ub[current_id] = -0.05;
    colname[current_id] = (char*) malloc(50*sizeof(char));
    sprintf(colname[current_id], "beta_%i\0", i);
    obj[current_id] = 0;
    current_id++;
  }
  // slack
  int tmp = current_id;
  int cell_num_tmp = 0;
  for (int i = 0; i < length; i++) {
    int max_i = ((i + max_length_in_obj < length)? (i + max_length_in_obj):(length-1));
    for (int j = i + min_length_in_obj; j <= max_i; j++) {
      //for (int j = i + 1; j <= max_i; j++) {
      lb[current_id] = 0;
      ub[current_id] = CPX_INFBOUND;
      colname[current_id] = (char*) malloc(50*sizeof(char));
      sprintf(colname[current_id], "e_%i,%i\0", i, j);
      obj[current_id] = 1;
      current_id++;
      cell_num_tmp += (4 + j - i);
    }
  }
  //cout << "#Var " << current_id - tmp << "\t" << e_size << endl;
  // insulator
  for (int i = 0; i < length; i++) {
    lb[current_id] = 0;
    ub[current_id] = CPX_INFBOUND;
    colname[current_id] = (char*) malloc(50*sizeof(char));
    sprintf(colname[current_id], "r_%i\0", i);
    obj[current_id] = 0;
    current_id++;
  }
  // Obj and the bound
  status = CPXnewcols (env, lp, var_num, obj, lb, ub, NULL, colname);
  if (status) {
    cout << "LINH: Fail to set up (bound and obj) the LP problem!" << endl;
    cout << "(start,end) = " << start << "\t" << end << endl;
  }
  //else 
  //	cout << "LINH: Add all variables successfully!" << endl;
  // Add constraints
  int row_num = 2*e_size;
  int cell_num = 2*cell_num_tmp;
  int* rmatbeg = (int*) malloc(row_num*sizeof(int));
  int* rmatind = (int*) malloc(cell_num*sizeof(int));
  double* rmatval = (double*) malloc(cell_num*sizeof(double));
  double* rhs = (double*) malloc(row_num*sizeof(double));
  char* sense = (char*) malloc(row_num*sizeof(char));
  char** rowname = (char**) malloc(row_num*sizeof(char*));
	
  int current_row_id = 0;
  int current_cell_id = 0;
  //cout << "Total: " << length << endl;
  for (int i = 0; i < length; i++) {
    int j_max = ((i + max_length_in_obj < length)? (i + max_length_in_obj):(length-1));
    //if (i > n - 10)
    //	cout << i << endl;
    //for (int j = i + 1; j <= j_max; j++) {
    for (int j = i + min_length_in_obj; j <= j_max; j++) {
      double logH = ((hic_mat[start + i][start + j] > 0)? log(hic_mat[start + i][start + j]) : log(zero_const));
      // (a_i + a_j)/2 + beta*log(d_ij) - sum(r_k) - e_ij <= log(H_ij)
      rmatbeg[current_row_id] = current_cell_id;
      sense[current_row_id] = 'L';
      rhs[current_row_id] = logH;
      rowname[current_row_id] = (char*) malloc(20*sizeof(char));
      sprintf(rowname[current_row_id], "ce1_%i,%i\0", i, j);
      rmatind[current_cell_id] = i;
      rmatval[current_cell_id] = 0.5;
      rmatind[current_cell_id + 1] = j;
      rmatval[current_cell_id + 1] = 0.5;
      rmatind[current_cell_id + 2] = ((j - i <= mid_length_in_obj)? alpha_size:(alpha_size + 1));
      rmatval[current_cell_id + 2] = log(j - i);
      rmatind[current_cell_id + 3] = alpha_size + beta_size + round(current_row_id/2);
      rmatval[current_cell_id + 3] = -1;
      for (int k = i + 1; k <= j; k++) {
        rmatind[current_cell_id + 3 + k - i] = alpha_size + beta_size + e_size + k;
        rmatval[current_cell_id + 3 + k - i] = -1; 
      }
      current_row_id++;
      current_cell_id += 4 + (j-i);   // 4: a_i, a_j, beta, e_ij; (j-i): sum(r_k) 
      // -(a_i + a_j)/2 - beta*log(d_ij) + sum(r_k) - e_ij <= -log(H_ij)
      rmatbeg[current_row_id] = current_cell_id;
      sense[current_row_id] = 'L';
      rhs[current_row_id] = -logH;
      rowname[current_row_id] = (char*) malloc(20*sizeof(char));
      sprintf(rowname[current_row_id], "ce2_%i,%i\0", i, j);
      rmatind[current_cell_id] = i;
      rmatval[current_cell_id] = -0.5;
      rmatind[current_cell_id + 1] = j;
      rmatval[current_cell_id + 1] = -0.5;
      rmatind[current_cell_id + 2] = ((j - i <= mid_length_in_obj)? alpha_size:(alpha_size + 1));
      rmatval[current_cell_id + 2] = -log(j - i);
      rmatind[current_cell_id + 3] = alpha_size + beta_size + round(current_row_id/2);
      rmatval[current_cell_id + 3] = -1;
      for (int k = i + 1; k <= j; k++) {
        rmatind[current_cell_id + 3 + k - i] = alpha_size + beta_size + e_size + k;
        rmatval[current_cell_id + 3 + k - i] = 1; 
      }
      current_row_id++;
      current_cell_id += 4 + (j-i);
    }
  }
  //cout << "Estimate: " << row_num << "\t" << cell_num << endl;
  //cout << "Iterate: " << current_row_id << "\t" << current_cell_id << endl;

  //status = CPXaddrows (env, lp, 0, row_num, cell_num, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);
  status = CPXaddrows (env, lp, 0, current_row_id, current_cell_id, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);
  if (status) {
    cout << "LINH: Fail to add constraints to LP problem!" << endl;
    //for (int i = 0; i < var_num; i++)
    //	cout << i << "\t" << lb[i] << "\t" << ub[i] << "\t" << obj[i] << "\t" << colname[i] << endl;
    //cout << "=============================" << endl;
    //for (int i = 0; i < current_row_id; i++)
    //	cout << i << "\t" << "rmatbeg: " << rmatbeg[i] << "\tsense: " << sense[i] << "\trhs: " << rhs[i] << "\t" << rowname[i] << endl;
    //for (int i = 0; i < current_cell_id; i++)
    //	cout << i << "\t" << "rmatind: " << rmatind[i] << "\t rmatval: " << rmatval[i] << endl;		
  }
  //else {
  //	cout << "LINH: Add all constraints successfully!" << endl;
  //}	
  // Solve the problem
  time_t t_start, t_end;
  t_start = clock();
  status = CPXlpopt (env, lp); 

  t_end = clock();
  double total_time = (double)(t_end - t_start)/(double)CLOCKS_PER_SEC;
  //if (print_solver_status)
  printf( "\t\t\tElapsed time : %0.3f \n", total_time );


  // Get the actual size of the problem
  int cur_numrows = CPXgetnumrows (env, lp);
  int cur_numcols = CPXgetnumcols (env, lp);
  double* x = (double*) malloc (cur_numcols*sizeof(double));
  double* slack = (double*) malloc (cur_numrows*sizeof(double));
  double* dj = (double*) malloc (cur_numcols*sizeof(double));
  double* pi = (double*) malloc (cur_numrows*sizeof(double));
  int solstat;
  double objval;
  // Get the solution
  //status = CPXwriteprob (env, lp, "linear_model.lp", NULL);
  status = CPXsolution (env, lp, &solstat, &objval, x, pi, slack, dj);
  if (status)
    cout << "LINH: Fail to solve the LP problem!" << endl;
  else {
    //if (print_solver_status)
    cout << "LINH: # Error term = " << e_size << ", final obj value = " << objval << endl;
    //for (int i = 0; i < var_num; i++)
    //	cout << colname[i] << "\t" << x[i] << endl;
  }
	
  // Export the parameter values
  //hic_data_model.start = start;
  //hic_data_model.end = end;
  cout << x[alpha_size] << "\t" << x[alpha_size + 1] << endl;
  resize_and_fill(hic_data_model.alpha, length, 0);
  resize_and_fill(hic_data_model.beta_1, length, 0);
  resize_and_fill(hic_data_model.beta_2, length, 0);
  resize_and_fill(hic_data_model.insulation, length, 0);
  for (int i = 0; i < length; i++) {
    hic_data_model.alpha[i] = x[i];
    hic_data_model.beta_1[i] = x[alpha_size];
    hic_data_model.beta_2[i] = x[alpha_size + 1];
    hic_data_model.insulation[i] = x[alpha_size + beta_size + e_size + i];
  }

  // Free the problem
  CPXfreeprob(env, &lp);
  // Free the environment
  CPXcloseCPLEX (&env);
  // Free all memory
  free(obj);
  free(ub);
  free(lb);
  for (int i = 0 ; i < var_num; i++)
    free(colname[i]);
  free(colname);
  free(x);
  free(slack);
  free(dj);
  free(pi);
  free(rmatbeg);
  free(rmatind);
  free(rmatval);
  free(rhs);
  free(sense);
  for (int i = 0; i < current_row_id; i++)
    free(rowname[i]);
  free(rowname);	
}

void fit_hic_data_model (const ValueMatrix& hic_mat, const ModelFittingParameter& model_fitting_parameter) {
  int start = model_fitting_parameter.start,
      end = model_fitting_parameter.end,
      seg_length = model_fitting_parameter.seg_length,
      min_seg_overlap = model_fitting_parameter.min_seg_overlap;
 
  if (start < 1)
    start = 0;
  if (end < 1)
    end = hic_mat.size() - 1;
  if (model_fitting_parameter.fitting_method.compare("full") == STR_EQ)
    seg_length = hic_mat.size();

  ValueList alpha, beta_1, beta_2, insulation;
  // Segment and fit
  int current_start = start;
  while (current_start <= end) {
    int current_end = ((current_start + seg_length >= end - 0.25*seg_length)? end:(current_start + seg_length)),
        current_length = current_end - current_start + 1;
    ModelFittingParameter para_tmp = model_fitting_parameter;
    para_tmp.start = current_start;
    para_tmp.end = current_end;
    HiCDataModel model_tmp;
    cout << "(start, end, current_start, current_end) = " << start << "\t" << end << "\t" << current_start << "\t" << current_end << endl;    
    full_fitting_hic_data_model(hic_mat, para_tmp, model_tmp);
    int break_point_pos = current_length;
    if (current_end < end) {
      double epsilon = 1e-10;
      int offset = current_length - min_seg_overlap;
      break_point_pos = -1;
      for (int i = offset; i >= 2; i--) {
        if (model_tmp.insulation[i] < epsilon && model_tmp.insulation[i - 1] < epsilon && model_tmp.insulation[i - 2] < epsilon) {
          break_point_pos = i - 1;
          break;
        }
      }
      if (break_point_pos < 0) {
        cout << "ERROR: Can not find the break point to segment the chromosome for fitting" << endl;
        break_point_pos = offset - 1;
      }
    }
    for (int i = 0; i < break_point_pos; i++) {
      alpha.push_back(model_tmp.alpha[i]);
      beta_1.push_back(model_tmp.beta_1[i]);
      beta_2.push_back(model_tmp.beta_2[i]);
      insulation.push_back(model_tmp.insulation[i]);
    }
    cout << "Break point " << break_point_pos << endl;
    current_start = current_start + break_point_pos;
  }
  // Print the model to file
  cout << "(alpha, beta_1, beta_2, insultion) = " 
       << alpha.size() << "\t" << beta_1.size() << "\t" << beta_2.size() << "\t" << insulation.size() << endl;
  ofstream out_file(model_fitting_parameter.output_filename.c_str());
  if (!out_file.is_open()) {
    cout << "ERROR: Can not open the file " << model_fitting_parameter.output_filename << endl;
    return;
  }
  out_file << "#" << model_fitting_parameter.output_comment;
  out_file << endl << model_fitting_parameter.resolution << "\t" << start << "\t" << end << "\t" << model_fitting_parameter.mid_length_in_obj;
  for (int i = 0; i < alpha.size(); i++)
    out_file << endl << alpha[i] << "\t" << beta_1[i] << "\t" << beta_2[i] << "\t" << insulation[i];
  out_file.close();
}

int read_Mundlos_file(string filename, int& resolution, int& min_start, int& max_end, ValueMatrix& hi_c_matrix) {
  IdList col_list = {0,1,2};
  StringMatrix str_mat;
  read_a_table_file (filename, ROW_DELIM, FIELD_DELIM, '#', 0, col_list, str_mat);
  string chr_name = "";
  resolution = min_start = max_end = -1;
  for (int pass = 1; pass <= 2; pass++) {
    if (pass == 2) {
      if ((max_end - min_start + 1) % resolution != 0) {
        cout << "(start, end, resolution) = (" << min_start << "," << max_end << "," << resolution << "), length is not a multiplier of the resolution!" << endl;
        return -1;
      }
      else {
        int matrix_size = (max_end - min_start + 1)/resolution;
        resize_and_fill(hi_c_matrix, matrix_size, 0);
      }
    }
    for (int i = 0; i < str_mat.size(); i++) {
      StringList tmp1 = split(str_mat[i][0], ','),
                 tmp2 = split(str_mat[i][1], ':'),
                 tmp3 = split(tmp2[1], '-');

      string chr1 = tmp1[0], chr2 = tmp2[0];
      int start_1 = atoi(tmp1[1].c_str()), end_1 = atoi(tmp1[2].c_str()),
          resolution_1 = end_1 - start_1 + 1,
          start_2 = atoi(tmp3[0].c_str()), end_2 = atoi(tmp3[1].c_str()),
          resolution_2 = end_2 - start_2 + 1;

      if (pass == 1) {
        // First bin
        if (chr_name.empty()) 
          chr_name = chr1;
        else if (chr_name.compare(chr1) != STR_EQ) {
          cout << "Line " << i + 1 << ", 1st chromosome name is different! " << chr1 << endl;
          return -1;
        }   
	if (resolution < 0) 
          resolution = resolution_1;
        else if (resolution != resolution_1) {
          cout << "Line " << i + 1 << ", 1st resolution is different! " << resolution_1 << endl;
          return -1;
        }
        if (min_start < 0 || min_start > start_1)
          min_start = start_1;
        if (max_end < 0 || max_end < end_1)
          max_end = end_1;

        // Second bin
        if (chr_name.empty()) 
          chr_name = chr2;
        else if (chr_name.compare(chr2) != STR_EQ) {
          cout << "Line " << i + 1 << ", 2nd chromosome name is different! " << chr2 << endl;
          return -1;
        }        
	if (resolution < 0) 
          resolution = resolution_2;
        else if (resolution != resolution_2) {
          cout << "Line " << i + 1 << ", 2nd resolution is different! " << resolution_2 << endl;
          return -1;
        }
        if (min_start < 0 || min_start > start_2)
          min_start = start_2;
        if (max_end < 0 || max_end < end_2)
          max_end = end_2;
      }
      else {
        int row = (start_1 - min_start)/resolution,
            col = (start_2 - min_start)/resolution;
        hi_c_matrix[row][col] = hi_c_matrix[col][row] = atof(str_mat[i][2].c_str());
      }
    }
  }
  //for (int i = 0; i < hi_c_matrix.size(); i++)
  //  for (int j = 0; j < i; j++)
  //    if (hi_c_matrix[j][i] < 0.1)
  //      cout << j << "\t" << i << "\t" << hi_c_matrix[j][i] << endl;
}

void read_sparse_hic_data_file (string filename, int resolution, ValueMatrix& hi_c_matrix) {
  IdList src, dest;
  ValueList val_list;
  int pos_max = -1;

  ifstream in_file(filename.c_str());
  if (in_file.is_open()) {		
    int line_num = 0;
    for (string row; getline(in_file, row, ROW_DELIM); ) {
      line_num++;
      istringstream ss(row);
      StringList record;
      for (string field; getline(ss, field, FIELD_DELIM); )
      record.push_back(field);
      if (!record.empty() && !record[0].empty() && record[0][0] != '#') {
        if (record.size() >= 3) {
          int x = atoi(record[0].c_str()), y = atoi(record[1].c_str());
          if (pos_max < x)
            pos_max = x;
          if (pos_max < y)
            pos_max = y;
          src.push_back(x);
          dest.push_back(y);
          val_list.push_back(atof(record[2].c_str()));				
        }
      }
    }
    cout << "File " << filename << ", #lines " << line_num << endl;
    in_file.close();
  }
  else {
    cout << "Can not open file " << filename << endl;
    return;
  }
  int size  = floor(pos_max/resolution) + 1;
  cout << "Matrix size " << size << endl;
  hi_c_matrix.resize(size);
  for (int i = 0; i < size; i++) {
    hi_c_matrix[i].resize(size);
    for (int j = 0; j < size; j++)
      hi_c_matrix[i][j] = 0;
  }
  for (int i = 0; i < src.size(); i++) {
    int x = floor(src[i]/resolution);
    int y = floor(dest[i]/resolution);
    if (hi_c_matrix[x][y] > 0 || hi_c_matrix[y][x] > 0)
      cout << "ERROR: matrix at " << x << "\t" << y << " has value " << hi_c_matrix[x][y] << endl;
    else {			
      hi_c_matrix[x][y] = val_list[i];
      hi_c_matrix[y][x] = val_list[i];
    }
  }
}

void read_full_hic_data_file (string filename, int header_row_num, int header_col_num, ValueMatrix& hi_c_matrix) {
  ifstream in_file(filename.c_str());
  if (in_file.is_open()) {		
    int line_num = 0;
  for (string row; getline(in_file, row, ROW_DELIM); ) {
    line_num++;
    if (line_num > header_row_num) {
      istringstream ss(row);
      StringList record;
      for (string field; getline(ss, field, FIELD_DELIM); )
      record.push_back(field);
      if (line_num == header_row_num + 1) 
        hi_c_matrix.resize(record.size() - header_col_num);
      if (record.size() != hi_c_matrix.size() + header_col_num) {
        cout << "ERROR: line " << line_num << ", #col = " << record.size() << ", the matrix size = " << hi_c_matrix.size() 
             << ", #header col = " << header_col_num << endl;
        return;
      }
      hi_c_matrix[line_num - header_row_num - 1].clear();
      for (int i = 0; i < hi_c_matrix.size(); i++)
        hi_c_matrix[line_num - header_row_num - 1].push_back(atoi(record[i + header_col_num].c_str()));				
    }
  }
  cout << "File " << filename << ", #lines " << line_num << endl;
  if (line_num != header_row_num + hi_c_matrix.size())
    cout << "ERROR: #header row = " << header_row_num << ", matrix size = " << hi_c_matrix.size() 
         << " with " << line_num << " lines!" << endl;
    in_file.close();
  }
  else {
    cout << "Can not open file " << filename << endl;
    return;
  }
}

void read_5C_file (string filename, int header_line_num, int matrix_size,  ValueMatrix& hi_c_matrix) {
  hi_c_matrix.resize(matrix_size);

  ifstream in_file(filename.c_str());
  if (in_file.is_open()) {		
    int line_num = 0;
    for (string row; getline(in_file, row, ROW_DELIM); ) {
      line_num++;
      istringstream ss(row);
      StringList record;
      for (string field; getline(ss, field, FIELD_DELIM); )
      record.push_back(field);
      if (line_num > header_line_num) {
        if (record.size() != matrix_size + 1) {
          cout << "ERROR: line " << line_num << " has " << record.size() << " elements while the matrix size = " << matrix_size << endl; 
          return;
        }
        hi_c_matrix[line_num - header_line_num - 1].clear();
        for (int i = 0; i < matrix_size; i++)
          hi_c_matrix[line_num - header_line_num - 1].push_back(atoi(record[i + 1].c_str()));			
      }
    }
    cout << "File " << filename << ", #lines " << line_num << endl;
    if (line_num != header_line_num + matrix_size)
      cout << "ERROR: Reading a 5C file with " << header_line_num << " header lines, matrix size = " << matrix_size 
           << " with " << line_num << " lines!" << endl;
    in_file.close();
  }
  else {
    cout << "Can not open file " << filename << endl;
    return;
  }
}

void read_a_table_file (string filename, char row_delim_char, char field_delim_char, char comment_char, int skipping_header_line_num, const IdList& col_list, StringMatrix& str_mat) {
  str_mat.clear();
  ifstream in_file(filename.c_str());
  if (in_file.is_open()) {		
    cout << "Reading file " << filename << endl;
    int line_num = 0;
    for (string row; getline(in_file, row, row_delim_char); ) {
      line_num++;
      if (line_num > skipping_header_line_num) {
        istringstream ss(row);
        StringList record;
        for (string field; getline(ss, field, field_delim_char); )
        record.push_back(field);
        if (!record.empty() && !record[0].empty() && record[0][0] != comment_char) {
          StringList tmp;
          for (int i = 0; i < col_list.size(); i++)
            if (record.size() > col_list[i])
              tmp.push_back(record[col_list[i]]);
            else
              cout << "\tline " << line_num << " has only " << record.size() << " fields while the required index is " << col_list[i] << endl;
          if (tmp.size() == col_list.size())
              str_mat.push_back(tmp);
        }
      }
    }
    cout << "\t#lines " << line_num << "\t#records " << str_mat.size() << endl;
    in_file.close();	
  }
  else {
    cout << "Can not open file " << filename << endl;
    return;
  }
}

void read_a_value_list (string filename, int col, ValueList& val_list) {
  IdList col_list(1, col);
  StringMatrix str_mat;
  read_a_table_file (filename, ROW_DELIM, FIELD_DELIM, '#', false, col_list, str_mat);
  val_list.clear();
  for (int i = 0; i < str_mat.size(); i++)
    val_list.push_back(strtod(str_mat[i][0].c_str(), NULL));
}

void write_a_value_list (string filename, ValueList& val_list) {
  ofstream out_file(filename.c_str());
  for (int i = 0; i < val_list.size(); i++)
    out_file << i << "\t" << val_list[i] << endl;
  out_file.close();
}

void sort_by_value (const ValueList& l, ValueList& sorted_index, ValueList& rank, bool is_asc) {
  ValueList x;
  for (int i = 0; i < l.size(); i++)
    x.push_back(l[i]*(is_asc? 1:-1));
  int n = x.size();
  sorted_index.resize(n);
  std::size_t tmp(0);
  std::generate(std::begin(sorted_index), std::end(sorted_index), [&]{ return tmp++; });
  std::sort(std::begin(sorted_index), std::end(sorted_index), [&](int i, int j) { // Sort min -> max
    return (x[i] < x[j]);
  });
  rank.resize(n);
  for (int i = 0; i < n; i++)
    rank[sorted_index[i]] = i + 1;
}

StringList split(const string&s, char delim) {
  stringstream ss(s);
  string item;
  StringList elems;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

string num_to_string (int n) {
  stringstream ss;
  ss << n;
  return ss.str();
}

double max(double a, double b) {
  return (a > b)? a:b;
}

double min(double a, double b) {
  return (a < b)? a:b;
}

void resize_and_fill(ValueList& x, int size, int val) {
  if (x.size() != size)
    x.resize(size);
  for (int i = 0; i < size; i++)
    x[i] = val;
}
void resize_and_fill(ValueMatrix& x, int size, int val) {
  if (x.size() != size)
    x.resize(size);
  for (int i = 0; i < size; i++)
    resize_and_fill(x[i], size, val);
}

