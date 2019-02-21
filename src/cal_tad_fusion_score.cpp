//#: Date         :05/05/2017
//#: Author       :Linh Huynh
//#: Version      :1.0.0 

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
                        const IdList& col_list, StringMatrix& str_mat); // index of 1st col = 0
void read_a_value_list (string filename, int col, ValueList& val_list); // index of 1st col = 0
void write_a_value_list (string filename, ValueList&);
void sort_by_value (const ValueList& l, ValueList& sorted_index, ValueList& rank, bool is_asc);

double estimate_fusion_score(const ValueList& alpha, const ValueList& beta, const ValueList& insulation, int del_start_bin, int del_end_bin, int window_size, double delta, bool debug) {
  int start_window = ((del_start_bin > window_size)? (del_start_bin - window_size):0),
      end_window = ((del_end_bin + window_size < (alpha.size() - 1))? (del_end_bin + window_size):(alpha.size() - 1));
  //cout << chr_length << "\t" << del_start_bin << "\t" << del_end_bin << endl;
  if (del_end_bin <= del_start_bin) {
    cout << "ERROR: Too short deletion, del_end_bin <= del_start_bin at " << del_start_bin << "\t" << del_end_bin << endl;
    return -1;
  }
  else {
    double tad_fusion_score = 0, del_insulation = 0;
    for (int i = del_start_bin; i <= del_end_bin; i++)
      del_insulation += insulation[i];
    double insulation_tmp_i = del_insulation;
    for (int i = del_start_bin - 1; i >= start_window; i--) {
      if (i < del_start_bin - 1)
        insulation_tmp_i += insulation[i+1];
      double all_insulation = insulation_tmp_i;
      for (int j = del_end_bin + 1; j <= end_window; j++) {
        all_insulation += insulation[j];
        double before_del_val = exp(((beta[i] + beta[j])/2)*log(j - i) - all_insulation),
               after_del_val = exp(((beta[i] + beta[j])/2)*log(j - i - (del_end_bin - del_start_bin + 1)) - all_insulation + del_insulation);
        if (before_del_val < delta && after_del_val > delta && alpha[i] >= 1 && alpha[j] >= 1)
          tad_fusion_score++;
        if (debug) {
          cout << i << "\t" << j << "\t"
               << alpha[i] << "\t" << alpha[j] << "\t"
               << beta[i] << "\t" << beta[j] << "\t"
               << all_insulation << "\t" << del_insulation << "\t"
               << before_del_val << "\t" << after_del_val << endl;
        }
      }
      //break;      
    }
    return tad_fusion_score;
  } 
}


int main (int argc, char** argv) {
  string del_filename, model_dir, out_filename;
  int chr_min = 1, chr_max = 23, min_del_length = 0, max_del_length = 1e7;
  int debug = 0, window_size, permutation_num = 0;
  double delta;

  // Parse the parameters
  if (argc % 2 == 0) {
    cout << "ERROR: Number of parameters must be even" << endl;
    return -1;
  }
  for (int i = 1; i < argc; i = i + 2) {
    string option(argv[i]);
    //cout << option << "\t" << argv[i+1] << endl;
    if (option.compare("-f") == STR_EQ)
      del_filename = argv[i+1];
    else if (option.compare("-md") == STR_EQ)
      model_dir = argv[i+1];
    else if (option.compare("-w") == STR_EQ)
      window_size = atoi(argv[i+1]);
    else if (option.compare("-d") == STR_EQ)
      delta = atof(argv[i+1]);
    else if (option.compare("-mnl") == STR_EQ)
      min_del_length = atoi(argv[i+1]);
    else if (option.compare("-mxl") == STR_EQ)
      max_del_length = atof(argv[i+1]);
    else if (option.compare("-chr") == STR_EQ)
      chr_min = chr_max = atoi(argv[i+1]);
    else if (option.compare("-per") == STR_EQ)
      permutation_num = atoi(argv[i+1]);
    else if (option.compare("-debug") == STR_EQ)
      debug = atoi(argv[i+1]);
    else if (option.compare("-o") == STR_EQ)
      out_filename = argv[i+1];
    else {
      cout << "ERROR: Can not recognize the option " << option << endl;
      return -1;
    }
  }
  // Ready for randomizing the permutation
  srand (time(NULL));
  // Read the deletion file
  StringMatrix del_str_mat;
  read_a_table_file (del_filename, ROW_DELIM, FIELD_DELIM, '#', 0, {0,1,2}, del_str_mat);
  ValueList tad_fusion_score_list(del_str_mat.size(), -1);
  ValueMatrix permutation_score_mat;
  for (int p = 0; p < permutation_num; p++)
    permutation_score_mat.push_back(tad_fusion_score_list);
 
  for (int chr = chr_min; chr <= chr_max; chr++) {
    string chr_str = "chr" + ((chr < 23)? num_to_string(chr):"X");
    cout << chr_str << endl;
    StringMatrix para_str_mat;
    read_a_table_file (model_dir + "/" + chr_str + ".model", ROW_DELIM, FIELD_DELIM,'#', 0, {0,1,2,3}, para_str_mat);
    double resolution = atof(para_str_mat[0][0].c_str());
    int chr_length = para_str_mat.size() - 1;  // skip the first line
    ValueList alpha(chr_length, 0), 
              beta_1(chr_length, 0), 
              beta_2(chr_length, 0),
              insulation(chr_length, 0); 
    for (int para = 0; para < chr_length; para++) {
      alpha[para] = atof(para_str_mat[para + 1][0].c_str());
      beta_1[para] = atof(para_str_mat[para + 1][1].c_str());
      beta_2[para] = atof(para_str_mat[para + 1][2].c_str());
      insulation[para] = atof(para_str_mat[para + 1][3].c_str());
    }
    for (int del = 0; del < del_str_mat.size(); del++) {
      if (del_str_mat[del][0].compare(chr_str) == STR_EQ) {
        int start_del = atoi(del_str_mat[del][1].c_str()),
            end_del = atoi(del_str_mat[del][2].c_str());
        if (end_del - start_del >= min_del_length && end_del - start_del <= max_del_length) {
          // Pass the length threshold, the score is zero. Otherwise, the score is -1 and will be filtered
          tad_fusion_score_list[del] = 0;
          int del_start_bin = floor(atoi(del_str_mat[del][1].c_str())/resolution),
              del_end_bin = floor(atoi(del_str_mat[del][2].c_str())/resolution);
          tad_fusion_score_list[del] = estimate_fusion_score(alpha, beta_2, insulation, del_start_bin, del_end_bin, window_size, delta, debug > 0);
          for (int p = 0; p < permutation_num; p++) {
            double r = (rand() % 10001)/10000.0;
            int del_length = del_end_bin - del_start_bin,
                random_start_bin = round(r*(alpha.size() - del_length - 1 - 2*window_size)) + window_size,
                random_end_bin = random_start_bin + del_length;
            permutation_score_mat[p][del] = estimate_fusion_score(alpha, beta_2, insulation, random_start_bin, random_end_bin, window_size, delta, false);
          } 
        }
      }
    }
  }
  ofstream out_file(out_filename.c_str());
  bool first_line = true;
  for (int del = 0; del < del_str_mat.size(); del++) {
    if (tad_fusion_score_list[del] >= 0) {
      if (!first_line)
        out_file << endl;
      first_line = false;
      out_file << del_str_mat[del][0] << "\t" 
               << del_str_mat[del][1] << "\t" 
               << del_str_mat[del][2] << "\t" 
               << tad_fusion_score_list[del];
    }
  }
  out_file.close();
  // print the permutation file
  if (permutation_num > 0) {
    string per_filename = out_filename + ".permutation";
    ofstream per_file(per_filename.c_str());
    for (int p = 0; p < permutation_num; p++) {
      if (p > 0)
        per_file << endl;
      for (int del = 0; del < permutation_score_mat[p].size(); del++) {
        if (del > 0)
          per_file << "\t";
        per_file << permutation_score_mat[p][del];
      }
    }
    per_file.close();
  }
  cout << "Complete!" << endl;
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

