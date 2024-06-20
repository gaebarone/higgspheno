#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

void concatenate_indices(vector<int> in, vector <int> & out){
  out.insert(out.end(), in.begin(), in.end());
}

#endif