#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

bool compare(string a,string b){
  return a+b<b+a;
}

class Solution {
  public: 
    string largestNumber(vector<int> &num) { 
      vector<string> str; 
      for(int i=0;i<num.size();i++){ 
        ostringstream ss;
        ss<<num[i];
        str.push_back(ss.str());
      } 
      ostringstream ss;
      sort(str.begin(),str.end(),compare);
      for(int i=num.size()-1;i>=0;i--){
        ss<<str[i];
      }
      if(ss.str()[0]=='0') return string("0");
      else return ss.str();
    }
};

int main(){
  std::vector<int> num;
  num.push_back(0);
  num.push_back(0);
  num.push_back(0);
  Solution* sol=new Solution();
  cout<<sol->largestNumber(num);
  return 0;
}
