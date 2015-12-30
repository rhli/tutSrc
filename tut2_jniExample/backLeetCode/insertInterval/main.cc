#include <iostream>
#include <algorithm>
#include <vector>
#include <climits>

using namespace std;

struct Interval {
    int start;
    int end;
    Interval() : start(0), end(0) {}
    Interval(int s, int e) : start(s), end(e) {}
};

class Solution {
    int searchStartInd(vector<int> &intList,int start,int end,int sVal){
        if(end==start+1) return end;
        int mid=(start+end)/2;
        int startVal=start==-1?INT_MIN:intList[start];
        int endVal=end==intList.size()?INT_MAX:intList[end];
        if(intList[mid]<sVal){
            return searchStartInd(intList,mid,end,sVal);
        }else{
            return searchStartInd(intList,start,mid,sVal);
        }
    }
    
    int searchEndInd(vector<int> &intList,int start,int end,int eVal){
        if(end==start+1) return start;
        int mid=(start+end)/2;
        int startVal=start==-1?INT_MIN:intList[start];
        int endVal=end==intList.size()?INT_MAX:intList[end];
        if(intList[mid]>eVal){
            return searchEndInd(intList,start,mid,eVal);
        }else{
            return searchEndInd(intList,mid,end,eVal);
        }
    }
    
public:
    vector<Interval> insert(vector<Interval> &intervals, Interval newInterval) {
        vector<int> startEnd;
        for(int i=0;i<intervals.size();i++){
            startEnd.push_back(intervals[i].start);
            startEnd.push_back(intervals[i].end);
        }
        int startInd=searchStartInd(startEnd,-1,startEnd.size(),newInterval.start);
        int endInd=searchEndInd(startEnd,-1,startEnd.size(),newInterval.end);
        int intervalStartInd=startInd/2;
        int intervalEndInd=endInd/2==intervals.size()?endInd/2-1:endInd/2;
        if(intervalStartInd>=intervals.size()){
          intervals.insert(intervals.begin()+intervals.size(),newInterval);
          return intervals;
        }
        if((endInd<0)||(intervals[0].start>newInterval.end)){
          intervals.insert(intervals.begin(),newInterval);
          return intervals;
        }
        cout<<intervalStartInd<<" "<<intervalEndInd<<endl;
        cout<<startInd<<" "<<endInd<<endl;
        Interval nIntval=Interval(min(intervals[intervalStartInd].start,newInterval.start),
            max(intervals[intervalEndInd].end,newInterval.end));
        for(int i=intervalEndInd;i>=intervalStartInd;i--){
            intervals.erase(intervals.begin()+i);
        }
        intervals.insert(intervals.begin()+intervalStartInd,nIntval);
        return intervals;
    }
};

int main(){
  Solution sol=Solution();
  vector<Interval> intervals;
  intervals.push_back(Interval(1,5));
  //intervals.push_back(Interval(7,12));
  sol.insert(intervals,Interval(0,0));
  for(int i=0;i<intervals.size();i++){
    cout<<i<<" "<<intervals[i].start<<" "<<intervals[i].end<<endl;
  }
  return 0;
}


