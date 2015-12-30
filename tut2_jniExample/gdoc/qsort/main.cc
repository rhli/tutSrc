#include <algorithm>
#include <climits>
#include <iostream>
#include <vector>

using namespace std;

void merge(vector<int>& nums, int s, int mid, int e) {
int flen = mid - s + 1;
int slen = e - mid;
vector<int> a = vector<int>(nums.begin() + s, nums.begin() + mid + 1);
vector<int> b = vector<int>(nums.begin() + mid + 1, nums.begin() + e + 1);
a.push_back(INT_MAX);
b.push_back(INT_MAX);

int id1 = 0, id2 = 0;
for (int i = 0; i < slen + flen; i ++) {
if (a[id1] < b[id2]) {
	nums[s + i] = a[id1 ++];
} else {
nums[s + i] = b[id2 ++];
}
}
}

void mergesort(vector<int>& nums, int s, int e) {
if (s >= e) return;
int mid = (s + e) / 2;
mergesort(nums, s, mid);
mergesort(nums, mid + 1, e);
merge(nums, s, mid, e);
}

void mergesort(vector<int>& nums) {
mergesort(nums, 0, nums.size() - 1);
}

int Partition(vector<int>& arr, int s, int e) {
int p = s - 1, i = s;
for (; i < e; i ++) {
if (arr[i] < arr[e]) swap(arr[++ p], arr[i]);
}
swap(arr[++ p], arr[i]);
return p;
}

void Qsort(vector<int>& arr, int s, int e) {
if (s >= e) return;
int m = Partition(arr, s, e);
Qsort(arr, s, m - 1);
Qsort(arr, m + 1, e);
}

void Qsort(vector<int>& arr){
Qsort(arr, 0, arr.size() - 1);
}

int main() {
  vector<int> input = vector<int>({4,1,8,5,2,7,6,3,10});
  mergesort(input);
  for (auto it : input) cout << it << endl;
  return 0;
}
