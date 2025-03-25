#include <bits/stdc++.h>
using namespace std;

const int MIN_LENGTH = 5;  // Minimum length of substrings to consider , to filter substrings that are too short to match the answer
const int BASE = 26, MOD = 1e9 + 7;  // Base and modulus for hash function

struct RepeatInfo {
    int size;
    int pos_in_ref;
    int repeat_count;
    bool inverse;
    string substring;
};

vector<RepeatInfo> results;

// Vectors to store precomputed hash values
vector<long long> pow_hash, prefix_hash;

// Precompute the hash values for the input string s
void precompute_hash(const string &s) {
    int n = s.size();
    pow_hash.resize(n + 1, 1);  // Initialize power hash array
    prefix_hash.resize(n + 1, 0);  // Initialize prefix hash array

    // Compute powers of BASE and prefix hashes for each position
    for (int i = 0; i < n; i++) {
        pow_hash[i + 1] = (pow_hash[i] * BASE) % MOD;
        prefix_hash[i + 1] = (prefix_hash[i] * BASE + (s[i] - 'A')) % MOD;
    }
}

// Return the hash of the substring from index l to r
long long get_hash(int l, int r) {
    return (prefix_hash[r] - prefix_hash[l] * pow_hash[r - l] % MOD + MOD) % MOD;
}


bool compare(int i, int j, const vector<int>& rank, int k, int n) {
    if (rank[i] != rank[j]) return rank[i] < rank[j];
    int ri = (i + k < n) ? rank[i + k] : -1;
    int rj = (j + k < n) ? rank[j + k] : -1;
    return ri < rj;
}

int partition(vector<int>& sa, int low, int high, const vector<int>& rank, int k, int n) {
    int pivot = sa[high];
    int i = low - 1;
    for (int j = low; j < high; j++) {
        // 使用原cmp逻辑比较sa[j]和pivot
        if (compare(sa[j], pivot, rank, k, n)) {
            i++;
            swap(sa[i], sa[j]);
        }
    }
    swap(sa[i + 1], sa[high]);
    return i + 1;
}

void custom_sort(vector<int>& sa, int low, int high, const vector<int>& rank, int k, int n) {
    if (low >= high) return;
    int pivot = partition(sa, low, high, rank, k, n);
    custom_sort(sa, low, pivot - 1, rank, k, n);
    custom_sort(sa, pivot + 1, high, rank, k, n);
}
// Build suffix array for the input string s
vector<int> build_suffix_array(const string &s) {
    int n = s.size();
    vector<int> sa(n), rank(n), temp(n);
    iota(sa.begin(), sa.end(), 0);
    copy(s.begin(), s.end(), rank.begin());

    for (int k = 1; k < n; k *= 2) {
        auto compare = [&](int i, int j) {
            if (rank[i] != rank[j]) return rank[i] < rank[j];
            int ri = (i + k < n) ? rank[i + k] : -1;
            int rj = (j + k < n) ? rank[j + k] : -1;
            return ri < rj;
        };

        // 手动快速排序替代标准库sort
        auto cmp_wrapper = [&](int i, int j) { return compare(i, j); };
        custom_sort(sa, 0, n - 1, rank, k, n);

        temp[sa[0]] = 0;
        for (int i = 1; i < n; i++) {
            temp[sa[i]] = temp[sa[i - 1]] + cmp_wrapper(sa[i - 1], sa[i]);
        }
        swap(rank, temp);
    }
    return sa;
}

// Build LCP (Longest Common Prefix) array for the input string s and its suffix array sa
vector<int> build_lcp(const string &s, const vector<int> &sa) {
    int n = s.size();
    vector<int> rank(n), lcp(n);
    for (int i = 0; i < n; i++) rank[sa[i]] = i;  // Build the rank array from the suffix array

    // Compute the LCP array using the rank array
    for (int i = 0, h = 0; i < n; i++) {
        if (rank[i] > 0) {
            int j = sa[rank[i] - 1];
            while (i + h < n && j + h < n && s[i + h] == s[j + h]) h++;  // Extend the common prefix
            lcp[rank[i]] = h;
            if (h > 0) h--;  // Decrease h for the next suffix
        }
    }
    return lcp;  // Return the LCP array
}

// KMP algorithm to find the first occurrence of the pattern in the reference string
int find_in_reference_kmp(const string &reference, const string &pattern) {
    int m = reference.size(), n = pattern.size();
    if (n == 0) return 0;  // If the pattern is empty, return 0

    vector<int> lps(n, 0);  // Longest Prefix Suffix (LPS) array for the pattern
    // Build the LPS array
    for (int i = 1, len = 0; i < n;) {
        if (pattern[i] == pattern[len]) {
            lps[i++] = ++len;
        } else if (len) {
            len = lps[len - 1];
        } else {
            lps[i++] = 0;
        }
    }

    // Perform the KMP search
    for (int i = 0, j = 0; i < m;) {
        if (reference[i] == pattern[j]) {
            i++; j++;
        }
        if (j == n) return i - j + n - 1;  // Pattern found, return position
        else if (i < m && reference[i] != pattern[j]) {
            if (j) j = lps[j - 1];  // Skip characters in the pattern
            else i++;  // Move to the next character in the reference
        }
    }
    return -1;  // Pattern not found
}

// Precompute repeated substrings and their occurrences in the query string
unordered_map<string, int> precompute_repetitions(const string &query) {
    precompute_hash(query);  // Precompute the hash values for the query string
    unordered_map<string, int> repetition_count;  // Map to store the repetition count of substrings
    int n = query.size();

    // Iterate over possible substring lengths from MIN_LENGTH to n / 2
    for (int len = MIN_LENGTH; len <= n / 2; len++) {
        unordered_map<long long, int> last_pos;  // Map to store the last position of each hash value
        // Iterate over all possible substrings of length 'len'
        for (int i = 0; i + len <= n; i++) {
            long long h = get_hash(i, i + len);  // Get the hash of the current substring
            // Check if the substring has been seen previously
            if (last_pos.count(h) && last_pos[h] == i - len) {
                repetition_count[query.substr(i, len)]++;  // Increment the repetition count
            }
            last_pos[h] = i;  // Update the last position of the current hash
        }
    }
    return repetition_count;  // Return the map of repeated substrings and their counts
}

// Function to find and print repeated substrings in both reference and query strings
void find_repeated_substrings(const string &reference, const string &query, bool inv) {
    string s = query + "$";
    vector<int> sa = build_suffix_array(s);
    vector<int> lcp = build_lcp(s, sa);
    auto repetition_count = precompute_repetitions(query);
    unordered_set<string> seen;

    for (int i = 1; i < lcp.size(); i++) {
        if (lcp[i] >= MIN_LENGTH) {
            string candidate = query.substr(sa[i], lcp[i]);
            if (seen.find(candidate) != seen.end()) continue;

            int occurrences = repetition_count.count(candidate) ? repetition_count[candidate] : 0;
            if (occurrences < 1) continue;

            int ref_pos = find_in_reference_kmp(reference, candidate);
            if (ref_pos != -1) {
                seen.insert(candidate);
                int reported_occurrences = occurrences + 1;
                int adjusted_pos = inv ? reference.length() - ref_pos - 1 + lcp[i] : ref_pos + 1;
                results.push_back({lcp[i], adjusted_pos, reported_occurrences, inv, candidate});
            }
        }
    }
}

int find_next_occurrence(string query, string substring, int start_pos) {
    size_t found_pos = query.find(substring, start_pos);
    return (found_pos != string::npos) ? found_pos : -1; // Return 1-based index, or -1 if not found
}

// Reverse the reference string and replace A/T/C/G with T/A/G/C
string rev_ref(string &str) {
    string ret = "";
    for (int i = str.length() - 1; i >= 0; --i) {
        char c = str[i];
        ret += (c == 'A' ? 'T' : c == 'T' ? 'A' : c == 'C' ? 'G' : 'C');
    }
    return ret;
}


int main() {
    string reference = "CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA";
    
    string ref2 = rev_ref(reference);
    string query =     "CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA";
    
    
    find_repeated_substrings(reference, query, 0);
    find_repeated_substrings(ref2, query, 1);
    
    for (const auto &res : results) {
        cout << "Repeat size: " << res.size << endl;
        cout << "POS in REF: " << res.pos_in_ref << endl;
        int repeat=(find_next_occurrence(query, res.substring, res.pos_in_ref)==res.pos_in_ref)? (res.repeat_count-1): res.repeat_count ;
        cout << "Repeat count: " << repeat<< endl;
        cout << "Inverse: " << (res.inverse ? "YES" : "NO") << endl;
        //cout <<"next_occurrence:"<<find_next_occurrence(query, res.substring, res.pos_in_ref)<<endl;
        //cout << "Substring: " << res.substring << endl;
        cout << "------------------------\n";
    }
    
    return 0;
}

